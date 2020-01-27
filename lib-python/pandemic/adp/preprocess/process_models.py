import os, shutil, collections
from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure
from bamboo.common.logs import Log
from bamboo.common.path import easy_directory

from giant.xray.crystal import CrystalSummary

def get_resolution_from_model_pdb_input(model=None, mtz_filename=None):
    return model.input.resolution()

def get_resolution_from_model_mtz(model=None, mtz_filename=None):
    return CrystalSummary.from_mtz(mtz_filename).high_res


class ProcessInputModelsTask:


    get_resolution = staticmethod(get_resolution_from_model_pdb_input)

    def __init__(self,
            output_directory,
            dataset_selection_params,
            weights_params,
            analysis_params,
            table_one_options,
            look_for_reflection_data,
            copy_reflection_data_to_output_folder,
            check_column_labels,
            verbose = False,
            log = None,
            ):
        if log is None: log = Log()
        has_reflection_data = False
        need_res_info = (dataset_selection_params.max_resolution is not None) or (weights_params.dataset_weights != 'one')
        table_one_cols = set(table_one_options.column_labels.split(',') + [table_one_options.r_free_label])

        adopt_init_args(self, locals())

    def run(self,
            models,
            ):

        self.log.subheading('Processing input models')

        # Use the first hierarchy as the reference
        self.reference_model = models[0]
        self.master_h = models[0].hierarchy.deep_copy()
        # Store input files in hash for convenience
        self.input_structures = collections.OrderedDict()
        self.input_reflection_data = collections.OrderedDict()

        errors = []

        for i_m, m in enumerate(models):
            # Check that all of the structures are the same
            if not self.master_h.is_similar_hierarchy(m.hierarchy):
                errors.append(Failure("Structures are not all the same. Model {}. File: {}".format(i_m, m.filename)))
                continue

            # Directory for this model
            m.directory = easy_directory(os.path.join(self.output_directory, m.tag))

            # filepaths for copying input files to output folder
            m_pdb = self.input_structures.setdefault(m.tag, os.path.join(m.directory, m.tag+'.input.pdb'))
            m_mtz = self.input_reflection_data.setdefault(m.tag, os.path.join(m.directory, m.tag+'.input.mtz'))

            # Write the input pdb to the output folder (without the header as this can affect R-factors calculations)
            m.hierarchy.write_pdb_file(
                    m_pdb,
                    crystal_symmetry=(m.crystal_symmetry if hasattr(m, 'crystal_symmetry') else None),
                    )
            assert os.path.exists(m_pdb), 'PDB does not exist: {}'.format(m_pdb)

            # Check MTZ exists if requested
            if self.look_for_reflection_data:
                # Look for same filename as the pdb
                mtz_file = m.filename.replace('.pdb','.mtz')
                # Check all or none of the datasets have data
                if (os.path.exists(mtz_file)):
                    # We have found reflection data -- update flags
                    if self.has_reflection_data is False:
                        # Update flag
                        self.has_reflection_data = True
                        # Update function to extract resolution from model
                        self.get_resolution = get_resolution_from_model_mtz
                    # Link the input mtz to the output folder
                    if self.copy_reflection_data_to_output_folder:
                        shutil.copy(mtz_file, m_mtz)
                    else:
                        rel_symlink(mtz_file, m_mtz)
                    assert os.path.exists(m_mtz), 'MTZ does not exist: {}'.format(m_mtz)
                elif (self.has_reflection_data is True):
                    errors.append(Sorry('MTZ files have been found for some datasets but not others.\n'\
                                        'To disable searching for MTZs, set input.look_for_reflection_data=False.\n'\
                                        'Missing File: {}'.format(mtz_file)))
                    continue

            # Extract the resolution of the dataset
            m.resolution = self.get_resolution(model=m, mtz_filename=m_mtz)

            # Check we have resolution information, and/or that correct columns are present in the MTZ file
            if self.has_reflection_data:
                # Extract crystal information from the MTZ
                cs = CrystalSummary.from_mtz(m_mtz)
                # Check for columns
                if self.check_column_labels and self.table_one_cols.difference(cs.column_labels):
                    errors.append(Failure("MTZ {} does not contain the correct columns.".format(m.filename) + \
                                          "\n\tLooking for: {}".format(','.join(self.table_one_cols)) + \
                                          "\n\tMTZ contains: {}".format(','.join(cs.column_labels)) + \
                                          "\n\tCan't find: {}".format(','.join(self.table_one_cols.difference(cs.column_labels))) + \
                                          "\nChange required columns with {}".format("table_ones_options.column_labels or table_ones_options.r_free_label")))

                # Check for high res information
                if (m.resolution is None) and self.need_res_info:
                    errors.append(Failure("MTZ does not contain resolution information: {}\n".format(m_mtz)))
                    continue
            else:
                # Check for high res information
                if (m.resolution is None) and self.need_res_info:
                    errors.append(Failure("PDB does not contain resolution information in the REMARK 3 records: {}\n".format(m_pdb) + \
                                          "This is normally present in structures that have come from refinement.\n" + \
                                          "You can also remove the high-resolution cutoff for dataset masking (optimisation.max_resolution=None),\n" + \
                                          "or use equal weighting for each dataset (optimisation.dataset_weights=one)"))
                    continue

        # Check for errors
        if errors:
            self.log.heading('Errors found while processing input files')
            for e in errors:
                self.log(str(e))
                self.log('')
            self.log.bar(False, True)
            raise Sorry("Some structures are invalid or missing information.")

        # Report dataset resolutions
        self.log('\nDataset resolutions (extracted from {} files):'.format('MTZ' if self.has_reflection_data else 'PDB'))
        for m in models:
            self.log('> {}: {}'.format(m.tag, m.resolution))

        # Get and print optimisation
        # Update parameter flags
        if not self.has_reflection_data:
            self.log.bar(True, False)
            self.log('Reflection data has not been found in the input folder(s) -- updating settings.')
            self.log('(Added reflection data must have the same filename as input structures')
            self.log('but with different extension, e.g. <filename>.pdb & <filename>.mtz)')
            self.log('Without reflection data, it is not possible to refine the fitted model or recalculate R-factors.')
            self.log('Setting analysis.refine_output_structures = False')
            self.analysis_params.refine_output_structures = False
            self.log('Setting analysis.calculate_r_factors = False')
            self.analysis_params.calculate_r_factors = False
            #self.log('Setting analysis.calculate_electron_density_metrics = False')
            #self.analysis_params.calculate_electron_density_metrics = False
            self.log.bar()



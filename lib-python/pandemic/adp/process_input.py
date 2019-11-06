import os, shutil, collections
import numpy

from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure
from bamboo.common.path import easy_directory
from bamboo.common.command import check_programs_are_available
from giant.dataset import CrystallographicModel, AtomicModel
from giant.xray.crystal import CrystalSummary

import math
EIGHT_PI_SQ = 8*math.pi*math.pi


def get_resolution_from_model_pdb_input(model=None, mtz_filename=None):
    return CrystalSummary.from_pdb(pdb_input=model.input).high_res

def get_resolution_from_model_mtz(model=None, mtz_filename=None):
    return CrystalSummary.from_mtz(mtz_filename).high_res

def validate_parameters(params, log=None):

    log.heading('Processing input parameters')

    if len(params.input.pdb) == 0:
        raise Sorry('No structures have been provided for analysis (input.pdb=[...])')

    if params.settings.debug:
        log('DEBUG is turned on -- setting verbose=True')
        params.settings.verbose = True

    # Output images
    if params.output.images.all:
        log('params.output.images.all = {}'.format(params.output.images.all))
        # ----------------
        new = 'all'
        log('Setting params.output.images.pymol to {}'.format(new))
        params.output.images.pymol = new
        # ----------------
        new = True
        log('Setting params.output.images.distributions to {}'.format(new))
        params.output.images.distributions = new
        # ----------------
        log.bar()
    
    # Special settings if only one model is loaded
    if (len(params.input.pdb) == 1):
        log.bar(True, False)
        log('One file provided for analysis -- updating settings')
        # TURN THIS OFF -- there may be advantages for fitting a simple hierarchical model to one structure (e.g. low resolution)
        #log('For one model, it is not neccessary to refine the fitted model or recalculate R-factors (or even look for reflection data)')
        #log('Setting analysis.refine_output_structures = False')
        #params.analysis.refine_output_structures = False
        #log('Setting analysis.calculate_r_factors = False')
        #params.analysis.calculate_r_factors = False
        #log('Setting analysis.calculate_electron_density_metrics = False')
        #params.analysis.calculate_electron_density_metrics = False
        #log('Setting input.look_for_reflection_data = False')
        #params.input.look_for_reflection_data = False
        log('Setting fitting.optimisation.dataset_weights = one')
        params.optimisation.weights.dataset_weights = 'one'
        log.bar()

    try:
        if params.analysis.refine_output_structures is True:
            if params.refinement.program == 'phenix':
                message = 'phenix.refine is required when analysis.refine_output_structures=True and refinement.program=phenix'
                check_programs_are_available(['phenix.refine'])
            elif params.refinement.program == 'refmac':
                message = 'refmac5 is required when analysis.refine_output_structures=True and refinement.program=refmac'
                check_programs_are_available(['refmac5'])
            else:
                raise Sorry('Must select refinement.program when analysis.refine_output_structures=True')
        if params.analysis.calculate_r_factors is True:
            message = 'phenix.table_one is required when analysis.calculate_r_factors is True'
            check_programs_are_available(['phenix.table_one'])
        #if params.analysis.calculate_electron_density_metrics:
        #    message = 'edstats (script name "edstats.pl") is required when analysis.calculate_electron_density_metrics is True'
        #    check_programs_are_available(['edstats.pl'])
        if params.output.images.pymol is not None:
            message = 'pymol is required when output.images.pymol is not set to "none"'
            check_programs_are_available(['pymol'])
    except Exception as e:
        #log(e)
        raise


class ModelLoader:


    def __init__(self,
            model_type = 'crystallographic',
            labelling = 'foldername',
            verbose = False,
            log = None,
            ):
        adopt_init_args(self, locals())

        # Store appropriate model labelling function
        from bamboo.common.path import foldername, filename
        if labelling == 'foldername':
            self.label_func = foldername
        elif labelling == 'filename':
            self.label_func = filename
        else:
            raise Exception('Invalid labelling function: {}'.format(labelling))

        # Store appropriate model class
        if model_type == "crystallographic":
            self.model_class = CrystallographicModel
        else:
            self.model_class = AtomicModel

    def __call__(self,
            pdb_files,
            ):

        log = self.log

        # Load input structures
        log.subheading('Building model list -- {} files'.format(len(pdb_files)))
        models = []
        for f in pdb_files:
            m = self.model_class.from_file(f)
            l = self.label_func(f)
            if not l:
                if len(pdb_files) == 1:
                    log('No label created for label function: {}'.format(self.labelling))
                    log('Trying to label by filename instead')
                    self.labelling = 'filename'
                    from bamboo.common.path import filename
                    self.label_func = filename
                    l = self.label_func(f)
                if not l:
                    raise Sorry('No label generated using labelling function "{}"\n\tLabel {}\n\tFile {}'.format(self.labelling, l, f))
            m.label(tag=l)
            models.append(m)

        # Check for duplicate labels
        all_labels = [m.tag for m in models]
        unq_labels = sorted(set(all_labels))
        if len(unq_labels) != len(models):
            counts = [(l, all_labels.count(l)) for l in unq_labels]
            dups = ['{} (found {} times)'.format(l,c) for l,c in counts if c>1]
            raise Sorry('Duplicate labels generated for models: \n\t{}'.format('\n\t'.join(dups)))

        # Sort models for convenience
        models = sorted(models, key=lambda m: m.tag)
        log('{} models loaded'.format(len(models)))

        return models


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


class SelectOptimisationDatasetsTask:


    def __init__(self,
            max_resolution,
            max_datasets,
            sort_datasets_by = 'resolution',
            random_seed = 0,
            verbose = False,
            log = None,
            ):
        adopt_init_args(self, locals())

    def run(self,
            models,
            resolutions = None,
            ):

        if resolutions is not None:
            assert len(models) == len(resolutions)

        if (self.max_resolution is None):
            opt_datasets = range(len(models))
        else:
            assert resolutions is not None, 'no resolutions provided'
            opt_datasets = numpy.where(resolutions < self.max_resolution)[0]

        if len(opt_datasets) == 0:
            raise Sorry('No datasets selected for optimisation (e.g. above resolution cutoff: {})'.format(self.max_resolution))

        if self.sort_datasets_by == 'resolution':
            assert resolutions is not None, 'no resolutions provided'
            opt_datasets = sorted(opt_datasets, key=lambda i: resolutions[i])
        elif self.sort_datasets_by == 'name':
            opt_datasets = sorted(opt_datasets, key=lambda i: self.models[i].tag)
        elif self.sort_datasets_by == 'random':
            self.log('Setting random seed: {}'.format(self.random_seed))
            numpy.random.seed(self.random_seed)
            opt_datasets = numpy.random.permutation(opt_datasets)

        self.log('After reordering:')
        for i_m in opt_datasets:
            self.log('\t{}: {}'.format(i_m, self.models[i_m].tag))

        # Limit the number of datasets for optimisation
        if (self.max_datasets is not None) and (len(opt_datasets) > self.max_datasets):
            self.log('\nLimiting list of datasets for TLS optimisation to {} datasets'.format(self.max_datasets))
            opt_datasets = opt_datasets[:self.max_datasets]

        assert len(selection) > 0, 'no datasets selected for optimisation with resolution cutoff: {}'.format(self.max_resolution)

        self.iselection = opt_datasets
        return iselection


class ExtractAndProcessModelUijsTask:


    def __init__(self,
            expected_disorder_model,
            verbose = False,
            log = None,
            ):
        if log is None: log = Log()
        actual_disorder_model = expected_disorder_model
        adopt_init_args(self, locals())

    @staticmethod
    def extract_uijs_from_models(
            models_atoms,
            overall_atom_mask = None,
            ):
        # Extract all atoms from all datasets
        assert len(models_atoms) > 0, 'No models have been used?!'
        assert len(models_atoms[0]) > 0, 'No atoms have been extracted from models?'
        # Extract uij and xyz from all datasets (and only select atoms we're interested in)
        model_uij = numpy.array([a.extract_uij() for a in models_atoms])
        # Filter on input mask
        if overall_atom_mask is not None:
            model_uij = model_uij[:,overall_atom_mask]
        return model_uij

    def determine_disorder_model_from_mask(self,
            isotropic_mask,
            ):
        # Initialise disorder model as expected
        actual_disorder_model = self.expected_disorder_model
        # Check that the isotropic mask is as expected
        if isotropic_mask.selection.all_eq(True):
            if (self.expected_disorder_model != 'isotropic'):
                raise Sorry('All atoms selected for fitting have isotropic ADPs and input.input_uij_model is currently set to "{}".\nOptions:'.format(self.expected_disorder_model) + \
                        '\n   1) Re-refine structures with anisotropic ADPs.' + \
                        '\n   2) Set input.input_uij_model=isotropic to allow use of isotropic ADPs.')
        elif isotropic_mask.selection.all_eq(False):
            if (self.expected_disorder_model == 'isotropic'):
                raise Sorry('input.input_uij_model is set to "{}" but all atoms selected have anisotropic uij.'.format(self.expected_disorder_model))
            if (self.expected_disorder_model == 'mixed'):
                self.log('Input disorder model is set to "{}", but all input atoms have anisotropic Uij.'.format(self.expected_disorder_model))
                actual_disorder_model = 'anisotropic'
                self.log('Updated disorder_model to "{}"'.format(actual_disorder_model))
        else:
            if (self.expected_disorder_model != 'mixed'):
                raise Sorry('Some atoms for fitting (but not all) have anisotropically-refined ADPs and input.input_uij_model is currently set to "{}".\nOptions:'.format(self.expected_disorder_model) + \
                        '\n   1) Re-refine structures so that all atoms selected have {} ADPs and set input.input_uij_model={}.'.format(self.expected_disorder_model, self.expected_disorder_model) + \
                        '\n   2) Change levels.overall_selection so that atoms with anisotropic/isotropic ADPs are not selected.' + \
                        '\n   3) Change input.input_uij_model to "mixed".')

        return actual_disorder_model

    @staticmethod
    def generate_uij_for_isotropic_atoms(
            models_atoms,
            model_uij,
            isotropic_mask,
            overall_atom_mask = None,
            ):

        # Create copy
        import copy
        model_uij = copy.deepcopy(model_uij)

        # Find atoms that are isotropic
        assert (model_uij[:,isotropic_mask.selection,:] == -1).all()
        # Set atoms to zeros
        model_uij[:,isotropic_mask.selection,:] = 0.0
        # Extract B-values and convert to Uijs
        model_b = numpy.array([a.extract_b()/EIGHT_PI_SQ for a in models_atoms])
        # Apply overall mask
        if overall_atom_mask is not None:
            model_b = model_b[:,overall_atom_mask]
        # Extract values for isotropic atoms
        model_b_iso = model_b[:,isotropic_mask.selection]
        # Set Uij values of isotropic atoms
        model_uij[:,isotropic_mask.selection,0] = model_b_iso
        model_uij[:,isotropic_mask.selection,1] = model_b_iso
        model_uij[:,isotropic_mask.selection,2] = model_b_iso

        return model_uij

    def run(self,
            models,
            overall_atom_mask = None,
            ):

        log = self.log

        log.subheading('Extracting Uijs from models')

        # Extract atoms from models
        models_atoms = [m.hierarchy.atoms() for m in models]

        model_uij = self.extract_uijs_from_models(
                models_atoms = models_atoms,
                overall_atom_mask = overall_atom_mask,
                )

        # Create a mask object for isotropic atoms
        from pandemic.adp.utils import UijIsotropicMask
        isotropic_mask = UijIsotropicMask.from_uij_array(model_uij, axis=(0,2))

        # Determine the disorder model
        actual_disorder_model = self.determine_disorder_model_from_mask(isotropic_mask)

        # Delete the isotropic mask if not needed
        #if isotropic_mask.selection.all_eq(False):
        #    isotropic_mask = None

        # Extract b-values for isotropic atoms and convert to isotropic Uij
        if not isotropic_mask.all_anisotropic:
            model_uij = self.generate_uij_for_isotropic_atoms(
                    models_atoms = models_atoms,
                    model_uij = model_uij,
                    isotropic_mask = isotropic_mask,
                    overall_atom_mask = overall_atom_mask,
                    )

        # All -1 values should now have been removed from input array
        assert not (model_uij==-1).all(axis=2).any()

        self.result = group_args(
                disorder_model = actual_disorder_model,
                model_uij = model_uij,
                isotropic_mask = isotropic_mask,
                )
        self.show_summary()
        return self.result

    def show_summary(self):
        log = self.log
        r = self.result

        log.bar()
        log('Expected disorder model was {}'.format(self.expected_disorder_model))
        log('Actual disorder model is {}'.format(r.disorder_model))
        log.bar()
        log('Anisotropic atoms: {}'.format(r.isotropic_mask.n_anisotropic))
        log('Isotropic atoms: {}'.format(r.isotropic_mask.n_isotropic))
        log.bar()


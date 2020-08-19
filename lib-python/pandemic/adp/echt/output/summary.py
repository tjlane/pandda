import giant.logs as lg
logger = lg.getLogger(__name__)

import os, glob, collections
from libtbx import adopt_init_args

import numpy
from scitbx.array_family import flex

from giant.paths import easy_directory
from giant.structure.pymol import auto_chain_images
from giant.structure.uij import \
        uij_to_b, \
        calculate_uij_anisotropy_ratio
from pandemic.adp.utils import show_file_dict
from pandemic.adp.hierarchy.utils import PartitionBordersFactory


class ModelSummaryOutput:

    def __init__(self,
        output_files,
        level_b_factor_statistics_table,
        ):

        adopt_init_args(self, locals())


class WriteEchtModelSummary:

    all_levels_uijs_profiles_png_prefix   = 'all-levels-uij-profile'
    all_levels_uijs_anisotropy_png_prefix = 'all-levels-uij-anisotropy'

    level_uijs_profile_png_prefix    = 'uij-profile-level_{}'
    level_uijs_anisotropy_png_prefix = 'uij-anisotropy-level_{}'

    target_uijs_pdb         = 'average_uijs_input.pdb'
    output_uijs_pdb         = 'average_uijs_output.pdb'

    level_uijs_pdb          = 'level_{}-all.pdb'
    level_uijs_rescaled_pdb = 'level_{}-all-rescaled.pdb'
    level_uijs_by_mode_pdb  = 'level_{}-mode_{}.pdb'

    level_uijs_pymol_by_chain_png_prefix   = 'level_{}'

    level_tls_matrices_csv   = 'tls_matrices_level_{:04d}.csv'
    level_tls_amplitudes_csv = 'tls_amplitudes_level_{:04d}.csv'
    level_tls_origins_csv    = 'tls_origins_level_{:04d}.csv'

    level_tls_amplitudes_png = 'tls-model-amplitudes-level-{}-group-{}.png'

    pymol_script_py = 'pymol_script.py'

    show_file_dict = show_file_dict

    def __init__(self,
        output_directory,
        pymol_images = None,
        distribution_images = False,
        ):
        adopt_init_args(self, locals())
        self.check_and_create_output_folders()

    def check_and_create_output_folders(self):
        """
        Checks output folders exist and creates if necessary.
        Run at beginning of calling function so that folders are updated relative to output_directory.
        """
        easy_directory(self.output_directory)
        self.csv_directory = easy_directory(os.path.join(self.output_directory, 'csvs'))
        self.pdb_directory = easy_directory(os.path.join(self.output_directory, 'structures'))
        if (self.pymol_images is not None):
            self.pymol_directory = easy_directory(os.path.join(self.output_directory, 'pymol_images'))

    def filepath(self, filename, directory=None):
        if directory is None:
            directory = self.output_directory
        return os.path.join(directory, filename)

    def __call__(self,
        overall_atom_mask,
        reference_model,
        model_object,
        isotropic_mask,
        uij_target,
        level_group_array,
        plotting_object,
        ):

        self.check_and_create_output_folders()

        self.plot = plotting_object

        # Ensure flex
        overall_atom_mask = flex.bool(overall_atom_mask)

        # Model values
        average_uijs_by_mode, \
        average_uijs_by_level, \
        atomic_uijs = self.extract_average_uij_values(model_object)

        # input/output values
        average_uijs_target = uij_target.mean(axis=0)
        average_uijs_output = average_uijs_by_level.sum(axis=0) + atomic_uijs

        # Structure-writer
        structure_factory = PartitionBordersFactory(master_h=reference_model.hierarchy)

        output_files = collections.OrderedDict()

        of = self.write_average_target_output_structures(
            average_uijs_target = average_uijs_target, # isotropic mask not necessary here because determined by input anyway...
            average_uijs_output = isotropic_mask(average_uijs_output),
            overall_atom_mask = overall_atom_mask,
            structure_factory = structure_factory,
            )
        self.show_file_dict(of)
        output_files.update(of)

        of = self.write_average_level_structures(
            average_mode_uijs = average_uijs_by_mode,
            average_level_uijs = average_uijs_by_level,
            atomic_level_uijs = atomic_uijs,
            overall_atom_mask = overall_atom_mask,
            isotropic_mask = isotropic_mask,
            structure_factory = structure_factory,
            model_object = model_object,
            )
        self.show_file_dict(of)
        output_files.update(of)

        of = self.write_average_structure_plots(
            level_group_array = level_group_array,
            average_mode_uijs = average_uijs_by_mode,
            average_level_uijs = average_uijs_by_level,
            atomic_level_uijs = atomic_uijs,
            average_target_uijs = average_uijs_target,
            overall_atom_mask = overall_atom_mask,
            structure_factory = structure_factory,
            model_object = model_object,
            )
        self.show_file_dict(of)
        output_files.update(of)

        of = self.write_tls_matrix_and_amplitude_csvs(
            model_object = model_object,
            )
        self.show_file_dict(of)
        output_files.update(of)

        of = self.make_pymol_images_of_average_structures(
            level_uijs_pdbs_dict = output_files['level_uijs_pdb'],
            model_object = model_object,
            reference_hierarchy = reference_model.hierarchy,
            overall_atom_mask = overall_atom_mask,
            )
        self.show_file_dict(of)
        output_files.update(of)

        of = {'pymol_script' : self.make_pymol_script_of_average_structures(file_dict=output_files)}
        self.show_file_dict(of)
        output_files.update(of)

        level_b_factor_statistics_table = self.calculate_various_average_b_statistics(
            level_b_values = [uij_to_b(u) for u in list(average_uijs_by_level)+[atomic_uijs] ],
            level_names = model_object.all_level_names,
            structure_factory = structure_factory,
            overall_atom_mask = overall_atom_mask,
        )

        self.plot = None

        self.result = ModelSummaryOutput(
            output_files = output_files,
            level_b_factor_statistics_table = level_b_factor_statistics_table,
        )

        return self.result

    def extract_average_uij_values(self,
        model_object,
        ):
        """Extract average values"""
        av_uijs_by_mode = model_object.uijs_tls().mean(axis=2)  # (n_tls, n_mode, *n_dataset*, n_atom, 6)
        av_uijs_by_level = av_uijs_by_mode.sum(axis=1)          # (n_tls, *n_mode*, n_atom, 6)
        atomic_uijs = numpy.array(model_object.adp_values)    # (n_atoms, 6)
        return av_uijs_by_mode, av_uijs_by_level, atomic_uijs

    def write_average_target_output_structures(self,
        average_uijs_target,
        average_uijs_output,
        overall_atom_mask,
        structure_factory,
        ):

        file_dict = collections.OrderedDict()

        uij = average_uijs_target
        m_h = structure_factory.custom_copy(
                uij=uij,
                iso=uij_to_b(uij),
                mask=overall_atom_mask,
                blank_copy=True,
                )
        m_f = self.filepath(self.target_uijs_pdb, self.pdb_directory)
        m_h.write_pdb_file(m_f)
        file_dict['target_uijs_pdb'] = m_f

        uij = average_uijs_output
        m_h = structure_factory.custom_copy(
                uij=uij,
                iso=uij_to_b(uij),
                mask=overall_atom_mask,
                blank_copy=True,
                )
        m_f = self.filepath(self.output_uijs_pdb, self.pdb_directory)
        m_h.write_pdb_file(m_f)
        file_dict['output_uijs_pdb'] = m_f

        return file_dict

    def write_average_level_structures(self,
        average_mode_uijs,
        average_level_uijs,
        atomic_level_uijs,
        overall_atom_mask,
        isotropic_mask,
        structure_factory,
        model_object,
        ):

        file_dict = collections.OrderedDict()
        # Iterate through the levels and plot TLS contributions for each mode for each model for each level
        for i_level, level_name in enumerate(model_object.tls_level_names):
            for i_tls in xrange(average_mode_uijs.shape[1]):
                uij = average_mode_uijs[i_level, i_tls]
                m_h = structure_factory.custom_copy(
                        uij=isotropic_mask(uij),
                        iso=uij_to_b(uij),
                        mask=overall_atom_mask,
                        blank_copy=True,
                        )
                m_f = self.filepath(self.level_uijs_by_mode_pdb.format(i_level+1, i_tls+1), self.pdb_directory)
                m_h.write_pdb_file(m_f)
                file_dict.setdefault('level_uijs_by_mode_pdb',collections.OrderedDict())[(level_name,i_tls)] = m_f
            # Write one for each level
            uij = average_level_uijs[i_level]
            m_h = structure_factory.custom_copy(
                    uij=isotropic_mask(uij),
                    iso=uij_to_b(uij),
                    mask=overall_atom_mask,
                    blank_copy=True,
                    )
            m_f = self.filepath(self.level_uijs_pdb.format(i_level+1), self.pdb_directory)
            m_h.write_pdb_file(m_f)
            file_dict.setdefault('level_uijs_pdb',collections.OrderedDict())[level_name] = m_f
        # Write atomic level
        uij = atomic_level_uijs
        m_h = structure_factory.custom_copy(
                uij=isotropic_mask(uij),
                iso=uij_to_b(uij),
                mask=overall_atom_mask,
                blank_copy=True,
                )
        m_f = self.filepath(self.level_uijs_pdb.format(model_object.n_levels), self.pdb_directory)
        m_h.write_pdb_file(m_f)
        file_dict.setdefault('level_uijs_pdb',collections.OrderedDict())[model_object.adp_level_name] = m_f
        return file_dict

    def write_average_structure_plots(self,
        level_group_array,
        average_mode_uijs,
        average_level_uijs,
        atomic_level_uijs,
        average_target_uijs,
        overall_atom_mask,
        structure_factory,
        model_object,
        ):

        atom_sel = overall_atom_mask
        n_modes = model_object.n_modes

        file_dict = collections.OrderedDict()

        # ------------------------
        # Write profiles for each tls level
        # ------------------------
        for i_level, level_name in enumerate(model_object.tls_level_names):
            # ------------------------
            # Get boundaries of groups in this level
            # ------------------------
            boundaries = structure_factory.partition_boundaries(
                atom_labels = level_group_array[i_level],
                mask = atom_sel,
                ).select(atom_sel)
            # Don't make boundaries if too many lines to plot
            if sum(boundaries.atoms().extract_b()) > 0.5*len(list(boundaries.residue_groups())):
                boundaries = None
            # ------------------------
            # Write stacked bar plot of TLS for this level
            # ------------------------
            uijs = average_mode_uijs[i_level]
            filenames_prefix = self.filepath(self.level_uijs_profile_png_prefix.format(i_level+1))
            output_hash = self.plot.stacked_bar(
                prefix        = filenames_prefix,
                hierarchies   = [structure_factory.custom_copy(iso=uij_to_b(u), mask=atom_sel).select(atom_sel) for u in uijs],
                legends       = ['TLS (Mode {})'.format(i+1) for i in xrange(n_modes)],
                title         = 'TLS contributions - Level {} ({})'.format(i_level+1, level_name),
                v_line_hierarchy = boundaries,
                colour_indices = [float(i_level)+(float(i_mode)/float(n_modes)) for i_mode in range(n_modes)],
                )
            if not output_hash:
                logger.warning('no plots have been generated! ({})'.format(filenames_prefix+'*'))
            else:
                # Store in output dictionary
                file_dict.setdefault('level_uijs_profiles_png',collections.OrderedDict())[level_name] = output_hash
            # ------------------------
            # Write graph of anisotropy
            # ------------------------
            ani = 1.0 - calculate_uij_anisotropy_ratio(uij=average_level_uijs[i_level])
            filenames_prefix = self.filepath(self.level_uijs_anisotropy_png_prefix.format(i_level+1))
            graph_title = ('Anisotropy of Level {} ({})'+ \
            '\nfully isotropic -> 0 (spheres) ' + \
            '\nfully anisotropic -> 1 (lines/disks)').format(i_level+1, level_name)
            output_hash = self.plot.stacked_bar(
                prefix        = filenames_prefix,
                hierarchies   = [structure_factory.custom_copy(iso=ani, mask=atom_sel).select(atom_sel)],
                legends       = ['Anisotropy'],
                title         = graph_title,
                y_lim         = (0.0,1.05),
                y_lab         = 'Anisotropy of Uij ($1 - \\frac{E_{min}}{E_{max}}$)',
                colours       = ['grey'],
                )
            if not output_hash:
                logger.warning('no plots have been generated! ({})'.format(filenames_prefix+'*'))
            else:
                # Store in output dictionary
                file_dict.setdefault('level_uijs_anisotropy_png',collections.OrderedDict())[level_name] = output_hash

        # ------------------------------------------------------------------------------------------------------------------------- #
        #                                                                                                                           #
        # ------------------------------------------------------------------------------------------------------------------------- #

        # ------------------------
        # Write profiles for atomic level
        # ------------------------
        uijs = atomic_level_uijs
        filenames_prefix = self.filepath(self.level_uijs_profile_png_prefix.format(model_object.n_levels))
        output_hash = self.plot.stacked_bar(
            prefix         = filenames_prefix,
            hierarchies    = [structure_factory.custom_copy(iso=uij_to_b(uijs), mask=atom_sel).select(atom_sel)],
            legends        = ['{!s:12}'.format(model_object.adp_level_name)],
            title          = 'Uij Profile of {} Level'.format(model_object.adp_level_name.title()),
            colour_indices = [model_object.n_levels-1],
            )
        if not output_hash:
            logger.warning('no plots have been generated! ({})'.format(filenames_prefix+'*'))
        else:
            # Store in output dictionary
            file_dict.setdefault('level_uijs_profiles_png',collections.OrderedDict())[model_object.adp_level_name] = output_hash

        # ------------------------
        # Write graph of anisotropy for atomic level
        # ------------------------
        ani = 1.0 - calculate_uij_anisotropy_ratio(uij=atomic_level_uijs)
        filenames_prefix = self.filepath(self.level_uijs_anisotropy_png_prefix.format(model_object.n_levels))
        graph_title = ('Anisotropy of Level {} ({})'+ \
        '\nfully isotropic -> 0 (spheres) ' + \
        '\nfully anisotropic -> 1 (lines/disks)').format(model_object.n_levels, model_object.adp_level_name)
        output_hash = self.plot.stacked_bar(
            prefix        = filenames_prefix,
            hierarchies   = [structure_factory.custom_copy(iso=ani, mask=atom_sel).select(atom_sel)],
            legends       = ['Anisotropy'],
            title         = graph_title,
            y_lim         = (0.0,1.05),
            y_lab         = 'Anisotropy of Uij ($1 - \\frac{E_{min}}{E_{max}}$)',
            colours       = ['grey'],
            )
        if not output_hash:
            logger.warning('no plots have been generated! ({})'.format(filenames_prefix+'*'))
        else:
            # Store in output dictionary
            file_dict.setdefault('level_uijs_anisotropy_png',collections.OrderedDict())[model_object.adp_level_name] = output_hash

        # ------------------------------------------------------------------------------------------------------------------------- #
        #                                                                                                                           #
        # ------------------------------------------------------------------------------------------------------------------------- #

        # ------------------------
        # Write stacked profiles for all levels
        # ------------------------
        uijs = numpy.append(average_level_uijs, [atomic_level_uijs], axis=0)
        assert uijs.shape == (model_object.n_levels, model_object.n_atoms, 6)
        filenames_prefix = self.filepath(self.all_levels_uijs_profiles_png_prefix)
        output_hash = self.plot.stacked_bar(
            prefix        = filenames_prefix,
            hierarchies   = [structure_factory.custom_copy(iso=uij_to_b(u), mask=atom_sel).select(atom_sel) for u in uijs],
            legends       = model_object.all_level_names,
            reference_hierarchy = structure_factory.custom_copy(iso=uij_to_b(average_target_uijs), mask=atom_sel).select(atom_sel),
            reference_legend    = 'Target',
            reference_functions = ['mean'],
            title         = 'All Level contributions',
            reverse_legend_order = False,
            colour_indices = range(model_object.n_levels),
            legend_kw_args = dict(ncol=3, bbox_to_anchor=(0.5, 0.0), loc=9, borderaxespad=0.),
            )
        if not output_hash:
            logger.warning('no plots have been generated! ({})'.format(filenames_prefix+'*'))
        else:
            # Store in output dictionary
            file_dict['all_levels_uijs_profiles_png'] = output_hash

        # ------------------------
        # Write graph of anisotropy
        # ------------------------
        ani = 1.0 - calculate_uij_anisotropy_ratio(uij=uijs.sum(axis=0))
        filenames_prefix = self.filepath(self.all_levels_uijs_anisotropy_png_prefix)
        graph_title = ('Anisotropy of complete model'+ \
        '\nfully isotropic -> 0 (spheres) ' + \
        '\nfully anisotropic -> 1 (lines/disks)')
        output_hash = self.plot.stacked_bar(
            prefix        = filenames_prefix,
            hierarchies   = [structure_factory.custom_copy(iso=ani, mask=atom_sel).select(atom_sel)],
            legends       = ['Anisotropy  '],
            title         = graph_title,
            y_lim         = (0.0,1.05),
            y_lab         = 'Anisotropy of Uij ($1 - \\frac{E_{min}}{E_{max}}$)',
            colours       = ['grey'],
            )
        if not output_hash:
            logger.warning('no plots have been generated! ({})'.format(filenames_prefix+'*'))
        else:
            # Store in output dictionary
            file_dict['all_levels_uijs_anisotropy_png'] = output_hash

        return file_dict

    def write_tls_matrix_and_amplitude_csvs(self,
        model_object,
        ):
        """Write the various TLS uijs to the master hierarchy structure"""

        import pandas

        file_dict = collections.OrderedDict()

        # ------------------------
        # Iterate through the levels
        # ------------------------
        for i_level, tls_objects in enumerate(model_object.tls_objects):
            level_name = model_object.tls_level_names[i_level]
            # ------------------------
            # Create output tables
            # ------------------------
            # TLS model components
            mat_csv = self.filepath(self.level_tls_matrices_csv.format(i_level+1), self.csv_directory)
            file_dict.setdefault('level_tls_matrices_csv',collections.OrderedDict())[level_name] = mat_csv
            mat_tab = pandas.DataFrame(columns=["group", "mode", "label",
                                                  "T11","T22","T33","T12","T13","T23",
                                                  "L11","L22","L33","L12","L13","L23",
                                                  "S11","S12","S13","S21","S22","S23","S31","S32","S33"])
            # TLS model amplitudes
            amp_csv = self.filepath(self.level_tls_amplitudes_csv.format(i_level+1), self.csv_directory)
            file_dict.setdefault('level_tls_amplitudes_csv',collections.OrderedDict())[level_name] = amp_csv
            amp_tab = pandas.DataFrame(columns=["group", "mode", "label"]+model_object.dataset_labels, dtype=object)
            # TLS origins
            org_csv = self.filepath(self.level_tls_origins_csv.format(i_level+1), self.csv_directory)
            file_dict.setdefault('level_tls_origins_csv',collections.OrderedDict())[level_name] = org_csv
            org_tab = pandas.DataFrame(columns=["group", "label"]+model_object.dataset_labels, dtype=object)

            # ------------------------
            # Iterate through the groups in this level
            # ------------------------
            for i_group, group_object in enumerate(tls_objects):

                tls_matrices = numpy.array([mode.matrices.get() for mode in group_object.tls_parameters])
                tls_amplitudes = numpy.array([mode.amplitudes.get() for mode in group_object.tls_parameters])
                tls_origins = numpy.array(map(str, group_object.origins.round(3)), dtype=str) # extract as simple list

                assert tls_matrices.shape == (model_object.n_modes, 21)
                assert tls_amplitudes.shape == (model_object.n_modes, model_object.n_datasets)
                assert tls_origins.shape == (model_object.n_datasets,)

                # Add to model and amplitudes tables
                for i_mode in xrange(model_object.n_modes):
                    headers = numpy.array([i_group+1, i_mode+1, group_object.label], dtype=object)
                    # Add model values to last row of table
                    mat_tab.loc[len(mat_tab.index)] = numpy.concatenate([headers, tls_matrices[i_mode]])
                    # Add amplitudes to last row of table
                    amp_tab.loc[len(amp_tab.index)] = numpy.concatenate([headers, tls_amplitudes[i_mode]])

                # Add to origins table
                headers = numpy.array([i_group+1, group_object.label], dtype=object)
                org_tab.loc[len(org_tab.index)] = numpy.concatenate([headers, tls_origins])

            # Write tables
            mat_tab.to_csv(mat_csv)
            amp_tab.to_csv(amp_csv)
            org_tab.to_csv(org_csv)

        return file_dict

    def make_pymol_images_of_average_structures(self,
        level_uijs_pdbs_dict,
        model_object,
        reference_hierarchy,
        overall_atom_mask,
        ):

        logger.subheading('Generating pymol images for levels in ECHT model')

        # Output file dict
        file_dict = collections.OrderedDict()

        # Record any failed/missing files
        missing_files = []

        # Write images for each chain (average structures)
        if self.pymol_images is not None:
            # Iterate all levels
            for i_level, level_name in enumerate(model_object.all_level_names):
                f_prefix = self.filepath(self.level_uijs_pymol_by_chain_png_prefix.format(i_level+1), self.pymol_directory)
                of = auto_chain_images(
                    structure_filename = level_uijs_pdbs_dict[level_name],
                    output_prefix = f_prefix,
                    style = 'lines+ellipsoids',
                    colours = 'bfactor',
                    width=1000, height=750,
                    )
                # Check output
                if (not of):
                    logger.warning('no images have been generated: {}...'.format(f_prefix))
                else:
                    # Append missing files
                    [missing_files.append(v) for v in of.values() if not os.path.exists(v)]
                    # Store in output dictionary
                    file_dict[level_name] = of

        if missing_files:
            logger.warning('\n'.join(['image does not exist: {}'.format(v) for v in missing_files]))

        return {'level_uijs_pymol_by_chain_png' : file_dict}

    def make_pymol_script_of_average_structures(self,
        file_dict,
        ):

        from giant.pymol_utils import PythonScript

        s = PythonScript(pretty_models=False, pretty_maps=False)

        s.change_into_directory_maybe(path=os.path.abspath(self.pdb_directory))

        for f in [file_dict.get('target_uijs_pdb',None)] + file_dict.get('level_uijs_pdb',{}).values():
            if f is None: continue
            obj = os.path.basename(f)
            s.load_pdb(
                f_name = os.path.relpath(os.path.abspath(f), start=self.pdb_directory),
                obj = obj,
                )

        for f in [file_dict.get('output_uijs_pdb',None)] + file_dict.get('level_uijs_by_mode_pdb',{}).values():
            obj = os.path.basename(f)
            s.load_pdb(
                f_name = os.path.relpath(os.path.abspath(f), start=self.pdb_directory),
                obj = obj,
                )
            s.disable(obj=obj)

        s.show_as(obj='all', style='lines')
        s.show(obj='all', style='ellipsoids')
        s.custom('spectrum', expression='b', selection='all')
        s.set('grid_mode', 1)

        s.orient(obj='all')

        filename = os.path.join(self.pdb_directory, self.pymol_script_py)
        s.write_script(filename)

        return filename

    def calculate_various_average_b_statistics(self,
        level_b_values,
        level_names,
        structure_factory,
        overall_atom_mask,
        ):

        total_table_dicts = []
        chain_table_dicts = []

        level_b_values = [flex.double(b) for b in level_b_values]

        assert len(level_names) == len(level_b_values)

        total_b_mean = sum([flex.mean(b) for b in level_b_values])
        chain_b_means = {}

        for l_name, l_bvals in zip(level_names, level_b_values):

            hierarchy = structure_factory.custom_copy(
                iso = l_bvals,
                mask = overall_atom_mask,
            ).select(
                overall_atom_mask,
            )

            atom_cache = hierarchy.atom_selection_cache()

            level_b_mean = flex.mean(hierarchy.atoms().extract_b())

            if (total_b_mean > 0.0):
                level_b_mean_perc = 100.0 * level_b_mean / total_b_mean
            else:
                level_b_mean_perc = 0.0

            total_table_dicts.append({
                'Level' : l_name,
                'Chain' : 'all',
                'Average B' : level_b_mean,
                'Average B (%)' : level_b_mean_perc,
            })

            chain_ids = sorted(set([c.id for c in hierarchy.chains()]))

            for c_id in chain_ids:

                chain_sel = atom_cache.selection('chain {}'.format(c_id))
                chain_atoms = hierarchy.atoms().select(chain_sel)
                chain_level_b_mean = flex.mean(chain_atoms.extract_b())

                chain_b_means[c_id] = chain_level_b_mean + chain_b_means.get(c_id, 0)

                chain_table_dicts.append({
                    'Level' : l_name,
                    'Chain' : c_id,
                    'Average B' : chain_level_b_mean,
                })

        for d in chain_table_dicts:

            chain_id = d['Chain']
            chain_b_mean = chain_b_means[chain_id]

            if chain_b_mean > 0.0:
                chain_level_b_mean_perc = 100.0 * d['Average B'] / chain_b_mean
            else:
                chain_level_b_mean_perc = 0.0

            d['Average B (%)'] = chain_level_b_mean_perc

        # Sort in place
        chain_table_dicts.sort(key=lambda d: d['Chain'])

        import pandas
        b_factor_statistics = pandas.DataFrame(
            data = (total_table_dicts + chain_table_dicts),
            columns = ['Level', 'Chain', 'Average B', 'Average B (%)'],
        )

        return b_factor_statistics

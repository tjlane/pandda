import os, glob, traceback
import operator, collections, math
import numpy
from libtbx import adopt_init_args, group_args
from scitbx.array_family import flex
from pandemic.adp.utils import StructureFactory, PartitionBordersFactory

from bamboo.common.path import easy_directory
from giant.structure.pymol import auto_residue_images, auto_chain_images, selection_images
from giant.structure.uij import uij_to_b, \
        calculate_uij_anisotropy_ratio, \
        scale_uij_to_target_by_selection

EIGHT_PI_SQ = 8*math.pi*math.pi


class MultiModelTLSHeaderFactory:


    def __init__(self,
        tls_group_array,
        tls_selection_strings,
        tls_objects,
        overall_selection,
        ):
        adopt_init_args(self, locals())

    def __call__(self,
        i_levels = None,
        ):

        import mmtbx.tls.tools
        from bamboo.common import ListStream

        if i_levels is None: i_levels = range(len(self.tls_objects))

        assert isinstance(i_levels, list)
        assert set(map(type, i_levels)) == {int}
        i_levels = sorted(i_levels)

        # Unique combinations of the tls groups across levels
        unique_tls_combinations = sorted(set(map(tuple,self.tls_group_array[i_levels].T.tolist())))

        # List of tlso objects for each dataset
        group_tlsos = []
        group_strings = []

        # Iterate through the different UNIQUE grouping combinations
        for i_groups in unique_tls_combinations:

            i_level_groups = [(i_l, i_g) for (i_l, i_g) in zip(i_levels, i_groups) if (i_g != -1)]

            if not i_level_groups: continue

            # Combine the selection strings for this combination
            group_selection_strings = [self.overall_selection]*(self.overall_selection is not None) + \
                                      [self.tls_selection_strings[i_l][i_g] for i_l, i_g in i_level_groups]
            combined_selection_string = '('+') and ('.join(group_selection_strings)+')'

            # Extract the tls parameters for each level
            tls_objects = [self.tls_objects[i_l][i_g] for i_l, i_g in i_level_groups]
            # Expand out the tls objects for each level
            tls_parameters_exp = [[ma.expand() for ma in obj.tls_parameters] for obj in tls_objects]
            #assert set(map(len,tls_parameters_exp)) ==
            # Calculate the total for each dataset for each level
            tls_matrices_levels = [[reduce(operator.add, ml) for ml in zip(*mll)] for mll in tls_parameters_exp]
            # Sum matrices over levels
            tls_matrices = [reduce(operator.add, ms) for ms in zip(*tls_matrices_levels)]

            # Extract origins for these groups and check they're all the same (!)
            all_group_origins = [obj.origins for obj in tls_objects]
            for o in all_group_origins:
                assert o.as_double().all_eq(all_group_origins[0].as_double())
            tls_origins = all_group_origins[0]
            assert len(tls_matrices) == len(tls_origins)

            # Create tlso objects
            tlsos = [mmtbx.tls.tools.tlso(t=tuple(m.get('T')), l=tuple(m.get('L')), s=tuple(m.get('S')), origin=tuple(o)) for m,o in zip(tls_matrices, tls_origins)]

            group_tlsos.append(tlsos)
            group_strings.append(combined_selection_string)

        # Create header for each dataset
        dataset_headers = []
        for tlsos in zip(*group_tlsos):
            # Store output in list stream
            l = ListStream()
            # Get the remark records
            mmtbx.tls.tools.remark_3_tls(tlsos=tlsos, selection_strings=group_strings, out=l)
            # And finally append
            dataset_headers.append(str(l))

        return dataset_headers


class MultiModelStructureWriter:


    def __init__(self,
        output_directory,
        models,
        atom_mask,
        isotropic_mask,
        log = None,
        ):
        if log is None: log = Log()

        # Make sure output directory exists
        output_directory = easy_directory(output_directory)

        atom_mask = flex.bool(atom_mask)

        if isotropic_mask is None:
            isotropic_mask = lambda x: x

        structure_factories = [StructureFactory(master_h=m.hierarchy) for m in models]

        # Extract useful integers
        n_datasets = len(models)
        n_atoms = sum(atom_mask)

        adopt_init_args(self, locals())

    def __call__(self,
        uij,
        iso = None,
        headers = None,
        model_suffix = '.pdb',
        blank_copy = False,
        ):

        n_datasets = self.n_datasets
        n_atoms = self.n_atoms

        # Validate AnoUijs
        uij = numpy.array(uij)
        assert uij.shape == (n_datasets, n_atoms, 6)

        # Validate IsoBs
        if iso is not None:
            iso = numpy.array(iso)
            assert iso.shape == (n_datasets, n_atoms)

        # Validate header lines
        if headers is not None:
            if isinstance(headers[0], str):
                headers = [headers]
            assert set(map(len,headers)) == {len(self.models)}
            open_append = True
        else:
            open_append = False

        # List of output pdbs
        pdbs = collections.OrderedDict()

        # Apply to each model and write out
        for i, mdl in enumerate(self.models):

            # Create model paths
            mdl_d = easy_directory(os.path.join(self.output_directory, mdl.tag))
            mdl_f = os.path.join(mdl_d, mdl.tag+model_suffix)

            # Store in output dict
            pdbs[mdl.tag] = mdl_f

            try:
                # Create new copy of hierarchy and extract atoms
                h = self.structure_factories[i].custom_copy(
                    uij = self.isotropic_mask(uij[i]),
                    iso = (iso[i] if (iso is not None) else None),
                    mask = self.atom_mask,
                    blank_copy = blank_copy,
                    )

                # Write headers
                if headers is not None:
                    with open(mdl_f, 'w') as fh:
                        [fh.write(l[i].strip()+'\n') for l in headers]

                # Write atoms
                crystal_symmetry = (mdl.crystal_symmetry if hasattr(mdl, 'crystal_symmetry') else None)
                h.write_pdb_file(
                        mdl_f,
                        open_append = open_append,
                        crystal_symmetry = crystal_symmetry,
                        )

            except Exception as e:
                log = self.log
                log.bar()
                log('ERROR: Failed to write structure for dataset {}: ({})'.format(mdl.tag, mdl_f))
                log(traceback.format_exc())
                log.bar()
                continue

        return pdbs


class WriteEchtParameterSummary:


    def __init__(self,
        output_directory,
        verbose = False,
        log = None,
        ):
        adopt_init_args(self, locals())

    def __call__(self,
        params,
        plotting_object,
        ):
        """Write out the composition of the penalty functions"""

        from bamboo.maths.functions import Sigmoid

        sb = params.optimisation.first_cycle.sigmoid_buffer

        penalties_meta = [
                ('barrier_penalty',
                    Sigmoid(
                        y_scale=sb.barrier_height,
                        x_width=sb.barrier_width,
                        x_offset=sb.barrier_offset),
                    'Penalty functions for eigenvalues\nof $\Delta$U=U$_{model}$-U$_{target}$',
                    'Eigenvalue of 8*$\pi^{2}$*$\Delta$U',
                    'Penalty Value',
                    (   sb.barrier_offset-10.0*sb.barrier_width,
                        sb.barrier_offset+10.0*sb.barrier_width),
                    8.*math.pi*math.pi,
                    ),
                ]

        file_dict = collections.OrderedDict()

        for key, func, title, x_lab, y_lab, x_lim, x_scale in penalties_meta:
            x_vals = numpy.linspace(x_lim[0], x_lim[1], 101)
            y_vals = func(x_vals)
            f_name = os.path.join(self.output_directory, key+'.png')
            file_dict[key] = f_name
            plotting_object.lineplot(
                x_vals = x_scale*x_vals,
                y_vals = y_vals,
                title = title,
                x_label = x_lab,
                y_label = y_lab,
                filename = f_name,
                )

        return file_dict


class WriteEchtModelSummary:


    all_levels_uijs_profiles_png = 'all-levels-uij-profile-chain_{}.png'
    all_levels_uijs_anisotropy_png = 'all-levels-uij-anisotropy-chain_{}.png'

    level_uijs_profile_png = 'uij-profile-level_{}-chain_{}.png'
    level_uijs_anisotropy_png = 'uij-anisotropy-level_{}-chain_{}.png'

    level_uijs_pdb          = 'level_{}-all.pdb'
    level_uijs_rescaled_pdb = 'level_{}-all-rescaled.pdb'
    level_uijs_by_mode_pdb  = 'level_{}-mode_{}.pdb'

    level_uijs_pymol_by_chain_png           = 'level_{}-chain_{}.png'
    level_uijs_pymol_by_group_png           = 'level_{}-group_{}.png'
    level_uijs_pymol_by_group_rescaled_png  = 'level_{}-rescaled-group_{}.png'

    level_tls_matrices_csv   = 'tls_matrices_level_{:04d}.csv'
    level_tls_amplitudes_csv = 'tls_amplitudes_level_{:04d}.csv'
    level_tls_amplitudes_png = 'tls-model-amplitudes-level-{}-group-{}.png'

    def __init__(self,
        output_directory,
        pymol_images = None,
        distribution_images = False,
        warnings = None,
        verbose = False,
        log = None,
        ):
        csv_directory = easy_directory(os.path.join(output_directory, 'csvs'))
        pdb_directory = easy_directory(os.path.join(output_directory, 'structures'))
        pymol_directory = easy_directory(os.path.join(output_directory, 'pymol_images'))
        amplitude_directory = easy_directory(os.path.join(output_directory, 'group_amplitude_distributions'))
        adopt_init_args(self, locals())

    def show(self, file_dict, indent=0):
        log = self.log
        s = '  '
        for k, v in file_dict.iteritems():
            if isinstance(v, dict):
                log(s*indent + '> {}'.format(k))
                self.show(v, indent+1)
            elif isinstance(v, str):
                log(s*indent + '> {}: {}'.format(k, v))
            else:
                log(s*indent + '> {}'.format(k))
                for vv in v:
                    log(s*(indent+1)+vv)

    def filepath(self, filename, directory=None):
        if directory is None: directory = self.output_directory
        return os.path.join(directory, filename)

    def __call__(self,
        overall_atom_mask,
        reference_model,
        model_object,
        isotropic_mask,
        uij_target,
        level_group_array,
        results_object,
        plotting_object,
        ):

        self.plot = plotting_object

        # Ensure flex
        overall_atom_mask = flex.bool(overall_atom_mask)

        # Model values
        average_uijs_by_mode, \
        average_uijs_by_level, \
        residual_uijs = self.extract_average_uij_values(model_object)

        # Target values
        average_uijs_target = uij_target.mean(axis=0)

        # Structure-writer
        structure_factory = PartitionBordersFactory(master_h=reference_model.hierarchy)

        output_files = collections.OrderedDict()

        of = self.write_average_structures(
            average_mode_uijs = average_uijs_by_mode,
            average_level_uijs = average_uijs_by_level,
            residual_level_uijs = residual_uijs,
            overall_atom_mask = overall_atom_mask,
            isotropic_mask = isotropic_mask,
            structure_factory = structure_factory,
            model_object = model_object,
            )
        self.show(of)
        output_files.update(of)

        of = self.write_average_structure_plots(
            level_group_array = level_group_array,
            average_mode_uijs = average_uijs_by_mode,
            average_level_uijs = average_uijs_by_level,
            residual_level_uijs = residual_uijs,
            average_target_uijs = average_uijs_target,
            overall_atom_mask = overall_atom_mask,
            structure_factory = structure_factory,
            model_object = model_object,
            )
        self.show(of)
        output_files.update(of)

        of = self.write_tls_matrix_and_amplitude_csvs(
            model_object = model_object,
            )
        self.show(of)
        output_files.update(of)

        of = self.make_pymol_images_of_average_structures(
            level_uijs_pdbs_dict = output_files['level_uijs_pdb'],
            level_uijs_rescaled_pdbs_dict = output_files['level_uijs_rescaled_pdb'],
            model_object = model_object,
            reference_hierarchy = reference_model.hierarchy,
            overall_atom_mask = overall_atom_mask,
            )
        self.show(of)
        output_files.update(of)

        self.plot = None

        # Show warnings
        self.warnings.flush()

        return output_files

    def extract_average_uij_values(self,
        model_object,
        ):
        """Extract average values"""
        av_uijs_by_mode = model_object.uijs_tls().mean(axis=2)  # (n_tls, n_mode, *n_dataset*, n_atom, 6)
        av_uijs_by_level = av_uijs_by_mode.sum(axis=1)          # (n_tls, *n_mode*, n_atom, 6)
        residual_uijs = numpy.array(model_object.adp_values)    # (n_atoms, 6)
        return av_uijs_by_mode, av_uijs_by_level, residual_uijs

    def write_average_structures(self,
        average_mode_uijs,
        average_level_uijs,
        residual_level_uijs,
        overall_atom_mask,
        isotropic_mask,
        structure_factory,
        model_object,
        ):

        if isotropic_mask is None:
            isotropic_mask = lambda x: x

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
            # Write a normalised version for each level
            m_h_scl = scale_uij_to_target_by_selection(
                    hierarchy   = m_h,
                    selections  = model_object.tls_selection_strings[i_level],
                    target      = 0.1,
                    )
            m_f = self.filepath(self.level_uijs_rescaled_pdb.format(i_level+1), self.pdb_directory)
            m_h_scl.write_pdb_file(m_f)
            file_dict.setdefault('level_uijs_rescaled_pdb',collections.OrderedDict())[level_name] = m_f
        # Write residual level
        uij = residual_level_uijs
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
        residual_level_uijs,
        average_target_uijs,
        overall_atom_mask,
        structure_factory,
        model_object,
        ):

        log = self.log

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
            filenames_glob   = self.filepath(self.level_uijs_profile_png.format(i_level+1, '*'))
            filenames_prefix = filenames_glob.replace('-chain_*.png','')
            self.plot.stacked_bar(prefix        = filenames_prefix,
                                  hierarchies   = [structure_factory.custom_copy(iso=uij_to_b(u), mask=atom_sel).select(atom_sel) for u in uijs],
                                  legends       = ['TLS (Mode {})'.format(i+1) for i in xrange(n_modes)],
                                  title         = 'TLS contributions - Level {} ({})'.format(i_level+1, level_name),
                                  v_line_hierarchy = boundaries,
                                  colour_indices = [float(i_level)+(float(i_mode)/float(n_modes)) for i_mode in range(n_modes)])
            output_images = glob.glob(filenames_glob)
            if not output_images:
                self.warnings.append('no plots have been generated! ({})'.format(filenames_glob))
            else:
                # Identify the chain for each image
                output_hash = collections.OrderedDict(sorted([(s.split('_')[-1].replace('.png',''), s) for s in output_images]))
                # Store in output dictionary
                file_dict.setdefault('level_uijs_profiles_png',collections.OrderedDict())[level_name] = output_hash
            # ------------------------
            # Write graph of anisotropy
            # ------------------------
            ani = 1.0 - calculate_uij_anisotropy_ratio(uij=average_level_uijs[i_level])
            filenames_glob   = self.filepath(self.level_uijs_anisotropy_png.format(i_level+1, '*'))
            filenames_prefix = filenames_glob.replace('-chain_*.png','')
            graph_title = ('Anisotropy of Level {} ({})'+ \
            '\nfully isotropic -> 0 (spheres) ' + \
            '\nfully anisotropic -> 1 (lines/disks)').format(i_level+1, level_name)
            self.plot.stacked_bar(prefix        = filenames_prefix,
                                  hierarchies   = [structure_factory.custom_copy(iso=ani, mask=atom_sel).select(atom_sel)],
                                  legends       = ['Anisotropy  '],
                                  title         = graph_title,
                                  y_lim         = (0.0,1.05),
                                  y_lab         = 'Anisotropy of Uij ($1 - \\frac{E_{min}}{E_{max}}$)',
                                  colours       = ['grey'])
            output_images = glob.glob(filenames_glob)
            if not output_images:
                self.warnings.append('no plots have been generated! ({})'.format(filenames_glob))
            else:
                # Identify the chain for each image
                output_hash = collections.OrderedDict(sorted([(s.split('_')[-1].replace('.png',''), s) for s in output_images]))
                # Store in output dictionary
                file_dict.setdefault('level_uijs_anisotropy_png',collections.OrderedDict())[level_name] = output_hash

        # ------------------------------------------------------------------------------------------------------------------------- #
        #                                                                                                                           #
        # ------------------------------------------------------------------------------------------------------------------------- #

        # ------------------------
        # Write profiles for residual level
        # ------------------------
        uijs = residual_level_uijs
        filenames_glob   = self.filepath(self.level_uijs_profile_png.format(model_object.n_levels, '*'))
        filenames_prefix = filenames_glob.replace('-chain_*.png','')
        self.plot.stacked_bar(prefix         = filenames_prefix,
                              hierarchies    = [structure_factory.custom_copy(iso=uij_to_b(uijs), mask=atom_sel).select(atom_sel)],
                              legends        = ['{!s:12}'.format(model_object.adp_level_name)],
                              title          = 'Uij Profile of {} Level'.format(model_object.adp_level_name.capitalize()),
                              colour_indices = [model_object.n_levels-1])
        output_images = glob.glob(filenames_glob)
        if not output_images:
            self.warnings.append('no plots have been generated! ({})'.format(filenames_glob))
        else:
            # Identify the chain for each image
            output_hash = collections.OrderedDict(sorted([(s.split('_')[-1].replace('.png',''), s) for s in output_images]))
            # Store in output dictionary
            file_dict.setdefault('level_uijs_profiles_png',collections.OrderedDict())[model_object.adp_level_name] = output_hash

        # ------------------------
        # Write graph of anisotropy for residual level
        # ------------------------
        ani = 1.0 - calculate_uij_anisotropy_ratio(uij=residual_level_uijs)
        filenames_glob   = self.filepath(self.level_uijs_anisotropy_png.format(model_object.n_levels, '*'))
        filenames_prefix = filenames_glob.replace('-chain_*.png','')
        graph_title = ('Anisotropy of Level {} ({})'+ \
        '\nfully isotropic -> 0 (spheres) ' + \
        '\nfully anisotropic -> 1 (lines/disks)').format(model_object.n_levels, model_object.adp_level_name)
        self.plot.stacked_bar(prefix        = filenames_prefix,
                              hierarchies   = [structure_factory.custom_copy(iso=ani, mask=atom_sel).select(atom_sel)],
                              legends       = ['Anisotropy  '],
                              title         = graph_title,
                              y_lim         = (0.0,1.05),
                              y_lab         = 'Anisotropy of Uij ($1 - \\frac{E_{min}}{E_{max}}$)',
                              colours       = ['grey'])
        output_images = glob.glob(filenames_glob)
        if not output_images:
            self.warnings.append('no plots have been generated! ({})'.format(filenames_glob))
        else:
            # Identify the chain for each image
            output_hash = collections.OrderedDict(sorted([(s.split('_')[-1].replace('.png',''), s) for s in output_images]))
            # Store in output dictionary
            file_dict.setdefault('level_uijs_anisotropy_png',collections.OrderedDict())[model_object.adp_level_name] = output_hash

        # ------------------------------------------------------------------------------------------------------------------------- #
        #                                                                                                                           #
        # ------------------------------------------------------------------------------------------------------------------------- #

        # ------------------------
        # Write stacked profiles for all levels
        # ------------------------
        uijs = numpy.append(average_level_uijs, [residual_level_uijs], axis=0)
        assert uijs.shape == (model_object.n_levels, model_object.n_atoms, 6)
        filenames_glob =  self.filepath(self.all_levels_uijs_profiles_png.format('*'))
        filenames_prefix = filenames_glob.replace('-chain_*.png','')
        self.plot.stacked_bar(prefix        = filenames_prefix,
                              hierarchies   = [structure_factory.custom_copy(iso=uij_to_b(u), mask=atom_sel).select(atom_sel) for u in uijs],
                              legends       = model_object.all_level_names,
                              reference_hierarchy = structure_factory.custom_copy(iso=uij_to_b(average_target_uijs), mask=atom_sel).select(atom_sel),
                              reference_legend    = 'Target',
                              title         = 'All Level contributions',
                              reverse_legend_order = True,
                              colour_indices = range(model_object.n_levels),
                              )
        output_images = glob.glob(filenames_glob)
        if not output_images:
            self.warnings.append('no plots have been generated! ({})'.format(filenames_glob))
        else:
            # Identify the chain for each image
            output_hash = collections.OrderedDict(sorted([(s.split('_')[-1].replace('.png',''), s) for s in output_images]))
            # Store in output dictionary
            file_dict['all_levels_uijs_profiles_png'] = output_hash

        # ------------------------
        # Write graph of anisotropy
        # ------------------------
        ani = 1.0 - calculate_uij_anisotropy_ratio(uij=uijs.sum(axis=0))
        filenames_glob   = self.filepath(self.all_levels_uijs_anisotropy_png.format('*'))
        filenames_prefix = filenames_glob.replace('-chain_*.png','')
        graph_title = ('Anisotropy of complete model'+ \
        '\nfully isotropic -> 0 (spheres) ' + \
        '\nfully anisotropic -> 1 (lines/disks)')
        self.plot.stacked_bar(prefix        = filenames_prefix,
                              hierarchies   = [structure_factory.custom_copy(iso=ani, mask=atom_sel).select(atom_sel)],
                              legends       = ['Anisotropy  '],
                              title         = graph_title,
                              y_lim         = (0.0,1.05),
                              y_lab         = 'Anisotropy of Uij ($1 - \\frac{E_{min}}{E_{max}}$)',
                              colours       = ['grey'])
        output_images = glob.glob(filenames_glob)
        if not output_images:
            self.warnings.append('no plots have been generated! ({})'.format(filenames_glob))
        else:
            # Identify the chain for each image
            output_hash = collections.OrderedDict(sorted([(s.split('_')[-1].replace('.png',''), s) for s in output_images]))
            # Store in output dictionary
            file_dict['all_levels_uijs_anisotropy_png'] = output_hash

        return file_dict

    def write_tls_matrix_and_amplitude_csvs(self,
        model_object,
        ):
        """Write the various TLS uijs to the master hierarchy structure"""

        import pandas

        file_dict = collections.OrderedDict()

        amp_tab_all = []
        # ------------------------
        # Iterate through the levels
        # ------------------------
        for i_level, tls_objects in enumerate(model_object.tls_objects):
            level_name = model_object.tls_level_names[i_level]
            # ------------------------
            # Create output tables
            # ------------------------
            # Table for TLS model components
            mat_csv = self.filepath(self.level_tls_matrices_csv.format(i_level+1), self.csv_directory)
            file_dict.setdefault('level_tls_matrices_csv',collections.OrderedDict())[level_name] = mat_csv
            mat_tab = pandas.DataFrame(columns=["group", "mode", "label",
                                                  "T11","T22","T33","T12","T13","T23",
                                                  "L11","L22","L33","L12","L13","L23",
                                                  "S11","S12","S13","S21","S22","S23","S31","S32","S33"])
            # Create amplitude table
            amp_csv = self.filepath(self.level_tls_amplitudes_csv.format(i_level+1), self.csv_directory)
            file_dict.setdefault('level_tls_amplitudes_csv',collections.OrderedDict())[level_name] = amp_csv
            amp_tab = pandas.DataFrame(columns=["group", "mode", "label"]+model_object.dataset_labels, dtype=object)
            amp_tab_all.append(amp_tab)

            # ------------------------
            # Iterate through the groups in this level
            # ------------------------
            for i_group, group_object in enumerate(tls_objects):

                tls_matrices = numpy.array([mode.matrices.get() for mode in group_object.tls_parameters])
                tls_amplitudes = numpy.array([mode.amplitudes.get() for mode in group_object.tls_parameters])

                assert tls_matrices.shape == (model_object.n_modes, 21)
                assert tls_amplitudes.shape == (model_object.n_modes, model_object.n_datasets)

                # Add to model and amplitudes tables
                for i_mode in xrange(model_object.n_modes):
                    headers = numpy.array([i_group+1, i_mode+1, group_object.label], dtype=object)
                    # Add model values to last row of table
                    mat_tab.loc[len(mat_tab.index)] = numpy.concatenate([headers, tls_matrices[i_mode]])
                    # Add amplitudes to last row of table
                    amp_tab.loc[len(amp_tab.index)] = numpy.concatenate([headers, tls_amplitudes[i_mode]])

                # Write histograms of amplitudes -- only for non-zero models
                if (tls_matrices.sum() > 0.0) and self.distribution_images:
                    hist_png = self.filepath(self.level_tls_amplitudes_png.format(i_level+1, i_group+1), self.amplitude_directory)
                    titles = ['Mode {}:'.format(i_m+1) for i_m in xrange(model_object.n_modes)]
                    x_vals = [tls_amplitudes[i_m,:]    for i_m in xrange(model_object.n_modes)]
                    self.plot.multi_histogram(filename  = hist_png,
                                              x_vals    = x_vals,
                                              titles    = titles,
                                              x_labs    = ['']*model_object.n_modes,
                                              rotate_x_labels = True,
                                              shape     = (tls_amplitudes.shape[0], 1),
                                              n_bins    = 30, x_lim=[0, None])
                    file_dict.setdefault('level_tls_amplitudes_png',collections.OrderedDict())[(level_name,group_object.label)] = amp_csv
            # Write tables
            mat_tab.to_csv(mat_csv)
            amp_tab.to_csv(amp_csv)

        self.show_group_amplitude_statistics(amp_tab_all)

        return file_dict

    def show_group_amplitude_statistics(self,
        group_amplitudes,
        ):

        log = self.log

        log.subheading('Summary of group amplitudes for levels in ECHT model')
        for i_level, ga in enumerate(group_amplitudes):
            log.bar()
            log('Amplitude statistics by dataset: Level {}'.format(i_level+1))
            log.bar(False, True)
            # Plot histogram of amplitudes for each dataset for this level
            ga_new = ga.copy().set_index(['group','mode','label'])
            ga_new.index = range(len(ga_new.index))
            ga_new = ga_new.transpose()
            for label, values in ga_new.iterrows():
                log.bar()
                # Plot histogram of amplitudes across all groups
                if len(values) > 5:
                    try:
                        from ascii_graph import Pyasciigraph
                        g=Pyasciigraph(float_format='{0:d}')
                        counts, bounds = numpy.histogram(a=values, bins=10)
                        graph_data = [('{:.3f}-{:.3f}'.format(bounds[i],bounds[i+1]), v) for i,v in enumerate(counts)]
                        for l in g.graph(label='\n> Histogram of amplitudes for dataset {}\n'.format(label), data=graph_data):
                            if l.startswith('#######'): continue
                            log(l.replace(u"\u2588", '=').replace('= ','> '))
                    except ImportError:
                        pass
                    except:
                        raise
                # Write summary statistics
                log('\n> Amplitude statistics for level {}, dataset {}\n'.format(i_level+1, label))
                mn = values.mean()
                md = values.median()
                sd = values.std()
                log('> Mean:   {:.3f} A^2 (B-factor {:.3f} A^2)'.format(mn, EIGHT_PI_SQ*mn))
                log('> Median: {:.3f} A^2 (B-factor {:.3f} A^2)'.format(md, EIGHT_PI_SQ*md))
                log('> Std:    {:.3f} A^2 (B-factor {:.3f} A^2)'.format(sd, EIGHT_PI_SQ*sd))
            log.bar(True, True)

    def make_pymol_images_of_average_structures(self,
        level_uijs_pdbs_dict,
        level_uijs_rescaled_pdbs_dict,
        model_object,
        reference_hierarchy,
        overall_atom_mask,
        ):

        self.log.subheading('Generating pymol images for levels in ECHT model')

        file_dict = collections.OrderedDict()

        # Write images for each chain (average structures)
        if self.pymol_images is not None:
            # Iterate all levels
            for i_level, level_name in enumerate(model_object.all_level_names):
                m_f = level_uijs_pdbs_dict[level_name]
                filenames_glob   = self.filepath(self.level_uijs_pymol_by_chain_png.format(i_level+1, '*'), self.pymol_directory)
                filenames_prefix = filenames_glob.replace('-chain_*.png','')
                auto_chain_images(structure_filename = m_f,
                                  output_prefix = filenames_prefix,
                                  style = 'lines+ellipsoids',
                                  colours = 'bfactor',
                                  width=1000, height=750)
                output_images = glob.glob(filenames_glob)
                if not output_images:
                    self.warnings.append('no plots have been generated! ({})'.format(filenames_glob))
                else:
                    # Identify the chain for each image
                    output_hash = collections.OrderedDict(sorted([(s.split('_')[-1].replace('.png',''), s) for s in output_images]))
                    # Store in output dictionary
                    file_dict.setdefault('level_uijs_pymol_by_chain_png',collections.OrderedDict())[level_name] = output_hash

        # Write images of each group (average and rescaled structures)
        if self.pymol_images == 'all':
            from giant.structure.formatting import PymolSelection
            # Extract atoms covered by global mask
            r_h_filt = reference_hierarchy.select(overall_atom_mask)
            r_h_cache = r_h_filt.atom_selection_cache()
            r_h_atoms = r_h_filt.atoms()
            # Iterate only tls levels
            for i_level, level_name in enumerate(model_object.tls_level_names):

                # Extract selections for groups
                m_f = level_uijs_pdbs_dict[level_name]
                selections_atoms = [r_h_atoms.select(r_h_cache.selection(s)) for s in model_object.tls_selection_strings[i_level]]
                selections_pymol = [PymolSelection.join_or(map(PymolSelection.format,atoms)) for atoms in selections_atoms]
                n_groups = len(selections_pymol)

                # Generate unscaled images
                filenames_glob =  self.filepath(self.level_uijs_pymol_by_group_png.format(i_level+1, '*'), self.pymol_directory)
                filenames_prefix = filenames_glob.replace('*.png','')
                selection_images(structure_filename = m_f,
                                 output_prefix = filenames_prefix,
                                 selections = selections_pymol,
                                 style = 'lines+ellipsoids',
                                 width = 250 if n_groups>50 else 1000,
                                 height= 250 if n_groups>50 else 750)
                output_images = glob.glob(filenames_glob)
                if not output_images:
                    self.warnings.append('no plots have been generated! ({})'.format(filenames_glob))
                else:
                    # Identify the chain for each image
                    output_hash = collections.OrderedDict(sorted([(int(s.split('_')[-1].replace('.png','')), s) for s in output_images]))
                    # Store in output dictionary
                    file_dict.setdefault('level_uijs_pymol_by_group_png',collections.OrderedDict())[level_name] = output_hash

                # Generate rescaled images
                m_f = level_uijs_rescaled_pdbs_dict[level_name]
                filenames_glob =  self.filepath(self.level_uijs_pymol_by_group_rescaled_png.format(i_level+1, '*'), self.pymol_directory)
                filenames_prefix = filenames_glob.replace('*.png','')
                selection_images(structure_filename = m_f,
                                 output_prefix = filenames_prefix,
                                 selections = selections_pymol,
                                 style = 'lines+ellipsoids',
                                 width = 250 if n_groups>50 else 1000,
                                 height= 250 if n_groups>50 else 750)
                output_images = glob.glob(filenames_glob)
                if not output_images:
                    self.warnings.append('no plots have been generated! ({})'.format(filenames_glob))
                else:
                    # Identify the chain for each image
                    output_hash = collections.OrderedDict(sorted([(int(s.split('_')[-1].replace('.png','')), s) for s in output_images]))
                    # Store in output dictionary
                    file_dict.setdefault('level_uijs_pymol_by_group_rescaled_png',collections.OrderedDict())[level_name] = output_hash

        return file_dict


class WriteEchtStructures:


    def __init__(self,
        output_directory,
        verbose = False,
        log = None
        ):
        adopt_init_args(self, locals())

    def show(self, label, pdbs, max=5):
        log = self.log
        log(label)
        for l, p in list(pdbs.iteritems())[:max]+[('', '...')]*(len(pdbs)>max):
            log('> {}: {}'.format(l, p))

    def __call__(self,
        level_group_array,
        model_object,
        isotropic_mask,
        models,
        overall_selection = None,
        overall_selection_bool = None,
        ):
        """Extract output and generate summaries"""

        log = self.log

        make_tls_headers = MultiModelTLSHeaderFactory(
            tls_group_array = level_group_array,
            tls_selection_strings = model_object.tls_selection_strings,
            tls_objects = model_object.tls_objects,
            overall_selection = overall_selection,
            )

        output_structures = MultiModelStructureWriter(
            output_directory = self.output_directory,
            models = models,
            atom_mask = overall_selection_bool,
            isotropic_mask = isotropic_mask,
            log = log,
            )

        #------------------------------------------------------------------------------#
        #---#                            Extract uijs                              #---#
        #------------------------------------------------------------------------------#

        uij_lvl = model_object.uijs()
        uij_all = uij_lvl.sum(axis=0)
        uij_tls = uij_lvl[:-1].sum(axis=0)

        #------------------------------------------------------------------------------#
        #---#                     Write output structures                          #---#
        #------------------------------------------------------------------------------#

        output_dict = {}

        # All TLS levels
        tls_headers = make_tls_headers(i_levels=None)

        # Fully parameterised structures
        log.subheading('Writing full uij models')
        # Full models (with original values for non-fitted atoms)
        pdbs = output_structures(
            uij = uij_all,
            iso = map(uij_to_b, uij_all),
            headers = tls_headers,
            model_suffix = '.all.pdb',
            blank_copy = False)
        self.show('All levels integrated into original structures', pdbs)
        output_dict['complete_structures'] = pdbs

        # Full models (zero values for non-fitted atoms)
        pdbs = output_structures(
            uij = uij_all,
            iso = map(uij_to_b, uij_all),
            headers = tls_headers,
            model_suffix = '.all-levels.pdb',
            blank_copy = True)
        self.show('All levels', pdbs)
        output_dict['all_components'] = pdbs

        # TLS-parameterised structures (total TLS)
        log.subheading('Writing all TLS contributions')
        pdbs = output_structures(
            uij = uij_tls,
            iso = map(uij_to_b, uij_tls),
            headers = tls_headers,
            model_suffix = '.all-tls-levels.pdb',
            blank_copy = True)
        self.show('All tls levels', pdbs)
        output_dict['all_tls_levels'] = pdbs

        # Level by level TLS-parameterised structures (single level contribution)
        log.subheading('Writing individual TLS levels')
        for i_level in xrange(model_object.n_tls_levels):
            uij_this = uij_lvl[i_level]
            pdbs = output_structures(
                uij = uij_this,
                iso = map(uij_to_b, uij_this),
                headers = make_tls_headers(i_levels=[i_level]),
                model_suffix = '.tls-level-{:04}.pdb'.format(i_level+1),
                blank_copy = True)
            self.show('Level {}'.format(i_level+1), pdbs)
            output_dict.setdefault('tls_levels', collections.OrderedDict())[i_level] = pdbs

        # Level by level TLS-parameterised structures (cumulative level contributions)
        log.subheading('Writing different combinations of TLS levels')
        for i_level, j_level in flex.nested_loop((model_object.n_tls_levels, model_object.n_tls_levels)):
            if i_level >= j_level: continue
            cuml_uij = uij_lvl[i_level:j_level+1].sum(axis=0)
            pdbs = output_structures(
                uij = cuml_uij,
                iso = map(uij_to_b, cuml_uij),
                headers = make_tls_headers(i_levels=range(i_level, j_level+1)),
                model_suffix = '.tls-level-{:04}-to-{:04}.pdb'.format(i_level+1,j_level+1),
                blank_copy = True)
            self.show('Levels {}-{}'.format(i_level+1, j_level+1), pdbs)
            output_dict.setdefault('tls_levels', collections.OrderedDict())[(i_level,j_level)] = pdbs

        return output_dict


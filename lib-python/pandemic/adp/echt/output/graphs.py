import giant.logs as lg
logger = lg.getLogger(__name__)

import collections
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import figaspect

from scipy.stats import kde


class WriteTotalBFactorProfile(object):

    output_key = 'all_levels_uijs_profiles_png'
    output_path_prefix = 'profile-all-levels'

    def __init__(self,
        output_directory,
        ):

        self.output_directory = output_directory

    def __call__(self,
        target_values,
        model_values,
        structure_factory,
        atom_selection,
        plotting_object,
        *args, **kwargs
        ):

        b_values = model_values.all_levels_b

        filenames_prefix = str(
            self.output_directory / self.output_path_prefix
            )

        of = plotting_object.stacked_bar(
            prefix = filenames_prefix,
            hierarchies = [
                structure_factory.custom_copy(
                    iso = b,
                    mask = atom_selection,
                    ).select(
                    atom_selection
                    )
                for b in b_values
                ],
            legends = model_values.all_level_names,
            reference_hierarchy = structure_factory.custom_copy(
                iso = target_values.total_b,
                mask = atom_selection,
                ).select(
                atom_selection
                ),
            reference_legend = 'Target',
            reference_functions = ['mean'],
            title = 'All Level contributions',
            reverse_legend_order = False,
            colour_indices = list(range(len(model_values.all_level_names))),
            legend_kw_args = dict(ncol=3, bbox_to_anchor=(0.5, 0.0), loc=9, borderaxespad=0.),
            )

        if not of:
            logger.warning(
                'no plots have been generated! ({})'.format(filenames_prefix+'*')
                )
            return {}

        return {self.output_key : of}


class WriteTLSLevelBFactorProfiles(object):

    output_key = 'level_uijs_profiles_png'
    output_path_template_prefix = 'profile-level_{level_num}'

    def __init__(self,
        output_directory,
        ):

        self.output_directory = output_directory

        assert '{level_num}' in self.output_path_template_prefix

    def __call__(self,
        level_group_array,
        model_values,
        structure_factory,
        atom_selection,
        plotting_object,
        *args, **kwargs
        ):

        output_files = collections.OrderedDict()

        for i_level, level_name in enumerate(model_values.tls_level_names):

            # Get boundaries of groups in this level
            boundaries = structure_factory.partition_boundaries(
                atom_labels = level_group_array[i_level],
                mask = atom_selection,
                ).select(
                atom_selection
                )

            # Don't make boundaries if too many lines to plot
            if sum(boundaries.atoms().extract_b()) > 0.5*len(list(boundaries.residue_groups())):
                boundaries = None

            b_values = model_values.tls_levels_modes_b[i_level]

            filenames_prefix = str(
                self.output_directory / self.output_path_template_prefix.format(
                    level_num = i_level+1,
                    )
                )

            of = plotting_object.stacked_bar(
                prefix = filenames_prefix,
                hierarchies = [
                    structure_factory.custom_copy(
                        iso = b,
                        mask = atom_selection,
                        ).select(
                        atom_selection
                        )
                    for b in b_values
                    ],
                legends = [
                    'TLS (Mode {mode_num})'.format(
                        mode_num = i+1,
                        )
                    for i in range(model_values.n_tls_modes)
                    ],
                title = 'TLS contributions - Level {} ({})'.format(i_level+1, level_name),
                v_line_hierarchy = boundaries,
                colour_indices = [
                    float(i_level)+(float(i_mode)/float(model_values.n_tls_modes))
                    for i_mode in range(model_values.n_tls_modes)
                    ],
                )

            if not of:
                logger.warning(
                    'no plots have been generated! ({})'.format(filenames_prefix+'*')
                    )
            else:
                output_files[level_name] = of

        return {self.output_key : output_files}


class WriteAtomicLevelBFactorProfiles(object):

    output_key = 'level_uijs_profiles_png'
    output_path_template_prefix = 'profile-level_{level_num}'

    def __init__(self,
        output_directory,
        ):

        self.output_directory = output_directory

        assert '{level_num}' in self.output_path_template_prefix

    def __call__(self,
        model_values,
        structure_factory,
        atom_selection,
        plotting_object,
        *args, **kwargs
        ):

        output_files = collections.OrderedDict()

        b_values = model_values.atomic_level_b

        filenames_prefix = str(
            self.output_directory / self.output_path_template_prefix.format(
                level_num = model_values.n_levels,
                )
            )

        of = plotting_object.stacked_bar(
            prefix = filenames_prefix,
            hierarchies = [
                structure_factory.custom_copy(
                    iso = b_values,
                    mask = atom_selection,
                    ).select(
                    atom_selection
                    )
                ],
            legends = [
                '{!s:12}'.format(model_values.adp_level_name),
                ],
            title = 'Uij Profile of {level_name} Level'.format(
                level_name = model_values.adp_level_name.title(),
                ),
            colour_indices = [model_values.n_levels-1],
            )

        if not of:
            logger.warning(
                'no plots have been generated! ({})'.format(filenames_prefix+'*')
                )
            return {}
        else:
            output_files[model_values.adp_level_name] = of

        return {self.output_key : output_files}


class WriteTotalAnisotropyProfile(object):

    output_key = 'all_levels_uijs_anisotropy_png'
    output_path_prefix = 'anisotropy-all-levels'

    def __init__(self,
        output_directory,
        ):

        self.output_directory = output_directory

    def __call__(self,
        model_values,
        structure_factory,
        atom_selection,
        plotting_object,
        *args, **kwargs
        ):

        from giant.structure.uij import calculate_uij_anisotropy_ratio

        anisotropy = 1.0 - calculate_uij_anisotropy_ratio(uij=model_values.total_uij)

        filenames_prefix = str(
            self.output_directory / self.output_path_prefix
            )

        graph_title = (
            'Anisotropy of complete model\n'
            'fully isotropic -> 0 (spheres)\n'
            'fully anisotropic -> 1 (lines/disks)'
            )

        of = plotting_object.stacked_bar(
            prefix = filenames_prefix,
            hierarchies = [
                structure_factory.custom_copy(
                    iso = anisotropy,
                    mask = atom_selection,
                    ).select(atom_selection)
                ],
            legends = ['Anisotropy'],
            title = graph_title,
            y_lim = (0.0,1.05),
            y_lab = 'Anisotropy of Uij ($1 - \\frac{E_{min}}{E_{max}}$)',
            colours = ['grey'],
            )

        if not of:
            logger.warning(
                'no plots have been generated! ({})'.format(filenames_prefix+'*')
                )
            return {}

        return {self.output_key : of}


class WriteLevelAnisotropyProfiles(object):

    output_key = 'level_uijs_anisotropy_png'
    output_path_template_prefix = 'anisotropy-level_{level_num}'

    def __init__(self,
        output_directory,
        ):

        self.output_directory = output_directory

        assert '{level_num}' in self.output_path_template_prefix

    def __call__(self,
        model_values,
        structure_factory,
        atom_selection,
        plotting_object,
        *args, **kwargs
        ):

        from giant.structure.uij import calculate_uij_anisotropy_ratio

        output_files = collections.OrderedDict()

        for i_level, level_name in enumerate(model_values.all_level_names):

            level_uijs = model_values.all_levels_uij[i_level]

            anisotropy = 1.0 - calculate_uij_anisotropy_ratio(uij=level_uijs)

            filenames_prefix = str(
                self.output_directory / self.output_path_template_prefix.format(
                    level_num = i_level+1,
                    )
                )

            graph_title = (
                'Anisotropy of Level {level_num} ({level_name})\n'
                'fully isotropic -> 0 (spheres)\n'
                'fully anisotropic -> 1 (lines/disks)'
                ).format(
                level_num = i_level+1,
                level_name = level_name,
                )

            of = plotting_object.stacked_bar(
                prefix = filenames_prefix,
                hierarchies = [
                    structure_factory.custom_copy(
                        iso = anisotropy,
                        mask = atom_selection,
                        ).select(atom_selection)
                    ],
                legends = ['Anisotropy'],
                title = graph_title,
                y_lim = (0.0,1.05),
                y_lab = 'Anisotropy of Uij ($1 - \\frac{E_{min}}{E_{max}}$)',
                colours = ['grey'],
                )

            if not of:
                logger.warning(
                    'no plots have been generated! ({})'.format(filenames_prefix+'*')
                    )
            else:
                output_files[level_name] = of

        return {self.output_key : output_files}


class WriteModelBFactorDistributions(object):

    output_key = 'b_factor_distributions'
    output_path = 'b-factor-distributions.png'

    def __init__(self,
        output_directory,
        ):

        self.output_directory = output_directory

    def __call__(self,
        model_values,
        structure_factory,
        atom_selection,
        plotting_object,
        *args, **kwargs
        ):

        level_b_values = np.array(model_values.all_levels_b)

        n_levels = len(level_b_values)
        frac_n_levels = 1.0 / float(n_levels)

        x_range_all = (0.0, max(1.0, np.max(level_b_values)))

        fig = plt.figure(figsize=figaspect(0.25))

        # axes (x, y, w, h)
        main_axis = fig.add_axes([
            0.05, 0.05,
            0.65, 0.95,
            ])

        handles = []

        for i_level, level_name in enumerate(model_values.all_level_names):

            this_axis = fig.add_axes([
                0.75, 1.0 - (float(i_level+0.7)*frac_n_levels),
                0.25, 0.7 * frac_n_levels,
                ])

            level_values = level_b_values[i_level]

            x_min = max(0.0, np.min(level_values))
            x_max = max(1.0, np.max(level_values))
            x_range = (x_min, x_max)

            colour = plotting_object.get_level_colours_arbitrary([i_level])[0]

            _ = self.make_density(
                axis = main_axis,
                values = level_values,
                colour = colour,
                label = level_name,
                x_range = x_range_all,
                )

            hdl = self.make_hist(
                axis = this_axis,
                values = level_values,
                colour = colour,
                label = level_name,
                x_range = x_range,
                )

            if hdl is not None:
                handles.append(hdl)


        filename = str(
            self.output_directory / self.output_path
            )

        main_axis.set_title('B-factor Distributions for all levels')
        main_axis.set_xlabel('B-factor ($\\AA$)')
        main_axis.set_ylabel('Number Density')

        lgd = fig.legend(
            handles = handles,
            fontsize = plotting_object.helper.labelsize(0),
            bbox_transform = fig.transFigure,
            ncol = 3, loc = 9,
            bbox_to_anchor = (0.5, -0.1),
            borderaxespad = 0.,
            )

        plotting_object.helper.write_and_close_fig(
            fig = fig,
            filename = filename,
            )

        return {self.output_key : filename}

    def make_hist(self,
        axis,
        values,
        colour,
        x_range,
        label = None,
        ):

        n_array, bin_array, handles = axis.hist(
            values,
            bins = np.arange(x_range[0], 1.0+x_range[1], 1.0),
            # range = x_range,
            density = False,
            color = colour,
            label = label,
            )

        axis.set_xlim(x_range)
        axis.set_ylim((0.0, None))

        return handles[0]

    def make_density(self,
        axis,
        values,
        colour,
        x_range,
        label = None,
        ):

        try:
            density = kde.gaussian_kde(values)
            x = np.linspace(x_range[0], x_range[1], 100)
            y = density(x)
        except Exception as e:
            return

        # do first so under line
        _ = axis.fill_between(
            x, y,
            alpha = 0.5,
            color = colour,
            )

        _ = axis.plot(
            x, y,
            color = 'k',
            label = None,
            linewidth = 3,
            )

        handles = axis.plot(
            x, y,
            color = colour,
            label = label,
            )

        axis.set_xlim(x_range)
        axis.set_ylim((0.0, None))

        return handles[0]


class WriteModelLevelBFactorDistributions(WriteModelBFactorDistributions):

    output_key = 'b_factor_distributions_level'
    output_path_template = 'b-factor-distribution-level_{level_num}.png'

    def __init__(self,
        output_directory,
        ):

        self.output_directory = output_directory

        assert '{level_num}' in self.output_path_template

    def __call__(self,
        model_values,
        structure_factory,
        atom_selection,
        plotting_object,
        *args, **kwargs
        ):

        level_b_values = np.array(model_values.all_levels_b)

        output_files = collections.OrderedDict()

        for i_level, level_name in enumerate(model_values.all_level_names):

            level_values = level_b_values[i_level]

            x_min = max(0.0, np.min(level_values))
            x_max = max(1.0, np.max(level_values))
            x_range = (x_min, x_max)

            colour = plotting_object.get_level_colours_arbitrary([i_level])[0]

            fig, (a1, a2) = plt.subplots(
                nrows = 1, ncols = 2,
                squeeze = True,
                gridspec_kw = {'width_ratios': [0.7,0.3]},
                figsize = figaspect(0.25),
                )

            a1.set_title(
                'B-factor Distributions for level {level_name}'.format(
                    level_name = level_name,
                    )
                )

            a1.set_xlabel('B-factor ($\\AA$)')
            a1.set_ylabel('Number Density')
            a2.set_xlabel('B-factor ($\\AA$)')
            a2.set_ylabel('Number')

            self.make_density(
                axis = a1,
                values = level_values,
                colour = colour,
                x_range = (0.0, x_max),
                )

            self.make_hist(
                axis = a2,
                values = level_values,
                colour = colour,
                x_range = (x_min, x_max),
                )

            filename = output_files.setdefault(
                level_name,
                str(
                    self.output_directory / self.output_path_template.format(
                        level_num = i_level+1,
                        )
                    ),
                )

            plotting_object.helper.write_and_close_fig(
                fig = fig,
                filename = filename,
                )

        return {self.output_key : output_files}


###


class WriteEchtModelGraphs(object):

    def __init__(self, output_directory):

        self.output_directory = output_directory

        self.functions = [
            WriteTotalBFactorProfile(
                output_directory = output_directory,
                ),
            WriteTLSLevelBFactorProfiles(
                output_directory = output_directory,
                ),
            WriteAtomicLevelBFactorProfiles(
                output_directory = output_directory,
                ),
            WriteTotalAnisotropyProfile(
                output_directory = output_directory,
                ),
            WriteLevelAnisotropyProfiles(
                output_directory = output_directory,
                ),
            WriteModelBFactorDistributions(
                output_directory = output_directory,
                ),
            WriteModelLevelBFactorDistributions(
                output_directory = output_directory,
                ),
            ]

    def __call__(self,
        target_values,
        level_group_array,
        model_values,
        structure_factory,
        atom_selection,
        plotting_object,
        ):

        from giant.utils import merge_dicts

        if not self.output_directory.exists():
            self.output_directory.mkdir(parents=True)

        output_files = collections.OrderedDict()

        for f in self.functions:

            try:

                of = f(
                    target_values = target_values,
                    level_group_array = level_group_array,
                    model_values = model_values,
                    structure_factory = structure_factory,
                    atom_selection = atom_selection,
                    plotting_object = plotting_object,
                    )

                merge_dicts(master_dict=output_files, merge_dict=of)

            except Exception as e:

                logger.warning('Error generating graph: '+str(e))
                continue

        return output_files


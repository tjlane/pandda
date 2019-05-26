import os, math, collections
import numpy, pandas

from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure

from bamboo.common.logs import Log
from bamboo.maths.functions import rms

from pandemic.adp.utils import show_file_dict

EIGHTPISQ = 8.0*math.pi*math.pi


class AnalyseResidualsTask:

    mean_model_fit_rmsd = 'Mean RMSD'
    median_model_fit_rmsd = 'Median RMSD'

    average_b_factor_column = "Average B-factor (masked atoms)"
    target_suffix = ' (Input)'
    fitted_suffix = ' (Fitted)'

    show_file_dict = show_file_dict

    def __init__(self,
        output_directory,
        plotting_object,
        verbose = False,
        log = None,
        ):
        if log is None: log = Log()

        adopt_init_args(self, locals())

    def run(self,
        uij_fitted,
        uij_target,
        level_names,
        dataset_labels,
        reference_hierarchy,
        ):

        if level_names is None:     level_names = ['Level {}'.format(i+1) for i in range(uij_fitted.shape[0])]
        if dataset_labels is None:  dataset_labels = ['Dataset {}'.format(i+1) for i in range(uij_fitted.shape[1])]

        results_table = pandas.DataFrame(index=dataset_labels)

        output_files = collections.OrderedDict()

        self.validate(
            uij_fitted = uij_fitted,
            uij_target = uij_target,
            level_names = level_names,
            dataset_labels = dataset_labels,
            reference_hierarchy = reference_hierarchy,
            )

        from pandemic.adp.utils import StructureFactory
        structure_factory = StructureFactory(master_h=reference_hierarchy)

        self.log.subheading('Assessing model fit and calculating summary statistics')

        self.rmsd_statistics(
            results_table = results_table,
            uij_target = uij_target,
            uij_fitted = uij_fitted.sum(axis=0),
            )

        self.b_factor_statistics(
            results_table = results_table,
            uij_target = uij_target,
            uij_fitted = uij_fitted.sum(axis=0),
            )

        of = self.fit_residual_plots(
            uij_target = uij_target,
            uij_fitted = uij_fitted.sum(axis=0),
            dataset_labels = dataset_labels,
            structure_factory = structure_factory,
            )
        output_files.update(of)

        of = self.residual_v_bfactor_plots(
            uij_target = uij_target,
            uij_fitted = uij_fitted.sum(axis=0),
            dataset_labels = dataset_labels,
            )
        output_files.update(of)

        of = self.correlation_plots(
            uij_target = uij_target,
            uij_fitted = uij_fitted,
            level_names = level_names,
            dataset_labels = dataset_labels,
            structure_factory = structure_factory,
            )
        output_files.update(of)

        self.show_file_dict(output_files)

        self.result = group_args(
            output_files = output_files,
            results_table = results_table,
            )

        return self.result

    def validate(self,
        uij_fitted,
        uij_target,
        level_names,
        dataset_labels,
        reference_hierarchy,
        ):

        assert uij_fitted.shape[0] == len(level_names)
        assert uij_fitted.shape[1] == len(dataset_labels)
        assert uij_fitted.shape[2] == reference_hierarchy.atoms_size()
        assert uij_fitted.shape[3] == 6
        assert uij_fitted.shape[1:] == uij_target.shape

    def rmsd_statistics(self,
        results_table,
        uij_target,
        uij_fitted,
        ):

        log = self.log

        # Calculate rmsd between input and fitted uijs
        from bamboo.maths.functions import rms
        uij_rmsd = rms(uij_target-uij_fitted, axis=2)

        # Extract mean/median dataset-by-dataset RMSDs
        mean_rmsds = numpy.mean(uij_rmsd, axis=1)
        median_rmsds = numpy.median(uij_rmsd, axis=1)

        results_table[self.mean_model_fit_rmsd] = mean_rmsds
        results_table[self.median_model_fit_rmsd] = median_rmsds

        log('\nMean and Median RMSDs to input Uijs:')
        for i_r, (tag, data) in enumerate(results_table.iterrows()):
            if i_r == 10:
                log('...')
                break
            log('> Dataset {:10}: {:6.3f} (mean), {:6.3f} (median)'.format(
                tag,
                data[self.mean_model_fit_rmsd],
                data[self.median_model_fit_rmsd],
                ))

    def b_factor_statistics(self,
        results_table,
        uij_target,
        uij_fitted,
        ):

        log = self.log

        # Calculate isotropic ADPs for input and fitted uijs
        from giant.structure.uij import uij_to_b
        b_target = numpy.array(map(uij_to_b, uij_target))
        b_fitted = numpy.array(map(uij_to_b, uij_fitted))

        # Calculate mean/median ADPs for each atom
        mean_b_target = numpy.mean(b_target, axis=1)
        mean_b_fitted = numpy.mean(b_fitted, axis=1)

        results_table[self.average_b_factor_column+self.target_suffix] = mean_b_target
        results_table[self.average_b_factor_column+self.fitted_suffix] = mean_b_fitted

        log('\nAverage B-factors for models (fitted atoms only)')
        for i_r, (tag, data) in enumerate(results_table.iterrows()):
            if i_r == 10:
                log('...')
                break
            log('Dataset {:10}: {:6.3f} (input) -> {:6.3f} (fitted)'.format(
                tag,
                data[self.average_b_factor_column+self.target_suffix],
                data[self.average_b_factor_column+self.fitted_suffix],
                ))

    def fit_residual_plots(self,
        uij_target,
        uij_fitted,
        dataset_labels,
        structure_factory,
        ):

        uij_diff = (uij_target-uij_fitted)
        atom_rmsds = rms(uij_diff, axis=-1)

        output_files = collections.OrderedDict()

        n_datasets, n_atoms = uij_diff.shape[0], uij_diff.shape[1]

        flierprops = dict(marker='D', markersize=2., markerfacecolor='g')

        #
        # plot by residue
        #

        atom_hierarchies = [structure_factory.custom_copy(iso=EIGHTPISQ*a) for a in atom_rmsds]
        of = self.plotting_object.multi_hierarchy_plot_by_residue(
            hierarchies = atom_hierarchies,
            plot_function = self.plotting_object.boxplot,
            plot_kw_args = {
                'title' : 'Fitting Residuals by Residue\n residual = $8\pi^2$rms($U_{model} - U_{target}$)',
                'x_label' : 'Residue',
                'y_label' : 'Residuals ($\AA^2$)',
                'legends' : dataset_labels if (len(dataset_labels) < 10) else None,
                'flierprops' : flierprops,
                },
            prefix = os.path.join(self.output_directory, 'residual_by_residue'),
            residue_values_function = None,
            y_array_values_function = self.plotting_object.array_concatenate,
            )
        output_files['residuals_by_residue'] = of

        #
        # plot by dataset
        #

        of = collections.OrderedDict()
        n_per_image = 10
        for i_min in range(0, len(dataset_labels), n_per_image):
            i_max = i_min + n_per_image
            filename = os.path.join(self.output_directory, 'residual_by_dataset_{}-{}.png'.format(i_min+1, i_max))
            self.plotting_object.boxplot(
                y_vals = list(EIGHTPISQ*atom_rmsds[i_min:i_max]),
                title = 'Fitting Residuals by Dataset\n residual = $8\pi^2$rms($U_{model} - U_{target}$)',
                x_label = 'Dataset',
                y_label = 'Residuals ($\AA^2$)',
                x_tick_labels = dataset_labels[i_min:i_max],
                flierprops = flierprops,
                rotate_x_labels = True,
                filename = filename,
                )
            of['{}-{}'.format(i_min+1,i_max)] = filename
        output_files['residuals_by_dataset'] = of

        return output_files

    def residual_v_bfactor_plots(self,
        uij_target,
        uij_fitted,
        dataset_labels,
        ):

        # Plot rmsds against B-factors
        atom_rmsds = EIGHTPISQ * rms(uij_target-uij_fitted, axis=-1)
        iso_target = EIGHTPISQ * uij_target[...,0:3].mean(axis=-1)

        filename = os.path.join(self.output_directory, 'residual_vs_bfactor.png')

        n = len(iso_target)

        alpha = 0.1 + 0.9*math.e ** (-(n-1)**2 / 100**2)

        self.plotting_object.scatter(
            x_vals_array = [a for a in iso_target],
            y_vals_array = [a for a in atom_rmsds],
            title = 'Fitting Residual vs Target B-factor\nresidual = $8\pi^2$rms($U_{model} - U_{target}$)',
            x_label = 'Target B-factor ($\AA^2$)',
            y_label = 'Residual ($\AA^2$)',
            alphas = [alpha] * n,
            s = 1.5,
            legends = dataset_labels if (n < 10) else None,
            filename = filename,
            )

        return {'residuals_vs_bfactor' : filename}

    def correlation_plots(self,
        uij_target,
        uij_fitted,
        level_names,
        dataset_labels,
        structure_factory,
        ):

        output_files = collections.OrderedDict()

        uij_diff = (uij_target - uij_fitted.sum(axis=0))

        # For extracting the cross-diagonal from the correlation matrix
        diag_idxs = (numpy.arange(uij_target.shape[1]), uij_target.shape[1]+numpy.arange(uij_target.shape[1]))

        for i_l, l_name in enumerate(level_names):

            # Extract cross-correlations between atoms
            uij_correlations = [numpy.corrcoef(uij_diff[i_d], uij_fitted[i_l, i_d])[diag_idxs[0], diag_idxs[1]] for i_d in range(uij_target.shape[0])]

            atom_hierarchies = [structure_factory.custom_copy(iso=c) for c in uij_correlations]

            output_files[l_name] = self.plotting_object.multi_hierarchy_plot_by_residue(
                hierarchies = atom_hierarchies,
                plot_function = self.plotting_object.lineplot,
                plot_kw_args = {
                    'title' : 'Residual v Target Correlations\n{} level'.format(l_name.capitalize()),
                    'x_label' : 'Residue',
                    'y_label' : 'Correlation',
                    'y_lim'   : (-1.0, 1.0),
                    'legends' : dataset_labels if (len(dataset_labels) < 10) else None,
                    'marker' : 'D',
                    'markersize' : 1.5,
                    'linewidth' : 1.0,
                    },
                prefix = os.path.join(self.output_directory, 'fitting_residual_correlations_level_{}'.format(i_l+1)),
                residue_values_function = numpy.mean,
                y_array_values_function = None,
                )

        return {'residual_correlations' : output_files}

    def as_html_summary(self):
        from pandemic.adp.analysis.residuals_html import AnalyseResidualHtmlSummary
        return AnalyseResidualHtmlSummary(self)



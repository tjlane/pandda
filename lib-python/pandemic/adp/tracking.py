import os, copy
import math, collections
import numpy, pandas

from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure

from bamboo.maths.functions import rms
from giant.structure.uij import uij_to_b
from pandemic.adp import constants


class PandemicConvergenceChecker:

    def __init__(self,
        max_rmsd_b = None,
        max_delta_b = None,
        delta_b_window_frac = 0.05,
        delta_b_window_min = 5,
        eps_b = 0.01,
        log = None,
        ):
        """
        Applies multiple checks to find if the model has converged.
        rmsd_b_max: rmsd of model and target must be less than this for convergence
        delta_b_max: change in all B-factors over the last delta_n_cycles must be less than this,
            where delta_n_cycles is less given by max(delta_cycle_min, delta_cycle_frac*n_cycle_total)
        eps_b: mean B-factor of model must be greater than this
        """

        # Store history of the class
        last_cycle_history = collections.OrderedDict()

        # Cycle number where the model became non-zero
        effective_n_start = None

        adopt_init_args(self, locals())

        self._new_last_cycle_object()

    def _new_last_cycle_object(self,
        uij = 0.0,
        non_zero = False,
        delta_b = None,
        mean_b = None,
        rmsd_b = None,
        checking_from = None,
        ):

        # Store convergence values for a specific point
        self.last_cycle = group_args(
            uij = uij,
            non_zero = non_zero,
            delta_b = delta_b,
            mean_b = mean_b,
            rmsd_b = rmsd_b,
            checking_from = checking_from,
        )

        return self.last_cycle

    def update(self,
        n_cycle,
        uij_level,
        rmsd_b,
        ):

        # Sanitise
        n_cycle = int(n_cycle)

        # average b-factor
        mean_b_lvl_atom = constants.EIGHTPISQ * uij_level[..., 0:3].mean(axis=-1) # average diagonal by atom
        mean_b = mean_b_lvl_atom.sum(axis=0).mean() # sum levels and average

        # Check if model is still zero and record cycle number if not
        non_zero_b = bool((mean_b_lvl_atom > self.eps_b).any())
        if (self.effective_n_start is None) and (non_zero_b is True):
            self.effective_n_start = n_cycle

        # Calculate delta B
        largest_delta_b = 0.0
        checking_from = None
        # First cycle where model is non-zero: compare to zero
        if (self.effective_n_start == n_cycle):
            largest_delta_b = mean_b_lvl_atom.max()
        # Compare over previous N cycles
        elif (self.effective_n_start is not None):
            # Calculate number of cycles to check the convergence over
            n_non_zero_cycles = max(1, n_cycle - self.effective_n_start)
            # Size of the window (minimum size or as fraction of cycles)
            n_check_start_delta = max(
                int(self.delta_b_window_min),
                int(numpy.ceil(self.delta_b_window_frac * n_non_zero_cycles)),
            )
            # Find the start cycle - must be at least the first cycle
            n_check_start = max(1, n_cycle - n_check_start_delta)
            checking_from = n_check_start

            # Calculate changes between selected previous cycles and current cycle
            for nn_cyc in range(n_check_start, n_cycle):
                # maximum change since last cycle
                largest_delta_b = max(
                    largest_delta_b,
                    constants.EIGHTPISQ * self._find_max_u_difference(
                        uij_1 = self.last_cycle_history[nn_cyc].uij,
                        uij_2 = uij_level,
                    ),
                )

        # Reinitialise last_cycle object for overriding
        self._new_last_cycle_object(
            uij = uij_level,
            non_zero = non_zero_b,
            delta_b = largest_delta_b,
            mean_b = mean_b,
            rmsd_b = rmsd_b,
            checking_from = checking_from,
        )

        # Store history
        self.last_cycle_history[n_cycle] = copy.copy(self.last_cycle)

    def show(self):

        lc = self.last_cycle

        s = 'Convergence Checker Summary:\n'
        s += '> Model is approximately zero: {}\n'.format(
            'no' if lc.non_zero else 'yes',
        )
        s += '> Average B-factor of model: {}\n'.format(
            lc.mean_b,
        )
        s += '> Maximum change over recent cycles (since {}): {} (B-factor)\n'.format(
            lc.checking_from,
            lc.delta_b,
        )
        s += '    (cutoff for convergence: {})\n'.format(
            self.max_delta_b,
        )
        s += '> RMSD from target and fitted Uij : {} (B-factor)\n'.format(
            lc.rmsd_b,
        )
        s += '    (cutoff for convergence: {})\n'.format(
            self.max_rmsd_b,
        )

        self.log(s)

    def is_converged(self):

        # Initalise to true -- any failed test will set to false
        converged = True

        lc = self.last_cycle

        # Check if model is still zero
        if (bool(lc.non_zero) is False):
            self.log('Model is zero -- not converged')
            converged = False

        if (self.max_rmsd_b is not None) and (lc.rmsd_b > self.max_rmsd_b):
            self.log('RMSD is above threshold -- not converged')
            converged = False

        # Check if the change in B is less than tolerance
        if (self.max_delta_b is not None) and (lc.delta_b > self.max_delta_b):
            self.log('Delta B is above threshold -- not converged')
            converged = False

        return converged

    def _find_max_u_difference(self, uij_1, uij_2):
        diff = (uij_2 - uij_1)
        rms_diffs = rms(diff, axis=-1)
        min_diff = rms_diffs.max()
        return min_diff


class PandemicTrackingObject:

    csv_name1 = 'tracking_levels.csv'
    png_base1 = 'tracking_levels_'
    csv_name2 = 'tracking_rmsds.csv'
    png_base2 = 'tracking_rmsds_'

    _convergence_checker_class = PandemicConvergenceChecker

    def __init__(self,
            output_directory,
            plotting_object,
            dataset_names,
            level_names,
            convergence_args = {},
            verbose = False,
            log = None,
            ):
        if log is None: log = Log()
        # Create table for tracking progress over cycles
        table = pandas.DataFrame(
            columns=[
                'cycle', 'step', 'level', 'level#', 'rmsd',
                'u_avg', 'b_avg', 'b_min', 'b_max',
                'u_avg (total)', 'b_avg (total)',
                ],
            )
        table_by_dataset = pandas.DataFrame(
            columns=['cycle', 'type', 'overall'] + dataset_names,
            )
        tracking_csv1 = os.path.join(output_directory, self.csv_name1)
        tracking_csv2 = os.path.join(output_directory, self.csv_name2)
        tracking_png1 = os.path.join(output_directory, self.png_base1+'1.png')
        tracking_png2 = os.path.join(output_directory, self.png_base1+'2.png')
        tracking_png3 = os.path.join(output_directory, self.png_base2+'1.png')
        n_cycle = 0

        # Class to check when the model has converged
        convergence_checker = self._convergence_checker_class(log=log, **convergence_args)

        adopt_init_args(self, locals(), exclude=('convergence_args',))

    def update(self,
        uij_target,
        uij_lvl,
        step,
        i_level=None,
        write_graphs=False,
        ):
        """Update the tracking table"""

        log = self.log
        log.subheading('Updating tracking...')

        # Extract uijs for all of the levels for all datasets
        uij_dst = uij_lvl.sum(axis=0)
        # Calculate the rms between fitted and input
        u_rmsd = rms(uij_target-uij_dst, axis=None)
        b_rmsd = constants.EIGHTPISQ*u_rmsd

        # Average over all datasets
        uij_tot = uij_dst.mean(axis=0)
        uij_lvl = uij_lvl.mean(axis=1)

        if not isinstance(i_level, list):
            i_level = [i_level]

        # Iterate through levels to be dumped into table
        for i_l in i_level:

            # Extract the Uij for the selected level(s)
            if isinstance(i_l, int):
                uij_sel = uij_lvl[i_l]
                level_number = i_l+1
                level_name = self.level_names[i_l]
            else:
                uij_sel = None
                level_number = None
                level_name = None

            # Calculate U-iso & B-iso for selected level
            if uij_sel is not None:
                # average values
                b_iso_sel = numpy.mean(uij_to_b(uij_sel))
                u_iso_sel = b_iso_sel / constants.EIGHTPISQ
                # min/max values
                b_min_sel = numpy.min(uij_to_b(uij_sel))
                b_max_sel = numpy.max(uij_to_b(uij_sel))
            else:
                # average values
                b_iso_sel = 0.0
                u_iso_sel = 0.0
                # min/max values
                b_min_sel = numpy.nan
                b_max_sel = numpy.nan

            # Calculate U-iso & B-iso for complete model
            b_iso_tot = numpy.mean(uij_to_b(uij_tot))
            u_iso_tot = b_iso_tot / constants.EIGHTPISQ

            # Create human-readable cycle number
            cycle_lab = self.n_cycle

            # Add to tracking table
            self.table.loc[len(self.table.index)] = {
                'cycle'  : cycle_lab,
                'step'   : step,
                'level#' : level_number,
                'level'  : level_name,
                'rmsd'   : round(b_rmsd,3),
                'u_avg' : round(u_iso_sel,3),
                'b_avg' : round(b_iso_sel,3),
                'b_min' : round(b_min_sel,3),
                'b_max' : round(b_max_sel,3),
                'u_avg (total)' : round(u_iso_tot,3),
                'b_avg (total)' : round(b_iso_tot,3),
                }

        log(self.table.loc[len(self.table)-len(i_level):].to_string())

        if (step == 'end'):

            # Store by-dataset RMSDs
            dataset_rmsds = [constants.EIGHTPISQ*rms(d, axis=None) for d in (uij_target-uij_dst)]
            self.table_by_dataset.loc[len(self.table_by_dataset)] = [cycle_lab, 'rmsd', b_rmsd] + list(dataset_rmsds)

            # Update convergence data
            self.convergence_checker.update(
                n_cycle = self.n_cycle,
                uij_level = uij_lvl,
                rmsd_b = b_rmsd,
            )

        # Dump to csv
        self.table.to_csv(self.tracking_csv1)
        self.table_by_dataset.to_csv(self.tracking_csv2)

        if write_graphs:
            self.write_graphs()

        return

    def write_graphs(self):

        # Make plots
        self.plotting_object.tracking_plots(
            table = self.table,
            filename = self.tracking_png1,
            )
        self.plotting_object.convergence_plots(
            table = self.table,
            filename = self.tracking_png2,
            )

        tmp_table = self.table_by_dataset[(self.table_by_dataset['type'] == 'rmsd')]
        self.plotting_object.lineplot(
            x_vals = list(tmp_table['cycle'].values),
            y_vals = list(tmp_table['overall'].values),
            title = 'Model fit across cycles',
            x_label = 'Cycle',
            y_label = 'RMSD (B-factor; $\AA^2$)',
            x_ticks = tmp_table['cycle'].values,
            legends = ['rmsd'],
            filename = self.tracking_png3,
            legend_kw_args = {'bbox_to_anchor':(1.0, -0.15), 'loc':1, 'borderaxespad':0.},
            marker = '.',
            markersize = 10,
            markeredgecolor = 'k',
            linewidth = 3,
            )

        return

    def is_converged(self):
        self.log.subheading('Checking convergence')
        self.convergence_checker.show()
        cvgd = self.convergence_checker.is_converged()
        return cvgd

    def as_html_summary(self):
        from pandemic.adp.html.tracking import TrackingHtmlSummary
        return TrackingHtmlSummary(self)


from pandemic.adp.plots import PandemicAdpPlotter
from matplotlib import pyplot


class TrackingPlotter(PandemicAdpPlotter):


    def tracking_plots(self,
        table,
        filename,
        number_to_plot=5,
        ):

        # trim the table to certain rows
        start_cycle = min(table['cycle'])
        n = max(table['cycle'])
        cycles_to_plot = range(start_cycle, n, int(n/number_to_plot)+1) + [n]
        cycles_to_plot_bool = table['cycle'].isin(cycles_to_plot)
        table = table[cycles_to_plot_bool]

        fig, axes = pyplot.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
        # Create list if only one plot
        axes = numpy.array(axes).flatten()

        # Group by cycle & step to allow stacking
        grouped = table.groupby(['cycle','step'], sort=False, as_index=False)

        n_total = len(grouped)

        grouped_reduced = grouped.max()
        grouped_reduced['x'] = range(len(grouped_reduced))

        prev_x = prev_r = prev_b = None

        for n_cycle, cycle_info in grouped_reduced.groupby('cycle', sort=False):

            x_vals = cycle_info['x'].values
            r_vals = cycle_info['rmsd'].values
            b_vals = cycle_info['b_avg (total)'].values

            # Create RMSD plot
            hdl1 = axes[0].plot(
                     x_vals, r_vals,
                     'bo-', label='rmsd',
                     lw=1, ms=max(1, min(3, 5-0.1*grouped.ngroups)),
                     )
            if prev_x is not None:
                axes[0].plot(
                         [prev_x, x_vals[0]], [prev_r, r_vals[0]],
                         'b:', label='rmsd',
                         lw=1, ms=max(1, min(3, 5-0.1*grouped.ngroups)),
                         )
            # Create an overall B-iso TOTAL line
            hdl2 = axes[1].plot(
                     x_vals, b_vals,
                     'ko-', label='total',
                     lw=1, ms=max(1, min(3, 5-0.1*grouped.ngroups)),
                     )
            if prev_x is not None:
                axes[1].plot(
                         [prev_x, x_vals[0]], [prev_b, b_vals[0]],
                         'k:', label='total',
                         lw=1, ms=max(1, min(3, 5-0.1*grouped.ngroups)),
                         )

            prev_x = x_vals[-1]
            prev_r = r_vals[-1]
            prev_b = b_vals[-1]

        # Other plot -- done as one

        x_vals = numpy.arange(0, grouped.ngroups)
        x_max = max(x_vals)                 # remove
        x_keys = [v[0] for v in grouped]    # remove
        x_labs = [v[1] for v in x_keys]     # remove

        # Create B-iso lines for each LEVEL
        colours = self.get_level_colours()
        # Bottoms of bar where stacking occurs
        y_cuml = numpy.zeros(len(x_vals))
        # handles for legend at bottom of image
        handles = []
        for lvl_no in sorted(set(table['level#'])):
            if lvl_no is None: continue
            assert isinstance(lvl_no, int)
            # Get values
            sel = (table['level#'] == lvl_no)
            sel_t = table[sel]
            lvl_name = sel_t['level'].values[0]

            # The indices of the x-axis positions
            i_x = [x_keys.index(v) for v in map(tuple,sel_t[['cycle','step']].values.tolist())]
            # Extract y_vals
            y_vals = sel_t['b_avg'].values

            # Plot
            hdl = axes[1].bar(
                         x = x_vals[i_x],
                         height = y_vals,
                         bottom = y_cuml[i_x],
                         width  = 0.8,
                         color  = colours[lvl_no-1],
                         edgecolor = 'k',
                         linewidth = 0.5,
                         align  = 'center',
                         label  = '{}: {}'.format(lvl_no, lvl_name),
                         )
            handles.append(hdl)
            # Add to cuml
            y_cuml[i_x] = y_cuml[i_x] + y_vals

        # Axis stuff
        ax = axes[0]
        ax.set_title('Hierarchical Model Fit')
        ax.set_xticks([])
        ax.xaxis.set_ticks_position('bottom')
        ax.set_ylabel('model fit\n($\AA^2$)')
        ax.set_ylim(bottom=0.0)

        # Axis stuff
        ax = axes[1]
        ax.xaxis.set_ticks_position('bottom')
        ax.set_title('B-factors of Hierarchical Model')
        ax.set_xlabel('Optimisation Stage/Cycle')
        ax.set_ylabel('Isotropic B\n($\AA^2$)')
        ax.set_xticks(x_vals)
        ax.set_xticklabels(x_labs, rotation=90, ha='center')
        ax.tick_params('x', labelsize=max(2, min(10, 14-0.15*len(grouped))))
        ax.set_xlim(left=-0.5, right=max(x_vals)+0.5)
        ax.set_ylim(bottom=0.0)

        # Add legend to first graph for both lines
        lgd0a = axes[0].legend(handles=hdl1, bbox_to_anchor=(1.02, 0.95), loc=2, borderaxespad=0.)
        lgd0b = axes[1].legend(handles=hdl2, bbox_to_anchor=(1.02, 0.95), loc=2, borderaxespad=0.)

        # Other legends
        ncol = 3
        #flip_h = []; [flip_h.extend(handles[i::ncol]) for i in range(ncol)]
        flip_h = handles # (used to reorder but now keep as are)
        lgd1 = axes[1].legend(
                handles=flip_h, ncol=ncol,
                bbox_to_anchor=(0.5, 0.0),
                bbox_transform=fig.transFigure,
                loc=9, borderaxespad=0.,
                )
        # Need to re-add the old legend!
        axes[1].add_artist(lgd0b)

        # BOTH AXES -- Add vertical lines between macro-cycles
        start_x = x_vals[[x_keys.index(v[:2]) for v in map(tuple,table[['cycle','step','level#']].values.tolist()) if v[1]=='start' and v[2]==1]]
        n_cycles = len(start_x)
        last_v = None
        delta = None
        for i, v in enumerate(start_x - 0.5):
            # Dashed lines to separate cycles
            if (v > 0) and (v < x_max):
                for ax in axes:
                    ax.axvline(x=v, linewidth=1, zorder=1, color='k', linestyle='--')
            # Text to label each cycle
            if (last_v is not None):
                delta = v - last_v
                axes[0].text(x=last_v+delta/2.0,
                             y=0.1*axes[0].get_ylim()[0] + 0.9*axes[0].get_ylim()[1],
                             s='cycle '*(n_cycles<6) +str(cycles_to_plot[i-1]), # This is plotting the previous point so need -1
                             horizontalalignment='center',
                             verticalalignment='top',
                            )
            last_v = v
        # Plot the last point (or do nothing for 1 cycle)
        if delta is not None:
            axes[0].text(x=min(v+delta/2.0, axes[0].get_xlim()[1]),
                         y=0.1*axes[0].get_ylim()[0] + 0.9*axes[0].get_ylim()[1],
                         s='cycle '*(n_cycles<6) +str(cycles_to_plot[i]), # This does not need a -1
                         horizontalalignment='center',
                         verticalalignment='top',
                        )

        fig.tight_layout()
        fig.savefig(filename,
                    bbox_extra_artists=[lgd0a,lgd0b,lgd1],
                    bbox_inches='tight',
                    dpi=200)
        pyplot.close(fig)

    def convergence_plots(self,
        table,
        filename,
        ):

        fig, axes = pyplot.subplots(nrows=1, ncols=2, sharex=True, sharey=False)
        axes = numpy.array(axes).flatten()

        # Extract rows with non-zero values
        #table = table[table['b_iso (level)']!=0.0]
        # Extract only end-of-cycle optimisation values (last step of each cycle)
        table = table[(table['step']=='end')]
        # Labels for each of the series to plot
        m_cyc = 0 if (len(table) == 0) else min(table['cycle'])
        labels = table[table['cycle']==m_cyc]['level'].values

        # Colours for each level
        colours = self.get_level_colours()

        ########################
        # FIRST AXIS
        ########################
        ax = axes[0]
        handles = []
        # Extract common list of x-values
        x_keys = sorted(set(table['cycle'].values))
        x_vals = numpy.array(x_keys)
        # Select labels to plot (up to maximum)
        max_labels = 1 + 10
        plot_every = max(1, 1+((len(x_vals)-1)//max_labels))
        x_labs = [x_vals[i] if (i%plot_every)==0 else '' for i in xrange(len(x_vals))]
        # Cumulative y-values for stacking
        y_cuml = numpy.zeros(len(x_vals))
        # Plot same values as stacked bars
        for l_name in labels:
            assert isinstance(l_name, str)
            # Extract relevant rows from table
            l_table = table[table['level']==l_name]
            l_no = l_table['level#'].values[0]
            # Extract plot vals
            i_x = [x_keys.index(v) for v in l_table['cycle'].values]
            y_vals = l_table['b_avg'].values
            # Plot stacked bar
            hdl = ax.bar(x      = x_vals[i_x],
                         height = y_vals,
                         bottom = y_cuml[i_x],
                         width  = 1.0,
                         color  = colours[l_no-1],
                         edgecolor = 'k',
                         linewidth = 0.5,
                         align  = 'center',
                         label  = '{}: {}'.format(l_no, l_name),
                         )
            handles.append(hdl)
            # Add to cuml
            y_cuml[i_x] = y_cuml[i_x] + y_vals

        # Create legend for axis
        ncol = 3
        #flip_h = []; [flip_h.extend(handles[i::ncol]) for i in range(ncol)]
        flip_h = handles # (used to reorder but now keep as are)
        lgd0 = ax.legend(
                handles=flip_h, ncol=ncol,
                bbox_to_anchor=(0.5, 0.0),
                bbox_transform=fig.transFigure,
                loc=9, borderaxespad=0.,
                )

        ax.xaxis.set_ticks_position('bottom')
        ax.set_xticks(x_vals)
        ax.set_xticklabels(x_labs)
        ax.tick_params('x', labelsize=max(6, min(10, 14-0.15*len(x_keys))))
        ax.set_xlabel('Optimisation Cycle')
        ax.set_ylabel('Average B-factor of Level')
        ax.set_ylim(bottom=0.0)

        ########################
        # SECOND AXIS
        ########################
        ax = axes[1]
        handles = []
        for l_name in labels:
            assert isinstance(l_name, str)
            # Extract relevant rows from table
            l_table = table[table['level']==l_name]
            l_no = l_table['level#'].values[0]
            # Extract plot vals
            x_vals = l_table['cycle'].values
            y_vals = l_table['b_avg'].values

            hd_ = ax.plot(x_vals, y_vals,
                    'ko-',
                    lw = 2,
                    ms = max(3, min(5, 7-0.1*len(x_vals))))
            hdl = ax.plot(x_vals, y_vals,
                    'o-',
                    lw = 1,
                    ms = max(1, min(3, 5-0.1*len(x_vals))),
                    color = colours[l_no-1],
                    label = '{}: {}'.format(l_no, l_name),
                    )
            handles.extend(hdl)

        # Axis stuff
        ax.xaxis.set_ticks_position('bottom')
        ax.tick_params('x', labelsize=max(6, min(10, 14-0.15*len(x_keys))))
        ax.set_xlabel('Optimisation Cycle')
        ax.set_ylabel('Average B-factor of Level')
        ax.set_ylim(bottom=0.0)

        # Create legend for axis
        #lgd1 = ax.legend(handles=handles, ncol=3, bbox_to_anchor=(0.00, -0.15), loc=9, borderaxespad=0.)

        t = fig.suptitle('Convergence of level b-factors',
                y = 1.00, verticalalignment='bottom')

        fig.tight_layout()
        fig.savefig(filename,
                    bbox_extra_artists=[t, lgd0],
                    bbox_inches='tight',
                    dpi=200)
        pyplot.close(fig)

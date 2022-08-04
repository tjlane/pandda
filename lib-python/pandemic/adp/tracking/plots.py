import os, collections
import numpy
from libtbx import adopt_init_args

from pandemic.adp import constants


class PandemicTrackingPlotter(object):

    _snapshots_png         = 'tracking_snapshots.png'
    _level_convergence_png = 'tracking_convergence.png'
    _rmsds_convergence_png = 'tracking_rmsds.png'
    _delta_u_png_prefix    = 'tracking_atoms'

    def __init__(self,
        parent,
        output_directory,
        plotting_object,
        ):

        output_files = collections.OrderedDict(
            snapshots = os.path.join(output_directory, self._snapshots_png),
            level_convergence = os.path.join(output_directory, self._level_convergence_png),
            rmsds_convergence = os.path.join(output_directory, self._rmsds_convergence_png),
            model_changes = collections.OrderedDict(),
            )

        adopt_init_args(self, locals())

    def __call__(self):

        self.cycle_snapshot_plot(
            table = self.parent.table,
            filename = self.output_files['snapshots'],
            )

        self.level_convergence_plot(
            table = self.parent.table,
            filename = self.output_files['level_convergence']
            )

        self.rmsds_convergence_plot(
            table = self.parent.table_by_dataset,
            filename = self.output_files['rmsds_convergence'],
            )

        uij_delta_eigenvalues = self.parent.uij_history.get_delta_eigenvalues()

        of = self.delta_u_plots(
            delta_eigenvalues = uij_delta_eigenvalues.mean(axis=1),
            structure_factory = self.parent.structure_factory,
            prefix = os.path.join(self.output_directory, self._delta_u_png_prefix),
            )
        self.output_files['model_changes'] = of

        return self.output_files

    def cycle_snapshot_plot(self,
        table,
        filename,
        number_to_plot=5,
        ):

        from matplotlib import pyplot

        # Useful functions
        helper = self.plotting_object.helper

        start_cycle = min(table['cycle'])
        n = max(table['cycle'])
        # trim the table to certain rows
        cycles_to_plot = list(range(start_cycle, n, int(n/number_to_plot)+1)) + [n]
        cycles_to_plot_bool = table['cycle'].isin(cycles_to_plot)
        table = table[cycles_to_plot_bool]

        # Group by cycle & step to allow stacking
        grouped = table.groupby(['cycle','step'], sort=False, as_index=False)
        n_total = len(grouped)
        grouped_reduced = grouped.max()
        grouped_reduced['x'] = list(range(len(grouped_reduced)))

        # Previous cycle variables for connecting lines
        prev_x = prev_r = prev_b = None

        # Setup plot and plot args - do in convoluted way to allow hiding of x-ticks
        # fig, axes = pyplot.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
        # axes = numpy.array(axes).flatten() # Creates list if only one plot
        fig = pyplot.figure()
        ax1 = fig.add_subplot(2,1,1)
        ax2 = fig.add_subplot(2,1,2, sharex=ax1)
        axes = [ax1, ax2]

        # Get ls/ms outside of loop
        lw = helper.lw(grouped.ngroups)
        ms = helper.ms(grouped.ngroups)

        # Iterate through the cycles
        for n_cycle, cycle_info in grouped_reduced.groupby('cycle', sort=False):

            x_vals = cycle_info['x'].values
            r_vals = cycle_info['rmsd'].values
            b_vals = cycle_info['b_avg (total)'].values

            # Create RMSD plot
            hdl1 = axes[0].plot(
                x_vals, r_vals, 'bo-',
                label = 'rmsd',
                linewidth = lw, markersize = ms,
                )
            if prev_x is not None:
                axes[0].plot(
                    [prev_x, x_vals[0]], [prev_r, r_vals[0]], 'b:',
                    label = 'rmsd',
                    linewidth = lw, markersize = ms,
                    )
            # Create an overall B-iso TOTAL line
            hdl2 = axes[1].plot(
                x_vals, b_vals, 'ko-',
                label = 'total',
                linewidth = lw, markersize = ms,
                )
            if prev_x is not None:
                axes[1].plot(
                    [prev_x, x_vals[0]], [prev_b, b_vals[0]], 'k:',
                    label = 'total',
                    linewidth = lw, markersize = ms,
                    )

            prev_x = x_vals[-1]
            prev_r = r_vals[-1]
            prev_b = b_vals[-1]

        # Other plot -- done as one

        x_vals = numpy.arange(0, grouped.ngroups)
        x_keys = [v[0] for v in grouped]    # remove
        x_labs = [v[1] for v in x_keys]     # remove

        # Create B-iso lines for each LEVEL
        colours = self.plotting_object.get_level_colours()
        # handles for figure legend (NOT axis legend)
        figure_legend_handles = []

        # Bottoms of bar where stacking occurs
        y_cuml = numpy.zeros(len(x_vals))
        # Iterate through the levels and plot cumulative bars
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
                x       = x_vals[i_x],
                height  = y_vals,
                bottom  = y_cuml[i_x],
                width   = 0.8,
                color   = colours[lvl_no-1],
                edgecolor = 'k',
                linewidth = 0.5,
                align   = 'center',
                label   = '{}: {}'.format(lvl_no, lvl_name),
                )
            figure_legend_handles.append(hdl)
            # Add to cuml
            y_cuml[i_x] = y_cuml[i_x] + y_vals

        #
        # Axis stuff
        #
        # 1
        ax = axes[0]
        # ax.xaxis.set_ticks_position('bottom')
        helper.set_axis_labels(
            axis = ax,
            title = 'Model Fit',
            x_label = None,
            y_label = 'Model fit\n($\\AA^2$)',
            )
        helper.hide_x_labels(axis=ax)
        #ax.set_ylim(bottom=0.0)
        ax.set_yscale('log')
        #
        # 2
        ax = axes[1]
        #ax.xaxis.set_ticks_position('bottom')
        helper.set_axis_labels(
            axis = ax,
            title = 'Model B-factors',
            x_label = 'Optimisation Stage/Cycle',
            y_label = 'Isotropic B\n($\\AA^2$)',
            )
        helper.make_x_ticks(
            axis = ax,
            x_ticks = x_vals,
            x_tick_labels = x_labs,
            n_labels = 20,
            )
        helper.rotate_x_tick_labels(axis=ax)
        ax.set_xlim(left=-0.5, right=max(x_vals)+0.5)
        ax.set_ylim(bottom=0.0)

        # Add legend to first graph for both lines
        lgd0a = axes[0].legend(handles=hdl1, bbox_to_anchor=(1.02, 0.95), loc=2, borderaxespad=0.)
        lgd0b = axes[1].legend(handles=hdl2, bbox_to_anchor=(1.02, 0.95), loc=2, borderaxespad=0.)

        # BOTH AXES -- Add vertical lines between macro-cycles
        start_x = x_vals[[x_keys.index(v[:2]) for v in map(tuple,table[['cycle','step','level#']].values.tolist()) if v[1]=='start' and v[2]==1]]
        n_cycles = len(start_x)
        x_max = max(x_vals)
        last_v = delta = None
        for i, v in enumerate(start_x - 0.5):
            # Dashed lines to separate cycles
            if (v > 0) and (v < x_max):
                for ax in axes:
                    ax.axvline(x=v, linewidth=1, zorder=1, color='k', linestyle='--')
            # Text to label each cycle
            if (last_v is not None):
                delta = v - last_v
                axes[0].text(
                    x = last_v+delta/2.0,
                    y = 0.1*axes[0].get_ylim()[0] + 0.9*axes[0].get_ylim()[1],
                    s = 'cycle '*(n_cycles<6) +str(cycles_to_plot[i-1]), # This is plotting the previous point so need -1
                    horizontalalignment = 'center',
                    verticalalignment   = 'top',
                    )
            last_v = v
        # Plot the last point (or do nothing for 1 cycle)
        if delta is not None:
            axes[0].text(
                x = min(v+delta/2.0, axes[0].get_xlim()[1]),
                y = 0.1*axes[0].get_ylim()[0] + 0.9*axes[0].get_ylim()[1],
                s = 'cycle '*(n_cycles<6) +str(cycles_to_plot[i]), # This does not need a -1
                horizontalalignment = 'center',
                verticalalignment   = 'top',
                )

        # Create legend
        lgd1 = fig.legend(
            handles=figure_legend_handles, ncol=3,
            bbox_to_anchor=(0.5, 0.0),
            bbox_transform=fig.transFigure,
            loc=9, borderaxespad=0.,
            )

        helper.write_and_close_fig(
            fig=fig,
            filename=filename,
            bbox_extra_artists=[lgd0a,lgd0b,lgd1],
            )

    def level_convergence_plot(self,
        table,
        filename,
        ):

        from matplotlib import pyplot

        # Useful functions
        helper = self.plotting_object.helper

        # fig, axes = pyplot.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
        # axes = numpy.array(axes).flatten() # Creates list if only one plot
        fig = pyplot.figure()
        ax1 = fig.add_subplot(2,1,1)
        ax2 = fig.add_subplot(2,2,3)
        ax3 = fig.add_subplot(2,2,4, sharex=ax2)
        axes = [ax1, ax2, ax3]

        # Extract only end-of-cycle optimisation values (last step of each cycle)
        table = table[(table['step']=='end')]
        # Labels for each of the series to plot
        m_cyc = 0 if (len(table) == 0) else min(table['cycle'])
        labels = table[table['cycle']==m_cyc]['level'].values
        # Extract common list of x-values
        x_keys = sorted(set(table['cycle'].values))
        x_vals = numpy.array(x_keys)

        # Colours for each level
        colours = self.plotting_object.get_level_colours()
        # List of handles for figure legend (NOT axis legend)
        figure_legend_handles = []

        ########################
        # FIRST AXIS
        ########################
        ax = axes[0]
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
            hdl = ax.bar(
                x       = x_vals[i_x],
                height  = y_vals,
                bottom  = y_cuml[i_x],
                width   = 1.0,
                color   = colours[l_no-1],
                edgecolor = 'k',
                linewidth = 0.5,
                align   = 'center',
                label   = '{}: {}'.format(l_no, l_name),
                )
            figure_legend_handles.append(hdl)
            # Add to cuml
            y_cuml[i_x] = y_cuml[i_x] + y_vals

        ########################
        # SECOND AXIS
        ########################
        ax = axes[1]
        for l_name in labels:
            assert isinstance(l_name, str)
            # Extract relevant rows from table
            l_table = table[table['level']==l_name]
            l_no = l_table['level#'].values[0]
            # Extract plot vals
            x_vals = l_table['cycle'].values
            y_vals = l_table['b_avg'].values
            nx = len(x_vals)

            hdl_ = ax.plot(
                x_vals, y_vals, 'ko-',
                linewidth = helper.lw(nx, 'chunky'),
                markersize = helper.ms(nx, 'chunky'),
                )
            hdl_ = ax.plot(
                x_vals, y_vals, 'o-',
                linewidth = helper.lw(nx, 'narrow'),
                markersize = helper.ms(nx, 'narrow'),
                color = colours[l_no-1],
                label = '{}: {}'.format(l_no, l_name),
                )

        ########################
        # THIRD AXIS
        ########################
        ax = axes[2]
        for l_name in labels:
            assert isinstance(l_name, str)
            # Extract relevant rows from table
            l_table = table[table['level']==l_name]
            l_no = l_table['level#'].values[0]
            # Extract plot vals
            x_vals = l_table['cycle'].values
            plot_vals = l_table['b_avg'].values
            y_vals = numpy.concatenate(([plot_vals[0]], plot_vals[1:]-plot_vals[:-1]))
            nx = len(x_vals)

            hdl_ = ax.plot(
                x_vals, y_vals, 'ko-',
                linewidth = helper.lw(nx, 'chunky'),
                markersize = helper.ms(nx, 'chunky'),
                )
            hdl_ = ax.plot(
                x_vals, y_vals, 'o-',
                linewidth = helper.lw(nx, 'narrow'),
                markersize = helper.ms(nx, 'narrow'),
                color = colours[l_no-1],
                label = '{}: {}'.format(l_no, l_name),
                )

        #
        # Axis stuff
        #
        # 1
        ax = axes[0]
        #ax.xaxis.set_ticks_position('bottom')
        helper.set_axis_labels(
            axis = ax,
            title = 'Total B-factors',
            x_label = 'Optimisation Cycle',
            y_label = 'B-factor ($\\AA^2$)',
            )
        helper.make_x_ticks(
            axis = ax,
            x_ticks = x_keys,
            n_labels = 20,
            )
        ax.set_ylim(bottom=0.0)
        #
        # 2
        ax = axes[1]
        helper.set_axis_labels(
            axis = ax,
            title = 'Level B-factors',
            x_label = 'Optimisation Cycle',
            y_label = 'B-factor ($\\AA^2$)',
            )
        helper.make_x_ticks(
            axis = ax,
            x_ticks = x_keys,
            n_labels = 5,
            )
        ax.set_ylim(bottom=0.0)
        #
        # 3
        ax = axes[2]
        helper.set_axis_labels(
            axis = ax,
            title ='Changes between cycles',
            x_label = 'Optimisation Cycle',
            y_label = None,
            )
        ax.set_yscale('symlog', linthresh=0.1)

        # Create legend
        lgd = fig.legend(
            handles=figure_legend_handles, ncol=3,
            bbox_to_anchor=(0.5, 0.0),
            bbox_transform=fig.transFigure,
            loc=9, borderaxespad=0.,
            )

        helper.write_and_close_fig(
            fig=fig,
            filename=filename,
            bbox_extra_artists=[lgd],
            )

    def rmsds_convergence_plot(self,
        table,
        filename,
        ):

        from matplotlib import pyplot

        # Useful functions
        helper = self.plotting_object.helper

        # Filter the given table
        tmp_table = table[(table['type'] == 'rmsd')]
        # Extract x- and y-values
        x_vals = tmp_table['cycle'].values
        y_vals = tmp_table['overall'].values
        # Calculate the differences between the values
        y_vals_delta = numpy.concatenate(([0.0], y_vals[1:]-y_vals[:-1]))
        # Convert to list?

        nx = len(x_vals)

        fig, (ax1, ax2) = pyplot.subplots(nrows=2, ncols=1)

        # t = fig.suptitle(
        #     'Model fit during optimisation',
        #     y = 1.00,
        #     verticalalignment='bottom',
        #     )

        _fig, _axis = self.plotting_object.lineplot(
            axis = ax1,
            x_vals = x_vals,
            y_vals = y_vals,
            title = 'RMSD to target values',
            x_label = 'Cycle',
            y_label = 'RMSD (B-factor; $\\AA^2$)',
            x_ticks = x_vals,
            legends = ['rmsd'],
            filename = None, # returns fig and axis
            legend_kw_args = {'bbox_to_anchor':(1.0, 1.0), 'loc':4, 'borderaxespad':0.5},
            background_line_type = 'chunky',
            )

        _fig, _axis = self.plotting_object.lineplot(
            axis = ax2,
            x_vals = x_vals,
            y_vals = y_vals_delta,
            title = 'RMSD change from previous cycle',
            x_label = 'Cycle',
            y_label = 'RMSD (B-factor; $\\AA^2$)',
            x_ticks = x_vals,
            legends = ['$\\Delta$ rmsd'],
            filename = None, # returns fig and axis
            legend_kw_args = {'bbox_to_anchor':(1.0, 1.0), 'loc':4, 'borderaxespad':0.5},
            background_line_type = 'chunky',
            )

        # Set log-scale when there are non-zero values
        if len(y_vals) and (numpy.min(y_vals) > 0.0):
            # Set one log-scale above and below
            y_min = 10**numpy.floor(numpy.log10(numpy.min(y_vals)))
            y_max = 10**numpy.ceil(numpy.log10(numpy.max(y_vals)))
            ax1.set_ylim((y_min, y_max))
            ax1.set_yscale('log')
            # Just set sym log for delta plot
            ax2.set_yscale('symlog', linthresh=0.001)

        helper.write_and_close_fig(
            fig = fig,
            filename = filename,
            )

        return

    def delta_u_plots(self,
        delta_eigenvalues,
        structure_factory,
        prefix,
        ):

        helper = self.parent.plotting_object.helper

        #delta_max_sel = numpy.argmax(numpy.abs(delta_eigenvalues), axis=-1)
        #delta_indices = tuple(numpy.indices(delta_eigenvalues.shape[:-1])) + (delta_max_sel, )
        #delta_b = constants.EIGHTPISQ * delta_eigenvalues[delta_indices]
        delta_b = constants.EIGHTPISQ * delta_eigenvalues.mean(axis=-1)

        hierarchies = [structure_factory.custom_copy(iso=db, blank_copy=False) for db in delta_b]

        of = self.plotting_object.multi_hierarchy_plot_by_residue(
            hierarchies = hierarchies,
            plot_function = self.plotting_object.lineplot,
            plot_kw_args = dict(
                title   = 'Changes over last cycle',
                x_label = 'Residue',
                y_label = 'B-factor Changes ($\\AA^2$)',
                legends = self.parent.level_names,
                marker  = 'o',
                legend_kw_args = dict(ncol=3, bbox_to_anchor=(0.5, 1.1), loc=8, borderaxespad=0.),
                ),
            prefix = prefix,
            residue_values_function = numpy.mean, #max?
            y_array_values_function = None,
            )

        return of

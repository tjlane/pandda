import matplotlib
matplotlib.interactive(0)
from matplotlib import pyplot

import numpy
from scitbx.math.distributions import normal_distribution

def mean_obs_scatter(f_name, mean_vals, obs_vals):
    fig = pyplot.figure()
    pyplot.title('UNSORTED MEAN v OBS SCATTER PLOT')
    pyplot.plot([-3, 10], [-3, 10], 'b--')
    pyplot.plot(mean_vals, obs_vals, 'go')
    pyplot.xlabel('MEAN MAP VALUES')
    pyplot.ylabel('OBSV MAP VALUES')
    # Apply tight layout to prevent overlaps
    pyplot.tight_layout()
    # Save
    pyplot.savefig(f_name)
    pyplot.close(fig)

def sorted_mean_obs_scatter(f_name, mean_vals, obs_vals):
    fig = pyplot.figure()
    pyplot.title('SORTED MEAN v OBS Q-Q PLOT')
    pyplot.plot([-3, 10], [-3, 10], 'b--')
    pyplot.plot(mean_vals, obs_vals, 'go-')
    pyplot.xlabel('SORTED MEAN MAP VALUES')
    pyplot.ylabel('SORTED OBSV MAP VALUES')
    # Apply tight layout to prevent overlaps
    pyplot.tight_layout()
    # Save
    pyplot.savefig(f_name)
    pyplot.close(fig)

def diff_mean_qqplot(f_name, map_off, map_unc, q_cut, obs_diff, quantile):
    fig = pyplot.figure()
    pyplot.title('DIFF-MEAN-MAP Q-Q PLOT')
    pyplot.plot([map_off-5*map_unc, map_off+5*map_unc], [-5, 5], 'b--')
    pyplot.plot([-1, 1], [q_cut, q_cut], 'k-.')
    pyplot.plot([-1, 1], [-q_cut, -q_cut], 'k-.')
    pyplot.plot(obs_diff, quantile, 'go-')
    pyplot.xlabel('DIFF-MEAN-MAP QUANTILES')
    pyplot.ylabel('THEORETICAL QUANTILES')
    pyplot.tight_layout()
    pyplot.savefig(f_name)
    pyplot.close(fig)

def med_mean_diff_hist(f_name, plot_vals):
    fig = pyplot.figure()
    pyplot.title('MEAN-MEDIAN DIFFERENCE HISTOGRAM')
    pyplot.hist(x=plot_vals, bins=30, normed=True)
    pyplot.xlabel('MAP VALUE')
    pyplot.ylabel('DENSITY')
    pyplot.tight_layout()
    pyplot.savefig(f_name)
    pyplot.close(fig)

def map_value_distribution(f_name, plot_vals, plot_normal):
    fig = pyplot.figure()
    pyplot.title('MAP VALUE DISTRIBUTION')
    pyplot.hist(x=plot_vals, bins=30, normed=True)
    if plot_normal:
        # Plot the distribution for N(0,1)
        nd_t = normal_distribution()
        theor_x = numpy.linspace(-5,5,101)
        theor_y = [nd_t.pdf(x) for x in theor_x]
        pyplot.plot(theor_x, theor_y, c='k', ls='--', marker='o')
        # Plot the distribution for the observed distribution
        nd_o = normal_distribution(mean=numpy.mean(plot_vals), sd=numpy.std(plot_vals))
        obs_x = numpy.linspace(-5,5,101)
        obs_y = [nd_o.pdf(x) for x in obs_x]
        pyplot.plot(obs_x, obs_y, c='g', ls='-', marker='o')
    pyplot.xlabel('MAP VALUE')
    pyplot.ylabel('DENSITY')
    pyplot.tight_layout()
    pyplot.savefig(f_name)
    pyplot.close(fig)

def qq_plot_against_normal(f_name, plot_vals):
    fig = pyplot.figure()
    pyplot.title('OBS v THEORETICAL Q-Q PLOT')
    expected_vals = normal_distribution().quantiles(len(plot_vals))
    pyplot.plot([min(expected_vals)-1, max(expected_vals)+1], [min(expected_vals)-1, max(expected_vals)+1], 'b--')
    pyplot.plot(sorted(plot_vals), expected_vals, 'go-')
    pyplot.xlabel('OBSERVED QUANTILES')
    pyplot.ylabel('THEORETICAL QUANTILES')
    pyplot.tight_layout()
    pyplot.savefig(f_name)
    pyplot.close(fig)

def write_occupancy_graph(occ_est_data, d_handler):

    # 2 subplots sharing x-axis
    fig, (axis_1_1, axis_2_1) = pyplot.subplots(2, sharex=True)
    # 1st Plot - 1st Y-Axis
    line_1_1, = axis_1_1.plot(occ_est_data.occ_values, occ_est_data.global_corr_vals, 'g--', label='GLOBAL')
    line_1_2, = axis_1_1.plot(occ_est_data.occ_values, occ_est_data.local_corr_vals, 'k--', label='EVENT')
    axis_1_1.set_xlabel('Event Occupancy')
    axis_1_1.set_ylabel('Correlation to Mean Map', color='k')
    axis_1_1.set_ylim((-1, 1))
    # 1st Plot - 2nd Y-Axis
    axis_1_2 = axis_1_1.twinx()
    line_1_3, = axis_1_2.plot(occ_est_data.occ_values, occ_est_data.corr_vals_diffs, 'b-', label='DISCREP')
    axis_1_2.set_ylabel('Correlation Differences', color='b')
    # Plot line at the maximum
    line_1_4, = axis_1_2.plot([occ_est_data.feature_occ_est,occ_est_data.feature_occ_est],[-1,1], 'k-')
    text_1_1 = axis_1_2.text(0.02+occ_est_data.feature_occ_est, 0.0, str(occ_est_data.feature_occ_est))
    # Joint legend
    axis_1_1.legend(handles=[line_1_1, line_1_2, line_1_3], loc=4)
    # Title
    axis_1_1.set_title('Estimating Occupancy by Correlation Value Differences')
    # 2nd Plot - 1st Y-Axis
    line_2_1, = axis_2_1.plot(occ_est_data.occ_values, occ_est_data.global_corr_grads, 'g--', label='GLOBAL')
    line_2_2, = axis_2_1.plot(occ_est_data.occ_values, occ_est_data.local_corr_grads, 'k--', label='EVENT')
    axis_2_1.set_xlabel('Event Occupancy')
    axis_2_1.set_ylabel('Correlation to Mean Map', color='k')
    axis_2_1.set_ylim((-1, 1))
    # 1st Plot - 2nd Y-Axis
    axis_2_2 = axis_2_1.twinx()
    line_2_3, = axis_2_2.plot(occ_est_data.occ_values, occ_est_data.corr_vals_diffs, 'b-', label='DISCREP')
    axis_2_2.set_ylabel('Correlation Differences', color='b')
    # Plot line at the maximum
    line_2_4, = axis_2_2.plot([occ_est_data.feature_occ_est,occ_est_data.feature_occ_est],[-1,1], 'k-')
    text_2_1 = axis_2_2.text(0.02+occ_est_data.feature_occ_est, 0.0, str(occ_est_data.feature_occ_est))
    # Joint legend
    axis_2_1.legend(handles=[line_2_1, line_2_2, line_2_3], loc=4)
    # Title
    axis_2_1.set_title('Estimating Occupancy by Correlation Gradient Differences')
    # Remove spacing between subplots
#        fig.subplots_adjust(hspace=0)
    pyplot.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    pyplot.setp([axis_2_1.get_xticklabels()], visible=True)
    # Apply tight layout to prevent overlaps
    #pyplot.tight_layout()
    # Save
    pyplot.savefig(d_handler.output_handler.get_file('occ_corr_png').format(event_key))
    pyplot.close(fig)

def multiple_bar_plot(f_name, plot_vals, colour_bool=None):

    if colour_bool:
        assert len(plot_vals) == len(colour_bool)
        assert map(len,plot_vals) == map(len,colour_bool)

    num_sites = len(plot_vals)
    fig, all_axes = pyplot.subplots(num_sites, sharex=True)

    for i_site, bar_vals in enumerate(plot_vals):
        s_axis = all_axes[i_site]

        num_vals = len(bar_vals)
        # Left sides of the bars
        bar_left = [x+0.6 for x in range(num_vals)]
        bar_hght = bar_vals
        # Colour of the bars
        if colour_bool: bar_colr = ['green' if b else 'red' for b in colour_bool[i_site]]
        else:           bar_colr = ['blue']*num_vals

        bar_bar = s_axis.bar(left=bar_left, height=bar_hght, width=0.8, color=bar_colr)

        s_axis.set_yticks([0, int(max(bar_hght)+0.5)])

    s_axis.set_xticks(range(1,max(map(len,plot_vals))+1))

    pyplot.savefig(f_name)
    pyplot.close(fig)



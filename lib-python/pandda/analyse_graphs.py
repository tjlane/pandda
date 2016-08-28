import os, sys

#################################
try:
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.interactive(0)
    from matplotlib import pyplot
    pyplot.style.use('ggplot')
except:
    pass
#################################

import numpy

from scitbx.array_family import flex

def filter_nans(x):
    return [v for v in x if not numpy.isnan(v)]

def map_value_distribution(f_name, plot_vals, plot_normal):
    """Plot histogram of values, with optional normal distribution"""
    from scitbx.math.distributions import normal_distribution
    fig = pyplot.figure()
    pyplot.title('Map Value Distribution')
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
    pyplot.xlabel('Map Value')
    pyplot.ylabel('Density')
    pyplot.tight_layout()
    pyplot.savefig(f_name)
    pyplot.close(fig)

def qq_plot_against_normal(f_name, plot_vals):
    """Sort and plot list of values against expected quantiles from a normal distribution"""
    from scitbx.math.distributions import normal_distribution
    fig = pyplot.figure()
    pyplot.title('Observed v Theoretical Q-Q Plot')
    expected_vals = normal_distribution().quantiles(len(plot_vals))
    pyplot.plot([min(expected_vals)-1, max(expected_vals)+1], [min(expected_vals)-1, max(expected_vals)+1], 'b--')
    pyplot.plot(sorted(plot_vals), expected_vals, 'go-')
    pyplot.xlabel('Observed Quantiles')
    pyplot.ylabel('Theoretical Quantiles')
    pyplot.tight_layout()
    pyplot.savefig(f_name)
    pyplot.close(fig)

def mean_obs_scatter(f_name, mean_vals, obs_vals):
    """Plot mean map against observed map values"""
    fig = pyplot.figure()
    pyplot.title('Unsorted Mean v Dataset Scatter Plot')
    pyplot.plot([-3, 10], [-3, 10], 'b--')
    pyplot.plot(mean_vals, obs_vals, 'go')
    pyplot.xlabel('Mean Map Values')
    pyplot.ylabel('Dataset Map Values')
    # Apply tight layout to prevent overlaps
    pyplot.tight_layout()
    # Save
    pyplot.savefig(f_name)
    pyplot.close(fig)

def sorted_mean_obs_scatter(f_name, mean_vals, obs_vals):
    """Plot sorted mean map against sorted observed map values"""
    fig = pyplot.figure()
    pyplot.title('Sorted Mean v Dataset Q-Q Scatter Plot')
    pyplot.plot([-3, 10], [-3, 10], 'b--')
    pyplot.plot(mean_vals, obs_vals, 'go')
    pyplot.xlabel('Sorted Mean Map Values')
    pyplot.ylabel('Sorted Dataset Map Values')
    # Apply tight layout to prevent overlaps
    pyplot.tight_layout()
    # Save
    pyplot.savefig(f_name)
    pyplot.close(fig)

def diff_mean_qqplot(f_name, map_off, map_unc, q_cut, obs_diff, quantile):
    """Plot diff-mean map against normal quantiles, with uncertainty lines"""
    fig = pyplot.figure()
    pyplot.title('Mean-Difference Map Q-Q Plot')
    pyplot.plot([map_off-5*map_unc, map_off+5*map_unc], [-5, 5], 'b--')
    pyplot.plot([-1, 1], [q_cut, q_cut], 'k-.')
    pyplot.plot([-1, 1], [-q_cut, -q_cut], 'k-.')
    pyplot.plot(obs_diff, quantile, 'go-')
    pyplot.xlabel('Difference from Mean Map Quantiles')
    pyplot.ylabel('Theoretical Quantiles')
    pyplot.tight_layout()
    pyplot.savefig(f_name)
    pyplot.close(fig)

def write_dataset_summary_graphs(pandda):
    """Plot dataset summary graphs of resolution, unit cell variation, etc"""

    pandda.log('----------------------------------->>>')
    pandda.log('Generating Summary Graphs')

    # Filter the dataset to non-rejected datasets
    non_rejected_dtags = [dh.tag for dh in pandda.datasets.mask(mask_name='rejected - total', invert=True)]

    d_info = pandda.tables.dataset_info.loc[non_rejected_dtags]
    n_bins = 30

    # ================================================>
    fig = pyplot.figure()
    pyplot.title('Resolution Histograms')
    pyplot.subplot(2, 1, 1)
    pyplot.hist(x=d_info['high_resolution'], bins=n_bins)
    pyplot.xlabel('High Resolution Limit (A)')
    pyplot.ylabel('Count')
    pyplot.subplot(2, 1, 2)
    pyplot.hist(x=d_info['low_resolution'], bins=n_bins)
    pyplot.xlabel('Low Resolution Limit (A)')
    pyplot.ylabel('Count')
    pyplot.tight_layout()
    pyplot.savefig(pandda.output_handler.get_file('d_resolutions'))
    pyplot.close(fig)
    # ================================================>
    try:
        fig = pyplot.figure()
        pyplot.title('R-Factor Histograms')
        pyplot.subplot(2, 1, 1)
        pyplot.hist(x=d_info['r_free'], bins=n_bins)
        pyplot.xlabel('R-Free')
        pyplot.ylabel('Count')
        pyplot.subplot(2, 1, 2)
        pyplot.hist(x=d_info['r_work'], bins=n_bins)
        pyplot.xlabel('R-Work')
        pyplot.ylabel('Count')
        pyplot.tight_layout()
        pyplot.savefig(pandda.output_handler.get_file('d_rfactors'))
        pyplot.close(fig)
    except:
        fig = pyplot.figure()
        pyplot.title('Failed to make R-Factor Histograms')
        pyplot.savefig(pandda.output_handler.get_file('d_rfactors'))
        pyplot.close(fig)
    # ================================================>
    fig = pyplot.figure()
    pyplot.title('RMSDs to Reference Structure Histogram')
    pyplot.hist(x=filter_nans(d_info['rmsd_to_reference']), bins=n_bins)
    pyplot.xlabel('RMSD (A)')
    pyplot.ylabel('Count')
    pyplot.tight_layout()
    pyplot.savefig(pandda.output_handler.get_file('d_global_rmsd_to_ref'))
    pyplot.close(fig)
    # ================================================>
    fig = pyplot.figure()
    pyplot.title('Unit Cell Axis Variation')
    pyplot.subplot(3, 1, 1)
    pyplot.hist(x=d_info['uc_a'], bins=n_bins)
    pyplot.xlabel('A (A)')
    pyplot.subplot(3, 1, 2)
    pyplot.hist(x=d_info['uc_b'], bins=n_bins)
    pyplot.xlabel('B (A)')
    pyplot.subplot(3, 1, 3)
    pyplot.hist(x=d_info['uc_c'], bins=n_bins)
    pyplot.xlabel('C (A)')
    pyplot.tight_layout()
    pyplot.savefig(pandda.output_handler.get_file('d_cell_axes'))
    pyplot.close(fig)
    # ================================================>
    fig = pyplot.figure()
    pyplot.title('Unit Cell Angle Variation')
    pyplot.subplot(3, 1, 1)
    pyplot.hist(x=d_info['uc_alpha'], bins=n_bins)
    pyplot.xlabel('Alpha')
    pyplot.subplot(3, 1, 2)
    pyplot.hist(x=d_info['uc_beta'], bins=n_bins)
    pyplot.xlabel('Beta')
    pyplot.subplot(3, 1, 3)
    pyplot.hist(x=d_info['uc_gamma'], bins=n_bins)
    pyplot.xlabel('Gamma')
    pyplot.tight_layout()
    pyplot.savefig(pandda.output_handler.get_file('d_cell_angles'))
    pyplot.close(fig)
    # ================================================>
    fig = pyplot.figure()
    pyplot.title('Unit Cell Volume Variation')
    pyplot.hist(x=d_info['uc_vol'], bins=n_bins)
    pyplot.xlabel('Volume (A^3)')
    pyplot.ylabel('Count')
    pyplot.tight_layout()
    pyplot.savefig(pandda.output_handler.get_file('d_cell_volumes'))
    pyplot.close(fig)

def write_map_analyser_graphs(pandda, map_analyser, analysis_mask_name):

    # Get the output directory to write the graphs into
    img_out_dir = os.path.join(pandda.output_handler.get_dir('analyses'), '{!s}A Maps'.format(map_analyser.meta.resolution))
    if not os.path.exists(img_out_dir): os.mkdir(img_out_dir)

    # Map Analyser Variables
    map_res = map_analyser.meta.resolution

    ########################################################

    # Statistical Map Values

    # Points to extract from the grid
    grid_mask = pandda.reference_grid().global_mask().outer_mask_indices()

    # Extract map values for plotting
    mean_map_vals = list(map_analyser.statistical_maps.mean_map.select(grid_mask))
    try:    medn_map_vals = list(map_analyser.statistical_maps.medn_map.select(grid_mask))
    except: medn_map_vals = mean_map_vals
    stds_map_vals = list(map_analyser.statistical_maps.stds_map.select(grid_mask))
    sadj_map_vals = list(map_analyser.statistical_maps.sadj_map.select(grid_mask))

    ########################################################

    # Dataset Parameters
    d_info = pandda.tables.dataset_info
    m_info = pandda.tables.dataset_map_info

    # All datasets
    high_res = [d_info['high_resolution'][mh.tag] for mh in map_analyser.dataset_maps.all()]
    low_res =  [d_info['low_resolution'][mh.tag]  for mh in map_analyser.dataset_maps.all()]
    rfree =    [d_info['r_free'][mh.tag]          for mh in map_analyser.dataset_maps.all()]
    rwork =    [d_info['r_work'][mh.tag]          for mh in map_analyser.dataset_maps.all()]

    # All datasets
    map_uncties = [m_info['map_uncertainty'][mh.tag] for mh in map_analyser.dataset_maps.all()]
    # Analysed datasets only
    z_map_mean  = [m_info['z_map_mean'][mh.tag] for mh in map_analyser.dataset_maps.mask(mask_name=analysis_mask_name)]
    z_map_std   = [m_info['z_map_std'][mh.tag]  for mh in map_analyser.dataset_maps.mask(mask_name=analysis_mask_name)]
    z_map_skew  = [m_info['z_map_skew'][mh.tag] for mh in map_analyser.dataset_maps.mask(mask_name=analysis_mask_name)]
    z_map_kurt  = [m_info['z_map_kurt'][mh.tag] for mh in map_analyser.dataset_maps.mask(mask_name=analysis_mask_name)]

    ########################################################

    n_bins = 30

    ########################################################

    pandda.log('=> Writing Statistical Map Distributions')

    ##################################
    # STATISTICAL MAP HISTOGRAMS
    ##################################
    fig = pyplot.figure()
    pyplot.title('Mean and Median Map Histograms')
    # MEAN MAP
    pyplot.subplot(2, 1, 1)
    pyplot.hist(x=mean_map_vals, bins=n_bins)
    pyplot.xlabel('Mean Map Values')
    pyplot.ylabel('Count')
    # MEDIAN MAP
    pyplot.subplot(2, 1, 2)
    pyplot.hist(x=medn_map_vals, bins=n_bins)
    pyplot.xlabel('Median Map Values')
    pyplot.ylabel('Count')
    # Apply tight layout to prevent overlaps
    pyplot.tight_layout()
    # Save both
    pyplot.savefig(os.path.join(img_out_dir, '{!s}A-mean_medn_map_hist.png'.format(map_res)))
    pyplot.close(fig)

    ##################################
    # MEAN Values v MEDIAN Values
    ##################################
    fig = pyplot.figure()
    pyplot.title('Mean v Median Map Values')
    pyplot.scatter(x=mean_map_vals, y=medn_map_vals)
    # Plot straight line between the min and max values
    min_val = min(mean_map_vals+medn_map_vals)
    max_val = max(mean_map_vals+medn_map_vals)
    pyplot.plot([min_val, max_val], [min_val, max_val], 'b--')
    # Axis labels
    pyplot.xlabel('Mean Map Value')
    pyplot.ylabel('Median Map Value')
    # Apply tight layout to prevent overlaps
    pyplot.tight_layout()
    # Save
    pyplot.savefig(os.path.join(img_out_dir, '{!s}A-mean_v_median_scatter.png'.format(map_res)))
    pyplot.close(fig)

    ##################################
    # MEAN-MEDIAN DIFFERENCE HISTOGRAM
    ##################################
    fig = pyplot.figure()
    pyplot.title('Mean-Median Difference Histogram')
    pyplot.hist(x=numpy.abs(flex.double(mean_map_vals)-flex.double(medn_map_vals)), bins=30, normed=True)
    pyplot.xlabel('Difference Map Value')
    pyplot.ylabel('Density')
    pyplot.tight_layout()
    pyplot.savefig(os.path.join(img_out_dir, '{!s}A-mean_median_diff_hist.png'.format(map_res)))
    pyplot.close(fig)

    ##################################
    # STATISTICAL MAP HISTOGRAMS
    ##################################
    fig = pyplot.figure()
    pyplot.title('Point Variation Map Histograms')
    # STANDARD DEVIATION MAPS
    pyplot.subplot(2, 1, 1)
    pyplot.hist(x=stds_map_vals, bins=n_bins)
    pyplot.xlabel('"Raw" Variation of Map Values')
    pyplot.ylabel('Count')
    # ADJUSTED STANDARD DEVIATION MAPS
    pyplot.subplot(2, 1, 2)
    pyplot.hist(x=sadj_map_vals, bins=n_bins)
    pyplot.xlabel('"Adjusted" Variation of Map Values')
    pyplot.ylabel('Count')
    # Apply tight layout to prevent overlaps
    pyplot.tight_layout()
    # Save both
    pyplot.savefig(os.path.join(img_out_dir, '{!s}A-stds_sadj_map_vals.png'.format(map_res)))
    pyplot.close(fig)

    ##################################
    # STD Values v ADJ STD Values
    ##################################
    fig = pyplot.figure()
    pyplot.title('Raw v Adjusted Map Variation')
    pyplot.scatter(x=stds_map_vals, y=sadj_map_vals)
    # Plot straight line between the min and max values
    min_val = min(stds_map_vals+sadj_map_vals)
    max_val = max(stds_map_vals+sadj_map_vals)
    pyplot.plot([min_val, max_val], [min_val, max_val], 'b--')
    # Axis labels
    pyplot.xlabel('"Raw" Variation of Map Values')
    pyplot.ylabel('"Adjusted" Variation of Map Values')
    # Apply tight layout to prevent overlaps
    pyplot.tight_layout()
    # Save
    pyplot.savefig(os.path.join(img_out_dir, '{!s}A-std_v_adj_std_scatter.png'.format(map_res)))
    pyplot.close(fig)

    ########################################################

    pandda.log('=> Writing Map Uncertainties')

    # MAP PARAMS
    fig = pyplot.figure()
    pyplot.title('Map Statistics Histrograms')
    # MAP UNCERTAINTIES
    pyplot.hist(x=map_uncties, bins=n_bins, range=(min(map_uncties)-0.1,max(map_uncties)+0.1))
    pyplot.xlabel('Dataset Map Uncertainties')
    pyplot.ylabel('Count')
    # Apply tight layout to prevent overlaps
    pyplot.tight_layout()
    # Save
    pyplot.savefig(os.path.join(img_out_dir, '{!s}A-d_map_uncertainties.png'.format(map_res)))
    pyplot.close(fig)

    ########################################################

    pandda.log('=> Scatter Plots')

    # MAP RESOLUTION V UNCERTAINTY
    fig = pyplot.figure()
    pyplot.title('High Resolution Limit v Map Uncertainty')
    pyplot.scatter(x=high_res, y=map_uncties)
    pyplot.xlabel('Dataset Resolution')
    pyplot.ylabel('Dataset Map Uncertainty')
    # Apply tight layout to prevent overlaps
    pyplot.tight_layout()
    # Save
    pyplot.savefig(os.path.join(img_out_dir, '{!s}A-resolution_v_uncertainty.png'.format(map_res)))
    pyplot.close(fig)

    # MAP RESOLUTION V UNCERTAINTY
    fig = pyplot.figure()
    pyplot.title('High Resolution Limit v R-Free')
    pyplot.scatter(x=high_res, y=rfree)
    pyplot.xlabel('Dataset Resolution')
    pyplot.ylabel('Dataset R-Free')
    # Apply tight layout to prevent overlaps
    pyplot.tight_layout()
    # Save
    pyplot.savefig(os.path.join(img_out_dir, '{!s}A-resolution_v_rfree.png'.format(map_res)))
    pyplot.close(fig)

    # RFREE V UNCERTAINTY
    fig = pyplot.figure()
    pyplot.title('R-Free v Uncertainty')
    pyplot.scatter(x=rfree, y=map_uncties)
    pyplot.xlabel('Dataset R-Free')
    pyplot.ylabel('Dataset Map Uncertainty')
    # Apply tight layout to prevent overlaps
    pyplot.tight_layout()
    # Save
    pyplot.savefig(os.path.join(img_out_dir, '{!s}A-rfree_v_uncertainty.png'.format(map_res)))
    pyplot.close(fig)

    ########################################################

    pandda.log('=> Z-Map Distribution')

    # R-FACTORS
    fig = pyplot.figure()
    pyplot.title('Z-Map Statistics Histograms')
    # RFree
    pyplot.subplot(2, 1, 1)
    pyplot.hist(x=z_map_mean, bins=n_bins, range=(min(z_map_mean)-0.1,max(z_map_mean)+0.1))
    pyplot.xlabel('Z-Map Mean')
    pyplot.ylabel('Count')
    # RWork
    pyplot.subplot(2, 1, 2)
    pyplot.hist(x=z_map_std, bins=n_bins, range=(0, max(z_map_std)+0.1))
    pyplot.xlabel('Z-Map Std')
    pyplot.ylabel('Count')
    # Apply tight layout to prevent overlaps
    pyplot.tight_layout()
    # Save both
    pyplot.savefig(os.path.join(img_out_dir, '{!s}A-z_map_statistics.png'.format(map_res)))
    pyplot.close(fig)

    # Z-MAP SKEW V UNCERTAINTY
    fig = pyplot.figure()
    pyplot.title('Z-Map Normality Plots')
    pyplot.scatter(x=z_map_skew, y=z_map_kurt)
    pyplot.xlabel('Skew')
    pyplot.ylabel('Kurtosis')
    # Apply tight layout to prevent overlaps
    pyplot.tight_layout()
    # Save
    pyplot.savefig(os.path.join(img_out_dir, '{!s}A-z_map_skew_v_kurtosis.png'.format(map_res)))
    pyplot.close(fig)

def write_occupancy_graph(x_values, global_values, local_values, diff_values, out_file):
    """Write output graph from occupancy estimation"""

    # Get the x-value of the maximum difference
    max_x = x_values[diff_values.index(max(diff_values))]

    fig, (axis_1_1, axis_2_1) = pyplot.subplots(2, sharex=True)
    # 1st Plot - 1st Y-Axis
    line_1_1, = axis_1_1.plot(x_values, global_values, 'g--', label='GLOBAL CORRELATION', linewidth=2)
    line_1_2, = axis_1_1.plot(x_values, local_values, 'k--', label='LOCAL CORRELATION', linewidth=2)
    axis_1_1.set_ylabel('CORRELATION\nTO GROUND STATE', color='k', size=16)
    axis_1_1.set_ylim((-1, 1))
    # 2nd Plot - 1st Y-Axis
    line_2_1, = axis_2_1.plot(x_values, diff_values, 'b-', label='DIFFERENCE', linewidth=2)
    axis_2_1.set_ylabel('CORRELATION\nDIFFERENCE', color='k', size=16)
    axis_2_1.set_xlabel('1-BDC', color='k', size=16)
    axis_2_1.set_ylim((min(diff_values)-0.2, max(diff_values)+0.2))
    # Plot line at the maximum
    line_2_2, = axis_2_1.plot([max_x,max_x],[-1,1], 'k-', linewidth=2)
    text_2_1 = axis_2_1.text(0.02+max_x, 0.0, 'BDC='+str(1-max_x), size=14)
    # Joint legend
    axis_1_1.legend(handles=[line_1_1, line_1_2, line_2_1], loc=4, fontsize=16)
#    # Title
#    axis_2_1.set_title('Estimating Occupancy by Correlation Gradient Differences')
    # Remove spacing between subplots
    pyplot.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    pyplot.setp([axis_2_1.get_xticklabels()], visible=True)
    # Apply tight layout to prevent overlaps
    pyplot.tight_layout()
    # Save
#    pyplot.savefig(d_handler.output_handler.get_file('occ_corr_png').format(event_key))
    pyplot.savefig(out_file)
    pyplot.close(fig)

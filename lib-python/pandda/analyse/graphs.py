import os, sys

#################################
try:
    import matplotlib
    matplotlib.interactive(False)
    from matplotlib import pyplot
    pyplot.style.use('ggplot')
except:
    pass
#################################

import numpy

from scitbx.array_family import flex

from bamboo.plot import bar, simple_histogram

def filter_nans(x):
    return [v for v in x if not numpy.isnan(v)]

def map_value_distribution(f_name, plot_vals, plot_normal=False):
    """Plot histogram of values, with optional normal distribution"""
    from scitbx.math.distributions import normal_distribution
    fig = pyplot.figure()
    pyplot.title('Distribution of map values')
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
    pyplot.xlabel('Map value')
    pyplot.ylabel('Density')
    fig.set_tight_layout(True)
    pyplot.savefig(f_name)
    pyplot.close(fig)

def qq_plot_against_normal(f_name, plot_vals):
    """Sort and plot list of values against expected quantiles from a normal distribution"""
    from scitbx.math.distributions import normal_distribution
    fig = pyplot.figure()
    pyplot.title('Q-Q plot for map values against normal distribution')
    expected_vals = normal_distribution().quantiles(len(plot_vals))
    pyplot.plot([min(expected_vals)-1, max(expected_vals)+1], [min(expected_vals)-1, max(expected_vals)+1], 'b--')
    pyplot.plot(sorted(plot_vals), expected_vals, 'go-')
    pyplot.xlabel('Observed quantiles')
    pyplot.ylabel('Theoretical quantiles')
    fig.set_tight_layout(True)
    pyplot.savefig(f_name)
    pyplot.close(fig)

def mean_obs_scatter(f_name, mean_vals, obs_vals):
    """Plot mean map against observed map values"""
    fig = pyplot.figure()
    pyplot.title('Mean map plotted against dataset map (unsorted)')
    pyplot.plot([-3, 10], [-3, 10], 'b--')
    pyplot.plot(mean_vals, obs_vals, 'go')
    pyplot.xlabel('Mean map value')
    pyplot.ylabel('Dataset map value')
    fig.set_tight_layout(True)
    pyplot.savefig(f_name)
    pyplot.close(fig)

def sorted_mean_obs_scatter(f_name, mean_vals, obs_vals):
    """Plot sorted mean map against sorted observed map values"""
    fig = pyplot.figure()
    pyplot.title('Mean map plotted against dataset map (sorted)')
    pyplot.plot([-3, 10], [-3, 10], 'b--')
    pyplot.plot(mean_vals, obs_vals, 'go')
    pyplot.xlabel('Mean map value')
    pyplot.ylabel('Dataset map value')
    fig.set_tight_layout(True)
    pyplot.savefig(f_name)
    pyplot.close(fig)

def uncertainty_qqplot(f_name, map_off, map_unc, q_cut, obs_diff, quantile):
    """Plot diff-mean map against normal quantiles, with uncertainty lines"""
    fig = pyplot.figure()
    pyplot.title('Q-Q plot for map values against normal distribution')
    pyplot.plot([map_off-5*map_unc, map_off+5*map_unc], [-5, 5], 'b--')
    pyplot.plot([-1, 1], [q_cut, q_cut], 'k-.')
    pyplot.plot([-1, 1], [-q_cut, -q_cut], 'k-.')
    pyplot.plot(obs_diff, quantile, 'go-')
    pyplot.xlabel('Difference from mean map')
    pyplot.ylabel('Theoretical Quantiles')
    fig.set_tight_layout(True)
    pyplot.savefig(f_name)
    pyplot.close(fig)

def write_occupancy_graph(f_name, x_values, global_values, local_values):
    """Write output graph from occupancy estimation"""

    # Get the x-value of the maximum difference
    diff_values = list(numpy.array(global_values) - numpy.array(local_values))
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
    # Remove spacing between subplots
    pyplot.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    pyplot.setp([axis_2_1.get_xticklabels()], visible=True)
    fig.set_tight_layout(True)
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
    fig.set_tight_layout(True)
    pyplot.savefig(pandda.file_manager.get_file('d_resolutions'))
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
        fig.set_tight_layout(True)
        pyplot.savefig(pandda.file_manager.get_file('d_rfactors'))
        pyplot.close(fig)
    except:
        fig = pyplot.figure()
        pyplot.title('Failed to make R-Factor Histograms')
        pyplot.savefig(pandda.file_manager.get_file('d_rfactors'))
        pyplot.close(fig)
    # ================================================>
    fig = pyplot.figure()
    pyplot.title('RMSDs to Reference Structure Histogram')
    pyplot.hist(x=filter_nans(d_info['rmsd_to_reference']), bins=n_bins)
    pyplot.xlabel('RMSD (A)')
    pyplot.ylabel('Count')
    fig.set_tight_layout(True)
    pyplot.savefig(pandda.file_manager.get_file('d_global_rmsd_to_ref'))
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
    fig.set_tight_layout(True)
    pyplot.savefig(pandda.file_manager.get_file('d_cell_axes'))
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
    fig.set_tight_layout(True)
    pyplot.savefig(pandda.file_manager.get_file('d_cell_angles'))
    pyplot.close(fig)
    # ================================================>
    fig = pyplot.figure()
    pyplot.title('Unit Cell Volume Variation')
    pyplot.hist(x=d_info['uc_vol'], bins=n_bins)
    pyplot.xlabel('Volume (A^3)')
    pyplot.ylabel('Count')
    fig.set_tight_layout(True)
    pyplot.savefig(pandda.file_manager.get_file('d_cell_volumes'))
    pyplot.close(fig)

def write_truncated_data_plots(pandda, resolution, miller_arrays):
    """Write summaries of the truncated data"""

    reslns = [m.d_min() for m in miller_arrays]
    min_res, max_res = min(reslns), max(reslns)
    pandda.log('After Truncation - Resolution Range: {!s}-{!s}'.format(min_res, max_res))

    # ================================================>
    # Resolution ranges of the truncated data
    simple_histogram(filename = pandda.file_manager.get_file('dataset_res_hist').format(resolution),
                     data     = reslns,
                     title    = 'Truncated dataset resolutions',
                     x_lab    = 'Resolution (A)',
                     n_bins   = 15)
    # ================================================>
    # Wilson plots of the truncated data
    v_min = v_max = 0
    # Make plot
    fig = pyplot.figure()
    pyplot.title('Wilson plot for truncated structure factors')
    for m in miller_arrays:
        # Convert to amplitudes
        ma = m.as_amplitude_array()
        # Create binning - just use the auto-binning for each dataset
        binner = ma.setup_binner(auto_binning=True)
        # Create the wilson plot
        binned = ma.wilson_plot(use_binning=True)
        # Extract bin centres and bin data
        bin_cent = binner.bin_centers(1)
        bin_data = binned.data[1:-1]
        assert len(bin_cent) == len(bin_data)
        pyplot.semilogy(bin_cent, bin_data, '-', linewidth=2, basey=numpy.e)
    pyplot.xlabel('1/resolution (A^-1)')
    pyplot.ylabel('Mean intensity')
    fig.set_tight_layout(True)
    pyplot.savefig(pandda.file_manager.get_file('dataset_wilson_plot').format(resolution))
    pyplot.close(fig)

def write_map_analyser_reference_dataset_graphs(pandda, map_analyser):

    # Plot the mean map against the reference map (unsorted)
    mean_obs_scatter(f_name      = pandda.file_manager.get_file('ref_v_mean_map_unsort').format(map_analyser.meta.resolution),
                     mean_vals   = map_analyser.statistical_maps.mean_map.as_sparse().data,
                     obs_vals    = pandda.datasets.reference().child.as_sparse().data)
    # Plot the mean map against the reference map (sorted)
    sorted_mean_obs_scatter(f_name    = pandda.file_manager.get_file('ref_v_mean_map_sort').format(map_analyser.meta.resolution),
                            mean_vals = sorted(map_analyser.statistical_maps.mean_map.as_sparse().data),
                            obs_vals  = sorted(pandda.datasets.reference().child.as_sparse().data))
    # Plot the reference map distribution
    map_value_distribution(f_name    = pandda.file_manager.get_file('ref_map_dist').format(map_analyser.meta.resolution),
                           plot_vals = pandda.datasets.reference().child.as_sparse().data)

def write_map_analyser_graphs(pandda, resolution, analysis_mask_name, building_mask_name):

    # Extract the statistical maps at this resolution
    statistical_maps = pandda.stat_maps.get(resolution)
    mean_map_vals = list(statistical_maps.mean_map.data)
    medn_map_vals = list(statistical_maps.medn_map.data)
    stds_map_vals = list(statistical_maps.stds_map.data)
    sadj_map_vals = list(statistical_maps.sadj_map.data)
    # Extract the dataset tags from the pandda object
    analysis_tags = [d.tag for d in pandda.datasets.mask(mask_name=analysis_mask_name)]
    building_tags = [d.tag for d in pandda.datasets.mask(mask_name=building_mask_name)]
    combined_tags = sorted(set(analysis_tags+building_tags))

    # Dataset Parameters
    d_info = pandda.tables.dataset_info
    m_info = pandda.tables.dataset_map_info

    # All datasets
    high_res = [d_info['high_resolution'][t] for t in combined_tags]
    low_res =  [d_info['low_resolution'][t]  for t in combined_tags]
    rfree =    [d_info['r_free'][t]          for t in combined_tags]
    rwork =    [d_info['r_work'][t]          for t in combined_tags]

    # All datasets
    map_uncties = [m_info['map_uncertainty'][t] for t in combined_tags]
    # Analysed datasets only
    z_map_mean  = [m_info['z_map_mean'][t] for t in analysis_tags]
    z_map_std   = [m_info['z_map_std'][t]  for t in analysis_tags]
    z_map_skew  = [m_info['z_map_skew'][t] for t in analysis_tags]
    z_map_kurt  = [m_info['z_map_kurt'][t] for t in analysis_tags]

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
    fig.set_tight_layout(True)
    pyplot.savefig(pandda.file_manager.get_file('map_mean_median_hist').format(resolution))
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
    fig.set_tight_layout(True)
    pyplot.savefig(pandda.file_manager.get_file('map_mean_median_scat').format(resolution))
    pyplot.close(fig)

    ##################################
    # MEAN-MEDIAN DIFFERENCE HISTOGRAM
    ##################################
    fig = pyplot.figure()
    pyplot.title('Mean-Median Difference Histogram')
    pyplot.hist(x=numpy.array(flex.double(mean_map_vals)-flex.double(medn_map_vals)), bins=30, normed=True)
    pyplot.xlabel('Difference Map Value')
    pyplot.ylabel('Density')
    fig.set_tight_layout(True)
    pyplot.savefig(pandda.file_manager.get_file('map_mean_median_diff').format(resolution))
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
    fig.set_tight_layout(True)
    pyplot.savefig(pandda.file_manager.get_file('map_stds_sadj_hist').format(resolution))
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
    fig.set_tight_layout(True)
    pyplot.savefig(pandda.file_manager.get_file('map_stds_sadj_scat').format(resolution))
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
    fig.set_tight_layout(True)
    pyplot.savefig(pandda.file_manager.get_file('dataset_unc_hist').format(resolution))
    pyplot.close(fig)

    ########################################################

    pandda.log('=> Scatter Plots')

    # MAP RESOLUTION V UNCERTAINTY
    fig = pyplot.figure()
    pyplot.title('High Resolution Limit v Map Uncertainty')
    pyplot.scatter(x=high_res, y=map_uncties)
    pyplot.xlabel('Dataset Resolution')
    pyplot.ylabel('Dataset Map Uncertainty')
    fig.set_tight_layout(True)
    pyplot.savefig(pandda.file_manager.get_file('dataset_res_unc_scat').format(resolution))
    pyplot.close(fig)

    # MAP RESOLUTION V UNCERTAINTY
    fig = pyplot.figure()
    pyplot.title('High Resolution Limit v R-Free')
    pyplot.scatter(x=high_res, y=rfree)
    pyplot.xlabel('Dataset Resolution')
    pyplot.ylabel('Dataset R-Free')
    fig.set_tight_layout(True)
    pyplot.savefig(pandda.file_manager.get_file('dataset_res_rfree_scat').format(resolution))
    pyplot.close(fig)

    # RFREE V UNCERTAINTY
    fig = pyplot.figure()
    pyplot.title('R-Free v Uncertainty')
    pyplot.scatter(x=rfree, y=map_uncties)
    pyplot.xlabel('Dataset R-Free')
    pyplot.ylabel('Dataset Map Uncertainty')
    fig.set_tight_layout(True)
    pyplot.savefig(pandda.file_manager.get_file('dataset_unc_rfree_scat').format(resolution))
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
    fig.set_tight_layout(True)
    pyplot.savefig(pandda.file_manager.get_file('zmap_mean_sadj_hist').format(resolution))
    pyplot.close(fig)

    # Z-MAP SKEW V UNCERTAINTY
    fig = pyplot.figure()
    pyplot.title('Z-Map Normality Plots')
    pyplot.scatter(x=z_map_skew, y=z_map_kurt)
    pyplot.xlabel('Skew')
    pyplot.ylabel('Kurtosis')
    fig.set_tight_layout(True)
    pyplot.savefig(pandda.file_manager.get_file('zmap_skew_kurt_scat').format(resolution))
    pyplot.close(fig)

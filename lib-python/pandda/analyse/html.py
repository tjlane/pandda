import os, glob

from bamboo.html import png2base64str
from pandda.html import PANDDA_HTML_ENV

def write_initial_html(pandda):
    # Get template to be filled in
    template = PANDDA_HTML_ENV.get_template('pandda_summary.html')
    # Output file
    out_file = pandda.file_manager.get_file(file_tag='initial_html')
    # Output directory (for relative symlinks)
    out_dir  = os.path.abspath(os.path.dirname(out_file))

    # ===========================================================>
    # Construct the data object to populate the template
    output_data = {}
    output_data['header'] = 'PANDDA Dataset Summaries'
    output_data['title'] = 'PANDDA Dataset Summaries'
    output_data['introduction'] = 'Summary of Added Datasets'
    # ===========================================================>
    # Summary bar
    output_data['summary_bar'] = []
    output_data['summary_bar'].append({'colour':'info',    'text':'Datasets Loaded:   {}'.format(pandda.datasets.size())                                            })
    output_data['summary_bar'].append({'colour':'success', 'text':'Datasets Accepted: {}'.format(pandda.datasets.size(mask_name='rejected - total', invert=True))   })
    output_data['summary_bar'].append({'colour':'danger',  'text':'Datasets Rejected: {}'.format(pandda.datasets.size(mask_name='rejected - total'))                })
    # ===========================================================>
    # Header Images
    output_data['top_images'] = []
    for file_tag, title, tooltip in [('d_resolutions',
                                      'Dataset Resolutions', ''),
                                     ('d_rfactors',
                                      'Dataset R-Factors', ''),
                                     ('d_unscaled_wilson_plots',
                                      'Unscaled Wilson plots', ''),
                                     ('d_scaled_wilson_plots',
                                      'Scaled Wilson plots', ''),
                                     ('d_unscaled_wilson_rmsds',
                                      'RMSDs to reference wilson plot (unscaled)', ''),
                                     ('d_scaled_wilson_rmsds',
                                      'RMSDs to reference wilson plot (scaled)', ''),
                                     ('d_global_rmsd_to_ref',
                                      'Dataset RMSD to Mean Structure', ''),
                                     ('d_cell_volumes',
                                      'Dataset Cell Volumes', ''),
                                     ('d_cell_axes',
                                      'Dataset Cell Axis Lengths', ''),
                                     ('d_cell_angles',
                                      'Dataset Cell Angles', '')]:
        f_name = pandda.file_manager.get_file(file_tag=file_tag)
        output_data['top_images'].append({ 'path': 'data:image/png;base64,{}'.format(png2base64str(path=f_name)),
                                           'title': title,
                                           'tooltip': tooltip  })
    # ===========================================================>
    # Write Output
    with open(out_file, 'w') as out_html:
        out_html.write(template.render(output_data))

def write_datasets_htmls(pandda, resolution, datasets):

    # Get template to be filled in
    template = PANDDA_HTML_ENV.get_template('pandda_summary.html')

    # Iterate through datasets and write graphs for each
    for dataset in datasets:
        # Output file
        out_file = dataset.file_manager.get_file(file_tag='dataset_html')

        # ===========================================================>
        # Construct the data object to populate the template
        output_data = {}
        output_data['header'] = 'Dataset {} Analysis Summary'.format(dataset.tag)
        output_data['title'] = 'Dataset {} Analysis Summary'.format(dataset.tag)
        output_data['introduction'] = 'Summary of dataset {} analysed @ {}A'.format(dataset.tag, resolution)
        # ===========================================================>
        # Dataset analysis images
        output_data['top_images'] = []
        for file_tag, title, tooltip in [('s_map_png',
                                          '<strong>Electron density distribution</strong>. ',
                                          'Tooltip. '\
                                          ''),
                                         ('d_mean_map_png',
                                          '<strong>Distribution of differences from the mean map</strong>. ',
                                          'Tooltip. '\
                                          ''),
                                         ('obs_qqplot_unsorted_png',
                                          '<strong>Dataset v Reference map values (unsorted)</strong>. ',
                                          'Tooltip. '\
                                          ''),
                                         ('obs_qqplot_sorted_png',
                                          '<strong>Dataset v Reference map values (sorted)</strong>. ',
                                          'Tooltip. '\
                                          ''),
                                         ('unc_qqplot_png',
                                          '<strong>Uncertainty Q-Q plot</strong>. ',
                                          'Tooltip. '\
                                          ''),
                                         ('z_map_qq_plot_png',
                                          '<strong>Z-map Q-Q plot</strong>. ',
                                          'Tooltip. '\
                                          ''),
#                                         ('z_map_naive_normalised_png',
#                                          '<strong>(Naive) Z-map distribution</strong>. ',
#                                          'Tooltip. '\
#                                          ''),
                                         ('z_map_corrected_normalised_png',
                                          '<strong>(Corrected) Z-map distribution</strong>. ',
                                          'Tooltip. '\
                                          ''),
                                         ('wilson_plot_png',
                                          '<strong>Wilson Plot. ',
                                          'Tooltip. '\
                                          '')]:
            f_name = dataset.file_manager.get_file(file_tag=file_tag)
            output_data['top_images'].append({ 'path': 'data:image/png;base64,{}'.format(png2base64str(path=f_name)),
                                               'title': title,
                                               'tooltip': tooltip  })
        # Any event maps (variable filenames)
        for f_name in sorted(glob.glob(dataset.file_manager.get_file(file_tag='bdc_est_png').format('*'))):
            output_data['top_images'].append({ 'path': 'data:image/png;base64,{}'.format(png2base64str(path=f_name)),
                                               'title': os.path.basename(os.path.splitext(f_name)[0]),
                                               'tooltip': '' })
        # ===========================================================>
        # Write Output
        with open(out_file, 'w') as out_html:
            out_html.write(template.render(output_data))

def write_map_analyser_html(pandda, resolution):

    # Get template to be filled in
    template = PANDDA_HTML_ENV.get_template('pandda_summary.html')
    # Output file
    out_file = pandda.file_manager.get_file(file_tag='map_html').format(resolution)
    # Output directory (for relative symlinks)
    out_dir  = os.path.abspath(os.path.dirname(out_file))

    # ===========================================================>
    # Construct the data object to populate the template
    output_data = {}
    output_data['header'] = 'PANDDA Statistical Map Summary'
    output_data['title'] = 'PANDDA Statistical Map Summary'
    output_data['introduction'] = 'Summary of Statistical Maps @ {}A'.format(resolution)
    # ===========================================================>
    # Summary bar
#    output_data['summary_bar'] = []
#    output_data['summary_bar'].append({'colour':'info',    'text':'Datasets Loaded:   {}'.format(pandda.datasets.size())                                            })
#    output_data['summary_bar'].append({'colour':'success', 'text':'Datasets Accepted: {}'.format(pandda.datasets.size(mask_name='rejected - total', invert=True))   })
#    output_data['summary_bar'].append({'colour':'danger',  'text':'Datasets Rejected: {}'.format(pandda.datasets.size(mask_name='rejected - total'))                })
    # ===========================================================>
    # Header Images
    output_data['top_images'] = []
    for file_tag, title, tooltip in [('truncated_res_hist',
                                      '<strong>Distribution of truncated dataset resolutions</strong>. ',
                                      'Histogram of the diffraction data resolution for each dataset after truncation to the same resolution. '\
                                      'The truncated data is used to generate the statistical maps at this resolution, so for best results, '\
                                      'the range on the x-axis should be very small. Large variation indicates that something has gone wrong during '\
                                      'data truncation.'),
                                     ('truncated_wilson_plot',
                                      '<strong>Wilson plots of truncated diffraction data</strong>.',
                                      'Plots of intensity vs resolution for the truncated diffraction data, showing how intensity falls off as a '\
                                      'function of resolution. For best results, all of the lines should be approximately overlaid on each other.'),
                                     ('ref_v_mean_map_unsort',
                                      '<strong>Comparison of electron density values in the reference and mean map (unsorted)</strong>.',
                                      'Scatter plot of equivalent electron density values in the reference map and the mean map (electron '\
                                      'density values sampled at the same point in both maps). The points in this plot should group along '\
                                      'the diagonal; the amount that the cloud deviates from the diagonal reflects the uncertainty in the '\
                                      'reference map and the variation in the electron density data.'),
                                     ('ref_v_mean_map_sort',
                                      '<strong>Comparison of electron density values in the reference and mean map (after sorting)</strong>.',
                                      'Scatter plot of sorted electron density values in the reference map and the mean map (largest in one '\
                                      'plotted against largest in other, and smallest in one plotted against smallest in the other, etc.), '\
                                      'resulting in a quantile-quantile plot. The results should be a line along the diagonal. Deviations from '\
                                      'the diagonal, except at the extremes of the axis may indicate problems in analysis.'),
                                     ('ref_map_dist',
                                      '<strong>Reference map distribution</strong>.',
                                      'Histogram of electron density values in the reference map. This should be the same shape as those for the '\
                                      'mean and median maps (a slightly skewed distribution).'),
                                     ('map_mean_median_hist',
                                      '<strong>Mean & median map distributions</strong>.',
                                      'Histogram of electron density values in the mean (average) and median electron density maps. These should '\
                                      'be the same shape as for the reference map.'),
                                     ('map_mean_median_scat',
                                      '<strong>Comparison of electron density values in the mean and median map</strong>.',
                                      'Scatter plot of electron density values at equivalent points in both maps.'),
                                     ('map_mean_median_diff',
                                      '<strong>Distribution of differences between the mean and median map</strong>.',
                                      'Histograms of median map subtracted from the mean map. Points should cluster around the diagonal'),
                                     ('map_stds_sadj_hist',
                                      '<strong>Distribution of values in the "raw" and "adjusted" variation maps</strong>.',
                                      'The "raw" variation at a point is the standard deviation of the observed electron density values at that point;'\
                                      'the "adjusted" variation accounts for the error associated with each electron density value. The graphs show '\
                                      'the histogram of these values for the whole map. The "adjusted" variation histogram should show smaller '\
                                      'values on average than the "raw" variation (the distribution should be shifted to the left)'),
                                     ('map_stds_sadj_scat',
                                      '<strong>Comparison of "raw" and "adjusted" variation for all points</strong>.',
                                      'Scatter plot of the different variation measures for all points. The "adjusted" variation should always be '\
                                      'lower than the "raw" variation, so all points on the plot should fall to the bottom-right of the diagonal.'),
                                     ('dataset_unc_hist',
                                      '<strong>Distribution of map uncertainties for analysed datasets</strong>.',
                                      'Histogram of datasets when compared to the mean map.'),
                                     ('dataset_res_unc_scat',
                                      '<strong>Comparison of dataset resolution and map uncertainty</strong>.',
                                      'Resolution of the full (untruncated) diffraction data for each dataset plotted against the uncertainty of the '\
                                      'map calculated at this resolution. There is rarely a strong correlation between these variables.'),
                                     ('dataset_res_rfree_scat',
                                      '<strong>Comparison of dataset resolution and R-free</strong>.',
                                      'Resolution of the full (untruncated) diffraction data for each dataset plotted against the R-free of the input model. '\
                                      'There is rarely a strong correlation between these variables.'),
                                     ('dataset_unc_rfree_scat',
                                      '<strong>Comparison of map uncertainty and dataset R-Free</strong>.',
                                      'Uncertainty of the map calculated at this resolution plotted against the R-free of the input model. There is rarely '\
                                      'a strong correlation between these variables.'),
                                     ('zmap_mean_sadj_hist',
                                      '<strong>Distribution of mean Z-values and standard deviation of Z-values for each dataset</strong>.',
                                      'The mean values should be approximately zero, and the standard deviations should be approximately one (although they '\
                                      'frequently are observed to be between 0.8-1.0.'),
                                     ('zmap_skew_kurt_scat',
                                      '<strong>Comparision of Z-map skew and kurtosis for each dataset</strong>.',
                                      'Scatter plot of skew and kurtosis for all z-maps. Skew should be approximately zero and kurtosis should be '\
                                      'approximately 3.0.')]:
        f_name = pandda.file_manager.get_file(file_tag=file_tag).format(resolution)
        output_data['top_images'].append({ 'path': 'data:image/png;base64,{}'.format(png2base64str(path=f_name)),
                                           'title': title,
                                           'tooltip': tooltip  })
    # ===========================================================>
    # Write Output
    with open(out_file, 'w') as out_html:
        out_html.write(template.render(output_data))

def write_analyse_html(pandda):

    # Get template to be filled in
    template = PANDDA_HTML_ENV.get_template('pandda_summary.html')
    # Output file
    out_file = pandda.file_manager.get_file(file_tag='analyse_html')
    # Output directory (for relative symlinks)
    out_dir  = os.path.abspath(os.path.dirname(out_file))

    # Extract the dataset info as a dictionary
    all_d_info = pandda.tables.dataset_info.transpose().to_dict()
    all_m_info = pandda.tables.dataset_map_info.transpose().to_dict()

    len_data = len(pandda.tables.dataset_info.index)
    num_interesting = len([d for d in pandda.datasets.mask(mask_name='rejected - total', invert=True) if d.events])
    num_analysed    = len([d for d in pandda.datasets.mask(mask_name='rejected - total', invert=True) if d.meta.analysed])
    num_rejected    = pandda.datasets.size(mask_name='rejected - total')
    num_not_interesting = num_analysed - num_interesting
    num_not_analysed    = len_data - num_analysed - num_rejected

    # ===========================================================>
    # Construct the data object to populate the template
    output_data = {}
    output_data['header'] = 'PANDDA Processing Output'
    output_data['title'] = 'PANDDA Processing Output'
    output_data['introduction'] = 'Summary of Processing of Datasets'
    # ===========================================================>
    # Summary bar
    output_data['summary_bar'] = []
    output_data['summary_bar'].append({'colour':'info',    'text':'Analysed: {}'.format(num_analysed)               })
    output_data['summary_bar'].append({'colour':'success', 'text':'Interesting: {}'.format(num_interesting)         })
    output_data['summary_bar'].append({'colour':'warning', 'text':'Not Interesting: {}'.format(num_not_interesting) })
    output_data['summary_bar'].append({'colour':'danger',  'text':'Not Analysed: {}'.format(num_not_analysed)       })
    # Progress Bars
    output_data['progress_bar'] = []
    output_data['progress_bar'].append({'title':'Dataset Summary', 'data':[]})
    output_data['progress_bar'][0]['data'].append({'text':'Interesting',     'colour':'success', 'size':100.0*num_interesting/len_data})
    output_data['progress_bar'][0]['data'].append({'text':'Not Interesting', 'colour':'info',    'size':100.0*num_not_interesting/len_data})
    output_data['progress_bar'][0]['data'].append({'text':'Not Analysed',    'colour':'warning', 'size':100.0*num_not_analysed/len_data})
    output_data['progress_bar'][0]['data'].append({'text':'Rejected',        'colour':'danger',  'size':100.0*num_rejected/len_data})
    # ===========================================================>
    # Header Images
    output_data['top_images'] = []
    if os.path.exists(pandda.file_manager.get_file(file_tag='pymol_sites_png_1')):
        output_data['top_images'].append({ 'path': 'data:image/png;base64,{}'.format(png2base64str(path=pandda.file_manager.get_file(file_tag='pymol_sites_png_1'))),
                                           'title': 'Identified Sites (Front)' })
    if os.path.exists(pandda.file_manager.get_file(file_tag='pymol_sites_png_2')):
        output_data['top_images'].append({ 'path': 'data:image/png;base64,{}'.format(png2base64str(path=pandda.file_manager.get_file(file_tag='pymol_sites_png_2'))),
                                           'title': 'Identified Sites (Back)' })
    for i_png, png in enumerate(sorted(glob.glob(pandda.file_manager.get_file(file_tag='analyse_site_graph_mult').format('*')))):
        output_data['top_images'].append({ 'path': 'data:image/png;base64,{}'.format(png2base64str(path=png)),
                                           'title': 'Identified Site Summary ({})'.format(i_png+1) })
    # ===========================================================>
    # Tables
    output_data['table'] = {}
    output_data['table']['column_headings'] = ['Dataset', 'Dataset Checks', 'RFree/RWork', 'Resolution', 'Reference RMSD', 'Map Uncertainty', 'Overall', 'Interesting Areas', 'Result']
    output_data['table']['rows'] = []
    # Add the datasets as rows
    for d in pandda.datasets.all():

        d_data = all_d_info[d.tag]
        m_data = all_m_info[d.tag]
        rmsd  = d_data['rmsd_to_reference']
        rfree = d_data['r_free']
        rwork = d_data['r_work']
        res_h = d_data['high_resolution']
        uncty = m_data['map_uncertainty']

        rejection_reason = d_data['rejection_reason']

        # colour choices - 'success', 'muted', 'danger'
        # icon choices   - 'ok', 'flag', 'remove'

        columns = []
        # ------------------------------>>>
        # Whether dataset was rejected
        # ------------------------------>>>
        if pandda.datasets.all_masks().get_value(name='rejected - total', id=d.tag) == True:
            columns.append({'colour':'danger',  'icon':'remove',    'message':rejection_reason})
        else:
            columns.append({'colour':'success', 'icon':'ok',        'message':'OK'})
        # ------------------------------>>>
        # Test for Refinement Success - some test on r-free
        # ------------------------------>>>
        if pandda.datasets.all_masks().get_value(name='rejected - rfree', id=d.tag) == True:
            columns.append({'colour':'danger',  'icon':'remove',    'message':'{}/{}'.format(rfree,rwork)})
        else:
            columns.append({'colour':'success', 'icon':'default',   'message':'{}/{}'.format(rfree,rwork)})
        # ------------------------------>>>
        # Resolution
        # ------------------------------>>>
        if 0:
            columns.append({'colour':'danger',  'icon':'remove',    'message':'{}'.format(res_h)})
        else:
            columns.append({'colour':'success', 'icon':'ok',        'message':'{}'.format(res_h)})
        # ------------------------------>>>
        # Test for Structure movement
        # ------------------------------>>>
        if pandda.datasets.all_masks().get_value(name='rejected - rmsd to reference', id=d.tag) == True:
            columns.append({'colour':'danger',  'icon':'remove',    'message':'{}'.format(rmsd)})
        else:
            columns.append({'colour':'success', 'icon':'ok',        'message':'{}'.format(rmsd)})
        # ------------------------------>>>
        # Uncertainty of Map
        # ------------------------------>>>
        if 0:
            columns.append({'colour':'danger',  'icon':'remove',    'message':'{}'.format(uncty)})
        else:
            columns.append({'colour':'success', 'icon':'ok',        'message':'{}'.format(uncty)})
        # ------------------------------>>>
        # Analysed
        # ------------------------------>>>
        if d.meta.analysed == True:
            columns.append({'colour':'success', 'icon':'ok',        'message':'Analysed'})
        else:
            columns.append({'colour':'danger',  'icon':'remove',    'message':'Not Analysed'})
        # ------------------------------>>>
        # Test for if it's interesting
        # ------------------------------>>>
        if d.events:
            columns.append({'colour':'success', 'icon':'ok',        'message':'{} Clusters'.format(len(d.events))})
        else:
            columns.append({'colour':'default', 'icon':'remove',    'message':''})

        row_message = 'Z-Blobs Found' if d.events else \
                      ''
        row_colour  = 'success' if row_message else \
                      'danger' if pandda.datasets.all_masks().get_value(name='rejected - total', id=d.tag) == True else \
                      'default'

        output_data['table']['rows'].append({'heading' : d.tag,
                                             'colour'  : row_colour,
                                             'message' : row_message,
                                             'columns' : columns})

    # ===========================================================>
    # Write Output
    with open(out_file, 'w') as out_html:
        out_html.write(template.render(output_data))

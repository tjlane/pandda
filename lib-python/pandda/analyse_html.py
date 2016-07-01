import os, glob

from bamboo.html import png2base64str
from pandda.html import PANDDA_HTML_ENV

def write_initial_html(pandda):
    # Get template to be filled in
    template = PANDDA_HTML_ENV.get_template('pandda_summary.html')
    # Output file
    out_file = pandda.output_handler.get_file(file_tag='initial_html')
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
    output_data['top_images'].append({ 'path': 'data:image/png;base64,{}'.format(png2base64str(path=pandda.output_handler.get_file(file_tag='d_resolutions'))),
                                       'title': 'Dataset Resolutions' })
    output_data['top_images'].append({ 'path': 'data:image/png;base64,{}'.format(png2base64str(path=pandda.output_handler.get_file(file_tag='d_rfactors'))),
                                       'title': 'Dataset R-Factors' })
    output_data['top_images'].append({ 'path': 'data:image/png;base64,{}'.format(png2base64str(path=pandda.output_handler.get_file(file_tag='d_global_rmsd_to_ref'))),
                                       'title': 'Dataset RMSD to Mean Structure' })
    output_data['top_images'].append({ 'path': 'data:image/png;base64,{}'.format(png2base64str(path=pandda.output_handler.get_file(file_tag='d_cell_volumes'))),
                                       'title': 'Dataset Cell Volumes' })
    output_data['top_images'].append({ 'path': 'data:image/png;base64,{}'.format(png2base64str(path=pandda.output_handler.get_file(file_tag='d_cell_axes'))),
                                       'title': 'Dataset Cell Axis Lengths' })
    output_data['top_images'].append({ 'path': 'data:image/png;base64,{}'.format(png2base64str(path=pandda.output_handler.get_file(file_tag='d_cell_angles'))),
                                       'title': 'Dataset Cell Angles' })
    # ===========================================================>
    # Write Output
    with open(out_file, 'w') as out_html:
        out_html.write(template.render(output_data))

def write_analyse_html(pandda):

    # Get template to be filled in
    template = PANDDA_HTML_ENV.get_template('pandda_summary.html')
    # Output file
    out_file = pandda.output_handler.get_file(file_tag='analyse_html')
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
    if os.path.exists(pandda.output_handler.get_file(file_tag='pymol_sites_png_1')):
        output_data['top_images'].append({ 'path': 'data:image/png;base64,{}'.format(png2base64str(path=pandda.output_handler.get_file(file_tag='pymol_sites_png_1'))),
                                           'title': 'Identified Sites (Front)' })
    if os.path.exists(pandda.output_handler.get_file(file_tag='pymol_sites_png_2')):
        output_data['top_images'].append({ 'path': 'data:image/png;base64,{}'.format(png2base64str(path=pandda.output_handler.get_file(file_tag='pymol_sites_png_2'))),
                                           'title': 'Identified Sites (Back)' })
    for i_png, png in enumerate(sorted(glob.glob(pandda.output_handler.get_file(file_tag='analyse_site_graph_mult').format('*')))):
        output_data['top_images'].append({ 'path': 'data:image/png;base64,{}'.format(png2base64str(path=png)),
                                           'title': 'Identified Site Summary ({})'.format(i_png+1) })
    # ===========================================================>
    # Tables
    output_data['table'] = {}
    output_data['table']['column_headings'] = ['Dataset', 'Crystal', 'Structure', 'RFree/RWork', 'Resolution', 'Reference RMSD', 'Map Uncertainty', 'Overall', 'Interesting Areas', 'Result']
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
        # Test for Bad Crystal
        # ------------------------------>>>
        if pandda.datasets.all_masks().get_mask_value(mask_name='rejected - crystal', entry_id=d.tag) == True:
            columns.append({'colour':'danger',  'icon':'remove',    'message':rejection_reason})
        else:
            columns.append({'colour':'success', 'icon':'ok',        'message':'OK'})
        # ------------------------------>>>
        # Test for Bad Structures
        # ------------------------------>>>
        if pandda.datasets.all_masks().get_mask_value(mask_name='rejected - structure', entry_id=d.tag) == True:
            columns.append({'colour':'danger',  'icon':'remove',    'message':rejection_reason})
        else:
            columns.append({'colour':'success', 'icon':'ok',        'message':'OK'})
        # ------------------------------>>>
        # Test for Refinement Success - some test on r-free
        # ------------------------------>>>
        if pandda.datasets.all_masks().get_mask_value(mask_name='bad crystal - rfree', entry_id=d.tag) == True:
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
        if pandda.datasets.all_masks().get_mask_value(mask_name='bad crystal - isomorphous structure', entry_id=d.tag) == True:
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
                      'danger' if pandda.datasets.all_masks().get_mask_value(mask_name='rejected - total', entry_id=d.tag) == True else \
                      'default'

        output_data['table']['rows'].append({'heading' : d.tag,
                                             'colour'  : row_colour,
                                             'message' : row_message,
                                             'columns' : columns})

    # ===========================================================>
    # Write Output
    with open(out_file, 'w') as out_html:
        out_html.write(template.render(output_data))

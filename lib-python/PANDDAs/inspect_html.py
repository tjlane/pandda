import os

from PANDDAs.settings import PANDDA_TOP, PANDDA_TEXT
from PANDDAs.html import PANDDA_HTML_ENV, path2url

def write_inspect_html(out_dir, inspector):

    # Get template to be filled in
    template = PANDDA_HTML_ENV.get_template('pandda_summary.html')

    all_data = inspector.log_table
    len_data = len(all_data.index)

    num_fitted = sum(all_data['Ligand Placed'])
    num_viewed = sum(all_data['Viewed'])
    num_empty  = num_viewed - num_fitted
    num_unviewed = len_data - num_viewed

    # ===========================================================>
    # Construct the data object to populate the template
    output_data = {'PANDDA_TOP' : path2url(PANDDA_TOP)}
    output_data['header'] = 'PANDDA Inspect Summary'
    output_data['title'] = 'PANDDA Inspect Summary'
    output_data['introduction'] = 'Summary of Inspection of Datasets'
    # ===========================================================>
    # Header Images
    output_data['top_images'] = []
    output_data['top_images'].append({ 'path': path2url(os.path.abspath(os.path.join(out_dir, 'results_summaries', 'pandda_inspect.png'))),
                                       'title': 'Identified Site Summary' })
    # ===========================================================>
    # Progress Bars
    output_data['progress_bar'] = []
    output_data['progress_bar'].append({'title':'Fitting Progress', 'data':[]})
    output_data['progress_bar'][0]['data'].append({'text':'Fitted',  'colour':'success', 'size':100.0*num_fitted/len_data})
    output_data['progress_bar'][0]['data'].append({'text':'Unviewed','colour':'warning', 'size':100.0*num_unviewed/len_data})
    output_data['progress_bar'][0]['data'].append({'text':'Empty',   'colour':'danger',  'size':100.0*num_empty/len_data})
    # ===========================================================>
    # Tables
    output_data['table'] = {}
    output_data['table']['column_headings'] = ['Event','Site','Est. Occupancy','Z-Peak','Map Resolution','Map Uncertainty','Interesting','Ligand Placed','Ligand Confidence','Comment','Viewed']
    output_data['table']['rows'] = []
    # Add the datasets as rows
    for i_d in range(len(all_data.index)):

        # colour choices - 'success', 'muted', 'danger'
        # icon choices   - 'ok', 'flag', 'remove'

        d_data = all_data.iloc[i_d].to_dict()
        d_tag, d_event = all_data.index[i_d]

        columns = []
        columns.append({'message':d_event})
        columns.append({'message':d_data['site_idx']})
        columns.append({'message':round(d_data['est_occupancy'],3)})
        columns.append({'message':round(d_data['z_peak'],3)})
        columns.append({'message':d_data['analysed_resolution']})
        columns.append({'message':round(d_data['map_uncertainty'],3)})

        if d_data['Interesting']:   columns.append({'colour':'success', 'icon':'ok',     'message':d_data['Interesting']})
        else:                       columns.append({'colour':'danger',  'icon':'remove', 'message':d_data['Interesting']})

        if d_data['Ligand Placed']: columns.append({'colour':'success', 'icon':'ok',     'message':d_data['Ligand Placed']})
        else:                       columns.append({'colour':'danger',  'icon':'remove', 'message':d_data['Ligand Placed']})

        columns.append({'message':d_data['Ligand Confidence']})
        columns.append({'message':d_data['Comment']})

        if d_data['Viewed']:        columns.append({'colour':'success', 'icon':'ok',     'message':d_data['Viewed']})
        else:                       columns.append({'colour':'danger',  'icon':'remove', 'message':d_data['Viewed']})

        row_message = 'Ligand' if d_data['Ligand Placed'] else \
                      ''
        row_colour  = 'success' if d_data['Ligand Placed'] else \
                      'danger' if not d_data['Viewed'] else \
                      'info'

        output_data['table']['rows'].append({'heading' : d_tag,
                                             'colour'  : row_colour,
                                             'message' : row_message,
                                             'columns' : columns})

    with open(os.path.join(out_dir, 'results_summaries', 'pandda_inspect.html'), 'w') as out_html:
        out_html.write(template.render(output_data))
import os

from PANDDAs.settings import PANDDA_TOP, PANDDA_TEXT
from PANDDAs.html import PANDDA_HTML_ENV, path2url

def write_inspect_html(out_dir, inspector):

    # Get template to be filled in
    template = PANDDA_HTML_ENV.get_template('pandda_summary.html')

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
    # Tables
    output_data['table'] = {}
    output_data['table']['column_headings'] = ['Event','Site','Est. Occupancy','Z-Peak','Map Resolution','Map Uncertainty','Interesting','Ligand Placed','Ligand Confidence','Comment','Viewed']
    output_data['table']['rows'] = []
    # Add the datasets as rows
    for i_d in range(len(inspector.log_table.index)):

        d_data = inspector.log_table.iloc[i_d].to_dict()
        d_tag, d_event = inspector.log_table.index[i_d]

        columns = []
        columns.append({'message':d_event})
        columns.append({'message':d_data['site_idx']})
        columns.append({'message':round(d_data['est_occupancy'],3)})
        columns.append({'message':round(d_data['z_peak'],3)})
        columns.append({'message':d_data['analysed_resolution']})
        columns.append({'message':round(d_data['map_uncertainty'],3)})

        if d_data['Interesting']:   columns.append({'flag':0,'message':d_data['Interesting']})
        else:                       columns.append({'flag':4,'message':d_data['Interesting']})

        if d_data['Ligand Placed']: columns.append({'flag':0,'message':d_data['Ligand Placed']})
        else:                       columns.append({'flag':4,'message':d_data['Ligand Placed']})

        columns.append({'message':d_data['Ligand Confidence']})
        columns.append({'message':d_data['Comment']})

        if d_data['Viewed']:        columns.append({'flag':0,'message':d_data['Viewed']})
        else:                       columns.append({'flag':4,'message':d_data['Viewed']})

        row_colour = 0 if d_data['Ligand Placed'] else 3

        output_data['table']['rows'].append({'heading':d_tag,
                                             'success':row_colour,
                                             'columns':columns})

    with open(os.path.join(out_dir, 'results_summaries', 'pandda_inspect.html'), 'w') as out_html:
        out_html.write(template.render(output_data))

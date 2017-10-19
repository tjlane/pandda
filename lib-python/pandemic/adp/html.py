import os, glob

from bamboo.html import png2base64str
from pandemic.html import PANDEMIC_HTML_ENV


def write_adp_summary(fitter):
    """Write an overall, and a level-by-level summary of the parameterised ADPs"""
    # Get the html template
    template = PANDDA_HTML_ENV.get_template('adp_summary.html')
    # Create output file path
    out_file = os.path.abspath(os.path.join(fitter.out_dir, 'results.html'))

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



    with open(out_file, 'w') as out_html:
        out_html.write(template.render(output_data))

import os
from libtbx import adopt_init_args
from pandemic.adp.html import HtmlSummary


class PostProcessingHtmlSummary(HtmlSummary):


    r_factor_keys = ['r_free', 'r_work', 'r_gap']

    def __init__(self,
        results_object,
        analysis_files,
        ):
        _json_plots = []
        adopt_init_args(self, locals())

    def json_plots(self):
        return self._json_plots

    def main_summary(self):

        tab = {'id'        : 'statistics',
               'alt_title' : 'R-factor Statistics',
               'title' : 'Model R-factor Summary',
               'contents'  : [],
               }

        tab['contents'].append(self.get_interactive_results_plot())
        tab['contents'].append(self.get_dataset_statistics_table())

        panels = self.get_rfactor_panels()
        if panels:
            tab['contents'].append({'title' : 'R-factor statistics'})
            tab['contents'].extend(panels)

        return [tab]

    def get_rfactor_panels(self):

        panels = []

        f = self.analysis_files.get('r_values')
        if f:
            images = []
            images.append({
                'width' : 6,
                'title' : 'R-Free Comparison',
                'image' : self.image(f.get('r_free')),
                'footnote': "Comparison of R-free values for the different input and output structures (and optionally reference R-values).",
                })
            images.append({
                'width' : 6,
                'title' : 'R-Work Comparison',
                'image' : self.image(f.get('r_work')),
                'footnote': "Comparison of R-work values for the different input and output structures (and optionally reference R-values).",
                })
            images.append({
                'width' : 6,
                'title' : 'R-Gap Comparison',
                'image' : self.image(f.get('r_gap')),
                'footnote' : "Comparison of R-work values for the different input and output structures (and optionally reference R-values).",
                })
            panels.append({
                'type' : 'panel',
                'title' : 'R-factor comparisons across structures',
                'contents' : images,
                })

        f = self.analysis_files.get('r_value_differences')
        if f:
            images = []
            for (s1, s2), img in f.iteritems():
                images.append({
                    'width' : 6,
                    'title' : 'R-factor differences between {} and {} structures'.format(s1.lower(), s2.lower()),
                    'image' : self.image(img),
                    'footnote' : "Changes in the R-factors from the {} model to the {} model (negative values mean lower values for the {} model).".format(s2.lower(), s1.lower(), s1.lower()),
                    })
            panels.append({
                'type' : 'panel',
                'title' : 'R-factor differences between structures',
                'contents' : images,
                })

        return panels

    def get_interactive_results_plot(self):

        # Extract data table for plot and table
        table = self.results_object.table.dropna(axis='columns', how='all')
        # Panel for the interactive plots
        json_plot = {
            'id'        : 'variable-plots-{}'.format(self.counter.next()),
            'json'      : table.T.to_json(orient='split'),
            'default_x' : 'High Resolution Limit',
            'default_y' : 'R-free Diff. (Fitted-Input)',
            }
        self._json_plots.append(json_plot)
        panel = {
            'type' : 'panel',
            'title' : 'Interactive Summary Graphs',
            'contents' : [{'id': json_plot['id']}],
            }

        return panel

    def get_dataset_statistics_table(self):

        table = self.results_object.table.dropna(axis='columns', how='all')
        panel = {
            'type'  : 'panel',
            'title' : 'Individual Dataset Statistics Table (click to expand/collapse)',
            'show'  : False,
            'contents'  : [
                {
                    'text': 'Data from output CSV',
                    'table': table.to_html(bold_rows=False, classes=['table table-striped table-hover datatable nowrap'])\
                               .replace('<th></th>','<th>Dataset</th>')\
                               .replace('border="1" ', ''),
                    },
                ],
            }

        return panel



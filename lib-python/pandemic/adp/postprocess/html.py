import os
from libtbx import adopt_init_args
from pandemic.adp.html import HtmlSummary, divs


class PostProcessingHtmlSummary(HtmlSummary):


    r_factor_keys = ['r_free', 'r_work', 'r_gap']

    def __init__(self,
        results_object,
        analysis_task,
        ):
        _json_plots = []
        adopt_init_args(self, locals())

    def json_plots(self):
        return self._json_plots

    def main_summary(self):

        tab = divs.Tab(
            id = 'statistics',
            title = 'Model R-factor Summary',
            alt_title = 'R-factor Statistics',
        )

        tab.extend(self.get_interactive_results_plot_panels())
        tab.extend(self.get_dataset_statistics_table_panels())
        tab.extend(self.get_rfactor_panels())

        if len(tab.contents) == 0:
            return []

        return [tab]

    def get_rfactor_panels(self):

        output = []

        analysis_files = self.analysis_task.result.output_files

        f = analysis_files.get('r_values')
        if f is not None:
            panel = divs.Panel(
                title = 'R-factor comparisons across structures',
                contents = [
                    divs.Block(
                        width = 6,
                        title = 'R-Free Comparison',
                        image = self.image(f.get('r_free')),
                        footnote = "Comparison of R-free values for the different input and output structures (and optionally reference R-values).",
                    ),
                    divs.Block(
                        width = 6,
                        title = 'R-Work Comparison',
                        image = self.image(f.get('r_work')),
                        footnote = "Comparison of R-work values for the different input and output structures (and optionally reference R-values).",
                    ),
                    divs.Block(
                        width = 6,
                        title = 'R-Gap Comparison',
                        image = self.image(f.get('r_gap')),
                        footnote = "Comparison of R-work values for the different input and output structures (and optionally reference R-values).",
                    ),
                ],
            )
            output.append(panel)

        f = analysis_files.get('r_value_differences')
        if f is not None:
            panel = divs.Panel(
                title = 'R-factor differences between structures',
            )
            output.append(panel)
            for (s1, s2), img in f.items():
                panel.append(
                    divs.Block(
                        width = 6,
                        title = 'R-factor differences between {} and {} structures'.format(
                            s1.lower(),
                            s2.lower(),
                        ),
                        image = self.image(img),
                        footnote = "Changes in the R-factors from the {} model to the {} model (negative values mean lower values for the {} model).".format(
                            s2.lower(),
                            s1.lower(),
                            s1.lower(),
                        ),
                    )
                )

        # Only add title if some data added
        if len(output) > 0:
            output.insert(0, divs.Block(title='R-factor statistics'))

        return output

    def get_interactive_results_plot_panels(self):

        output = []

        # Extract data for plotting
        table = self.results_object.table.dropna(axis='columns', how='all')
        if len(table.columns) > 0:
            # JSON data for the interactive plots
            json_plot = {
                'id'        : 'variable-plots-{}'.format(self.counter.next()),
                'json'      : table.T.to_json(orient='split'),
                'default_x' : 'High Resolution Limit',
                'default_y' : 'R-free Diff. (Fitted-Input)',
                }
            self._json_plots.append(json_plot)
            # Panel for the interactive plots to be loaded into
            output.append(
                divs.Panel(
                    title = 'Interactive Summary Graphs',
                    contents = [divs.Div(id=json_plot['id'])],
                )
            )

        return output

    def get_dataset_statistics_table_panels(self):

        output = []

        # Extract table columns with data
        table = self.results_object.table.dropna(axis='columns', how='all')
        if len(table.columns) > 0:
            output.append(
                divs.Panel(
                    title = 'Individual Dataset Statistics Table (click to expand/collapse)',
                    show = False,
                    contents = [
                        divs.Block(
                            text = 'Data from output CSV',
                            table = table.to_html(
                                bold_rows = False,
                                classes = ['table table-striped table-hover datatable nowrap'],
                            ) \
                            .replace('<th></th>','<th>Dataset</th>') \
                            .replace('border="1" ', ''),
                        ),
                    ],
                )
            )

        return output

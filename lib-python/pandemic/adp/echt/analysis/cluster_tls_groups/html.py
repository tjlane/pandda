import os
import pandas
from libtbx import adopt_init_args
from pandemic.adp.html import HtmlSummary, divs


class ClusterTLSGroupsTaskHtmlSummary(HtmlSummary):


    def __init__(self, task):
        adopt_init_args(self, locals())

    def main_summary(self):

        output = divs.Tab(alt_title='TLS Group Clustering')

        txt = """
        > TLS Group Clustering
        Clustering of the TLS groups in each level to identify possible collective motions composed of multiple groups.
        """
        output.extend(self.format_summary(txt, classes=["square-corners-top"]))
        output.extend(self.make_intro_panel())
        output.extend(self.make_tab_set())

        return [output]

    def make_intro_panel(self):

        return self.format_summary(
            self.task.cluster.description() \
                .replace("->","&#8594;") \
                .replace("Angstrom", "&#8491;"),
            width = 6,
            )

    def make_tab_set(self):

        tab_set = divs.TabSet(title='Clustering Results')

        for level_name, clustering_result in self.task.result.clustering_results.items():
            clustering_files = self.task.result.output_files[level_name]

            cluster_summary = divs.Tab(
                title = 'Clustering results for {} Level'.format(level_name.title()),
                alt_title = level_name.title(),
                contents = self.cluster_result_html_summary_contents(
                    cluster_result = clustering_result,
                    cluster_files = clustering_files,
                ),
            )
            tab_set.append(cluster_summary)
        tab_set.set_active()

        return [tab_set]

    def cluster_result_html_summary_contents(self, cluster_result, cluster_files):

        clustering_function = self.task.cluster

        contents = []

        txt = """
        > New groupings for different B-factor thresholds
        At each threshold, modules are created from TLS groups that can reproduce the B-factors of other TLS groups above that threshold.
        If one module contains another module, then only the larger module is shown. If two modules overlap, but not completely, both are shown.

        > Note
        This analysis does not account for the shape of the groups and so will likely generate groups that do not represent physical motions.
        The similarity between two groups does not mean that they move together, only that they have similar/overlapping disorder patterns.
        Manual inspection and filtering of the groups is crucial.
        """
        contents.append(divs.Block(width=4, contents=self.format_summary(txt)))

        pymol_txt = """
        Command to view in PyMOL:
        <pre>pymol {script}</pre>
        The colours in the graph below are the same as the colours in the pymol session.
        """.format(
            script=cluster_files.get('pymol_script'),
        )

        pymol_block = divs.Block(width=8)
        pymol_block.extend(self.format_summary(pymol_txt))
        pymol_block.append(divs.Block(image=self.image(cluster_files.get('modules_png'))))
        contents.append(pymol_block)

        txt = """
        > Parameter files for new groupings
        Input parameters (*.eff files) for these groups can be found in: {out_dir}
        """.format(
            out_dir=os.path.dirname(cluster_files.get('modules_png','')),
        )
        contents.extend(self.format_summary(txt))

        comparison_table = pandas.DataFrame(cluster_result.comparison_matrix.round(1))
        comparison_table.index = ['Group {}'.format(i) for i in range(1,len(comparison_table)+1)]
        comparison_table.columns = comparison_table.index.values
        panel = divs.Panel(
            title = 'Similarity Matrix',
            contents = [
                divs.Block(
                    text = '{} between TLS Groups'.format(clustering_function.metric_type.title()),
                    table = comparison_table.to_html(
                        bold_rows = False,
                        justify = 'center',
                        classes = ['table table-hover datatable nowrap text-center']) \
                            .replace('<th>', '<th class="text-center">') \
                            .replace('border="1" ', ''),
                ),
            ],
        )
        contents.append(panel)

        connectivity_table = pandas.DataFrame(cluster_result.connectivity.astype(int))
        connectivity_table.index = ['Group {}'.format(i) for i in range(1,len(connectivity_table)+1)]
        connectivity_table.columns = connectivity_table.index.values
        panel = divs.Panel(
            title = 'Connectivity Matrix',
            show = False,
            contents = [
                divs.Block(
                    text = 'Connectivity between TLS Groups',
                    table = connectivity_table.to_html(
                        bold_rows = False,
                        justify = 'center',
                        classes = ['table table-hover datatable nowrap text-center']) \
                            .replace('<th>', '<th class="text-center">') \
                            .replace('border="1" ', ''),
                ),
            ],
        )
        contents.append(panel)

        return contents

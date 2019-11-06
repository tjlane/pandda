import pandas
from libtbx import adopt_init_args
from pandemic.adp.html import HtmlSummary


class ClusterTLSGroupsTaskHtmlSummary(HtmlSummary):


    def __init__(self, task):
        adopt_init_args(self, locals())

    def main_summary(self):

        output = {
            'title' : 'TLS Group Clustering',
            'alt_title' : 'Clustering',
            'contents' : [],
        }

        output['contents'].extend(self.make_intro_panel())
        output['contents'].extend(self.make_tab_set())

        return [output]

    def make_intro_panel(self):

        return self.format_summary(self.task.cluster.description())

    def make_tab_set(self):

        tabs = []
        for level_name, clustering_result in self.task.result.clustering_results.iteritems():
            clustering_files = self.task.result.output_files[level_name]
            cluster_summary = self.cluster_result_html_summary(
                cluster_result = clustering_result,
                cluster_files = clustering_files,
                )
            cluster_summary['title'] = 'Clustering results for Level "{}"'.format(level_name)
            cluster_summary['alt_title'] = ' '.join([s.capitalize() for s in level_name.split(' ')])
            tabs.append(cluster_summary)

        # Make first active
        tabs[0]['active'] = True

        return [{
            'type' : 'tabs',
            'title' : 'Clustering Results',
            'contents' : tabs,
        }]

    def cluster_result_html_summary(self, cluster_result, cluster_files):

        clustering_function = self.task.cluster

        out = {'contents'  : []}

        out['contents'] += [{
            'title' : 'New groupings for different thresholds',
            'width' : 8,
            'image' : self.image(cluster_files.get('modules_png')),
            'footnote' : 'to visualise, run: <pre>pymol {script}</pre>'.format(script=cluster_files.get('pymol_script')),
        }]

        out['contents'] += [{
            'type' : 'alert',
            'colour' : 'info',
            'width' : 6,
            'title' : 'Connectivity Matrix',
            'text' : 'Connectivity between TLS Groups',
            'table' : pandas.DataFrame(cluster_result.connectivity.astype(int)).to_html(
                bold_rows=False,
                classes=['table table-striped table-hover datatable nowrap'])\
                               .replace('border="1" ', ''),
        }]

        out['contents'] += [{
            'type' : 'alert',
            'colour' : 'info',
            'width' : 6,
            'title' : 'Similarity Matrix',
            'text' : '{} between TLS Groups'.format(clustering_function.metric_type.capitalize()),
            'table' : pandas.DataFrame(cluster_result.comparison_matrix.round(1)).to_html(
                bold_rows=False,
                classes=['table table-striped table-hover datatable nowrap'])\
                               .replace('border="1" ', ''),
        }]

        return out

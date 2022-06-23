import giant.logs as lg
logger = lg.getLogger(__name__)

import copy, collections

from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure

from pandemic.adp.echt.analysis.cluster_tls_groups import ClusterTLSGroupsTask
from pandemic.adp.echt.analysis.amplitudes import AnalyseTLSAmplitudesTask

import pathlib as pl


class AnalyseEchtModelTask(object):


    def __init__(self,
        output_directory,
        master_phil,
        analysis_parameters = None,
        ):

        output_directory = pl.Path(output_directory)

        cluster_tls_groups_task = ClusterTLSGroupsTask(
            output_directory = (output_directory / 'tls_group_clustering'),
            metric = analysis_parameters.tls_group_clustering.metric,
            parameters = analysis_parameters.tls_group_clustering,
            )

        analyse_tls_amplitudes_task = AnalyseTLSAmplitudesTask(
            output_directory = (output_directory / 'dataset_clustering'),
            )

        adopt_init_args(self, locals())

    def run(self,
        model_object,
        model_summary_output,
        ):

        model_files = model_summary_output.output_files

        # Add function to clustering to allow writing of levels
        # (requires model object which is not available at initialisation)
        from pandemic.adp.hierarchy.custom_levels import MakeNewCustomLevelEffFilesFromIndices
        self.cluster_tls_groups_task.write_levels_function = MakeNewCustomLevelEffFilesFromIndices.from_model_object(
            model_object = model_object,
            master_phil = self.master_phil,
            custom_level_scope_name = 'model.custom_level',
            )

        output_files = collections.OrderedDict()

        ######################################

        clustering_results = self.cluster_tls_groups_task.run(
            model_object = model_object,
            model_files = model_files,
            )
        output_files['clustering'] = clustering_results.output_files

        ######################################

        amplitude_analysis = self.analyse_tls_amplitudes_task.run(
            )

        ######################################

        self.result = group_args(
            clustering_results = clustering_results,
            output_files = output_files,
            )

        return self.result

    def as_html_summary(self):

        from pandemic.adp.html import (
            HtmlSummaryCollator,
            )

        return HtmlSummaryCollator(
            title = 'Analysis of fitted ECHT model',
            alt_title = 'ECHT Analysis',
            summaries = [
                self.cluster_tls_groups_task.as_html_summary(),
                #self.analyse_tls_amplitudes_task.as_html_summary(),
                ],
            )


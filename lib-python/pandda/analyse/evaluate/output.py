import giant.logs as lg
logger = lg.getLogger(__name__)

from pandda.utils import (
    merge_dicts,
    )


class PanddaDatasetEvaluatorOutputter(object):
    """Output dataset maps and graphs (inside processing loop)"""

    def __init__(self,
        functions,
        ):

        self.functions = (
            functions
            )

    def __call__(self,
        datasets,
        datasets_map_dict,
        datasets_results,
        statistical_model,
        map_resolution,
        ):

        output_files = {}

        for func in self.functions:

            of = func(
                datasets = datasets,
                datasets_map_dict = datasets_map_dict,
                datasets_results = datasets_results,
                statistical_model = statistical_model,
                map_resolution = map_resolution,
                )

            merge_dicts(
                master_dict = output_files,
                merge_dict = of,
                )

        return output_files


class GetPanddaDatasetEvaluatorOutputter(object):

    def __init__(self,
        get_dataset_map_writer,
        output_dir,
        dataset_dir,
        processor = None,
        output_requires_events = True,
        ):

        from .mtzs import (
            GetMakePanddaEvaluationMtzs,
            )

        from .graphs import (
            GetMakePanddaEvaluationGraphs,
            )

        self.getters = [
            GetMakePanddaEvaluationMtzs(
                get_dataset_map_writer = get_dataset_map_writer,
                output_dir = output_dir,
                dataset_dir = dataset_dir, # only maps put in the dataset folders
                processor = processor,
                output_requires_events = output_requires_events,
                ),
            GetMakePanddaEvaluationGraphs(
                output_dir = output_dir,
                dataset_dir = None, # forces redirect to output dir
                processor = processor,
                ),
        ]

    def __call__(self, label):

        return PanddaDatasetEvaluatorOutputter(
            functions = [
                g(label=label)
                for g in self.getters
                ],
            )


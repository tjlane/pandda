import giant.logs as lg
logger = lg.getLogger(__name__)

from pandda.utils import (
    merge_dicts,
    )

from .maps import (
    GetMakePanddaEvaluationMaps,
    )

from .graphs import (
    GetMakePanddaEvaluationGraphs,
    )


class PanddaDatasetEvaluatorOutputter:
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
        ):

        output_files = {}

        for func in self.functions:

            of = func(
                datasets = datasets,
                datasets_map_dict = datasets_map_dict,
                datasets_results = datasets_results,
                statistical_model = statistical_model,
                )

            merge_dicts(
                master_dict = output_files,
                merge_dict = of,
                )

        return output_files


class GetPanddaDatasetEvaluatorOutputter:

    def __init__(self,
        get_dataset_map_writer,
        output_dir,
        dataset_dir,
        processor = None,
        ):

        self.getters = [
            GetMakePanddaEvaluationMaps(
                get_dataset_map_writer = get_dataset_map_writer,
                output_dir = output_dir,
                dataset_dir = dataset_dir, # only maps put in the dataset folders
                processor = processor,
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


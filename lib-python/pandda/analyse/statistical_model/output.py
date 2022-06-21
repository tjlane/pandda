import giant.logs as lg
logger = lg.getLogger(__name__)

from pandda.utils import (
    merge_dicts,
    )

from .graphs import (
    GetMakePanddaStatisticalModelGraphs,
    )


class PanddaStatisticalModelOutputter(object):
    """Output statistical model maps and graphs (inside processing loop)"""

    def __init__(self,
        functions,
        ):

        self.functions = (
            functions
            )

    def __call__(self,
        statistical_model,
        **kwargs
        ):

        output_files = {}

        for func in self.functions:

            of = func(
                statistical_model = statistical_model,
                **kwargs
                )

            merge_dicts(
                master_dict = output_files,
                merge_dict = of,
                )

        return output_files


class GetPanddaStatisticalModelOutputter(object):

    def __init__(self,
        output_dir,
        processor = None,
        ):

        self.getters = [
            GetMakePanddaStatisticalModelGraphs(
                output_dir = output_dir,
                processor = processor,
                ),
        ]

    def __call__(self, label):

        return PanddaStatisticalModelOutputter(
            functions = [
                g(label=label)
                for g in self.getters
                ],
            )



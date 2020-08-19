import giant.logs as lg
logger = lg.getLogger(__name__)

from giant.mulch.statistics import (
    ExtractAndClassifyDatasetStatistics,
    )

from giant.mulch.statistics.xray import (
    ExtractBasicXrayStatistics,
    ExtractWilsonStatistics, 
    ClassifyWilsonStatistics,
    )

from giant.mulch.statistics.scaling import (
    ExtractScalingStatistics, 
    ClassifyScalingStatistics,
    )


class PanddaExtractDatasetStatistics(ExtractAndClassifyDatasetStatistics):

    def __init__(self,
        max_scaling_z_score,
        max_wilson_z_score,
        ):

        extracters = [
            ExtractBasicXrayStatistics(),
            ExtractWilsonStatistics(),
            ExtractScalingStatistics(),
        ]

        classifiers = [
            ClassifyWilsonStatistics(
                outlier_partition = "not_train",
                max_z_score = max_wilson_z_score,
                ),
            ClassifyScalingStatistics(
                outlier_partition = "not_train",
                max_z_score = max_scaling_z_score,
                ),
        ]

        super(PanddaExtractDatasetStatistics, self).__init__(
            extracters = extracters,
            classifiers = classifiers,
            )

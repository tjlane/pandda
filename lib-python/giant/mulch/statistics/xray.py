import collections

from giant.mulch.statistics.classifiers import (
    ZScoreClassifier,
    )


class ExtractBasicXrayStatistics(object):

    def __init__(self):
        pass

    def __call__(self, dataset, **kw_args):

        d_dict = self.extract(dataset)

        return d_dict

    def extract(self, dataset):

        d = collections.OrderedDict([
            (
                "r_free",
                dataset.model.crystal.r_free,
                ),
            (
                "r_work",
                dataset.model.crystal.r_work,
                ),
            (
                "high_resolution",
                dataset.model.crystal.resolution_high,
                ),
            (
                "low_resolution",
                dataset.model.crystal.resolution_low,
                ),
            (
                "r_free/r_work",
                (dataset.model.crystal.r_free / dataset.model.crystal.r_work),
                ),
            (
                "space_group",
                str(dataset.model.crystal.space_group.info()),
                ),
            (
                "unit_cell",
                dataset.model.crystal.unit_cell.parameters(),
                ),
            ])

        return d


class ExtractWilsonStatistics(object):

    def __init__(self):
        pass

    def __call__(self, 
        dataset, 
        get_miller_array, 
        get_scaling_object, 
        **kw_args
        ):

        # Default return None
        wilson_b_input = wilson_b_scaled = None

        # Get input miller array
        miller_array = get_miller_array(dataset)
        # Get wilson B of input array
        wilson_b_input = self.get_wilson_b_factor(miller_array)

        # Get scaling function
        scale_miller_array = get_scaling_object(dataset)

        # Scale if provided
        if (scale_miller_array is not None):
            # Get scaled array
            scaled_miller_array = scale_miller_array(miller_array)
            # Calculate scaled wilson B
            wilson_b_scaled = self.get_wilson_b_factor(scaled_miller_array)

        return collections.OrderedDict([
            ('wilson_b_input', wilson_b_input),
            ('wilson_b_scaled', wilson_b_scaled),
            ])

    def get_wilson_b_factor(self, miller_array):
        from giant.xray.data import estimate_wilson_b_factor
        return estimate_wilson_b_factor(miller_array=miller_array)


class ClassifyWilsonStatistics(ZScoreClassifier):
    
    default_columns = ["wilson_b_scaled"]

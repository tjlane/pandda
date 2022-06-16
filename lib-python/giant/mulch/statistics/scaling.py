import collections
import numpy as np

from .classifiers import ZScoreClassifier


class ExtractScalingStatistics(object):

    resolution_divider = 4.0

    def __init__(self):
        pass

    def __call__(self, dataset, get_scaling_object, **kw_args):

        scaling = get_scaling_object(dataset)

        if (not scaling):
            return {}

        s_info = self.extract_scaling_statistics(
            scaling_result = scaling.scaling_result,
            )

        return s_info

    def extract_scaling_statistics(self, scaling_result):

        scaling = scaling_result

        # Selections for high resolution and low resolution
        high_res_sel = scaling.x_values > 1.0 / (self.resolution_divider ** 2)
        low_res_sel = scaling.x_values <= 1.0 / (self.resolution_divider ** 2)

        return collections.OrderedDict([
            ##########
            (
                "b_factor_scaling",
                np.round(scaling.info["scaling_b_factor"], 3),
                ),
            ##########
            # (
            #     "scaled_wilson_rmsd_all",
            #     np.round(scaling.info["scaled_rmsd"], 3),
            #     ),
            (
                "scaled_wilson_rmsd",
                np.round(scaling.rmsd_to_ref(values=scaling.out_values), 3),
                ),
            (
                "scaled_wilson_rmsd_high",
                np.round(scaling.rmsd_to_ref(values=scaling.out_values, sel=high_res_sel), 3),
                ),
            (
                "scaled_wilson_rmsd_low",
                np.round(scaling.rmsd_to_ref(values=scaling.out_values, sel=low_res_sel), 3),
                ),
            (
                "scaled_wilson_ln_rmsd",
                np.round(scaling.info["scaled_ln_rmsd"], 3),
                ),
            (
                "scaled_wilson_ln_dev",
                np.round(scaling.info["scaled_ln_dev"], 3),
                ),
            ##########
            # (
            #     "unscaled_wilson_rmsd_all",
            #     np.round(scaling.unscaled_rmsd, 3),
            #     ),
            (
                "unscaled_wilson_rmsd",
                np.round(scaling.rmsd_to_ref(values=scaling.scl_values), 3),
                ),
            (
                "unscaled_wilson_rmsd_high",
                np.round(scaling.rmsd_to_ref(values=scaling.scl_values, sel=high_res_sel), 3),
                ),
            (
                "unscaled_wilson_rmsd_low",
                np.round(scaling.rmsd_to_ref(values=scaling.scl_values, sel=low_res_sel), 3),
                ),
            (
                "unscaled_wilson_ln_rmsd",
                np.round(scaling.info["unscaled_ln_rmsd"], 3),
                ),
            (
                "unscaled_wilson_ln_dev",
                np.round(scaling.info["unscaled_ln_dev"], 3),
                ),
            ##########
            ])


class ClassifyScalingStatistics(ZScoreClassifier):

    default_columns = [
        "scaled_wilson_rmsd",
        "scaled_wilson_rmsd_high",
        "scaled_wilson_rmsd_low",
        "scaled_wilson_ln_rmsd",
        "scaled_wilson_ln_dev",
        ]

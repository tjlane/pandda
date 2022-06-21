
from pandda.utils import (
    merge_dicts,
    make_sorted_dict,
    )

from giant.mulch.statistics import (
    ExtractDatasetStatistics,
    )

from giant.mulch.statistics.xray import (
    ExtractBasicXrayStatistics,
    )


class ExtractPanddaDatasetOutputInfo(object):

    def __init__(self):

        self.extract = ExtractDatasetStatistics(
            extracters = [
                ExtractBasicXrayStatistics(),
                ],
            )

    def __call__(self, 
        mcd, 
        shell_dicts,
        dataset_dicts,
        ):

        # Input data from extract is
        #    { dataset_key : { variable_key : value } }
        # Want dicts of 
        #    { variable_key : { dataset_key : value } }
        # pandas dataframe is return the desired way
        dataset_statistics = self.extract(
            mcd = mcd,
            )
        
        merge_dicts(
            master_dict = dataset_dicts,
            merge_dict = make_sorted_dict(
                dataset_statistics.as_dataframe().to_dict()
                ),
            )

        ###

        merge_dicts(
            master_dict = dataset_dicts,
            merge_dict = make_sorted_dict(
                self.unpack_shell_dicts(
                    shell_dicts = shell_dicts,
                    all_dataset_keys = list(mcd.datasets.keys()),
                    )
                ),
            )

        return dataset_dicts

    def unpack_shell_dicts(self, 
        shell_dicts,
        all_dataset_keys,
        ):

        shell_info = {
            'train' : {k:False for k in all_dataset_keys},
            'test' : {k:False for k in all_dataset_keys},
            }

        for s_dict in shell_dicts:

            for k in s_dict['train']:
                shell_info['train'][k] = True

            for k in s_dict['test']:
                shell_info['test'][k] = True

        return shell_info
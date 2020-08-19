import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy as np


class DatasetKeySelector:

    def __init__(self, dataset_keys):

        self.dataset_keys = dataset_keys

    def __call__(self, datasets):

        overlap_keys = set(self.dataset_keys).intersection(datasets.keys())

        selected = {}

        for d_key in overlap_keys:

            selected[d_key] = datasets[d_key]

        return selected


class SortedDatasetSelector:

    def __init__(self,
        max_datasets,
        sort_datasets_func,
        ):

        self.max_datasets = max_datasets
        self.sort_datasets_func = sort_datasets_func

    def __call__(self, datasets):

        if (self.max_datasets is None):
            return datasets

        sorted_dtag_dataset_tuples = sorted(
            list(datasets.items()),
            key = lambda d_tuple: self.sort_datasets_func(d_tuple[1]),
            )

        sorted_datasets = dict(
            sorted_dtag_dataset_tuples[:self.max_datasets]
            )

        return sorted_datasets


class RandomDatasetSelector: 

    def __init__(self,
        max_datasets,
        seed = 616,
        replace = False,
        ):

        self.max_datasets = max_datasets

        self.generator = np.random.RandomState(seed)
        self.replace = False

    def __call__(self, datasets):

        if (self.max_datasets is None):
            return datasets

        if (self.replace is False) and (len(datasets) <= self.max_datasets): 
            return datasets

        all_keys = sorted(datasets.keys())

        sample_keys = self.generator.choice(
            all_keys, 
            size = self.max_datasets,
            replace = self.replace,
            )

        sample_datasets = {
            dkey : datasets[dkey]
            for dkey in sample_keys
        }

        return sample_datasets


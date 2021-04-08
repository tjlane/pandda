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
        dataset_sorter,
        ):

        self.max_datasets = max_datasets

        self.dataset_sorter = dataset_sorter

    def __call__(self, datasets):

        if (self.max_datasets is None):
            return datasets

        sorted_dtag_dataset_tuples = self.dataset_sorter(datasets)

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


class SoulmateSelector:

    def __init__(self,
        dataset_sorter,
        ):

        self.dataset_sorter = dataset_sorter

        self.chosen_dkey = None

    def __call__(self, datasets):

        if self.chosen_dkey is None:

            sorted_dtag_dataset_tuples = self.dataset_sorter(datasets)

            self.chosen_dkey = sorted_dtag_dataset_tuples[0][0]

        chosen_dataset = datasets.get(self.chosen_dkey, None)

        if chosen_dataset is None:
            return {}

        return {self.chosen_dkey : chosen_dataset}

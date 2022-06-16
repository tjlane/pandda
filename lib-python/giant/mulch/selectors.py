import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy as np


class DatasetKeySelector(object):

    def __init__(self, dataset_keys):

        self.dataset_keys = dataset_keys

    def __call__(self, datasets):

        overlap_keys = set(self.dataset_keys).intersection(list(datasets.keys()))

        selected = {}

        for d_key in overlap_keys:

            selected[d_key] = datasets[d_key]

        return selected


class SortedDatasetSelector(object):


    def __init__(self,
        max_datasets,
        sort_datasets,
        ):
        """
        max_datasets: int
        sort_dataset_items:
            function that takes dictionary of (key, dataset) pairs and returns sorted list of (key, dataset) tuples
        """

        self.max_datasets = max_datasets

        self.sort_datasets = sort_datasets

    def __call__(self, datasets):
        """
        datasets: dictionary of key-dataset pairs
        """

        if (self.max_datasets is None):
            return datasets

        sorted_dataset_items = self.sort_datasets(datasets)

        selected_datasets = dict(
            sorted_dataset_items[:self.max_datasets]
            )

        return selected_datasets


class RandomDatasetSelector(object):

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


class SoulmateSelector(object):

    def __init__(self,
        sort_datasets,
        ):
        """
        sort_datasets: DatasetSorter-like class that returns a sorted list of (key, dataset) tuples
        """

        self.sort_datasets = sort_datasets

        self.chosen_dkey = None

    def __call__(self, datasets):

        if self.chosen_dkey is None:

            sorted_dataset_items = self.sort_datasets(datasets)

            self.chosen_dkey = sorted_dataset_items[0][0]

        chosen_dataset = datasets.get(self.chosen_dkey, None)

        if chosen_dataset is None:
            return {}

        return {self.chosen_dkey : chosen_dataset}

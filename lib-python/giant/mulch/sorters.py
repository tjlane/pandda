import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy as np


class DatasetSorter(object):

    def __init__(self, reverse=False):

        self.reverse = reverse

    def __call__(self, datasets):

        if hasattr(datasets, 'items'):
            sorted_datasets = sorted(
                list(datasets.items()),
                key = lambda d: self.sort_function(d[1]),
                reverse = self.reverse,
                )
        else:
            sorted_datasets = sorted(
                list(datasets),
                key = lambda d: self.sort_function(d),
                reverse = self.reverse,
                )

        return sorted_datasets
    
    def sort_function(self, dataset):

        raise NotImplemented()


class HighResolutionSorter(DatasetSorter):

    def sort_function(self, dataset):

        return dataset.data.crystal.resolution_high


class ArbitraryDatasetSorter(DatasetSorter):

    def __init__(self, sort_function, reverse=False):

        self._sort_function = sort_function
        self.reverse = reverse

    def sort_function(self, dataset):
        return self._sort_function(dataset)

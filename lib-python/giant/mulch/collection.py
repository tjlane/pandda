import giant.logs as lg
logger = lg.getLogger(__name__)

import copy 

from giant.mulch.dataset import CrystallographicDataset


class MultiCrystalDataset(object):

    name = "MultiCrystalDataset"

    def __init__(self, datasets=None, sort_keys=True):

        if (datasets is None):
            datasets = {}

        self.datasets = datasets
        self.dataset_keys = list(datasets.keys()) # can be reordered as desired

        if (sort_keys is True):
            self.dataset_keys.sort()

        self.populate()

    def __iter__(self):
        # TODO make this also yield dataset tag?
        for dtag in self.dataset_keys:
            yield self.datasets[dtag]

    def __str__(self):

        datasets_strs = [
            'Dataset {dkey}:\n\t{dataset}'.format(
                dkey = dkey, 
                dataset = (
                    str(self.datasets[dkey])
                    ).strip('\n').replace('\n','\n\t'),
                )
            for dkey in self.dataset_keys
        ]

        s_ = (
            'Collection Type: {name}\n'
            '| Number of Datasets: {n_datasets}\n'
            '| Datasets:\n'
            '|\t{datasets}\n'
            '`---->'
            ).format(
            name = self.name,
            n_datasets = self.n_datasets(),
            datasets = (
                '\n'.join(datasets_strs)
                ).strip('\n').replace('\n','\n|\t'),
            )

        return s_

    def populate(self):

        # Move to slots
        self.resolution_high = None
        self.resolution_low = None

        self.reference_dataset = None

        if (self.datasets is None) or (len(self.datasets) == 0):
            return self

        self.resolution_high = min([
            d.model.crystal.resolution_high
            for dtag, d
            in self.datasets.items()
        ])

        self.resolution_low = max([
            d.model.crystal.resolution_high
            for dtag, d
            in self.datasets.items()
        ])

        return self

    def n_datasets(self):
        return len(self.datasets)

    def new_from_datasets(self, datasets):
        """Create a new collections of datasets inheriting class from the parent MultiCrystalDataset"""

        clone = MultiCrystalDataset(
            datasets = datasets,
        )

        clone.set_reference_dataset(
            copy.deepcopy(self.reference_dataset)
            )

        return clone

    def partition(self, partitioner):
        """
        Partition the dataset into different multi-crystal datasets.

        "partitioner" must return a dict of partitions, e.g.
        {
            'partition1' : dictofdatasets,
            'partition2' : dictofdatasets,
        }
        """

        partitioned = partitioner(self.datasets)

        assert hasattr(partitioned, 'keys'), "partitioners must return dict-like object"

        result = {
            p_key : self.new_from_datasets(partitioned[p_key])
            for p_key in partitioned.keys()
        }

        return result

    def set_reference_dataset(self, dataset):

        self.reference_dataset = dataset

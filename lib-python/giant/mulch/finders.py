import giant.logs as lg
logger = lg.getLogger(__name__)


class DatasetFinder(object):
    """Returns the key to the dataset identified by a supplied function"""

    def __init__(self):

        self.name = None

    def __call__(self, mcd):
        
        dataset_keys, datasets = list(zip(*list(mcd.datasets.items())))

        # Get the rankings of the datasets
        dataset_ranks = self.get_ranks(datasets)

        # Get the index of the highest rankes dataset (lowest value)
        i_dst = dataset_ranks.index(min(dataset_ranks))

        return dataset_keys[i_dst]

    def get_ranks(self, datasets):
        raise Exception('not implemented in base class')


class DummyDatasetFinder(DatasetFinder):

    def __init__(self):

        super(DummyDatasetFinder, self).__init__()

        self.name = "DummyDatasetFinder"

    def get_ranks(self, datasets):
        return list(range(len(datasets)))


class HighestResolutionFinder(DatasetFinder):

    def __init__(self):

        super(HighestResolutionFinder, self).__init__()

        self.name = "HighestResolutionFinder"

    @classmethod
    def get_resolution(cls, dataset):
        return dataset.data.mtz_object().max_min_resolution()[1]

    def get_ranks(self, datasets):
        return [self.get_resolution(d) for d in datasets]


class LowestRfreeFinder(DatasetFinder):

    def __init__(self):

        super(LowestRfreeFinder, self).__init__()

        self.name = "LowestRfreeFinder"

    @classmethod
    def get_rfree(cls, dataset):
        return dataset.model.input.get_r_rfree_sigma().r_free

    def get_ranks(self, datasets):
        return [self.get_rfree(d) for d in datasets]

import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy
from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure


class SelectOptimisationDatasetsTask(object):


    def __init__(self,
            max_resolution,
            max_datasets,
            sort_datasets_by = 'resolution',
            random_seed = 0,
            ):
        adopt_init_args(self, locals())

    def run(self,
            models,
            resolutions = None,
            ):

        if resolutions is not None:
            assert len(models) == len(resolutions)

        if (self.max_resolution is None):
            opt_datasets = list(range(len(models)))
        else:
            assert resolutions is not None, 'no resolutions provided'
            opt_datasets = numpy.where(resolutions < self.max_resolution)[0]

        if len(opt_datasets) == 0:
            raise Sorry('No datasets selected for optimisation (e.g. above resolution cutoff: {})'.format(self.max_resolution))

        if self.sort_datasets_by == 'resolution':
            assert resolutions is not None, 'no resolutions provided'
            opt_datasets = sorted(opt_datasets, key=lambda i: resolutions[i])
        elif self.sort_datasets_by == 'name':
            opt_datasets = sorted(opt_datasets, key=lambda i: self.models[i].tag)
        elif self.sort_datasets_by == 'random':
            logger.debug('Setting random seed: {}'.format(self.random_seed))
            numpy.random.seed(self.random_seed)
            opt_datasets = numpy.random.permutation(opt_datasets)

        logger('After reordering:')
        for i_m in opt_datasets:
            logger('\t{}: {}'.format(i_m, self.models[i_m].tag))

        # Limit the number of datasets for optimisation
        if (self.max_datasets is not None) and (len(opt_datasets) > self.max_datasets):
            logger('\nLimiting list of datasets for TLS optimisation to {} datasets'.format(self.max_datasets))
            opt_datasets = opt_datasets[:self.max_datasets]

        assert len(selection) > 0, 'no datasets selected for optimisation with resolution cutoff: {}'.format(self.max_resolution)

        self.iselection = opt_datasets
        return iselection



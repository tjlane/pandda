import giant.logs as lg
logger = lg.getLogger(__name__)


class DatasetFilter(object):

    name = "DatasetFilter"

    def __init__(self):

        self.rejections = {}

    def __call__(self, mcd, **kw_args):

        remove_datasets = []

        # Filter datasets
        for dtag, dataset in mcd.datasets.items():

            valid, reason = self.is_valid(dataset, **kw_args)

            if (valid is False):
                self.rejections[dtag] = reason
                remove_datasets.append(dtag)

        # Create a new dataset
        new_datasets = {
            dtag: d
            for dtag, d in mcd.datasets.items()
            if (dtag not in remove_datasets)
        }
        new_dataset = mcd.new_from_datasets(
            datasets = new_datasets,
        )

        return new_dataset

    def __str__(self):

        parameters = self.parameters()
        rejections = self.rejections 

        parameter_strings = [
            '{k}: {v}'.format(
                k=k, 
                v=str(v).strip().replace('\n','\n\t'),
                )
            for k,v in parameters.items()
            ]

        rejection_strings = [
            '{k}: {v}'.format(
                k=k, 
                v=str(v).strip().replace('\n','\n\t'),
                )
            for k,v in sorted(rejections.items())
        ]

        s_ = (
            'Filter Type: {name}\n'
            '| Parameters: \n'
            '|\t{parameters}\n'
            '| Rejected Datasets: \n'
            '|\t{rejections}\n'
            '`---->'
            ).format(
            name = self.name,
            parameters = (
                '\n'.join(parameter_strings).strip().replace('\n','\n|\t')
                if parameter_strings else "None"
                ),
            rejections = (
                '\n'.join(rejection_strings).strip().replace('\n','\n|\t')
                if rejection_strings else "None"
                ),
            )

        return s_.strip()

    def parameters(self):
        return {}

    def is_valid(self, dataset):
        raise Exception('not implemented')


class ManualDatasetFilter(DatasetFilter):

    name = "ManualDatasetFilter"

    def __init__(self, rejections):
        """Filter datasets by a predetermined set of rejections"""

        super(ManualDatasetFilter, self).__init__()

        if not hasattr(rejections, "keys"): 
            rejections = dict.fromkeys(rejections, None)

        self.rejections.update(rejections)

    def __call__(self, mcd):

        remove_datasets = list(self.rejections.keys())

        # Create a new dataset
        new_datasets = {
            dtag: d
            for dtag, d in mcd.datasets.items()
            if (dtag not in remove_datasets)
        }
        new_dataset = mcd.new_from_datasets(
            datasets = new_datasets,
        )

        return new_dataset


class InclusiveDatasetTagFilter(DatasetFilter):

    name = "InclusiveDatasetTagFilter"

    def __init__(self, dataset_tags):

        super(InclusiveDatasetTagFilter, self).__init__()

        self.dataset_tags = dataset_tags

    def parameters(self):

        p_ = {
            'dataset_tags' : self.dataset_tags,
        }

        return p_

    def get_dataset_tag(self, dataset):

        if dataset.tag is None: 
            raise ValueError(
                "This function can only be used with labelled datasets with a 'tag' attribute. "
                "Datasets can be labelled using the label method."
                )

        return dataset.tag

    def is_valid(self, dataset, **kwargs):

        if self.get_dataset_tag(dataset) in self.dataset_tags:
            return True, None
        else:
            return False, "Not in the list of allowed datasets"


class ExclusiveDatasetTagFilter(DatasetFilter):

    name = "ExclusiveDatasetTagFilter"
    
    def __init__(self, dataset_tags):

        super(ExclusiveDatasetTagFilter, self).__init__()

        self.dataset_tags = dataset_tags

    def parameters(self):

        p_ = {
            'dataset_tags' : self.dataset_tags,
        }

        return p_

    def get_dataset_tag(self, dataset):

        if dataset.tag is None: 
            raise ValueError(
                "This function can only be used with labelled datasets with a 'tag' attribute. "
                "Datasets can be labelled using the label method."
                )

        return dataset.tag

    def is_valid(self, dataset, **kwargs):

        if self.get_dataset_tag(dataset) in self.dataset_tags:
            return False, "Dataset in list of excluded datasets"
        else:
            return True, None


class HighResolutionFilter(DatasetFilter):
    
    name = "HighResolutionFilter"

    def __init__(self, high_resolution_cutoff):
        """
        Reject datasets with a high resolution limit worse than
        the provided value.

        min_high_resolution = 5
          Accepted: 4
          Rejected: 6
        """

        super(HighResolutionFilter, self).__init__()

        self.high_resolution_cutoff = high_resolution_cutoff

    def parameters(self):

        p_ = {
            'high_resolution_cutoff' : self.high_resolution_cutoff,
        }

        return p_

    def get_high_resolution(self, dataset):
        return dataset.data.mtz_object().max_min_resolution()[1]

    def is_valid(self, dataset, **kwargs):

        if self.get_high_resolution(dataset) > self.high_resolution_cutoff:
            return False, "High resolution limit above cutoff"

        return True, None


class RValueFilter(DatasetFilter):

    name = "RValueFilter"

    def __init__(self, max_rfree=None, max_rwork=None):

        super(RValueFilter, self).__init__()

        self.max_rfree = max_rfree
        self.max_rwork = max_rwork

    def parameters(self):

        p_ = {
            'max_rfree' : self.max_rfree,
            'max_rwork' : self.max_rwork,
        }

        return p_

    @classmethod
    def get_rfree(cls, dataset):
        return dataset.model.input.get_r_rfree_sigma().r_free

    @classmethod
    def get_rwork(cls, dataset):
        return dataset.model.input.get_r_rfree_sigma().r_work

    def is_valid(self, dataset, **kwargs):

        if (self.max_rfree is not None):
            if self.get_rfree(dataset) > self.max_rfree:
                return False, "R-free is too high"

        if (self.max_rwork is not None):
            if self.get_rwork(dataset) > self.max_rwork:
                return False, "R-work is too high"

        return True, None


class SpaceGroupFilter(DatasetFilter):

    name = "SpaceGroupFilter"

    def __init__(self, reference_dataset=None):

        super(SpaceGroupFilter, self).__init__()

        self.reference_dataset = reference_dataset

    def parameters(self):

        p_ = {
            'reference_dataset' : self.reference_dataset,
        }

        return p_

    @classmethod
    def get_spacegroup(cls, dataset):
        return dataset.model.crystal.space_group.info().symbol_and_number()

    def is_valid(self, dataset, reference_dataset=None, **kwargs):

        if (self.reference_dataset is not None):
            reference_dataset = self.reference_dataset

        ref_sg = self.get_spacegroup(reference_dataset)
        dat_sg = self.get_spacegroup(dataset)

        if ref_sg != dat_sg:
            return False, "Different space group ({} != {})".format(dat_sg, ref_sg)

        return True, None


class IdenticalHierarchyFilter(DatasetFilter):

    name = "IdenticalHierarchyFilter"

    def __init__(self, reference_dataset=None, atom_selection_string=None):

        super(IdenticalHierarchyFilter, self).__init__()

        self.reference_dataset = reference_dataset
        self.atom_selection_string = atom_selection_string

    def parameters(self):

        p_ = {
            'reference_dataset' : reference_dataset,
            'atom_selection_string' : self.atom_selection_string,
        }

        return p_

    def is_valid(self, dataset, reference_dataset=None, **kwargs):

        if (self.reference_dataset is not None):
            reference_dataset = self.reference_dataset

        if not dataset.model.hierarchy.is_similar_hierarchy(
            reference_dataset.model.hierarchy
            ):

            return False, "Different structural compositions"

        return True, None


##


class DatasetFilterGroup(DatasetFilter):
    
    name = "DatasetFilterGroup"

    def __init__(self, filters=None):

        super(DatasetFilterGroup, self).__init__()

        if (filters is None):
            filters = []

        self.filters = filters

    def __call__(self, mcd, **kwargs):

        for f in self.filters:
            mcd = f(mcd=mcd, **kwargs)
            self.rejections.update(f.rejections)

        return mcd

    def __str__(self):

        parameters = self.parameters()
        rejections = self.rejections 
        filters = self.filters

        parameter_strings = [
            '{k}: {v}'.format(
                k=k, 
                v=str(v).strip().replace('\n','\n\t'),
                )
            for k,v in parameters.items()
            ]

        filter_strings = [
            "> Filter {i}\n\t{s}".format(
                i = i+1, 
                s = str(f).strip().replace('\n','\n\t'),
                )
            for i, f in enumerate(filters)
            ]

        rejection_strings = [
            '{k}: {v}'.format(
                k=k, 
                v=str(v).strip().replace('\n','\n\t'),
                )
            for k,v in sorted(rejections.items())
        ]

        s_ = (
            'Filter Type: {name}\n'
            '| Parameters: \n'
            '|\t{parameters}\n'
            '| Filters: \n'
            '|\t{filters}\n'
            '| Rejected Datasets: \n'
            '|\t{rejections}\n'
            '`---->'
            ).format(
            name = self.name,
            parameters = (
                '\n'.join(parameter_strings).strip().replace('\n','\n|\t')
                if parameter_strings else "None"
                ),
            filters = (
                '\n'.join(filter_strings).strip().replace('\n','\n|\t')
                if filter_strings else "None"
                ),
            rejections = (
                '\n'.join(rejection_strings).strip().replace('\n','\n|\t')
                if rejection_strings else "None"
                ),
            )

        return s_.strip()


class DefaultDatasetFilter(DatasetFilterGroup):

    name = "DefaultDatasetFilter"

    def __init__(self,
        same_space_group_only = True,
        similar_models_only = False,
        max_rfree = None,
        max_rwork = None,
        reference_dataset = None,
        ):

        super(DefaultDatasetFilter, self).__init__()

        self.same_space_group_only = bool(same_space_group_only)
        self.similar_models_only = bool(similar_models_only)
        self.max_rfree = max_rfree
        self.max_rwork = max_rwork
        self.reference_dataset = reference_dataset

        # Create the list of filters

        if (self.same_space_group_only is True):
            self.filters.append(
                SpaceGroupFilter(
                    reference_dataset = reference_dataset,
                    )
                )

        if (self.similar_models_only is True):
            self.filters.append(
                IdenticalHierarchyFilter(
                    reference_dataset = reference_dataset,
                    )
                )

        if [self.max_rwork, self.max_rfree].count(None) < 2:
            self.filters.append(
                RValueFilter(
                    max_rfree = max_rfree,
                    max_rwork = max_rwork,
                    )
                )

    def parameters(self):

        p_ = {
            'similar_models_only' : self.similar_models_only,
            'same_space_group_only' : self.same_space_group_only,
            'max_rfree' : self.max_rfree,
            'max_rwork' : self.max_rwork,
            'reference_dataset' : self.reference_dataset,
        }

        return p_


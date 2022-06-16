import giant.logs as lg
logger = lg.getLogger(__name__)

import copy, collections

from giant.mulch.selectors import (
    DatasetKeySelector,
    )

from giant.utils import (
    pretty_format_list,
    )


class UpdateablePartitioner(object):

    name = "DatasetPartitioner"

    def __init__(self):

        self.sets = dict()

    def __str__(self):

        parameters = self.parameters()

        parameter_strings = [
            '{k}: {v}'.format(
                k=k, 
                v=str(v).strip().replace('\n','\n\t'),
                ) 
            for k,v in parameters.items()
            ]

        sets = self.as_dict()

        set_strings = [
            (
                pretty_format_list(
                    values = v,
                    key = k,
                    ).strip().replace('\n','\n\t')
                if (len(v) > 0) 
                else ("{k} : None".format(k=k))
                )
            for k, v in sets.items()
        ]

        s_ = (
            'Partitioner Type: {name}\n'
            '|\n'
            '| Parameters:\n'
            '|\t{parameters}\n'
            '| Sets:\n'
            '|\t{sets}\n'
            '`---->\n'
            ).format(
            name = self.name,
            parameters = (
                '\n'.join(parameter_strings).strip().replace('\n','\n|\t')
                if parameter_strings 
                else "None"
                ),
            sets = (
                '\n'.join(set_strings).strip().replace('\n','\n|\t')
                if set_strings 
                else "None"
                ),
            )

        return s_.strip()

    def parameters(self):
        r_ = {}
        return r_

    def update(self, partitions_dict):

        for p_key, p_tags in partitions_dict.items():

            p_set = self.sets.get(p_key)

            if (p_set is None): 
                self.sets[p_key] = set(p_tags)
            else:
                p_set.update(p_tags)

    def remove(self, partition, values): 

        if partition not in self.sets: 
            raise ValueError(
                'partition "{p}" does not exist'.format(partition)
                )

        p_set = self.sets[partition]

        if p_set is None: 
            return 

        p_set.difference_update(values)

    def as_dict(self):

        return {
            p_key : (
                list(p_set) 
                if (p_set is not None)
                else []
                )
            for p_key, p_set
            in self.sets.items()
        }


class DatasetKeyPartitioner(UpdateablePartitioner):

    name = "DatasetKeyPartitioner"

    def __init__(self, **kw_args):

        self.sets = {
            p_key : set(dataset_keys)
            for p_key, dataset_keys 
            in kw_args.items()
        }

    def __call__(self, datasets):

        partitions = {}

        for p_key, p_set in self.sets.items():

            selector = DatasetKeySelector(
                dataset_keys = list(p_set),
                )

            partitions[p_key] = selector(datasets)

        return partitions


class TestTrainPartitioner(UpdateablePartitioner):

    name = "TestTrainPartitioner"

    def __init__(self,
        test = None,
        train = None,
        not_test = None,
        not_train = None,
        test_selector = None,
        train_selector = None,
        ):

        self.sets = dict(
            test = self._sanitise(test),
            train = self._sanitise(train),
            not_test = self._sanitise(not_test),
            not_train = self._sanitise(not_train),
        )

        self.test_selector = test_selector
        self.train_selector = train_selector

        self.validate()

    def __call__(self, datasets):

        return {
            'test' : self.get_test_datasets(datasets),
            'train' : self.get_train_datasets(datasets),
        }

    def _sanitise(self, name_list):

        if (name_list is None):
            return None

        return set(name_list)

    def validate(self):

        if (self.sets["train"] is not None) and (self.sets["not_train"] is not None):

            overlap = self.sets["train"].intersection(self.sets["not_train"])

            if len(overlap) > 0:
                raise ValueError(
                    "Keys present in set 'train' and set 'not_train': {}".format(
                        list(overlap)
                        )
                    )

        if (self.sets["test"] is not None) and (self.sets["not_test"] is not None):

            overlap = self.sets["test"].intersection(self.sets["not_test"])
            
            if len(overlap) > 0:
                raise ValueError(
                    "Keys present in set 'test' and set 'not_test': {}".format(
                        list(overlap)
                        )
                    )

    def get_train_datasets(self, datasets):

        if (self.sets["not_train"] is not None):
            datasets = {dtag: d for dtag, d in datasets.items() if (dtag not in self.sets["not_train"])}

        if (self.sets["train"] is not None):
            datasets = {dtag: d for dtag, d in datasets.items() if (dtag in self.sets["train"])}

        if (self.train_selector is not None):
            datasets = self.train_selector(datasets)

        return datasets

    def get_test_datasets(self, datasets):

        if (self.sets["not_test"] is not None):
            datasets = {dtag: d for dtag, d in datasets.items() if (dtag not in self.sets["not_test"])}

        if (self.sets["test"] is not None):
            datasets = {dtag: d for dtag, d in datasets.items() if (dtag in self.sets["test"])}

        if (self.test_selector is not None):
            datasets = self.test_selector(datasets)

        return datasets


class ResolutionShellTestTrainPartitioner(object):

    name = "ResolutionShellTestTrainPartitioner"

    def __init__(self,
        shell_thickness = 0.5,
        high_resolution = 0.0,
        low_resolution = 999.0,
        test_train_partitioner = None,
        min_train_datasets = None,
        ):

        self.shell_thickness = shell_thickness
        self.high_resolution = high_resolution
        self.low_resolution = low_resolution
        self.partition_datasets = test_train_partitioner
        self.min_train_datasets = min_train_datasets

    def parameters(self):

        r_ = {
            'shell_thickness' : self.shell_thickness,
            'high_resolution' : self.high_resolution,
            'low_resolution' : self.low_resolution,
            'partition_datasets' : self.partition_datasets,
            'min_train_datasets' : self.min_train_datasets,
        }

        return r_

    def __str__(self):

        parameters = self.parameters()

        test_train_partitioner = parameters.pop('partition_datasets')

        parameter_strings = [
            '{k}: {v}'.format(
                k=k, 
                v=str(v).strip().replace('\n','\n\t'),
                ) 
            for k,v in parameters.items()
            ]

        s_ = (
            'Partitioner Type: {name}\n'
            '|\n'
            '| Parameters:\n'
            '|\t{parameters}\n'
            '|\n'
            '| Test/Train Partitioner:\n'
            '|\t{partitioner}\n'
            '`---->\n'
            ).format(
            name = self.name,
            parameters = (
                '\n'.join(parameter_strings).strip().replace('\n','\n|\t')
                ),
            partitioner = (
                str(test_train_partitioner).strip().replace('\n','\n|\t')
                ),
            )

        return s_.strip()

    def __call__(self, datasets):

        if len(datasets) < self.min_train_datasets:
            raise ValueError("Number of datasets is less than min_train_datasets")

        # Create a copy so we can modify it!
        partition_datasets = copy.deepcopy(self.partition_datasets)

        # Loop termination - limited by params or dataset resolution
        max_resolution = min(
            self.low_resolution,
            max([self.get_resolution(d) for d in datasets.values()])
        )
        # Ensure not less than high resolution
        max_resolution = max(max_resolution, self.high_resolution)

        # Initialise loop variables - initialise at high resolution
        r_outer = self.high_resolution
        r_inner = self.high_resolution

        # Initialise output dict
        output_shells = []

        while r_inner <= max_resolution:

            r_inner = (r_inner + self.shell_thickness)

            # Filter datasets by resolution
            sel_datasets = {
                dtag: d 
                for dtag, d in datasets.items() 
                if (self.get_resolution(d) < r_inner)
                }

            # Any datasets? 
            if (not sel_datasets):
                continue

            #####################
            # Extract these functions to a partition validator! 
            #
            # Get test and train datasets
            dataset_partitions = partition_datasets(sel_datasets)
            #
            if not self.validate_partitions(dataset_partitions):
                continue

            # Get the dataset keys selected for testing
            test_dtags = list(dataset_partitions['test'].keys())

            # Remove them from future test sets (add to excluded datasets)
            partition_datasets.update( # add to negative set
                partitions_dict = {'not_test': test_dtags},
                )

            # Label for this shell
            shell_bounds = (r_outer, r_inner) # () "{:.2f} - {:.2f}A".format(...)

            # Add to output dict
            output_shells.append(
                (shell_bounds, dataset_partitions)
                )

            # Update for next shell
            r_outer = r_inner

        if not output_shells:
            raise ValueError('No shell partitions generated')

        return collections.OrderedDict(output_shells)

    @classmethod
    def get_resolution(cls, dataset):
        return dataset.data.crystal.resolution_high

    def validate_partitions(self, p_dict):

        test_keys = list(p_dict['test'].keys())
        train_keys = list(p_dict['train'].keys())

        # Skip to next if no test datasets
        if (len(test_keys) == 0) or (len(train_keys) == 0):
            return False

        # Same dataset only in test and train? can't compare to itself! 
        if (len(train_keys) == 1) and (train_keys == test_keys):
            return False

        # Skip if not enough training datasets
        if (len(train_keys) < self.min_train_datasets):
            return False

        #
        # MODIFIERS OF THE PARTITIONS
        #

        # One dataset in train and this dataset also in test? Remove from test.
        if (len(train_keys) == 1) and (train_keys[0] in test_keys): 

            p_dict['test'].pop(train_keys[0])

        return True

import giant.logs as lg
logger = lg.getLogger(__name__)

import pandas
import collections


class DatasetStatistics(object):

    def __init__(self, data):
        """
        Output data-holder class from ExtractDatasetStatistics class.
        Manages generation of panda dataframes, generating partition objects, etc.
        """

        self.data = data

    def __call__(self):
        return self.as_dataframe()

    def get_row_names(self):
        return sorted(self.data.keys())

    def get_column_names(self):
        """Extract the superset of all data keys"""
        keys_dict = collections.OrderedDict()
        for dtag, d_dict in self.data.items():
            keys_dict.update(d_dict)
        return list(keys_dict.keys())

    def as_dataframe(self):
        index = self.get_row_names()
        columns = self.get_column_names()
        return pandas.DataFrame(
            index = index,
            columns = columns,
            data = [self.data[k] for k in index],
            )

    def to_csv(self, path):
        dataframe = self.as_dataframe()
        dataframe.to_csv(path)


class ExtractDatasetStatistics(object):

    def __init__(self,
        extracters,
        ):

        self.extracters = extracters

    def __call__(self, mcd, **kw_args):

        statistics_dicts = collections.OrderedDict()

        for dtag, dataset in mcd.datasets.items():

            dataset_dict = statistics_dicts.setdefault(
                dtag,
                collections.OrderedDict(),
                )

            for extract in self.extracters:

                dataset_dict.update(
                    extract(
                        dataset = dataset,
                        **kw_args
                        )
                    )

        statistics = DatasetStatistics(
            data = statistics_dicts,
            )

        return statistics
        

class ClassifyDatasetStatistics(object):

    def __init__(self,
        classifiers,
        combination_approach = "union",
        ):

        self.classifiers = classifiers

        if combination_approach == "union":
            self.combine_sets = lambda main_set, new_set: set(main_set).union(new_set)
        elif combination_approach == "intersection":
            self.combine_sets = lambda main_set, new_set: set(main_set).intersection(new_set)
        else:
            raise ValueError("invalid combination_approach: {}".format(combination_approach))

    def __call__(self, 
        dataframe,
        ):

        classifications = collections.OrderedDict()

        for classify in self.classifiers:

            # Get the classification dict (partition name -> list of ids)
            c_partitions = classify(dataframe)

            # Add to current sets using
            for p_key, p_set in c_partitions.items():

                current_set = classifications.get(p_key, None)

                if (current_set is None):
                    classifications[p_key] = p_set
                else:
                    classifications[p_key] = self.combine_sets(current_set, p_set)

        return classifications


class ExtractAndClassifyDatasetStatistics(object):

    def __init__(self, 
        extracters, 
        classifiers,
        classifier_combination_approach = "union",
        ):

        self.extract = ExtractDatasetStatistics(
            extracters = extracters,
            )
        self.classify = ClassifyDatasetStatistics(
            classifiers = classifiers,
            combination_approach = classifier_combination_approach,
            )

    def __call__(self, mcd, **kw_args):

        self.statistics = self.extract(
            mcd = mcd, 
            **kw_args
            )

        self.dataframe = self.statistics.as_dataframe()

        self.classifications = self.classify(
            dataframe = self.dataframe,
            )

        return self

    def __str__(self):

        s = (
            "Dataset Statistics: \n{dataframe}\n\n"
            "Classifications: \n\t{classifications}\n\n"
            ).format(
            dataframe = str(self.dataframe),
            classifications = '\n\t'.join([
                "{}: {}".format(p,d) 
                for p, d in self.classifications.items()
                ]),
            )
        
        return s

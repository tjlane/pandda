import numpy as np

from giant.common.statistics import modified_z_scores


class ZScoreClassifier(object):

    default_columns = []

    def __init__(self,
        columns = None,
        outlier_partition = "outliers",
        column_suffix = "_z",
        max_z_score = 5.0,
        use_absolute_z_score = True,
        z_score_type = "normal",
        ):

        assert z_score_type in ["normal", "modified"]

        if (columns is None):
            columns = []

        # Concatenate columns
        columns = list(columns) + list(self.default_columns)
        
        if len(columns) == 0:
            raise ValueError("no columns supplied and no default_columns defined.")

        self.columns = columns
        self.outlier_partition = (
            str(outlier_partition) 
            if (outlier_partition is not None) 
            else outlier_partition
            )
        self.column_suffix = (
            str(column_suffix) 
            if (column_suffix is not None) 
            else column_suffix
            )
        self.max_z_score = max_z_score
        self.use_absolute_z_score = bool(use_absolute_z_score)
        self.z_score_type = z_score_type

    def __call__(self, dataframe):

        outliers = set()

        missing_columns = set(self.columns).difference(dataframe.columns)
        if len(missing_columns) > 0:
            raise ValueError(
                "Columns '{columns}' are missing from dataframe: \n{dataframe}".format(
                    columns = "', '".join(map(str,missing_columns)),
                    dataframe = str(dataframe),
                    )
                )

        for c in self.columns:
            # Extract column data
            series = dataframe[c]
            # Convert to zscores
            zscores = self.calculate_z_scores(series.values)
            # Store in input dataframe
            if (self.column_suffix is not None):
                dataframe[c+self.column_suffix] = zscores
            # Exclude negative Z-scores? 
            if (self.use_absolute_z_score is True):
                zscores = np.abs(zscores)
            # Identify values above cutoff
            exclude = (
                zscores > self.max_z_score
                )
            # Create union with existing set
            outliers = outliers.union(
                list(
                    series.index[exclude].values
                    )
                )

        return {
            self.outlier_partition : sorted(outliers),
            }

    def calculate_z_scores(self, values):

        values = np.array(values)

        if self.z_score_type == "normal":
            return (values - np.mean(values)) / np.std(values)

        elif self.z_score_type == "modified":
            return modified_z_scores(values)

        raise ValueError("invalid z_score_type: {}".format(self.z_score_type))

import pathlib as pl
import pandas as pd


class MakePanddaDatasetSummaryTable(object):

    columns = [
        "high_resolution",
        "low_resolution",
        "r_free",
        "r_work",
        "test",
        "train",
        "not_test",
        "not_train",
        "space_group",
        "unit_cell",
        "rejection_reason",
        # !!! others,
        ]

    def __init__(self, output_path=None):
        self.output_path = output_path

    def __call__(self, dataset_info_dicts, output_path=None):

        if (output_path is None):
            output_path = self.output_path

        if not pl.Path(output_path).parent.exists():
            pl.Path(output_path).parent.mkdir(parents=True)

        df = pd.DataFrame(
            data = dataset_info_dicts,
            columns = self.columns,
            )

        df.index.name = 'dtag'

        df = df.sort_values(
            by = 'dtag',
            ascending = True,
            )

        if (output_path is not None):
            df.to_csv(
                path_or_buf = str(output_path),
                )

        return df

import giant.logs as lg
logger = lg.getLogger(__name__)

import pathlib as pl
import pandas as pd

from .tables import (
    MakePanddaDatasetSummaryTable,
    )

from .graphs import (
    MakePanddaDatasetSummaryGraphs,
    )

from .html import (
    MakePanddaDatasetSummaryHtml,
    )

from pandda.utils import (
    show_dict,
    )


class WritePanddaDatasetSummary(object):

    def __init__(self,
        output_dir,
        dataset_dir = None,
        write_table = None,
        write_graphs = None,
        write_html = None,
        # write_json = None,
        ):

        output_dir = pl.Path(output_dir)

        summary_dir = (output_dir / "dataset_summary")

        if (write_table is None):
            write_table = MakePanddaDatasetSummaryTable(
                output_path = str(
                    summary_dir / "input_datasets_info.csv"
                    ),
                )

        if (write_graphs is None):
            write_graphs = MakePanddaDatasetSummaryGraphs(
                output_dir = (
                    summary_dir
                    ),
                )

        if (write_html is None):
            write_html = MakePanddaDatasetSummaryHtml(
                output_directory = str(
                    output_dir / "html"
                    ),
                )

        # if (write_json is None):
        #     write_json = MakePanddaInputSummaryJson(
        #         output_directory = str(
        #             )
        #         )

        self.write_table = write_table
        self.write_graphs = write_graphs
        self.write_html = write_html

    def __call__(self, 
        mcd, 
        dataset_statistics,
        dataloader,
        filters,
        partitions,
        data_getters,
        reference_dataset,
        get_reference_data,
        ):

        output_files = {}

        # Get the dataset information as {columns: {index: data}}
        all_dataset_info = dataset_statistics.dataframe.to_dict()

        self.update_from_filters(
            all_dataset_info = all_dataset_info,
            all_dataset_keys = sorted(mcd.datasets.keys()),
            filters = filters,
            )

        self.update_from_partitions(
            all_dataset_info = all_dataset_info,
            all_dataset_keys = sorted(mcd.datasets.keys()),
            partitions = partitions,
            )

        ###

        logger.subheading('Writing dataset table')
        
        logger(
            self.write_table(
                dataset_info_dicts = all_dataset_info,
                )
            )

        output_files['tables'] = {
            'dataset_info' : str(self.write_table.output_path)
            }

        #

        logger.subheading('Writing dataset graphs')

        output_files['graphs'] = self.write_graphs(
            datasets = mcd.datasets,
            dataset_dicts = all_dataset_info,
            data_getters = data_getters,
            reference_miller_array = get_reference_data(
                reference_dataset,
                ),
            )

        #
        
        output_info = self.extract_html_info(
            mcd = mcd,
            dataloader = dataloader,
            dataset_info = all_dataset_info,
            )

        self.write_html(
            datasets = mcd.datasets,
            output_info = output_info,
            output_graphs = output_files['graphs'],
            output_tables = output_files['tables'],
            )

        output_files['html'] = {
            'dataset_html' : str(self.write_html.output_path)
            }

        #

        logger.subheading('New Output Files')
        show_dict(output_files, logger=logger)

        ###

        return output_files

    def update_from_filters(self, 
        all_dataset_info, 
        all_dataset_keys,
        filters,
        ):

        if (filters is None):
            return

        d_dict = all_dataset_info.setdefault('rejection_reason', {})

        # for dkey in all_dataset_keys: 
        #     if (dkey not in d_dict):
        #         d_dict[dkey] = 'none'

        for f in filters:

            if hasattr(f, "as_filter"):
                f = f.as_filter()

            for dkey, reason in f.rejections.items():

                d_dict[dkey] = reason

    def update_from_partitions(self,
        all_dataset_info,
        all_dataset_keys,
        partitions,
        ):
        
        if (partitions is None):
            return

        for p_name, p_dkeys in partitions.items():

            assert p_name not in all_dataset_info # if this is removed, need to merge dictionaries rather than override

            p_dict = {
                dkey : (dkey in p_dkeys)
                for dkey in all_dataset_keys
            }

            all_dataset_info[p_name] = p_dict

    def extract_html_info(self, 
        mcd,
        dataloader,
        dataset_info,
        ):

        output_info = {
            'rejected_datasets' : dataset_info.get('rejection_reason', {}),
            'loaded_datasets' : list(mcd.datasets.keys()),
        }

        output_info.update(
            dataloader.info()
            )

        return output_info

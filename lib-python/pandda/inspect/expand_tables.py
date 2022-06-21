import pathlib as pl
import pandas as pd

from pandda.analyse.output.tables import (
    MakePanddaEventTable,
    MakePanddaSiteTable,
    )


class ExpandEventTable(object):

    def __init__(self,
        dataset_info_csv,
        dataset_tags,
        output_path = None,
        keep_original_events = False,
        remove_datasets_with_events = False,
        ):
        
        self.dataset_info_csv = (
            dataset_info_csv
            )

        self.dataset_tags = (
            dataset_tags
            )

        self.make_pandda_event_table = MakePanddaEventTable(
            output_path = (
                str(output_path)
                if output_path is not None
                else None
                ),
            )

        self.keep_original_events = (
            keep_original_events
            )

        self.remove_datasets_with_events = (
            remove_datasets_with_events
            )

    def __call__(self, 
        event_table = None, 
        ):

        dataset_table = pd.read_csv(
            str(self.dataset_info_csv), 
            sep = ',', 
            dtype = {'dtag': str},
            )

        dummy_event_dicts = self.make_event_dicts_from_dataset_table(
            dataset_dicts = list(
                dataset_table.to_dict('records')
                ),
            )

        #

        if self.keep_original_events is False:

            event_dicts = []

        else:

            assert event_table is not None

            event_dicts = list(
                event_table.to_records()
                )

            if self.remove_datasets_with_events is True:

                event_dtags = sorted(
                    set([
                        e['dtag']
                        for e in event_dicts
                        ])
                    )

                dummy_event_dicts = [
                    d
                    for d in dummy_event_dicts
                    if d['dtag'] not in event_dtags
                ]

        #

        df = self.make_pandda_event_table(
            event_dicts = (
                event_dicts + dummy_event_dicts
                ),
            )

        return df

    def make_event_dicts_from_dataset_table(self, dataset_dicts):

        event_dicts = []

        defaults = {
            "event_num" : 1,
            "site_num" : 1,
            "event_fraction" : 1.0,
            "bdc" : 0.0,
        }

        for d in dataset_dicts: 

            if d['dtag'] not in self.dataset_tags:
                continue

            e = {
                c : d.get(
                    c, 
                    defaults.get(c, None),
                    ) 
                for c in self.make_pandda_event_table.columns
            }

            event_dicts.append(e)

        return event_dicts


class MakeDummySiteTable(object):

    def __init__(self,
        output_path = None,
        ):

        self.make_pandda_site_table = MakePanddaSiteTable(
            output_path = str(output_path),
            )

    def __call__(self):

        site_dicts = [
            {
            "site_num" : 1,
            "n_events" : None,
            "max_value" : None,
            "xyz_centroid" : None,
            "xyz_extent" : None,
            }
        ]

        df = self.make_pandda_site_table(
            site_dicts = (
                site_dicts
                ),
            )

        return df
import giant.logs as lg
logger = lg.getLogger(__name__)

import copy, collections

import pandas as pd
import pathlib as pl

from pandda.utils import (
    show_dict,
    )


class MakePanddaDatasetTable(object):

    columns = [
        "map_uncertainty",
        "analysed_resolution",
        "test",
        "train",
        "high_resolution",
        "low_resolution",
        "r_free",
        "r_work",
        "space_group",
        "unit_cell",
        ]

    def __init__(self, output_path=None):
        self.output_path = output_path

    def __call__(self, dataset_dicts, output_path=None):

        if (output_path is None):
            output_path = self.output_path

        if not pl.Path(output_path).parent.exists():
            pl.Path(output_path).parent.mkdir(parents=True)

        dataset_dicts = self.get_columns(dataset_dicts)

        df = pd.DataFrame.from_dict(
            data = dataset_dicts,
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

    def get_columns(self, dataset_dicts):

        dataset_dicts = collections.OrderedDict([
            (
                col, 
                dataset_dicts.get(col, {})
                )
            for col in self.columns
            ])

        return dataset_dicts
        

class MakePanddaEventTable(object):

    columns = [
        "dtag",
        "event_num",
        "site_num",
        "event_fraction",
        "bdc",
        "z_peak",
        "z_mean",
        "cluster_size",
        "x",
        "y",
        "z",
        "analysed_resolution",
        "map_uncertainty",
        "global_correlation",
        "local_correlation",
        "r_work",
        "r_free",
        ]

    def __init__(self, output_path=None):
        self.output_path = output_path

    def __call__(self, event_dicts, output_path=None):

        if (output_path is None): 
            output_path = self.output_path

        event_dicts = self.populate_dicts_maybe(
            event_dicts = event_dicts,
            )

        df = pd.DataFrame(
            data = event_dicts,
            columns = self.columns,
            )

        df['interesting'] = False

        df = df.set_index(["dtag", "event_num"])

        # df = df.sort_values(
        #     by = ['site_idx','z_peak'], 
        #     ascending = [1,0],
        #     )

        if (output_path is not None): 
            df.to_csv(
                path_or_buf = str(output_path),
                )

        return df

    def populate_dicts_maybe(self, event_dicts):

        event_dicts = copy.deepcopy(event_dicts)

        for e in event_dicts:

            xyz = e.get('xyz_centroid', (None, None, None))

            if xyz is None: 
                xyz = (None, None, None)

            e.setdefault('x', xyz[0])
            e.setdefault('y', xyz[1])
            e.setdefault('z', xyz[2])

        return event_dicts


class MakePanddaSiteTable(object):

    columns = [
        "site_num",
        "n_events",
        "max_value",
        "xyz_centroid",
        "xyz_extent",
        ]

    def __init__(self, output_path=None):
        self.output_path = output_path

    def __call__(self, site_dicts, output_path=None):

        if (output_path is None):
            output_path = self.output_path

        site_dicts = self.populate_dicts_maybe(
            site_dicts = site_dicts,
            )

        df = pd.DataFrame(
            data = site_dicts,
            columns = self.columns,
            )

        df = df.set_index("site_num")

        df = df.sort_values(
            by = 'site_num',
            ascending = True,
            )

        if (output_path is not None):
            df.to_csv(
                path_or_buf = str(output_path),
                )

        return df

    def populate_dicts_maybe(self, site_dicts):

        site_dicts = copy.deepcopy(site_dicts)

        for s in site_dicts:

            xyz = s.get('xyz_centroid', (None, None, None))

            if xyz is None: 
                xyz = (None, None, None)

            s.setdefault('x', xyz[0])
            s.setdefault('y', xyz[1])
            s.setdefault('z', xyz[2])

        return site_dicts


class MakePanddaResultsTables(object):

    def __init__(self,
        output_directory,
        make_dataset_table = None,
        make_events_table = None,
        make_sites_table = None,
        ):

        output_directory = pl.Path(output_directory)

        if (make_dataset_table is None):
            make_dataset_table = MakePanddaDatasetTable(
                output_path = str(
                    output_directory / "pandda_analyse_datasets.csv"
                    ),
                )

        if (make_events_table is None):
            make_events_table = MakePanddaEventTable(
                output_path = str(
                    output_directory / "pandda_analyse_events.csv"
                    ),
                )

        if (make_sites_table is None):
            make_sites_table = MakePanddaSiteTable(
                output_path = str(
                    output_directory / "pandda_analyse_sites.csv"
                    ),
                )

        self.make_dataset_table = make_dataset_table
        self.make_events_table = make_events_table
        self.make_sites_table = make_sites_table

    def __call__(self,
        dataset_dicts,
        event_dicts,
        site_dicts,
        ):

        ###
        # Create output dicts from events
        #
        all_event_dicts = self.populate_event_dicts(
            dataset_dicts = dataset_dicts,
            event_dicts = event_dicts,
            )
        #
        ###

        ###
        # Write output and return files
        #
        logger.subheading('Output Dataset Table')
        logger(
            self.make_dataset_table(
                dataset_dicts = dataset_dicts,
                )
            )
        #
        logger.subheading('Output Events Table')
        logger(
            self.make_events_table(
                event_dicts = all_event_dicts,
                )
            )
        #
        logger.subheading('Output Sites Table')
        logger(
            self.make_sites_table(
                site_dicts = site_dicts,
                )
            )
        #
        ###

        output_files = {
            'dataset_table' : self.make_dataset_table.output_path,
            'events_table' : self.make_events_table.output_path,
            'sites_table' : self.make_sites_table.output_path,
        }

        return output_files

    def populate_event_dicts(self, 
        dataset_dicts,
        event_dicts,
        ):

        all_event_dicts = []
        
        for d_event in event_dicts:

            d = collections.OrderedDict(d_event)

            dtag = d['dtag']
            for key, value_dict in dataset_dicts.items(): 
                d[key] = value_dict.get(dtag)

            all_event_dicts.append(d)

        return all_event_dicts

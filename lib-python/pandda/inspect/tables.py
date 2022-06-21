import giant.logs as lg
logger = lg.getLogger(__name__)

import os
import pandas as pd
import pathlib as pl

from .exceptions import MissingFile

def validate_path(path):

    if not os.path.exists(str(path)):
        raise MissingFile(
            'File does not exist: {}'.format(str(path))
            )


class TableHandler(object):

    require_not_empty = True

    def __init__(self, 
        input_csv_path,
        output_csv_path,
        ):

        validate_path(input_csv_path)
        
        self.input_path = pl.Path(input_csv_path)
        self.output_path = pl.Path(output_csv_path)

        table = self.read_csv(self.input_path)

        if self.require_not_empty is True:
            if len(table) == 0: 
                logger(str(table))
                raise SystemExit(
                    '\n\nInput csv contains no data!\n\n'
                    )

        self.initialise(table)

        if self.output_path.exists():

            self.update(
                master = table,
                other = self.read_csv(self.output_path),
                )
        
        self.table = table

        self.write()

    def get(self, index, col):

        return self.table[col].iat[index]

    def set(self, index, col, value):

        self.table[col].iat[index] = value

    def write(self):

        logger('Writing output csv: {}'.format(str(self.output_path)))

        self.table.to_csv(str(self.output_path))


class PanddaEventTableHandler(TableHandler):

    def read_csv(self, path):
            
        table = pd.read_csv(str(path), sep=',', dtype={'dtag': str})

        table = table.set_index(['dtag', 'event_num'])

        if self.require_not_empty is True:
            if len(table) == 0: 
                logger(str(table))
                raise SystemExit(
                    '\n\nNo events to inspect: Input events csv contains no data!\n\n'
                    )

        return table

    def initialise(self, table):

        table['Interesting'] = False
        table['Ligand Placed'] = False
        table['Ligand Confidence'] = 'Low'
        table['Comment'] = 'None'
        table['Viewed'] = False
        table['Ligand Names'] = ''

    def update(self, master, other):

        master.update(
            other[
                [
                    'Interesting', 
                    'Ligand Placed', 
                    'Ligand Confidence', 
                    'Comment', 
                    'Viewed',
                    'Ligand Names',
                    ]
                ]
            )


class PanddaSiteTableHandler(TableHandler):

    def read_csv(self, path):
            
        table = pd.read_csv(str(path), sep=',')

        table = table.set_index('site_num')

        if self.require_not_empty is True:
            if len(table) == 0: 
                logger(str(table))
                raise SystemExit(
                    '\n\nNo sites: Input sites csv contains no data!\n\n'
                    )

        return table

    def initialise(self, table):

        table['Name'] = None
        table['Comment'] = None

    def update(self, master, other):

        master.update(
            other[
                [
                    'Name', 
                    'Comment', 
                    ]
                ]
            )
        

class PanddaInspectTableHandler(object): 

    def __init__(self,
        files_dict,
        ):

        self.events = PanddaEventTableHandler(
            input_csv_path = files_dict['input_events'],
            output_csv_path = files_dict['output_events'],
            )

        self.sites = PanddaSiteTableHandler(
            input_csv_path = files_dict['input_sites'],
            output_csv_path = files_dict['output_sites'],
            )

    def write(self):
        self.events.write()
        self.sites.write()

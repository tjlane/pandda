import giant.paths as lg
logger = lg.getLogger(__name__)

import os
import re
import glob
import shutil
import pathlib as p

from itertools import chain
from collections import OrderedDict

class GetPanddaFileTree:

    def __init__(self,
        output_dir,
        ):

        self.output_dir = output_dir

    def __call__(self,
        dataset,
        shells,
        ):

        from giant.paths import Path, File, Dir, MultiIndexedFile, MultiIndexedDir

        tree = Dir(
            name = "{}".format(p.Path(self.output_dir).name),
            root = p.Path(self.output_dir).parent,
            children = {
                "processed_datasets": Dir(
                    name = "processed_datasets",
                    children = {
                        dtag: Dir(
                            name = dtag,
                            children = {
                                "z_map": File(name="{}-z_map.ccp4".format(dtag)),
                                "event_map": MultiIndexedFile(name="{}-event_{}_1-BDC_{}_map.ccp4"),
                                "autobuilt": Dir(
                                    name = "autobuilt",
                                    children = {
                                        "dummy": File(name="dummy"),
                                        "stripped_protein": File("stripped_protein.pdb"),
                                        "event": MultiIndexedDir(name="event_{}_{}"),
                                        "initial_build": MultiIndexedFile(name="initial_build_{}_{}.pdb"),
                                    },
                                ),
                                "initial_data": File(name="{}-pandda-input.mtz".format(dtag)),
                                "initial_model": File(name="{}-pandda-input.pdb".format(dtag)),
                                "ligand_files": Dir(
                                    name="ligand_files",
                                    children={"dummy": File(name="dummy")},
                                ),
                            },
                        ) for dtag, dst
                        in dataset.datasets.items()
                    },
                ),
                "shells": Dir(
                    name = "shells",
                    children = {
                        shell_num: Dir(
                            name = str(shell_num),
                            children = {
                                "mean_map": File(name="mean_map.ccp4"),
                                "event_table": File(name="event_table.csv"),
                            }
                        ) for shell_num, shell_dataset
                        in shells.items()
                    }
                ),
                "event_table": File(name="event_table.csv"),
                "analyses": Dir(
                    name = "analyses",
                    children = {
                        "pandda_analyse_events": File(name="pandda_analyse_events.csv"),
                        "pandda_analyse_sites": File(name="pandda_analyse_sites.csv"),
                        "html_summaries": Dir(
                            name = "html_summaries",
                            children = {
                                "dummy": File(name="dummy"),
                            },
                        ),
                    },
                ),
                "ligand_files": Dir(name="ligand_files",
                                    children={"dummy": File(name="dummy")},
                                    ),
                "modelled_structures": Dir(name="modelled_structures",
                                           children={"dummy": File(name="dummy")},
                                           ),
                }
                   )

        return tree






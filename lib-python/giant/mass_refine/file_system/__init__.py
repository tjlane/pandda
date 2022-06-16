import giant.logs as lg
logger = lg.getLogger(__name__)

import pathlib as pl

from .datasets import (
    DatasetsFileSystem,
    )

from .models import (
    ModelsFileSystem,
    )


class MassRefineFileSystem(object):

    def __init__(self, work_folder):

        self.work_folder = pl.Path(
            work_folder
            )

        self.datasets = DatasetsFileSystem(
            self.work_folder / 'datasets',
            )

        self.models = ModelsFileSystem(
            self.work_folder / 'models',
            )

        self.config_path = (
            self.work_folder / 'config.json'
            )

    def update(self):

        if self.exists(): 
            self.validate()
            self.rebuild()
        else:
            self.initialise()
        
        return self

    def exists(self):
        return self.work_folder.exists()

    def initialise(self):

        assert not self.exists()

        self.work_folder.mkdir(parents=True)
        
        self.datasets.initialise()
        self.models.initialise()

    def validate(self):
        return
        self.datasets.validate()
        self.models.validate()

    def rebuild(self):
        pass



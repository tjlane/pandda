import giant.logs as lg
logger = lg.getLogger(__name__)

import shutil
import pathlib as pl

from .base import (
    FolderWithMeta,
    )

from .refinements import (
    RefinementsFileSystem,
    )


class MakeDatasetFolder(object):

    def __init__(self, dataset_folder):

        self.df = dataset_folder

    def __call__(self, 
        data_path, 
        meta_dict = None,
        ):

        assert (
            not self.df.root_dir.exists()
            ), 'dataset_folder {} already exists'.format(
            str(self.df.root_dir),
            )
        
        self.df.root_dir.mkdir()

        self.copy_data(
            data_path
            )

        self.write_meta(
            meta_dict
            )

        self.df.refinements.initialise()

    def copy_data(self, path):

        shutil.copy(
            str(path),
            str(self.df.data_path),
            )

    def write_meta(self, meta_dict):

        if meta_dict is None: 
            meta_dict = {}

        self.df.write_meta(meta_dict)


class DatasetFolder(FolderWithMeta): 

    name = "Dataset"

    def __init__(self, root_dir):

        self.root_dir = pl.Path(root_dir)

        self.data_path = (
            self.root_dir / 'imported.mtz'
            )

        self.refinements = RefinementsFileSystem(
            root_dir = self.root_dir / 'refinements',
            data_folder = self,
            )

        self.meta_path = (
            self.root_dir / 'meta.json'
            )

    def __str__(self):

        s_ = (
            'Object Type: {name}\n'
            '| Location: {path}\n'
            '| Data: {data}\n'
            '| Meta: \n'
            '| \t{meta_string}\n'
            '`---->'
            ).format(
            name = self.name,
            path = self.root_dir,
            data = self.data(),
            meta_string = '\n'.join([
                "{k:10s} : {v}".format(
                    k=k, v=v,
                    )
                for k, v in sorted(
                    self.get_meta().items()
                    )
                ]).replace('\n','\n| \t'),
            )

        return s_.strip()

    def initialise(self, 
        data_path,
        meta_dict = None,
        ):

        make = MakeDatasetFolder(
            dataset_folder = self,
            )

        make(
            data_path = data_path,
            meta_dict = meta_dict,
            )
        
        return self

    def data(self): 

        return self.data_path


class DatasetsFileSystem(object):

    name = "Datasets"

    def __init__(self, root_dir):

        self.root_dir = pl.Path(root_dir)

    def __str__(self):

        s_ = (
            'Object Type: {name}\n'
            '| Location: {path}\n'
            '| Datasets: \n'
            '| \t{datasets_string}\n'
            '`---->'
            ).format(
            name = self.name,
            path = self.root_dir,
            datasets_string = '\n'.join([
                str(v) for k,v in sorted(
                    self.get_all_as_dict().items()
                    )
                ]).replace('\n','\n| \t'),
            )

        return s_.strip()

    def _path_from_label(self, label):

        return (
            self.root_dir / str(label)
            )

    def initialise(self):

        assert not self.root_dir.exists(), 'already initialised'

        self.root_dir.mkdir()

        return self

    def validate(self):

        return self

    def create_dataset(self, 
        label, 
        data_path,
        meta_dict = None,
        ):

        dir_path = self._path_from_label(
            label,
            )

        assert not dir_path.exists()

        meta_dict = (
            meta_dict if meta_dict is not None else dict()
            )

        meta_dict.setdefault('label', label)

        return DatasetFolder(
            root_dir = dir_path,
            ).initialise(
            data_path = data_path,
            meta_dict = meta_dict,
            )

    def get(self, label):

        dir_path = self._path_from_label(label)

        assert dir_path.exists()

        return DatasetFolder(dir_path)

    def get_all_as_dict(self): 

        return {
            d.stem : self.get(d.stem)
            for d in self.root_dir.glob('*')
            if (
                self.exists(d.stem)
                )
        }

    def exists(self, label):

        dir_path = self._path_from_label(label)

        return dir_path.exists() and dir_path.is_dir()

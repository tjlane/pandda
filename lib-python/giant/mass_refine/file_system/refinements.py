import giant.logs as lg
logger = lg.getLogger(__name__)

import shutil
import pathlib as pl

from .base import (
    FolderWithMeta,
    DirectoryPathGenerator,
    )


class MakeRefinementFolder(object):

    def __init__(self, refinement_folder):

        self.f = refinement_folder

    def __call__(self,
        meta_dict,
        ):

        assert (
            not self.f.root_dir.exists()
            ), 'refinement_folder {} already exists'.format(
            str(self.f.root_dir),
            )

        self.f.root_dir.mkdir()

        self.write_meta(
            meta_dict
            )

    def write_meta(self, meta_dict):

        if meta_dict is None: 
            meta_dict = {}

        self.f.write_meta(meta_dict)


class RefinementFolderPathGenerator(DirectoryPathGenerator):

    prefix = "refine"


class RefinementFolder(FolderWithMeta):

    name = "Refinement"

    def __init__(self, 
        root_dir,
        input_data_folder = None,
        input_model_folder = None,
        ):

        self.root_dir = pl.Path(root_dir)

        self.model_path = (
            self.root_dir / 'output.pdb'
            )

        self.data_path = (
            self.root_dir / 'output.mtz'
            )

        self.input_data = (
            input_data_folder
            )

        self.input_model = (
            input_model_folder
            )

        self.meta_path = (
            self.root_dir / 'meta.json'
            )

    def __str__(self):

        s_ = (
            'Object Type: {name}\n'
            '| Location: {path}\n'
            '| Refined model: {output_model}\n'
            '| Refined data: {output_data}\n'
            '| Input model: \n'
            '| \t{input_model}\n'
            '| Input data: \n'
            '| \t{input_data}\n'
            '| Meta: \n'
            '| \t{meta_string}\n'
            '`---->'
            ).format(
            name = self.name,
            path = self.root_dir,
            input_model = str(
                self.input_model_folder
                ).replace('\n','\n| \t'),
            input_data = str(
                self.input_data_folder
                ).replace('\n','\n| \t'),
            output_model = self.output_model,
            output_data = self.output_data,
            meta_string = '\n'.join([
                "{k:>10s} : {v}".format(
                    k=k, v=v,
                    )
                for k, v in sorted(
                    self.get_meta().items()
                    )
                ]).replace('\n','\n| \t'),
            )

        return s_.strip()

    def initialise(self,
        meta_dict,
        ):

        make = MakeRefinementFolder(
            refinement_folder = self,
            )

        make(
            meta_dict = meta_dict,
            )

        return self

    def validate(self):

        pass

    def model(self):

        if self.model_path.exists():
            return self.model_path
        else:
            return None

    def data(self):

        if self.data_path.exists():
            return self.data_path
        else: 
            return None


class RefinementsFileSystem(object): 

    name = "Refinements"

    def __init__(self, 
        root_dir, 
        data_folder = None,
        ):

        self.root_dir = pl.Path(root_dir)

        self.path_generator = RefinementFolderPathGenerator(
            root_dir = self.root_dir,
            )

        self.data_folder = (
            data_folder
            )

    def __str__(self):

        s_ = (
            'Object Type: {name}\n'
            '| Location: {path}\n'
            '| Refinements: \n'
            '| \t{refinements_string}\n'
            '`---->'
            ).format(
            name = self.name,
            path = self.root_dir,
            refinements_string = '\n'.join([
                str(v) for k,v in sorted(
                    self.get_all_by_number().items()
                    )
                ]).replace('\n','\n| \t'),
            )

        return s_.strip()

    def initialise(self):

        assert not self.root_dir.exists(), 'already initialised'

        self.root_dir.mkdir()

        return self

    def exists(self):

        return self.root_dir.exists()

    def validate(self):

        return self

    def create_refinement_folder(self, 
        label,
        model_folder = None,
        meta_dict = None,
        ):

        dir_path = self.path_generator.get_next(
            suffix = label,
            )

        assert not dir_path.exists()

        meta_dict = (
            meta_dict if meta_dict is not None else dict()
            )
        
        meta_dict.setdefault('label', label)

        return RefinementFolder(
            root_dir = dir_path,
            input_data_folder = self.data_folder,
            input_model_folder = model_folder,
            ).initialise(
            meta_dict = meta_dict,
            )

    def get(self, integer=None, label=None):

        if integer is not None: 

            dir_path_dict = self.path_generator.get_all_as_dict_by_int()
            dir_path = dir_path_dict[integer]

        elif label is not None: 

            dir_path_dict = self.path_generator.get_all_as_dict_by_name()
            dir_path = dir_path_dict[label]

        else: 

            raise Exception('must provide integer or label')

        return RefinementFolder(
            root_dir = dir_path,
            input_data_folder = self.data_folder,
            )

    def get_all_by_number(self):

        return {
            k : RefinementFolder(
                root_dir = p,
                input_data_folder = self.data_folder,
                )
            for k, p in 
            self.path_generator.get_all_as_dict_by_int().items()
        }

    def get_all_by_label(self):

        return {
            k : RefinementFolder(
                root_dir = p,
                input_data_folder = self.data_folder,
                )
            for k, p in 
            self.path_generator.get_all_as_dict_by_name().items()
        }

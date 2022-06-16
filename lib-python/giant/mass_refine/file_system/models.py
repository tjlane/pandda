import giant.logs as lg
logger = lg.getLogger(__name__)

import shutil
import pathlib as pl

from .base import (
    FolderWithMeta,
    DirectoryPathGenerator,
    )


class MakeModelFolder(object):

    def __init__(self, model_folder):

        self.f = model_folder

    def __call__(self, 
        model_path,
        refinement_parameters_paths = None,
        cif_restraints_paths = None,
        meta_dict = None,
        ):

        assert (
            not self.f.root_dir.exists()
            ), 'model_folder {} already exists'.format(
            str(self.f.root_dir),
            )

        self.f.root_dir.mkdir()

        self.copy_model(
            model_path,
            )

        self.copy_refinement_parameters(
            refinement_parameters_paths,
            )

        self.copy_cif_restraints(
            cif_restraints_paths,
            )

        self.write_meta(
            meta_dict,
            )

    def copy_model(self, path):

        shutil.copy(
            str(path),
            str(self.f.model_path),
            )

    def copy_refinement_parameters(self, paths):

        if (paths is None) or not (paths): 
            return

        for i, p in enumerate(paths): 

            ipath = pl.Path(p)

            opath = self.f.root_dir / (
                "refinement-{i:02d}.params".format(
                    i = i+1,
                    )
                )

            shutil.copy(
                str(ipath),
                str(opath),
                )

    def copy_cif_restraints(self, paths):

        if (paths is None) or (not paths): 
            return 

        self.f.cif_dir.mkdir()

        for i, p in enumerate(paths): 

            ipath = pl.Path(p)

            opath = self.f.cif_dir / (
                "{i:02d}-{p}".format(i=i+1, p=ipath.name)
                )

            shutil.copy(
                str(ipath),
                str(opath),
                )

    def write_meta(self, meta_dict):

        if meta_dict is None: 
            meta_dict = {}

        self.f.write_meta(meta_dict)


class ModelFolderPathGenerator(DirectoryPathGenerator):

    prefix = "model"


class ModelFolder(FolderWithMeta): 

    name = "Model"

    def __init__(self, root_dir):

        self.root_dir = pl.Path(root_dir)

        self.meta_path = (
            self.root_dir / 'meta.json'
            )

        self.model_path = (
            self.root_dir / 'model.pdb'
            )

        self.cif_dir = (
            self.root_dir / 'cif_files'
            )

    def __str__(self):

        s_ = (
            'Object Type: {name}\n'
            '| Location: {path}\n'
            '| Model: {model}\n'
            '| Info: \n'
            '| \t{meta_string}\n'
            '`---->'
            ).format(
            name = self.name,
            path = self.root_dir,
            model = self.model(),
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
        model_path,
        refinement_parameters_paths = None,
        cif_restraints_paths = None,
        meta_dict = None, 
        ):

        make = MakeModelFolder(
            model_folder = self,
            )

        make(
            model_path = model_path,
            refinement_parameters_paths = refinement_parameters_paths,
            cif_restraints_paths = cif_restraints_paths,
            meta_dict = meta_dict,
            )

        return self

    def validate(self):

        pass
        
    def model(self):
        
        if self.model_path.exists():
            return self.model_path
        else: 
            raise IOError(
                'model does not exist: {}'.format(
                    self.model_path
                    )
                )

    def refinement_parameters(self):

        ref_params = sorted(
            self.root_dir.glob(
                "*.params"
                )
            )

        if len(ref_params) == 0: 
            return None

        return ref_params

    def cif_restraints(self):

        if not self.cif_dir.exists():
            return None
            
        cif_files = sorted(
            self.cif_dir.glob(
                "*.cif"
                )
            )

        if len(cif_files) == 0: 
            return None

        return cif_files


class ModelsFileSystem(object):

    name = "Models"

    def __init__(self, root_dir):

        self.root_dir = pl.Path(root_dir)

        self.path_generator = ModelFolderPathGenerator(
            root_dir = self.root_dir,
            )

    def __str__(self):

        s_ = (
            'Object Type: {name}\n'
            '| Location: {path}\n'
            '| Models: \n'
            '| \t{models_string}\n'
            '`---->'
            ).format(
            name = self.name,
            path = self.root_dir,
            models_string = '\n'.join([
                str(r) for r in self.get_all()
                ]).replace('\n','\n| \t'),
            )

        return s_.strip()

    def initialise(self):

        assert not self.root_dir.exists(), 'already initialised'

        self.root_dir.mkdir()

        return self
        
    def validate(self):

        return self

    def import_model(self, 
        label, 
        model_path,
        refinement_parameters_paths = None,    
        cif_restraints_paths = None,
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

        return ModelFolder(
            root_dir = dir_path,
            ).initialise(
            model_path = model_path,
            refinement_parameters_paths = refinement_parameters_paths,
            cif_restraints_paths = cif_restraints_paths,
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

        return ModelFolder(dir_path)

    def get_all_by_number(self):

        return {
            k : ModelFolder(p)
            for k, p in 
            self.path_generator.get_all_as_dict_by_int().items()
        }

    def get_all_by_label(self):

        return {
            k : ModelFolder(p)
            for k, p in 
            self.path_generator.get_all_as_dict_by_name().items()
        }

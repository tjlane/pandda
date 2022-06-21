import giant.logs as lg 
logger = lg.getLogger(__name__)

import collections
import pathlib as pl

from pandda.utils import (
    show_dict,
    )


class CopyDatasetPath(object):

    out_dir_template = '{dtag}'
    out_path_template = None

    def __init__(self, mode="link"):
        self.mode = mode

    def __call__(self, dataset, target_dir):
        fmt_dict = self.get_dict(dataset)
        src_path = self.get_source_path(dataset)
        out_dir = pl.Path(target_dir) / self.get_output_directory(**fmt_dict)
        if not out_dir.exists():
            out_dir.mkdir(parents=True)
        out_path = (
            out_dir / self.get_output_filename(**fmt_dict)
            )
        self.copy_file(
            source_path = src_path,
            output_path = out_path,
            )
        return str(out_path)

    def get_tag(self, dataset):
        return dataset.tag

    def get_dict(self, dataset):
        return {
            'dtag': self.get_tag(dataset),
            }

    def get_source_path(self, dataset):
        return None

    def get_output_directory(self, **kw_args):
        return pl.Path(
            self.out_dir_template.format(
                **kw_args
                )
            )

    def get_output_filename(self, **kw_args):
        return pl.Path(
            self.out_path_template.format(
                **kw_args
                )
            )

    def copy_file(self, source_path, output_path):
        if (self.mode == "link"):
            from giant.paths import rel_symlink
            rel_symlink(
                str(source_path),
                str(output_path),
                )
        elif (self.mode == "copy"):
            import shutil
            shutil.copy(
                str(source_path),
                str(output_path),
                )
        else:
            raise NotImplemented('mode "{}" not valid'.formats(self.mode))


class CopyInputStructure(CopyDatasetPath):

    out_path_template = "{dtag}-pandda-input.pdb"

    def get_source_path(self, dataset):
        return pl.Path(dataset.model.filename)


class CopyInputData(CopyDatasetPath):

    out_path_template = "{dtag}-pandda-input.mtz"

    def get_source_path(self, dataset):
        return pl.Path(dataset.data.filename)


class FindAndCopyLigands(CopyDatasetPath):

    out_dir_template = "{dtag}/ligand_files"
    out_path_template = "{dtag}-pandda-ligand-{num:03d}.cif"

    def __init__(self,
        mode = "copy",
        regex = "*.cif",
        rel_search_path = None,
        find_and_copy_matching_pdbs = True,
        rename_files = False,
        ):

        super(FindAndCopyLigands, self).__init__(
            mode = mode,
            )

        self.regex = regex
        self.rel_search_path = rel_search_path
        self.find_and_copy_matching_pdbs = find_and_copy_matching_pdbs
        self.rename_files = rename_files

    def __call__(self, dataset, target_dir):

        fmt_dict = self.get_dict(dataset)

        output_files = []

        out_dir = pl.Path(target_dir) / self.get_output_directory(**fmt_dict)
        if not out_dir.exists():
            out_dir.mkdir(parents=True)

        for i, src_path in enumerate(self.get_ligand_paths(dataset)):

            if (self.rename_files is True):
                out_path = out_dir / self.get_output_filename(num=i+1, **fmt_dict)
            else:
                out_path = out_dir / src_path.name

            self.copy_file(
                source_path = src_path,
                output_path = out_path,
                )

            output_files.append(
                str(out_path)
                )

            if (self.find_and_copy_matching_pdbs is True):
                src_pdb = src_path.with_suffix('.pdb')
                if src_pdb.exists():
                    out_pdb = out_path.with_suffix('.pdb')
                    self.copy_file(
                        source_path = src_pdb,
                        output_path = out_pdb,
                        )

        return output_files

    def get_ligand_paths(self, dataset):
        search_dir = pl.Path(dataset.data.filename).parent
        if (self.rel_search_path is not None):
            search_dir = (search_dir / self.rel_search_path)
        return search_dir.glob(self.regex)


#########


class CopyPanddaInputFiles(object):

    def __init__(self, output_dir):

        self.copy_structure = CopyInputStructure()
        self.copy_data = CopyInputData()
        self.copy_ligands = FindAndCopyLigands()

        self.output_dir = output_dir

    def __call__(self, mcd):

        output_files = collections.OrderedDict()

        for dtag, dataset in mcd.datasets.items():

            output_files[dtag] = self.copy_dataset_files(
                dataset = dataset,
                )

        logger.subheading('Files copied to output folder')
        show_dict(output_files, logger=logger)

        return output_files

    def copy_dataset_files(self, dataset):

        output_files = collections.OrderedDict()

        output_files['structure'] = self.copy_structure(
            dataset = dataset,
            target_dir = self.output_dir,
            )

        output_files['data'] = self.copy_data(
            dataset = dataset,
            target_dir = self.output_dir,
            )

        output_files['ligands'] = self.copy_ligands(
            dataset = dataset,
            target_dir = self.output_dir,
            )
        
        return output_files

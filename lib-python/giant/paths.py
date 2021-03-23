import os
import re
import glob
import shutil

import pathlib as pl

from itertools import chain
from collections import OrderedDict

def filename(f):
    return os.path.basename(os.path.splitext(str(f))[0])

def foldername(f):
    return os.path.basename(os.path.dirname(str(f)))

def easy_directory(directory):
    """Checks a directory exists and creates it if not"""
    directory = str(directory)
    if not os.path.exists(directory):
        up_dir = os.path.dirname(directory)
        if (up_dir != '') and (not os.path.exists(up_dir)):
            easy_directory(up_dir)
        if not os.path.exists(directory):
            os.mkdir(directory)
    return directory

def resolve_glob(path, n=None):
    import glob
    glob_list = glob.glob(path)
    if (n is not None) and (len(glob_list) != n):
        raise ValueError(
            'glob did not return expected number of results ({}): {}'.format(n,glob_list)
        )
    if n == 1:
        return glob_list[0]
    return sorted(glob_list)

def is_available(program):
    """Return whether the program is available in the path"""
    from distutils.spawn import find_executable
    if find_executable(program):
        return True
    return False

def not_available(programs):
    """Return list of programs that are not available in the path"""
    missing = []
    for p in programs:
        if not is_available(p):
            missing.append(p)
    return missing

def rel_symlink(orig, link):
    """Make a relative symlink from link to orig"""
    orig = str(orig)
    link = str(link)
    assert not os.path.exists(link), 'Link already exists: {!s}'.format(link)
    orig = os.path.abspath(orig)
    link = os.path.abspath(link)
    assert not link.endswith('/'), 'LINK CANNOT END WITH /'
    os.symlink(os.path.relpath(orig, start=os.path.dirname(link)), link)

def splice_ext(path, new, position=-1):
    dirname, basename = os.path.split(path)
    split = basename.split(os.path.extsep)
    assert abs(position) <= len(split), 'invalid position selected ({}) to insert "{}" into {}'.format(position, new, split)
    spliced = split[:position] + [new] + split[position:]
    joined = os.path.extsep.join(spliced)
    return os.path.join(dirname, joined)

def compress_file(filename, delete_original=True):
    """Compress a file with gzip"""
    import gzip
    zip_file = filename + '.gz'
    with open(filename, 'r') as fh:
        f = gzip.open(zip_file, 'wb')
        f.write(fh.read())
        f.close()
    if delete_original:
        os.remove(filename)
    return zip_file

def decompress_all_files_in_directory(path, recursive=True, extension='.gz', delete_original=True):
    """Decompress in-place all files in a directory with a particular extension"""
    assert extension in ['.gz'], 'other extensions not yet implemented'
    assert os.path.exists(path), 'path does not exist: {}'.format(path)
    for cur_dir, _, cur_files in os.walk(path):
        for item in cur_files:
            if not item.endswith(extension): continue
            f_zip = os.path.join(cur_dir, item)
            f_out = os.path.splitext(f_zip)[0]
            assert not os.path.exists(f_out), 'output file already exists: {}'.format(f_out)
            if extension == '.gz':
                import gzip
                with gzip.open(f_zip) as fh_zip:
                    with open(f_out, 'w') as fh_out:
                        fh_out.write(fh_zip.read())
            else:
                pass
            assert os.path.exists(f_out), 'output file does not exist: {}'.format(f_out)
            if delete_original is True:
                os.remove(f_zip)
            #zip_obj = zipfile.ZipFile(item_path, 'r')
            #zip_obj.extractall(output_path)
            #zip_obj.close()
        # Breaking in top directory == non-recursive
        if recursive is False:
            break


class Tree:

    def __init__(self, root, structure):

        self.root = root
        self.structure = structure

        self.build_dirs_recursive(root, structure)

    def __call__(self, res):

        paths = self.get_paths_recursive(self.structure, res)

        root_path = p.Path(self.root)

        # append root
        paths = [str(root_path / path) for path in paths]

        return paths

    def get_paths_recursive(self, struc, regexes):

        if len(regexes) == 1:
            return [key for key in struc if re.match(regexes[0], key)]
        else:
            # get all paths that match the rest of the re
            list_of_lists_of_paths = [
                [p.Path(key) / x for x in self.get_paths_recursive(struc[key], regexes[1:])]
                for key in struc
                if re.match(regexes[0], key)
            ]

            # unpack
            list_of_paths = list(chain.from_iterable(list_of_lists_of_paths))

            # return
            return list_of_paths

    def build_dirs_recursive(self, root, structure):

        root_path = p.Path(root)

        for key in structure:
            # Get paths to items
            path = root_path / key

            # recurse if not a file
            if structure[key] is not None:
                try:
                    os.mkdir(str(path))
                    self.build_dirs_recursive(path, structure[key])
                except Exception as e:
                    print(e)

    def update(self, new):

        # Update structure
        self.structure.update(new)

        # Build new
        self.build_dirs_recursive(self.root, new)


class Trace:

    def __init__(
        self,
        path,
        dirct,
        ):

        self.path = path
        self.dirct = dirct

    def __getitem__(
        self,
        item,
        ):
        return Trace(
            self.path / self.dirct.children[item].name,
            self.dirct.children[item],
        )

    def __call__(self, *args):
        return self.dirct(
            path = self.path,
            *args
        )

    def make(self, overwrite=True):
        self.dirct.make(
            path = self.path,
            overwrite = overwrite,
        )


class Dir:

    def __init__(
        self,
        name = None,
        children = None,
        root = None,
        ):

        self.name = name

        if root:
            self.path = p.Path(root) / name
        else:
            self.path = p.Path(name)

        self.children = children

    def __getitem__(self, item):
        return Trace(
            self.path / self.children[item].name,
            self.children[item],
        )

    def __call__(self, path=None):
        if not path:
            path = self.path
        return path

    def make(self, path=None, overwrite=True):
        if not path:
            path = self.path

        if overwrite:
            shutil.rmtree(
                str(path),
                ignore_errors=True,
            )

        os.mkdir(str(path))

        for child_name, child in self.children.items():
            Trace(
                path / child.name,
                child,
            ).make(
                overwrite = overwrite,
            )


class File:

    def __init__(
        self,
        name = None,
        root = None,
        ):

        self.name = name

        if root:
            self.path = p.Path(root) / self.name
        else:
            self.path = p.Path(self.name)

    def __call__(self, path=None):
        if not path:
            path = self.path
        return path

    def make(self, path=None, overwrite=True):
        return


class IndexedFile:

    def __init__(
        self,
        name = None,
        root = None,
        ):
        self.name = name

        if root:
            self.path = p.Path(root) / self.name
        else:
            self.path = p.Path(self.name)

    def __call__(self, index, path=None):
        if not path:
            path = self.path
        return p.Path(str(path).format(index))

    def make(self, path=None, overwrite=True):
        return


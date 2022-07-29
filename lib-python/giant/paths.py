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
    with open(filename, 'rb') as fh:
        f = gzip.open(zip_file, 'wb')
        f.write(
            fh.read()
            )
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
                with gzip.open(f_zip, 'rb') as fh_zip:
                    with open(f_out, 'wb') as fh_out:
                        fh_out.write(
                            fh_zip.read()
                            )
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

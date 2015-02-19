#! /usr/local/python/python2.7.3-64bit/bin/python

# The NEW File objects

import os, datetime

class output_file_object(object):
    def __init__(self, rootdir):
        """Creates and stores output files to be neatly created and retrieved"""

        self.topdir = rootdir
        self.output_dirs = {'root':rootdir}
        self.output_files = {}

        # Check root exists
        assert os.path.exists(rootdir)

    def add_dir(self, dir_name, dir_tag, top_dir_tag=None):
        """Store a directory `dir_name` under `dir_tag` in directory under `top_dir_tag`"""
        # Check that it hasn't been created already
        assert dir_tag not in self.output_dirs.keys()
        # Check the directory that it's being added to exists
        if top_dir_tag == None: top_dir_tag = 'root'
        top_dir_name = self.output_dirs[top_dir_tag]
        # Create dirname and store
        self.output_dirs[dir_tag] = os.path.join(top_dir_name, dir_name)
        # Create if it doesn't exist
        if not os.path.exists(self.output_dirs[dir_tag]):
            os.mkdir(self.output_dirs[dir_tag])
            assert os.path.exists(self.output_dirs[dir_tag])

    def get_dir(self, dir_tag):
        """Retrieve a dirname by it's dir_tag"""
        assert dir_tag in self.output_dirs.keys()
        return self.output_dirs[dir_tag]

    def add_file(self, file_name, file_tag, dir_tag=None):
        """Store a filename `file_name` under `file_tag` in directory under `dir_tag`"""
        # Check that it hasn't been created already
        assert file_tag not in self.output_files.keys()
        # Check the directory that it's beeing added to exists
        if dir_tag == None: dir_tag = 'root'
        dir_name = self.output_dirs[dir_tag]
        # Create filename and store
        self.output_files[file_tag] = os.path.join(dir_name, file_name)

    def get_file(self, file_tag):
        """Retrieve a filename by it's file_tag"""
        assert file_tag in self.output_files.keys()
        return self.output_files[file_tag]

class fileObj(object):
    def __init__(self, file):
        """File information object"""

        self.input = file

        self.path = os.path.abspath(file)
        self.tag = None
        self.dir, self.name = os.path.split(self.path)
        self.base, self.ext = os.path.splitext(self.name)

        stats = os.lstat(file)

        self.created  = datetime.datetime.fromtimestamp(stats.st_ctime)
        self.modified = datetime.datetime.fromtimestamp(stats.st_mtime)
        self.accessed = datetime.datetime.fromtimestamp(stats.st_atime)

    def __call__(self):
        return self.path

def easy_directory(directory):
    """Checks a directory exists and creates it if not"""
    if not os.path.exists(directory):
        os.mkdir(directory)
    return directory

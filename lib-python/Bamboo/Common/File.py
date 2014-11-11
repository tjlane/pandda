#! /usr/local/python/python2.7.3-64bit/bin/python

# The NEW File objects

import os, datetime

class fileObj(object):
    """File information object"""

    def __init__(self, file):

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


#! /usr/local/python/python2.7.3-64bit/bin/python

import os, datetime

class modelHistory(object):
    """Class for containing the modelling history of a Model"""

    def __init__(self, name='modelHistory'):

        self.name = name
        self.history = []
        self.original = None

    def __getitem__(self, index):
        return self.history[index]

    def __str__(self):
        return 'modelHistory: {!s}'.format(self.name)

    def __iter__(self):
        return iter(self.history)

    def set_original_file(self, origfile):
        """Marks a file as the first file"""
        self.original = origfile

    def add_file(self, newfile, index=0):
        """Appends a new model to the beginning of the history list"""
        self.history.insert(index, newfile)

    def get_current_file(self):
        """Returns the newest model (top of the list)"""
        return self.history[0]


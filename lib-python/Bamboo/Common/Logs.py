#! /usr/local/python/python2.7.3-64bit/bin/python

import os, datetime

class ligandRecord(object):

    def __str__(self):
        return 'LigandRecord: {!s}'.format(self.name)

    def __call__(self):
        pass

    def __init__(self, name='ligand', smile=None, builder=None):

        # Name of Ligand
        self.name = name
        # Program used to generate the ligand model
        self.builder = builder
        # Initialise
        self.smile = smile
        # Ligand Model (Generated from smile, random pose)
        self.unfitted = None
        self.restraints = None
        # Correctly oriented and positioned ligand file (relative to protein)
        self.fitted = None
        # Fitted ligand merged with protein structure PDB
        self.merged = None
        # Refined Merged
        self.refined = None
        # LogFiles
        self.buildinglog = None
        self.fittinglog = None
        self.refininglog = None

class fileHistory(object):

    def __init__(self, name=None, path=None):

        # Name of object
        self.name = name
        # File history
        self.history = ModelHistory('{!s} History'.format(self.name))
        # Current Models (most recent)
        self.current = path
        self.previous = None
        self.history.setOriginalFile(self.current)

    def __call__(self):
        return self.current

class pdbRecord(fileHistory):

    def __init__(self, name='PDBRecords', path=None):

        super(pdbRecord, self).__init__(name, path)

        # Reference PDB
        self.reference = None
        self.apo = None

    def __setattr__(self, name, value):

        if value == None:
            self.__dict__[name] = value
        elif name == 'current':
            if isinstance(value,str):
                self.__dict__['previous'] = self.__dict__[name]
                self.__dict__[name] = PDBFileObj(value)
                self.history.addFile(self.__dict__[name])
            else:
                raise ValueError('File Path must be of type str')
        elif name == 'previous':
            print('Cannot set self.previous - not allowed')
        elif name in ['reference']:
            self.__dict__[name] = PDBFileObj(value)
            if not self.current:
                self.__dict__['current'] = self.__dict__[name]
                self.history.addFile(self.__dict__[name])
            else:
                self.history.addFile(self.__dict__[name], index=-1)
        elif name in ['apo']:
            self.__dict__[name] = PDBFileObj(value)
            if not self.current:
                # No files - add as current
                self.__dict__['current'] = self.__dict__[name]
                self.history.addFile(self.__dict__[name])
            elif len(self.history.history)==1 and self.__dict__['reference']:
                # Only reference - add as current
                self.__dict__['current'] = self.__dict__[name]
                self.history.addFile(self.__dict__[name])
            elif self.current and not self.__dict__['reference']:
                # Files present, but not reference - put as earliest file
                self.history.addFile(self.__dict__[name], index=-1)
            else:
                # Multiple present, as is reference - set as second file
                self.history.addFile(self.__dict__[name], index=-2)
        else:
            self.__dict__[name] = value

        # TODO add RDKIT autoparsing?

class mtzRecord(fileHistory):

    def __init__(self, name='MTZRecords', path=None):

        super(mtzRecord, self).__init__(name, path)

        # Raw (No Phases) & Apo (Initial Density)
        self.raw = None
        self.apo = None

    def __setattr__(self, name, value):

        if value == None:
            self.__dict__[name] = value
        elif name == 'current':
            if isinstance(value,str):
                self.__dict__['previous'] = self.__dict__[name]
                self.__dict__[name] = MTZFileObj(value)
                self.history.addFile(self.__dict__[name])
            else:
                raise ValueError('File Path must be of type str')
        elif name == 'previous':
            print('Cannot set self.previous - not allowed')
        elif name in ['raw']:
            self.__dict__[name] = MTZFileObj(value)
            if not self.current:
                # No files - add as current
                self.__dict__['current'] = self.__dict__[name]
                self.history.addFile(self.__dict__[name])
            else:
                self.history.addFile(self.__dict__[name], index=-1)
        elif name in ['apo']:
            self.__dict__[name] = PDBFileObj(value)
            if not self.current:
                self.__dict__['current'] = self.__dict__[name]
                self.history.addFile(self.__dict__[name])
            elif len(self.history.history)==1 and self.__dict__['raw']:
                # Only raw - add as current
                self.__dict__['current'] = self.__dict__[name]
                self.history.addFile(self.__dict__[name])
            elif self.current and not self.__dict__['raw']:
                # Files present, but not raw - put as earliest file
                self.history.addFile(self.__dict__[name], index=-1)
            else:
                # Multiple present, as is raw - set as second file
                self.history.addFile(self.__dict__[name], index=-2)
#        elif name in ['raw','apo']:
#            self.__dict__[name] = MTZFileObj(value)
#            if not self.current:
#                self.__dict__['current'] = self.__dict__[name]
#            self.history.addFile(self.__dict__[name])
        else:
            self.__dict__[name] = value

class modelHistory(object):
    """Class for containing the modelling history of a Model"""

    def __init__(self, name='ModelHistory'):

        self.name = name
        self.history = []
        self.original = None

    def __getitem__(self, index):
        return self.history[index]

    def __str__(self):
        return 'ModelHistory: {!s}'.format(self.name)

    def __iter__(self):
        return iter(self.history)

    def setOriginalFile(self, origfile):
        """Marks a file as the first file"""
        self.original = origfile

    def addFile(self, newfile, index=0):
        """Appends a new model to the beginning of the history list"""
        self.history.insert(index, newfile)

    def getCurrentFile(self):
        """Returns the newest model (top of the list)"""
        return self.history[0]

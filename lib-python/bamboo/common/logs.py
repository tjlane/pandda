import os, sys, datetime

class Log(object):
    def __init__(self, log_file=None, stdout=sys.stdout, verbose=False):
        """Log Object for writing to logs...?"""
        # File that the log will be written to
        self.log_file = log_file
        self.stdout = stdout
        self.verbose = verbose

    def __call__(self, message, show=False, hide=False):
        """Log message to file, and mirror to stdout if verbose or force_print (hide overrules show)"""

        # Format message
        if not isinstance(message, str):
            message = str(message)

        # Print to stdout
        if (not hide) and (show or self.verbose):
            self.show(message)

        if self.log_file:
            # Remove \r from message as this spoils the log (^Ms)
            message = message.replace('\r','')
            # Write to file
            self.write(message=message+'\n', mode='a')

    def show(self, message):
        print(message)

    def write(self, message, mode='a'):
        with open(self.log_file, mode) as fh: fh.write(message)

    def read_all(self):
        return open(self.log_file, 'r').read()

##################################################################
#
#   TODO Redo all code below this line - and move to utils TODO
#
##################################################################

class LigandRecord(object):

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

class ModelHistory(object):
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

    def set_original_file(self, origfile):
        """Marks a file as the first file"""
        self.original = origfile

    def add_file(self, newfile, index=0):
        """Appends a new model to the beginning of the history list"""
        self.history.insert(index, newfile)

    def get_current_file(self):
        """Returns the newest model (top of the list)"""
        return self.history[0]

class FileHistory(object):

    def __init__(self, name=None, path=None):

        # Name of object
        self.name = name
        # File history
        self.history = ModelHistory('{!s} History'.format(self.name))
        # Current Models (most recent)
        self.current = path
        self.previous = None
        self.history.set_original_file(self.current)

    def __call__(self):
        return self.current


#! /usr/local/python/python2.7.3-64bit/bin/python

# Contains General File Utility Functions

import os, sys, shutil

#:::::::::::::::::::::::::::::::::#
# ############################### #
# ### File Utility Functions  ### #
# ############################### #
#:::::::::::::::::::::::::::::::::#

def delete_temporary_directory(tempdir):
    """Delete A Temporary Directory"""

    # Remove Temporary Files
    if tempdir.startswith('/tmp/'):
        shutil.rmtree(tempdir)
        return 0
    else:
        FlagError("Will Not Delete This Directory - Not in '/tmp' - Remove Manually : "+tempdir)
        return 1

def list_directory_contents(directory, templatename=None, templatestyle=None):
    """List Folders in a Directory"""

    # List the contents of the directory
    contents = os.listdir(directory)
    # Decide which are directories and which to keep
    if   templatename and templatestyle=='End':
        subdirs = [dir for dir in contents if os.isdir(os.path.join(directory,dir)) and dir.endswith(templatename)]
    elif templatename and templatestyle=='Start':
        subdirs = [dir for dir in contents if os.isdir(os.path.join(directory,dir)) and dir.startswith(templatename)]
    elif templatename and templatestyle=='Contains':
        subdirs = [dir for dir in contents if os.isdir(os.path.join(directory,dir)) and templatename in dir]
    else:
        subdirs = [dir for dir in contents if os.isdir(os.path.join(directory,dir))]
    # Construct Full Paths
    subdirs = [os.path.join(directory,dir) for dir in subdirs]

    return subdirs


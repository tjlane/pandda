#! /usr/local/python/python2.7.3-64bit/bin/python

import os, sys

def save_compounds_to_temporary_files(residues, prefix='Temp_Residues'):
    """Save a List of Residues to individual PDB files"""

    import tempfile

    # Create Temporary Directory for Files
    tempdir = tempfile.mkdtemp()
    # List of Output Files
    outputfiles = []
    # Save all of the residues to Temporary Files
    for i in xrange(0,len(residues)):
        # Generate Filename
        tempfile = os.path.join(tempdir,prefix+'-Compound_'+str(i)+'.pdb')
        # Append to OutputFiles
        outputfiles.append(tempfile)
        # Save File
        WriteStructureToFile(residues[i],tempfile)

    return tempdir, outputfiles


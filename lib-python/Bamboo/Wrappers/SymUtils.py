#! /usr/local/python/python2.7.3-64bit/bin/python

import os, sys, shutil

from Bamboo.Common.Command import CommandManager
from Bamboo.Utils.Rdkit.Mol import check_pdb_readable

def generate_symmetry_mates(pdbin, pdbout, sgno, cell):
    """Takes the input pdb and generates the unit cell from the point group symmetry"""

    if not pdbout.endswith('.pdb'):
        pdbout = pdbout+'.pdb'
    if not os.path.exists(pdbin):
        raise IOError('pdbin does not exist! {!s}'.format(pdbin))
    if os.path.exists(pdbout):
        raise Exception('pdbout already exists! {!s}'.format(pdbout))

    assert isinstance(sgno, int), 'SPACEGROUP MUST BE AN INTEGER! {!s}'.format(sgno)
    assert isinstance(cell, list), 'CELL MUST BE A LIST! {!s}'.format(cell)

    # Initialise Commander
    PDBSET = CommandManager('pdbset')
    # Set Command Arguments
    PDBSET.SetArguments('XYZIN',os.path.abspath(pdbin),'XYZOUT',os.path.abspath(pdbout))
    # Set inputs
    PDBSET.SetInput(['SYMGEN {!s}'.format(sgno),'CELL {!s}'.format(' '.join(map(str,cell)))])
    # Add Termination
    PDBSET.AppendInput('END')
    # Run!
    PDBSET.Run()

    if not os.path.exists(pdbout):
        raise ExternalProgramError('PDBSET has failed to generate SYMMETRY mates. {!s}\nCOM: {!s}\nOUT: {!s}\nERR: {!s}'.format(pdbin,PDBSET.command,PDBSET.out,PDBSET.err))

    return pdbout

def map_to_reference_using_symmetry(refpdb, movpdb, pdbout, conrad=5):
    """Transforms `movpdb` to the closest symmetry site to `refpdb` using csymmatch - Symmetry info must be contained in the header of the pdb file"""

    if not pdbout.endswith('.pdb'):
        pdbout = pdbout+'.pdb'
    if not os.path.exists(refpdb):
        raise IOError('refpdb does not exist! {!s}'.format(refpdb))
    if not os.path.exists(movpdb):
        raise IOError('movpdb does not exist! {!s}'.format(movpdb))
    if os.path.exists(pdbout):
        raise Exception('pdbout already exists! {!s}'.format(pdbout))

    # Initialise Commander
    CSYMMATCH = CommandManager('csymmatch')
    # Set Command Arguments
    CSYMMATCH.SetArguments('-pdbin-ref',os.path.abspath(refpdb),'-pdbin',os.path.abspath(movpdb),'-pdbout',os.path.abspath(pdbout),'-connectivity-radius',str(conrad))
    # Run!
    CSYMMATCH.Run()

    if not os.path.exists(pdbout):
        raise ExternalProgramError('CSYMMATCH has failed to map {!s} to {!s}.\nERR: {!s}'.format(movpdb,refpdb,CSYMMATCH.err))

    # Check output is valid (raises LigandCheckError if it fails)
    check_pdb_readable(pdbout)

    return pdbout

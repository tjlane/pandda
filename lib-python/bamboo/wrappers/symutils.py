import os, sys, shutil

from bamboo.common.command import CommandManager
from bamboo.rdkit_utils.mol import check_pdb_readable

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
    PDBSET.add_command_line_arguments('XYZIN',os.path.abspath(pdbin),'XYZOUT',os.path.abspath(pdbout))
    # Set inputs
    PDBSET.add_standard_input(['SYMGEN {!s}'.format(sgno),'CELL {!s}'.format(' '.join(map(str,cell))), 'END'])
    # run!
    PDBSET.run()

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
    CSYMMATCH.add_command_line_arguments('-pdbin-ref',os.path.abspath(refpdb),'-pdbin',os.path.abspath(movpdb),'-pdbout',os.path.abspath(pdbout),'-connectivity-radius',str(conrad))
    # run!
    CSYMMATCH.run()

    if not os.path.exists(pdbout):
        raise ExternalProgramError('CSYMMATCH has failed to map {!s} to {!s}.\nERR: {!s}'.format(movpdb,refpdb,CSYMMATCH.err))

    # Check output is valid (raises LigandCheckError if it fails)
    check_pdb_readable(pdbout)

    return pdbout

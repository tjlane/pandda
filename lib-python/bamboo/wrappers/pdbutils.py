import os, sys, shutil

from bamboo.common.command import CommandManager
from bamboo.macro.utils import get_residue_labels
from bamboo.constants import DEFAULT_OUTPUT_CHAIN, DEFAULT_OUTPUT_RESNUM

def isolate_residue(inpdb, outpdb, resname):
    """Extract the residues identified by resname using pdbcur - i.e. 'UNL'"""

    PDBCUR = CommandManager('pdbcur')
    PDBCUR.add_command_line_arguments('XYZIN',inpdb,'XYZOUT',outpdb)
    PDBCUR.add_standard_input(['lvresidue /*/*/({!s})'.format(resname),'END'])
    PDBCUR.run()

    return PDBCUR

def remove_residue(inpdb, outpdb, resname):
    """Delete the residues identified by resname using pdbcur - i.e. 'UNL'"""

    PDBCUR = CommandManager('pdbcur')
    PDBCUR.add_command_line_arguments('XYZIN',inpdb,'XYZOUT',outpdb)
    PDBCUR.add_standard_input(['delresidue /*/*/({!s})'.format(resname),'END'])
    PDBCUR.run()

    return PDBCUR

def change_residue_chain_and_number_to_default(inpdb):
    """Takes the residue in inpdb and moves it to chain X, residue 666 - IN PLACE (keeps original copy as inpdb+'.origfile')"""

    # Get the current labels
    res_id = get_residue_labels(inpdb)
    assert len(res_id)!=0, 'NO RESIDUES IN FILE! {!s}'.format(inpdb)
    assert len(res_id)==1, 'MORE THAN ONE RESIDUE IS PRESENT IN THE FILE! {!s}'.format(inpdb)
    curr_resname, curr_chain, curr_resnum, curr_inscode = res_id[0]
    if not curr_chain: curr_chain='\'\''
    new_chain, new_resnum = DEFAULT_OUTPUT_CHAIN[0], DEFAULT_OUTPUT_RESNUM[0]

    # Move file to prepare for new file
    origfile = inpdb + '.origfile'
    shutil.move(inpdb, origfile)

    PDBSET = CommandManager('pdbset')
    PDBSET.add_command_line_arguments('XYZIN',origfile,'XYZOUT',inpdb)
    PDBSET.add_standard_input(['RENUMBER {!s} CHAIN {!s} TO {!s}'.format(new_resnum, curr_chain, new_chain), 'END'])
    PDBSET.run()

    return PDBSET

def isolate_residue_by_res_id(inpdb, outpdb, chain, resnum, model='*', inscode=''):
    """Isolate the residues identified by residue ids using pdbcur - i.e. '0/A/54.A/'"""

    if inscode:
        selection = '/{!s}/{!s}/{!s}.{!s}'.format(model,chain,resnum,inscode)
    else:
        selection = '/{!s}/{!s}/{!s}'.format(model,chain,resnum)

    PDBCUR = CommandManager('pdbcur')
    PDBCUR.add_command_line_arguments('XYZIN',inpdb,'XYZOUT',outpdb)
    PDBCUR.add_standard_input(['lvresidue {!s}'.format(selection),'END'])
    PDBCUR.run()

    return PDBCUR

def remove_residue_by_res_id(inpdb, outpdb, chain, resnum, model='*', inscode='', removeSolvent=False):
    """Remove the residues identified by res info using pdbcur - i.e. 'UNL'"""

    if inscode:
        selection = '/{!s}/{!s}/{!s}.{!s}'.format(model,chain,resnum,inscode)
    else:
        selection = '/{!s}/{!s}/{!s}'.format(model,chain,resnum)

    std_input = ['delresidue {!s}'.format(selection)]+(removeSolvent)*['delsolvent']+['END']

    PDBCUR = CommandManager('pdbcur')
    PDBCUR.add_command_line_arguments('XYZIN',inpdb,'XYZOUT',outpdb)
    PDBCUR.add_standard_input(std_input)
    PDBCUR.run()

    return PDBCUR

def merge_pdb_files(pdb1, pdb2, pdbout):
    """Merge two PDB Files using CCP4s pdb_merge"""

    # Initialise Commander
    MERGER = CommandManager('pdb_merge')
    # Set command arguments
    MERGER.add_command_line_arguments('xyzin1', pdb1, 'xyzin2', pdb2, 'xyzout', pdbout)
    # Set inputs
    MERGER.add_standard_input('END')
    # run!
    MERGER.run()

    return MERGER

def reset_pdb_file(pdbin, pdbout):
    """Resets the B-Factors in a file and removes anisotropy"""

    if not pdbout.endswith('.pdb'):
        pdbout = pdbout+'.pdb'
    if not os.path.exists(pdbin):
        raise IOError('pdbin does not exist! {!s}'.format(pdbin))
    if os.path.exists(pdbout):
        raise Exception('pdbout already exists! {!s}'.format(pdbout))

    pdbtemp = pdbout.replace('.pdb','.temp.pdb')

    # Initialise Commander
    PDBCUR = CommandManager('pdbcur')
    # Set Command Arguments
    PDBCUR.add_command_line_arguments('XYZIN',pdbin,'XYZOUT',pdbtemp)
    # Set inputs
    PDBCUR.add_standard_input(['NOANISOU','DELSOLVENT','END'])
    # run!
    PDBCUR.run()

    if not os.path.exists(pdbtemp):
        raise ExternalProgramError('PDBCUR has failed to remove anisotropy and delete solvent. {!s}\nOUT: {!s}\nERR: {!s}'.format(pdbin, PDBCUR.out, PDBCUR.err))

    # Initialise Commander
    PDBSET = CommandManager('pdbset')
    # Set Command Arguments
    PDBSET.add_command_line_arguments('XYZIN',pdbtemp,'XYZOUT',pdbout)
    # Set inputs
    PDBSET.add_standard_input(['BFACTOR','END'])
    # run!
    PDBSET.run()

    if not os.path.exists(pdbout):
        raise ExternalProgramError('PDBSET has failed to reset B-factors. {!s}\nOUT: {!s}\nERR: {!s}'.format(pdbtemp, PDBSET.out, PDBSET.err))
    else:
        os.remove(pdbtemp)

    return pdbout

def remove_modres_records(pdb):
    """greps out lines beginning with MODRES from a file"""

    assert os.path.exists(pdb), 'input file must exist'

    origpdb = pdb + '.origfile'
    assert not os.path.exists(origpdb), 'temporary file must not exist'

    shutil.move(pdb, origpdb)

    cmd = "grep -v '^MODRES' {!s} > {!s}".format(origpdb, pdb)
    os.system(cmd)

    if not os.path.exists(pdb):
        raise ExternalProgramError('Failed to remove modres records from {!s}'.format(pdb))

    return

def create_alpha_carbon_backbone(pdbin, pdbout):
    """Takes a pdb files and removes eveything except for the alpha carbons"""

    if not pdbout.endswith('.pdb'):
        pdbout = pdbout+'.pdb'
    if not os.path.exists(pdbin):
        raise IOError('pdbin does not exist! {!s}'.format(pdbin))
    if os.path.exists(pdbout):
        raise Exception('pdbout already exists! {!s}'.format(pdbout))

    # Initialise Commander
    PDBCUR = CommandManager('pdbcur')
    # Set Command Arguments
    PDBCUR.add_command_line_arguments('XYZIN',pdbin,'XYZOUT',pdbout)
    # Set inputs
    PDBCUR.add_standard_input(['lvatom "CA[C]:*"','END'])
    # run!
    PDBCUR.run()

    if not os.path.exists(pdbout):
        raise ExternalProgramError('PDBCUR has failed to create carbon backbone. {!s}\nOUT: {!s}\nERR: {!s}'.format(pdbin, PDBCUR.out, PDBCUR.err))

    return pdbout

def create_cryst_line(pdbin, pdbout, sg, cell):
    """Adds a cryst line to pdbin"""

    if not pdbout.endswith('.pdb'):
        pdbout = pdbout+'.pdb'
    if not os.path.exists(pdbin):
        raise IOError('pdbin does not exist! {!s}'.format(pdbin))
    if os.path.exists(pdbout):
        raise Exception('pdbout already exists! {!s}'.format(pdbout))

    # Initialise Commander
    PDBSET = CommandManager('pdbset')
    # Set Command Arguments
    PDBSET.add_command_line_arguments('XYZIN',os.path.abspath(pdbin),'XYZOUT',os.path.abspath(pdbout))
    # Set Stdin
    PDBSET.add_standard_input(['SPACEGROUP {!s}'.format(sg),'CELL {!s}'.format(' '.join(map(str,cell)))])
    # run!
    PDBSET.run()

    if not os.path.exists(pdbout):
        raise ExternalProgramError('PDBSET has failed to create cryst line for {!s}\nOUT: {!s}\nERR: {!s}'.format(pdbin, PDBSET.out, PDBSET.err))

    return PDBSET

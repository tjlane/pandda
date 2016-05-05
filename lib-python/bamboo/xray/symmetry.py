import os, sys, shutil

from rdkit import Chem

from bamboo.utils.mtz import get_mtz_summary_dict
from bamboo.wrappers.symutils import generate_symmetry_mates
from bamboo.rdkit_utils.contacts.utils import order_structures_by_minimum_distance_to_reference, order_structures_by_number_of_contacts_to_reference

def generate_full_unit_cell(pdbin, pdbout, mtz=None, sgno=None, cell=None):
    """Map pdbin to the whole of the unit cell using symmetry operations"""

    if not os.path.exists(pdbin):
        raise IOError('pdbin does not exist! {!s}'.format(pdbin))
    if not pdbout.endswith('.pdb'):
        pdbout += '.pdb'
    if os.path.exists(pdbout):
        raise Exception('File already exists: {!s}'.format(pdbout))

    if (sgno != None) and (cell != None):
        pass
    elif (mtz != None):
        # Get the spacegroup number from the mtzfile
        mtzsummary = get_mtz_summary_dict(mtz)
        sgno = mtzsummary['spacegroupno']
        cell = mtzsummary['cell']
    else:
        raise Exception('mtz, or sgno AND cell, must be given!')

    # Generate the symmetry equivalent molecules
    pdbout = generate_symmetry_mates(pdbin=pdbin, pdbout=pdbout, sgno=sgno, cell=cell)

    return pdbout

def generate_symmetry_equivalent_ligands(pdbin, pdbouttemplate, mtz=None, sgno=None, cell=None, delete_structures=True):
    """Takes the ligand and a spacegroup number and creates symmetry copies. Returns a list of the new files
        Only designed to be used on ligands - do not use for large structures!"""

    if not os.path.exists(pdbin):
        raise IOError('pdbin does not exist! {!s}'.format(pdbin))
    if pdbouttemplate.endswith('.pdb'):
        pdbouttemplate = pdbouttemplate.replace('.pdb','')
    pdbout_individual = pdbouttemplate+'.{!s}.pdb'

    if (sgno != None) and (cell != None):
        pass
    elif (mtz != None):
        # Get the spacegroup number from the mtzfile
        mtzsummary = get_mtz_summary_dict(mtz)
        sgno = mtzsummary['spacegroupno']
        cell = mtzsummary['cell']
    else:
        raise Exception('SGNO and Cell, or MTZ must be given!')

    # Check that it can be read
    test_mol = Chem.MolFromPDBFile(pdbin)

    if not test_mol:
        raise RDkitReadError('"generate_symmetry_equivalent_ligands": Could not read input molecule {!s}'.format(pdbin))

    # Generate the symmetry equivalent molecules
    pdbtemp = generate_symmetry_mates(pdbin=pdbin, pdbout=pdbouttemplate+'.allsym.pdb', sgno=sgno, cell=cell)

    if not os.path.exists(pdbtemp):
        raise IOError('No Symmetry mates have been generated for {!s}'.format(pdbin))

    # Now use rdkit to split them into separate molecules
    all_mols = Chem.MolFromPDBFile(pdbtemp)

    # This will fail if the model is on an axis of symmetry (molecule maps on top of itself) - return the original model
    if not all_mols:
        # Copy the input model to the output
        out_filename = pdbout_individual.format(0)
        shutil.copy(pdbin, out_filename)
        # Delete the temporary file and return
        if delete_structures:
            os.remove(pdbtemp)
        return [out_filename]

    ind_mols = Chem.GetMolFrags(all_mols, asMols=True)

    # Write them to individual files
    outfiles = []
    for i, mol in enumerate(ind_mols):
        out_filename = pdbout_individual.format(i)
        Chem.MolToPDBFile(mol, out_filename)
        outfiles.append(out_filename)

    # Remove the tempfile
    if delete_structures: os.remove(pdbtemp)

    return outfiles


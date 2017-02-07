
import pandas
import iotbx.pdb

from scitbx.array_family import flex
from giant.maths.geometry import is_within

import iotbx.pdb.hierarchy as iotbx_pdbh

###############################################################################
###                                                                         ###
###                                LABELS                                   ###
###                                                                         ###
###############################################################################

def make_label(obj):
    """Return the relevant label for a supplied hierarchy/atom object"""
    if isinstance(obj, iotbx_pdbh.residue_group):
        return label_from_residue_group(obj)
    elif isinstance(obj, iotbx_pdbh.atom_group):
        return label_from_atom_group(obj)
    elif isinstance(obj, iotbx_pdbh.conformer):
        return label_from_conformer(obj)
    elif isinstance(obj, iotbx_pdbh.residue):
        return label_from_residue(obj)
    elif isinstance(obj, iotbx_pdbh.atom):
        return label_from_atom(obj)
    elif isinstance(obj, iotbx_pdbh.atom_with_labels):
        return label_from_atom_with_labels(obj)
    else:
        raise Exception('Invalid object type provided: {}'.format(type(obj)))

def label_from_residue_group(rg):
    """Return (chain_id, resid)"""
    ch = rg.parent()
    return (ch.id, rg.resid())

def label_from_atom_group(ag):
    """Return (chain_id, resid, altloc)"""
    rg = ag.parent()
    ch = rg.parent()
    return (ch.id, rg.resid(), ag.altloc)

def label_from_conformer(cnf):
    """Return (chain_id, resid, altloc). Must only have one residue."""
    ch = cnf.parent()
    res = cnf.only_residue()
    return (ch.id, res.resid(), cnf.altloc)

def label_from_residue(res):
    """Return (chain_id, resid, altloc)."""
    cnf = res.parent()
    ch = cnf.parent()
    return (ch.id, res.resid(), cnf.altloc)

def label_from_atom(at):
    """Return (chain_id, resid, altloc)"""
    ag = at.parent()
    rg = ag.parent()
    ch = rg.parent()
    return (ch.id, rg.resid(), ag.altloc)

def label_from_atom_with_labels(at):
    """Return (chain_id, resid, altloc)"""
    return (at.chain_id, at.resid(), at.altloc)

###############################################################################
###                                                                         ###
###                         CONVENIENCE FUNCTIONS                           ###
###                                                                         ###
###############################################################################

def get_atom_pairs(residue_1, residue_2, fetch_labels=True):
    atom_pairs = []
    a1 = residue_1.atoms()
    a2 = residue_2.atoms()
    a1_d = a1.build_dict()
    a2_d = a2.build_dict()
    for a_name in a1.extract_name():
        if a_name not in a2_d: continue
        if fetch_labels:
            atom_pairs.append((a1_d[a_name].fetch_labels(),a2_d[a_name].fetch_labels()))
        else:
            atom_pairs.append((a1_d[a_name],a2_d[a_name]))
    return atom_pairs

def calculate_residue_group_occupancy(residue_group):
    """
    Calculate the occupancy of a whole residue, allowing for partial occupancy of partial conformers.
    If conformers are present, the maximum occupancies of each conformer are summed
    """
    rg = residue_group
    if [c.altloc for c in rg.conformers()] == ['']:
        res_occ = max(rg.atoms().extract_occ())
    else:
        res_occ = sum([max(c.atoms().extract_occ()) for c in rg.conformers() if c.altloc])
    if res_occ > 1:
        # Currently errors silently by returning None
        return None
    return res_occ

def calculate_residue_group_rmsd(residue_group_1, residue_group_2):
    """Calculate rmsd between two hierarchically-identical residue groups (only coordinates, occupancies, and b-factors allowed to be different)"""
    rg_1 = residue_group_1; rg_2 = residue_group_2
    ats_1 = rg_1.atoms();   ats_2 = rg_2.atoms()
    # Check that the residues are equivalent
    assert ats_1.size()                       == ats_2.size()
    assert list(ats_1.extract_name()   )      == list(ats_2.extract_name()   )
    assert list(ats_1.extract_element())      == list(ats_2.extract_element())
    assert [a.parent().altloc for a in ats_1] == [a.parent().altloc for a in ats_2]
    rmsd = (ats_1.extract_xyz()-ats_2.extract_xyz()).rms_length()
    return rmsd

def calculate_rmsd(atoms_1, atoms_2, sort=True, truncate_to_common_set=True, remove_H=True):
    """
    Calculate the RMSD between two sets of atoms
    - sort                   - sort atoms prior to rmsd calculation (so that equivalent atoms are compared)
    - truncate_to_common_set - only calculate rmsd over the common_set of atoms, return None if different atoms are provided. Setting this will set sort to true.
    """

    # Need to sort by atoms to check if atoms are the same in both sets
    if truncate_to_common_set:
        sort=True

    # Sort atoms by name if required
    if sort:
        # Extract atom names
        names_1 = list(atoms_1.extract_name())
        names_2 = list(atoms_2.extract_name())
        # Get name overlap between atom sets
        common_atoms = list(set(names_1).intersection(names_2))
        if (not truncate_to_common_set):
            if not (len(common_atoms)==len(atoms_1)==len(atoms_2)):
                return None
        # Get reorderings of atoms
        sort_1 = flex.size_t([names_1.index(an) for an in common_atoms])
        sort_2 = flex.size_t([names_2.index(an) for an in common_atoms])
        # Reorder the atoms
        atoms_1 = atoms_1.select(sort_1)
        atoms_2 = atoms_2.select(sort_2)

    # Remove hydrogen atoms
    if remove_H:
        sel_h = (atoms_1.extract_element() != ' H')
        atoms_1 = atoms_1.select(sel_h)
        atoms_2 = atoms_2.select(sel_h)
    # Check selection working as it should
    if not atoms_1.extract_name().all_eq(atoms_2.extract_name()):
        return None
    # Calculate RMSD and return
    return flex.mean((atoms_1.extract_xyz() - atoms_2.extract_xyz()).dot())**0.5


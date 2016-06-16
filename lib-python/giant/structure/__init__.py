
import pandas
import iotbx.pdb

from scitbx.array_family import flex
from giant.maths.geometry import is_within

import iotbx.pdb.hierarchy as iotbx_pdbh

#---------------------------------------------------->

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

#---------------------------------------------------->

def label_from_residue_group(rg):
    """Return (chain_id, resid)"""
    ch = rg.parent()
    return (ch.id, rg.resid())

def label_from_atom_group(ag):
    """Return (chain_id, resid, altloc)"""
    rg = ag.parent()
    ch = rg.parent()
    return (ch.id, rg.resid(), ag.altloc)

#---------------------------------------------------->

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

#---------------------------------------------------->

def label_from_atom(at):
    """Return (chain_id, resid, altloc)"""
    ag = at.parent()
    rg = ag.parent()
    ch = rg.parent()
    return (ch.id, rg.resid(), ag.altloc)

def label_from_atom_with_labels(at):
    """Return (chain_id, resid, altloc)"""
    return (at.chain_id, at.resid(), at.altloc)

#---------------------------------------------------->

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

def calculate_residue_group_bfactor_ratio(residue_group, hierarchy, data_table=None, rg_label=None, column_suffix=''):
    """Calculate bfactor quality metrics of the residue to surrounding residues"""

    rg_sel = residue_group
    # Set defaults
    if rg_label is None:    rg_label = (rg_sel.unique_resnames()[0]+'-'+rg_sel.parent().id+'-'+rg_sel.resseq+rg_sel.icode).replace(' ','')
    if data_table is None:  data_table = pandas.DataFrame(index=[rg_label], column=[])
    # Check validity
    if len(rg_sel.unique_resnames()) != 1: raise Exception(rg_label+': More than one residue name associated with residue group -- cannot process')

    # Extract rg_sel objects
    rg_sel_ags    = rg_sel.atom_groups()
    rg_sel_atoms  = rg_sel.atoms()
    rg_sel_coords = rg_sel_atoms.extract_xyz()
    # Select nearby atom_group for scoring
    near_ags = [ag for ag in hierarchy.atom_groups() if (is_within(4, rg_sel_coords, ag.atoms().extract_xyz()) and (ag not in rg_sel_ags))]
    if near_ags:
        # Extract atoms from nearby groups
        near_ats = iotbx.pdb.hierarchy.af_shared_atom()
        [near_ats.extend(ag.detached_copy().atoms()) for ag in near_ags]
        # Calculate B-factors of the residue
        res_mean_b = flex.mean_weighted(rg_sel_atoms.extract_b(), rg_sel_atoms.extract_occ())
        data_table.set_value(   index = rg_label,
                                col   = 'Average B-factor (Residue)'+column_suffix,
                                value = res_mean_b )
        # Calculate B-factors of the surrounding atoms
        sch_mean_b = flex.mean_weighted(near_ats.extract_b(), near_ats.extract_occ())
        data_table.set_value(   index = rg_label,
                                col   = 'Average B-factor (Surroundings)'+column_suffix,
                                value = sch_mean_b )
        # Store the ratio of the b-factors
        data_table.set_value(   index = rg_label,
                                col   = 'Surroundings B-factor Ratio'+column_suffix,
                                value = res_mean_b/sch_mean_b )

    return data_table


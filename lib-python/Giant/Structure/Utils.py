import os, sys, copy
import numpy
import iotbx.pdb
from scitbx.array_family import flex

def normalise_occupancies(hierarchy, exclude_conformers=[], max_occ=1.0, min_occ=0.0, in_place=False):
    """Normalise the occupancies of a hierarchy so that the occupancies for a residue sum to 1.0"""

    # Edit original copy or create new one?
    if (not in_place): hierarchy = hierarchy.deep_copy()
    # Iterate through the output structure, and normalise the occupancies if necessary
    for rg in hierarchy.residue_groups():
        # No conformers - don't normalise
        if not rg.have_conformers(): continue

        # Get the groups to change and the groups to keep constant
        ag_chnge = [ag for ag in rg.atom_groups() if ag.altloc and (ag.altloc not in exclude_conformers)]
        if not ag_chnge: continue
        ag_const = [ag for ag in rg.atom_groups() if ag.altloc and (ag.altloc in exclude_conformers)]
        # Calculate the total occupancy of the groups
        occ_chnge = sum([max(ag.atoms().extract_occ()) for ag in ag_chnge])
        occ_const = sum([max(ag.atoms().extract_occ()) for ag in ag_const])

        if (occ_const+occ_chnge < min_occ):
            # Occupancy too low - normalise to minimum value
            new_occ_chnge = min_occ - occ_const
        elif (occ_const+occ_chnge > max_occ):
            # Occupancy too high - normalise to maximum value
            new_occ_chnge = max_occ - occ_const
        else:
            # Not selected for normalising
            continue

        # Normalise residue occupancies
        assert new_occ_chnge > 0.0, 'Occupancy of INVARIABLE AGs is more than 1.0: {}'.format(occ_const)
        [ag.atoms().set_occ(ag.atoms().extract_occ()*new_occ_chnge/occ_chnge) for ag in ag_chnge]

    return hierarchy

def old_normalise_occupancies(hierarchy, in_place=False):
    """Normalise the occupancies of a hierarchy so that the occupancies for a residue sum to 1.0"""

    # Edit original copy or create new one?
    if (not in_place): hierarchy = hierarchy.deep_copy()
    # Iterate through the output structure, and normalise the occupancies if necessary
    for rg in hierarchy.residue_groups():
        # Calculate the total occupancy for this residue group
        occ_total = sum([max(ag.atoms().extract_occ()) for ag in rg.atom_groups()])
        # If occupancy is greater than 1, normalise by this total
        if occ_total > 1.0:
            [ag.atoms().set_occ(ag.atoms().extract_occ()/occ_total) for ag in rg.atom_groups()]
            assert sum([max(ag.atoms().extract_occ()) for ag in rg.atom_groups()]) <= 1.1, [max(ag.atoms().extract_occ()) for ag in rg.atom_groups()]
    return hierarchy

def set_conformer_occupancy(hierarchy, conf_id, conf_occ, in_place=False):
    """Normalise the occupancies of a hierarchy so that the occupancies for a residue sum to 1.0"""

    # Edit original copy or create new one?
    if (not in_place): hierarchy = hierarchy.deep_copy()
    # Iterate through the output structure, and set the occupancies where selected
    for rg in hierarchy.residue_groups():
        # Calculate the total occupancy for this residue group
        for ag in rg.atom_groups():
            if ag.altloc == conf_id: ag.atoms().set_occ(flex.double([conf_occ]*len(ag.atoms())))
    return hierarchy

def find_unused_chain_ids(hierarchies):
    unused_chain_ids = iotbx.pdb.systematic_chain_ids()
    # Extract current chain_ids
    current_chain_ids=[]; [current_chain_ids.extend([ch.id for ch in h.chains()]) for h in hierarchies]
    current_chain_ids=list(set(current_chain_ids))
    # Remove present chains
    for curr_id in current_chain_ids:   unused_chain_ids.remove(curr_id)
    return unused_chain_ids

def resolve_residue_id_clashes(ref_hierarchy, mov_hierarchy, in_place=False):
    """Move residues in mov_hierarchiy to new chains if they have the same resid as a residue in ref_hierarchy but different resnames"""

    # Edit original copy or create new one?
    if (not in_place):
        ref_hierarchy = ref_hierarchy.deep_copy()
        mov_hierarchy = mov_hierarchy.deep_copy()

    # New chain to add clashing groups to
    new_chain = None
    new_chain_ids = find_unused_chain_ids(hierarchies=[ref_hierarchy, mov_hierarchy])

    # Find the residues with clashing resids
    residues_to_update = []
    for rg_mov in mov_hierarchy.residue_groups():
        # Extract label for this residue
        resid = rg_mov.resid(); chainid = rg_mov.parent().id
        rg_ref = [rg for rg in ref_hierarchy.residue_groups() if rg.resid()==resid and rg.parent().id==chainid]
        # MORE THAN ONE MATCHING RESIDUE -- ERROR?
        if len(rg_ref) > 1:     raise Exception('MORE THAN ONE MATCHING')
        # PRESENT ONLY IN MOVING
        elif len(rg_ref) == 0:  continue
        # PRESENT IN BOTH
        rg_ref = rg_ref[0]
        # Check to see if the residue is the same type as in the reference structure
        if map(str.strip, rg_ref.unique_resnames()) == map(str.strip,rg_mov.unique_resnames()):
            # Same residue -- that's fine
            continue
        else:
            print 'DIFFERENT RESIDUES WITH SAME ID - CHANGING CHAINS: {} != {}'.format(list(rg_ref.unique_resnames()), list(rg_mov.unique_resnames()))
            residues_to_update.append(rg_mov)

    # Go through and transfer the residue groups to new chains
    for rg_mov in residues_to_update:
        # Remove from old chain
        old_chain = rg_mov.parent()
        old_chain.remove_residue_group(rg_mov)
        # See if there is a residue with this id already present in the new chain
        if (not new_chain) or (rg_mov.resid() in new_chain.get_residue_ids()):
            new_chain = iotbx.pdb.hierarchy.chain(id=new_chain_ids.pop(0))
            # Add new chain to the same model as the old chain
            old_chain.parent().append_chain(new_chain)
        # Add to the new chain
        new_chain.append_residue_group(rg_mov)

    return mov_hierarchy

def find_next_conformer_idx(hierarchy, all_ids=iotbx.pdb.systematic_chain_ids()):
    # Current altlocs in the ref structure
    current_conf_ids = sorted(list(hierarchy.altloc_indices()))
    # Find the next one that doesn't exist
    for i_id, new_conf_id in enumerate(all_ids):
        if new_conf_id not in current_conf_ids:
            return i_id

def merge_hierarchies(ref_hierarchy, mov_hierarchy):
    """Merge the hierarchies - creating new altlocs where necessary"""

    # Create new hierarchies to prevent editing the originals
    new_ref = ref_hierarchy.deep_copy()
    new_mov = mov_hierarchy.deep_copy()

    ######################################################################
    # FINDING NEW CONFORMER IDS
    ######################################################################

    # Find the next unused atlloc index to be used for the
    all_conformer_ids = iotbx.pdb.systematic_chain_ids()
    new_conformer_for_ref_idx = find_next_conformer_idx(hierarchy=ref_hierarchy, all_ids=all_conformer_ids)
    new_conformer_for_ref = all_conformer_ids[new_conformer_for_ref_idx]
    new_conformer_for_mov_idx = new_conformer_for_ref_idx + 1
    new_conformer_for_mov = all_conformer_ids[new_conformer_for_mov_idx]

    ######################################################################
    # IDENTIFYING WHICH RESIDUES ARE DIFFERENT BETWEEN MOVING AND REFERENCE STRUCTURES
    ######################################################################

    # Iterate through the mov conformation and find which residues do NOT need merging (have not changed)
    ref_only = []; ref_mov_diff = []; ref_mov_same = []
    # Iterate through residues
    for rg_ref in new_ref.residue_groups():
        # Check that the residue groups only have one residue name in each
        assert len(rg_ref.unique_resnames()) == 1
        # Extract label for this residue
        resid = rg_ref.resid()
        chainid = rg_ref.parent().id
        # Extract same rg for mov
        rg_mov = [rg for rg in new_mov.residue_groups() if rg.resid()==resid and rg.parent().id==chainid]

        # MORE THAN ONE MATCHING RESIDUE -- ERROR?
        if len(rg_mov) > 1:     raise Exception('MORE THAN ONE MATCHING')
        # PRESENT ONLY IN REFERENCE
        elif len(rg_mov) == 0:  ref_only.append((chainid, resid))
        # PRESENT IN BOTH
        elif len(rg_mov) == 1:
            # Check to see if the residue has moved
            rg_mov = rg_mov[0]
            # Check that the residue name is the same
            assert map(str.strip, rg_ref.unique_resnames()) == map(str.strip,rg_mov.unique_resnames()), 'Cannot Merge - Different Residues!: {} != {}'.format(list(rg_ref.unique_resnames()), list(rg_mov.unique_resnames()))
            # Extract atoms and check to see if the same
            rg_ref_ats = rg_ref.atoms()
            rg_mov_ats = rg_mov.atoms()
            # Check to see if the same length
            if len(rg_ref_ats) != len(rg_mov_ats):
                # There's been a change in the structure
                ref_mov_diff.append((chainid, resid))
            else:
                # Extract coordinates for the atoms
                rg_ref_xyz = rg_ref_ats.extract_xyz()
                rg_mov_xyz = rg_mov_ats.extract_xyz()
                # Calculate the rmsd of the atoms (ASSUMING IN THE SAME ORDER)
                rmsd = (rg_ref_xyz - rg_mov_xyz).norm()
                # Check to see if the atoms have moved
                if rmsd != 0.0:
                    # There's been a change in the structure
                    ref_mov_diff.append((chainid, resid))
                else:
                    # Residues haven't moved - no need to transfer this residue
                    ref_mov_same.append((chainid, resid))
    print '=====================>>>'
    print '{!s} RESIDUES CONSERVED (SAME IN REFERENCE AND MOVING)'.format(len(ref_mov_same))
    print '{!s} RESIDUES MOVED (DIFFERENT IN REFERENCE AND MOVING)'.format(len(ref_mov_diff))
    print '{!s} RESIDUES UNIQUE TO REFERENCE CONFORMATION'.format(len(ref_only))

    ######################################################################
    # UPDATING CONFORMER LABELS FOR REFERENCE CONFORMERS
    ######################################################################

    # Record the number of conformers introduced into the ref structure
    conf_introduced = 0
    # Iterate through and update ref conformers
    for rg_ref in new_ref.residue_groups():
        # Extract label for this residue
        resid = rg_ref.resid()
        chainid = rg_ref.parent().id
        # Major and mov are the same - Don't need to change the altloc
        if (chainid, resid) in ref_mov_same:    continue
        # Major only or difference - Change the altloc for this residue to the defaults for the ref (if no conformers present)
        elif not rg_ref.have_conformers():
            assert len(rg_ref.atom_groups()) == 1
            rg_ref.atom_groups()[0].altloc = new_conformer_for_ref
            conf_introduced += 1
        # Conformers already present - don't need to do anything for the ref conformation
        else:
#            print 'CONFORMERS PRESENT (NOT CHANGING): {!s}'.format((chainid, resid))
#            for ag in rg_ref.atom_groups():
#                print ag.altloc if ag.altloc else '-', '->', ag.altloc if ag.altloc else '-'
            pass

    print '=====================>>>'
    print '{!s} CONFORMER IDS UPDATED IN REFERENCE STRUCTURE'.format(conf_introduced)

    ######################################################################
    # UPDATING CONFORMER LABELS FOR MOVING CONFORMERS
    ######################################################################

    # If no conformer, change to default
    # If conformers, change to default + e.g. 'A'

    # Record the number of conformers introduced/incremented into the mov structure
    conf_introduced = 0; conf_incremented = 0
    # Iterate through and update mov conformers
    for rg_mov in new_mov.residue_groups():
        # Extract label for this residue
        resid = rg_mov.resid()
        chainid = rg_mov.parent().id
        # If residue has not changed - skip as will not be transferred anyway
        if (chainid, resid) in ref_mov_same:    continue
        # No conformers - set to the default
        if not rg_mov.have_conformers():
            assert len(rg_mov.atom_groups()) == 1
            rg_mov.atom_groups()[0].altloc = new_conformer_for_mov
            conf_introduced += 1
        # Conformers already present - increment the conformer id by 1
        else:
#            print 'CONFORMERS PRESENT: {!s}'.format((chainid, resid))
            for ag_min in rg_mov.atom_groups():
                conf_id = ag_min.altloc
                # No conformer yet - set to defaults for mov
                if conf_id == '':   new_conf_id = new_conformer_for_mov
                # Existing conformers - increment conformer letter
                else:
                    new_conf_idx = new_conformer_for_mov_idx + all_conformer_ids.index(conf_id) + 1
                    new_conf_id  = all_conformer_ids[new_conf_idx]
                # Update the conformation
#                print conf_id if conf_id else '-', '->', new_conf_id
                ag_min.altloc = new_conf_id
            conf_incremented += 1

    print '=====================>>>'
    print '{!s} CONFORMER IDS CREATED IN MOVING STRUCTURE'.format(conf_introduced)
    print '{!s} CONFORMER IDS UPDATED IN MOVING STRUCTURE'.format(conf_incremented)

    ######################################################################
    # TRANSFERRING RESIDUES FROM MOVING CONFORMATION
    ######################################################################

    # Create another copy of the ref conformation
    final_struct = new_ref.deep_copy()

    # Iterate through chains
    for ch_min in new_mov.chains():

        # Find the matching chain in the ref conformation
        ch_fin = [c for c in final_struct.chains() if c.id==ch_min.id]
        # MORE THAN ONE MATCHING RESIDUE -- ERROR?
        if len(ch_fin) > 1:     raise Exception('MORE THAN ONE MATCHING')
        # If not a matching chain, add to a new chain
        elif len(ch_fin) == 0:  final_struct.models()[0].append_chain(ch_min.detached_copy())
        # One matching chain, so transfer the residue groups
        elif len(ch_fin) == 1:
            ch_fin = ch_fin[0]
            # Iterate through residue groups
            for rg_mov in ch_min.residue_groups():
                # Extract label for this residue
                resid = rg_mov.resid()
                chainid = rg_mov.parent().id
                # If residue has not changed - do not transfer
                if (chainid, resid) in ref_mov_same:    continue
                # Find the matching residue group in the ref conformation
                rg_fin = [rg for rg in ch_fin.residue_groups() if resid==rg.resid()]
                # MORE THAN ONE MATCHING RESIDUE -- ERROR?
                if len(rg_fin) > 1:     raise Exception('MORE THAN ONE MATCHING')
                # If not a matching residue group, add to a new residue group
                elif len(rg_fin)==0:    ch_fin.append_residue_group(rg_mov.detached_copy())
                # One matching residue, so transfer the atom groups (don't have to worry about clashes as altlocs have already been changed)
                elif len(rg_fin)==1:
                    rg_fin = rg_fin[0]
                    # Transfer the atom groups
                    for min_ag in rg_mov.atom_groups(): rg_fin.append_atom_group(min_ag.detached_copy())

    return final_struct


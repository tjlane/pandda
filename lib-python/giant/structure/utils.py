import os, sys, copy
import numpy
import iotbx.pdb
from scitbx.array_family import flex

def find_unused_chain_ids(hierarchies):
    unused_chain_ids = iotbx.pdb.systematic_chain_ids()
    # Extract current chain_ids
    current_chain_ids=[]; [current_chain_ids.extend([ch.id for ch in h.chains()]) for h in hierarchies]
    current_chain_ids=list(set(current_chain_ids))
    # Remove present chains
    for curr_id in current_chain_ids:   unused_chain_ids.remove(curr_id)
    return unused_chain_ids

def find_next_conformer_idx(hierarchy, all_ids=iotbx.pdb.systematic_chain_ids()):
    # Current altlocs in the ref structure
    current_conf_ids = sorted(list(hierarchy.altloc_indices()))
    # Find the next one that doesn't exist
    for i_id, new_conf_id in enumerate(all_ids):
        if new_conf_id not in current_conf_ids:
            return i_id

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

def compare_hierarchies(reference_hierarchy, query_hierarchy):
    """Find which residues are the same, different, or unique to each of the hierarchies"""

    # Can't be bothered to type those names all the time
    ref_h = reference_hierarchy; qry_h = query_hierarchy
    # Iterate through the reference hierarchy and find which residues have changed/are unique/are the same
    ref_only = []; ref_diff = []; ref_same = []

    # Iterate through residues
    for rg_ref in ref_h.residue_groups():
        # Check that the residue groups only have one residue name in each
        assert len(rg_ref.unique_resnames()) == 1, "Can't handle residues with more than one residue name"
        # Extract label for this residue
        resid = rg_ref.resid(); chainid = rg_ref.parent().id
        # Extract same rg for qry
        rg_qry = [rg for rg in qry_h.residue_groups() if rg.resid()==resid and rg.parent().id==chainid]

        # MORE THAN ONE MATCHING RESIDUE -- ERROR?
        if len(rg_qry) > 1:
            raise Exception('MORE THAN ONE MATCHING')
        # PRESENT ONLY IN REFERENCE
        elif len(rg_qry) == 0:
            ref_only.append((chainid, resid))
        # PRESENT IN BOTH
        elif len(rg_qry) == 1:
            # Check to see if the residue has moved
            rg_qry = rg_qry[0]
            # Check that the residue name is the same
            assert map(str.strip, rg_ref.unique_resnames()) == map(str.strip,rg_qry.unique_resnames()), 'Cannot Merge - Different Residues!: {} != {}'.format(list(rg_ref.unique_resnames()), list(rg_qry.unique_resnames()))
            # Extract atoms and check to see if the same
            rg_ref_ats = rg_ref.atoms()
            rg_qry_ats = rg_qry.atoms()
            # Check to see if the same length
            if len(rg_ref_ats) != len(rg_qry_ats):
                # There's been a change in the structure
                ref_diff.append((chainid, resid))
            else:
                # Extract coordinates for the atoms
                rg_ref_xyz = rg_ref_ats.extract_xyz()
                rg_qry_xyz = rg_qry_ats.extract_xyz()
                # Calculate the rmsd of the atoms (ASSUMING IN THE SAME ORDER)
                rmsd = (rg_ref_xyz - rg_qry_xyz).norm()
                # Check to see if the atoms have moved
                if rmsd != 0.0:
                    # There's been a change in the structure
                    ref_diff.append((chainid, resid))
                else:
                    # Same between the two structures
                    ref_same.append((chainid, resid))

    return ref_only, ref_diff, ref_same

def create_pure_alt_conf_from_proper_alt_conf(residue_group):
    """Take 'main conf' atoms and add to 'pure conformer' atom_groups"""

    main_ags = [ag for ag in residue_group.atom_groups() if ag.altloc=='']
    conf_ags = [ag for ag in residue_group.atom_groups() if ag.altloc!='']
    assert len(main_ags)==1, "Must be one atom_group in residue_group with altloc==''"
    assert len(conf_ags)!=0, "Must be at least one alternate conformer present"

    # Get a blank copy of the residue_group
    blank_rg = residue_group.detached_copy()
    for ag in blank_rg.atom_groups(): blank_rg.remove_atom_group(ag)
    assert not blank_rg.atom_groups(), blank_rg.atom_groups()

    # Atoms to be transferred to the other ags
    main_ag = main_ags[0]
    # Iterate through and add atoms from main_conf to each of the other atom_groups
    for conf_ag in conf_ags:
        new_conf_ag = conf_ag.detached_copy()
        new_main_ag = main_ag.detached_copy()
        # Set the occupancy of the main_conf atoms to be transferred
        max_occ = max(new_conf_ag.atoms().extract_occ())
        new_main_ag.atoms().set_occ(flex.double([max_occ]*new_main_ag.atoms().size()))
        # Copy atoms from the main conf to the alt conf
        for i_at, at in enumerate(new_main_ag.atoms()):
            new_conf_ag.insert_atom(i_at, at.detached_copy())
        # Add pure conf ag to the blank rg
        blank_rg.append_atom_group(new_conf_ag)
    return blank_rg

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

    ref_only, ref_diff, ref_same = compare_hierarchies( reference_hierarchy = new_ref,
                                                        query_hierarchy     = new_mov   )
    print '=====================>>>'
    print '{!s} RESIDUES CONSERVED (SAME IN REFERENCE AND MOVING)'.format(len(ref_same))
    print '{!s} RESIDUES MOVED (DIFFERENT IN REFERENCE AND MOVING)'.format(len(ref_diff))
    print '{!s} RESIDUES UNIQUE TO REFERENCE CONFORMATION'.format(len(ref_only))

    ######################################################################
    # UPDATING CONFORMER LABELS FOR REFERENCE CONFORMERS
    ######################################################################

    # Record the number of conformers introduced into the ref structure
    conf_introduced = 0; conf_distributed = 0
    # Iterate through and update ref conformers
    for rg_ref in new_ref.residue_groups():
        # Extract label for this residue
        resid = rg_ref.resid(); chainid = rg_ref.parent().id
        # Major and mov are the same - Don't need to change the altloc
        if (chainid, resid) in ref_same:    continue
        # Major only or difference - Change the altloc for this residue to the defaults for the ref (if no conformers present)
        elif not rg_ref.have_conformers():
            assert len(rg_ref.atom_groups()) == 1
            rg_ref.atom_groups()[0].altloc = new_conformer_for_ref
            conf_introduced += 1
        # Conformers already present - make sure they are pure not proper
        elif '' in [ag.altloc for ag in rg_ref.atom_groups()]:
            new_rg = create_pure_alt_conf_from_proper_alt_conf(residue_group=rg_ref)
            # Find the rg's place in the parent chain
            ch = rg_ref.parent()
            rg_idx = ch.find_residue_group_index(rg_ref)
            # Remove the old rg and insert the new rg in it's place
            ch.remove_residue_group(rg_idx)
            ch.insert_residue_group(rg_idx, new_rg)
            conf_distributed += 1

    print '=====================>>>'
    print '{!s} CONFORMER IDS CREATED IN REFERENCE STRUCTURE'.format(conf_introduced)
    print '{!s} CONFORMER IDS UPDATED IN REFERENCE STRUCTURE'.format(conf_distributed)

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
        resid = rg_mov.resid(); chainid = rg_mov.parent().id
        # If residue has not changed - skip as will not be transferred anyway
        if (chainid, resid) in ref_same:    continue
        # Transfer but no conformers - set to the default
        elif not rg_mov.have_conformers():
            assert len(rg_mov.atom_groups()) == 1
            rg_mov.atom_groups()[0].altloc = new_conformer_for_mov
            conf_introduced += 1
        # Conformers already present
        else:
            # Increment conformer letter
            for ag_min in rg_mov.atom_groups():
                if ag_min.altloc:
                    new_conf_idx = new_conformer_for_mov_idx + all_conformer_ids.index(ag_min.altloc) + 1
                    ag_min.altloc = all_conformer_ids[new_conf_idx]
            # Make sure they are pure not proper
            if '' in [ag.altloc for ag in rg_mov.atom_groups()]:
                new_rg = create_pure_alt_conf_from_proper_alt_conf(residue_group=rg_mov)
                # Find the rg's place in the parent chain
                ch = rg_mov.parent()
                rg_idx = ch.find_residue_group_index(rg_mov)
                # Remove the old rg and insert the new rg
                ch.remove_residue_group(rg_idx)
                ch.insert_residue_group(rg_idx, new_rg)
            conf_incremented += 1

    print '=====================>>>'
    print '{!s} CONFORMER IDS CREATED IN MOVING STRUCTURE'.format(conf_introduced)
    print '{!s} CONFORMER IDS UPDATED IN MOVING STRUCTURE'.format(conf_incremented)

    ######################################################################
    # TRANSFERRING RESIDUES FROM MOVING CONFORMATION
    ######################################################################

    # Create another copy of the ref conformation
    final_struct = new_ref.deep_copy()

    ref_chain_ids = [c.id for c in new_ref.chains()]

    # Iterate through chains
    for ch_min in new_mov.chains():

        # If not a matching chain, add to a new chain
        if ch_min.id not in ref_chain_ids:
            final_struct.models()[0].append_chain(ch_min.detached_copy())
        # Iterate through residue groups
        else:
            for rg_mov in ch_min.residue_groups():
                # Extract label for this residue
                resid = rg_mov.resid()
                chainid = rg_mov.parent().id
                # If residue has not changed - do not transfer
                if (chainid, resid) in ref_same:    continue
                # Find the matching residue group in the ref conformation
                rg_fin = [rg for rg in final_struct.residue_groups() if (resid==rg.resid() and chainid==rg.parent().id)]
                # MORE THAN ONE MATCHING RESIDUE -- ERROR?
                if len(rg_fin) > 1:     raise Exception('MORE THAN ONE MATCHING')
                # If not a matching residue group, add to a new residue group
                elif len(rg_fin)==0:
                    # Get the class name of the residue so that we can identify an appropriate chain object
                    res_class = [iotbx.pdb.common_residue_names_get_class(n) for n in rg_mov.unique_resnames()][0]
                    # Get the chains that already contain a residue with the same class
                    ch_fin = [c for c in final_struct.chains() if (res_class in c.get_residue_names_and_classes()[1].keys())]
                    # If can't find one, then just select all
                    if not ch_fin: ch_fin = final_struct.chains()
                    # Append the residue group to the first identified chain
                    ch_fin[0].append_residue_group(rg_mov.detached_copy())
                # One matching residue, so transfer the atom groups (don't have to worry about clashes as altlocs have already been changed)
                elif len(rg_fin)==1:
                    rg_fin = rg_fin[0]
                    # Transfer the atom groups
                    for min_ag in rg_mov.atom_groups(): rg_fin.append_atom_group(min_ag.detached_copy())

    return final_struct


import os, sys, copy

import libtbx.phil

import numpy

import iotbx.pdb

from scitbx.array_family import flex

############################################################################

systematic_letters = iotbx.pdb.systematic_chain_ids()

############################################################################

master_phil = libtbx.phil.parse("""
major = None
    .help = 'The major conformation of the protein (normally the unbound or reference structure)'
    .type = str

minor = None
    .help = 'The minor conformation of the protein (normally the bound or "interesting" structure)'
    .type = str

output = None
    .help = 'output pdb file'
    .type = str

verbose = False
    .type = bool
""")

############################################################################

def run(params):

    ######################################################################
    print '===========================================>>>'
    print 'READING INPUT FILES'
    ######################################################################

    # Read in the ligand file and set each residue to the requested conformer
    maj_obj = iotbx.pdb.hierarchy.input(params.major)
    min_obj = iotbx.pdb.hierarchy.input(params.minor)

    ######################################################################
    print '===========================================>>>'
    print 'AUTO-GENERATING PARAMETERS'
    ######################################################################

    # Find the next conformer that's not used in the major structure
    break_iter = False
    # Current altlocs in the major structure
    current_conf_ids = list(maj_obj.hierarchy.altloc_indices())
    # Iterate through and select new conformer ids
    for i_id, new_conf_id in enumerate(systematic_letters):

        # STEP 2 (LOOK BELOW FOR STEP 1)
        # Find new id to re-assign to conformers in the minor structure
        # (Once new id has been selected to assign to major structure residues)
        if break_iter:
            params.__inject__('new_conformer_for_minor', new_conf_id)
            params.__inject__('new_conformer_for_minor_idx', i_id)
            print 'Conformers in the minor structure will be re-assigned starting from CONFORMER: {!s}'.format(params.new_conformer_for_minor)
            break

        # STEP 1
        # Find new id to assign to unassigned conformers in the major structure
        if new_conf_id not in current_conf_ids:
            params.__inject__('new_conformer_for_major', new_conf_id)
            print 'Unassigned conformers in the major structure will be assigned to CONFORMER: {!s}'.format(params.new_conformer_for_major)
            break_iter = True

    ######################################################################
    print '===========================================>>>'
    print 'VALIDATING PARAMETERS'
    ######################################################################

    # Check that ... something
    maj_obj.hierarchy.only_model()
    min_obj.hierarchy.only_model()

    ######################################################################
    print '===========================================>>>'
    print 'PREPARING THE STRUCTURES FOR MERGING'
    ######################################################################

    # Create a new copy of the structures
    new_major = maj_obj.hierarchy.deep_copy()
    new_minor = min_obj.hierarchy.deep_copy()

    # Iterate through the minor conformation and find which residues do NOT need merging (have not changed)
    major_only = []
    major_minor_diff = []
    major_minor_same = []
    # Iterate through residues
    for rg_maj in new_major.residue_groups():
        # Extract label for this residue
        resid = rg_maj.resid()
        chainid = rg_maj.parent().id
        # Extract same rg for minor
        rg_min = [rg for rg in new_minor.residue_groups() if rg.resid()==resid and rg.parent().id==chainid]

        # MORE THAN ONE MATCHING RESIDUE -- ERROR?
        if len(rg_min) > 1:
            raise Exception('MORE THAN ONE MATCHING')
        # PRESENT ONLY IN MAJOR
        elif len(rg_min) == 0:
            # No need to worry about merging
            major_only.append((chainid, resid))
        # PRESENT IN BOTH
        elif len(rg_min) == 1:
            # Check to see if the residue has moved
            rg_min = rg_min[0]
            # Extract atoms and check to see if the same
            rg_maj_ats = rg_maj.atoms()
            rg_min_ats = rg_min.atoms()
            # Check to see if the same length
            if len(rg_maj_ats) != len(rg_min_ats):
                # There's been a change in the structure
                major_minor_diff.append((chainid, resid))
            else:
                # Extract corrdinates for the atoms
                rg_maj_xyz = rg_maj_ats.extract_xyz()
                rg_min_xyz = rg_min_ats.extract_xyz()
                # Calculate the rmsd of the atoms
                rmsd = (rg_maj_xyz - rg_min_xyz).norm()
                # Check to see if the atoms have moved
                if rmsd != 0.0:
                    # There's been a change in the structure
                    major_minor_diff.append((chainid, resid))
                else:
                    # Residues haven't moved - no need to transfer this residue
                    major_minor_same.append((chainid, resid))
    print '=====================>>>'
    print '{!s} RESIDUES CONSERVED'.format(len(major_minor_same))
    print '{!s} RESIDUES MOVED'.format(len(major_minor_diff))
    print '{!s} RESIDUES UNIQUE TO MAJOR CONFORMATION'.format(len(major_only))

    ######################################################################
    print '===========================================>>>'
    print 'UPDATING CONFORMER LABELS FOR MAJOR CONFORMERS'
    ######################################################################

    # Record the number of conformers introduced into the major structure
    conf_introduced = 0
    # Iterate through and update major conformers
    for rg_maj in new_major.residue_groups():
        # Extract label for this residue
        resid = rg_maj.resid()
        chainid = rg_maj.parent().id
        # Major and minor are the same - Don't need to change the altloc
        if (chainid, resid) in major_minor_same:
            continue
        # Major only or difference - Change the altloc for this residue to the defaults for the major (if no conformers present)
        elif not rg_maj.have_conformers():
            assert len(rg_maj.atom_groups()) == 1
            rg_maj.atom_groups()[0].altloc = params.new_conformer_for_major
            conf_introduced += 1
        # Conformers already present - don't need to do anything for the major conformation
        else:
            print 'CONFORMERS PRESENT: {!s}'.format((chainid, resid))
            for ag in rg_maj.atom_groups():
                print ag.altloc if ag.altloc else '-', '->', ag.altloc if ag.altloc else '-'
            pass

    print '=====================>>>'
    print '{!s} CONFORMER IDS UPDATED IN MAJOR STRUCTURE'.format(conf_introduced)

    ######################################################################
    print '===========================================>>>'
    print 'UPDATING CONFORMER LABELS FOR MINOR CONFORMERS'
    ######################################################################

    # If no conformer, change to default
    # If conformers, change to default + e.g. 'A'

    # Record the number of conformers introduced into the minor structure
    conf_introduced = 0
    # Record the number of conformers incremented in the minor structure
    conf_incremented = 0
    # Iterate through and update minor conformers
    for rg_min in new_minor.residue_groups():
        # Extract label for this residue
        resid = rg_min.resid()
        chainid = rg_min.parent().id
        # If residue has not changed - skip as will not be transferred anyway
        if (chainid, resid) in major_minor_same:
            continue
        # No conformers - set to the default
        if not rg_min.have_conformers():
            assert len(rg_min.atom_groups()) == 1
            rg_min.atom_groups()[0].altloc = params.new_conformer_for_minor
            conf_introduced += 1
        # Conformers already present - increment the conformer id by 1
        else:
            print 'CONFORMERS PRESENT: {!s}'.format((chainid, resid))
            for ag_min in rg_min.atom_groups():
                conf_id = ag_min.altloc
                # No conformer yet - set to defaults for minor
                if conf_id == '':
                    new_conf_id = params.new_conformer_for_minor
                # Existing conformers - increment conformer letter
                else:
                    new_conf_idx = params.new_conformer_for_minor_idx + systematic_letters.index(conf_id) + 1
                    new_conf_id  = systematic_letters[new_conf_idx]
                # Update the conformation
                print conf_id if conf_id else '-', '->', new_conf_id
                ag_min.altloc = new_conf_id
            conf_incremented += 1

    print '=====================>>>'
    print '{!s} CONFORMER IDS CREATED IN MINOR STRUCTURE'.format(conf_introduced)
    print '{!s} CONFORMER IDS UPDATED IN MINOR STRUCTURE'.format(conf_incremented)

    ######################################################################
    print '===========================================>>>'
    print 'TRANSFERRING RESIDUES FROM MINOR CONFORMATION'
    ######################################################################

    # Create another copy of the major conformation
    final_struct = new_major.deep_copy()

    # Iterate through chains
    for ch_min in new_minor.chains():

        # Find the matching chain in the major conformation
        ch_fin = [c for c in final_struct.chains() if c.id==ch_min.id]
        # MORE THAN ONE MATCHING RESIDUE -- ERROR?
        if len(ch_fin) > 1:
            raise Exception('MORE THAN ONE MATCHING')
        # If not a matching chain, add to a new chain
        elif len(ch_fin) == 0:
            final_struct.models()[0].append_chain(ch_min.detached_copy())
        # One matching chain, so transfer the residue groups
        elif len(ch_fin) == 1:
            ch_fin = ch_fin[0]
            # Iterate through residue groups
            for rg_min in ch_min.residue_groups():
                # Extract label for this residue
                resid = rg_min.resid()
                chainid = rg_min.parent().id
                # If residue has not changed - do not transfer
                if (chainid, resid) in major_minor_same:
                    continue
                # Find the matching residue group in the major conformation
                rg_fin = [rg for rg in ch_fin.residue_groups() if resid==rg.resid()]
                # MORE THAN ONE MATCHING RESIDUE -- ERROR?
                if len(rg_fin) > 1:
                    raise Exception('MORE THAN ONE MATCHING')
                # If not a matching residue group, add to a new residue group
                elif len(rg_fin)==0:
                    ch_fin.append_residue_group(rg_min.detached_copy())
                # One matching residue, so transfer the atom groups
                elif len(rg_fin)==1:
                    rg_fin = rg_fin[0]
                    # Transfer the atom groups
                    for min_ag in rg_min.atom_groups():
                        rg_fin.append_atom_group(min_ag.detached_copy())

    ######################################################################
    print '===========================================>>>'
    print 'WRITING OUTPUT STRUCTURE'
    ######################################################################

    # Update the atoms numbering
    final_struct.atoms_reset_serial()
    # Write output file
    final_struct.write_pdb_file(params.output)

    return


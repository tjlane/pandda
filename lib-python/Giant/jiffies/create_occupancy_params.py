#!/usr/bin/env pandda.python

import os, sys, copy

import numpy

import iotbx.pdb

import libtbx.phil
from libtbx.utils import Sorry, null_out
from scitbx.array_family import flex

from PANDDAs.jiffies import parse_phil_args
from Giant.Maths.geometry import pairwise_dists
from Giant.Stats.Cluster import find_connected_groups

master_phil = libtbx.phil.parse("""
pdb = None
    .help = 'The major conformation of the protein (normally the unbound or reference structure)'
    .type = str
lig = LIG,UNL
    .help = 'Residues to be occupancy refined (comma separated list of residue identifiers, i.e. lig=LIG or lig=LIG,UNL)'
    .type = str

phenix_occ_out = 'phenix_occ.params'
    .help = 'Output occupancy coupling parameter file for phenix (allowing refinement of the water model in superposition with the ligand)'
    .type = path
refmac_occ_out = 'refmac_occ.params'
    .help = 'Output occupancy coupling parameter file for phenix (allowing refinement of the water model in superposition with the ligand)'
    .type = path

group_dist = 5
    .type = float
    .help = 'Radius around ligands to include in grouped occupancies (must have same conformer as ligand)'

clash_dist = 2
    .type = float
    .help = 'The distance at which we consider the ligand to clash with a nearby (solvent) atom. Clashing atoms will be placed into another refinement group'

overwrite = True
    .type = bool
verbose = True
    .type = bool

""")

def generate_phenix_occupancy_params(occupancy_groups):
    """Using pairs of sets of atom groups, generate occupancy groups for phenix"""

    occ_params_template = """refinement {{\n  refine {{\n    occupancies {{\n{!s}\n    }}\n  }}\n}}"""
    constrained_group_template = """      constrained_group {{\n{!s}\n      }}"""
    selection_template = '        selection = {!s}'
    res_sel_template = '(chain {!s} and resid {!s} and altid {!s})'

    constrained_groups = []
    for (g1, g2) in occupancy_groups:
        if g1:
            # Selections for the first group
            g1_sel = [res_sel_template.format(ag.parent().parent().id, ag.parent().resid(), ag.altloc) for ag in g1]
            g1_str = selection_template.format(' or \\\n                    '.join(g1_sel))
        if g2:
            # Selections for the second group
            g2_sel = [res_sel_template.format(ag.parent().parent().id, ag.parent().resid(), ag.altloc) for ag in g2]
            g2_str = selection_template.format(' or \\\n                    '.join(g2_sel))
        if (g1 and g2): constrained_group_str = constrained_group_template.format('\n'.join([g1_str, g2_str]))
        elif g1:        constrained_group_str = constrained_group_template.format(g1_str)
        elif g2:        constrained_group_str = constrained_group_template.format(g2_str)
        constrained_groups.append(constrained_group_str)
    occ_params = occ_params_template.format('\n'.join(constrained_groups))
    return occ_params

def generate_refmac_occupancy_params(occupancy_groups):
    """Using pairs of sets of atom groups, generate occupancy groups for refmac"""

    selection_template = 'occupancy group id {!s} chain {!s} residue {!s} alt {!s}'
    exclude_template = 'occupancy group alts incomplete {!s} {!s}'
    final_line = 'occupancy refine'

    out_lines = []
    for i_g, (g1, g2) in enumerate(occupancy_groups):
        # Selections for the first group
        if g1: out_lines.extend([selection_template.format(2*i_g+1, ag.parent().parent().id, ag.parent().resseq, ag.altloc) for ag in g1])
        if g2: out_lines.extend([selection_template.format(2*i_g+2, ag.parent().parent().id, ag.parent().resseq, ag.altloc) for ag in g2])
        if (g1 and g2): out_lines.append(exclude_template.format(2*i_g+1, 2*i_g+2))
        elif g1:        out_lines.append(exclude_template.format(2*i_g+1, ''))
        elif g2:        out_lines.append(exclude_template.format(2*i_g+2, ''))
    out_lines.append(final_line)
    occ_params = '\n'.join(out_lines)
    return occ_params

def run(params):

    ######################################################################
    # VALIDATE INPUT
    ######################################################################

    assert params.pdb, 'No PDB File Provided'
    assert params.lig, 'No Ligand Identifier Provided (e.g. LIG or LIG,UNL)'

    assert params.phenix_occ_out or params.refmac_occ_out, 'Must specify at least one of refmac or phenix output files'

    if params.phenix_occ_out and os.path.exists(params.phenix_occ_out):
        if params.overwrite: os.remove(params.phenix_occ_out)
        else: raise Exception('File already exists: {}'.format(params.phenix_occ_out))
    if params.refmac_occ_out and os.path.exists(params.refmac_occ_out):
        if params.overwrite: os.remove(params.refmac_occ_out)
        else: raise Exception('File already exists: {}'.format(params.refmac_occ_out))

    ######################################################################
    # READ IN INPUT FILES
    ######################################################################

    # Read in the ligand file and set each residue to the requested conformer
    pdb_obj = iotbx.pdb.hierarchy.input(params.pdb)

    ######################################################################
    # Iterate through and create refinement groups
    ######################################################################

    # Ligand resname identifiers
    lig_resns = params.lig.split(',')
    if params.verbose: print 'Looking for ligands with resname {!s}'.format(' or '.join(lig_resns))

    # Residues matching the resn that we're looking for
    matching_resn = [ag for ag in pdb_obj.hierarchy.atom_groups() if ag.resname in lig_resns]
    if not matching_resn: raise Sorry('No Matching Ligands found')
    # Allow us to exclude resn to avoid using twice
    use_resn = [True]*len(matching_resn)

    occupancy_groups = []

    for i_lig_ag, lig_ag in enumerate(matching_resn):

        if params.verbose:
            print '============================================>'
            print 'Creating Occupancy Group:', lig_ag.id_str()

        # If excluded - skip (has probably been used already)
        if not use_resn[i_lig_ag]: continue

        ######################################################################
        # SELECT RESIDUES OF THE SAME CONFORMER AROUND THE LIGAND
        ######################################################################

        # Find nearby residues with the same conformer
        altloc_ag = [lig_ag] + [ag for ag in pdb_obj.hierarchy.atom_groups() if (lig_ag.altloc == ag.altloc) and (lig_ag.id_str() != ag.id_str())]

        if len(altloc_ag) > 1:
            # Calculate the distances between these residues
            connection_matrix = numpy.array([[(pairwise_dists(ag1.atoms().extract_xyz(), ag2.atoms().extract_xyz()).min() < params.group_dist) for ag2 in altloc_ag] for ag1 in altloc_ag])
            # Clusters these residues by cutoff distance params.group_dist
            clusters = find_connected_groups(connection_matrix=connection_matrix)

            # Select the cluster containing the ligand (the first residue)
            n_lig_clust = clusters[0]
            # Select the atom groups in this group for the selection
            ligand_group_ag = [ag for i_ag, ag in enumerate(altloc_ag) if clusters[i_ag] == n_lig_clust]
        else:
            ligand_group_ag = altloc_ag

        # Check to see if any of the ligands in matching_resn have been found and remove them
        if params.verbose: print '============================================>'
        for ag in ligand_group_ag:
            if ag in matching_resn:
                if params.verbose: print 'FOUND GROUP FOR:', ag.id_str()
                use_resn[matching_resn.index(ag)] = False

        ######################################################################
        # FIND RESIDUES THAT CLASH WITH THE GROUP
        ######################################################################

        # Find nearby residues with a different conformer (excluding those with more than two conformers as the constraints are too difficult)
        diff_altloc_ag = [ag for ag in pdb_obj.hierarchy.atom_groups() if (ag.altloc) and (lig_ag.altloc != ag.altloc) and (len(ag.parent().conformers())<=2)]
        # Remove those that have multiple conformers other than the ligand conformer (i.e. 'A','B','C' where ligand is conformer 'C')
        clash_ag = [ag for ag in diff_altloc_ag if [1 for gr_ag in ligand_group_ag if (pairwise_dists(gr_ag.atoms().extract_xyz(), ag.atoms().extract_xyz()).min() < params.clash_dist)]]

        ######################################################################
        # APPEND TO OUTPUT
        ######################################################################

        occupancy_groups.append((ligand_group_ag, clash_ag))

    ######################################################################
    # GENERATE OCCUPANCY RESTRAINTS FOR REFINEMENT
    ######################################################################

    ##########################################
    # REFMAC
    ##########################################
    if params.refmac_occ_out:
        if params.verbose:
            print '============================================>'
            print 'CREATING REFMAC OCCUPANCY REFINEMENT PARAMETERS'
        occ_params = generate_refmac_occupancy_params(occupancy_groups=occupancy_groups)
        print '============================================>'
        print 'REFMAC Occupancy Refinement Parameter File Output'
        print '============================================>'
        print occ_params
        with open(params.refmac_occ_out, 'w') as fh: fh.write(occ_params)

    ##########################################
    # PHENIX PARAMETERS
    ##########################################
    if params.phenix_occ_out:
        if params.verbose:
            print '============================================>'
            print 'CREATING PHENIX OCCUPANCY REFINEMENT PARAMETERS'
        occ_params = generate_phenix_occupancy_params(occupancy_groups=occupancy_groups)
        print '============================================>'
        print 'PHENIX Occupancy Refinement Parameter File Output'
        print '============================================>'
        print occ_params
        with open(params.phenix_occ_out, 'w') as fh: fh.write(occ_params)


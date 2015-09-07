#!/usr/bin/env pandda.python

import os, sys, copy

import numpy

import iotbx.pdb

import libtbx.phil
from libtbx.utils import Sorry, null_out
from scitbx.array_family import flex

from PANDDAs.jiffies import parse_phil_args
from Giant.Maths import pairwise_dists
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

verbose = True
    .type = bool

""")

#out_phil = libtbx.phil.parse("""
#refinement
#{
#    refine
#    {
#        occupancies
#        {
#            constrained_group
#                .multiple = True
#            {
#                selection = None
#                    .multiple = True
#                    .type = str
#            }
#        }
#    }
#}
#""")

def run(params):

    ######################################################################
    # VALIDATE INPUT
    ######################################################################

    assert params.pdb, 'No PDB File Provided'
    assert params.lig, 'No Ligand Identifier Provided (e.g. LIG or LIG,UNL)'

    assert params.phenix_occ_out or params.refmac_occ_out, 'Must specify refmac or phenix output file'

    if params.phenix_occ_out:
        assert not os.path.exists(params.phenix_occ_out)
    if params.refmac_occ_out:
        assert not os.path.exists(params.refmac_occ_out)

    ######################################################################
    # READ IN INPUT FILES
    ######################################################################

    # Read in the ligand file and set each residue to the requested conformer
    pdb_obj = iotbx.pdb.hierarchy.input(params.pdb)

#    ######################################################################
#    # CREATE OUTPUT PHIL READY TO BE POPULATED
#    ######################################################################
#
#    # Format the output phil (for Phenix Refine)
#    out_params = out_phil.extract()
#    # Extract the constrained_group scope object
#    constraint_template = out_params.refinement.refine.occupancies.constrained_group[0]
#    # Clear the constrained_group list
#    out_params.refinement.refine.occupancies.constrained_group = []
#    # Clear the selection list
#    constraint_template.selection = []

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

        # If excluded - skip
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

        # Find nearby residues with a different conformer
        diff_altloc_ag = [ag for ag in pdb_obj.hierarchy.atom_groups() if (ag.altloc) and (lig_ag.altloc != ag.altloc)]
        clash_ag = [ag for ag in diff_altloc_ag if [1 for gr_ag in ligand_group_ag if (pairwise_dists(gr_ag.atoms().extract_xyz(), ag.atoms().extract_xyz()).min() < params.clash_dist)]]

        ######################################################################
        # APPEND TO OUTPUT
        ######################################################################

        occupancy_groups.append((ligand_group_ag, clash_ag))

    ######################################################################
    # GENERATE OCCUPANCY RESTRAINTS FOR REFINEMENT
    ######################################################################

    # REFMAC

    if params.refmac_occ_out:

        if params.verbose:
            print '============================================>'
            print 'CREATING REFMAC OCCUPANCY REFINEMENT PARAMETERS'

        selection_template = 'occupancy group id {!s} chain {!s} residue {!s} alt {!s}'
        exclude_template = 'occupancy group alts incomplete {!s} {!s}'
        final_line = 'occupancy refine'

        out_lines = []
        for i_g, (g1, g2) in enumerate(occupancy_groups):

            # Selections for the first group
            out_lines.extend([selection_template.format(2*i_g+1, ag.parent().parent().id, ag.parent().resseq, ag.altloc) for ag in g1])
            out_lines.extend([selection_template.format(2*i_g+2, ag.parent().parent().id, ag.parent().resseq, ag.altloc) for ag in g2])
            out_lines.append(exclude_template.format(2*i_g+1, 2*i_g+2))

        out_lines.append(final_line)

        # Combine into final template
        occ_params = '\n'.join(out_lines)
        print '============================================>'
        print 'REFMAC Occupancy Refinement Parameter File Output'
        print '============================================>'
        print occ_params
        # Write
        with open(params.refmac_occ_out, 'w') as fh:
            fh.write(occ_params)

    # PHENIX

    if params.phenix_occ_out:

        if params.verbose:
            print '============================================>'
            print 'CREATING PHENIX OCCUPANCY REFINEMENT PARAMETERS'

        occ_params_template = """refinement {{\n  refine {{\n    occupancies {{\n{!s}\n    }}\n  }}\n}}"""
        constrained_group_template = """      constrained_group {{\n{!s}\n      }}"""
        selection_template = '        selection = {!s}'
        #res_sel_template = '(chain "{!s}" and resid "{!s}" and altid "{!s}")'
        res_sel_template = '(chain {!s} and resid {!s} and altid {!s})'

        constrained_groups = []
        for (g1, g2) in occupancy_groups:
            # Selections for the first group
            g1_sel = [res_sel_template.format(ag.parent().parent().id, ag.parent().resid(), ag.altloc) for ag in g1]
            g1_str = selection_template.format(' or \\\n                    '.join(g1_sel))
            # Selections for the second group
            g2_sel = [res_sel_template.format(ag.parent().parent().id, ag.parent().resid(), ag.altloc) for ag in g2]
            g2_str = selection_template.format(' or \\\n                    '.join(g2_sel))
            # Create the constrained group strings
            constrained_group_str = constrained_group_template.format('\n'.join([g1_str, g2_str]))
            constrained_groups.append(constrained_group_str)
        # Combine into final template
        occ_params = occ_params_template.format('\n'.join(constrained_groups))
        print '============================================>'
        print 'PHENIX Occupancy Refinement Parameter File Output'
        print '============================================>'
        print occ_params
        # Write
        with open(params.phenix_occ_out, 'w') as fh:
            fh.write(occ_params)


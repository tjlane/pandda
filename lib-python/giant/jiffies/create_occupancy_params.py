import os, sys, copy

import numpy

import iotbx.pdb

import libtbx.phil
from libtbx.utils import Sorry, null_out
from scitbx.array_family import flex

from giant.maths.geometry import pairwise_dists, is_within
from giant.stats.cluster import find_connected_groups, generate_group_idxs

#######################################

blank_arg_prepend = {'.pdb' : 'pdb='}

master_phil = libtbx.phil.parse("""
pdb = None
    .help = 'The major conformation of the protein (normally the unbound or reference structure)'
    .type = str
lig = LIG,UNL,DRG,FRG
    .help = 'Residues to be occupancy refined (comma separated list of residue identifiers, i.e. lig=LIG or lig=LIG,UNL)'
    .type = str

phenix_occ_out = 'occ_phenix.params'
    .help = 'Output occupancy coupling parameter file for phenix (allowing refinement of the water model in superposition with the ligand)'
    .type = path
refmac_occ_out = 'occ_refmac.params'
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

#######################################

def generate_phenix_occupancy_params(occupancy_groups):
    """Using pairs of sets of atom groups, generate occupancy groups for phenix"""

    occ_params_template = """refinement {{\n  refine {{\n    occupancies {{\n{!s}\n    }}\n  }}\n}}"""
    constrained_group_template = """      constrained_group {{\n{!s}\n      }}"""
    selection_template = '        selection = {!s}'
    res_sel_template = '(chain {!s} and resid {!s} and altid "{!s}")'

    constrained_groups = []
    for groups in occupancy_groups:
        all_g_strs = []
        for g in groups:
            # Selections for each group
            g_sel = sorted([res_sel_template.format(ag.parent().parent().id, ag.parent().resid(), ag.altloc) for ag in g])
            g_str = selection_template.format(' or \\\n                    '.join(g_sel))
            all_g_strs.append(g_str)
        constrained_group_str = constrained_group_template.format('\n'.join(all_g_strs))
        constrained_groups.append(constrained_group_str)
    occ_params = occ_params_template.format('\n'.join(constrained_groups))
    return occ_params

def generate_refmac_occupancy_params(occupancy_groups):
    """Using pairs of sets of atom groups, generate occupancy groups for refmac"""

    selection_template = 'occupancy group id {!s} chain {!s} residue {!s} alt {!s}'
    exclude_template = 'occupancy group alts incomplete {!s}'
    final_line = 'occupancy refine'

    out_lines = []
    g_idx = 1
    for group in occupancy_groups:
        # Selections for each group
        for i_g, g in enumerate(group):
            out_lines.extend(sorted([selection_template.format(g_idx+i_g, ag.parent().parent().id, ag.parent().resseq, ag.altloc) for ag in g]))
        out_lines.append(exclude_template.format(' '.join(map(str,range(g_idx, g_idx+len(group))))))
        # Increment the group number
        g_idx += len(group)

    out_lines.append(final_line)
    occ_params = '\n'.join(out_lines)
    return occ_params

#######################################

def run(params):

    ######################################################################
    # VALIDATE INPUT
    ######################################################################

    assert params.pdb, 'No PDB File Provided'
    assert params.lig, 'No Ligand Identifier Provided (e.g. LIG or LIG,UNL)'

#    assert params.phenix_occ_out or params.refmac_occ_out, 'Must specify at least one of refmac or phenix output files'

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
    matching_resn = [ag for ag in pdb_obj.hierarchy.atom_groups() if (ag.resname in lig_resns) and ag.altloc]
    if not matching_resn: raise Sorry('No Matching Ligands found')
    # Allow us to exclude resn to avoid using twice
    use_resn = [True]*len(matching_resn)

    occupancy_groups = []

    for i_lig_ag, lig_ag in enumerate(matching_resn):

        # If excluded - skip (has probably been used already)
        if not use_resn[i_lig_ag]: continue

        if params.verbose:
            print '============================================>'
            print 'Creating Occupancy Group:', lig_ag.id_str(), lig_ag.altloc

        ######################################################################
        # SELECT RESIDUES OF THE SAME CONFORMER AROUND THE LIGAND
        ######################################################################

        # Find residues with the same conformer
        altloc_ag = [lig_ag] + [ag for ag in pdb_obj.hierarchy.atom_groups() if (lig_ag.altloc == ag.altloc) and (lig_ag.id_str() != ag.id_str())]
        # Cluster the residues and chose those near the residue of interest
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

        ######################################################################
        # FIND RESIDUES THAT CLASH WITH THE GROUP
        ######################################################################

        # Find atom groups with different conformers to the residue of interest
        diff_altloc_ag = [ag for ag in pdb_obj.hierarchy.atom_groups() if (ag.altloc) and (lig_ag.altloc != ag.altloc)]
        # Find atom groups that clash with the residue of interest
        clash_ag = [ag for ag in diff_altloc_ag if [1 for gr_ag in ligand_group_ag if is_within(params.clash_dist, gr_ag.atoms().extract_xyz(), ag.atoms().extract_xyz())]]

        # Check to see if any of the ags in matching_resn have been found and remove them
        if params.verbose: print '============================================>'
        for ag in ligand_group_ag+clash_ag:
            if ag in matching_resn:
                if params.verbose: print 'Creating Occupancy Group for:', ag.id_str(), ag.altloc
                use_resn[matching_resn.index(ag)] = False

        ######################################################################
        # GROUP ALL OF THESE RESIDUES INTO OCCUPANCY GROUPS
        ######################################################################

        # Expand atom groups to full residue groups
        all_rgs = set([ag.parent() for ag in ligand_group_ag+clash_ag])

        # Extract all atom groups from these residue groups
        all_ags = []; [all_ags.extend(rg.atom_groups()) for rg in all_rgs]

        # Pull out the matching resn ags again
        matching_ags = [ag for ag in all_ags if (ag.resname in lig_resns)]
        matching_alt = set([ag.altloc for ag in matching_ags])
        print 'Looking for residues with altlocs: {}'.format(','.join(sorted(matching_alt)))

        # Filter out residue groups
        filtered_rgs = [rg for rg in all_rgs if (len([ag.altloc for ag in rg.atom_groups()]) == 1) or matching_alt.intersection([ag.altloc for ag in rg.atom_groups()])]
        filtered_ags = []; [filtered_ags.extend(rg.atom_groups()) for rg in filtered_rgs]
        # Remove blank altlocs
        filtered_ags = [ag for ag in filtered_ags if ag.altloc]
        print '{} residue groups, {} atom groups'.format(len(filtered_rgs), len(filtered_ags))

        filtered_alt = [ag.altloc for ag in filtered_ags]
        groups = []
        for a, g in generate_group_idxs(filtered_alt):
            groups.append([filtered_ags[i] for i in g])

        ######################################################################
        # APPEND TO OUTPUT
        ######################################################################

        occupancy_groups.append(groups)

    ######################################################################
    # GENERATE OCCUPANCY RESTRAINTS FOR REFINEMENT
    ######################################################################

    ##########################################
    # REFMAC
    ##########################################
    if params.verbose:
        print '============================================>'
        print 'CREATING REFMAC OCCUPANCY REFINEMENT PARAMETERS'
    occ_params = generate_refmac_occupancy_params(occupancy_groups=occupancy_groups)
    if params.verbose:
        print '============================================>'
        print 'REFMAC Occupancy Refinement Parameter File Output'
        print '============================================>'
        print occ_params
    if params.refmac_occ_out:
        with open(params.refmac_occ_out, 'w') as fh: fh.write(occ_params)

    ##########################################
    # PHENIX PARAMETERS
    ##########################################
    if params.verbose:
        print '============================================>'
        print 'CREATING PHENIX OCCUPANCY REFINEMENT PARAMETERS'
    occ_params = generate_phenix_occupancy_params(occupancy_groups=occupancy_groups)
    if params.verbose:
        print '============================================>'
        print 'PHENIX Occupancy Refinement Parameter File Output'
        print '============================================>'
        print occ_params
    if params.phenix_occ_out:
        with open(params.phenix_occ_out, 'w') as fh: fh.write(occ_params)

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)

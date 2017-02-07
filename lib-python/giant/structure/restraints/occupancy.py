import numpy

from libtbx.utils import Sorry
from scitbx.array_family import flex

from giant.maths.geometry import is_within
from giant.stats.cluster import find_connected_groups, generate_group_idxs
from giant.structure.formatting import GenericSelection
from giant.structure.select import non_h

def filter_by_distance(atom_groups, xyz, cutoff):
    return [ag for ag in atom_groups if is_within(cutoff, ag.atoms().extract_xyz(), xyz)]
def filter_by_altloc(atom_groups, altloc):
    return [ag for ag in atom_groups if ag.altloc==altloc]

def cluster_atom_groups(atom_groups, cutoff):
    if len(atom_groups) == 1:
        return [0]
    connection_matrix = numpy.array([[is_within(cutoff, ag1.atoms().extract_xyz(), ag2.atoms().extract_xyz()) for ag2 in atom_groups] for ag1 in atom_groups])
    clusters = find_connected_groups(connection_matrix=connection_matrix)
    return clusters

def simple_occupancy_groups(hierarchy, include_single_conformer_groups=False, verbose=False):
    """Given a selection, return the default occupancy groups"""

    occupancy_groups = []
    # Iterate through the default occupancy groups
    for g in hierarchy.occupancy_groups_simple():
        # Skip single groups
        if (len(g) == 1) and (len(g[0]) == 1):
            if verbose:
                print 'Not making simple restraints for single-atom groups:', ','.join([GenericSelection.format(a) for a in hierarchy.select(flex.size_t(g[0])).atoms()])
            continue
        if (len(g) == 1) and (not include_single_conformer_groups):
            if verbose:
                print 'Not making simple restraints for single-conformer groups:\n\t', '\n\t'.join([GenericSelection.format(a) for a in hierarchy.select(flex.size_t(g[0])).atom_groups()])
            continue
        selections = []
        for sel in g:
            ags = [GenericSelection.to_dict(ag) for ag in hierarchy.select(flex.size_t(sel)).atom_groups()]
            selections.append(ags)
        occupancy_groups.append(selections)
    return occupancy_groups


def overlapping_occupancy_groups(hierarchy, resnames, group_dist, overlap_dist, verbose=False):

    # Remove hydrogens to prevent ridiculous amounts of restraints
    hierarchy = non_h(hierarchy)
    # Extract all altlocs and ags with altlocs
    all_altlocs = list(hierarchy.altloc_indices())
    all_alt_ags = [ag for ag in hierarchy.atom_groups() if ag.altloc]

    # Record for each altloc
    # - atom groups for each altloc
    # - assigment of each ag to a cluster of ags
    cluster_dict = {}

    if verbose:
        print '-------------------------------------->'
        print ''
        print 'Generating groups of nearby alternate conformers (cutoff {}A)'.format(group_dist)
        print ''

    for altloc in all_altlocs:
        if altloc=='': continue
        # Select atom groups with this altloc
        altloc_ags = filter_by_altloc(all_alt_ags, altloc)
        # Cluster the atom groups
        altloc_clusters = cluster_atom_groups(altloc_ags, cutoff=group_dist)
        # Dictionary mapping altlocs to ags to clusters
        cluster_dict[altloc] = (altloc_ags, altloc_clusters)
        if verbose:
            print '- altloc {}: {} residues clustered into {} clusters'.format(altloc, len(altloc_ags), len(set(altloc_clusters)))
    if verbose:
        print ''

    # Find atom_groups with the selected resnames
    seed_ags = [ag for ag in all_alt_ags if (ag.resname in resnames)]
    # List of 2-length tuples (containing constrained pairs)
    constrain_groups = []

    # Loop until all atom groups have been used
    while seed_ags:
        # Pick the first residue to focus on
        focus_ag = seed_ags.pop(0)

        # Find which cluster this ag is in
        altloc_ags, altloc_clusters = cluster_dict[focus_ag.altloc]
        focus_clust = altloc_clusters[altloc_ags.index(focus_ag)]
        # Extract all ags in this cluster
        group_ags = [ag for i,ag in enumerate(altloc_ags) if altloc_clusters[i] == focus_clust]
        group_xyz = group_ags[0].atoms().extract_xyz()
        for ag in group_ags[1:]:
            group_xyz = group_xyz.concatenate(ag.atoms().extract_xyz())

        if verbose:
            print '-------------------------------------->'
            print ''
            print 'Creating occupancy group based on: {}'.format(GenericSelection.to_str(focus_ag))
            print '- this residue is part of alternate conformer {}'.format(focus_ag.altloc)
            print '- there are {} atom groups in this cluster'.format(len(group_ags))
            print ''
            print 'Looking for overlapping groups of residues with different alternate conformers:'
            print ''

        tmp_constrain_groups = []
        for altloc in all_altlocs:
            # Skip blank altloc or the selected altloc
            if altloc=='' or altloc==focus_ag.altloc:
                continue
            # Find all ags for this altloc that overlap with the selected cluster
            altloc_ags, altloc_clusters = cluster_dict[altloc]
            overlap_ags = filter_by_distance(atom_groups=altloc_ags, xyz=group_xyz, cutoff=overlap_dist)
            overlap_clusts = sorted(set([altloc_clusters[altloc_ags.index(ag)] for ag in overlap_ags]))

            if verbose:
                print '- altloc {}: overlaps with {} group(s) of residues'.format(altloc, len(overlap_clusts))

            for cluster in overlap_clusts:
                tmp_constrain_groups.append(((focus_ag.altloc, focus_clust),(altloc, cluster)))

        # Add to the complete list
        if tmp_constrain_groups:
            if verbose:
                print ''
                print 'Occupancy groups for this residue'
                print '- {} overlapping group(s) found'.format(len(tmp_constrain_groups))
                print '- creating {} occupancy group constraints'.format(len(tmp_constrain_groups))
            constrain_groups.extend(tmp_constrain_groups)
        else:
            if verbose:
                print '...no overlapping groups found.'
                print '- not creating any occupancy groups for this residue'
        if verbose:
            print ''

        # Remove any used seed groups
        [seed_ags.remove(ag) for ag in group_ags if ag in seed_ags]

    # Filter duplicated restraint groups
    tmp = []
    for g in map(sorted,constrain_groups):
        if g not in tmp: tmp.append(g)
    constrain_groups = tmp

    # Format to generic residue selections
    occupancy_groups = []
    for g in constrain_groups:
        ags_1 = [GenericSelection.to_dict(ag) for ag,c in zip(*cluster_dict[g[0][0]]) if c==g[0][1]]
        ags_2 = [GenericSelection.to_dict(ag) for ag,c in zip(*cluster_dict[g[1][0]]) if c==g[1][1]]
        occupancy_groups.append([ags_1,ags_2])

    return occupancy_groups

#def _find_occupancy_groups(hierarchy, resnames, group_dist, clash_dist, verbose=False):
#
#    # Residues matching the resn that we're looking for
#    if not matching_resn: raise Sorry('No Matching Ligands found')
#    # Allow us to exclude resn to avoid using twice
#    use_resn = [True]*len(matching_resn)
#
#    occupancy_groups = []
#
#    for i_query_ag, query_ag in enumerate(matching_resn):
#
#        # If excluded - skip (has probably been used already)
#        if not use_resn[i_query_ag]: continue
#
#        if verbose:
#            print '============================================>'
#            print 'Creating Occupancy Group:', query_ag.id_str(), query_ag.altloc
#
#        ######################################################################
#        # SELECT RESIDUES OF THE SAME CONFORMER AROUND THE LIGAND
#        ######################################################################
#
#        # Find residues with the same conformer
#        altloc_ag = [query_ag] + [ag for ag in hierarchy.atom_groups() if (query_ag.altloc == ag.altloc) and (query_ag.id_str() != ag.id_str())]
#        # Cluster the residues and chose those near the residue of interest
#        if len(altloc_ag) > 1:
#
#            # Select the cluster containing the ligand (the first residue)
#            n_query_clust = clusters[0]
#            # Select the atom groups in this group for the selection
#            ligand_group_ag = [ag for i_ag, ag in enumerate(altloc_ag) if clusters[i_ag] == n_query_clust]
#        else:
#            ligand_group_ag = altloc_ag
#
#        ######################################################################
#        # FIND RESIDUES THAT CLASH WITH THE GROUP
#        ######################################################################
#
#        # Find atom groups with different conformers to the residue of interest
#        diff_altloc_ag = [ag for ag in hierarchy.atom_groups() if (ag.altloc) and (query_ag.altloc != ag.altloc)]
#        # Find atom groups that clash with the residue of interest
#        clash_ag = [ag for ag in diff_altloc_ag if [1 for gr_ag in ligand_group_ag if is_within(clash_dist, gr_ag.atoms().extract_xyz(), ag.atoms().extract_xyz())]]
#
#        # Check to see if any of the ags in matching_resn have been found and remove them
#        if verbose: print '============================================>'
#        for ag in ligand_group_ag+clash_ag:
#            if ag in matching_resn:
#                if verbose: print 'Creating Occupancy Group for:', ag.id_str(), ag.altloc
#                use_resn[matching_resn.index(ag)] = False
#
#        ######################################################################
#        # GROUP ALL OF THESE RESIDUES INTO OCCUPANCY GROUPS
#        ######################################################################
#
#        # Expand atom groups to full residue groups
#        all_rgs = set([ag.parent() for ag in ligand_group_ag+clash_ag])
#
#        # Extract all atom groups from these residue groups
#        all_ags = []; [all_ags.extend(rg.atom_groups()) for rg in all_rgs]
#
#        # Pull out the matching resn ags again
#        matching_ags = [ag for ag in all_ags if (ag.resname in resnames)]
#        matching_alt = set([ag.altloc for ag in matching_ags])
#        print 'Looking for residues with altlocs: {}'.format(','.join(sorted(matching_alt)))
#
#        # Filter out residue groups
#        filtered_rgs = [rg for rg in all_rgs if (len([ag.altloc for ag in rg.atom_groups()]) == 1) or matching_alt.intersection([ag.altloc for ag in rg.atom_groups()])]
#        filtered_ags = []; [filtered_ags.extend(rg.atom_groups()) for rg in filtered_rgs]
#        # Remove blank altlocs
#        filtered_ags = [ag for ag in filtered_ags if ag.altloc]
#        print '{} residue groups, {} atom groups'.format(len(filtered_rgs), len(filtered_ags))
#
#        filtered_alt = [ag.altloc for ag in filtered_ags]
#        groups = []
#        for a, g in generate_group_idxs(filtered_alt):
#            groups.append([filtered_ags[i] for i in g])
#
#        ######################################################################
#        # APPEND TO OUTPUT
#        ######################################################################
#
#        occupancy_groups.append(groups)
#
#    return occupancy_groups

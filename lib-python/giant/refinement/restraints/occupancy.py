import giant.logs as lg
logger = lg.getLogger(__name__)

import copy

import numpy as np
from scitbx.array_family import flex

from giant.common.clustering import (
    find_connected_groups
    )

from giant.structure.common import (
    GetInterestingResnames,
    )

from giant.structure.formatting import (
    Labeller,
    GenericSelection,
    )

from .base import (
    OccupancyGroup,
    OccupancyRestraint,
    OccupancyRestraintList,
    RestraintsCollection,
    )

#####

class _BaseRestraintMaker(object):

    def get_selected_atoms(self, hierarchy):

        from giant.structure.select import non_h
        hierarchy = non_h(hierarchy)

        if self.atom_selection is not None:
            asc = hierarchy.atom_selection_cache()
            selection = asc.selection(self.atom_selection)
            hierarchy = hierarchy.select(selection)

        return hierarchy


class MakeSimpleOccupancyRestraints(_BaseRestraintMaker):

    name = "MakeSimpleOccupancyRestraints"

    def __init__(self,
        single_atom_groups = True,
        single_conformer_groups = True,
        atom_selection = None,
        ):

        self.single_atom_groups = bool(single_atom_groups)
        self.single_conformer_groups = bool(single_conformer_groups)

        self.atom_selection = (
            str(atom_selection)
            if atom_selection is not None
            else None
            )

    def __str__(self):

        s_ = (
            'Task: {name}\n'
            '| single_atom_groups: {single_atom_groups}\n'
            '| single_conformer_groups: {single_conformer_groups}\n'
            '| atom_selection: {atom_selection}\n'
            '`---->'
            ).format(
                name = self.name,
                single_atom_groups = str(self.single_atom_groups),
                single_conformer_groups = str(self.single_conformer_groups),
                atom_selection = str(self.atom_selection),
                )

        return s_.strip()

    def __call__(self, hierarchy):

        occupancy_restraints = []

        hierarchy = self.get_selected_atoms(
            hierarchy = hierarchy,
            )

        ag_set = set()

        # Iterate through the default occupancy groups
        for g in hierarchy.occupancy_groups_simple():

            if (not self.single_atom_groups) and (len(g) == 1) and (len(g[0]) == 1):
                logger.debug(
                    '\n\t'.join([
                        'Not making simple restraints for single-atom groups:',
                        '\n\t'.join([
                            Labeller.format(ag)
                            for ag in hierarchy.select(
                                flex.size_t(g[0])
                                ).atom_groups()
                            ]),
                        ])
                    )
                continue

            if (not self.single_conformer_groups) and (len(g) == 1):
                logger.debug(
                    '\n\t'.join([
                        'Not making simple restraints for single-conformer groups:',
                        '\n\t'.join([
                            Labeller.format(ag)
                            for ag in hierarchy.select(
                                flex.size_t(g[0])
                                ).atom_groups()
                            ]),
                        ])
                    )
                continue

            # All atoms selected by this group
            combined_selection = (
                flex.size_t(np.sort(np.concatenate(g)))
                )

            ag_labels = [
                Labeller.format(ag)
                for ag in hierarchy.select(
                    combined_selection
                    ).atom_groups()
                ]

            # There are duplicates... prevent duplicates
            if ag_set.intersection(ag_labels):
                continue

            # Add to set
            ag_set.update(ag_labels)

            # List of lists
            occupancy_restraints.append(
                OccupancyRestraint(
                    occupancy_groups = [
                        OccupancyGroup(
                            objects = [
                                GenericSelection.to_dict(ag)
                                for ag in hierarchy.select(
                                    flex.size_t(sel)
                                    ).atom_groups()
                                ],
                            )
                        for sel in g
                        ],
                    complete = (len(g) > 1),
                    )
                )

        return RestraintsCollection(
            occupancy_restraints = occupancy_restraints,
            )


class MakeMultiStateOccupancyRestraints(_BaseRestraintMaker):

    name = "MakeMultiStateOccupancyRestraints"

    def __init__(self,
        group_distance_cutoff = 6.0,
        overlap_distance_cutoff = 6.0,
        ignore_common_molecules = True,
        include_resnames_list = None,
        ignore_resnames_list = None,
        set_group_completeness_to = None,
        atom_selection = None,
        make_intra_conformer_distance_restraints = None,
        ):

        # if we have exclude altlocs, need to set group_completeness to False?
        # if exclude_altlocs is None: exclude_altlocs = []
        # if exclude_altlocs == ['']: exclude_altlocs = []

        self.group_distance_cutoff = float(
            group_distance_cutoff
            )

        self.overlap_distance_cutoff = float(
            overlap_distance_cutoff
            )

        self.get_interesting_resnames = GetInterestingResnames(
            ignore_common_molecules = ignore_common_molecules,
            include_resnames_list = include_resnames_list,
            ignore_resnames_list = ignore_resnames_list,
            )

        self.set_group_completeness_to = (
            bool(set_group_completeness_to)
            if set_group_completeness_to is not None
            else None
            )

        self.atom_selection = (
            str(atom_selection)
            if atom_selection is not None
            else None
            )

        self.make_intra_conformer_distance_restraints = (
            make_intra_conformer_distance_restraints
            if make_intra_conformer_distance_restraints is not None
            else None
            )

    def __str__(self):

        s_ = (
            'Task: {name}\n'
            '| group_distance_cutoff: {group_distance_cutoff}\n'
            '| overlap_distance_cutoff: {overlap_distance_cutoff}\n'
            '| get_interesting_resnames: \n'
            '| \t{get_interesting_resnames}\n'
            '| set_group_completeness_to: {set_group_completeness_to}\n'
            '| atom_selection: {atom_selection}\n'
            '`---->'
            ).format(
                name = self.name,
                group_distance_cutoff = str(self.group_distance_cutoff),
                overlap_distance_cutoff = str(self.overlap_distance_cutoff),
                get_interesting_resnames = str(self.get_interesting_resnames).replace('\n','\n| \t'),
                set_group_completeness_to = str(self.set_group_completeness_to),
                atom_selection = str(self.atom_selection),
                )

        return s_.strip()

    def __call__(self, hierarchy):

        hierarchy = self.get_selected_atoms(
            hierarchy = hierarchy,
            )

        rc = RestraintsCollection()

        # Get interesting residues
        interesting_resnames = self.get_interesting_resnames(hierarchy)

        if (not interesting_resnames):
            return rc # return null object

        # Group alternate conformers into clusters (by conformer)
        altloc_clusters_dict = self.cluster_altloc_atoms(
            hierarchy = hierarchy,
            )

        # Group the clusters by finding all overlapping clusters
        altloc_groups = self.find_overlapping_altloc_groups(
            hierarchy = hierarchy,
            altloc_groups_dict = altloc_clusters_dict,
            )

        # Remove groups that don't have a residue of interest (unless None)
        if (interesting_resnames is not None):
            altloc_groups = self.filter_groups(
                hierarchy = hierarchy,
                altloc_groups = altloc_groups,
                filter_resnames = interesting_resnames,
                )

        # Convert groups into restraints
        altloc_occ_restraints = self.make_occupancy_restraints(
            hierarchy = hierarchy,
            altloc_groups = altloc_groups,
            )

        rc.add(altloc_occ_restraints)

        # Make distance restraints within the groups (optional)
        altloc_dist_restraints = self.get_intra_group_distance_restraints(
            hierarchy = hierarchy,
            altloc_groups = altloc_groups,
            )

        rc.add(altloc_dist_restraints)

        return rc

    def cluster_altloc_atoms(self, hierarchy):

        logger.debug('** Clustering alternate conformer atom groups **')

        h_atoms = hierarchy.atoms()

        # Get the indices of each conformer in the structure
        alt_indices = hierarchy.get_conformer_indices() # (0,0,0,1,1,2,2,0,0,0,0,...)

        # Get all the altlocs in the structure
        all_altlocs = list(hierarchy.altloc_indices()) # ['','A','B'...]

        # clusters for each altloc {altloc: [indices1, indices2,...], ...}
        altloc_cluster_dict = {}
        #
        for i_altloc, altloc in enumerate(all_altlocs):

            if altloc.strip() == '':
                continue

            # Get the altloc atoms
            alt_selection = (alt_indices == i_altloc) # flex bool

            # Get the clusters relative to the altloc atoms
            cluster_atom_selections_rel = self.cluster_atoms(
                atoms = h_atoms.select(alt_selection),
                distance_cutoff = self.group_distance_cutoff,
                )

            # Remake the clusters in terms of the full hierarchy
            alt_selection_idx = np.where(np.array(alt_selection))[0]
            cluster_atom_selections = [
                alt_selection_idx[cluster_selection_rel]
                for cluster_selection_rel in cluster_atom_selections_rel
            ]

            # store clusters by altloc
            altloc_cluster_dict[altloc] = (
                cluster_atom_selections
                )

        for k, v in sorted(altloc_cluster_dict.items()):
            logger.debug(
                '\n** Clustered atoms, altloc {k} **\n'.format(k=k)
                )
            for i, vv in enumerate(v):
                logger.debug(
                    'Group {i}\n{s}'.format(
                        i = i+1,
                        s = hierarchy.select(flex.size_t(vv)).as_str(),
                        )
                    )

        return altloc_cluster_dict

    def find_overlapping_altloc_groups(self, hierarchy, altloc_groups_dict):

        logger.debug('** identifying overlapping altloc groups **')

        h_atoms = hierarchy.atoms()

        # flatten dict to 1d list
        all_altloc_group_atoms = []
        #
        for altloc, group_selections in sorted(altloc_groups_dict.items()):

            for group_idx, group_sel in enumerate(group_selections):

                all_altloc_group_atoms.append(
                    (
                        altloc,
                        group_idx,
                        h_atoms.select(
                            flex.size_t(group_sel)
                            )
                        )
                    )

        # How many groups do we have
        n_groups = len(all_altloc_group_atoms)

        # cluster the groups
        distance_matrix = np.empty(
            (n_groups, n_groups),
            dtype = float,
            )
        distance_matrix.fill(-1.)
        np.fill_diagonal(distance_matrix, 0.0)

        # Do in for loops as shouldn't be /too/ many groups
        for i_group_1, (alt_1, i_idx_1, atoms_1) in enumerate(all_altloc_group_atoms):

            for i_group_2, (alt_2, i_idx_2, atoms_2) in enumerate(all_altloc_group_atoms):

                # Triangular
                if (i_group_2 <= i_group_1):
                    continue

                # Don't cluster groups with the same altloc (should be done in "group_atoms" step)
                if (alt_1 == alt_2):
                    continue

                min_distance = self.shortest_distance_between_atom_sets(
                    atoms_1 = atoms_1,
                    atoms_2 = atoms_2,
                    )

                distance_matrix[i_group_1, i_group_2] = min_distance
                distance_matrix[i_group_2, i_group_1] = min_distance

        connection_matrix = (
            distance_matrix >= 0.0
            ) & (
            distance_matrix < self.overlap_distance_cutoff
            )

        # for r in connection_matrix:
        #     logger.debug(np.where(r)[0])

        cluster_indices = np.array(
            find_connected_groups(
                connection_matrix = connection_matrix,
                ),
            dtype = int,
            )

        assert cluster_indices.shape == (n_groups,)

        cluster_selections = [
            np.where(
                cluster_indices == idx,
                )[0]
            for idx in sorted(
                set(cluster_indices)
                )
            ]

        # Format groupings for output
        group_dicts = []
        #
        for cluster_sel in cluster_selections:

            # Get selections and sort by altloc again
            group_dict_tmp = {}
            for i in cluster_sel:
                altloc, idx, _ = all_altloc_group_atoms[i]
                group_dict_tmp.setdefault(altloc, []).append(
                    altloc_groups_dict[altloc][idx]
                    )

            group_dicts.append(
                {
                    k: np.sort(np.concatenate(v))
                    for k, v in group_dict_tmp.items()
                    }
                )

        return group_dicts

    def filter_groups(self, hierarchy, altloc_groups, filter_resnames):

        filter_resnames = set(filter_resnames)

        filtered_groups = []

        for a_group_dict in altloc_groups:

            # Get all the atoms in this group across altlocs
            full_selection = flex.size_t(np.sort(np.concatenate(list(a_group_dict.values()))))

            h = hierarchy.select(full_selection)

            this_resnames = list(h.overall_counts().resnames.keys())

            overlap = filter_resnames.intersection(this_resnames)

            if len(overlap) > 0:
                filtered_groups.append(a_group_dict)

        return filtered_groups

    def make_occupancy_restraints(self, hierarchy, altloc_groups):

        occupancy_restraints = []

        for a_group_dict in altloc_groups:

            occupancy_restraints.append(
                OccupancyRestraint(
                    occupancy_groups = [
                        OccupancyGroup(
                            label = 'Conformer {}'.format(altloc),
                            objects = [
                                GenericSelection.to_dict(o)
                                for o in hierarchy.select(
                                    flex.size_t(atom_selection),
                                    ).atom_groups()
                                ],
                            )
                        for altloc, atom_selection in sorted(
                            a_group_dict.items()
                            )
                        ],
                    )
                )

        occ_r_list = OccupancyRestraintList(
            occupancy_restraints = occupancy_restraints,
            )

        if self.set_group_completeness_to is not None:
            occ_r_list.set_complete(
                self.set_group_completeness_to
                )

        return RestraintsCollection(
            occupancy_restraints = occ_r_list,
            )

    def get_intra_group_distance_restraints(self, hierarchy, altloc_groups):

        logger.debug('** making intra-group distance restraints **')

        rc = RestraintsCollection()

        if self.make_intra_conformer_distance_restraints is None:
            return rc

        alt_indices = hierarchy.get_conformer_indices() # (0,0,0,0,1,1,2,2,0,0,0...)
        all_altlocs = list(hierarchy.altloc_indices()) # ['','A','B']
        alt_blank_sel = (
            alt_indices == all_altlocs.index('')
            ) # flex.bool

        for a_group_dict in altloc_groups:

            atom_selection = copy.deepcopy(alt_blank_sel)

            # Get selection for this group and all main conf atoms
            for altloc, alt_atom_selection in sorted(a_group_dict.items()):

                atom_selection.set_selected(
                    flex.size_t(alt_atom_selection), True,
                    ) # in_place

            rc_group = self.make_intra_conformer_distance_restraints(
                hierarchy = hierarchy.select(atom_selection),
                )

            rc.add(rc_group)

        return rc

    def cluster_atoms(self, atoms, distance_cutoff):

        xyz = np.array(atoms.extract_xyz())

        n_atoms = len(xyz)

        # manual return for one atom
        if (n_atoms == 1):
            return [np.array([0])]

        xyz_grid = xyz.reshape(
            (n_atoms, 1, 3)
            ).repeat(n_atoms, axis=1)

        xyz_grid_diff_sq = np.power(xyz_grid - xyz_grid.transpose((1,0,2)), 2).sum(axis=2)

        connection_matrix = xyz_grid_diff_sq < (distance_cutoff ** 2)

        cluster_indices = np.array(
            find_connected_groups(
                connection_matrix = connection_matrix,
                ),
            dtype = int,
            )

        assert cluster_indices.shape == (n_atoms,)

        cluster_selections = [
            np.array(
                (cluster_indices == idx),
                dtype = bool,
                )
            for idx in sorted(
                set(cluster_indices)
                )
            ]

        # return as bool selections -- could've equally returned as lists of indices

        return cluster_selections

    def shortest_distance_between_atom_sets(self, atoms_1, atoms_2):

        n_atoms_1 = len(atoms_1)
        n_atoms_2 = len(atoms_2)

        xyz_1 = np.array(atoms_1.extract_xyz())
        xyz_2 = np.array(atoms_2.extract_xyz())

        xyz_1_grid = xyz_1.reshape((n_atoms_1, 1, 3)).repeat(n_atoms_2, axis=1)
        xyz_2_grid = xyz_2.reshape((n_atoms_2, 1, 3)).repeat(n_atoms_1, axis=1)

        xyz_diffs_sq = np.power(xyz_1_grid - xyz_2_grid.transpose((1,0,2)), 2).sum(axis=2)

        xyz_diffs_min = xyz_diffs_sq.min() ** 0.5

        # logger.debug(
        #     '\n'.join(
        #         ['Atoms1']+[a.format_atom_record() for a in atoms_1] +
        #         ['Atoms2']+[a.format_atom_record() for a in atoms_2] +
        #         ['Min distance: {}'.format(xyz_diffs_min)]
        #         )
        #     )

        # from IPython import embed; embed(); exist()
        return xyz_diffs_min


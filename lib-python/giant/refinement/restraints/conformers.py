import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy as np

from .base import (
    RestraintsCollection,
    DistanceRestraint,
    )


class _BaseRestraintMaker(object):

    def get_selected_atoms(self, hierarchy):

        if self.exclude_hydrogens is True:
            from giant.structure.select import non_h
            hierarchy = non_h(hierarchy)

        if self.atom_selection is not None:
            asc = hierarchy.atom_selection_cache()
            selection = asc.selection(self.atom_selection)
            hierarchy = hierarchy.select(selection)

        return hierarchy


class MakeIntraConformerRestraints(_BaseRestraintMaker):

    name = "MakeIntraConformerRestraints"

    def __init__(self,
        min_distance_cutoff = 1.6,
        max_distance_cutoff = 4.2,
        distance_restraint_sigma = 0.1,
        #torsion_restraint_sigma = None, TODO
        select_altlocs = None,
        atom_selection = None,
        exclude_hydrogens = True,
        ):

        self.min_distance_cutoff = float(min_distance_cutoff)
        self.max_distance_cutoff = float(max_distance_cutoff)

        self.distance_restraint_sigma = float(distance_restraint_sigma)

        self.select_altlocs = (
            list(select_altlocs)
            if select_altlocs is not None
            else None
            )

        self.atom_selection = (
            str(atom_selection)
            if atom_selection is not None
            else None
            )

        self.exclude_hydrogens = bool(exclude_hydrogens)

        assert self.min_distance_cutoff > 0.0
        assert self.max_distance_cutoff > self.min_distance_cutoff

        self.min_distance_cutoff_sq = (min_distance_cutoff ** 2)
        self.max_distance_cutoff_sq = (max_distance_cutoff ** 2)

    def __str__(self):

        s_ = (
            'Task: {name}\n'
            '| min_distance_cutoff: {min_distance_cutoff}\n'
            '| max_distance_cutoff: {max_distance_cutoff}\n'
            '| distance_restraint_sigma: {distance_restraint_sigma}\n'
            '| select_altlocs: {select_altlocs}\n'
            '| atom_selection: {atom_selection}\n'
            '| exclude_hydrogens: {exclude_hydrogens}\n'
            '`---->'
            ).format(
                name = self.name,
                min_distance_cutoff = self.min_distance_cutoff,
                max_distance_cutoff = self.max_distance_cutoff,
                distance_restraint_sigma = self.distance_restraint_sigma,
                select_altlocs = str(self.select_altlocs),
                atom_selection = str(self.atom_selection),
                exclude_hydrogens = str(self.exclude_hydrogens),
                )

        return s_.strip()

    def __call__(self, hierarchy):

        hierarchy = self.get_selected_atoms(
            hierarchy = hierarchy,
            )

        rc = RestraintsCollection()

        for (altloc, conformer_atoms, other_atoms) in self.iterate_conformer_sets(hierarchy):

            atom_pairs = self.iterate_atom_pairs_within_cutoff(
                conformer_atoms = conformer_atoms,
                other_atoms = other_atoms,
                )

            restraints = self.make_atom_pair_restraints(
                atom_pairs = atom_pairs,
                )

            rc.add(restraints)

        return rc

    def iterate_conformer_sets(self, hierarchy):

        h_atoms = hierarchy.atoms()

        # Get the indices of each conformer in the structure
        alt_indices = hierarchy.get_conformer_indices()

        # Get all the altlocs in the structure
        all_altlocs = list(hierarchy.altloc_indices())

        # Subselect altlocs?
        sel_altlocs = (
            self.select_altlocs
            if self.select_altlocs is not None
            else all_altlocs
            )

        # Get blank altloc atom selection
        alt_blank_sel = (
            alt_indices == all_altlocs.index('')
            ) # flex.bool

        for i_alt, alt in enumerate(all_altlocs):

            if alt.strip() == '':
                continue

            # Skip skipped altlocs
            if (alt not in sel_altlocs):
                continue

            # Make selection for the altloc atoms
            alt_sel = (alt_indices == i_alt) # flex.bool

            # Make selection for altloc atoms and main conf atoms
            combined_sel = (alt_indices == i_alt).set_selected(alt_blank_sel, True)

            yield (
                alt,
                h_atoms.select(alt_sel),
                h_atoms.select(combined_sel),
                )

    def iterate_atom_pairs_within_cutoff(self, conformer_atoms, other_atoms):

        conformer_xyz = np.array(conformer_atoms.extract_xyz())
        other_xyz = np.array(other_atoms.extract_xyz())

        pairwise_distances_sq = np.zeros(
            (len(conformer_xyz), len(other_xyz)),
            dtype = float,
            )

        for i_conformer, conf_xyz in enumerate(conformer_xyz):

            # calculate inter-atom distances
            pairwise_distances_sq[i_conformer] = (
                np.power(other_xyz - conf_xyz, 2).sum(axis=1)
                )

        # Truncate large distances
        pairwise_distances_sq[pairwise_distances_sq > self.max_distance_cutoff_sq] = -1

        # Select distances above cutoff
        atom_pair_indices = np.array(
            np.where(
                pairwise_distances_sq >= self.min_distance_cutoff_sq
                )
            ).T

        atom_check_hash = {}

        for i_at_1, i_at_2 in atom_pair_indices:

            a1 = conformer_atoms[i_at_1]
            a2 = other_atoms[i_at_2]

            hash_lab = tuple(sorted([a1.i_seq,a2.i_seq]))

            # avoid duplicates
            if atom_check_hash.get(hash_lab):
                continue

            atom_check_hash[hash_lab] = 1

            yield (
                a1.fetch_labels(),
                a2.fetch_labels(),
                )

    def make_atom_pair_restraints(self, atom_pairs):

        distance_restraints = []

        for a1, a2 in atom_pairs:

            # Skip atoms in the same residue - TODO

            distance_restraints.append(
                DistanceRestraint(
                    atom1 = a1,
                    atom2 = a2,
                    length = a1.distance(a2),
                    sigma = self.distance_restraint_sigma,
                    )
                )

        return RestraintsCollection(
            distance_restraints = distance_restraints,
            )


class MakeDuplicateConformerRestraints(_BaseRestraintMaker):

    name = "MakeDuplicateConformerRestraints"

    def __init__(self,
        rmsd_cutoff = 0.1,
        distance_restraint_sigma = 0.02,
        atom_selection = None,
        exclude_hydrogens = True,
        ):

        self.rmsd_cutoff = float(rmsd_cutoff)
        self.distance_restraint_sigma = float(distance_restraint_sigma)

        self.atom_selection = (
            str(atom_selection)
            if atom_selection is not None
            else None
            )

        self.exclude_hydrogens = bool(exclude_hydrogens)


    def __str__(self):

        s_ = (
            'Task: {name}\n'
            '| rmsd_cutoff: {rmsd_cutoff}\n'
            '| distance_restraint_sigma: {distance_restraint_sigma}\n'
            '| atom_selection: {atom_selection}\n'
            '| exclude_hydrogens: {exclude_hydrogens}\n'
            '`---->'
            ).format(
                name = self.name,
                rmsd_cutoff = self.rmsd_cutoff,
                distance_restraint_sigma = self.distance_restraint_sigma,
                atom_selection = str(self.atom_selection),
                exclude_hydrogens = str(self.exclude_hydrogens),
                )

        return s_.strip()

    def __call__(self,
        hierarchy,
        ):

        hierarchy = self.get_selected_atoms(
            hierarchy = hierarchy,
            )

        rc = RestraintsCollection()

        altconf_dict = self.get_altconf_dict(
            hierarchy = hierarchy,
            )

        duplicate_residues = self.get_duplicate_residues(
            altconf_dict = altconf_dict,
            )

        restraints = self.make_atom_restraints(
            residue_pairs = duplicate_residues,
            )

        rc.add(restraints)

        return rc

    def get_altconf_dict(self, hierarchy):

        blank_idx = list(hierarchy.altloc_indices()).index('')

        not_blank_sel = (hierarchy.get_conformer_indices() != blank_idx)

        sel_hierarchy = hierarchy.select(not_blank_sel)

        rg_dict = {}

        for rg in sel_hierarchy.residue_groups():

            # Need a unique label
            label = rg.id_str()

            assert label not in rg_dict

            rg_dict[label] = {
                cf.altloc : (
                    cf.only_residue(),
                    dict(zip(
                        cf.only_residue().atoms().extract_name(),
                        cf.only_residue().atoms().extract_xyz(),
                        ))
                    )
                for cf in rg.conformers()
            }

        return rg_dict

    def get_duplicate_residues(self, altconf_dict):

        duplicate_list = []

        for rg_lab, rg_dict in sorted(altconf_dict.items()):

            rg_items = sorted(rg_dict.items())

            for i1, (altloc1, (residue1, xyz_dict1)) in enumerate(rg_items):

                for i2, (altloc2, (residue2, xyz_dict2)) in enumerate(rg_items):

                    if i2 <= i1:
                        continue

                    if residue1.resname != residue2.resname:
                        continue

                    rmsd = self.calculate_rmsd(
                        xyz_dict_1 = xyz_dict1,
                        xyz_dict_2 = xyz_dict2,
                        )

                    if rmsd is None:
                        continue

                    elif rmsd < self.rmsd_cutoff:
                        duplicate_list.append(
                            (residue1, residue2)
                            )

        return duplicate_list

    def make_atom_restraints(self, residue_pairs):

        distance_restraints = []

        for r1, r2 in residue_pairs:

            ats1 = r1.atoms().build_dict()
            ats2 = r2.atoms().build_dict()

            for k in sorted(
                set(ats1.keys()).intersection(ats2.keys())
                ):

                distance_restraints.append(
                    DistanceRestraint(
                        atom1 = ats1[k],
                        atom2 = ats2[k],
                        length = 0.,
                        sigma = self.distance_restraint_sigma,
                        )
                    )

        return RestraintsCollection(
            distance_restraints = distance_restraints,
            )

    def calculate_rmsd(self, xyz_dict_1, xyz_dict_2):

        common_keys = set(xyz_dict_1.keys()).intersection(xyz_dict_2.keys())

        if not common_keys:
            return None

        xyz_1 = np.array([xyz_dict_1[k] for k in common_keys])
        xyz_2 = np.array([xyz_dict_2[k] for k in common_keys])

        d = np.power(xyz_1 - xyz_2, 2).sum() ** 0.5

        return d



import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy as np

from .base import (
    RestraintsCollection,
    DistanceRestraint,
    )

from giant.structure.formatting import (
    Labeller,
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

    _backbone_atom_set = set(['CA','C','N','O'])
    _single_bond_elems = ['C','N','O','S','H','BR','CL','F'] # upper case
    _multi_bond_resnames = ['HOH']

    def __init__(self,
        min_distance_cutoff = 0.1,
        max_distance_cutoff = 4.0,
        distance_restraint_sigma = 0.1,
        #torsion_restraint_sigma = None, TODO
        select_altlocs = None,
        atom_selection = None,
        exclude_hydrogens = True,
        filter_c_x_pairs = True,
        filter_c_c_pairs = True,
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

        self.filter_c_x_pairs = bool(filter_c_x_pairs)
        self.filter_c_c_pairs = bool(filter_c_c_pairs)

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

        for (altloc, conformer_h, other_h) in self.iterate_conformer_sets(hierarchy):

            logger.debug(
                (
                    '** Generating restraints for altloc {altloc} **\n'
                    '\tAlt. conf. atoms: {n_conf}\n'
                    '\tMain+Alt conf atoms: {n_main}\n'
                    ).format(
                    altloc = altloc,
                    n_conf = conformer_h.atoms().size(),
                    n_main = other_h.atoms().size(),
                    )
                )

            atom_pairs = self.iterate_atom_pairs_within_cutoff(
                conformer_hierarchy = conformer_h,
                other_hierarchy = other_h,
                )

            atom_pairs = self.prune_atom_pairs(
                atom_pairs = atom_pairs,
                )

            restraints = self.make_atom_pair_restraints(
                atom_pairs = atom_pairs,
                label = 'Restraints for altloc {alt}'.format(alt=altloc),
                )

            logger.debug(
                str(restraints)
                )

            rc.add(restraints)

        return rc

    def iterate_conformer_sets(self, hierarchy):

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
            # (doing it this way keeps atom ordering -- desirable but inefficient?)
            combined_sel = (alt_indices == i_alt).set_selected(alt_blank_sel, True) # flex.bool

            yield (
                alt,
                hierarchy.select(alt_sel),
                hierarchy.select(combined_sel),
                )

    def iterate_atom_pairs_within_cutoff(self, conformer_hierarchy, other_hierarchy):

        conformer_xyz = np.array(
            conformer_hierarchy.atoms().extract_xyz()
            )
        other_xyz = np.array(
            other_hierarchy.atoms().extract_xyz()
            )

        pairwise_distances_sq = np.zeros(
            (len(conformer_xyz), len(other_xyz)),
            dtype = float,
            )

        for i_conformer, conf_xyz in enumerate(conformer_xyz):

            # calculate inter-atom distances -- only for i > j
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

        # Sort the atom_pair_indices by the distance
        # so that the shortest distances are prioritised
        sorted_atom_pair_indices = sorted(
            atom_pair_indices.tolist(),
            key = lambda i_j: pairwise_distances_sq[i_j[0],i_j[1]],
            )

        atom_check_hash = {}

        # This is annoying -- could probably be optimised further
        conformer_atoms_with_labels = list(
            conformer_hierarchy.atoms_with_labels()
            )
        other_atoms_with_labels = list(
            other_hierarchy.atoms_with_labels()
            )

        for i_at_conf, i_at_other in sorted_atom_pair_indices:

            a1 = conformer_atoms_with_labels[i_at_conf]
            a2 = other_atoms_with_labels[i_at_other]

            # sort atoms by serial (unique)
            a1, a2 = sorted(
                [a1, a2],
                key = lambda a: a.serial_as_int(),
                )

            hash_lab = (
                a1.serial_as_int(),
                a2.serial_as_int(),
                )

            logger.debug(
                    'Atom pair: \n\t{l1}\n\t{l2}'.format(
                        l1 = Labeller.format(a1),
                        l2 = Labeller.format(a2),
                        )
                    )

            # avoid exact duplicates
            if atom_check_hash.get(hash_lab):

                logger.debug('Already done, skipping.\n')

                continue

            # Record this pair as "seen before"
            atom_check_hash[hash_lab] = 1

            # Check not in the same residue, etc
            if not self.is_valid_atom_pair(a1,a2):

                logger.debug('Not valid atom pair, skipping.\n')

                continue

            logger.debug('Keeping atom pair.\n')

            yield (a1, a2)

    def is_valid_atom_pair(self, atom1, atom2):

        if atom1.chain_id != atom2.chain_id:
            logger.debug('USE: Different chains')
            return True

        # do not allow same residue
        if atom1.resid() == atom2.resid():
            logger.debug('SKIP: Same atom group')
            return False

        # now quick checks

        if atom1.hetero or atom2.hetero:
            logger.debug('USE: At least one atom is HETATM')
            return True

        if self.filter_c_c_pairs is True:
            if (atom1.element.strip().upper() == 'C') and (atom2.element.strip().upper() == 'C'):
                logger.debug('SKIP: Skipping C-C pairs')
                return False

        # only allow interactions between adjacent sidechains

        if abs(
            atom1.resseq_as_int() - atom2.resseq_as_int()
            ) < 2:

            if self._backbone_atom_set.intersection(
                [
                    atom1.name.strip(),
                    atom2.name.strip(),
                    ]
                ):
                logger.debug('SKIP: Main chains of adjacent residues')
                return False

        return True

    def prune_atom_pairs(self, atom_pairs):

        # sort the atom pairs by residue-interactions
        sorting_hash = {}

        for a1, a2 in atom_pairs:

            assert a1.serial_as_int() < a2.serial_as_int(), "must be preordered"

            # Sort by residue interactions
            hash_key = (
                a1.chain_id+a1.resid(),
                a2.chain_id+a2.resid(),
                )

            hash_tuple = (
                a1,
                a2,
                a1.distance(a2),
                )

            sorting_hash.setdefault(
                hash_key, []
                ).append(
                hash_tuple
                )

        # Now return one one interaction for each atom for each residue pairing

        for k, atom_pairs in sorted(sorting_hash.items()):

            logger.debug(
                'Pruning atom pairs between "{}" and "{}"'.format(
                    *k
                    )
                )

            sorted_atom_pairs = sorted(
                atom_pairs,
                key = lambda t: t[2],
                )

            logger.debug(
                'Current atom pairs (sorted by distance):\n{}'.format(
                    '\n'.join([
                        (
                        '\tAtom1: {a1}\n'
                        '\tAtom2: {a2}\n'
                        '\tDistance: {d}'
                        ).format(
                        a1 = Labeller.format(ap[0]),
                        a2 = Labeller.format(ap[1]),
                        d = ap[2],
                        )
                        for ap in atom_pairs
                        ])
                    )
                )

            #

            use_this_atom_hash = {}

            for a1, a2, d in sorted_atom_pairs:

                use_this_pair = True

                # Check if either atom is marked for not using

                for a in [a1, a2]:

                    a_id = a.id_str()

                    use_this_atom = use_this_atom_hash.get(a_id, True)

                    #
                    # Flagged not to be used
                    #
                    if use_this_atom is False:

                        # TODO Add exception here for atoms that aren't carbon that link to "exotics"?
                        # Could be longer-range interactions that are required to keep geometry for heavier atoms

                        logger.debug(
                            'Ignoring atom pair:\n\t{a1}\n\t{a2}\n\t{m}'.format(
                                a1 = Labeller.format(a1),
                                a2 = Labeller.format(a2),
                                m = '(Atom {} is marked as already used.)'.format(a_id),
                                )
                            )

                        use_this_pair = False
                    #
                    # Previously seen, check we want to use it to be used again
                    # (even if it's not used in this pair).
                    #
                    elif (use_this_atom is None):
                        #
                        # some atoms e.g. water can have as many bonds as it likes to the same residue
                        # or at least until it starts to link to carbon atoms, then mark as done
                        # stops too many restraints between carbon atoms and e.g. waters
                        #
                        if 'C' in [a1.element.strip().upper(), a2.element.strip().upper()]:

                            # Flag this atom not to be used again
                            use_this_atom_hash[a_id] = False

                            # Set pair flag
                            use_this_pair = False

                            logger.debug(
                                'Ignoring atom pair:\n\t{a1}\n\t{a2}\n\t{m}'.format(
                                    a1 = Labeller.format(a1),
                                    a2 = Labeller.format(a2),
                                    m = '(Not making multiple restraints including C-X interactions.)'.format(a_id),
                                    )
                                )

                if (use_this_pair is False):
                    continue

                #

                logger.debug(
                    'Keeping atom pair:\n\t{a1}\n\t{a2}\n\t{m}'.format(
                        a1 = Labeller.format(a1),
                        a2 = Labeller.format(a2),
                        m = 'Distance = {}'.format(a1.distance(a2)),
                        )
                    )

                # Mark each atom are now used

                for a in [a1, a2]:

                    a_id = a.id_str()

                    # Assign check flag if unassigned
                    if a_id not in use_this_atom_hash:
                        #
                        # Mark as used, but able to be considered
                        #
                        use_this_atom_hash[a_id] = None
                        #

                yield (a1, a2)

            # makes log easier to read
            logger.debug('')

    def make_atom_pair_restraints(self, atom_pairs, label=None):

        distance_restraints = []

        for a1, a2 in list(atom_pairs):

            distance_restraints.append(
                DistanceRestraint(
                    atom1 = a1,
                    atom2 = a2,
                    length = a1.distance(a2),
                    sigma = self.distance_restraint_sigma,
                    )
                )

        return RestraintsCollection(
            label = label,
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



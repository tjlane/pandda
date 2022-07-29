import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy as np

from itertools import cycle

import iotbx.pdb
from scitbx.array_family import flex

from giant.common.geometry import (
    rmsd_coordinates,
    )

from ..formatting import (
    labeller,
    )

from ..common import (
    GetInterestingResnames,
    )

from ..occupancy import (
    ResetOccupancies,
    )

def rg_label(rg):
    return (rg.resseq_as_int(), rg.icode)

def hierarchy_to_residue_group_dict(hierarchy):

    h_dict = {}

    for ch in hierarchy.only_model().chains():

        ch_dict = h_dict.setdefault(
            ch.id, {},
            )

        for rg in ch.residue_groups():

            ch_dict.setdefault(
                rg_label(rg), [],
                ).append(
                rg
                )

    return h_dict


class ResolveHierarchyConflicts(object):
    """
    Move residues in mov_hierarchy to new residue numbers
    if they have the same resid as a residue in fixed_hierarchy
    but different resnames
    """

    def __init__(self, in_place=False):

        self.in_place = bool(in_place)

    def __call__(self, fixed_hierarchy, moving_hierarchy):

        if self.in_place is False:
            moving_hierarchy = moving_hierarchy.deep_copy()

        moving_hierarchy = self.resolve_residue_conflicts(
            fixed_hierarchy = fixed_hierarchy,
            moving_hierarchy = moving_hierarchy,
            )

        return moving_hierarchy

    def resolve_residue_conflicts(self, fixed_hierarchy, moving_hierarchy):

        # Sort all residues (by chain then id) for the fixed hierarchy - chains > residue_ids > residue_groups (objects)
        fixed_dict = hierarchy_to_residue_group_dict(fixed_hierarchy)

        # Find the residues with clashing resids
        residues_to_update = []
        #
        for ch_mov in moving_hierarchy.only_model().chains():

            for rg_mov in ch_mov.residue_groups():

                # Extract equivalent group for this residue
                rg_fixs = fixed_dict.get(
                    ch_mov.id, {},
                    ).get(
                    rg_label(rg_mov), None,
                    )

                if rg_fixs is None:
                    continue

                if len(rg_fixs) > 1:
                    # not sure how this arises...
                    err_str = (
                        "A residue is present more than once in the output hierarchy? (but not the same residue group?)\n"
                        "> Moving residue: \n\t{rg_mov_str}\n"
                        "> Target residues: \n\t{rg_ref_strs}\n"
                        "Each residue can only be present once unless labelled with appropriate alternate conformations.\n"
                        ).format(
                        rg_mov_str = (
                            '{label} ({resnames})'.format(
                                label = labeller(rg_mov),
                                resnames = '/'.join(rg_mov.unique_resnames()),
                                )
                            ),
                        rg_ref_strs = '\n\t'.join(
                            [
                                '{label} ({resnames})'.format(#
                                    label = labeller(r),
                                    resnames = '/'.join(r.unique_resnames()),
                                    )
                                for r in rg_fixs
                                ]
                            ),
                        )
                    raise Exception(err_str)

                # Get the single rg
                rg_fix = rg_fixs[0]

                # Check to see if the residue is the same type as in the reference structure
                if not set(rg_fix.unique_resnames()).symmetric_difference(rg_mov.unique_resnames()):

                    # TODO allow this if the resnames match OR if it's a protein residue? (allows mutations/reactions?)

                    continue

                else:

                    # Will need to be changed
                    logger.debug(
                        'Different residues with same id: {label}, {r1} & {r2}'.format(
                            label = labeller(rg_mov),
                            r1 = '/'.join(rg_fix.unique_resnames()),
                            r2 = '/'.join(rg_mov.unique_resnames()),
                            )
                        )

                    residues_to_update.append(rg_mov)

        # Nothing to do -- return
        if len(residues_to_update) == 0:
            return moving_hierarchy

        # Find the next unused residue number for each chain
        resnum_dict = self.get_next_unused_resnums(
            hierarchies = [fixed_hierarchy, moving_hierarchy],
            )

        # Go through and increment residue groups ressnums to unused values
        for rg_mov in residues_to_update:

            ch_id = rg_mov.parent().id

            # Get next resnum for the chain
            next_resnum = resnum_dict[ch_id]

            # Update the dict for any future rgs
            resnum_dict[ch_id] = (next_resnum + 1)

            rg_mov.resseq = '{:4d}'.format(next_resnum)
            rg_mov.icode = ' '

        return moving_hierarchy

    def get_next_unused_resnums(self, hierarchies):

        resnum_dict = {}

        for h in hierarchies:

            for ch in h.chains():

                resseqs = [
                    rg.resseq_as_int() for rg in ch.residue_groups()
                    ]

                resnum_dict[ch.id] = max(
                    resnum_dict.get(ch.id, 1),
                    max(resseqs) + 1,
                    )

        return resnum_dict


class ExpandToFullMultiConformer(object):

    def __init__(self, in_place=False):

        self.in_place = bool(in_place)

        self.protein_amino_acid_set = set(
            iotbx.pdb.common_residue_names_amino_acid +
            iotbx.pdb.common_residue_names_modified_amino_acid
            )

    def __call__(self, hierarchy):

        if self.in_place is False:
            hierarchy = hierarchy.deep_copy()

        hierarchy = self.expand(
            hierarchy = hierarchy,
            )

        return hierarchy

    def expand(self, hierarchy):

        full_altloc_set = sorted([
            a for a in hierarchy.altloc_indices() if a
            ])

        # If no altlocs, expand all to "A"
        if len(full_altloc_set) == 0:
            logger.debug(
                'No altlocs in structure: expanding all residues to conformer "A"'
                )
            full_altloc_set = ['A']

        logger.debug(
            'Expanding all (appropriate) residues to have altlocs {}'.format(
                str(full_altloc_set)
                )
            )

        # Iterate through and expand each residue group to have all conformers
        for chain in hierarchy.chains():

            for residue_group in chain.residue_groups():

                # If has conformers but has blank altloc atoms (add blank ag to all other ags)
                if (
                    residue_group.have_conformers() and
                    residue_group.move_blank_altloc_atom_groups_to_front()
                    ):

                    logger.debug(
                        '{label} - expanding to pure conformer (altlocs {altlocs})'.format(
                            label = labeller(residue_group),
                            altlocs = str([a.altloc for a in residue_group.atom_groups()]),
                            )
                        )

                    # Convert all residue_groups to pure alt-conf
                    self.convert_proper_alt_conf_to_pure_alt_conf(
                        residue_group = residue_group,
                        )

                # Can go to next if all conformers are present for this residue group
                current_set = {
                    a.altloc for a in residue_group.atom_groups()
                    }

                if not current_set.symmetric_difference(full_altloc_set):
                    continue

                # should check here if the whole chain is alternate conformers then should not expand?

                # Expand incomplete conformer sets
                if (not residue_group.have_conformers()) or self.protein_amino_acid_set.intersection(residue_group.unique_resnames()):

                    # Only want to expand conformers for protein atoms (which should be present in all conformers)
                    # or where the residue group is only present in one conformation (single conformer water)
                    # but DO NOT want to expand waters in conformer A to A,B,C etc...


                    logger.debug(
                        '{label} - populating missing conformers (current altlocs {current}, target set {target})'.format(
                            label = labeller(residue_group),
                            current = str(current_set),
                            target = str(full_altloc_set),
                            )
                        )

                    # Populate missing conformers (from the other conformers)
                    self.populate_missing_conformers(
                        residue_group = residue_group,
                        full_altloc_set = full_altloc_set,
                        )

                    assert full_altloc_set == sorted(
                        [a.altloc for a in residue_group.atom_groups()]
                        )

                    logger.debug(
                        '{label} - updated conformer list: (current altlocs {current}, target set {target})'.format(
                            label = labeller(residue_group),
                            current = str([a.altloc for a in residue_group.atom_groups()]),
                            target = str(full_altloc_set),
                            )
                        )

        return hierarchy

    def convert_proper_alt_conf_to_pure_alt_conf(self, residue_group):

        main_ags = [ag for ag in residue_group.atom_groups() if ag.altloc=='']
        conf_ags = [ag for ag in residue_group.atom_groups() if ag.altloc!='']

        assert len(main_ags)==1, "Must be one atom_group in residue_group with altloc==''"
        assert len(conf_ags)!=0, "Must be at least one alternate conformer present"

        # Atom group to be merged into the other ags
        main_ag = main_ags[0]

        # Remove from the residue group
        residue_group.remove_atom_group(main_ag)

        # Iterate through and add atoms from main_conf to each of the other atom_groups
        for conf_ag in conf_ags:

            new_main_ag = main_ag.detached_copy()

            # Set the occupancy of the main_conf atoms to be transferred
            max_occ = max(conf_ag.atoms().extract_occ())
            new_main_ag.atoms().set_occ(flex.double(new_main_ag.atoms().size(), max_occ))

            # Change the altloc - very important
            new_main_ag.altloc = conf_ag.altloc

            # Merge the atom groups
            residue_group.merge_atom_groups(
                primary = conf_ag,
                secondary = new_main_ag,
                )

        return residue_group

    def populate_missing_conformers(self, residue_group, full_altloc_set):
        """
        Expand the alternate conformations of the residue_group to full_set, using the existing conformations
        - no conformers are deleted if existing conformers are not in full_altloc_set
        """

        # pool_ags -> original multi-conf atom groups
        # pool_alts -> altlocs of the original ags
        # occupancy_corrections -> multipliers to apply to the ags based on how many times
        #                           they will be copied to maintain total occupancy

        # Create pool of atom groups to use for the new conformers
        if residue_group.have_conformers():

            # Use any atom_group with an altloc - sort in decreasing occupancy
            pool_ags = sorted(
                [ag for ag in residue_group.atom_groups() if ag.altloc],
                key = lambda x: max(x.atoms().extract_occ()),
                reverse = True,
                )

            pool_alts = [ag.altloc for ag in pool_ags]

            # Occupancy multipliers for each atom group dependant on how much they'll be duplicated
            n_current_alt = len(pool_alts)
            n_target_alt = (
                n_current_alt + len(set(full_altloc_set).difference(pool_alts))
                )
            occupancy_corrections = [
                1.0 / (
                    (n_target_alt//n_current_alt) +
                    ((n_target_alt%n_current_alt)>k)
                    )
                for k in range(n_current_alt)
                ]

        else:

            # Only one conformation
            pool_ags = [
                residue_group.only_atom_group()
                ]

            pool_alts = [] # none

            # Occupancy multipliers for each atom group dependant on how much they'll be duplicated
            occupancy_corrections = [
                1.0/len(full_altloc_set)
                ]

            # Remove the blank conformer (to add it again with different altlocs later)
            assert pool_ags[0].altloc == ''
            residue_group.remove_atom_group(pool_ags[0])

        logger.debug(
            '{label} - occupancy corrections: {occ_cor}'.format(
                label = labeller(residue_group),
                occ_cor = str(occupancy_corrections),
                )
            )

        # Apply occupancy multipliers
        for mult, ag in zip(occupancy_corrections,pool_ags):
            ag.atoms().set_occ(
                ag.atoms().extract_occ() * mult
                )

        # Create cycle to iterate through ags as needed
        ag_cycle = cycle(pool_ags)

        # Iterate through and create alternate conformers as required
        for altloc in full_altloc_set:

            # Check if altloc already present
            if altloc in pool_alts:
                continue

            # Create copy and add to residue_group
            new_ag = next(ag_cycle).detached_copy()
            new_ag.altloc = altloc
            residue_group.append_atom_group(
                new_ag
                )

        return residue_group


class _AltlocReassigner(object):

    def reassign_altlocs(self, hierarchy, altloc_hash):

        logger.debug(
            "Updating altlocs: \n\t{}".format(
                '\n\t'.join([
                    "{a} -> {b}".format(
                        a = k,
                        b = v,
                        )
                    for k,v in sorted(altloc_hash.items())
                    ])
                )
            )

        for atom_group in hierarchy.atom_groups():

            if atom_group.altloc == '':
                continue

            logger.debug(
                '{label} - updating altloc: {a} -> {b}'.format(
                    label = labeller(atom_group),
                    a = atom_group.altloc,
                    b = altloc_hash[atom_group.altloc],
                    )
                )

            atom_group.altloc = altloc_hash[atom_group.altloc]


class IncrementConformerAltlocs(_AltlocReassigner):

    def __init__(self, in_place=False):

        self.in_place = bool(in_place)

        self.conformer_ids = iotbx.pdb.systematic_chain_ids()

    def __call__(self, fixed_hierarchy, moving_hierarchy):

        if self.in_place is False:
            moving_hierarchy = moving_hierarchy.deep_copy()

        n_shift = self.get_next_conformer_idx(
            hierarchy = fixed_hierarchy,
            )

        new_altlocs_hash = self.get_altloc_hash(
            hierarchy = moving_hierarchy,
            i_start = n_shift,
            )

        self.reassign_altlocs(
            hierarchy = moving_hierarchy,
            altloc_hash = new_altlocs_hash,
            )

        return moving_hierarchy

    def get_next_conformer_idx(self, hierarchy):

        cur_altlocs = [
            a for a in hierarchy.altloc_indices() if a
            ]

        if len(cur_altlocs) == 0:
            raise Exception('Hierarchy must have alternate conformations')

        cur_altloc_idxs = [
            self.conformer_ids.index(c)
            for c in cur_altlocs
            ]

        return max(cur_altloc_idxs) + 1

    def get_altloc_hash(self, hierarchy, i_start):

        current_altlocs = sorted([a for a in hierarchy.altloc_indices() if a])

        new_altlocs = self.conformer_ids[i_start:i_start+len(current_altlocs)]

        altlocs_hash = dict(zip(current_altlocs, new_altlocs))

        return altlocs_hash


class SanitiseAltlocs(_AltlocReassigner):

    def __init__(self, in_place=False):

        self.in_place = bool(in_place)

        self.conformer_ids = iotbx.pdb.systematic_chain_ids()

    def __call__(self, hierarchy):

        if self.in_place is False:
            hierarchy = hierarchy.deep_copy()

        new_altlocs_hash = self.get_altloc_hash(
            hierarchy = hierarchy,
            )

        self.reassign_altlocs(
            hierarchy = hierarchy,
            altloc_hash = new_altlocs_hash,
            )

        return hierarchy

    def get_altloc_hash(self, hierarchy):

        current_altlocs = sorted([a for a in hierarchy.altloc_indices() if a])

        new_altlocs = self.conformer_ids[0:len(current_altlocs)]

        altlocs_hash = dict(zip(current_altlocs, new_altlocs))

        return altlocs_hash


class JoinHierarchies(object):
    """
    Transfer atom_groups from donor_hierarchy to matching
    residue_groups in acceptor_hierarchy, creating new chains
    and residue groups only where necessary.
    """

    def __init__(self, in_place=False):

        self.in_place = bool(in_place)

    def __call__(self, acceptor_hierarchy, donor_hierarchy):

        if self.in_place is False:
            acceptor_hierarchy = acceptor_hierarchy.deep_copy()

        acceptor_hierarchy = self.transfer_residue_groups(
            acceptor_hierarchy = acceptor_hierarchy,
            donor_hierarchy = donor_hierarchy,
            )

        return acceptor_hierarchy

    def transfer_residue_groups(self, acceptor_hierarchy, donor_hierarchy):

        acceptor_model = acceptor_hierarchy.only_model()

        acceptor_dict = hierarchy_to_residue_group_dict(
            hierarchy = acceptor_hierarchy,
            )

        # Dictionary to link matching chains (allows multiple chain As to be linked uniquely to multiple chain As)
        link_dict = {}

        # Residues that don't have a matching partner in the old hierarchy
        tricky_rgs = []

        # Iterate through donor chains
        for donor_ch in donor_hierarchy.only_model().chains():

            # If chain not in hierarchy, simply copy across
            if acceptor_dict.get(donor_ch.id, None) is None:

                logger.debug(
                    'Transferring whole chain to acceptor hierarchy: {label}'.format(
                        label = labeller(donor_ch),
                        )
                    )

                acceptor_model.append_chain(
                    donor_ch.detached_copy()
                    )

                continue

            # Chain present, copy by residue_group
            for donor_rg in donor_ch.residue_groups():

                # Find equivalent residue groups in the other hierarchy
                acceptor_rgs = acceptor_dict.get(
                    donor_ch.id
                    ).get(
                    rg_label(donor_rg), [],
                    )

                if len(acceptor_rgs) == 0:

                    # Have the possibility of multiple chains with the same id, so at the moment, store for later
                    tricky_rgs.append(donor_rg)

                    continue

                acceptor_rg = acceptor_rgs[0]

                if len(acceptor_rgs) > 1:

                    # Should only be one... just send warning...
                    logger.warning(
                        (
                            'More than one residue group in acceptor hierarchy '
                            'with the same residue_id and chain_id\n\t{labels}.\n'
                            'Adding to the first matching residue: {label1}.'
                            ).format(
                                labels = '\n\t'.join([
                                    labeller(rg)
                                    for rg in acceptor_rgs
                                    ]),
                                label1 = labeller(acceptor_rg),
                                )
                        )

                # Record the links between these chains (maybe useful later!)
                link_dict.setdefault(
                    donor_ch,
                    acceptor_rg.parent(),
                    )

                # Transfer atom groups to this residue_group
                logger.debug(
                    'Transferring atom groups: {label_donor} > {label_acceptor}'.format(
                        label_donor = labeller(donor_rg),
                        label_acceptor = labeller(acceptor_rg),
                        )
                    )

                for donor_ag in donor_rg.atom_groups():

                    acceptor_rg.append_atom_group(
                        donor_ag.detached_copy()
                        )

        # Transfer residues that have chain matches but don't have residue
        # matches in the acceptor structures. Sort by whether protein, etc,
        # but then just make a new chain.
        #
        for donor_rg in tricky_rgs:

            donor_ch = donor_rg.parent()
            donor_ch_id = donor_ch.id

            # Try to get chain from link_dict
            acceptor_ch = link_dict.get(
                donor_ch, # This actually uses the chain to hash
                None,
                )

            # If the chain isn't linked, find or make a new chain:
            if acceptor_ch is None:

                # If there's only one chain with the same ID, choose this one
                possible_chains = [
                    c for c in acceptor_model.chains()
                    if (c.id == donor_ch_id)
                    ]

                if len(possible_chains) == 1:

                    acceptor_ch = possible_chains[0]

            # Catch the rest here
            if acceptor_ch is None:

                acceptor_ch = iotbx.pdb.hierarchy.chain(
                    id = donor_ch_id,
                    )

                acceptor_model.append_chain(
                    acceptor_ch,
                    )

                link_dict.setdefault(
                    donor_ch,
                    acceptor_ch,
                    )

            # Append residue_group to chain
            logger.debug(
                'Transferring residue group {rg} to chain {ch} in acceptor hierarchy'.format(
                    rg = labeller(donor_rg),
                    ch = labeller(acceptor_ch),
                    )
                )

            acceptor_ch.append_residue_group(
                donor_rg.detached_copy()
                )

        return acceptor_hierarchy


class PruneRedundantConformers(object):
    """
    Remove alternate conformers of residues if residue has conformers
    of required_altlocs and all conformers are within rmsd_cutoff
    """

    def __init__(self,
        rmsd_cutoff = 0.05,
        in_place = True,
        ):

        self.rmsd_cutoff = float(rmsd_cutoff)

        self.in_place = bool(in_place)

    def __call__(self, hierarchy, required_altlocs=None):

        if self.in_place is False:
            hierarchy = hierarchy.deep_copy()

        # This is required
        hierarchy.sort_atoms_in_place()

        if required_altlocs is None:
            required_altlocs = [
                a for a in hierarchy.altloc_indices() if a.strip()
                ]

        hierarchy = self.prune(
            hierarchy = hierarchy,
            required_altlocs = required_altlocs,
            rmsd_cutoff = self.rmsd_cutoff,
            )

        return hierarchy

    def prune(self, hierarchy, required_altlocs, rmsd_cutoff):

        required_altlocs = set(required_altlocs)

        for residue_group in hierarchy.residue_groups():

            # Skip if no conformers
            if not residue_group.have_conformers():
                continue

            # Get the blank and non-blank altloc atom_groups
            if residue_group.move_blank_altloc_atom_groups_to_front() != 0:

                main_ag = residue_group.atom_groups()[0]
                alt_ags = residue_group.atom_groups()[1:]

                assert main_ag.altloc == ''
                assert len(alt_ags) > 0

            else:

                main_ag = None
                alt_ags = residue_group.atom_groups()

            # Check no misplaced main conf
            assert '' not in [ag.altloc for ag in alt_ags]

            # Check if ALL required altlocs are present (skip if not)
            if bool(
                required_altlocs.difference(
                    [ag.altloc for ag in alt_ags]
                    )
                ):
                continue

            # Check if all pair of conformers are within rmsd cutoff
            prune = True
            #
            for i, ag_1 in enumerate(alt_ags):

                atoms_1 = ag_1.atoms()

                for j, ag_2 in enumerate(alt_ags):

                    if j<=i:
                        continue

                    atoms_2 = ag_2.atoms()

                    if not self.are_mergable(atoms_1=atoms_1, atoms_2=atoms_2):
                        prune = False
                        break

                    rmsd = rmsd_coordinates(
                        atoms_1.extract_xyz(),
                        atoms_2.extract_xyz(),
                        )

                    logger.debug(
                        'Residue {label}, alt {a} - alt {b}: rmsd {rmsd}'.format(
                            label = labeller.format(residue_group),
                            a = i,
                            b = j,
                            rmsd = rmsd,
                            )
                        )

                    if (rmsd > rmsd_cutoff):

                        logger.debug(
                            '> Not pruning {label}'.format(
                                label = labeller(residue_group),
                                )
                            )

                        prune = False

                        break

                if prune is False:
                    break

            if prune is False:
                continue

            # All rmsds below cutoff - prune!
            #
            logger.debug(
                '> Pruning {label}: altlocs {alts} -> [""]'.format(
                    label = labeller(residue_group),
                    alts = str([ag.altloc for ag in alt_ags]),
                    )
                )
            #
            if (main_ag is not None):

                # Merge one alt group with the main atom_group

                new_main_ag = alt_ags[0].detached_copy()
                new_main_ag.altloc = ''

                self.rescale_occupancies(
                    atoms = new_main_ag.atoms(),
                    max_occ = max(main_ag.atoms().extract_occ()),
                    )

                residue_group.merge_atom_groups(
                    primary = main_ag,
                    secondary = new_main_ag,
                    )

            else:

                # Remove one atom_group and set altloc to ''

                new_main_ag = alt_ags.pop(0)
                new_main_ag.altloc = ''

                self.rescale_occupancies(
                    atoms = new_main_ag.atoms(),
                    max_occ = sum(
                        [max(ag.atoms().extract_occ()) for ag in [new_main_ag]+alt_ags]
                        )
                    )

            # Remove all remaining alternate groups
            [residue_group.remove_atom_group(ag) for ag in alt_ags]

            assert len(residue_group.atom_groups())==1

        return hierarchy

    @staticmethod
    def are_mergable(atoms_1, atoms_2):

        if atoms_1.size() != atoms_2.size():
            return False

        if not atoms_1.extract_name().all_eq(atoms_2.extract_name()):
            return False

        return True

    @staticmethod
    def rescale_occupancies(atoms, max_occ=1.0):
        """Scale the maximum occupancy of a group of atoms to max_occ"""

        occ = atoms.extract_occ()

        assert min(occ) >= 0.0, 'occupancies cannot be negative!'

        occ_mult = max_occ/max(occ)#

        atoms.set_occ(occ*occ_mult)

        return atoms


class MakeMultiStateModel(object):

    def __init__(self,
        prune_rmsd_cutoff = 0.05,
        in_place = False,
        ):

        self.in_place = bool(in_place)

        self.resolve_hierarchy_conflicts = ResolveHierarchyConflicts(
            in_place = True,
            )

        self.expand_to_full_multi_conformer = ExpandToFullMultiConformer(
            in_place = True,
            )

        self.increment_altlocs = IncrementConformerAltlocs(
            in_place = True,
            )

        self.join_hierarchies = JoinHierarchies(
            in_place = True,
            )

        self.prune_redundant_conformers = PruneRedundantConformers(
            rmsd_cutoff = prune_rmsd_cutoff,
            in_place = True,
            )

    def __call__(self,
        hierarchies,
        ):

        # Need to create separate list as will be popping elements
        hierarchies = list(hierarchies)

        if self.in_place is False:
            hierarchies = [
                h.deep_copy() for h in hierarchies
                ]

        logger(
            '* Preparing input structures *'
            )

        for h in hierarchies:
            h.sort_atoms_in_place()

        main_hierarchy = hierarchies.pop(0)

        logger(
            (
                'Taken the first hierarchy as the main hierarchy\n'
                '{n} other models provided.'
                ).format(
                n = len(hierarchies),
                )
            )

        # main_hierarchy  = fixed / acceptor
        # other_hierarchy = moving / donor

        logger(
            '* Expanding main hierarchy to multi-conformer model *'
            )

        self.expand_to_full_multi_conformer(
            hierarchy = main_hierarchy,
            )

        for i_h, other_hierarchy in enumerate(hierarchies):

            logger(
                '* Merging model {} of {} into main hierarchy *'.format(
                    i_h+1, len(hierarchies),
                    )
                )

            logger(
                'Expanding to multi-conformer model'
                )

            self.expand_to_full_multi_conformer(
                hierarchy = other_hierarchy,
                )

            logger(
                'Resolving any hierarchy residue conflicts'
                )

            self.resolve_hierarchy_conflicts(
                fixed_hierarchy = main_hierarchy,
                moving_hierarchy = other_hierarchy,
                )

            logger(
                'Incrementing altlocs of the secondary hierarchy'
                )

            self.increment_altlocs(
                fixed_hierarchy = main_hierarchy,
                moving_hierarchy = other_hierarchy,
                )

            logger(
                'Merging secondary hierarchy into primary hierarchy'
                )

            self.join_hierarchies(
                acceptor_hierarchy = main_hierarchy,
                donor_hierarchy = other_hierarchy,
                )

        logger(
            '* Pruning redundant conformations from the output hierarchy *'
            )

        self.prune_redundant_conformers(
            main_hierarchy,
            )

        return main_hierarchy


class SplitHierarchyByConformer(object):

    def __init__(self):

        pass

    def __call__(self, hierarchy):

        output_dict = {}

        for conf_id in sorted(hierarchy.altloc_indices()):

            # Skip main conf
            if not conf_id.strip():
                continue

            label = str(conf_id)

            output_dict[label] = self.select(
                hierarchy = hierarchy,
                conformers = [conf_id],
                )

        return output_dict

    def select(self, hierarchy, conformers):

        ac = hierarchy.atom_selection_cache()

        sel_string = ' or '.join(
            [
                'altid "{}"'.format(c)
                for c in sorted(set([' ']+list(conformers)))
                ]
            )

        sel_bool = ac.selection(sel_string)

        sel_hierarchy = hierarchy.select(
            sel_bool,
            copy_atoms = True,
            )

        return sel_hierarchy


class SplitHierarchyByConformerGroup(SplitHierarchyByConformer):

    name = "SplitHierarchyByConformerGroup"

    def __init__(self, conformer_id_sets):

        self.conformer_id_sets = [list(s) for s in conformer_id_sets]

    def __str__(self):

        s_ = (
            'Task: {name}\n'
            '| conformer_id_sets: \n'
            '| \t{conformer_id_sets}\n'
            '`---->'
            ).format(
                name = self.name,
                conformer_id_sets = '\n'.join(
                    map(str,self.conformer_id_sets)
                    ).replace('\n','\n| \t'),
                )

        return s_.strip()

    def __call__(self, hierarchy):

        logger(str(self))

        output_dict = {}

        for conf_ids in self.conformer_id_sets:

            label = ''.join(sorted(conf_ids))

            output_dict[label] = self.select(
                hierarchy = hierarchy,
                conformers = conf_ids,
                )

        return output_dict


class SplitHierarchyByResidueNames(SplitHierarchyByConformer):

    name = "SplitHierarchyByResidueNames"

    def __init__(self,
        ignore_common_molecules = True,
        include_resnames_list = None,
        ignore_resnames_list = None,
        atom_selection = None,
        selected_name = 'selected',
        unselected_name = 'unselected',
        combine_split_states = False,
        ):

        self.get_interesting_resnames = GetInterestingResnames(
            ignore_common_molecules = ignore_common_molecules,
            include_resnames_list = include_resnames_list,
            ignore_resnames_list = ignore_resnames_list,
            )

        self.atom_selection = (
            atom_selection
            if atom_selection is not None
            else None
            )

        self.selected_name = (
            str(selected_name)
            )

        self.unselected_name = (
            str(unselected_name)
            )

        self.combine_split_states = bool(
            combine_split_states
            )

        assert (
            self.selected_name != self.unselected_name
            ), "Selected and unselected names cannot be the same"

    def __str__(self):

        s_ = (
            'Task: {name}\n'
            '| get_interesting_resnames: \n'
            '| \t{get_interesting_resnames}\n'
            '| atom_selection: {atom_selection}\n'
            '| selected_name: {selected_name}\n'
            '| unselected_name: {unselected_name}\n'
            '`---->'
            ).format(
                name = self.name,
                get_interesting_resnames = str(self.get_interesting_resnames).replace('\n','\n| \t'),
                atom_selection = str(self.atom_selection),
                selected_name = str(self.selected_name),
                unselected_name = str(self.unselected_name),
                )

        return s_.strip()

    def __call__(self, hierarchy):

        logger(str(self))

        altloc_dict = self.get_conformer_dict(
            hierarchy = hierarchy,
            )

        logger(
            '\n'+'\n'.join([
                "{k} : {conf_ids}".format(
                    k = k,
                    conf_ids = str(conf_ids),
                    )
                for k, conf_ids in sorted(altloc_dict.items())
                ])
            )

        output_dict = {}

        for conf_key, conf_ids in sorted(altloc_dict.items()):

            output_dict[conf_key] = self.select(
                hierarchy = hierarchy,
                conformers = conf_ids,
                )

        return output_dict

    def get_conformer_dict(self, hierarchy):

        altlocs_dict = {}

        if self.atom_selection is not None:

            sel_altlocs_dict = self.get_altlocs_dict_from_selected_atoms(
                hierarchy = hierarchy,
                )

        else:

            sel_altlocs_dict = self.get_altlocs_dict_from_interesting_resnames(
                hierarchy = hierarchy,
                )

        altlocs_dict.update(
            sel_altlocs_dict
            )

        other_altlocs = self.get_remaining_altlocs_dict(
            hierarchy = hierarchy,
            altlocs = (
                sorted(set(list(
                    np.concatenate(list(sel_altlocs_dict.values()))
                    )))
                if sel_altlocs_dict else []
                ),
            )

        altlocs_dict.update(
            other_altlocs
            )

        return altlocs_dict

    def get_altlocs_dict_from_selected_atoms(self, hierarchy):

        assert self.atom_selection is not None

        sel_altlocs_groups = self.get_selection_altlocs_groups_by_residue(
            hierarchy = hierarchy,
            atom_selection = self.atom_selection,
            )

        if len(sel_altlocs_groups) == 0:
            logger.warning(
                (
                    'No altlocs selected in hierarchy by residue name.\n'
                    '{task}\n'
                    'Atom selection: {selection}.'
                    ).format(
                    task = str(self),
                    selection = str(self.atom_selection),
                    )
                )

        return {
            self.selected_name+'-'+''.join(sel_a) : sel_a
            for sel_a in sel_altlocs_groups
            }

    def get_altlocs_dict_from_interesting_resnames(self, hierarchy):

        resnames = self.get_interesting_resnames(hierarchy)

        logger.debug('Selected resnames: {s}'.format(s=str(resnames)))

        atom_selection = ' or '.join(
            ['resname "{r}"'.format(r=r) for r in resnames]
            )

        sel_altlocs_groups = self.get_selection_altlocs_groups_by_residue(
            hierarchy = hierarchy,
            atom_selection = atom_selection,
            )

        if len(sel_altlocs_groups) == 0:
            logger.warning(
                (
                    '{task}\n'
                    '** No altlocs selected in hierarchy by residue name. **\n'
                    '** Identified residue names: {selection}. **'
                    ).format(
                    task = str(self),
                    selection = str(atom_selection),
                    )
                )

        return {
            self.selected_name+'-'+''.join(sel_a) : sel_a
            for sel_a in sel_altlocs_groups
            }

    def get_remaining_altlocs_dict(self, hierarchy, altlocs):

        altlocs2 = set(hierarchy.altloc_indices())
        altlocs2.discard('')
        altlocs2.difference_update(altlocs)
        altlocs2 = tuple(sorted(altlocs2))

        return {
            self.unselected_name+'-'+''.join(altlocs2) : altlocs2,
            }

    def get_selection_altlocs_groups(self, hierarchy, atom_selection):

        asc = hierarchy.atom_selection_cache()
        selection = asc.selection(atom_selection)
        hierarchy = hierarchy.select(selection)

        sel_altlocs = sorted([
            a for a in hierarchy.altloc_indices() if a.strip()
            ])

        return [sel_altlocs]

    def get_selection_altlocs_groups_by_residue(self, hierarchy, atom_selection):

        asc = hierarchy.atom_selection_cache()
        selection = asc.selection(atom_selection)
        hierarchy = hierarchy.select(selection)

        sel_altlocs_groups = set()

        for rg in hierarchy.residue_groups():

            n_blank = rg.move_blank_altloc_atom_groups_to_front()

            # if n_blank == 0:
            #     continue

            altlocs = [
                ag.altloc for ag in rg.atom_groups()[n_blank:]
                ]

            logger(
                "Residues {resnames} : Conformers {conformers}".format(
                    resnames = '/'.join(rg.unique_resnames()),
                    conformers = ','.join(altlocs),
                    )
                )

            sel_altlocs_groups.add(
                tuple(sorted(altlocs))
                )

        if (self.combine_split_states is True):

            logger(
                'Pooling all conformers into one group: \n\t{string}'.format(
                    string = '\n\t'.join([
                        "{v}".format(
                            v = str(v),
                            )
                        for v in sorted(sel_altlocs_groups)
                        ])
                    )
                )

            combined_set = set()
            for g in sel_altlocs_groups:
                combined_set.update(g)

            sel_altlocs_groups = set()
            sel_altlocs_groups.add(
                tuple(sorted(combined_set))
                )

            logger(
                'Output groups: {}'.format(
                    str(sorted(sel_altlocs_groups))
                    )
                )

        return sorted(sel_altlocs_groups)


class SplitMultiStateModel(object):

    def __init__(self,
        split_hierarchy,
        prune_duplicates_rmsd = 0.05,
        reset_altlocs = False,
        reset_occupancies = False,
        ):

        self.split_hierarchy = (
            split_hierarchy
            )

        self.prune_redundant_conformers = (
            PruneRedundantConformers(
                rmsd_cutoff = prune_duplicates_rmsd,
                in_place = True,
                )
            if (prune_duplicates_rmsd is not None)
            else None
            )

        self.sanitise_altlocs = (
            SanitiseAltlocs(
                in_place = True,
                )
            if (reset_altlocs is True)
            else None
            )

        self.reset_occupancies = (
            ResetOccupancies(
                in_place = True,
                )
            if (reset_occupancies is True)
            else None
            )

    def __call__(self, hierarchy):

        hierarchy_dict = self.split_hierarchy(hierarchy)

        logger(
            (
                "\n"
                "Split hierarchy:\n"
                "\t{output}\n"
                ).format(
                    output = '\n\t'.join([
                        (
                            "Hierarchy label: {name}\n"
                            "Altlocs: {altlocs}"
                            ).format(
                            name = k,
                            altlocs = str(list(h.altloc_indices()))
                            ).replace('\n','\n\t')
                        for k, h in sorted(hierarchy_dict.items())
                        ])
                )
            )

        for key, split_h in hierarchy_dict.items():

            if self.prune_redundant_conformers is not None:
                self.prune_redundant_conformers(split_h)

            if self.sanitise_altlocs is not None:
                self.sanitise_altlocs(split_h)

            if self.reset_occupancies is not None:
                self.reset_occupancies(split_h)

        return hierarchy_dict



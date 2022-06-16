import giant.logs as lg
logger = lg.getLogger(__name__)

import collections
import numpy as np

from scitbx.array_family import flex

from giant.common.geometry import (
    pairwise_distances,
    )

def occupancy_weighted_average_b_factor(atoms):
    from scitbx.array_family import flex
    return flex.mean_weighted(atoms.extract_b(), atoms.extract_occ())

BFactorRatioInfo = collections.namedtuple(
    typename = 'BFactorRatioInfo',
    field_names = [
        'selected_atoms',
        'selected_average_b',
        'neighbours_atoms',
        'neighbours_average_b',
        'b_factor_ratio',
    ],
)


class CalculateSurroundingsBFactorRatio(object):

    def __init__(self,
        distance_cutoff = 6.0,
        include_het_atoms_in_surroundings = False,
        ):

        self.distance_cutoff = float(distance_cutoff)

        self.include_het_atoms_in_surroundings = bool(include_het_atoms_in_surroundings)

        assert self.distance_cutoff > 0.

    def __call__(self,
        selected_atoms,
        hierarchy,
        ):

        all_atoms = hierarchy.atoms()
        sel_atoms = selected_atoms

        if (self.include_het_atoms_in_surroundings is False):
            all_atoms = all_atoms.select(
                hierarchy.atom_selection_cache().selection('not hetatm')
                )

        neighbour_selection = self.find_neighbours(
            sel_atoms,
            all_atoms,
            )

        neighbour_atoms = (
            all_atoms.select(neighbour_selection)
            )

        selected_b = occupancy_weighted_average_b_factor(
            sel_atoms,
            )

        neighbour_b = occupancy_weighted_average_b_factor(
            neighbour_atoms,
            )

        b_factor_ratio = None

        if neighbour_b != 0.0:

            b_factor_ratio = (
                selected_b / neighbour_b
                )

        return BFactorRatioInfo(
            selected_atoms = selected_atoms,
            selected_average_b = selected_b,
            neighbours_atoms = neighbour_atoms,
            neighbours_average_b = neighbour_b,
            b_factor_ratio = b_factor_ratio,
        )

    def find_neighbours(self, selected_atoms, all_atoms):

        sel_xyz = np.array(selected_atoms.extract_xyz())
        all_xyz = np.array(all_atoms.extract_xyz())

        distances = pairwise_distances(
            points_1 = sel_xyz,
            points_2 = all_xyz,
            )

        within_cutoff = (distances < self.distance_cutoff)

        self_selection = (distances == 0.0)
        within_cutoff[self_selection] = False

        # flatten to 1d selection
        neighbour_selection = (
            within_cutoff.any(axis=0)
            )

        assert len(neighbour_selection) == len(all_atoms)

        return flex.bool(neighbour_selection)



import giant.logs as lg
logger = lg.getLogger(__name__)

from scitbx.array_family import flex

def set_occupancy(atoms, occupancy):

    atoms.set_occ(
        flex.double(atoms.size(), occupancy)
        )


class CalculateResidueGroupOccupancy(object):

    def __init__(self):

        pass

    def __call__(self, residue_group):

        main_occ = self.get_main_conf_occupancy(residue_group)

        if main_occ != -1:

            return main_occ

        alt_occ = self.get_alt_conf_occupancy(residue_group)

        return alt_occ

    @staticmethod
    def get_main_conf_occupancy(residue_group):

        n_blank = residue_group.move_blank_altloc_atom_groups_to_front()

        if n_blank == 0:
            return -1

        main_ag_occs = [
            max(
                ag.atoms().extract_occ()
                )
            for ag in residue_group.atom_groups()[:n_blank]
            ]

        return max(main_ag_occs)

    @staticmethod
    def get_alt_conf_occupancy(residue_group):

        n_blank = residue_group.move_blank_altloc_atom_groups_to_front()

        alt_ag_occs = [
            max(
                ag.atoms().extract_occ()
                )
            for ag in residue_group.atom_groups()[n_blank:]
            ]

        if len(alt_ag_occs) == 0:
            return -1

        return sum(alt_ag_occs)


class ScaleOccupancies(object):

    def __init__(self, in_place=True):

        self.in_place = in_place

    def __call__(self, atoms, multiplier):

        assert multiplier >= 0.0

        if self.in_place is False:
            atoms = atoms.deep_copy()

        atoms.set_occ(
            atoms.extract_occ() * float(multiplier)
        )

        return atoms


class ResetOccupancies(object):

    def __init__(self,
        default_occupancy = 1.0,
        ):

        self.default_occupancy = float(default_occupancy)

    def __call__(self, hierarchy):

        altloc_indices = list(hierarchy.altloc_indices())

        conformer_indices = hierarchy.get_conformer_indices()

        n_altlocs = len(
            [a for a in altloc_indices if a.strip()]
            )

        atoms = hierarchy.atoms()

        ##

        if '' in altloc_indices:

            i_main = altloc_indices.index('')

            main_sel = (conformer_indices == i_main)
            alt_sel = (conformer_indices != i_main)

        else:

            main_sel = None
            alt_sel = flex.bool(atoms.size(), True)

        ##

        if main_sel is not None:

            main_atoms = atoms.select(main_sel)

            set_occupancy(
                atoms = main_atoms,
                occupancy = self.default_occupancy,
                )

        if n_altlocs > 0:

            alt_atoms = atoms.select(alt_sel)

            alt_occupancy = (
                self.default_occupancy / float(n_altlocs)
                )

            set_occupancy(
                atoms = alt_atoms,
                occupancy = alt_occupancy,
                )

        return hierarchy


class SanitiseOccupancies(object):

    def __init__(self,
        min_occupancy = 0.0,
        max_occupancy = 1.0,
        ):

        self.min_occupancy = float(min_occupancy)
        self.max_occupancy = float(max_occupancy)

        self.scale_occupancy = ScaleOccupancies()

    def __call__(self, hierarchy):

        self.scale_main_to_max(hierarchy)

        self.enforce_limits(hierarchy)

        return hierarchy

    def scale_main_to_max(self, hierarchy):

        altloc_indices = list(hierarchy.altloc_indices())

        if '' not in altloc_indices:
            return

        i_main = altloc_indices.index('')

        main_sel = (hierarchy.get_conformer_indices() == i_main)

        all_atoms = hierarchy.atoms()
        main_atoms = all_atoms.select(main_sel)

        max_main_occ = max(main_atoms.extract_occ())

        scale_occ = (
            self.max_occupancy / max_main_occ
            )

        logger.debug(
            (
                'Current maximum occupancy {cur_max}\n'
                'Scaling occupancies by factor of {scale} to {new_max}'
                ).format(
                cur_max = max_main_occ,
                scale = scale_occ,
                new_max = self.max_occupancy,
                )
            )

        self.scale_occupancy(
            atoms = all_atoms,
            multiplier = scale_occ,
            )

    def enforce_limits(self, hierarchy):

        atoms = hierarchy.atoms()
        occupancies = atoms.extract_occ()

        sel_min = (occupancies < self.min_occupancy)
        sel_max = (occupancies > self.max_occupancy)

        if not sel_min.all_eq(False):
            set_occupancy(
                atoms = atoms.select(sel_min),
                occupancy = self.min_occupancy,
                )

        if not sel_max.all_eq(False):
            set_occupancy(
                atoms = atoms.select(sel_max),
                occupancy = self.max_occupancy,
                )


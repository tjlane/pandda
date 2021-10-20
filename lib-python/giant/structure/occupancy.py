import giant.logs as lg
logger = lg.getLogger(__name__)

from scitbx.array_family import flex
from giant.structure.common import split_main_and_alt_conf_atom_groups
from giant.structure.formatting import Labeller

def set_occupancy(atoms, occupancy):

    atoms.set_occ(
        flex.double(atoms.size(), occupancy)
        )


class ScaleOccupancies:

    def __init__(self, in_place=True):

        self.in_place = in_place

    def __call__(self, atoms, multiplier):

        assert multiplier > 0.0

        if self.in_place is False:
            atoms = atoms.deep_copy()

        atoms.set_occ(
            atoms.extract_occ() * float(multiplier)
        )

        return atoms


class ResetOccupancies:

    def __init__(self,
        default_occupancy = 1.0,
        ):

        self.default_occupancy = float(default_occupancy)

    def __call__(self, hierarchy):

        altloc_indices = list(hierarchy.altloc_indices())

        n_altlocs = len(
            [a for a in altloc_indices if a.strip()]
            )

        i_main = altloc_indices.index('')

        conformer_indices = hierarchy.get_conformer_indices()

        main_sel = (conformer_indices == i_main)
        alt_sel = (conformer_indices != i_main)

        atoms = hierarchy.atoms()

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


class SanitiseOccupancies:

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


def calculate_residue_group_occupancy(residue_group):
    """
    Extract the total occupancy of a residue, allowing for alternate conformers.
    - If conformers are present, the occupancies of each conformer are summed.
    - If multiple occupancies are present for a conformer, the maximum is taken for each.
    """
    rg = residue_group
    # Extract main and alt confs
    main_ag, alt_ags = split_main_and_alt_conf_atom_groups(rg)
    # Extract main-conf occupancy
    if main_ag is not None:
        main_occ = max(main_ag.atoms().extract_occ())
    else:
        main_occ = 0.0
    # Extract alt-conf atom groups
    if alt_ags is not None:
        ag_occs = [max(ag.atoms().extract_occ()) for ag in alt_ags]
        alt_occ = sum(ag_occs)
    else:
        alt_occ = 0.0
    # Residue occupancy is max of both
    res_occ = max([main_occ, alt_occ])
    return res_occ



def sanitise_occupancies(hierarchy, fixed_conformers=None, min_occ=0.0, max_occ=1.0, in_place=False, verbose=False):
    """Sanitise the occupancies of a hierarchy so that the occupancies for a residue sum to 1.0"""
    assert (min_occ >= 0.0) and (max_occ <= 1.0)
    if fixed_conformers is None: fixed_conformers = []
    if (not in_place): hierarchy = hierarchy.deep_copy()
    # Iterate through the output structure, and normalise the occupancies if necessary
    for rg in hierarchy.residue_groups():
        # Calculate occupancy of the residue group
        rg_occ = calculate_residue_group_occupancy(residue_group=rg)
        # If occupancy in range, continue
        if min_occ <= rg_occ <= max_occ:
            continue
        if verbose: print 'Occupancy of residue {} is {} -- sanitising'.format(Labeller.format(rg), rg_occ)
        # Extract main-conf and alt-conf ags
        main_ag, alt_ags = split_main_and_alt_conf_atom_groups(rg)
        # Sanitise main conf
        if main_ag is not None:
            if verbose:
                print '------------------>'
                print 'Sanitising main-conf atom group:\n\t{}'.format(Labeller.format(main_ag))
                print 'Current occupancy: {}'.format(max(main_ag.atoms().extract_occ()))
            sanitise_atom_group_occupancies_in_place(main_ag, min_occ=min_occ, max_occ=max_occ)
            if verbose:
                print 'New occupancy:     {}'.format(max(main_ag.atoms().extract_occ()))
        # Sanitise alt confs
        if alt_ags is not None:
            # Get the groups to change and the groups to keep constant
            ag_chnge = [ag for ag in alt_ags if (ag.altloc not in fixed_conformers)]
            ag_const = [ag for ag in alt_ags if (ag.altloc in fixed_conformers)]
            # Calculate the total occupancy of the groups
            occ_chnge = sum([max(ag.atoms().extract_occ()) for ag in ag_chnge])
            occ_const = sum([max(ag.atoms().extract_occ()) for ag in ag_const])
            if occ_const > max_occ: raise Exception('Occupancy of fixed atom groups ({}) is already greater than maximum ({})'.format(occ_const, max_occ))
            # Normalise the occupancies of the changing groups
            if verbose:
                print '------------------>'
                print 'Sanitising alt-conf atom groups:'
                print 'Fixed conformers:\n\t{}'.format('\n\t'.join([Labeller.format(ag) for ag in ag_const]))
                print 'Other conformers:\n\t{}'.format('\n\t'.join([Labeller.format(ag) for ag in ag_chnge]))
                print 'Total occupancy (fixed): {}'.format(occ_const)
                print 'Individual occupancies:  {}'.format(', '.join(map(str,[max(ag.atoms().extract_occ()) for ag in ag_const])))
                print 'Total occupancy (other): {}'.format(occ_chnge)
                print 'Individual occupancies:  {}'.format(', '.join(map(str,[max(ag.atoms().extract_occ()) for ag in ag_chnge])))
            sanitise_multiple_atom_group_occupancies_in_place(atom_groups = ag_chnge,
                                                              min_occ = max(0.0, min_occ-occ_const),
                                                              max_occ = min(1.0, max_occ-occ_const))
            if verbose:
                print 'New total occupancy (other): {}'.format(sum([max(ag.atoms().extract_occ()) for ag in ag_chnge]))
                print 'Individual occupancies:  {}'.format(', '.join(map(str,[max(ag.atoms().extract_occ()) for ag in ag_chnge])))

    return hierarchy

def sanitise_atom_group_occupancies_in_place(atom_group, min_occ=0.0, max_occ=1.0):
    """Normalise an atom-group so total occupancy between min_occ & max_occ"""
    ag = atom_group
    ag_occ = max(ag.atoms().extract_occ())
    if ag_occ == 0.0:
        ag.atoms().set_occ(flex.double(ag.atoms().size(), min_occ))
    elif ag_occ < min_occ:
        normalise_occupancies(ag.atoms(), max_occ=min_occ)
    elif ag_occ > max_occ:
        normalise_occupancies(ag.atoms(), max_occ=max_occ)
    return

def sanitise_multiple_atom_group_occupancies_in_place(atom_groups, min_occ=0.0, max_occ=1.0):
    """Normalise a group of atom-groups so total occupancy between min_occ & max_occ"""
    tot_occ = sum([max(ag.atoms().extract_occ()) for ag in atom_groups])
    # If all are zero -- divide min_occ between them
    if tot_occ <= 0.0:
        div_occ = min_occ/len(atom_groups)
        for ag in atom_groups:
            ag.atoms().set_occ(flex.double(ag.atoms().size(), div_occ))
        return
    # Too large or too small -- calculate scaling factor
    if tot_occ > max_occ:
        occ_mult = max_occ/tot_occ
    elif tot_occ < min_occ:
        occ_mult = min_occ/tot_occ
    else:
        return
    for ag in atom_groups:
        ag.atoms().set_occ(occ_mult*ag.atoms().extract_occ())
    return

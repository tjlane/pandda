import giant.logs as lg
logger = lg.getLogger(__name__)


class BfactorStatistics(object):
    """B-factor statstics for a structure"""

    def __init__(self,
        all_atoms = None,
        protein = None,
        backbone = None,
        sidechain = None,
        ):

        self.statistics = dict(
            all_atoms = all_atoms,
            protein   = protein,
            backbone  = backbone,
            sidechain = sidechain,
        )

    def to_z_scores(self, b_values, normalise_by='backbone'):
        """Convert a set of B-factors to Z-scores"""

        assert normalise_by in self.statistics.keys()

        stats = self.statistics.get(normalise_by)

        if stats.biased_standard_deviation == 0.0:
            # No variation - return list of zeroes
            return (b_values - b_values)

        return (b_values - stats.mean) / stats.biased_standard_deviation

    @classmethod
    def from_pdb(cls, pdb_file=None, pdb_hierarchy=None):
        """Calculate the b-factor statistics of a model"""

        assert [pdb_file, pdb_hierarchy].count(None)==1, 'Provide pdb_input OR pdb_hierarchy'
        if (pdb_file is not None):
            import iotbx.pdb
            pdb_hierarchy = iotbx.pdb.hierarchy.input(pdb_file).hierarchy

        cache = pdb_hierarchy.atom_selection_cache()

        from giant.structure.select import non_h, protein, backbone, sidechains

        all_b       = non_h(hierarchy=pdb_hierarchy, cache=cache, copy=True).atoms().extract_b()
        protein_b   = protein(hierarchy=pdb_hierarchy, cache=cache, copy=True).atoms().extract_b()
        backbone_b  = backbone(hierarchy=pdb_hierarchy, cache=cache, copy=True).atoms().extract_b()
        sidechain_b = sidechains(hierarchy=pdb_hierarchy, cache=cache, copy=True).atoms().extract_b()

        from scitbx.math import basic_statistics

        return cls(
            all_atoms = basic_statistics(all_b),
            protein = basic_statistics(protein_b),
            backbone = basic_statistics(backbone_b),
            sidechain = basic_statistics(sidechain_b),
        )

    def __str__(self):

        s = ""

        for k in ['all_atoms', 'protein', 'backbone', 'sidechain']:
            out = lg.ListStream()
            stats = self.statistics[k]
            stats.show(f=out)
            s += (
                '"{}" B-factor summary statistics'.format(k) +
                ':\n\n' +
                str(out) +
                '\n'
            )

        return s.strip()


def normalise_b_factors_to_z_scores(
    pdb_file = None,
    pdb_hierarchy = None,
    normalise_by = 'backbone',
    ):
    """Calculate the b-factor statistics of the model"""

    assert [pdb_file, pdb_hierarchy].count(None)==1,'Provide pdb_input OR pdb_hierarchy'

    if (pdb_hierarchy is None):
        import iotbx.pdb
        pdb_hierarchy = iotbx.pdb.hierarchy.input(pdb_file).hierarchy

    bfs = BfactorStatistics.from_pdb(
        pdb_hierarchy = pdb_hierarchy,
    )

    new_b = bfs.to_z_scores(
        b_values = pdb_hierarchy.atoms().extract_b(),
        normalise_by = normalise_by,
    )

    output_h = pdb_hierarchy.deep_copy()
    output_h.atoms().set_b(new_b)

    return output_h

def occupancy_weighted_average_b_factor(atoms):
    from scitbx.array_family import flex
    return flex.mean_weighted(atoms.extract_b(), atoms.extract_occ())

import collections
BFactorRatioInfo = collections.namedtuple(
    typename = 'BFactorRatioInfo',
    field_names = [
        'selection_name',
        'surroundings_names',
        'selection_av_bfactor',
        'surroundings_av_bfactor',
        'b_factor_ratio',
    ],
)

def calculate_residue_group_bfactor_ratio(
    residue_group,
    hierarchy,
    distance_cutoff = 4.0,
    ):
    """Calculate bfactor quality metrics of the residue to surrounding residues"""

    import iotbx.pdb

    from giant.maths.geometry import is_within
    from giant.structure.formatting import ShortLabeller

    rg_sel = residue_group

    # Check validity
    if len(rg_sel.unique_resnames()) != 1:
        raise ValueError('More than one residue name associated with residue group -- cannot process')

    if (distance_cutoff < 0.0):
        raise ValueError('distance_cutoff must be >= 0.0')

    # Extract rg_sel objects
    rg_sel_ags    = rg_sel.atom_groups()
    rg_sel_atoms  = rg_sel.atoms()
    rg_sel_coords = rg_sel_atoms.extract_xyz()

    # Select nearby atom_group for scoring
    near_ags = [
        ag for ag in hierarchy.atom_groups()
        if (
            is_within(
                distance_cutoff,
                rg_sel_coords,
                ag.atoms().extract_xyz(),
            )
            and (ag not in rg_sel_ags)
        )
    ]

    if (not near_ags):
        return None

    # Extract atoms from nearby groups
    near_ats = iotbx.pdb.hierarchy.af_shared_atom()
    [near_ats.extend(ag.detached_copy().atoms()) for ag in near_ags]

    # Extract the names of the selection & nearby residues
    selection_name = ShortLabeller.format(rg_sel)
    surround_names = sorted(set([ShortLabeller.format(ag) for ag in near_ags]))

    # Calculate B-factors of the residue & surroundings
    sel_mean_b = occupancy_weighted_average_b_factor(
        atoms = rg_sel_atoms,
    )
    sur_mean_b = occupancy_weighted_average_b_factor(
        atoms = near_ats,
    )

    # Store the ratio of the b-factors
    b_factor_ratio = None
    if sur_mean_b != 0.0:
        b_factor_ratio = sel_mean_b / sur_mean_b

    return BFactorRatioInfo(
        selection_name = selection_name,
        surroundings_names = surround_names,
        selection_av_bfactor = sel_mean_b,
        surroundings_av_bfactor = sur_mean_b,
        b_factor_ratio = b_factor_ratio,
    )


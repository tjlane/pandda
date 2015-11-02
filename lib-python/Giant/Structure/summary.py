import iotbx.pdb

from scitbx.math import basic_statistics

from Bamboo.Density.Edstats import edstats
from Giant.Structure.b_factors import b_factor_statistics, normalise_b_factors
from Giant.Structure.html import write_html_summary

class structureSummary(object):

    def __init__(self, pdb_input=None, pdb_hierarchy=None, pdb_file=None, mtz_file=None):
        """Summarise some characteristics of a structure"""

        assert [pdb_file, pdb_input, pdb_hierarchy].count(None)==2,'Provide pdb_file or pdb_input OR pdb_hierarchy'
        if pdb_file:  pdb_input     = iotbx.pdb.input(file_name = pdb_file)
        if pdb_input: pdb_hierarchy = pdb_input.construct_hierarchy()
        # Structure + Data
        self.mtz_file       = mtz_file
        self.pdb_input      = pdb_input
        self.pdb_hierarchy  = pdb_hierarchy
        # B-factors
        self.b_factors = b_factor_statistics(pdb_hierarchy=pdb_hierarchy)
        # Edstats
        if mtz_file: self.edstats = edstats(mtz_file=mtz_file, pdb_file=pdb_file)
        else:        self.edstats = None

    def normalise_b_factors(self, root=None):
        if root is None: root = self.pdb_hierarchy
        return normalise_b_factors(pdb_hierarchy=root, b_factor_summary=self.b_factors)

    def _get_atom_group_edstats_scores(self, atom_group):
        resname = atom_group.resname
        chain   = atom_group.parent().parent().id
        res_i   = atom_group.parent().resseq_as_int()
#        conf_id = atom_group.altloc if atom_group.altloc else ' '
        conf_id = ' '
        return self.edstats.scores[(resname, chain, res_i, conf_id)]

    def get_atom_group_summary(self, atom_group):
        return atomGroupSummary( atom_group = atom_group,
                                 edstats_scores = self._get_atom_group_edstats_scores(atom_group = atom_group),
                                 global_b_factor_statistics = self.b_factors )

    def to_html(self, fname):
        write_html_summary(fname=fname, atom_group_summaries=[self.get_atom_group_summary(ag) for ag in self.pdb_hierarchy.atom_groups()])

class atomGroupSummary(object):

    def __init__(self, atom_group, edstats_scores, global_b_factor_statistics):

        self.atom_group = atom_group
        self.edstats_scores = edstats_scores
        self.b_factors              = atom_group.atoms().extract_b()
        self.b_factors_statistics   = basic_statistics(atom_group.atoms().extract_b())
        self.b_factors_z            = global_b_factor_statistics.to_z_score(b_vals=atom_group.atoms().extract_b())
        self.b_factors_z_statistics = basic_statistics(self.b_factors_z)
#        self._global_b_factor_statistics = global_b_factor_statistics

    def to_html(self, fname):
        write_html_summary(fname=fname, atom_group_summaries=[self])

#class atomGroupSummary(object):
#
#    def __init__(self, atom_group):
#
#        self.atom_group = atom_group
#
#        # Meta about atom group
#        self.residue_class = iotbx.pdb.common_residue_names_get_class(atom_group.resname)
#
#        # Occupancy summary
#        occ_stats = atom_group.atoms().extract_occ().min_max_mean()
#        self.occ_min = occ_stats.min
#        self.occ_max = occ_stats.max
#        self.occ_mean = occ_stats.mean
#
#        # B-factor summary
#        b_stats = atom_group.atoms().extract_b().min_max_mean()
#        self.b_min = b_stats.min
#        self.b_max = b_stats.max
#        self.b_mean = b_stats.mean
#
#        # Coordinate summary
#        xyz_coords = atom_group.atoms().extract_xyz()
#        self.xyz_min = xyz_coords.min()
#        self.xyz_max = xyz_coords.max()
#        self.xyz_mean = xyz_coords.mean()



from Bamboo.Density.Edstats import edstats
from Giant.Structure.b_factors import b_factor_statistics, normalise_b_factors

class structureSummary(object):

    def __init__(self, pdb_input=None, pdb_hierarchy=None, mtz_file=None):
       """Summarise some characteristics of a structure"""

        assert [pdb_input, pdb_hierarchy].count(None)==1,'Provide pdb_input OR pdb_hierarchy'
        if pdb_input: pdb_hierarchy = pdb_input.construct_hierarchy()

        self.mtz_file       = mtz_file
        self.pdb_input      = pdb_input
        self.pdb_hierarchy  = pdb_hierarchy

        self.b_factors = b_factor_statistics(pdb_hierarchy=pdb_hierarchy)

        if mtz_file: self._score_against_density()

    def _score_against_density(self, mtz_file):
        self.edstats = edstats(mtz_file=mtz_file, pdb_file=pdb_file)

    def normalise_b_factors(self):
        return normalise_b_factors(pdb_hierarchy=self.pdb_hierarchy, b_factor_summary=self.b_factors):

    def to_html(self, fname):
        pass




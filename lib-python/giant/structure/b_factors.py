from scitbx.math import basic_statistics

class b_factor_summary(object):
    def __init__(self, all=None, protein=None, backbone=None, sidechain=None):
        self.all       = all
        self.protein   = protein
        self.backbone  = backbone
        self.sidechain = sidechain

    def to_z_score(self, b_vals, method='all'):
        """Convert the atoms B-factors to Z-scores"""
        assert method in ['all','backbone','sidechain']
        if method == 'all':       stats = self.all
        if method == 'protein':   stats = self.protein
        if method == 'backbone':  stats = self.backbone
        if method == 'sidechain': stats = self.sidechain
        return (b_vals - stats.mean)/stats.biased_standard_deviation

    def show(self):
        if self.all:
            print '================>'
            print 'All Atoms:'
            print self.all.show()
        if self.protein:
            print '================>'
            print 'Protein Atoms:'
            print self.protein.show()
        if self.backbone:
            print '================>'
            print 'Backbone Atoms:'
            print self.backbone.show()
        if self.sidechain:
            print '================>'
            print 'Sidechain Atoms:'
            print self.sidechain.show()

def b_factor_statistics(pdb_input=None, pdb_hierarchy=None):
    """Calculate the b-factor statistics of the model"""

    assert [pdb_input, pdb_hierarchy].count(None)==1,'Provide pdb_input OR pdb_hierarchy'
    if pdb_input: pdb_hierarchy = pdb_input.construct_hierarchy()

    cache = pdb_hierarchy.atom_selection_cache()
    sel_protein   = cache.selection('pepnames')
    sel_backbone  = cache.selection('pepnames and (name C or name CA or name N or name O)')
    sel_sidechain = cache.selection('pepnames and not (name C or name CA or name N or name O)')

    all_b = pdb_hierarchy.atoms().extract_b()
    protein_b   = pdb_hierarchy.select(sel_protein).atoms().extract_b()
    backbone_b  = pdb_hierarchy.select(sel_backbone).atoms().extract_b()
    sidechain_b = pdb_hierarchy.select(sel_sidechain).atoms().extract_b()

    summary = b_factor_summary( all       = basic_statistics(all_b),
                                protein   = basic_statistics(protein_b),
                                backbone  = basic_statistics(backbone_b),
                                sidechain = basic_statistics(sidechain_b) )
    return summary

def normalise_b_factors(pdb_input=None, pdb_hierarchy=None, b_factor_summary=None):
    """Calculate the b-factor statistics of the model"""

    assert [pdb_input, pdb_hierarchy].count(None)==1,'Provide pdb_input OR pdb_hierarchy'
    if pdb_input: pdb_hierarchy = pdb_input.construct_hierarchy()

    if not b_factor_summary: b_factor_summary = b_factor_statistics(pdb_hierarchy=pdb_hierarchy)
    new_b = b_factor_summary.to_z_score(b_vals=pdb_hierarchy.atoms().extract_b(), method='backbone')
    output_h = pdb_hierarchy.deep_copy()
    output_h.atoms().set_b(new_b)
    return output_h


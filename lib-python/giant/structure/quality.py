
class SelectionSummary(object):

    def __init__(self, ref_hierarchy, sub_hierarchy=None, selection=None):
        """Summarise a selection of a hierarchy relative to the whole hierarchy"""

        assert [selection, sub_hierarchy].count(None)==1,'Provide selection OR sub_hierarchy'

        if selection:
            cache = ref_hierarchy.atom_selection_cache()
            lig_sel = cache.selection(ligand_selection)
            sub_hierarchy = hierarchy.select(lig_sel)

        self.hierarchy_ref = ref_hierarchy
        self.hierarchy_sub = sub_hierarchy








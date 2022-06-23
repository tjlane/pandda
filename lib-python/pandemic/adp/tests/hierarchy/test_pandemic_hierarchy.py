from pytest import approx, raises, mark

from libtbx import group_args
import iotbx.pdb

from pandemic.resources.structure_test_snippets import pdb_test_structure_atoms

def test_CreateHierarchicalModelTask():
    """Integration test"""

    import numpy

    from pandemic.adp.hierarchy import CreateHierarchicalModelTask

    hierarchy = iotbx.pdb.hierarchy.input(pdb_string=pdb_test_structure_atoms[0]).hierarchy

    hm = CreateHierarchicalModelTask(
        auto_levels = ['chain', 'ss'],
        custom_levels = None,
        overall_selection = 'resseq 1900:1920',
        cbeta_in_backbone = True,
        assign_het_residues_to_nearest_ss_groups = False,
        assign_het_residues_to_nearest_custom_groups = False,
        )

    hm.run(hierarchy)

    assert hm.result.level_labels == ['chain', 'sec. struct.']
    #
    assert hm.result.level_group_array.shape == (2, 171)
    #
    assert (hm.result.level_group_array[0] == 0).all()
    #
    assert set(hm.result.level_group_array[1]) == {0,1,2}
    assert numpy.where(hm.result.level_group_array[1]==0)[0].tolist() == list(range(0,49))
    assert numpy.where(hm.result.level_group_array[1]==1)[0].tolist() == list(range(49,82))
    assert numpy.where(hm.result.level_group_array[1]==2)[0].tolist() == list(range(82,171))
    #
    assert len(hm.result.level_group_selection_strings) == 2
    assert hm.result.level_group_selection_strings[0] == ["chain 'A'"]
    assert hm.result.level_group_selection_strings[1] == ["chain 'A' and resid 1900  through 1905", "chain 'A' and resid 1906  through 1909", "chain 'A' and resid 1910  through 1920"]
    #
    assert hm.result.level_group_tree.links == {0: {0: {1: [0, 1, 2]}}, 1: {0: {}, 1: {}, 2: {}}}
    #
    assert hm.result.overall_atom_mask.shape == (428,)
    assert numpy.where(hm.result.overall_atom_mask)[0].tolist() == list(range(0,171))

def test_translate_phenix_selections_to_pymol_selections_simple():

    from pandemic.adp.hierarchy.summary import translate_phenix_selections_to_pymol_selections_simple

    for phenix_str, pymol_str in [
            (
                ['chain B'],
                ['chain B'],
                ),
            (
                ['chain C and resid 10 through 21'],
                ['chain C and resi 10:21'],
                ),
            (
                ['chain B and resid 121B through 130'],
                ['chain B and resi 121B:130'],
                ),
            (
                ['chain A and resseq 15:20'],
                [None],
                ),
            (
                ['chain A', 'chains B'],
                ['chain A', None],
                ),
            ]:

        t = translate_phenix_selections_to_pymol_selections_simple(phenix_str)
        assert t == pymol_str


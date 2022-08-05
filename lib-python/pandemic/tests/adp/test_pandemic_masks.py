from pytest import approx, raises
import numpy
from scitbx.array_family import flex

from pandemic.resources.structure_test_snippets import pdb_test_structure_short as pdb_str

def test_uij_isotropic_mask():

    from iotbx.pdb import hierarchy
    pdb_h = hierarchy.input(pdb_string=pdb_str).hierarchy

    from pandemic.adp.uijs.masks import UijIsotropicMask

    sel = pdb_h.atom_selection_cache().selection('resseq 12')
    mask = UijIsotropicMask(sel)
    uijs = mask(pdb_h.atoms().extract_uij())
    uijs_man = [((u[0]+u[1]+u[2])/3.,)*3+(0.,)*3 for u in pdb_h.select(sel).atoms().extract_uij()]
    assert list(uijs.select(sel)) == uijs_man
    assert list(uijs.select(sel==False)) == list(pdb_h.select(sel==False).atoms().extract_uij())

    mask_all = UijIsotropicMask(sel)
    mask_cut = mask_all[pdb_h.atom_selection_cache().selection('resseq 12:13')]
    assert mask_cut.selection.all_eq(flex.bool([True]*9 + [False]*8))

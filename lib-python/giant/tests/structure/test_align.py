import giant.logs as lg
logger = lg.getLogger(__name__)

import os

from pytest import (
    approx,
    mark,
    fixture,
    )

from iotbx.pdb.atom_selection import AtomSelectionError

@fixture
def example_pdbs(tmpdir):
    from giant.resources.test_data import get_test_data
    pdb_mov, pdb_ref = sorted(
        get_test_data(
            str(tmpdir),
            n = 2,
        )
    )
    return pdb_mov, pdb_ref

@mark.xfail(raises=AtomSelectionError, reason="library clash?")
def test_align_structures_rigid(example_pdbs):

    pdb_mov, pdb_ref = example_pdbs

    import iotbx.pdb.hierarchy
    h_mov = iotbx.pdb.hierarchy.input(pdb_mov).hierarchy
    h_ref = iotbx.pdb.hierarchy.input(pdb_ref).hierarchy

    from giant.structure.align import align_structures_rigid

    alignment = align_structures_rigid(
        mov_hierarchy = h_mov,
        ref_hierarchy = h_ref,
    )

    from scitbx.array_family import flex
    test_coords = flex.vec3_double([(4., 6., -2.)])

    assert alignment.nat2ref(test_coords)[0] == approx((4.037401, 5.862414, -2.025546))
    assert alignment.ref2nat(test_coords)[0] == approx((3.962991, 6.137686, -1.974420))

@mark.xfail(raises=AtomSelectionError, reason="library clash?")
def test_align_structures_flexible(example_pdbs):

    pdb_mov, pdb_ref = example_pdbs

    import iotbx.pdb.hierarchy
    h_mov = iotbx.pdb.hierarchy.input(pdb_mov).hierarchy
    h_ref = iotbx.pdb.hierarchy.input(pdb_ref).hierarchy

    from giant.structure.align import align_structures_flexible

    alignment = align_structures_flexible(
        mov_hierarchy = h_mov,
        ref_hierarchy = h_ref,
        altlocs = ['','A'],
        cutoff_radius = 15,
        sequence_identity_threshold = 0.95,
        one_to_one_mapping = True,
        require_hierarchies_identical = True,
    )

    from scitbx.array_family import flex
    test_coords = flex.vec3_double([(4., 6., -2.)])

    assert alignment.nat2ref(test_coords)[0] == approx((3.999301, 5.855366, -1.975495))
    assert alignment.ref2nat(test_coords)[0] == approx((4.001194, 6.144656, -2.024354))

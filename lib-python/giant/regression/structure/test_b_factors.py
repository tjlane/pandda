import giant.logs as lg
logger = lg.getLogger(__name__)

import os, numpy

from pytest import approx, mark, fixture

from iotbx.pdb.atom_selection import AtomSelectionError

@fixture
def example_pdb(tmpdir):
    from giant.resources.test_data import get_test_data
    pdb_file = get_test_data(str(tmpdir), n=1)
    return pdb_file

@mark.xfail(raises=AtomSelectionError, reason="library clash?")
def test_b_factor_statistics(example_pdb):

    pdb_file = example_pdb

    from giant.structure.b_factors import BfactorStatistics
    bfs = BfactorStatistics.from_pdb(
        pdb_file = pdb_file,
    )

    s = bfs.statistics['all_atoms']
    assert s.mean == approx(35.574075)
    assert s.biased_standard_deviation == approx(13.600492)
    s = bfs.statistics['protein']
    assert s.mean == approx(32.654358)
    assert s.biased_standard_deviation == approx(11.714119)
    s = bfs.statistics['backbone']
    assert s.mean == approx(30.351616)
    assert s.biased_standard_deviation == approx(9.774830)
    s = bfs.statistics['sidechain']
    assert s.mean == approx(34.884990)
    assert s.biased_standard_deviation == approx(12.941274)

    b_factors = numpy.array([10, 20, 30])

    z = bfs.to_z_scores(b_factors, normalise_by='all_atoms')
    assert z == approx([-1.88037869, -1.14511116, -0.40984362])

    z = bfs.to_z_scores(b_factors, normalise_by='protein')
    assert z == approx([-1.93393612, -1.08026546, -0.22659479])

    z = bfs.to_z_scores(b_factors, normalise_by='backbone')
    assert z == approx([-2.08204293, -1.05900727, -0.03597161])

    z = bfs.to_z_scores(b_factors, normalise_by='sidechain')
    assert z == approx([-1.92291647, -1.15019504, -0.37747361])

    str(bfs)

@mark.xfail(raises=AtomSelectionError, reason="library clash?")
def test_b_factor_z_scores(example_pdb):

    pdb_file = example_pdb

    from giant.structure.b_factors import normalise_b_factors_to_z_scores
    zh = normalise_b_factors_to_z_scores(
        pdb_file = pdb_file,
    )

    assert list(zh.atoms().extract_b()[5:10]) == approx(
        [
            -0.16794321004786722,
            0.2985610521354835,
            0.0919078482735604,
            0.6545774627490928,
            0.39268033310229933,
        ]
    )

def test_b_factor_ratios(example_pdb):

    pdb_file = example_pdb

    import iotbx.pdb.hierarchy
    h = iotbx.pdb.hierarchy.input(pdb_file).hierarchy

    rg = list(h.residue_groups())[20]

    from giant.structure.b_factors import calculate_residue_group_bfactor_ratio
    out = calculate_residue_group_bfactor_ratio(
        residue_group = rg,
        hierarchy = h,
        #distance_cutoff,
    )

    assert out.selection_name          == "A-1876"
    assert out.surroundings_names      == [
        'A-1872-LEU-()',
        'A-1873-CYS-()',
        'A-1874-SER-()',
        'A-1875-MET-()',
        'A-1877-LEU-()',
        'A-1878-THR-()',
        'A-1879-GLU-()',
        'A-1880-MET-()',
        'A-1961-PHE-()',
        'A-1965-TRP-()',
        'A-1968-THR-()',
        'A-1969-PHE-()',
    ]
    assert out.selection_av_bfactor    == approx(26.05375)
    assert out.surroundings_av_bfactor == approx(32.66340)
    assert out.b_factor_ratio          == approx(0.797644)


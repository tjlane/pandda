import giant.logs as lg
logger = lg.getLogger(__name__)

import os

from pytest import mark, fixture

from giant.paths import is_available

only_if_mtzdmp_is_available = mark.skipif(
    not is_available('mtzdmp'),
    reason = "mtzdmp is not available",
)

@fixture
def example_mtz_file(tmpdir):
    from giant.resources.test_data import get_test_data
    pdb = get_test_data(str(tmpdir), n=1)
    mtz = pdb.replace('.pdb', '.mtz')
    assert os.path.exists(mtz)
    return mtz

def test_mtz_header():

    from giant.io.mtz import MtzHeaderData
    header = MtzHeaderData(
        resolution_low = 99.9,
        resolution_high = 1.0,
        spacegroup = 'P1',
        spacegroup_no = 1,
        cell = (2., 3., 4., 90., 56., 87.,),
    )

    assert header.resolution_low == 99.9
    assert header.resolution_high == 1.0
    assert header.spacegroup == 'P1'
    assert header.spacegroup_no == 1
    assert header.cell == (2., 3., 4., 90., 56., 87.,)

    logger(str(header))
    logger(header.summary())

def test_mtz_column_labels():

    from giant.io.mtz import MtzColumnLabels
    labels = MtzColumnLabels(
        f = 'F',
        sigf = 'SIGF',
        i = 'I',
        sigi = 'SIGI',
        free = 'FREER',
        fom = 'FOM',
        f_calc = 'FC',
        phi_calc = 'PHIC',
        f_comp = 'FWT',
        phi_comp = 'PHWT',
        f_diff = 'DELFWT',
        phi_diff = 'DELPHWT',
    )

    assert labels.f == 'F'
    assert labels.sigf == 'SIGF'
    assert labels.i == 'I'
    assert labels.sigi == 'SIGI'
    assert labels.free == 'FREER'
    assert labels.fom == 'FOM'
    assert labels.f_calc == 'FC'
    assert labels.phi_calc == 'PHIC'
    assert labels.f_comp == 'FWT'
    assert labels.phi_comp == 'PHWT'
    assert labels.f_diff == 'DELFWT'
    assert labels.phi_diff == 'DELPHWT'

    logger(str(labels))
    logger(labels.summary())

@only_if_mtzdmp_is_available
def test_mtz_summary(example_mtz_file):

    from giant.io.mtz import MtzSummary
    s = MtzSummary(example_mtz_file)

    assert s.header.summary().strip() == """
        Resolution Range: 62.616 - 1.786 A
        Spacegroup: C 2 2 21 (No. 20)
        Cell: (82.128, 96.768, 57.955, 90.0, 90.0, 90.0)
        """.strip()

    assert s.labels.summary().strip() == """
        MTZ Column labels:
            F: F (SIGF)
            I: None (None)
            FREE: FreeR_flag
            FOM: FOM
            CALC: FC / PHIC
            COMP: 2FOFCWT / PH2FOFCWT
            DIFF: FOFCWT / PHFOFCWT
        """.strip()

@only_if_mtzdmp_is_available
def test_mtz_summary_dict(example_mtz_file):

    from giant.io.mtz import get_mtz_summary_dict
    summary_dict = get_mtz_summary_dict(example_mtz_file)

    expected_dict  = {
        'cell': [82.128, 96.768, 57.955, 90.0, 90.0, 90.0],
        'spacegroup': 'C 2 2 21',
        'spacegroupno': 20,
        'reslow': 62.616,
        'reshigh': 1.786,
        'numreflections': 22289,
        'numcols': 17,
        'coltypes': ['H', 'H', 'H', 'I', 'F', 'Q', 'F', 'P', 'F', 'P', 'F', 'P', 'F', 'P', 'W', 'F', 'P'],
        'colheadings': ['H', 'K', 'L', 'FreeR_flag', 'F', 'SIGF', 'FC', 'PHIC', 'FC_ALL', 'PHIC_ALL', '2FOFCWT', 'PH2FOFCWT', 'FOFCWT', 'PHFOFCWT', 'FOM', 'FC_ALL_LS', 'PHIC_ALL_LS'],
        'coldatasets': ['0', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1'],
    }

    assert summary_dict == expected_dict

def test_assign_labels():

    from giant.io.mtz import assign_column_labels
    labels_dict = assign_column_labels(
        col_labs = ['H', 'K', 'L', 'FreeR_flag', 'F', 'SIGF', 'FC', 'PHIC', 'FC_ALL', 'PHIC_ALL', '2FOFCWT', 'PH2FOFCWT', 'FOFCWT', 'PHFOFCWT', 'FOM', 'FC_ALL_LS', 'PHIC_ALL_LS'],
        col_types = ['H', 'H', 'H', 'I', 'F', 'Q', 'F', 'P', 'F', 'P', 'F', 'P', 'F', 'P', 'W', 'F', 'P'],
    )

    expected_dict = {
        'f_labels': ['F'],
        'sigf_labels': ['SIGF'],
        'i_labels': [],
        'sigi_labels': [],
        'r_free_labels': ['FreeR_flag'],
        'fom_labels' : ['FOM'],
        'f_calc_labels': ['FC', 'FC_ALL', 'FC_ALL_LS'],
        'phi_calc_labels': ['PHIC', 'PHIC_ALL', 'PHIC_ALL_LS'],
        'f_comp_labels': ['2FOFCWT'],
        'phi_comp_labels': ['PH2FOFCWT'],
        'f_diff_labels': ['FOFCWT'],
        'phi_diff_labels': ['PHFOFCWT'],
        'unassigned_labels': ['H', 'K', 'L'],
        'assigned_labels': ['FreeR_flag', 'F', 'SIGF', 'FC', 'PHIC', 'FC_ALL', 'PHIC_ALL', '2FOFCWT', 'PH2FOFCWT', 'FOFCWT', 'PHFOFCWT', 'FOM', 'FC_ALL_LS', 'PHIC_ALL_LS'],
    }

    assert labels_dict == expected_dict

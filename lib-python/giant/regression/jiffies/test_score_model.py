import giant.logs as lg
logger = lg.getLogger(__name__)

import os

from pytest import mark, fixture

from giant.regression.jiffies import run_jiffy, check_files_exist

from giant.paths import is_available

only_if_edstats_is_available = mark.xfail(
    not is_available('edstats.pl'),
    reason = "edstats.pl is not available",
)

from giant.jiffies import score_model

@fixture
def example_data(tmpdir):

    from giant.resources.test_data import get_test_data
    pdbs = get_test_data(
        str(tmpdir),
        n = 3,
    )
    return pdbs

@only_if_edstats_is_available
def test_score_model_single(example_data, tmpdir):

    pdb1 = example_data[0]

    args = [
        pdb1,
        "resname=EDO",
        "mode=separate+combined",
        "out_dir={}".format(str(tmpdir)),
    ]

    run_jiffy(
        args = args,
        module = score_model,
    )

    check_files_exist(
        [
            os.path.join(str(tmpdir), 'residue_plots', 'BAZ2BA-x430.dimple-B-1.png'),
            os.path.join(str(tmpdir), 'residue_scores.csv'),
            os.path.join(str(tmpdir), 'residue_scores.html'),
        ]
    )

def test_score_model_multiple(example_data, tmpdir):

    pdb1, pdb2 = example_data[0:2]

    args = [
        pdb1,
        pdb2,
        "resname=EDO",
        "mode=separate+combined",
        "out_dir={}".format(str(tmpdir)),
    ]

    run_jiffy(
        args = args,
        module = score_model,
    )

    check_files_exist(
        [
            os.path.join(str(tmpdir), 'residue_plots', 'BAZ2BA-x430.dimple-B-1.png'),
            os.path.join(str(tmpdir), 'residue_plots', 'BAZ2BA-x431.dimple-B-1.png'),
            os.path.join(str(tmpdir), 'residue_plots', 'compare-residue-B-1.png'),
            os.path.join(str(tmpdir), 'residue_scores.csv'),
            os.path.join(str(tmpdir), 'residue_scores.html'),
        ]
    )

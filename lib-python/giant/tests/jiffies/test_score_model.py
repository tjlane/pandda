import giant.logs as lg
logger = lg.getLogger(__name__)

import os
import pathlib as pl

from pytest import mark, fixture

from giant.tests.jiffies import (
    run_module,
    check_files,
    )

from giant.paths import is_available
only_if_edstats_is_available = mark.xfail(
    not is_available('edstats.pl'),
    reason = "edstats.pl is not available",
)

def run_score_model(args):

    import giant.jiffies.score_model

    run_module(
        module = giant.jiffies.score_model,
        args = args,
        )

@fixture
def score_model_data(tmp_path):

    from giant.resources.test_data import get_test_data

    pdbs = get_test_data(tmp_path, 3)

    return pdbs

@only_if_edstats_is_available
def test_score_model_single(score_model_data, tmp_path):

    pdb1 = score_model_data[0]

    out_dir = (tmp_path / 'score_model')

    args = [
        pdb1,
        "include_resname=EDO",
        "out_dir={}".format(str(out_dir)),
    ]

    run_score_model(args=args)

    check_files(
        files_list = [
            p.name for p in out_dir.glob('*')
            ],
        check_list = [
            "residue_scores_BAZ2BA-x430_B-   1.png",
            "residue_scores.csv",
            "score_distributions.png",
            "score_model.log",
            ],
        )

def test_score_model_multiple(score_model_data, tmp_path):

    pdb1, pdb2 = score_model_data[0:2]

    out_dir = (tmp_path / 'score_model')

    args = [
        pdb1,
        pdb2,
        "include_resname=EDO",
        "out_dir={}".format(str(out_dir)),
    ]

    run_score_model(args=args)

    check_files(
        files_list = [
            p.name for p in out_dir.glob('*')
            ],
        check_list = [
            "residue_scores_BAZ2BA-x430_B-   1.png",
            "residue_scores_BAZ2BA-x431_B-   1.png",
            "residue_scores.csv",
            "score_distributions.png",
            "score_model.log",
            ],
        )

def test_score_model_multiple_from_directories(score_model_data, tmp_path):

    pdb_dirs = [
        str(
            pl.Path(p).parent
            )
        for p in score_model_data
        ]

    out_dir = (tmp_path / 'score_model')

    args = pdb_dirs + [
        "include_resname=EDO",
        "out_dir={}".format(str(out_dir)),
        "pdb_style=*.dimple.pdb",
    ]

    run_score_model(args=args)

    check_files(
        files_list = [
            p.name for p in out_dir.glob('*')
            ],
        check_list = [
            "residue_scores_BAZ2BA-x430_B-   1.png",
            "residue_scores_BAZ2BA-x431_B-   1.png",
            "residue_scores_BAZ2BA-x432_B-   1.png",
            "residue_scores.csv",
            "score_distributions.png",
            "score_model.log",
            ],
        )

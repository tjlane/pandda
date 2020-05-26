import giant.logs as lg
logger = lg.getLogger(__name__)

import os

from giant.paths import is_available
from giant.regression.jiffies import run_jiffy, check_files_exist

from giant.jiffies import cluster_datasets

from pytest import fixture
@fixture
def example_data(tmpdir):

    from giant.resources.test_data import get_test_data
    pdbs = get_test_data(
        str(tmpdir),
        n = 5,
    )
    return pdbs

def test_cluster_datasets_pdbs(example_data, tmpdir):

    pdb_files = example_data

    out_dir = os.path.join(str(tmpdir), 'clustered_datasets')
    log_file = os.path.join(out_dir, 'clustering.log')

    args = pdb_files + [
        'out_dir={}'.format(out_dir),
        'cutoff=0.004',
    ]

    run_jiffy(
        args = args,
        module = cluster_datasets,
    )

    check_files_exist(
        [
            log_file,
            os.path.join(out_dir, 'dendrograms/sg-C_2_2_21-all.png'),
            os.path.join(out_dir, 'dendrograms/sg-C_2_2_21-cluster-1.png'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-1'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-1/BAZ2BA-x430.dimple/BAZ2BA-x430.dimple.pdb'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-1/BAZ2BA-x430.dimple/BAZ2BA-x430.dimple.mtz'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-1/BAZ2BA-x431.dimple/BAZ2BA-x431.dimple.mtz'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-1/BAZ2BA-x431.dimple/BAZ2BA-x431.dimple.pdb'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-1/BAZ2BA-x432.dimple/BAZ2BA-x432.dimple.pdb'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-1/BAZ2BA-x432.dimple/BAZ2BA-x432.dimple.mtz'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-1/BAZ2BA-x433.dimple/BAZ2BA-x433.dimple.mtz'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-1/BAZ2BA-x433.dimple/BAZ2BA-x433.dimple.pdb'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-2/BAZ2BA-x434.dimple/BAZ2BA-x434.dimple.mtz'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-2/BAZ2BA-x434.dimple/BAZ2BA-x434.dimple.pdb'),
        ]
    )

def test_cluster_datasets_mtzs(example_data, tmpdir):

    pdb_files = example_data
    mtz_files = [p.replace('.pdb', '.mtz') for p in pdb_files]

    out_dir = os.path.join(str(tmpdir), 'clustered_datasets')
    log_file = os.path.join(out_dir, 'clustering.log')

    args = mtz_files + [
        'out_dir={}'.format(out_dir),
        'cutoff=0.004',
    ]

    run_jiffy(
        args = args,
        module = cluster_datasets,
    )

    check_files_exist(
        [
            log_file,
            os.path.join(out_dir, 'dendrograms/sg-C_2_2_21-all.png'),
            os.path.join(out_dir, 'dendrograms/sg-C_2_2_21-cluster-1.png'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-1'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-1/BAZ2BA-x430.dimple/BAZ2BA-x430.dimple.pdb'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-1/BAZ2BA-x430.dimple/BAZ2BA-x430.dimple.mtz'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-1/BAZ2BA-x431.dimple/BAZ2BA-x431.dimple.mtz'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-1/BAZ2BA-x431.dimple/BAZ2BA-x431.dimple.pdb'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-1/BAZ2BA-x432.dimple/BAZ2BA-x432.dimple.pdb'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-1/BAZ2BA-x432.dimple/BAZ2BA-x432.dimple.mtz'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-1/BAZ2BA-x433.dimple/BAZ2BA-x433.dimple.mtz'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-1/BAZ2BA-x433.dimple/BAZ2BA-x433.dimple.pdb'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-2/BAZ2BA-x434.dimple/BAZ2BA-x434.dimple.mtz'),
            os.path.join(out_dir, 'sg-C_2_2_21-cluster-2/BAZ2BA-x434.dimple/BAZ2BA-x434.dimple.pdb'),
        ]
    )

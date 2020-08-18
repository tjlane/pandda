import giant.logs as lg
logger = lg.getLogger(__name__)

import os, shutil

from giant.paths import is_available, splice_ext
from pytest import mark

from giant.regression.jiffies import run_jiffy, check_files_exist

from giant.regression.jiffies.multi_conf_test_data import \
        LIGCODE, \
        MERGED_PDB_STR, \
        GROUND_PDB_STR, \
        BOUND_PDB_STR

def test_split_conformations(tmpdir):

    import giant.jiffies.split_conformations as split_c

    in_dir = tmpdir

    # Write input string to PDB file
    input_pdb = os.path.join(str(in_dir), 'input.pdb')
    with open(input_pdb, 'w') as fh:
        fh.write(MERGED_PDB_STR)

    dp = split_c.master_phil.extract()

    output_log = os.path.join(str(in_dir), 'split_conf.log')
    output_bound = os.path.join(
        str(in_dir),
        splice_ext(
            splice_ext(
                input_pdb,
                dp.output.suffix_prefix,
            ),
            dp.options.by_residue_name.selected_name,
        ),
    )
    output_ground = os.path.join(
        str(in_dir),
        splice_ext(
            splice_ext(
                input_pdb,
                dp.output.suffix_prefix,
            ),
            dp.options.by_residue_name.unselected_name,
        ),
    )

    split_args = [
        input_pdb,
        'output.log={}'.format(output_log),
        'by_residue_name.resname={}'.format(LIGCODE),
    ]

    expected_output_files = [
        output_log,
        output_ground,
        output_bound,
    ]

    run_jiffy(
        args = split_args,
        module = split_c,
    )

    check_files_exist(expected_output_files)

    with open(output_ground) as fh:
        assert fh.read().strip() == GROUND_PDB_STR.strip()

    with open(output_bound) as fh:
        assert fh.read().strip() == BOUND_PDB_STR.strip()

def test_merge_conformations(tmpdir):

    import giant.jiffies.merge_conformations as merge_c

    in_dir = tmpdir

    # Write input string to PDB file
    input_ground = os.path.join(str(in_dir), 'input-ground.pdb')
    with open(input_ground, 'w') as fh:
        fh.write(GROUND_PDB_STR)

    # Write input string to PDB file
    input_bound = os.path.join(str(in_dir), 'input-bound.pdb')
    with open(input_bound, 'w') as fh:
        fh.write(BOUND_PDB_STR)

    output_pdb = os.path.join(str(in_dir), 'merged.pdb')
    output_log = os.path.join(str(in_dir), 'merged.log')
    output_restraints_phenix = os.path.join(str(in_dir), 'phenix.params')
    output_restraints_refmac = os.path.join(str(in_dir), 'refmac.params')

    merge_args = [
        "major={}".format(input_ground),
        "minor={}".format(input_bound),
        "output.pdb={}".format(output_pdb),
        "output.log={}".format(output_log),
        "make_restraints=False",
        #"restraints.output.phenix={}".format(output_restraints_phenix),
        #"restraints.output.refmac={}".format(output_restraints_refmac),
    ]

    expected_output_files = [
        output_pdb,
        output_log,
        #output_restraints_phenix,
        #output_restraints_refmac,
    ]

    run_jiffy(
        args = merge_args,
        module = merge_c,
    )

    check_files_exist(expected_output_files)

    with open(output_pdb) as fh:
        assert fh.read().strip() == MERGED_PDB_STR.strip()


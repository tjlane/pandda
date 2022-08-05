import giant.logs as lg
logger = lg.getLogger(__name__)

from pytest import mark

import os, copy

from giant.tests.jiffies import (
    run_module,
    check_files,
    )

def run_make_restraints(args):

    import giant.jiffies.make_restraints

    return run_module(
        module = giant.jiffies.make_restraints,
        args = args,
        )

@mark.slow
def test_make_restraints_two_states(resources, tmp_path):

    input_pdb = resources.merged_apo_glc(tmp_path/"input.pdb")

    out_root = (tmp_path / "restraints")

    run_make_restraints(
        args=[
            "input.pdb={}".format(input_pdb),
            "output_root={}".format(out_root),
            ],
        )

    out_phenix = out_root.with_suffix('.phenix.params')
    out_refmac = out_root.with_suffix('.refmac.params')

    check_files(
        files_list = [
            p.name for p in tmp_path.glob('*')
            ],
        check_list = [
            input_pdb.name,
            out_phenix.name,
            out_refmac.name,
            out_root.with_suffix('.log').name,
            ],
        )

    resources.assert_same(
        resources.restraints_apo_glc_refmac(),
        out_refmac,
        )

    resources.assert_same(
        resources.restraints_apo_glc_phenix(),
        out_phenix,
        )

@mark.slow
def test_make_restraints_three_states(resources, tmp_path):

    input_pdb = resources.merged_apo_glc_glo(tmp_path/"input.pdb")

    out_root = (tmp_path / "restraints")

    run_make_restraints(
        args=[
            "input.pdb={}".format(input_pdb),
            "output_root={}".format(out_root),
            ],
        )

    out_phenix = out_root.with_suffix('.phenix.params')
    out_refmac = out_root.with_suffix('.refmac.params')

    check_files(
        files_list = [
            p.name for p in tmp_path.glob('*')
            ],
        check_list = [
            input_pdb.name,
            out_phenix.name,
            out_refmac.name,
            out_root.with_suffix('.log').name,
            ],
        )

    resources.assert_same(
        resources.restraints_apo_glc_glo_refmac(),
        out_refmac,
        )

    resources.assert_same(
        resources.restraints_apo_glc_glo_phenix(),
        out_phenix,
        )

@mark.slow
def test_make_restraints_four_states(resources, tmp_path):

    input_pdb = resources.merged_apo_glc_glo_z9n(tmp_path/"input.pdb")

    out_root = (tmp_path / "restraints")

    run_make_restraints(
        args=[
            "input.pdb={}".format(input_pdb),
            "output_root={}".format(out_root),
            ],
        )

    out_phenix = out_root.with_suffix('.phenix.params')
    out_refmac = out_root.with_suffix('.refmac.params')

    check_files(
        files_list = [
            p.name for p in tmp_path.glob('*')
            ],
        check_list = [
            input_pdb.name,
            out_phenix.name,
            out_refmac.name,
            out_root.with_suffix('.log').name,
            ],
        )

    resources.assert_same(
        resources.restraints_apo_glc_glo_z9n_refmac(),
        out_refmac,
        )

    resources.assert_same(
        resources.restraints_apo_glc_glo_z9n_phenix(),
        out_phenix,
        )



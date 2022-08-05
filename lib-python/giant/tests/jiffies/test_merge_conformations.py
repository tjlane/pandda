import giant.logs as lg
logger = lg.getLogger(__name__)

import os, copy
import pytest

from pytest import mark

from giant.tests.jiffies import (
    run_module,
    check_files,
    )

def run_merge_conformations(args):

    import giant.jiffies.merge_conformations

    return run_module(
        module = giant.jiffies.merge_conformations,
        args = args,
        )

@mark.slow
def test_merge_conformations_two_states(resources, tmp_path):

    input_apo = resources.input_apo(tmp_path/"state1.pdb")
    input_glc = resources.input_glc(tmp_path/"state2.pdb")

    output = (tmp_path / "merged.pdb")

    run_merge_conformations(
        args=[
            "input.pdb={}".format(input_apo),
            "input.pdb={}".format(input_glc),
            "output.pdb={}".format(output),
            "restraints.output.output_root={}".format(
                output.with_name(output.stem),
                ),
            ],
        )

    check_files(
        files_list = [
            p.name for p in tmp_path.glob('*')
            ],
        check_list = [
            input_apo.name, 
            input_glc.name, 
            output.name, 
            output.with_suffix('.log').name,
            output.stem+'.phenix.params',
            output.stem+'.refmac.params',
            ],
        )

    resources.assert_same(
        resources.merged_apo_glc(), 
        output,
        )

    resources.assert_same(
        resources.restraints_apo_glc_phenix(),
        output.with_suffix('.phenix.params'),
        )

    resources.assert_same(
        resources.restraints_apo_glc_refmac(),
        output.with_suffix('.refmac.params'),
        )

@mark.slow
def test_merge_conformations_two_states_no_restraints(resources, tmp_path):

    input_apo = resources.input_apo(tmp_path/"state1.pdb")
    input_glc = resources.input_glc(tmp_path/"state2.pdb")

    output = (tmp_path / "merged.pdb")

    run_merge_conformations(
        args=[
            "input.pdb={}".format(input_apo),
            "input.pdb={}".format(input_glc),
            "output.pdb={}".format(output),
            "make_restraints=False",
            ],
        )

    check_files(
        files_list = [
            p.name for p in tmp_path.glob('*')
            ],
        check_list = [
            input_apo.name, 
            input_glc.name, 
            output.name, 
            output.with_suffix('.log').name,
            ],
        )

    resources.assert_same(
        resources.merged_apo_glc(), 
        output,
        )

@mark.slow
def test_merge_conformations_three_states_no_restraints(resources, tmp_path):

    input_apo = resources.input_apo(tmp_path/"state1.pdb")
    input_glc = resources.input_glc(tmp_path/"state2.pdb")
    input_glo = resources.input_glo(tmp_path/"state3.pdb")

    output = (tmp_path / "merged.pdb")

    run_merge_conformations(
        args=[
            "input.pdb={}".format(input_apo),
            "input.pdb={}".format(input_glc),
            "input.pdb={}".format(input_glo),
            "output.pdb={}".format(output),
            "make_restraints=False",
            ],
        )

    check_files(
        files_list = [
            p.name for p in tmp_path.glob('*')
            ],
        check_list = [
            input_apo.name, 
            input_glc.name, 
            input_glo.name, 
            output.name, 
            output.with_suffix('.log').name,
            ],
        )

    resources.assert_same(
        resources.merged_apo_glc_glo(), 
        output,
        )

@mark.slow
def test_merge_conformations_four_states_no_restraints(resources, tmp_path):

    input_apo = resources.input_apo(tmp_path/"state1.pdb")
    input_glc = resources.input_glc(tmp_path/"state2.pdb")
    input_glo = resources.input_glo(tmp_path/"state3.pdb")
    input_z9n = resources.input_z9n(tmp_path/"state4.pdb")

    output = (tmp_path / "merged.pdb")

    run_merge_conformations(
        args=[
            "input.pdb={}".format(input_apo),
            "input.pdb={}".format(input_glc),
            "input.pdb={}".format(input_glo),
            "input.pdb={}".format(input_z9n),
            "output.pdb={}".format(output),
            "make_restraints=False",
            ],
        )

    check_files(
        files_list = [
            p.name for p in tmp_path.glob('*')
            ],
        check_list = [
            input_apo.name, 
            input_glc.name, 
            input_glo.name, 
            input_z9n.name, 
            output.name, 
            output.with_suffix('.log').name,
            ],
        )

    resources.assert_same(
        resources.merged_apo_glc_glo_z9n(), 
        output,
        )

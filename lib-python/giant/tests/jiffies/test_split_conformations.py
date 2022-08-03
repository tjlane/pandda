import os, copy
import pytest

from giant.tests.jiffies import (
    run_module,
    check_files,
    )

def run_split_conformations(args):

    import giant.jiffies.split_conformations

    run_module(
        module = giant.jiffies.split_conformations,
        args = args,
        )

def test_split_two_states_defaults(resources, tmp_path):

    input_pdb = resources.merged_apo_glc(tmp_path/"merged.pdb")

    run_split_conformations(
        args=[
            "input.pdb={}".format(input_pdb),
            ],
        )

    output_ground = tmp_path / (
        input_pdb.stem+'-split-ground-state-ABCDEF.pdb'
        )

    check_files(
        files_list = [
            p.name for p in tmp_path.glob('*')
            ],
        check_list = [
            input_pdb.name,
            input_pdb.stem+'-split_conformations.log',
            output_ground.name,
            ],
        )

    resources.assert_same(
        resources.split_apo_glc_to_apo_glc(),
        output_ground,
        )

def test_split_two_states_by_conformer(resources, tmp_path):

    input_pdb = resources.merged_apo_glc(tmp_path/"merged.pdb")

    run_split_conformations(
        args=[
            "input.pdb={}".format(input_pdb),
            "mode=by_conformer",
            ],
        )

    outputs = [
        tmp_path / (
            input_pdb.stem+'-split-{}.pdb'.format(l)
            )
        for l in 'ABCDEF'
        ]

    check_files(
        files_list = [
            p.name for p in tmp_path.glob('*')
            ],
        check_list = [
            input_pdb.name,
            input_pdb.stem+'-split_conformations.log',
            ] + [
            o.name for o in outputs
            ],
        )

    resources.assert_same(
        resources.split_apo_glc_to_A(),
        outputs[0],
        )

    resources.assert_same(
        resources.split_apo_glc_to_B(),
        outputs[1],
        )

    resources.assert_same(
        resources.split_apo_glc_to_C(),
        outputs[2],
        )

    resources.assert_same(
        resources.split_apo_glc_to_D(),
        outputs[3],
        )

    resources.assert_same(
        resources.split_apo_glc_to_E(),
        outputs[4],
        )

    resources.assert_same(
        resources.split_apo_glc_to_F(),
        outputs[5],
        )

def test_split_two_states_by_residue_name_ignore_common_false(resources, tmp_path):

    input_pdb = resources.merged_apo_glc(tmp_path/"merged.pdb")

    run_split_conformations(
        args=[
            "input.pdb={}".format(input_pdb),
            "mode=by_residue_name",
            "ignore_common_molecules=False",
            ],
        )

    output_ground = tmp_path / (
        input_pdb.stem+'-split-ground-state-ABC.pdb'
        )

    output_bound = tmp_path / (
        input_pdb.stem+'-split-bound-state-DEF.pdb'
        )

    check_files(
        files_list = [
            p.name for p in tmp_path.glob('*')
            ],
        check_list = [
            input_pdb.name,
            input_pdb.stem+'-split_conformations.log',
            output_ground.name,
            output_bound.name,
            ],
        )

    resources.assert_same(
        resources.split_apo_glc_to_apo(),
        output_ground,
        )

    resources.assert_same(
        resources.split_apo_glc_to_glc(),
        output_bound,
        )

def test_split_two_states_by_residue_name_include_resname(resources, tmp_path):

    input_pdb = resources.merged_apo_glc(tmp_path/"merged.pdb")

    run_split_conformations(
        args=[
            "input.pdb={}".format(input_pdb),
            "mode=by_residue_name",
            "include_resname=GLC",
            ],
        )

    output_ground = tmp_path / (
        input_pdb.stem+'-split-ground-state-ABC.pdb'
        )

    output_bound = tmp_path / (
        input_pdb.stem+'-split-bound-state-DEF.pdb'
        )

    check_files(
        files_list = [
            p.name for p in tmp_path.glob('*')
            ],
        check_list = [
            input_pdb.name,
            input_pdb.stem+'-split_conformations.log',
            output_ground.name,
            output_bound.name,
            ],
        )

    resources.assert_same(
        resources.split_apo_glc_to_apo(),
        output_ground,
        )

    resources.assert_same(
        resources.split_apo_glc_to_glc(),
        output_bound,
        )

def test_split_two_states_by_residue_name_ignore_resname(resources, tmp_path):

    input_pdb = resources.merged_apo_glc(tmp_path/"merged.pdb")

    run_split_conformations(
        args=[
            "input.pdb={}".format(input_pdb),
            "mode=by_residue_name",
            "ignore_common_molecules=False",
            "ignore_resname=GLC",
            ],
        )

    output_ground = tmp_path / (
        input_pdb.stem+'-split-ground-state-ABCDEF.pdb'
        )

    check_files(
        files_list = [
            p.name for p in tmp_path.glob('*')
            ],
        check_list = [
            input_pdb.name,
            input_pdb.stem+'-split_conformations.log',
            output_ground.name,
            ],
        )

    resources.assert_same(
        resources.split_apo_glc_to_apo_glc(),
        output_ground,
        )

def test_split_two_states_by_residue_name_atom_selection(resources, tmp_path):

    input_pdb = resources.merged_apo_glc(tmp_path/"merged.pdb")

    run_split_conformations(
        args=[
            "input.pdb={}".format(input_pdb),
            "mode=by_residue_name",
            "by_residue_name.atom_selection='resname GLC'",
            ],
        )

    output_ground = tmp_path / (
        input_pdb.stem+'-split-ground-state-ABC.pdb'
        )

    output_bound = tmp_path / (
        input_pdb.stem+'-split-bound-state-DEF.pdb'
        )

    check_files(
        files_list = [
            p.name for p in tmp_path.glob('*')
            ],
        check_list = [
            input_pdb.name,
            input_pdb.stem+'-split_conformations.log',
            output_ground.name,
            output_bound.name,
            ],
        )

    resources.assert_same(
        resources.split_apo_glc_to_apo(),
        output_ground,
        )

    resources.assert_same(
        resources.split_apo_glc_to_glc(),
        output_bound,
        )

def test_split_two_states_by_residue_name_output_names(resources, tmp_path):

    input_pdb = resources.merged_apo_glc(tmp_path/"merged.pdb")

    run_split_conformations(
        args=[
            "input.pdb={}".format(input_pdb),
            "mode=by_residue_name",
            "include_resname=GLC",
            'unselected_name=alpha',
            'selected_name=beta',
            ],
        )

    output_ground = tmp_path / (
        input_pdb.stem+'-split-alpha-ABC.pdb'
        )

    output_bound = tmp_path / (
        input_pdb.stem+'-split-beta-DEF.pdb'
        )

    check_files(
        files_list = [
            p.name for p in tmp_path.glob('*')
            ],
        check_list = [
            input_pdb.name,
            input_pdb.stem+'-split_conformations.log',
            output_ground.name,
            output_bound.name,
            ],
        )

    resources.assert_same(
        resources.split_apo_glc_to_apo(),
        output_ground,
        )

    resources.assert_same(
        resources.split_apo_glc_to_glc(),
        output_bound,
        )

def test_split_two_states_by_residue_name_no_prune(resources, tmp_path):

    input_pdb = resources.merged_apo_glc(tmp_path/"merged.pdb")

    run_split_conformations(
        args=[
            "input.pdb={}".format(input_pdb),
            "mode=by_residue_name",
            "include_resname=GLC",
            "prune_duplicates=False",
            ],
        )

    output_ground = tmp_path / (
        input_pdb.stem+'-split-ground-state-ABC.pdb'
        )

    output_bound = tmp_path / (
        input_pdb.stem+'-split-bound-state-DEF.pdb'
        )

    check_files(
        files_list = [
            p.name for p in tmp_path.glob('*')
            ],
        check_list = [
            input_pdb.name,
            input_pdb.stem+'-split_conformations.log',
            output_ground.name,
            output_bound.name,
            ],
        )

    resources.assert_same(
        resources.split_apo_glc_to_apo_no_prune(),
        output_ground,
        )

    resources.assert_same(
        resources.split_apo_glc_to_glc_no_prune(),
        output_bound,
        )

def test_split_two_states_by_residue_name_no_reset_altlocs(resources, tmp_path):

    input_pdb = resources.merged_apo_glc(tmp_path/"merged.pdb")

    run_split_conformations(
        args=[
            "input.pdb={}".format(input_pdb),
            "mode=by_residue_name",
            "include_resname=GLC",
            "reset_altlocs=False",
            ],
        )

    output_ground = tmp_path / (
        input_pdb.stem+'-split-ground-state-ABC.pdb'
        )

    output_bound = tmp_path / (
        input_pdb.stem+'-split-bound-state-DEF.pdb'
        )

    check_files(
        files_list = [
            p.name for p in tmp_path.glob('*')
            ],
        check_list = [
            input_pdb.name,
            input_pdb.stem+'-split_conformations.log',
            output_ground.name,
            output_bound.name,
            ],
        )

    resources.assert_same(
        resources.split_apo_glc_to_apo_no_reset_altlocs(),
        output_ground,
        )

    resources.assert_same(
        resources.split_apo_glc_to_glc_no_reset_altlocs(),
        output_bound,
        )

def test_split_three_states_defaults(resources, tmp_path):

    input_pdb = resources.merged_apo_glc_glo(tmp_path/"merged.pdb")

    run_split_conformations(
        args=[
            "input.pdb={}".format(input_pdb),
            ],
        )

    # only GLO will go into the bound-split

    output_ground = tmp_path / (
        input_pdb.stem+'-split-ground-state-ABCDEF.pdb'
        )

    output_bound = tmp_path / (
        input_pdb.stem+'-split-bound-state-GHI.pdb'
        )

    check_files(
        files_list = [
            p.name for p in tmp_path.glob('*')
            ],
        check_list = [
            input_pdb.name,
            input_pdb.stem+'-split_conformations.log',
            output_ground.name,
            output_bound.name,
            ],
        )

    resources.assert_same(
        resources.split_apo_glc_glo_to_apo_glc(),
        output_ground,
        )

    resources.assert_same(
        resources.split_apo_glc_glo_to_glo(),
        output_bound,
        )

def test_split_three_states_by_residue_name_combine(resources, tmp_path):

    input_pdb = resources.merged_apo_glc_glo(tmp_path/"merged.pdb")

    run_split_conformations(
        args=[
            "input.pdb={}".format(input_pdb),
            "mode=by_residue_name",
            "include_resname=GLC", # GLO should be included by default
            "combine_bound_states=True",
            ],
        )

    output_ground = tmp_path / (
        input_pdb.stem+'-split-ground-state-ABC.pdb'
        )

    output_bound = tmp_path / (
        input_pdb.stem+'-split-bound-state-DEFGHI.pdb'
        )

    check_files(
        files_list = [
            p.name for p in tmp_path.glob('*')
            ],
        check_list = [
            input_pdb.name,
            input_pdb.stem+'-split_conformations.log',
            output_ground.name,
            output_bound.name,
            ],
        )

    resources.assert_same(
        resources.split_apo_glc_glo_to_apo(),
        output_ground,
        )

    resources.assert_same(
        resources.split_apo_glc_glo_to_glc_glo(),
        output_bound,
        )

def test_split_three_states_by_conformer_group_1(resources, tmp_path):

    input_pdb = resources.merged_apo_glc_glo(tmp_path/"merged.pdb")

    run_split_conformations(
        args=[
            "input.pdb={}".format(input_pdb),
            "mode=by_conformer_group",
            "conformer=ABC",
            "conformer=DEFGHI",
            ],
        )

    output_1 = tmp_path / (
        input_pdb.stem+'-split-ABC.pdb'
        )

    output_2 = tmp_path / (
        input_pdb.stem+'-split-DEFGHI.pdb'
        )

    check_files(
        files_list = [
            p.name for p in tmp_path.glob('*')
            ],
        check_list = [
            input_pdb.name,
            input_pdb.stem+'-split_conformations.log',
            output_1.name,
            output_2.name,
            ],
        )

    resources.assert_same(
        resources.split_apo_glc_glo_to_apo(),
        output_1,
        )

    resources.assert_same(
        resources.split_apo_glc_glo_to_glc_glo(),
        output_2,
        )

def test_split_three_states_by_conformer_group_2(resources, tmp_path):

    input_pdb = resources.merged_apo_glc_glo(tmp_path/"merged.pdb")

    run_split_conformations(
        args=[
            "input.pdb={}".format(input_pdb),
            "mode=by_conformer_group",
            "conformer=ABC",
            "conformer=DEF",
            "conformer=GHI",
            "conformer=ABCDEF",
            "conformer=DEFGHI",
            ],
        )

    output_1 = tmp_path / (
        input_pdb.stem+'-split-ABC.pdb'
        )

    output_2 = tmp_path / (
        input_pdb.stem+'-split-DEF.pdb'
        )

    output_3 = tmp_path / (
        input_pdb.stem+'-split-GHI.pdb'
        )

    output_4 = tmp_path / (
        input_pdb.stem+'-split-ABCDEF.pdb'
        )

    output_5 = tmp_path / (
        input_pdb.stem+'-split-DEFGHI.pdb'
        )

    check_files(
        files_list = [
            p.name for p in tmp_path.glob('*')
            ],
        check_list = [
            input_pdb.name,
            input_pdb.stem+'-split_conformations.log',
            output_1.name,
            output_2.name,
            output_3.name,
            output_4.name,
            output_5.name,
            ],
        )

    resources.assert_same(
        resources.split_apo_glc_glo_to_apo(),
        output_1,
        )

    resources.assert_same(
        resources.split_apo_glc_glo_to_glc(),
        output_2,
        )

    resources.assert_same(
        resources.split_apo_glc_glo_to_glo(),
        output_3,
        )

    resources.assert_same(
        resources.split_apo_glc_glo_to_apo_glc(),
        output_4,
        )

    resources.assert_same(
        resources.split_apo_glc_glo_to_glc_glo(),
        output_5,
        )


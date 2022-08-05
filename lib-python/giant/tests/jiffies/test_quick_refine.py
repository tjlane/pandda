import os

from pytest import mark, xfail

from giant.paths import is_available

from giant.tests.jiffies import (
    run_module,
    )

fail_if_phenix_refine_unavailable = mark.xfail(
    not is_available('phenix.refine'),
    reason = "phenix.refine is not available",
)

fail_if_refmac_unavailable = mark.xfail(
    not is_available('refmac5'),
    reason = "refmac5 is not available",
)

def run_quick_refine(args):

    import giant.jiffies.quick_refine

    run_module(
        module = giant.jiffies.quick_refine,
        args = args,
    )

def run_quick_refine_simple(dir_path, other_args, run_no=1):

    from giant.resources.test_data import get_test_data

    pdb_file = get_test_data(dir_path, n=1)

    mtz_file = pdb_file.replace('.pdb', '.mtz')

    dir_prefix = str(dir_path / 'refine_test')

    link_prefix = str(dir_path / 'link_test')

    test_args = [
        str(pdb_file),
        str(mtz_file),
        'dir_prefix='+dir_prefix,
        'link_prefix='+link_prefix,
        ]

    for i in range(run_no):
        run_quick_refine(test_args+other_args)

@mark.slow
@fail_if_refmac_unavailable
def test_jiffies_quick_refine_refmac(tmp_path):

    params_str = (
        "NCYC 1"
        )

    params_file = (tmp_path / 'refmac.params')

    with open(str(params_file), 'w') as fh:
        fh.write(params_str)

    other_args = [
        'program=refmac',
        str(params_file),
        ]

    run_quick_refine_simple(
        tmp_path,
        other_args,
        run_no=2,
        ) # Run twice -- shouldn't crash

@mark.slow
@fail_if_phenix_refine_unavailable
def test_jiffies_quick_refine_phenix(tmp_path):

    params_str = """
    refinement {
      refine {
        strategy = individual_adp
      }
      main {
        number_of_macro_cycles = 1
      }
    }
    """

    params_file = (tmp_path / 'phenix.params')

    with open(str(params_file), 'w') as fh:
        fh.write(params_str)

    other_args = [
        'program=phenix',
        str(params_file),
        ]

    run_quick_refine_simple(
        tmp_path,
        other_args,
        )


import os

from pytest import mark, xfail

from giant.paths import is_available

fail_if_phenix_refine_unavailable = mark.xfail(
    not is_available('phenix.refine'),
    reason = "phenix.refine is not available",
)

fail_if_refmac_unavailable = mark.xfail(
    not is_available('refmac5'),
    reason = "refmac5 is not available",
)

def run_quick_refine(args):

    from giant.regression.jiffies import run_jiffy
    import giant.jiffies.quick_refine as test_mod

    run_jiffy(
        args = args,
        module = test_mod,
    )

def run_quick_refine_simple(directory, other_args, run_no=1):

    from giant.resources.test_data import get_test_data

    pdb_file = get_test_data(directory, n=1)
    mtz_file = pdb_file.replace('.pdb', '.mtz')
    dir_prefix = os.path.join(directory, 'refine_test')
    link_prefix = os.path.join(directory, 'link_test')

    test_args = [
        pdb_file,
        mtz_file,
        'dir_prefix='+dir_prefix,
        'link_prefix='+link_prefix,
        ]

    for i in range(run_no):
        run_quick_refine(test_args+other_args)

#@mark.slow
@fail_if_refmac_unavailable
def test_jiffies_quick_refine_refmac(tmpdir):

    params_str = """
    NCYC 1
    """
    params_file = os.path.join(str(tmpdir), 'refmac.params')
    with open(params_file, 'w') as fh:
        fh.write(params_str)

    other_args = ['program=refmac', params_file]
    run_quick_refine_simple(str(tmpdir), other_args, run_no=2) # Run twice -- shouldn't crash

@mark.slow
@fail_if_phenix_refine_unavailable
def test_jiffies_quick_refine_phenix(tmpdir):

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
    params_file = os.path.join(str(tmpdir), 'phenix.params')
    with open(params_file, 'w') as fh:
        fh.write(params_str)

    other_args = ['program=phenix', params_file]

    run_quick_refine_simple(str(tmpdir), other_args)


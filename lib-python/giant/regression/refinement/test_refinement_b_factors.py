import os

from pytest import mark

@mark.slow
def test_b_factor_factory(tmpdir):

    from giant.resources.test_data import get_test_data
    pdb_file = get_test_data(str(tmpdir), n=1)
    mtz_file = pdb_file.replace('.pdb', '.mtz')

    assert os.path.exists(pdb_file)
    assert os.path.exists(mtz_file)

    from giant.refinement.b_factors import BFactorRefinementFactory
    factory = BFactorRefinementFactory(
        pdb_file = pdb_file,
        mtz_file = mtz_file,
        out_dir = str(tmpdir),
        )

    # Make quick
    factory._n_cycles = 1

    factory.refine_b_factors(mode='isotropic')
    factory.refine_b_factors(mode='tls')
    factory.refine_b_factors(mode='anisotropic')

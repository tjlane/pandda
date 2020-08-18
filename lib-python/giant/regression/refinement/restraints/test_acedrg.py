from giant.paths import is_available

from pytest import mark

@mark.slow
def test_acedrg(tmpdir):

    # Do not error if can't find acedrg
    if not is_available('acedrg'):
        return None

    out_prefix = tmpdir.mkdir('acedrg1') / 'ligand'

    from giant.refinement.restraints.acedrg import generate_restraints

    out_pdb, out_cif = generate_restraints(
        smiles = 'OCCO',
        prefix = str(out_prefix),
    )


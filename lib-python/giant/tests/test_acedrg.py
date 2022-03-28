from giant.paths import is_available

from pytest import mark

@mark.slow
def test_acedrg(tmp_path):

    # Do not error if can't find acedrg
    if not is_available('acedrg'):
        return None

    out_prefix = (tmp_path / 'ligand')

    from giant.refinement.restraints.acedrg import generate_restraints

    out_pdb, out_cif = generate_restraints(
        smiles = 'OCCO',
        prefix = str(out_prefix),
    )


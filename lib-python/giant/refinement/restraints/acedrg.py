import giant.logs as lg
logger = lg.getLogger(__name__)

import os

from giant.dispatcher import Dispatcher

def generate_restraints(smiles, name='LIG', prefix='ligand', verbose=False):
    """Generate pdb and cif files from smiles string"""

    assert len(name) == 3

    # Common fixes
    smiles = smiles.replace('CL', 'Cl')
    smiles = smiles.replace('BR', 'Br')

    out_pdb = prefix+'.pdb'
    out_cif = prefix+'.cif'
    out_log = prefix+'-acedrg.log'

    if os.path.exists(out_pdb):
        raise IOError('Output PDB file already exists: {}'.format(out_pdb))
    if os.path.exists(out_cif):
        raise IOError('Output CIF file already exists: {}'.format(out_cif))

    # Run acedrg
    acedrg = Dispatcher('acedrg')

    acedrg.extend_args([
        '--smi={}'.format(smiles),
        '-r', name,
        '-o', prefix,
    ])

    if (verbose is True):
        logger(acedrg.as_string())

    acedrg.run()
    acedrg.write_output(out_log)

    if not (os.path.exists(out_pdb) and os.path.exists(out_cif)):
        logger(str(acedrg.result.stdout))
        logger(str(acedrg.result.stderr))
        raise Exception('acedrg failed during ligand generation')

    return out_pdb, out_cif


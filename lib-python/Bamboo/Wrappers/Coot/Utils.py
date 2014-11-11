#! /usr/local/python/python2.7.3-64bit/bin/python

# Import all coot scripts (all are strings starting with COOT_)
from Bamboo.Common.Command import CommandManager
from Bamboo.Wrappers.Coot.Scripts import *
from Bamboo.Utils.MTZ import mtzFile

def validate_coot_script(script):
    """Performs checks on a coot script"""

    contents = open(script, 'r').read().strip()

    assert contents.endswith('coot_real_exit(1)'), 'Coot script does not EXIT! {!s}'.format(script)

    return

def run_coot(script, graphical=False, noguano=True):
    """Runs coot with the provided script"""

    # Check that the script is valid (e.g. causes coot to exit at the end)
    validate_coot_script(script)

    # Stop coot droppings?
    coot_flags = ['--no-guano']*(noguano)
    # Run with/without graphics?
    if graphical:
        coot_flags.extend(['-script',script])
    else:
        coot_flags.extend(['--no-graphics','-s',script])

    # Initialise
    COOT = CommandManager('coot')
    # Load arguments
    COOT.SetArguments(*coot_flags)
    # Run!
    COOT.Run()

    return COOT

def coot_real_space_refine(pdbin, mtzin, pdbout, cifin=None, scriptname=None):
    """Uses coot's regularize-zone function to real-space-refine a ligand"""

    # Generate scriptname
    if not scriptname:
        scriptname = pdbout + '.coot_script.py'

    # Get mtz summary
    mtz_obj = mtzFile(mtzin)

    coot_script = []

    coot_script += ["mol1 = " + COOT_load_pdb.replace('<pdbin>',pdbin)]
    coot_script += ["mtz1 = " + COOT_load_mtz.replace('<mtzin>',mtzin)]

    coot_script += ["print('MOLLLL: ' + str(mol1))"]
    coot_script += ["print('MTZZZZ: ' + str(mtz1))"]

    coot_script += [COOT_exit]

    with open(scriptname, 'w') as scrip:
        scrip.write('\n'.join(coot_script))

    COOT = run_coot(scriptname)

    return COOT


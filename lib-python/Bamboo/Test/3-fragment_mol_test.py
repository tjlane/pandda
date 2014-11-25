#! /usr/local/python/python2.7.3-64bit/bin/python

import os

from Bamboo.Rdkit.Bonds.Fragment import break_and_rename_mol_to_file

def test_fragmenting_method():
    try:
        whole_file = '/home/npearce/bin/PANDDAs/Bamboo/Test/FRAGMENTTESTFILES/fragment_ligand_whole.pdb'
        fragd_file = '/home/npearce/bin/PANDDAs/Bamboo/Test/FRAGMENTTESTFILES/fragment_ligand_fragged.pdb'
        refce_file = '/home/npearce/bin/PANDDAs/Bamboo/Test/FRAGMENTTESTFILES/fragment_ligand_fragged_ref.pdb'
        if os.path.exists(fragd_file):
            os.remove(fragd_file)
        # Fragment mol (should look like reference mol)
        outfile=break_and_rename_mol_to_file(whole_file, fragd_file)
        ref_block = open(refce_file,'r').read()
        new_block = open(fragd_file,'r').read()
        # Assert
        assert ref_block == new_block
        # Return
        print('SUCCESS: FRAGMENT MOLECULE')
        return 0
    except:
        print('FAILURE: FRAGMENT MOLECULE')
        raise
        return 1

if __name__ == '__main__':

    failure = 0
    failure += test_fragmenting_method()

    print ' ==> Test 3: {!s} failures.'.format(failure)

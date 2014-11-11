#! /usr/local/python/python2.7.3-64bit/bin/python

import os

from Bamboo.Macro.Molecule import MacroMol

def test_structure_reading_1():
    try:
        whole_file = '/home/npearce/bin/PANDDAs/Bamboo/Test/FRAGMENTTESTFILES/fragment_ligand_whole.pdb'
        mol = MacroMol(whole_file)
        r = mol[0][0]
        # Assert
        assert len(mol.getChains())==1, 'FAILURE: test_structure_reading_1'
        assert len(mol.getResidues())==1, 'FAILURE: test_structure_reading_1'
        assert r.get_centroid() == (26.66569565217391, 50.656086956521726, 0.2699565217391303), 'FAILURE: test_structure_reading_1'
        assert r.get_mean_bfactor() == 38.31304347826087, 'FAILURE: test_structure_reading_1'
        assert r.get_mean_occupancy() == 1.0, 'FAILURE: test_structure_reading_1'
        # Return
        print('SUCCESS: RESIDUE METHODS 1')
        return 0
    except:
        print('FAILURE: RESIDUE METHODS 2')
        raise
        return 1

def test_structure_reading_2():
    try:
        whole_file = '/home/npearce/bin/PANDDAs/Bamboo/Test/STRUCTURETESTFILES/varied_occupancy_ligand_test.pdb'
        mol = MacroMol(whole_file)
        r = mol[0][0]
        # Assert
        assert len(mol.getChains())==1
        assert len(mol.getResidues())==1
        assert r.get_centroid() == (26.66569565217391, 50.656086956521726, 0.2699565217391303)
        assert r.get_mean_bfactor() == 25.00
        assert r.get_mean_occupancy() == 0.39999999999999997
        # Return
        print('SUCCESS: RESIDUE METHODS 2')
        return 0
    except:
        print('FAILURE: RESIDUE METHODS 2')
        raise
        return 1

if __name__ == '__main__':

    failure = 0
    failure += test_structure_reading_1()
    failure += test_structure_reading_2()

    print ' ==> Test 4: {!s} failures.'.format(failure)

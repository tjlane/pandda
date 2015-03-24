#! /usr/local/python/python2.7.3-64bit/bin/python

import os, sys

from Bamboo.Analysis.Angles import calculate_dihedral_angle_rmsd_from_file
from Bamboo.Analysis.Coords import calculate_rmsd, calculate_centroid_difference

class BambooTestError(Exception):
        pass

def CalcRMSDTest():
    try:
        original_lig = '/home/npearce/bin/PANDDAs/Bamboo/Test/RMSDTESTFILES/lig-orig.pdb'
        shifted_lig = '/home/npearce/bin/PANDDAs/Bamboo/Test/RMSDTESTFILES/lig-shift-1.pdb'
        # Assert
        assert calculate_rmsd(original_lig, original_lig) == 0.0
        assert calculate_rmsd(shifted_lig, shifted_lig) == 0.0
        assert calculate_rmsd(shifted_lig, original_lig) == 1.0
        # Return
        print('SUCCESS: CALCULATE RMSD')
        return 0
    except:
        print('FAILURE: CALCULATE RMSD')
        return 1

def CalcCentroidTest():
    try:
        original_lig = '/home/npearce/bin/PANDDAs/Bamboo/Test/RMSDTESTFILES/lig-orig.pdb'
        shifted_lig = '/home/npearce/bin/PANDDAs/Bamboo/Test/RMSDTESTFILES/lig-shift-1.pdb'
        # Assert
        assert calculate_centroid_difference(original_lig, original_lig) == 0.0
        assert calculate_centroid_difference(shifted_lig, shifted_lig) == 0.0
        assert calculate_centroid_difference(shifted_lig, original_lig) == 0.9999999999999987
        # Return
        print('SUCCESS: CALCULATE CENTROID')
        return 0
    except:
        print('FAILURE: CALCULATE CENTROID')
        return 1

def CalcAngleRMSDTest():
    try:
        # Test 1 (should be ~0)
        original_lig = '/home/npearce/bin/PANDDAs/Bamboo/Test/RMSDTESTFILES/lig-orig.pdb'
        shifted_lig = '/home/npearce/bin/PANDDAs/Bamboo/Test/RMSDTESTFILES/lig-shift-1.pdb'
        # Assert
        assert calculate_dihedral_angle_rmsd_from_file(original_lig, shifted_lig) < 1e-10
        # Test 2 (should be ~0)
        conf1 = '/home/npearce/bin/PANDDAs/Bamboo/Test/RMSDTESTFILES/lig-conf-1.pdb'
        conf2 = '/home/npearce/bin/PANDDAs/Bamboo/Test/RMSDTESTFILES/lig-conf-2.pdb'
        # Assert
        assert round(calculate_dihedral_angle_rmsd_from_file(conf1, conf2)) == 59
        # Return
        print('SUCCESS: CALCULATE ANGLE RMSD')
        return 0
    except:
        print('FAILURE: CALCULATE ANGLE RMSD')
        return 1

if __name__ == '__main__':
    failure = 0
    failure += CalcRMSDTest()
    failure += CalcCentroidTest()
    failure += CalcAngleRMSDTest()

    print ' ==> Test 2: {!s} failures.'.format(failure)

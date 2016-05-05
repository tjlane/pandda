import math

from bamboo.rdkit_utils.mol import check_pdb_readable
from bamboo.rdkit_utils.bonds.dihedral import calculate_dihedral_angle_differences

def calculate_dihedral_angle_rmsd_from_file(file1, file2):
    """Calculate angle rmsd between structures"""
    return calculate_dihedral_angle_rmsd(check_pdb_readable(file1),check_pdb_readable(file2))

def calculate_dihedral_angle_rmsd(mol1, mol2):
    """Calculate angle rmsd between structures"""

    try:
        # Get the atoms and the differences
        atoms, diffs = calculate_dihedral_angle_differences(mol1, mol2)
    except EqualityError as Err:
        print Err
        return None

    # List of rmsds for different ways of matching
    rmsds = []
    # Iterate through the different ways of matching the atoms
    for diff in diffs:
        if not diff:
            continue
        else:
            sqdiff = [d**2 for d in diff]
            # Calculate the RMSD
            rmsds.append(math.sqrt(sum(sqdiff)/len(sqdiff)))
    # Return the minimum rmsd
    if rmsds:
        return min(rmsds)
    else:
        return None


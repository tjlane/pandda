#! /usr/local/python/python2.7.3-64bit/bin/python

# Contains functions to filter BioPython Structure Objects

# Biopython
from Bio.PDB.Polypeptide import aa3, aa1

#:::::::::::::::::::::::::::::::::#
# ############################### #
# ###      PDB Functions      ### #
# ############################### #
#:::::::::::::::::::::::::::::::::#

# 3 Letter Codes of Amino Acids : from Bio.PDB.Polypeptide import aa3
# 1 Letter Codes of Amino Acids : from Bio.PDB.Polypeptide import aa1

#------------ Biopython PDB Functions ------------

def GetNonAminoAcids(structure, KeepWaters=False):
    """Get all of the Non-Amino-Acids in a Structure Object"""

    # Get 'Residues' from the Structure
    residues = structure.get_residues()
    # Create Filter
    filter = aa3 + ['HOH']*(not KeepWaters)
    # Filter out Amino Acids
    filtered = [res for res in residues if res.get_resname() not in filter]

    return filtered

def FilterResiduesByList(structure, list, KeepObjectsInList=False):
    """Filter Residues in List from the Structure (List needs to be 3-letter codes)"""

    # Get Residues
    residues = structure.get_residues()
    # Keep List / Remove List?
    if KeepObjectsInList:
        # Keep Residues in the List
        filtered = [res for res in residues if res.get_resname() in list]
    else:
        # Filter out Residues in List
        filtered = [res for res in residues if res.get_resname() not in list]

    return filtered


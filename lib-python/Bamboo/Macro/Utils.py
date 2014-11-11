#! /usr/local/python/python2.7.3-64bit/bin/python

from Bamboo.Macro.Molecule import MacroMol

def get_residue_labels(pdbpath):
    """Takes the PDB file and returns a list of the residues labels in the form (res.resname, res.chain, res.resnum, res.inscode) (e.g. for edstats)"""

    labels = []
    # Create a structure object
    structure = MacroMol(pdbpath)
    # Iterate through and create labels
    for res in structure.getResidues():
        labels.append(res.get_res_id())

    return labels

def get_mean_occupancy(pdbpath, type='aminos'):
    """Gets the mean occupancy of the structure in pdbpath. Changing `type` changes the residues that are used to calculate the occupancy. Choices are ['all','aminos','ions','waters']"""

    m = MacroMol(pdbpath)

    occupancies = []

    if type=='all':
        residues = m.getResidues()
    elif type=='aminos':
        residues = m.getAminos()
    elif type=='waters':
        residues = m.getWaters()
    elif type=='ions':
        residues = m.getIons()

    [occupancies.append(res.get_mean_occupancy()) for res in residues]

    mean_occ = sum(occupancies)/len(occupancies)

    return mean_occ

def get_residue_occupancies(pdbpath):
    """Gets the occupancies of individual residues. Returns a dict"""

    m = MacroMol(pdbpath)

    occupancies = {}

    for res in m.getResidues():
        occupancies[res.get_res_id()] = res.get_mean_occupancy()

    return occupancies


from bamboo.macro.molecule import MacroMol, MacroMolError
from bamboo.rdkit_utils.smile import find_structure_matches

def get_residue_labels(pdb_file):
    """Takes the PDB file and returns a list of the residues labels in the form (res.resname, res.chain, res.resnum, res.inscode) (e.g. for edstats)"""

    labels = []
    # Create a structure object
    structure = MacroMol(pdb_file)
    # Iterate through and create labels
    for res in structure.get_residues():
        labels.append(res.get_res_id())

    return labels

def get_mean_occupancy(pdb_file, type='aminos'):
    """Gets the mean occupancy of the structure in pdb_file. Changing `type` changes the residues that are used to calculate the occupancy. Choices are ['all','aminos','ions','waters']"""

    m = MacroMol(pdb_file)

    occupancies = []

    if type=='all':
        residues = m.get_residues()
    elif type=='aminos':
        residues = m.get_aminos()
    elif type=='waters':
        residues = m.get_waters()
    elif type=='ions':
        residues = m.get_ions()

    [occupancies.append(res.get_mean_occupancy()) for res in residues]

    mean_occ = sum(occupancies)/len(occupancies)

    return mean_occ

def get_residue_occupancies(pdb_file):
    """Gets the occupancies of individual residues. Returns a dict"""

    m = MacroMol(pdb_file)

    occupancies = {}

    for res in m.get_residues():
        occupancies[res.get_res_id()] = res.get_mean_occupancy()

    return occupancies

def extract_compound_smiles(pdbin, find_boring_compounds=False):
    """Takes a PDB file and returns smiles strings of all of the compounds in the file"""

    try:
        m = MacroMol(pdbin)
    except MacroMolError:
        raise

    cpds = m.get_unknowns()

    # ADD FILTERING STEP!
    if find_boring_compounds:
        cpds.extend(m.get_solvent())

    [c.get_smiles() for c in cpds]
    return [c.smile for c in cpds if c.smile]


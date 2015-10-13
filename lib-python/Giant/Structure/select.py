import iotbx.pdb

def get_backbone_atoms(atom_group):
    """Get the backbone atoms for the atom_group"""

    # Get the right interpreter for this residue
    intr = iotbx.pdb.protein_atom_name_interpreters[atom_group.resname]
    # Match the atom names
    matches = intr.match_atom_names(atom_group.atoms().extract_name())
    # Select the backbone atoms
    c = matches.expected['C']

    return






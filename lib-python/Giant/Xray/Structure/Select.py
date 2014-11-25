from scitbx.array_family import flex
#from mmtbx.command_line.super import extract_sequence_and_sites

def get_calpha_sites(input_obj, structure_obj):
    """Gets the coordinates of alpha carbons from the structure"""
    return flex.vec3_double([a.xyz for a in get_atom_selection(input_obj, atom_names=[' CA '], atom=True, hetero=False)])

def get_backbone_sites(input_obj, structure_obj):
    """Gets the coordinates of backbone atoms"""
    return flex.vec3_double([a.xyz for a in get_atom_selection(input_obj, atom_names=[' C  ',' N  ',' CA ',' O  '], atom=True, hetero=False)])

def get_atom_selection(input_obj, chain_names=[], res_names=[], atom_names=[], atom=True, hetero=False):
    """Iterate through atoms in input_obj and select atoms matching the criteria"""

    all_atoms = input_obj.atoms_with_labels()

    # Select Hetero
    if (atom == True) and (hetero == True):
        filt_atoms = all_atoms
    elif (atom == True) and (hetero == False):
        filt_atoms = [a for a in all_atoms if a.hetero == False]
    elif (atom == False) and (hetero == True):
        filt_atoms = [a for a in all_atoms if a.hetero == True]
    else:
        raise Exception('No Atoms Selected')

    if chain_names:
        filt_atoms = [a for a in filt_atoms if a.chain_id in chain_names]

    if res_names:
        filt_atoms = [a for a in filt_atoms if a.resname in res_names]

    if atom_names:
        filt_atoms = [a for a in filt_atoms if a.name in atom_names]

    return filt_atoms

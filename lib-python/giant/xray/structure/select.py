from scitbx.array_family import flex
#from mmtbx.command_line.super import extract_sequence_and_sites

def get_calpha_sites(input_hierarchy):
    """Gets the coordinates of alpha carbons from the structure"""
    return flex.vec3_double([a.xyz for a in get_atom_selection(input_hierarchy, atom_names=[' CA '], atom=True, hetero=False)])

def get_backbone_sites(input_hierarchy, cbeta=True):
    """Gets the coordinates of backbone atoms"""
    if cbeta: atoms = [' C  ',' N  ',' CA ',' O  ',' CB ']
    else:     atoms = [' C  ',' N  ',' CA ',' O  ']
    return flex.vec3_double([a.xyz for a in get_atom_selection(input_hierarchy, atom_names=[' C  ',' N  ',' CA ',' O  '], atom=True, hetero=False)])

def get_atom_selection(input_hierarchy, chain_names=None, res_names=None, atom_names=None, atom=True, hetero=False, proteinonly=True):
    """Iterate through atoms in input_hierarchy and select atoms matching the criteria"""

    if chain_names is None: chain_names = []
    if res_names   is None: res_names   = []
    if atom_names  is None: atom_names  = []

    filt_atoms = input_hierarchy.atoms_with_labels()

    # Select Hetero
    if (atom == True) and (hetero == True):
        filt_atoms = filt_atoms
    elif (atom == True) and (hetero == False):
        filt_atoms = [a for a in filt_atoms if a.hetero == False]
    elif (atom == False) and (hetero == True):
        filt_atoms = [a for a in filt_atoms if a.hetero == True]
    else:
        raise Exception('No Atoms Selected')

    # Get protein chains
    if proteinonly:
        prot_chains = [ch.id for ch in input_hierarchy.chains() if ch.is_protein()]
        filt_atoms = [a for a in filt_atoms if a.chain_id in prot_chains]

    # Get selected chains
    if chain_names:
        filt_atoms = [a for a in filt_atoms if a.chain_id in chain_names]

    # Get selected residues
    if res_names:
        filt_atoms = [a for a in filt_atoms if a.resname in res_names]

    # Get particular atoms
    if atom_names:
        filt_atoms = [a for a in filt_atoms if a.name in atom_names]

    return filt_atoms

from scitbx.array_family import flex

########################################################################################

def non_h(hierarchy, cache=None, copy=True):
    if not cache: cache=hierarchy.atom_selection_cache()
    sel = cache.selection('(not element H)')
    return hierarchy.select(sel, copy_atoms=copy)

def protein(hierarchy, cache=None, copy=True):
    if not cache: cache=hierarchy.atom_selection_cache()
    sel = cache.selection('(not element H) and pepnames')
    return hierarchy.select(sel, copy_atoms=copy)

def calphas(hierarchy, cache=None, copy=True):
    if not cache: cache=hierarchy.atom_selection_cache()
    sel = cache.selection('(not element H) and pepnames and (name CA)')
    return hierarchy.select(sel, copy_atoms=copy)

def backbone(hierarchy, cache=None, copy=True):
    if not cache: cache=hierarchy.atom_selection_cache()
    sel = cache.selection('(not element H) and pepnames and (name C or name CA or name N or name O)')
    return hierarchy.select(sel, copy_atoms=copy)

def sidechains(hierarchy, cache=None, copy=True):
    if not cache: cache=hierarchy.atom_selection_cache()
    sel = cache.selection('(not element H) and pepnames and not (name C or name CA or name N or name O)')
    return hierarchy.select(sel, copy_atoms=copy)

########################################################################################

def extract_backbone_atoms(residue):
    """Extract N, CA, C triplets for a residue"""
    residue_int = residue.residue_name_plus_atom_names_interpreter().atom_name_interpretation.expected
    residue_ats = residue.atoms().build_dict()
    n  = residue_ats[residue_int['N'][0]]
    ca = residue_ats[residue_int['CA'][0]]
    c  = residue_ats[residue_int['C'][0]]
    return (n, ca, c)

########################################################################################

#def get_backbone_atoms(atom_group):
#    """Get the backbone atoms for the atom_group"""
#
#    # Get the right interpreter for this residue
#    intr = iotbx.pdb.protein_atom_name_interpreters[atom_group.resname]
#    # Match the atom names
#    matches = intr.match_atom_names(atom_group.atoms().extract_name())
#    # Select the backbone atoms
#    c = matches.expected['C']
#
#    return

########################################################################################

# Deprecated
def get_calpha_sites(input_hierarchy):
    """Gets the coordinates of alpha carbons from the structure"""
    return flex.vec3_double([a.xyz for a in get_atom_selection(input_hierarchy, atom_names=[' CA '], atom=True, hetero=False)])

# Deprecated
def get_backbone_sites(input_hierarchy, cbeta=True):
    """Gets the coordinates of backbone atoms"""
    if cbeta: atoms = [' C  ',' N  ',' CA ',' O  ',' CB ']
    else:     atoms = [' C  ',' N  ',' CA ',' O  ']
    return flex.vec3_double([a.xyz for a in get_atom_selection(input_hierarchy, atom_names=[' C  ',' N  ',' CA ',' O  '], atom=True, hetero=False)])

########################################################################################

# Deprecated
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


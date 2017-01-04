
import iotbx.pdb
from scitbx.array_family import flex

from giant.structure.sequence import align_sequences_default

####################################################################################
###                             HIERARCHY FUNCTIONS                              ###
####################################################################################

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

####################################################################################
###                              CHAIN FUNCTIONS                                 ###
####################################################################################

def _truncate_by_idx(chn, i_rg_sel):
    """Truncate a chn by a list of indices"""
    assert len(i_rg_sel)>0, 'No residues selected for trunction!'
    new_h = iotbx.pdb.hierarchy.new_hierarchy_from_chain(chn)
    cache = new_h.atom_selection_cache()
    rg_all = list(new_h.residue_groups())
    resid_sel = ' or '.join(['resid {}'.format(rg_all[i].resid()) for i in i_rg_sel])
    return new_h.select(cache.selection(resid_sel)).only_chain()

def common_residues(chn_1, chn_2):
    """
    Truncates input chains to the common set of residues in chn_1 and chn_2 after sequence alignment.
    Returns both truncated chains.
    """
    # Apply default quick alignment
    alignment = align_sequences_default(seq_a=chn_1.as_sequence(),
                                        seq_b=chn_2.as_sequence())
    # Flags for which residues to use
    m_seq_1, m_seq_2 = alignment.exact_match_selections()
    assert len(m_seq_1) == len(m_seq_2), 'Something has gone wrong: these should be the same length!'
    assert (max(m_seq_1)<len(alignment.a)) and (max(m_seq_2)<len(alignment.b)), 'Something has gone wrong: selecting residue index greater than chain length'
    # Truncate down to the identical selections
    out_c_1 = _truncate_by_idx(chn_1, m_seq_1)
    out_c_2 = _truncate_by_idx(chn_2, m_seq_2)
    return out_c_1, out_c_2

def extract_conformer(chn, conf=['','A']):
    h = iotbx.pdb.hierarchy.new_hierarchy_from_chain(chn)
    c = h.atom_selection_cache()
    s = c.selection(' or '.join(['altid "{:1}"'.format(i) for i in conf]))
    return h.select(s, copy_atoms=False).only_chain()

def complete_backbone(chn, conf=['','A']):
    """
    Truncate a chn to only residues with a complete set of backbone atoms.
    Return chn object composing of residues with a complete set of backbone atoms
    """
    # If string given, put in list
    if isinstance(conf, str): conf = [conf]
    # If nothing given, default to blank conformer
    if not conf: conf = ['']
    # Check that there's a maximum of one non-blank conformer ID
    assert (conf == ['']) or ( (len(conf)-conf.count(''))==1 ), 'Supply up to one non-blank conformer (e.g. ["","A"], [""] or ["A"]): {}'.format(conf)
    # Iterate through and record indices of complete residues
    i_sel = []
    for i_rg,rg in enumerate(chn.residue_groups()):
        confs = [c for c in rg.conformers() if c.altloc in conf]
        if not confs: continue
        assert len(confs) == 1, 'Something has gone wrong here. This should never be longer than 1.'
        # Extract the residue (as this has atom name interpretation)
        res = confs[0].only_residue()
        # Check all backbone atoms are there and append
        if check_backbone_atoms(res): i_sel.append(i_rg)
    return _truncate_by_idx(chn, i_sel)

####################################################################################
###                             RESIDUE FUNCTIONS                                ###
####################################################################################

def check_backbone_atoms(residue):
    """Returns True if N,CA,C,O are present, else False"""
    return not bool(residue.residue_name_plus_atom_names_interpreter() \
                        .atom_name_interpretation.missing_atom_names().intersection(['N','CA','C','O']))

def extract_atom(residue, atom='CA'):
    i = residue.residue_name_plus_atom_names_interpreter().atom_name_interpretation.expected
    d = residue.atoms().build_dict()
    return d[i[atom][0]]

def extract_backbone_atoms(residue):
    """Extract N, CA, C triplets for a residue"""
    residue_int = residue.residue_name_plus_atom_names_interpreter().atom_name_interpretation.expected
    residue_ats = residue.atoms().build_dict()
    n  = residue_ats[residue_int['N'][0]]
    ca = residue_ats[residue_int['CA'][0]]
    c  = residue_ats[residue_int['C'][0]]
    return (n, ca, c)

####################################################################################
###                              DISTANCE FUNCTIONS                              ###
####################################################################################

def find_nearest_atoms(atoms, query):
    """Find the nearest atom for each coordinate in query"""
    label_indices = find_closest_points(reference_points=atoms.extract_xyz(), query=query)
    return [atoms[i] for i in indices]

def find_closest_points(points, query):
    """Find and return the index of the closest point in points for each coordinate in query"""
    tree = scipy.spatial.KDTree(data=points)
    nn_dists, nn_groups = tree.query(query)
    return nn_groups

####################################################################################
###                           DEPRECATED FUNCTIONS                               ###
####################################################################################

# Deprecated
def get_calpha_sites(input_hierarchy):
    """Gets the coordinates of alpha carbons from the structure"""
    print 'This function is deprecated and will be deleted'
    return flex.vec3_double([a.xyz for a in get_atom_selection(input_hierarchy, atom_names=[' CA '], atom=True, hetero=False)])
# Deprecated
def get_backbone_sites(input_hierarchy, cbeta=True):
    """Gets the coordinates of backbone atoms"""
    print 'This function is deprecated and will be deleted'
    if cbeta: atoms = [' C  ',' N  ',' CA ',' O  ',' CB ']
    else:     atoms = [' C  ',' N  ',' CA ',' O  ']
    return flex.vec3_double([a.xyz for a in get_atom_selection(input_hierarchy, atom_names=[' C  ',' N  ',' CA ',' O  '], atom=True, hetero=False)])
# Deprecated
def get_atom_selection(input_hierarchy, chain_names=None, res_names=None, atom_names=None, atom=True, hetero=False, proteinonly=True):
    """Iterate through atoms in input_hierarchy and select atoms matching the criteria"""
    print 'This function is deprecated and will be deleted'
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


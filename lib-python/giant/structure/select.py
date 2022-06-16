import giant.logs as lg
logger = lg.getLogger(__name__)

import re

import numpy
import scipy.spatial

import iotbx.pdb
from scitbx.array_family import flex
from libtbx.utils import null_out, Sorry, Failure

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

def backbone(hierarchy, cache=None, copy=True, cbeta=False):
    if not cache: cache=hierarchy.atom_selection_cache()
    if cbeta:   sel = cache.selection('pepnames and (name C or name CA or name N or name O or name CB)')
    else:       sel = cache.selection('pepnames and (name C or name CA or name N or name O)')
    return hierarchy.select(sel, copy_atoms=copy)

def sidechains(hierarchy, cache=None, copy=True, cbeta=True):
    if not cache: cache=hierarchy.atom_selection_cache()
    if cbeta:   sel = cache.selection('(not element H) and pepnames and not (name C or name CA or name N or name O)')
    else:       sel = cache.selection('(not element H) and pepnames and not (name C or name CA or name N or name O or name CB)')
    return hierarchy.select(sel, copy_atoms=copy)

def non_water(hierarchy, cache=None, copy=True):
    if not cache: cache=hierarchy.atom_selection_cache()
    sel = cache.selection('(not resname HOH)')
    return hierarchy.select(sel, copy_atoms=copy)

####################################################################################

def default_secondary_structure_selections(hierarchy):
    hierarchy = protein(hierarchy, copy=True)
    hierarchy.reset_atom_i_seqs()
    from mmtbx.secondary_structure import dssp
    result = dssp.dssp(hierarchy, log=null_out(), out=null_out()).get_annotation().as_atom_selections()
    if len(result)==0:
        raise Sorry('No secondary structure elements were identified by dssp')
    return result

def default_secondary_structure_selections_filled(hierarchy):
    """Return secondary structure selections and fill gaps with new selections"""

    hierarchy = protein(hierarchy, copy=True)

    logger.debug('>>> Creating secondary structure selections for input hierarchy <<<')

    # Get automatic secondary structure elements
    auto_sel = default_secondary_structure_selections(hierarchy=hierarchy)

    logger.debug('\n>> Default Selections (from dssp):\n\t{}'.format('\n\t'.join(auto_sel)))

    # Extract chain and residue numbers
    sel_regex = re.compile("chain '(.*?)' and resid (.*?). through (.*)")
    auto_sel_proc = [sel_regex.findall(s)[0] for s in auto_sel]

    # Output template and list
    output_sel_all = []

    # Extract the chain IDs
    chain_ids = numpy.unique(list(zip(*auto_sel_proc))[0])

    # Iterate through chains
    for c_id in chain_ids:

        # Extract residue start and end numbers
        c_sels = sorted([list(map(int, s[1:])) for s in auto_sel_proc if s[0]==c_id])
        # Extract chain
        c_obj  = [c for c in hierarchy.chains() if (c.is_protein() and c.id==c_id)]
        assert len(c_obj) == 1
        c_obj = c_obj[0]
        # Get the first and last residues of the chain
        c_start = min([r.resseq_as_int() for r in c_obj.residue_groups()])
        c_end   = max([r.resseq_as_int() for r in c_obj.residue_groups()])
        # Create boolean selection for residues in the chain
        n_res = c_end-c_start+1
        # Create residue mask with group numbers; mark missing residues with -1
        t_sel = -1*numpy.ones(n_res, dtype=int)
        t_sel[[r.resseq_as_int()-c_start for r in c_obj.residue_groups()]] = 0

        # log strings
        logger.debug('\n'.join((
            '>> Processing chain {}'.format(c_id),
            '>> Sorted selections for chain {}:\n\t{}'.format(c_id, '\n\t'.join(['{:>4d} -> {:>4d}'.format(*s) for s in c_sels])),
            '> Chain start: {}\n> Chain end: {}'.format(c_start, c_end),
            '> {} missing residues in chain {}'.format(sum(t_sel==-1), c_id),
            '\n>> Removing overlapping groups:',
        )))

        # First pass: Iterate through and place group numbers into mask
        for i,g in enumerate(c_sels):

            # Convert to selections
            n = i+1
            g_sel = numpy.zeros_like(t_sel, dtype=bool)
            g_sel[list(range(g[0]-c_start,g[1]-c_start+1))] = True
            # Multiply the unassigned atoms with this group mask
            new_g_sel = (t_sel==0)*g_sel

            # log strings
            logger.debug('\n'.join((
                'Group: {}'.format(str(g)),
                'Group selection: {} residues'.format(sum(g_sel)),
                'Filtered Group selection: {} residues'.format(sum(new_g_sel)),
                '(after removing residues that are already in other groups)',
            )))

            # Skip if no residues left
            if sum(new_g_sel) == 0:
                logger.debug('No residues for this group that are not already in another group')
                continue

            # Store as this group
            t_sel[new_g_sel] = n

        # Second pass: Process groups & create new groups for gaps
        g_start = 0
        f_sels = []
        logger.debug('>> Creating new non-overlapping groups:')
        for i_this, v_this in enumerate(t_sel):
            # Skip missing residues
            if v_this == -1:
                continue
            # If not in a group, start a group
            if g_start is None:
                g_start = i_this
            # When in a group, check if next is still in group
            v_next = t_sel[i_this+1] if (i_this+1)<len(t_sel) else None
            if v_this != v_next:
                f_sels.append([g_start+c_start, i_this+c_start])
                g_start = None
        logger.debug('Filtered & Filled selections for chain {}:\n\t{}'.format(c_id, '\n\t'.join(['{:>4d} -> {:>4d}'.format(*s) for s in f_sels])))

        # Third pass: Merge small groups
        o_sels = []
        logger.debug('>> Merging small groups with neighbours where possible:')
        for i, (start, end) in enumerate(f_sels):
            g_size = end - start + 1
            if g_size < 3:
                logger.debug('Group ({},{}) is less than three residues ({} residues)'.format(start,end,g_size))
                g_before = int((not i==0)
                        and (f_sels[i-1] is not None)
                        and (f_sels[i-1][1]==(start-1)))
                g_after = int((not i+1==len(f_sels))
                        and (f_sels[i+1] is not None)
                        and (f_sels[i+1][0]==(end+1)))
                # Decide how/if to split group between neighbours
                if (g_before or g_after):
                    logger.debug('Splitting group between neighbouring groups')
                    n_after = g_after * int(float(g_size)/float(g_before+g_after))
                    n_before = g_size - n_after
                    assert n_before + n_after == g_size
                    if n_before:
                        logger.debug('Adding {} residues to previous group'.format(n_before))
                        f_sels[i-1][1] += n_before
                    if n_after:
                        logger.debug('Adding {} residues to next group'.format(n_after))
                        f_sels[i+1][0] -= n_after
                    # Remove the group
                    f_sels[i] = None
                elif g_size == 1:
                    logger.debug('Group of one residue with no neighbours -- removing group')
                    # Remove the group
                    f_sels[i] = None
                else:
                    logger.debug('No neighbouring groups -- leaving as group of two residues')

        # Remove the none groups
        o_sels = [(c_id, v[0], v[1]) for v in f_sels if v is not None]
        logger.debug('Merged selections for chain {}:\n\t{}'.format(c_id, '\n\t'.join(['{:>4d} -> {:>4d}'.format(*s[1:]) for s in o_sels])))

        # Append to overall list
        output_sel_all.extend(o_sels)

    # Sort output for fun
    sel_template = "chain '{}' and resid {:d}  through {:d}"
    output_sel_all = sorted(output_sel_all)
    output_sel = [sel_template.format(*v) for v in output_sel_all]
    logger.debug('Processed Selections:\n\t{}'.format('\n\t'.join(output_sel)))

    return output_sel

####################################################################################

def get_select_function(selection):
    if selection is None: return lambda x: x
    functions = ['non_h','protein','calphas','backbone','sidechains','non_water']
    assert selection in functions, 'Function not found. Choices: {}'.format(','.join(functions))
    func = eval(selection)
    return func

####################################################################################

def sel_altloc(hierarchy, altlocs=['','A'], cache=None, copy=True):
    if not cache: cache=hierarchy.atom_selection_cache()
    sel = cache.selection('(not element H) and ({})'.format(' or '.join(['altid "{:1}"'.format(c) for c in altlocs])))
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

def chain_conformer(chn, altlocs=['','A']):
    if isinstance(altlocs, str): altlocs = [altlocs]
    if not altlocs: altlocs = ['']
    assert (altlocs == ['']) or ( (len(altlocs)-altlocs.count(''))==1 ), 'Supply up to one non-blank altloc (e.g. ["","A"], [""] or ["A"]): {}'.format(altlocs)
    return sel_altloc(hierarchy=iotbx.pdb.hierarchy.new_hierarchy_from_chain(chn), altlocs=altlocs).only_chain()

def common_residues(chn_1, chn_2):
    """
    Truncates input chains to the common set of residues in chn_1 and chn_2 after sequence alignment.
    Returns both truncated chains.
    """

    # Apply default quick alignment
    from giant.structure.sequence import align_sequences_default
    alignment = align_sequences_default(seq_a=chn_1.as_sequence(), seq_b=chn_2.as_sequence())
    # Flags for which residues to use
    m_seq_1, m_seq_2 = alignment.exact_match_selections()
    assert len(m_seq_1) == len(m_seq_2),          'Something has gone wrong: these should be the same length!'
    assert len(alignment.a) == len(alignment.b),  'Something has gone wrong: these should be the same length!'
    assert max(m_seq_1)<len(chn_1.as_sequence()), 'Something has gone wrong: selecting residue index greater than chain length'
    assert max(m_seq_2)<len(chn_2.as_sequence()), 'Something has gone wrong: selecting residue index greater than chain length'
    # Truncate down to the identical selections
    out_c_1 = _truncate_by_idx(chn_1, m_seq_1)
    out_c_2 = _truncate_by_idx(chn_2, m_seq_2)
    return out_c_1, out_c_2

def complete_backbone(chn, altlocs=['','A']):
    """
    Truncate a chn to only residues with a complete set of backbone atoms.
    Return chn object composing of residues with a complete set of backbone atoms
    """
    # Extract only the requested conformers
    chn = chain_conformer(chn, altlocs=altlocs)
    # Iterate through and record indices of complete residues
    i_sel = []
    for i_rg,rg in enumerate(chn.residue_groups()):
        confs = [c for c in rg.conformers() if c.altloc in altlocs]
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
    intrp = residue.residue_name_plus_atom_names_interpreter().atom_name_interpretation
    if not intrp: return False
    return not bool(intrp.missing_atom_names().intersection(['N','CA','C','O']))

def extract_backbone_atoms(residue, atoms=('N','CA','C')):
    """Extract N, CA, C triplets for a residue"""
    residue_int = residue.residue_name_plus_atom_names_interpreter().atom_name_interpretation.expected
    residue_ats = residue.atoms().build_dict()
    return tuple([residue_ats[residue_int[a][0]] for a in atoms])

def extract_atom(residue, atom='CA'):
    i = residue.residue_name_plus_atom_names_interpreter().atom_name_interpretation.expected
    d = residue.atoms().build_dict()
    return d[i[atom][0]]

####################################################################################
###                              DISTANCE FUNCTIONS                              ###
####################################################################################

def find_nearest_atoms(atoms, query):
    """Find the nearest atom for each coordinate in query"""
    if isinstance(atoms, list): xyz = [a.xyz for a in atoms]
    else:                       xyz = atoms.extract_xyz()
    label_indices = find_closest_points(points=xyz, query=query)
    return [atoms[i] for i in label_indices]

def find_closest_points(points, query):
    """Find and return the index of the closest point in points for each coordinate in query"""
    tree = scipy.spatial.KDTree(data=points)
    nn_dists, nn_groups = tree.query(query)
    return nn_groups



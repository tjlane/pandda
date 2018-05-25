import re

import numpy
import scipy.spatial

import iotbx.pdb
from scitbx.array_family import flex
from libtbx.utils import null_out

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

def non_water(hierarchy, cache=None, copy=True):
    if not cache: cache=hierarchy.atom_selection_cache()
    sel = cache.selection('(not resname HOH)')
    return hierarchy.select(sel, copy_atoms=copy)

####################################################################################

def default_secondary_structure_selections(hierarchy):
    from mmtbx.secondary_structure import dssp
    return dssp.dssp(hierarchy, log=null_out(), out=null_out()).get_annotation().as_atom_selections()

def default_secondary_structure_selections_filled(hierarchy):
    """Return secondary structure selections and fill gaps with new selections"""

    # Get automatic secondary structure elements
    auto_sel = default_secondary_structure_selections(hierarchy=hierarchy)
    # Extract chain and residue numbers
    sel_regex = re.compile("chain '(.*?)' and resid (.*?). through (.*)")
    auto_sel_proc = [sel_regex.findall(s)[0] for s in auto_sel]

    # Output template and list
    sel_template = "chain '{}' and resid {:d}  through {:d} "
    output_sel_all = []

    # Extract the chain IDs
    chain_ids = numpy.unique(zip(*auto_sel_proc)[0])

    # Iterate through chains
    for c_id in chain_ids:
        # Output list for this chain
        output_sel = []
        # Extract residue start and end numbers
        c_sels = sorted([map(int, s[1:]) for s in auto_sel_proc if s[0]==c_id])
        # Extract chain
        c_obj  = [c for c in hierarchy.chains() if (c.is_protein() and c.id==c_id)]
        assert len(c_obj) == 1
        c_obj = c_obj[0]
        # Get the first and last residues of the chain
        c_start = min([r.resseq_as_int() for r in c_obj.residue_groups()])
        c_end   = max([r.resseq_as_int() for r in c_obj.residue_groups()])
        # Get first residue in a ss group
        first_start, first_end = c_sels[0]
        # Find number of unused residues at the beginning of the chain
        unused_front = first_start - c_start
        # If one residue, add to first group
        if unused_front == 1:
            c_sels[0][0] -= 1
        # Else add own selection for these residues
        elif unused_front > 1:
            output_sel.append((c_id, c_start, first_start-1))
        # Process groups for middle of the chain
        for i in xrange(len(c_sels)):
            # Get the start and end of this group
            this_start, this_end = c_sels[i]
            # Perform sanity checks
            assert this_start <= this_end, 'group with negative size?! {} - {}'.format(this_start, this_end)
            # Get boolean for whether this is the final group of the chain
            final_group = (i+1 == len(c_sels))
            # More sanity checks
            if not final_group:
                # Check it is compatible with the next group
                next_start, next_end = c_sels[i+1]
                assert next_start <= next_end,   'group with negative size?! {} - {}'.format(next_start, next_end)
                assert this_end   <  next_end,   'end of next group before end of this group?! {} - {}'.format(this_end, next_end)
                assert this_start <= next_start, 'start of next group before start of this group?! {} - {}'.format(this_start, next_start)
            # Special case: group of one residue -- do not allow
            if (this_start == this_end):
                # delete this group from the input list
                c_sels[i] = None
                # If the last group, just ignore it and break
                if final_group:
                    break
                # Shuffle the group back one to simulate being still in the group before this one
                this_start -= 1
                this_end -= 1
                if output_sel:
                    assert output_sel[-1][2] == this_start, 'Last residue of last group is incorrect: {} != {}'.format(output_sel[-1][2], this_start)
            # Check for gap/overlap to next group
            if not final_group:
                # Get the start of the next group
                next_start, next_end = c_sels[i+1]
                # Check for overlap with next group
                if this_end >= next_start:
                    overlap = (this_end + 1 - next_start)
                    remove_this = int(numpy.ceil(overlap/2.0))
                    remove_next = overlap - remove_this
                    # Check that truncation of this group will leave at least two residues
                    this_remain = (this_end + 2 - this_start) - remove_this
                    # If one or no residues remain for this group
                    if this_remain <= 1:
                        # make this group not exist
                        c_sels[i] = None
                        # merge the two groups (set next start to this start)
                        c_sels[i+1][0] = this_start
                        break
                    # Check that at least two residues will remain in this group
                    assert this_remain >= 2
                    # Truncate this group
                    this_end -= remove_this
                    # Move the start of the next group (and end if necessary)
                    c_sels[i+1] = (next_start+remove_next, max(next_start+remove_next, next_end))
                    # Need to refresh these variables
                    next_start, next_end = c_sels[i+1]
                    assert this_end + 1 == next_start, 'error: this end {}; next start: {}'.format(this_end, next_start)
                # Check to see if the gap between this and the next group is one residue
                elif (next_start - this_end) == 2:
                    # Gap is one residue -- Move the beginning of the next residue
                    c_sels[i+1][0] -= 1
            # Create selection for this group if it has more than two residues
            assert this_end - this_start >= 0, 'error width of group less than one residue: {} - {}'.format(this_start, this_end)
            if (this_end - this_start) >= 1:
                # Recreate selection for this group
                output_sel.append((c_id, this_start, this_end))
            # Add selection for gap to next selection if it exists (and that more than one residue in gap)
            if (not final_group) and ((next_start - this_end) > 2):
                output_sel.append((c_id, this_end+1, next_start-1))

        # Deal with residues at the end of the chain
        last_chain, last_start, last_end = output_sel[-1]
        remain_residues = (c_end - last_end)
        # One residue -- add to last group
        if remain_residues == 1:
            output_sel[-1] = (c_id, last_start, last_end+1)
        # Mroe -- add new group
        elif remain_residues > 1:
            output_sel.append((c_id, last_end+1, c_end))

        # VALIDATE the output for this chain
        assert output_sel[0][1]  == c_start
        assert output_sel[-1][2] == c_end
        for i in xrange(len(output_sel)-1):
            this_chain, this_start, this_end = output_sel[i]
            next_chain, next_start, next_end = output_sel[i+1]
            assert this_chain   == next_chain   # Right chain!
            assert (this_end+1) == next_start   # Groups are contiguous
            assert (this_end - this_start) > 0  # At least two residues

        # Append to overall list
        output_sel_all.extend(output_sel)

    # Sort output for fun
    output_sel_all = sorted(output_sel_all)

    return [sel_template.format(*v) for v in output_sel_all]

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


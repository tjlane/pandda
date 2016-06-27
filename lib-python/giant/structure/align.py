
import itertools
import numpy

import iotbx.pdb
import iotbx.pdb.amino_acid_codes
import mmtbx.alignment

from scitbx.math import superpose
from scitbx.array_family import flex

from giant.structure import make_label
from giant.structure.select import protein

def nearby_coords_bool(query, coords, cutoff):
    """Find all points in coords within cutoff of query. Return boolean selection"""
    assert isinstance(coords, flex.vec3_double)
    return (coords-query).dot() < (cutoff**2)

def perform_flexible_alignment(mov_hierarchy, ref_hierarchy, cutoff_radius=10, require_non_protein_identical=True):
    """Perform a flexible alignment on the two hierarchies. Will return a list of alignment metrics, one for each residue."""

    try:
        mov_hierarchy.only_model()
        ref_hierarchy.only_model()
    except:
        raise Exception('structures for alignment can only have one model!')

    # Filter the structures to only the protein (don't care about waters)
    if not require_non_protein_identical:
        # Extract caches
        ref_cache = ref_hierarchy.atom_selection_cache()
        mov_cache = mov_hierarchy.atom_selection_cache()
        # Strip structures to protein only
        ref_hierarchy = ref_hierarchy.select(ref_cache.iselection('pepnames'))
        mov_hierarchy = mov_hierarchy.select(mov_cache.iselection('pepnames'))

    # Check the structures are identical
    assert mov_hierarchy.is_similar_hierarchy(ref_hierarchy), 'HIERARCHIES MUST BE SIMILAR'

    # Extract caches
    ref_cache = ref_hierarchy.atom_selection_cache()
    mov_cache = mov_hierarchy.atom_selection_cache()

    # Caches for the reference structure
    ref_backbone_hier  = ref_hierarchy.select(ref_cache.iselection('pepnames and (name CA or name C or name O or name N)'))
    ref_backbone_cache = ref_backbone_hier.atom_selection_cache()
    # Caches for the query structure
    mov_backbone_hier  = mov_hierarchy.select(mov_cache.iselection('pepnames and (name CA or name C or name O or name N)'))
    mov_backbone_cache = mov_backbone_hier.atom_selection_cache()

    # Tranformation dictionaries
    rts = {}

    # Calculate a transformation matrix for each calpha
    for calpha_at in ref_hierarchy.select(ref_cache.iselection('pepnames and name CA')).atoms_with_labels():
        # Get selection within cutoff distance
        align_sel = nearby_coords_bool( query  = calpha_at.xyz,
                                        coords = ref_backbone_hier.atoms().extract_xyz(),
                                        cutoff = cutoff_radius  )
        # Extract full atom groups where one of the atoms has been selected
        align_ags = ref_backbone_hier.select(align_sel).atom_groups()
        align_str = ' or '.join(['(chain "{}" and resid "{}" and altloc "{}")'.format(ag.parent().parent().id, ag.parent().resid(), ag.altloc if ag.altloc else ' ') for ag in align_ags])
        # Extract atoms to be aligned
        ref_align_ats = ref_backbone_hier.select(ref_backbone_cache.selection(align_str)).atoms()
        mov_align_ats = mov_backbone_hier.select(mov_backbone_cache.selection(align_str)).atoms()
        # Check that atoms are identical in the two hierarchies
        assert len(ref_align_ats) == len(mov_align_ats), 'atoms to be aligned are not the same length'
        assert len(ref_align_ats) > 0, len(ref_align_ats)
        assert len(ref_align_ats) >= sum(align_sel), 'Aligned atoms fewer than select atoms: {} < {}'.format(len(ref_align_ats), sum(align_sel))
        for ref_at, mov_at in zip(ref_align_ats, mov_align_ats):
            assert ref_at.pdb_label_columns() == mov_at.pdb_label_columns(), 'ATOMS MUST BE IDENTICAL'
        # Calculate the alignment for this residue
        rt = superpose.least_squares_fit(   reference_sites = ref_align_ats.extract_xyz(),
                                            other_sites     = mov_align_ats.extract_xyz()   ).rt()
        # Save the rotation matrix
        rts[make_label(calpha_at)] = rt
    return rts

def find_nearest_calphas(hierarchy, coordinates):
    """Find the nearest calpha in hierarchy for each coordinate in coordinates"""
    import scipy.spatial
    # Extract calphas and associated labels
    hierarchy = hierarchy.select(hierarchy.atom_selection_cache().selection('pepnames and name CA'))
    calpha_sites, calpha_labs = zip(*[(a.xyz, make_label(a)) for a in hierarchy.atoms_with_labels()])
    # Build a tree to find the smallest distance
    tree = scipy.spatial.KDTree(data=calpha_sites)
    # Find the nearest neighbours
    nn_dists, nn_groups = tree.query(coordinates)
    # Extract the associated labels
    mappings = [calpha_labs[i] for i in nn_groups]
    return mappings

def transform_coordinates_with_flexible_alignment(alignments, coordinates, mappings, inverse=False):
    """Transform coordinates by associated alignments associated with mappings values"""

    assert len(coordinates) == len(mappings)
    if not isinstance(coordinates, flex.vec3_double):
        coordinates = flex.vec3_double(coordinates)
    # Sort the indices by the mapping values
    num_tot = len(coordinates)
    sorted_idxs     = flex.size_t(sorted(range(num_tot), key=lambda i: mappings[i]))
    sorted_coords   = coordinates.select(sorted_idxs)
    sorted_mappings = [mappings[i] for i in sorted_idxs]

    # Initialise output array
    out_coords = numpy.zeros(len(coordinates), dtype=[('x',float),('y',float),('z',float)])

    # Iterate through the coords in groups and transform
    for lab, lab_idxs in itertools.groupby(range(num_tot), key=lambda i: sorted_mappings[i]):
        lab_idxs = flex.size_t(lab_idxs)
#        print 'Using RT for {}'.format(lab)
#        print 'on {} points'.format(len(lab_idxs))
#        print 'from idx {} to {}'.format(lab_idxs[0], lab_idxs[-1])
        # Extract coordinates for this block
        lab_coords   = sorted_coords.select(lab_idxs)
        lab_mappings = [sorted_mappings[i] for i in lab_idxs]
        orig_idxs    = sorted_idxs.select(lab_idxs)
        assert max(lab_idxs)-min(lab_idxs) == len(lab_idxs)-1
        assert len(set(lab_mappings)) == 1
        # Extract RT matrix
        rt = alignments[lab]
        if inverse: rt = rt.inverse()
        # Transform and put back
        rt_coords = rt * lab_coords
        out_coords.put(orig_idxs, rt_coords)

    return flex.vec3_double(out_coords)

def extract_sites_for_alignment(chain_obj):
    """Extract sequence and sites of c-alphas - adapted from mmtbx.command_line.super"""

    seq = []
    sites = flex.vec3_double()
    use_sites = flex.bool()
    for resi in chain_obj.conformers()[0].residues():
        if (   iotbx.pdb.common_residue_names_get_class(name=resi.resname) != "common_amino_acid"):
            continue
        resn = resi.resname
        single = iotbx.pdb.amino_acid_codes.one_letter_given_three_letter[resn]
        seq.append(single)
        use = False
        xyz = (0,0,0)
        for atom in resi.atoms():
            if (atom.name == " CA "):
              xyz = atom.xyz
              use = True
              break
        sites.append(xyz)
        use_sites.append(use)
    return "".join(seq), sites, use_sites

def align_hierarchies(mov_hier, ref_hier):
    """Extract c-alpha sites from the structures and align"""
    return align_chains(mov_chain=protein(mov_hier, copy=True).models()[0].only_chain(),
                        ref_chain=protein(ref_hier, copy=True).models()[0].only_chain() )

def align_chains(mov_chain, ref_chain):
    """Takes two chains and aligns them - return rt_mx"""

    mov_seq, mov_sites, mov_flags = extract_sites_for_alignment(mov_chain)
    ref_seq, ref_sites, ref_flags = extract_sites_for_alignment(ref_chain)

    align_obj = mmtbx.alignment.align(
        seq_a=ref_seq,
        seq_b=mov_seq,
        gap_opening_penalty = 20,
        gap_extension_penalty = 2,
        similarity_function = 'blosum50',
        style = 'local')

    # Extract the alignment
    alignment = align_obj.extract_alignment()
    # List of matches - '|' for exact match, '*' for good match
    matches = alignment.matches()
    equal = matches.count("|")
    similar = matches.count("*")
    total = len(alignment.a) - alignment.a.count("-")
    alignment.pretty_print(
        matches=matches,
        block_size=50,
        n_block=1,
        top_name="fixed",
        bottom_name="moving")

    # Create list of selected sites
    ref_sites_sel = flex.vec3_double()
    mov_sites_sel = flex.vec3_double()
    for ia,ib,m in zip(alignment.i_seqs_a, alignment.i_seqs_b, matches):
        if (m not in ["|", "*"]): continue
        # Check that the sites are flagged to be used
        if (ref_flags[ia] and mov_flags[ib]):
            # Append sites to list to align
            ref_sites_sel.append(ref_sites[ia])
            mov_sites_sel.append(mov_sites[ib])

    if (ref_sites_sel.size() == 0):
      raise Exception("No matching C-alpha atoms.")

    # Create LSQ
    lsq_fit = superpose.least_squares_fit(
        reference_sites = ref_sites_sel,
        other_sites     = mov_sites_sel)

    return lsq_fit

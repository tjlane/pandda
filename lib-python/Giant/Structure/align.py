
from scitbx.math import superpose
from scitbx.array_family import flex

import iotbx.pdb
import mmtbx.alignment

def perform_flexible_alignment(mov_hierarchy, ref_hierarchy, cutoff_radius=10):
    """Perform a flexible alignment on the two hierarchies. Will return a list of alignment metrics, one for each residue."""

    assert mov_hierarchy.is_similar_hierarchy(ref_hierarchy), 'HIERARCHIES MUST BE SIMILAR'

    ref_cache = ref_hierarchy.atom_selection_cache()
    calpha_idxs = ref_cache.iselection('pepnames and name CA')
    backbone_idxs = ref_cache.iselection('pepnames and (name CA or name C or name O or name N)')

    # Extract Xray Structures for the reference structure
    ref_xray = ref_hierarchy.extract_xray_structure()

    calpha_labels = [at.id_str() for at in ref_hierarchy.select(calpha_idxs).atoms()]

    # Tranformation dictionaries
    output_fits = {}

    # Iterate through the calphas
    for calpha_at in ref_hierarchy.select(calpha_idxs).atoms_with_labels():
        # Pull out the backbone atoms for this residue
        res_idxs = ref_cache.sel_resid(calpha_at.resid()).iselection()
        res_backbone_idxs = res_idxs.intersection(backbone_idxs)
        res_backbone_ats = ref_hierarchy.select(res_backbone_idxs).atoms()

        # Select atoms within cutoff_radius of the backbone of this residue
        nearby_scats = ref_xray.select(ref_xray.selection_within(cutoff_radius, ref_xray.label_selection(*[at.id_str() for at in res_backbone_ats]))).scatterers()
        # Convert back to hierarchy and select the backbone atoms that have been found
        nearby_labels = [s.label for s in nearby_scats]
        nearby_ats = [at for at in ref_hierarchy.select(backbone_idxs).atoms_with_labels() if at.id_str() in nearby_labels]
        # Build a list of residue ids where a backbone atom is within range
        nearby_resids = sorted(list(set([at.resid() for at in nearby_ats])))

        # Build a list of atom indices where a backbone atom is within range - for the reference structure
        ref_selected_idxs = flex.size_t(sorted([idx for rid in nearby_resids for idx in ref_cache.sel_resid(rid).iselection()]))
        # Get the atoms and sites for these atoms
        ref_selected_ats = ref_hierarchy.select(ref_selected_idxs).atoms()

        # Build a list of atom indices where a backbone atom is within range - for the moving structure
        mov_selected_idxs = flex.size_t(sorted([idx for rid in nearby_resids for idx in ref_cache.sel_resid(rid).iselection()]))
        # Get the atoms and sites for these for the
        mov_selected_ats = mov_hierarchy.select(mov_selected_idxs).atoms()

        # Check that
        for ref_at, mov_at in zip(ref_selected_ats, mov_selected_ats):
            assert ref_at.pdb_label_columns() == mov_at.pdb_label_columns(), 'ATOMS MUST BE IDENTICAL'

        # Extract the actual coordinates of the sites
        ref_sites = ref_selected_ats.extract_xyz()
        mov_sites = mov_selected_ats.extract_xyz()

        # Calculate the alignment for this residue
        res_rt = superpose.least_squares_fit(reference_sites=ref_sites, other_sites=mov_sites).rt()

        output_fits[(calpha_at.chain_id, calpha_at.resid())] = res_rt

    return output_fits

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

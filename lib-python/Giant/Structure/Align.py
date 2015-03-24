
from scitbx.math import superpose
from scitbx.array_family import flex

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
        ref_selected_idxs = flex.size_t([idx for rid in nearby_resids for idx in ref_cache.sel_resid(rid).iselection()])
        # Get the atoms and sites for these for the
        ref_selected_ats = ref_hierarchy.select(ref_selected_idxs).atoms()

        # Build a list of atom indices where a backbone atom is within range - for the moving structure
        mov_selected_idxs = flex.size_t([idx for rid in nearby_resids for idx in ref_cache.sel_resid(rid).iselection()])
        # Get the atoms and sites for these for the
        mov_selected_ats = mov_hierarchy.select(mov_selected_idxs).atoms()

        # Check that
        for ref_at, mov_at in zip(ref_selected_ats, mov_selected_ats):
            assert ref_at.pdb_label_columns() == mov_at.pdb_label_columns(), 'ATOMS MUST BE IDENTICAL'

        # Extract the actual coordinates of the sites
        ref_sites = ref_selected_ats.extract_xyz()
        mov_sites = mov_selected_ats.extract_xyz()

        # Calculate the alignment for this residue
        res_rt = superpose.least_squares_fit(reference_sites=ref_sites, other_sites=mov_sites)

        output_fits[(calpha_at.chain_id, calpha_at.resid())] = res_rt

    return output_fits

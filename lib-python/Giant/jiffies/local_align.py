import os, sys, copy, re

import libtbx.phil

import numpy

import iotbx.pdb
import scipy.spatial

from scitbx.array_family import flex
from Giant.Structure.align import perform_flexible_alignment

blank_arg_prepend = 'pdb='

master_phil = libtbx.phil.parse("""
input {
    reference_pdb = None
        .type = path
    pdb = None
        .type = path
        .multiple = True
}
alignment {
    cutoff_radius = 10.0
        .type = float
}
output {
    prefix = aligned_
        .type = path
}
""")

def run(params):

    assert params.input.reference_pdb is not None, 'No reference_pdb provided'

    # =================================================>
    # Load structure
    # =================================================>
    ref_hierarchy = iotbx.pdb.hierarchy.input(params.input.reference_pdb).hierarchy
    ref_hierarchy_trimmed = ref_hierarchy.select(ref_hierarchy.atom_selection_cache().selection('pepnames and (name CA or name C or name O or name N)'))

    for mov_pdb in params.input.pdb:

        # =================================================>
        # Load structure
        # =================================================>
        mov_hierarchy = iotbx.pdb.hierarchy.input(mov_pdb).hierarchy
        mov_hierarchy_trimmed = mov_hierarchy.select(mov_hierarchy.atom_selection_cache().selection('pepnames and (name CA or name C or name O or name N)'))

        # =================================================>
        # Calculate alignments
        # =================================================>
        alignment_mxs = perform_flexible_alignment( mov_hierarchy = ref_hierarchy_trimmed,
                                                    ref_hierarchy = mov_hierarchy_trimmed,
                                                    cutoff_radius = params.alignment.cutoff_radius    )

        # =================================================>
        # Build clustering tree
        # =================================================>
        mov_calpha_hierarchy = mov_hierarchy.select(mov_hierarchy.atom_selection_cache().selection('pepnames and name CA'))
        mov_atom_sites, mov_atom_labels = zip(*[(a.xyz, (a.chain_id, a.resid())) for a in mov_calpha_hierarchy.atoms_with_labels()])
        mov_tree = scipy.spatial.KDTree(data = mov_atom_sites)

        # =================================================>
        # Find calpha mappings to atom sites
        # =================================================>
        nn_dists, nn_groups = mov_tree.query(mov_hierarchy.atoms().extract_xyz())
        point_mappings = [mov_atom_labels[i] for i in nn_groups]

        # =================================================>
        # Transform coordinates
        # =================================================>
        mov_coords = mov_hierarchy.atoms().extract_xyz()
        # Get the set of labels to transform in groups
        lab_set = sorted(list(set(point_mappings)))
        # Initialise output array
        rt_coords = numpy.zeros(len(mov_coords), dtype=[('x',float),('y',float),('z',float)])
        for r_lab in lab_set:
            # Extract the idxs and points with this label in the mapping
            lab_idxs, lab_coords = zip(*[(i, mov_coords[i]) for i, i_lab in enumerate(point_mappings) if r_lab==i_lab])
            # Transform all of these points at once
            lab_rt_coords = alignment_mxs[r_lab].inverse() * flex.vec3_double(lab_coords)
            # Populate the array at the appropriate place
            rt_coords.put(lab_idxs, lab_rt_coords)
        aligned_coords = flex.vec3_double(rt_coords)

        # =================================================>
        # Output structure
        # =================================================>
        # Copy the input structure
        new_hierarchy = mov_hierarchy.deep_copy()
        new_hierarchy.atoms().set_xyz(aligned_coords)
        new_hierarchy.write_pdb_file(file_name=os.path.join(os.path.dirname(mov_pdb), params.output.prefix+os.path.basename(mov_pdb)))


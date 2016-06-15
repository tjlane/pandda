import os, sys, copy, re
import time

import libtbx.phil

import numpy

import iotbx.pdb
import scipy.spatial

from scitbx.array_family import flex
from giant.structure.align import perform_flexible_alignment, find_nearest_calphas, transform_coordinates_with_flexible_alignment

#######################################

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

#######################################

def run(params):

    assert params.input.reference_pdb is not None, 'No reference_pdb provided'

    # =================================================>
    # Load structure
    # =================================================>
    ref_hierarchy = iotbx.pdb.hierarchy.input(params.input.reference_pdb).hierarchy

    t_start = time.time()

    for mov_pdb in params.input.pdb:

        print 'ALIGNING: {}'.format(mov_pdb)

        # =================================================>
        # Load structure
        # =================================================>
        mov_hierarchy = iotbx.pdb.hierarchy.input(mov_pdb).hierarchy
#        mov_hierarchy_trimmed = mov_hierarchy.select(mov_hierarchy.atom_selection_cache().selection('pepnames and (name CA or name C or name O or name N)'))

        # =================================================>
        # Calculate alignments
        # =================================================>
        alignments = perform_flexible_alignment( mov_hierarchy = ref_hierarchy,
                                                 ref_hierarchy = mov_hierarchy,
                                                 cutoff_radius = params.alignment.cutoff_radius    )

        # =================================================>
        # Find which alignments should be used for coords
        # =================================================>
        mappings = find_nearest_calphas( hierarchy   = mov_hierarchy,
                                         coordinates = mov_hierarchy.atoms().extract_xyz()   )

        # =================================================>
        # Transform coordinates
        # =================================================>
        aligned_coords = transform_coordinates_with_flexible_alignment( alignments  = alignments,
                                                                        coordinates = mov_hierarchy.atoms().extract_xyz(),
                                                                        mappings    = mappings,
                                                                        inverse     = True   )

        # =================================================>
        # Output structure
        # =================================================>
        # Copy the input structure
        new_hierarchy = mov_hierarchy.deep_copy()
        new_hierarchy.atoms().set_xyz(aligned_coords)
        new_hierarchy.write_pdb_file(file_name=os.path.join(os.path.dirname(mov_pdb), params.output.prefix+os.path.basename(mov_pdb)))

    t_end = time.time()
    print 'Total Runtime: {!s}'.format(time.strftime("%H hours:%M minutes:%S seconds", time.gmtime(t_end-t_start)))

#######################################

if __name__ == '__main__':
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)

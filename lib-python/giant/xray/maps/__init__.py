import numpy

def scale_map_to_reference(ref_map, map, mask_idxs):
    """Scale given map to reference over the point_mask"""

    # Select values for the masked points
    ref_vals = ref_map.select(mask_idxs)
    map_vals =     map.select(mask_idxs)

    # Fit ref = a*map + b
    # a, b
    scale, offset = numpy.polyfit(x=map_vals, y=ref_vals, deg=1)
#    print '\nScale: {}\tOffset: {}'.format(scale, offset)

    return (map*scale) + offset

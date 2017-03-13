import numpy

import iotbx.ccp4_map
import cctbx.maptbx

from scitbx.array_family import flex

from pandda.analyse import graphs as analyse_graphs

def write_array_as_map(grid, array, f_name):
    """Take array on the reference grid and write to map"""
    iotbx.ccp4_map.write_ccp4_map(  file_name   = f_name,
                                    unit_cell   = grid.unit_cell(),
                                    space_group = grid.space_group(),
                                    map_data    = array,
                                    labels      = flex.std_string(['Map from pandda'])     )

def write_bool_as_map(grid, array, f_name):
    """Write boolean array as electron density map"""
    # Create binary map data
    array = array.astype(int)
    # Convert to flex and reshape
    array = flex.double(array.tolist())
    array.reshape(grid.indexer())
    # Write as array
    write_array_as_map(grid=grid, array=array, f_name=f_name)

def write_indices_as_map(grid, indices, f_name):
    """Take list of indices and write as a mask"""
    # Create bool map data
    array = numpy.zeros(grid.grid_size_1d(), dtype=bool)
    array.put(indices, True)
    # Write as array
    write_bool_as_map(grid=grid, array=array, f_name=f_name)

def rotate_map(grid, d_handler, map_data, align_on_grid_point=None):
    """Apply an RT matrix to an array on the reference grid"""

    if (align_on_grid_point is not None):
        # For the local alignment transformation
        rt_lab = grid.partition().query_by_grid_points([align_on_grid_point])[0]
        print('=> Aligning Event Map to: Chain {}, Residue {}'.format(rt_lab[0], rt_lab[1].strip()))
        rt = d_handler.local_alignment_transforms()[rt_lab]
    else:
        # For the global alignment transformation
        rt = d_handler.global_alignment_transform()

    return cctbx.maptbx.rotate_translate_map(   unit_cell          = grid.unit_cell(),
                                                map_data           = map_data,
                                                rotation_matrix    = rt.r.elems,
                                                translation_vector = rt.t.elems    )


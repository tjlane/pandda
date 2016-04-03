import iotbx.ccp4_map
import cctbx.maptbx

from scitbx.array_family import flex

from pandda import analyse_graphs

def write_map_value_distribution(map_vals, output_file, plot_indices=None, plot_normal=False):
    """Write out the value distribution for a map"""
    try:
        if plot_indices: plot_vals = [map_vals[i] for i in plot_indices]
        else:            plot_vals = list(map_vals)
        analyse_graphs.map_value_distribution(f_name=output_file, plot_vals=plot_vals, plot_normal=plot_normal)
    except: pass

def write_qq_plot_against_normal(map_vals, output_file, plot_indices=None):
    """Plot the values in map_vals against those expected from a normal distribution"""
    try:
        if plot_indices: plot_vals = [map_vals[i] for i in plot_indices]
        else:            plot_vals = list(map_vals)
        analyse_graphs.qq_plot_against_normal(f_name=output_file, plot_vals=plot_vals)
    except: pass

def write_array_to_map(output_file, map_data, grid):
    """Take array on the reference grid and write to map"""
    iotbx.ccp4_map.write_ccp4_map(  file_name   = output_file,
                                    unit_cell   = grid.unit_cell(),
                                    space_group = grid.space_group(),
                                    map_data    = map_data,
                                    labels      = flex.std_string(['Map from pandda'])     )

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


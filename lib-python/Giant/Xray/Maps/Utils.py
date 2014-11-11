import mmtbx.f_model as model_handler
import iotbx.map_tools as map_tools

from cctbx import maptbx
from libtbx.math_utils import ifloor, iceil
from iotbx import crystal_symmetry_from_any
from scitbx.array_family import flex

def get_fft_map_from_f_obs_and_structure(f_obs, xray_structure, map_type='Fobs', d_min=None):
    """Create map from structure factors from array, phases from structure"""

    fmodel = model_handler.manager(f_obs=f_obs, xray_structure = xray_structure)
    density_map = fmodel.electron_density_map(update_f_part1=False)
    return density_map.map_coefficients(map_type).fft_map(d_min=d_min, symmetry_flags=maptbx.use_space_group_symmetry)

def write_1d_array_as_p1_map(file_name, map_data, grid_size, grid_spacing):
    """Take map data and write out to ccp4 map file - assumes the grid starts from the origin and the grid is cartesian!"""

    # Calculate the dimensions of the p1 box needed to enclose the map
    if isinstance(grid_spacing, tuple) or isinstance(grid_spacing, list):
        # (d1,d2,d3) -> non-isotropic grid spacing
        assert len(grid_spacing)==3, 'Grid spacing must be either a number or a tuple of length 3'
        dims_of_p1_box = tuple(flex.double(grid_size)*flex.double(grid_spacing))
    else:
        # d -> isotropic grid spacing
        dims_of_p1_box = tuple(flex.double(grid_size)*grid_spacing)
#    print 'Size of P1 Box:', dims_of_p1_box
    p1_box = generate_p1_box_from_size(dims_of_p1_box)
    # Copy and reshape
    map_data_to_write = map_data.deep_copy()
    map_data_to_write.reshape(flex.grid(grid_size))
    # Write map
    map_tools.write_ccp4_map(sites_cart=None, unit_cell=p1_box.unit_cell(), map_data=map_data_to_write, n_real=grid_size, file_name=file_name, buffer=0)

def generate_p1_box_from_sites(sites_cart, buffer=10.0) :
    xyz_max = sites_cart.max()
    a = xyz_max[0] + buffer
    b = xyz_max[1] + buffer
    c = xyz_max[2] + buffer
    return generate_p1_box_from_size((a,b,c))

def generate_p1_box_from_size(max_size):
    a,b,c = max_size
    combined = "%.3f,%.3f,%.3f,90,90,90,P1" % (a, b, c)
    symm = crystal_symmetry_from_any.from_string(combined)
    return symm


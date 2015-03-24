import os, sys, glob

import iotbx.pdb as pdb_reader
import iotbx.map_tools as map_tools

import numpy
from scitbx.math import superpose
from scitbx.array_family import flex

from cctbx import maptbx
from iotbx.reflection_file_utils import extract_miller_array_from_file
from libtbx.math_utils import ifloor, iceil

from Giant.Xray.Miller.Utils import apply_symmetry_to_miller_array, check_miller_arrays_compatible_for_scaling
from Giant.Xray.Miller.Utils import scale_miller_array_to_reference, scale_amplitude_arrays_via_intensities
from Giant.Xray.Maps.Utils import get_fft_map_from_f_obs_and_structure, generate_p1_box_from_size
from Giant.Grid.Utils import calculate_sampling_distance, get_bounding_box_for_structure, calculate_grid_size, create_cartesian_grid
from Giant.Xray.Structure.Select import get_calpha_sites
from Giant.Stats.Moments import skew, kurtosis
from Giant.Stats.Normalise import normalise_array_to_z_scores

# ============================================================================>

#####
# Read Reference Data
#####

# ============================================================================>

input_dir = '/home/npearce/Analyses/FragmentSoaks/BAZ2BA-Feb13-Soak'

ref_pdb = os.path.join(input_dir,'reference.pdb')
ref_mtz = os.path.join(input_dir,'reference.mtz')
assert os.path.exists(ref_pdb)
assert os.path.exists(ref_mtz)

ref_input = pdb_reader.input(source_info=None,lines=open(ref_pdb,'r').read())
ref_struc = ref_input.xray_structure_simple()

ref_calpha_sites = get_calpha_sites(input_obj=ref_input, structure_obj=ref_struc)
print 'Num C-Alpha Sites: ', ref_calpha_sites.size()

# Get the observed data for the reference structure
ref_miller = extract_miller_array_from_file(ref_mtz, 'F,SIGF', log=open(os.devnull,'w'))
ref_miller.show_summary()

# Get the map data for the reference structure
ref_map = get_fft_map_from_f_obs_and_structure(ref_miller, ref_struc)
ref_map.apply_sigma_scaling()
ref_unit_cell = ref_map.unit_cell()
ref_space_group = ref_map.space_group()

# Create map handler for the reference structure
ref_basic_map = maptbx.basic_map(
        maptbx.basic_map_unit_cell_flag(),
        ref_map.real_map(),
        ref_map.real_map().focus(),
        ref_unit_cell.orthogonalization_matrix(),
        maptbx.out_of_bounds_clamp(0).as_handle(),
        ref_unit_cell)

print 'Size of Map: ', ref_map.n_real(), ref_map.n_grid_points()

print '===================================>>>'

# ============================================================================>

#####
# Create and map sites to the unit cell
#####

# ============================================================================>

# Set the target resolution for the other datasets # XXX Change this to be dynamic or otherwise defined in the next iterations of the script XXX
ref_res = ref_miller.d_min()
print 'Reference Resolution', ref_res

ref_struc_extent = get_bounding_box_for_structure(ref_struc)
print 'Min/Max Sites', ref_struc_extent

# Calculate the sampling distance between points
res_factor = 0.25
grid_spacing = calculate_sampling_distance(ref_res, res_factor)

# Create a grid (cartesian) of points for sampling
box_cart_size, grid_size, sample_points_cart_ref = create_cartesian_grid(min_carts=(0,0,0), max_carts=ref_struc_extent[1], grid_spacing=grid_spacing)
print 'Size of Ref Structure Bounding Box', box_cart_size
print 'Size of Ref Sample Grid', grid_size

# Calculate the distance between the grid points
grid_point_volume = grid_spacing**3
print 'Grid Point Volume', grid_point_volume

# Pull the reference structure map values at the points
sample_map_vals_ref = ref_basic_map.get_cart_values(sample_points_cart_ref)

# Number of grid points
num_grid_points = grid_size[0]*grid_size[1]*grid_size[2]

assert num_grid_points == len(sample_map_vals_ref)

print 'First Point: ', sample_points_cart_ref[0]
print 'Last Point: ', sample_points_cart_ref[-1]

# ============================================================================>

#####
# Save map of the reference unit cell in P1
#####

# ============================================================================>

# P1 BOX THE SIZE OF THE UNIT CELL
dims_of_p1_box = ref_unit_cell.orthogonalize((1,1,1))
print 'Size of P1 Box:', dims_of_p1_box
p1_box = generate_p1_box_from_size(dims_of_p1_box)

# P1 BOX THE SIZE OF THE PROTEIN
dims_of_p1_box = box_cart_size
print 'Size of P1 Box:', dims_of_p1_box
p1_box = generate_p1_box_from_size(dims_of_p1_box)

# P1 BOX FROM ORIGIN TO ACTUAL EXTENT OF GRID (using the fact that the grid is slightly larger than the box)
dims_of_p1_box = tuple(flex.double(grid_size)*grid_spacing)
print 'Size of P1 Box:', dims_of_p1_box
p1_box = generate_p1_box_from_size(dims_of_p1_box)

sample_map_vals_ref_copy = sample_map_vals_ref.deep_copy()
sample_map_vals_ref_copy.reshape(flex.grid(grid_size))

map_tools.write_ccp4_map(sites_cart=sample_points_cart_ref, unit_cell=p1_box.unit_cell(), map_data=sample_map_vals_ref, n_real=grid_size, file_name=ref_mtz.replace('.mtz','.sampled.ccp4'), buffer=0)

#import iotbx.ccp4_map
#iotbx.ccp4_map.write_ccp4_map(unit_cell=p1_box.unit_cell(), map_data=sample_map_vals_ref, file_name=ref_mtz.replace('.mtz','.sampled.ccp4'), buffer=0)













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
from Giant.Xray.Maps.Grid import calculate_sampling_distance, get_bounding_box_for_structure, calculate_grid_size, create_cartesian_grid
from Giant.Xray.Structure.Select import get_calpha_sites
from Giant.Stats.Moments import skew, kurtosis
from Giant.Stats.Normalise import normalise_array_to_z_scores

# ============================================================================>

#####
# Output Files
#####

# ============================================================================>

output_dir = '/work/SP100A/25-FS-Maybridge-Poised-lib'

output_points = os.path.join(output_dir, 'sample_points.csv')
output_map_point_summary = os.path.join(output_dir, 'map_point_summary.csv')
output_z_values = os.path.join(output_dir, 'z_values.csv')
output_z_values_summary = os.path.join(output_dir, 'z_values_summary.csv')
if os.path.exists(output_points):
    os.remove(output_points)
if os.path.exists(output_map_point_summary):
    os.remove(output_map_point_summary)
if os.path.exists(output_z_values):
    os.remove(output_z_values)
if os.path.exists(output_z_values_summary):
    os.remove(output_z_values_summary)

# ============================================================================>

#####
# Global Variables
#####

# ============================================================================>

map_type = ['Fobs','2mFo-DFc'][1]

# ============================================================================>

#####
# Read Reference Data
#####

# ============================================================================>

input_dir = '/work/SP100A/25-FS-Maybridge-Poised-lib'

ref_pdb = os.path.join(input_dir,'Reference','reference.pdb')
ref_mtz = os.path.join(input_dir,'Reference','reference.mtz')
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
ref_map = get_fft_map_from_f_obs_and_structure(f_obs=ref_miller, xray_structure=ref_struc, map_type=map_type)
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

include_origin_in_grid = True

# Create a grid (cartesian) of points for sampling
if include_origin_in_grid:
    box_cart_size, grid_size, sample_points_cart_ref = create_cartesian_grid(min_carts=(0,0,0), max_carts=ref_struc_extent[1], grid_spacing=grid_spacing)
else:
    box_cart_size, grid_size, sample_points_cart_ref = create_cartesian_grid(min_carts=ref_struc_extent[0], max_carts=ref_struc_extent[1], grid_spacing=grid_spacing)
print 'Size of Sampling Box', box_cart_size
print 'Size of Sample Grid', grid_size

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

# Write the output points
with open(output_points, 'w') as points_file:
    print 'Writing Grid Sample Points File...'
    # Write (x,y,z) tuples on each line
    for i in xrange(num_grid_points):
        if i%1000==0: print '\r>>', round(100.0*(i+1)/num_grid_points,0), '%',; sys.stdout.flush()
        output_line = [str(round(x,3)) for x in sample_points_cart_ref[i]]
        points_file.write(', '.join(output_line)+'\n')
    print ''

# Write output map (just for fun)
print 'Writing Maps!'
dims_of_ref_map_p1_box = tuple(flex.double(grid_size)*grid_spacing)
print 'Size of P1 Box:', dims_of_ref_map_p1_box
ref_map_p1_box = generate_p1_box_from_size(dims_of_ref_map_p1_box)
# Copy and reshape
map_data_to_write = sample_map_vals_ref.deep_copy()
map_data_to_write.reshape(flex.grid(grid_size))
map_name_to_write = ref_mtz.replace('.mtz','.sampled.ccp4')
# Write map
map_tools.write_ccp4_map(sites_cart=sample_points_cart_ref, unit_cell=ref_map_p1_box.unit_cell(), map_data=map_data_to_write, n_real=grid_size, file_name=map_name_to_write, buffer=0)
print 'done!'

# ============================================================================>

#####
# Create list of data
#####

# ============================================================================>

raw_pdbs = sorted(glob.glob(os.path.join(input_dir,'Processing/SP100A-*/refine.pdb')))
pdbs = []
mtzs = []

raw_file_pairs = []

for pdb in raw_pdbs:
    mtz = pdb.replace('.pdb','.mtz')
    if os.path.exists(mtz):
        raw_file_pairs.append((pdb,mtz))

if not raw_file_pairs: raise SystemExit('No Files Found!')

output_transformations = []
output_inv_transformations = []

output_maps = []
output_z_maps = []

num_raw_datasets = len(raw_file_pairs)
file_pairs = []

print 'Number of Datasets: ', num_raw_datasets

# ============================================================================>

#####
# Scale and Process All Data
#####

# ============================================================================>

for d_num, (pdb, mtz) in enumerate(raw_file_pairs):

    # Separator
    print '===================================>>>'

    print 'Dataset', d_num+1

    print pdb
    print mtz

    # Read New Data
    new_miller = extract_miller_array_from_file(mtz, 'F-obs,SIGF-obs', log=open(os.devnull,'w'))

    # Scale New Data
    new_miller_sym = apply_symmetry_to_miller_array(chg_miller=new_miller, ref_miller=ref_miller)
    new_miller_sym_scale = scale_amplitude_arrays_via_intensities(chg_arr_a=new_miller_sym, ref_arr_a=ref_miller, n_bins=20)
    new_miller_scale = apply_symmetry_to_miller_array(chg_miller=new_miller_sym_scale, ref_miller=new_miller)
#    new_miller_scale.show_summary()

    new_ref_correlation = new_miller_sym_scale.correlation(ref_miller).coefficient()
    print 'Scaled-Ref Correlation: ', new_ref_correlation
    new_new_correlation = new_miller_scale.correlation(new_miller).coefficient()
    print 'Scaled-Orig Correlation:', new_new_correlation

    if new_new_correlation < 0.9:
        print 'Removing Dataset!'
        continue

    # Read PDB File
    new_input = pdb_reader.input(source_info=None,lines=open(pdb,'r').read())
    new_struc = new_input.xray_structure_simple()

    # Convert To Map
    new_map = get_fft_map_from_f_obs_and_structure(f_obs=new_miller_scale, xray_structure=new_struc, map_type=map_type, d_min=ref_res)
    new_map.apply_sigma_scaling()
    new_unit_cell = new_map.unit_cell()
    new_space_group = new_map.space_group()

    print 'Map Resolution: ', new_map.d_min()
    print 'Size of Map: ', new_map.n_real(), new_map.n_grid_points()

    # Align Maps
    new_calpha_sites = get_calpha_sites(input_obj=new_input, structure_obj=new_struc)
    print 'Num C-Alpha Sites: ', new_calpha_sites.size()

    try:
        # Sample at the sites in cart_sampling_sites
        alignment_fit = superpose.least_squares_fit(reference_sites=ref_calpha_sites, other_sites=new_calpha_sites)
    except RuntimeError:
        continue

#    print 'R:\t', '\n\t'.join([' '.join(map(str,l)) for l in alignment_fit.r.as_list_of_lists()])
#    print 'T:\t', '\n\t'.join([' '.join(map(str,l)) for l in alignment_fit.t.as_list_of_lists()])

    # New -> Ref mappings
    new_to_ref_rt = alignment_fit.rt()
    # Ref -> New mappings
    ref_to_new_rt = alignment_fit.rt().inverse()

    # Transform the sites to sample (ref -> new transform)
    sample_points_cart_new = ref_to_new_rt*sample_points_cart_ref

    # Check that the right matrices are being used
    print 'Initial C-alpha RMSD', new_calpha_sites.rms_difference(ref_calpha_sites)
    print 'Using ref_to_new to transform ref', new_calpha_sites.rms_difference(ref_to_new_rt*ref_calpha_sites)

    # Create map handler for the new structure
    new_basic_map = maptbx.basic_map(
        maptbx.basic_map_unit_cell_flag(),
        new_map.real_map(),
        new_map.real_map().focus(),
        new_unit_cell.orthogonalization_matrix(),
        maptbx.out_of_bounds_clamp(0).as_handle(),
        new_unit_cell)

    sample_map_vals_new = new_basic_map.get_cart_values(sample_points_cart_new)

    # Write output map (just for fun)
    print 'Writing Maps!'
    # Copy and reshape
    map_data_to_write = sample_map_vals_new.deep_copy()
    map_data_to_write.reshape(flex.grid(grid_size))
    map_name_to_write = mtz.replace('.mtz','.sampled.ccp4')
    # Write map
    map_tools.write_ccp4_map(sites_cart=sample_points_cart_new, unit_cell=ref_map_p1_box.unit_cell(), map_data=map_data_to_write, n_real=grid_size, file_name=map_name_to_write, buffer=0)
    print 'done!'

    # Save things
    output_maps.append(sample_map_vals_new)
    output_transformations.append(alignment_fit.rt())
    output_inv_transformations.append(alignment_fit.rt().inverse())
    file_pairs.append((pdb, mtz))

#    # Superpose data
#    new_map_superposed = maptbx.superpose_maps(
#        unit_cell_1        = new_map.unit_cell(),
#        unit_cell_2        = ref_map.unit_cell(),
#        map_data_1         = new_map.real_map_unpadded(),
#        n_real_2           = ref_map.n_real(),
#        rotation_matrix    = alignment_fit.r.elems,
#        translation_vector = alignment_fit.t.elems)
#
#    print 'Size of Map: ', new_map_superposed.last(), new_map_superposed.size()
#

print '===================================>>>'

# ============================================================================>

#####
# ANALYSIS
#####

# ============================================================================>

num_datasets = len(file_pairs)

# Create array of the output data
print 'Combining output maps into array...',; sys.stdout.flush()
all_data = numpy.array(output_maps)
print 'done'

# Correlations between all of the datasets
#all_data_coeffs = numpy.corrcoef(all_data)
#print 'Correlations less than 0.9:', sum(all_data_coeffs < 0.9)
#print 'Average Cross-Correlation:', ' '.join([str(round(x, 3)) for x in (all_data_coeffs.sum(0)-1)/(all_data_coeffs.shape[1]-1)])

# ============================================================================>

# Calculate the mean and standard deviation for each grid point
# 1st moment - Mean
sample_map_vals_mean = all_data.mean(axis=0)
assert len(sample_map_vals_mean) == num_grid_points
print 'Means:', numpy.array([(round(x, 3)) for x in sample_map_vals_mean][0:5]), '(First 5 of {!s})'.format(len(sample_map_vals_mean))
# 2nd moment - Std
sample_map_vals_stds = all_data.std(axis=0)
assert len(sample_map_vals_stds) == num_grid_points
print 'Stds: ', numpy.array([(round(x, 3)) for x in sample_map_vals_stds][0:5])
# 3rd moment - Skew
sample_map_vals_skew = numpy.array([skew(all_data[:,i]) for i in xrange(num_grid_points)])
assert len(sample_map_vals_skew) == num_grid_points
print 'Skew: ', numpy.array([(round(x, 3)) for x in sample_map_vals_skew][0:5])
# 4th moment - Kurtosis
sample_map_vals_kurt = numpy.array([kurtosis(all_data[:,i]) for i in xrange(num_grid_points)])
assert len(sample_map_vals_kurt) == num_grid_points
print 'Kurt: ', numpy.array([(round(x, 3)) for x in sample_map_vals_kurt][0:5])

# ============================================================================>

# Writing out the means and the other moments for the map
with open(output_map_point_summary, 'a') as map_point_summary_file:
    # Write (mean,std,skew,kurtosis) tuples on each line
    print 'Writing Map Point Summary To File...'
    for i in xrange(num_grid_points):
        if i%1000==0: print '\r>>', round(100.0*(i+1)/num_grid_points,0), '%',; sys.stdout.flush()
        output_summary_line = [round(sample_map_vals_mean[i],3), round(sample_map_vals_stds[i],3),
                                round(sample_map_vals_skew[i],3), round(sample_map_vals_kurt[i],3)]
        map_point_summary_file.write(', '.join(map(str,output_summary_line))+'\n')
    print ''

# ============================================================================>

# Convert the map values in all_data to Z-scores based on the means and stds for the grid points
for i_data, data in enumerate(all_data):

    print '===================================>>>'

    z_array = normalise_array_to_z_scores(input_array=data, element_means=sample_map_vals_mean, element_stds=sample_map_vals_stds)

    assert i+1 == num_grid_points

    pdb, mtz = file_pairs[i_data]

    z_mean = round(numpy.mean(z_array),3)
    z_std  = round(numpy.std(z_array),3)
    z_skew = round(skew(z_array),3)
    z_kurt = round(kurtosis(z_array),3)
    z_max  = round(z_array.max(),3)

    print 'Dataset: ', i_data+1
    print 'PDB:', pdb
    print 'MTZ:', mtz
    print 'Mean:   ', z_mean
    print 'StdDev: ', z_std
    print 'Skew:   ', z_skew
    print 'Kurt:   ', z_kurt
    print 'MAX:    ', z_max

    # Append data to total array
    output_z_maps.append(z_array)

    # Write output map (just for fun)
    print 'Writing Maps!'
    # Copy and reshape
    map_data_to_write = flex.double(z_array).deep_copy()
    map_data_to_write.reshape(flex.grid(grid_size))
    map_name_to_write = mtz.replace('.mtz','.zvalues.ccp4')
    # Map the reference grid to the new unit cell
    map_points_trans = output_inv_transformations[i_data]*sample_points_cart_ref
    # Write map
    map_tools.write_ccp4_map(sites_cart=map_points_trans, unit_cell=ref_map_p1_box.unit_cell(), map_data=map_data_to_write, n_real=grid_size, file_name=map_name_to_write, buffer=0)
    print 'done!'

    # Write out data summary
    with open(output_z_values_summary, 'a') as z_values_summary_file:
        print 'Writing Z-value Summary To File...',; sys.stdout.flush()
        output_summary_line = [z_max, z_mean, z_std, z_skew, z_kurt, pdb, mtz]
        z_values_summary_file.write(', '.join(map(str,output_summary_line))+'\n')
        print 'done'

#    # Write out the data
#    with open(output_z_values, 'a') as z_values_file:
#        print 'Writing Z-values To File...',
#        output_line = [round(z,3) for z in z_array] + [pdb, mtz]
#        z_values_file.write(', '.join(map(str,output_line))+'\n')
#        print 'done'

# ============================================================================>

print '===================================>>>'

# Create array for all z-values
output_z_maps = numpy.array(output_z_maps)

# Write out the z-values for each map individually relative to the other maps
with open(output_z_values, 'a') as z_values_file:

    print 'Writing ALL Z-values To File... ({!s} Datasets, {!s} Grid Points)'.format(num_datasets, num_grid_points)

    # Column Headings and Filenames
    output_line = ['Dataset {!s}'.format(i+1) for i in xrange(num_datasets)]
    z_values_file.write(', '.join(output_line)+'\n')
    output_line = [file_pairs[i][0] for i in xrange(num_datasets)]
    z_values_file.write(', '.join(output_line)+'\n')
    output_line = [file_pairs[i][1] for i in xrange(num_datasets)]
    z_values_file.write(', '.join(output_line)+'\n')

    # Z-values at each point for each of the datasets
    for i in xrange(num_grid_points):
        if i%1000==0: print '\r>>', round(100.0*(i+1)/num_grid_points,0), '%',; sys.stdout.flush()
        output_line = map(str,output_z_maps[:,i].round(3))
        z_values_file.write(', '.join(output_line)+'\n')
    print ''

# ============================================================================>


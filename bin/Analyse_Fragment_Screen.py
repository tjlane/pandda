import os, sys, glob, time
import numpy

from PANDDAs.Main import multi_dataset_analyser, spherical_mask, protein_mask

# ============================================================================>
#####
# Reference Dataset
#####
# ============================================================================>

input_dir = '/home/npearce/Analyses/FragmentSoaks/BAZ2BA-Feb13-Soak'
ref_pdb = os.path.join(input_dir,'reference.pdb')
ref_mtz = os.path.join(input_dir,'reference.mtz')

# ============================================================================>
#####
# Create list of data
#####
# ============================================================================>

raw_pdbs = sorted(glob.glob(os.path.join(input_dir,'ProcessedFragmentSoak/BAZ2BA-*/1-apo/apo-BAZ2BA-*-refmac.pdb')))
raw_file_pairs = []

for pdb in raw_pdbs:
    mtz = pdb.replace('.pdb','.mtz')
    if os.path.exists(mtz):
        raw_file_pairs.append((pdb,mtz))

#raw_file_pairs = raw_file_pairs[0:5]

if not raw_file_pairs: raise SystemExit('No Files Found!')

num_raw_datasets = len(raw_file_pairs)
file_pairs = []

print 'Number of Datasets: ', num_raw_datasets

# ============================================================================>
# ============================================================================>
# TODO Change the name of the main processor object

pandda = multi_dataset_analyser(outdir=input_dir)
pandda.set_map_type(map_type='2mFo-DFc')

# ============================================================================>
#####
# Set Reference Dataset
#####
# Select the reference dataset
# TODO Choose the reference dataset better
# ============================================================================>

pandda.load_reference_dataset(ref_pdb=ref_pdb, ref_mtz=ref_mtz)

# ============================================================================>
#####
# Create Sample Grid
#####
# Choose the resolution limit for the maps
# Create reference grid based on the reference structure
#
# TODO Choose the resolution limit better
# TODO Add a buffer to the grid
# ============================================================================>

pandda.set_cut_resolution(d_min=pandda.get_reference_dataset().get_dataset_resolution())

if pandda.get_reference_grid() is None:
    pandda.create_reference_grid(res_factor=0.25, include_origin=True, buffer=None)

pandda.extract_reference_map_values()

# ============================================================================>
#####
# Create Local and Global Masks
#####
# Local mask used for forming groups of points around a grid point
# Global mask used for removing points in the bulk solvent regions
# ============================================================================>

if pandda.get_reference_grid().get_local_mask() is None:
    print('===================================>>>')
    print('Generating Local Mask')
    local_mask = spherical_mask(grid_spacing=pandda.get_reference_grid().get_grid_spacing(), distance_cutoff=3, grid_jump=None)
    pandda.get_reference_grid().set_local_mask(local_mask)

if pandda.get_reference_grid().get_global_mask() is None:
    print('===================================>>>')
    print('Generating Global Mask')
    global_mask = protein_mask(cart_sites=pandda.get_reference_dataset().get_calpha_sites(), grid_spacing=pandda.get_reference_grid().get_grid_spacing(), distance_cutoff=6)
    pandda.get_reference_grid().set_global_mask(global_mask)

if pandda.get_reference_grid().get_masked_grid_points() is None:
    pandda.mask_resampled_reference_grid()

# ============================================================================>
#####
# Scale and Process All Data
#####
# Read in all of the files
# Scale and align maps to reference
# TODO Revisit Scaling
# ============================================================================>

if not (pandda.get_used_files() and pandda.get_used_datasets()):
    pandda.add_input_files(raw_file_pairs)
    pandda.load_input_datasets()
    pandda.scale_and_align_datasets()

# This always needs to be done, even if the objects are retrieved from a pickle
pandda.load_map_handlers()

if not (pandda.get_used_files() and pandda.get_used_datasets()):
    pandda.filter_datasets()

if not pandda.get_maps():
    pandda.extract_map_values()

# ============================================================================>
#####
# ANALYSIS
#####
# Calculate the moments of the distributions at the grid points
# Use the means and the stds to convert the maps to z-maps
# Use the local mask to look for groups of significant z-values
# ============================================================================>

if not (pandda.get_mean_map() and pandda.get_stds_map() and pandda.get_skew_map() and pandda.get_kurt_map()):
    pandda.calculate_map_statistics()

if not pandda.get_z_maps():
    pandda.normalise_maps_to_z_maps()

if not pandda.get_modified_z_maps():
    pandda.post_process_z_maps()

# ============================================================================>

pandda.pickle_the_pandda()

# ============================================================================>

# ============================================================================>

# ============================================================================>

# ============================================================================>

hits = {}
hit_count = {}
for z_cutoff in [3,5,7,8]:
    hits[z_cutoff] = pandda.extract_modz_values_and_coords(z_cutoff=z_cutoff)
    print('{!s} hits found for z={!s}'.format(len(hits[z_cutoff]), z_cutoff))

    if hits[z_cutoff]:
        bin_count = numpy.bincount([h[0] for h in hits[z_cutoff]])
        freq_count = sorted(zip(numpy.nonzero(bin_count)[0],bin_count[numpy.nonzero(bin_count)[0]]), key=lambda tup: tup[1])
        dataset_count = [(num, pandda.get_used_files()[d_num][0]) for d_num, num in freq_count]
        hit_count[z_cutoff] = dataset_count


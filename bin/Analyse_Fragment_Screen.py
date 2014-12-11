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
#####
# Initialise and Settings
#####
# ============================================================================>
# TODO Change the name of the main processor object
# TODO Change the way the datasets are scaled
# TODO Choose the resolution limit better
# TODO Add a buffer to the grid

pandda = multi_dataset_analyser(outdir=input_dir, keep_in_memory=False)
pandda.set_map_type(map_type='2mFo-DFc')
pandda.set_map_scaling(scaling='sigma')
pandda.set_cut_resolution(d_min=2)
pandda.set_border_padding(6)

# ============================================================================>
#####
# Set Reference Dataset
#####
# Select the reference dataset
# TODO Choose the reference dataset better
# ============================================================================>

pandda.load_reference_dataset(ref_pdb=ref_pdb, ref_mtz=ref_mtz)
pandda.load_reference_map_handler()

# ============================================================================>
#####
# Create Sample Grid
#####
# Choose the resolution limit for the maps
# Create reference grid based on the reference structure
#
# ============================================================================>

if pandda.get_reference_grid() is None:
    pandda.create_reference_grid(res_factor=0.25, include_origin=True, buffer=pandda.get_border_padding())

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
    local_mask = spherical_mask(grid_spacing=pandda.get_reference_grid().get_grid_spacing(), distance_cutoff=1.2, grid_jump=2)
    pandda.get_reference_grid().set_local_mask(local_mask)

if pandda.get_reference_grid().get_global_mask() is None:
    print('===================================>>>')
    print('Generating Global Mask')
    global_mask = protein_mask(cart_sites=pandda.get_reference_dataset().get_backbone_sites(), grid_spacing=pandda.get_reference_grid().get_grid_spacing(), max_dist=pandda.get_border_padding(), min_dist=1.5)
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
pandda.collect_dataset_variation_statistics()

if not (pandda.get_used_files() and pandda.get_used_datasets()):
    pandda.filter_datasets()

if not pandda.get_maps():
    pandda.extract_map_values()

# ============================================================================>
#####
# DATA PROCESSING
#####
# Calculate the moments of the distributions at the grid points
# Use the means and the stds to convert the maps to z-maps
# Use the local mask to look for groups of significant z-values
# ============================================================================>

if not (pandda.get_mean_map() and pandda.get_stds_map() and pandda.get_skew_map() and pandda.get_kurt_map()):
    pandda.calculate_map_statistics()

if not pandda.get_z_maps():
    pandda.normalise_maps_to_z_maps()
    pandda.update_pandda_size(tag='After Creating Z-Maps')

renormalise = False
if renormalise and not pandda.get_renormalised_z_maps() and not pandda.get_modified_z_maps():
    pandda.renormalise_z_maps()

if not pandda.get_modified_z_maps():
    pandda.post_process_z_maps(local_mask_function=None, renormalise=renormalise)
    pandda.update_pandda_size(tag='After Post-Processing Z-Maps')

# ============================================================================>
#####
# STORING FOR REUSE
#####
# Pickle all of the large arrays so they can be reloaded
# ============================================================================>

pandda.pickle_the_pandda()

pandda.update_pandda_size(tag='After Pickling...')

# ============================================================================>
#####
# FINISHED ------------------------------>>>
#####
# ============================================================================>

from scitbx.array_family import flex

# XXX Setting Analysis Variables XXX
# Minimum size of cluster
min_cluster_size = 3
# Z cutoff for maps
z_cutoff = 3
# Cutoff for separation of clusters
cluster_cutoff = 5
# Clustering Methods
cluster_criterion='distance'
cluster_metric='euclidean'
cluster_method='average'

output_cluster_nums_results = []
output_cluster_nums_file = os.path.join(input_dir, 'cluster_nums_results.csv')

output_cluster_size_results = []
output_cluster_size_file = os.path.join(input_dir, 'cluster_size_results.csv')

output_cluster_hits_results = []
output_cluster_hits_file = os.path.join(input_dir, 'cluster_roc_results.csv')

# Cutoff for how close the "discovered" points have to be to the centroid of the ligand to be correct
max_peak_dist = 5

xtal_dict = dict([(i, file_pair[0][84:88]) for i, file_pair in enumerate(pandda.get_used_files())])
solution_file = os.path.join(input_dir, 'BAZ2BA.summary')
solutions = [line.strip().split(', ') for line in open(solution_file, 'r').readlines()]

correct_loc_dict = {}
[correct_loc_dict.setdefault(xtal,[]) for xtal, model, smile, x, y, z in solutions]
[correct_loc_dict[xtal].append((int(x), int(y), int(z))) for xtal, model, smile, x, y, z in solutions]
#[correct_loc_dict.setdefault(xtal,[]) for xtal, model, smile, x, y, z in solutions if smile != 'OCCO']
#[correct_loc_dict[xtal].append((int(x), int(y), int(z))) for xtal, model, smile, x, y, z in solutions if smile != 'OCCO']
# XXX -------------------------- XXX

# ============================================================================>
#####
# ANALYSIS
#####
# Analyse the processed data
# ============================================================================>

print('===================================>>>')
print('Getting Clusters of Z-scores')

for aaaaaaaaahhhhhh in [1]:
    print('===================================>>>')
    hits = pandda.cluster_modz_values(z_cutoff=z_cutoff, cluster_cutoff=cluster_cutoff,
            cluster_criterion=cluster_criterion, cluster_metric=cluster_metric, cluster_method=cluster_method)

    print('Clustering Cutoff: {!s}'.format(cluster_cutoff))
    cluster_num = [(k, xtal_dict[k], len(hits[k])) for k in hits.keys() if hits[k]]
    print('\n'.join(map(str,cluster_num)))

    # XXX Pulling out cluster nums for each dataset
    output_cluster_nums_results = cluster_num
    # XXX

    cluster_total = sum([a[2] for a in cluster_num])
    print('Total Datasets with Clusters: {!s}'.format(len(cluster_num)))
    print('Total Clusters: {!s}'.format(cluster_total))

    print('===================================>>>')
    print('Sorting Clusters of Z-scores')

    processed_clusters = []

    # Process the output clusters
    for d_num in sorted(hits.keys()):
        dataset_hits = hits[d_num]
        if dataset_hits == None:
            continue
        for c_num in sorted(dataset_hits.keys()):
            cluster_hits = dataset_hits[c_num]
            # Calculate Cluster Statistics
            c_size = len(cluster_hits)
            c_mean = numpy.mean([h[0] for h in cluster_hits])
            c_cent = tuple(map(int,numpy.mean([h[1] for h in cluster_hits], axis=0).tolist()))
            c_zmax = max([h[0] for h in cluster_hits])
            c_peak = tuple(map(int,cluster_hits[[h[0] for h in cluster_hits].index(c_zmax)][1]))
            if c_size > min_cluster_size:
#                print '{!s:3}'.format(d_num), '-', '{!s:3}'.format(c_num), '-', '{!s:5}'.format(c_size),
#                print '-', '{!s:>6}'.format(round(c_mean,3)), '{!s:>15}'.format(c_cent),
#                print '-', '{!s:>6}'.format(round(c_zmax,3)), '{!s:>15}'.format(c_peak)
                processed_clusters.append((d_num, c_num, c_size, c_mean, c_cent, c_zmax, c_peak))
    #                                    0       1       2       3       4       5       6

    # XXX Pulling out sizes etc for each cluster
    output_cluster_size_results = processed_clusters
    # XXX

    # Sort the clusters by size
    sorted_by_size = sorted(processed_clusters, key=lambda tup: tup[2], reverse=True)

    print('===================================>>>')
    print 'Sorted by Size:'
    print('===================================>>>')
    print '\n'.join(map(str,[(xtal_dict[x[0]],x[4]) for i,x in enumerate(sorted_by_size)])[0:20])
    print '...'

    # Sort the clusters by mean of cluster
    sorted_by_mean = sorted(processed_clusters, key=lambda tup: tup[3], reverse=True)

    print('===================================>>>')
    print 'Sorted by Mean Value:'
    print('===================================>>>')
    print '\n'.join(map(str,[(xtal_dict[x[0]],x[4]) for i,x in enumerate(sorted_by_mean)])[0:20])
    print '...'

    # Sort the clusters by peak of cluster
    sorted_by_peak = sorted(processed_clusters, key=lambda tup: tup[5], reverse=True)

    print('===================================>>>')
    print 'Sorted by Peak Map Value:'
    print('===================================>>>')
    print '\n'.join(map(str,[(xtal_dict[x[0]],x[6]) for i,x in enumerate(sorted_by_peak)])[0:20])
    print '...'

    # Iterate through list and check if point is near to correct point in dictionary

    for iteration, sorted_list in enumerate([sorted_by_size, sorted_by_mean, sorted_by_peak]):
        iter_type = ['size_c','mean_z','peak_z'][iteration]
        iter_point = [4,4,6][iteration]
        rank_val = [2,3,5][iteration]
        print('===================================>>>')
        print('Looking for Hits: {!s}'.format(iter_type))
        # Iterate through the sorted lists of dataset points
        for rank, list_item in enumerate(sorted_list):
            # Unpack cluster information
            d_num, c_num, c_size, c_mean, c_cent, c_zmax, c_peak = list_item
            # Check to see if list_item is in a dataset with a bound ligand
            if xtal_dict[d_num] not in correct_loc_dict.keys():
                is_correct = 0
            else:
                # Check if the selected point is near a bound ligand
                dists = []
                for cart_point in correct_loc_dict[xtal_dict[d_num]]:
                    dists.append((flex.double(cart_point) - flex.double(list_item[iter_point])).norm())
                if [1 for d in dists if d < max_peak_dist]:
                    is_correct = 1
                else:
                    is_correct = 0

            # XXX Pulling out correctness of clusters - Store the (sorting_type, xtal, sorting_value, correctness, cluster_size, cluster_zmean, cluster_zmax)
            output_cluster_hits_results.append((iter_type, xtal_dict[d_num], list_item[rank_val], is_correct, c_size, c_mean, c_zmax))
            # XXX

if 0:
    print('===================================>>>')
    print('Writing Cluster Summaries and ROC Data')

    # Write the cluster results to file
    with open(output_cluster_nums_file, 'w') as output_fh:
        for output_tup in output_cluster_nums_results:
            output_line = ', '.join(map(str, output_tup))
            output_fh.write(output_line+'\n')

    with open(output_cluster_size_file, 'w') as output_fh:
        for output_tup in output_cluster_size_results:
            output_line = ', '.join(map(str, output_tup))
            output_fh.write(output_line+'\n')

    with open(output_cluster_hits_file, 'w') as output_fh:
        for output_tup in output_cluster_hits_results:
            output_line = ', '.join(map(str, output_tup))
            output_fh.write(output_line+'\n')

# ============================================================================>
# ======================================>
# Manual Analyses - These only need processing once
# ======================================>
# ============================================================================>

if 0:
    print('===================================>>>')
    print('Getting Grid Points distributions')

    grid_indexer = pandda.get_reference_grid().get_grid_indexer()

    with open(os.path.join(pandda.outdir,'point_dist.txt'), 'w') as fh:
        fh.write('grid_point, dataset, map_val\n')
        # Write each map value for each grid point, labelled with the dataset
        for grid_point in [(25,80,60),(25,85,60),(25,80,65),(25,85,65)]:
            map_vals = [ma[grid_indexer(grid_point)] for ma in pandda.get_maps()]
            for d, map_val in enumerate(map_vals):
                out_list = [grid_point,  d, map_val]
                out_line = ', '.join(map(str,out_list)) + '\n'
                fh.write(out_line)

# ============================================================================>

if 0:
    print('===================================>>>')
    print('Writing Datasets Summaries')

    data_summary = pandda.get_dataset_variation_summary()

    with open(os.path.join(pandda.outdir,'dataset_summary.txt'), 'w') as fh:
        fh.write('res_low, res_high, a, b, c, alpha, beta, gamma, cell volume, rmsd to reference structure\n')
        # Write out parameters for each dataset
        for d_num in range(203):
            out_list = list(data_summary._resolution_pairs[d_num]) + list(data_summary._cell_params[d_num]) + [data_summary._cell_vols[d_num], data_summary._rmsds[d_num]]
            out_line = ', '.join(map(str,out_list)) + '\n'
            fh.write(out_line)

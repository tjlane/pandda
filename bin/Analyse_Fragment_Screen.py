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

pandda = multi_dataset_analyser(outdir=input_dir, keep_maps_in_memory=False)
pandda.set_map_type(map_type='2mFo-DFc')
pandda.set_map_scaling(scaling='none')
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
# Load, Scale and Filter All Data
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
# TODO Make these load dynamically...
if not pandda.has_maps():
    pandda.load_map_handlers()

# Analyses the crystallographic and structural variability of the datasets
pandda.calculate_mean_structure_and_protein_masks(deviation_cutoff=0.5)
pandda.collect_dataset_variation_statistics()

# Filters and Clusters the datasets using the information gathered above...
if not (pandda.get_used_files() and pandda.get_used_datasets()):
    pandda.filter_datasets()

if not pandda.has_maps():
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

if not pandda.has_z_maps():
    pandda.normalise_maps_to_z_maps()
    pandda.update_pandda_size(tag='After Creating Z-Maps')

if not pandda.has_modified_z_maps():
    pandda.post_process_z_maps(local_mask_function=None)
    pandda.update_pandda_size(tag='After Post-Processing Z-Maps')

# ============================================================================>
#####
# STORING FOR REUSE
#####
# Pickle all of the large arrays so they can be reloaded
# ============================================================================>

pandda.pickle_the_pandda()

pandda.update_pandda_size(tag='After Processing')

# ============================================================================>
#####
# SUMMARIES ------------------------------>>>
#####
# ============================================================================>

# TODO Write Summaries to be written the log/output file...

# ============================================================================>
#####
# FINISHED ------------------------------>>>
#####
# ============================================================================>

from scitbx.array_family import flex
from Giant.Stats.Cluster import cluster_data, combine_clusters

# XXX Setting Analysis Variables XXX

# Minimum size of cluster
min_cluster_size = 2
# Z cutoff for maps
z_cutoff = 3
# Cutoff for separation of clusters (sqrt((2x)**2 + (2y)**2 + (2z)**2))
cluster_cutoff = 1.5 * numpy.sqrt(3*4*pandda.get_reference_grid().get_grid_spacing()**2)
# Clustering Methods
cluster_criterion='distance'
cluster_metric='euclidean'
cluster_method='average'
# Cutoff for how close the "discovered" points have to be to the centroid of the ligand to be correct
max_peak_dist = 5

output_cluster_nums_results = []
output_cluster_nums_file = os.path.join(input_dir, 'cluster_per_dataset.csv')

output_cluster_summ_results = []
output_cluster_summ_file = os.path.join(input_dir, 'cluster_individual_summaries.csv')

output_cluster_hits_results = []
output_cluster_hits_file = os.path.join(input_dir, 'cluster_roc_results.csv')

dnum_to_xtal_dict = dict([(i, file_pair[0][84:88]) for i, file_pair in enumerate(pandda.get_used_files())])
xtal_to_dnum_dict = dict([(file_pair[0][84:88], i) for i, file_pair in enumerate(pandda.get_used_files())])
solution_file = os.path.join(input_dir, 'BAZ2BA.summary')
solutions = [line.strip().split(', ') for line in open(solution_file, 'r').readlines()]

correct_ligands = {}
[correct_ligands.setdefault(xtal_to_dnum_dict[xtal],[]) for xtal, model, smile, x, y, z in solutions if xtal in xtal_to_dnum_dict.keys()]
[correct_ligands[xtal_to_dnum_dict[xtal]].append([(int(x), int(y), int(z)), smile]) for xtal, model, smile, x, y, z in solutions if xtal in xtal_to_dnum_dict.keys()]

# XXX -------------------------- XXX

# ============================================================================>
#####
# ANALYSIS
#####
# Analyse the processed data
# TODO Create a `hit processor` class
# ============================================================================>

initial_results_dir = os.path.join(pandda.outdir, 'initial_hits')
if not os.path.exists(initial_results_dir): os.mkdir(initial_results_dir)

# ============================================================================>

if 1:

    print('===================================>>>')
    print('Getting Clusters of Z-scores')

    print('===================================>>>')
    hits = pandda.cluster_modz_values(z_cutoff=z_cutoff, cluster_cutoff=cluster_cutoff,
            cluster_criterion=cluster_criterion, cluster_metric=cluster_metric, cluster_method=cluster_method)

    print('===================================>>>')
    print('Clustering Cutoff: {!s}'.format(cluster_cutoff))

    cluster_num = [(k, dnum_to_xtal_dict[k], len(hits[k])) for k in hits.keys() if hits[k]]
    print('Total Datasets with Clusters: {!s}'.format(len(cluster_num)))
    cluster_total = sum([a[2] for a in cluster_num])
    print('Total Clusters: {!s}'.format(cluster_total))

    #print('\n'.join(map(str,cluster_num)))

    # XXX Pulling out cluster nums for each dataset
    output_cluster_nums_results = cluster_num
    # XXX

# ============================================================================>

if 1:

    print('===================================>>>')
    print('Processing Clusters of Z-scores')

    # Use the new cluster objects
    d_clusters = [cluster_data(hits[d_num]) if hits[d_num] else None for d_num in xrange(203)]

    cluster_summaries = []
    filtered_d_clusters = []

    for d_num, clust in enumerate(d_clusters):

        if clust == None:
            filtered_d_clusters.append(None)
            continue

        if 1:
            # Check how far apart the clusters are spaced
            dists = []
            c_num = len(clust.get_centroids())
            for cent_a in clust.get_centroids():
                for cent_b in clust.get_centroids():
                    if cent_a == cent_b: continue
                    dists.append((flex.double(cent_a) - flex.double(cent_b)).norm()*pandda.get_reference_grid().get_grid_spacing())
            if not dists: print 'Clusters, Dataset {:4}:'.format(d_num), '\tNum: {:4}'.format(c_num)
            else:         print 'Clusters, Dataset {:4}:'.format(d_num), '\tNum: {:4}'.format(c_num), '\tMin Spacing: {:5.3}'.format(min(dists)), '\tMax Spacing: {:5.3}'.format(max(dists))

        # Pull out a cluster summary
        for c_num, c_key in enumerate(clust.get_keys()):
            cluster_summaries.append([d_num, c_key] + clust.get_sizes([c_num]) + clust.get_means([c_num]) + clust.get_maxima([c_num]))

        # Filter out the small clusters
        new_clust = clust.create_new_cluster_from_mask(mask=[1 if s>=min_cluster_size else 0 for s in clust.get_sizes()])
        filtered_d_clusters.append(new_clust)

        if new_clust == None:
            continue

        # Link the maps to the output directory
        hit_subdir = os.path.join(initial_results_dir, 'Dataset-{:04}'.format(d_num))
        if not os.path.exists(hit_subdir):
            os.mkdir(hit_subdir)
            os.symlink(pandda.get_used_files()[d_num][0], os.path.join(hit_subdir, 'initial_model.pdb'.format(d_num)))
            os.symlink(pandda.get_used_files()[d_num][1], os.path.join(hit_subdir, 'initial_model.mtz'.format(d_num)))

        # Translate GP to array index
        grid_indexer = pandda.get_reference_grid().get_grid_indexer()

        # Create maps of the high z-value points (significant points)
        highz_points = []; [highz_points.extend(pg) for pg in new_clust.get_points()]
        highz_map_array = numpy.zeros(pandda.get_reference_grid().get_grid_size_1d(), dtype=int)
        highz_map_array.put(map(grid_indexer, highz_points), range(1,len(highz_points)+1))
        d_highz_map_file = os.path.join(hit_subdir, '{!s}-high_z_vals.ccp4'.format(dnum_to_xtal_dict[d_num]))
        pandda.write_array_to_map(d_highz_map_file, flex.double(highz_map_array.tolist()))

        # Create maps of the clusters (map values are how many points are in the cluster)
        clust_map_array = numpy.zeros(pandda.get_reference_grid().get_grid_size_1d(), dtype=int)
        clust_map_array.put(map(grid_indexer, new_clust.get_peaks()), new_clust.get_sizes())
        d_cluster_map_file = os.path.join(hit_subdir, '{!s}-cluster_map.ccp4'.format(dnum_to_xtal_dict[d_num]))
        pandda.write_array_to_map(d_cluster_map_file, flex.double(clust_map_array.tolist()))

    # Combine all of the clusters into one cluster object - clusters will be marked by dataset
    combined_cluster = combine_clusters(filtered_d_clusters)

    # XXX hits/ROC information for the clusters
    output_cluster_summ_results = cluster_summaries
    # XXX

    # TODO Look for common points in the datasets?

if 1:

    # Iterate through the correct ligands and check that there is a cluster near them
    print('===================================>>>')
    print('TRAINING - Checking Hits have been identified')

    missed_datasets = 0
    missed_ligands = 0
    missed_glycols = 0
    minimum_dists = []

    for d_num in sorted(correct_ligands.keys()):
        print 'Checking for the ligands in Dataset {!s}'.format(d_num)

        # Get the list of ligands that we should be detecting
        dataset_ligands = correct_ligands[d_num]

        # Get the clusters of identified points for this dataset
        clust = filtered_d_clusters[d_num]

        for loc_cart_tuple, smile in dataset_ligands:

            loc_grid_tuple = tuple(flex.double(loc_cart_tuple)/pandda.get_reference_grid().get_grid_spacing())
            print '\tLooking for {!s} at {!s}'.format(smile, loc_grid_tuple)

            # Check to see if there are any identified points in this dataset
            if clust==None:
                print "\t\tTHIS DATASET WASN'T IDENTIFIED AS INTERESTING - MISSED IT"
                missed_datasets += 1
                continue

            # Calculate the distances between ligands and the identified clusters
            dists = []
            found_it = False
            for cart_point in clust.get_peaks():
                dists.append((flex.double(cart_point) - flex.double(loc_grid_tuple)).norm())
                #print '\t\tHow about {!s}?'.format(cart_point)
                if dists[-1] <  max_peak_dist/pandda.get_reference_grid().get_grid_spacing():
                    if found_it==True:
                        pass
                    else:
                        print '\t\t\tFOUND IT!'
                        found_it = True

            # Check to see if the ligands has been found
            if not found_it:
                print '\t\t\tMISSED IT...'
                if smile == 'OCCO':
                    missed_glycols += 1
                else:
                    missed_ligands += 1

            if smile == 'OCCO':
                pass
            else:
                minimum_dists.append(min(dists)*pandda.get_reference_grid().get_grid_spacing())

    print('----------------------------------->>>')
    print('MISSED DATASETS: {!s}'.format(missed_datasets))
    print('MISSED LIGANDS: {!s}'.format(missed_ligands))
    print('MISSED GLYCOLS: {!s}'.format(missed_glycols))
    print('----------------------------------->>>')
    print('MINIMUM DISTANCES TO LIGAND (A): \n{!s}'.format(map(int,minimum_dists)))
    print('----------------------------------->>>')

if 1:

    print('===================================>>>')
    print('Sorting Clusters')

    # Sort the clusters by size
    size_sorted_indices = combined_cluster.sort(sorting_data='sizes', sorting_function=None, decreasing=True)
    # Sort the clusters by mean
    mean_sorted_indices = combined_cluster.sort(sorting_data='values', sorting_function=numpy.mean, decreasing=True)
    # Sort the clusters by max
    max_sorted_indices = combined_cluster.sort(sorting_data='values', sorting_function=max, decreasing=True)

    print('===================================>>>')
    print('TRAINING - Calculating ROC Curves')

    # Output ROC information
    roc_testing = []

    for test_num, (c_indices, rank_vals) in enumerate([(size_sorted_indices, combined_cluster.get_sizes(size_sorted_indices)),
                                                       (mean_sorted_indices, combined_cluster.get_means(mean_sorted_indices)),
                                                       (max_sorted_indices,  combined_cluster.get_maxima(max_sorted_indices))]):

        test_type = ['c_size','z_mean','z_peak'][test_num]

        print('===================================>>>')
        print 'Ranking by: {!s}'.format(test_type)
        print('===================================>>>')

        # List to see when a ligand is identified - avoids multiple identifications
        identified_ligands = []

        # Calculate ROC Curves etc - Iterate through clusters and check if near ligand
        for c_rank, c_index in enumerate(c_indices):
            # Pull out information to identify the cluster
            d_num, c_num = combined_cluster.get_keys([c_index])[0]
            c_centroid = combined_cluster.get_centroids([c_index])[0]
            print 'Checking Dataset {:4}, Cluster {:4}, Rank: {:4}, Val: {:6.2f}'.format(d_num, c_num, c_rank, rank_vals[c_rank])

            is_a_ligand = 0
            skip = False

            # Check if there is a ligand in the dataset
            if d_num not in correct_ligands.keys():
#                print '\t0 - No Ligand in this dataset'
                pass
            else:
                # Get the ligands present in the dataset
                dataset_ligands = correct_ligands[d_num]

                dists = []

                for l_num, (loc_cart_tuple, smile) in enumerate(dataset_ligands):
                    # Don't bother for Glycols
                    if smile == 'OCCO':
                        continue

                    # Convert to grid index for distance comparison
                    loc_grid_tuple = tuple(flex.double(loc_cart_tuple)/pandda.get_reference_grid().get_grid_spacing())
                    # Calculate distance to ligand
                    dists.append((flex.double(loc_grid_tuple) - flex.double(c_centroid)).norm())

                    # Check if correct
                    if dists[-1] < max_peak_dist/pandda.get_reference_grid().get_grid_spacing():
                        is_a_ligand = 1

                    # If it's a ligand, check to see if it's been previously identified
                    if is_a_ligand:
                        if (d_num, l_num) in identified_ligands:
                            print '\t. - Ligand already found: Dataset {!s}'.format(d_num)
                            skip = True
                        else:
                            print "\t1 - Ligand found: Dataset {!s}".format(d_num)
                            identified_ligands.append((d_num, l_num))

            if not skip:
                roc_testing.append([test_type, dnum_to_xtal_dict[d_num], d_num, c_num, is_a_ligand, c_rank, rank_vals[c_rank]])

        print 'Ligands Identified:', identified_ligands

        # XXX hits/ROC information for the clusters
        output_cluster_hits_results = roc_testing
        # XXX

# ============================================================================>

if 0:
    print('===================================>>>')
    print('Writing Cluster Summaries and ROC Data')

    for o_file, o_contents in [(output_cluster_nums_file, output_cluster_nums_results),
                               (output_cluster_summ_file, output_cluster_summ_results),
                               (output_cluster_hits_file, output_cluster_hits_results)]:

        # Write the results to file
        with open(o_file, 'w') as output_fh:
            print 'Writing {!s}'.format(o_file)
            for output_tup in o_contents:
                output_line = ', '.join(map(str, output_tup))
                output_fh.write(output_line+'\n')

# ============================================================================>
# ======================================>
# Manual Analyses - These only need processing once
# ======================================>
# ============================================================================>
analyses_dir = os.path.join(pandda.outdir, 'manual_analyses')
if not os.path.exists(analyses_dir): os.mkdir(analyses_dir)
# ============================================================================>

if 0:
    print('===================================>>>')
    print('Getting Grid Points distributions')

    grid_indexer = pandda.get_reference_grid().get_grid_indexer()

    with open(os.path.join(analyses_dir,'point_dist.csv'), 'w') as fh:
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
    print('Writing Dataset Crystal Summaries')

    data_summary = pandda.get_dataset_variation_summary()

    with open(os.path.join(analyses_dir,'dataset_xtal_summary.csv'), 'w') as fh:
        fh.write('res_low, res_high, a, b, c, alpha, beta, gamma, cell volume, rmsd to reference structure\n')
        # Write out parameters for each dataset
        for d_num in range(203):
            out_list = list(data_summary._resolution_pairs[d_num]) + list(data_summary._cell_params[d_num]) + [data_summary._cell_vols[d_num], data_summary._rmsds[d_num]]
            out_line = ', '.join(map(str,out_list)) + '\n'
            fh.write(out_line)

# ============================================================================>

if 0:
    print('===================================>>>')
    print('Writing Dataset Map Summaries')

    data_summary = pandda.get_dataset_variation_summary()

    with open(os.path.join(analyses_dir,'dataset_map_summary.csv'), 'w') as fh:
        fh.write('headers\n')
        # Write out parameters for each dataset
        for d_num in range(203):
            out_list = 'a'
            out_line = ', '.join(map(str,out_list)) + '\n'
            fh.write(out_line)

# ============================================================================>

if 0:
    print('===================================>>>')
    print('Calculating Deviations of C-alphas between structures')

    rms = lambda vals: numpy.sqrt(numpy.mean(numpy.abs(vals)**2))
    norm = lambda vals: numpy.sqrt(numpy.sum(numpy.abs(vals)**2))

    # Pull all c-alpha sites for each structure
    all_sites = numpy.array([d.transform_points_to_reference(d.get_calpha_sites()) for d in pandda.get_all_datasets()])
    # Calculate the mean x,y,z for each c-alpha
    mean_sites = numpy.mean(all_sites, axis=0)
    # Differences from the mean
    diff_sites = all_sites - mean_sites
    # Euclidean norms of the distances moved
    diff_norms = numpy.apply_along_axis(norm, axis=2, arr=diff_sites)

    with open(os.path.join(analyses_dir,'calpha_variation.csv'), 'w') as fh:
        for row in diff_norms:
            out_list = row.round(3).tolist()
            out_line = ', '.join(map(str,out_list)) + '\n'
            fh.write(out_line)

    print 'Largest deviation from the mean site: {!s}'.format(diff_norms.max())
    print 'Average deviation from the mean site: {!s}'.format(diff_norms.mean())

# ============================================================================>

if 0:
    print('===================================>>>')
    print('Clustering the Refined Structures')

    distance_matrix = []
    for d1 in pandda.get_all_datasets():
       distance_matrix.append([d1.transform_points_to_reference(d1.get_calpha_sites()).rms_difference(d2.transform_points_to_reference(d2.get_calpha_sites())) for d2 in pandda.get_all_datasets()])

    distance_matrix = numpy.array(distance_matrix)

    with open(os.path.join(analyses_dir,'calpha_distance_matrix.csv'), 'w') as fh:
        for row in distance_matrix:
            out_list = row.round(3).tolist()
            out_line = ', '.join(map(str,out_list)) + '\n'
            fh.write(out_line)

# ============================================================================>

if 1:
    print('===================================>>>')
    print('...')


# # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # #

###from cctbx import uctbx, maptbx, masks
###
#### Have a look at flood-filling the modified z-maps
###
#### Unit Cell for periodicity - also an object called a non-crystallographic-unit-cell
#### Get the size of the cartesian grid in angstroms
###uc = uctbx.unit_cell('{!s} {!s} {!s} 90 90 90'.format(*pandda.get_reference_grid().get_cart_max()))
###
#### Gridding object
#### Get the size of the cartesian grid in grid points
###gridding = maptbx.crystal_gridding(unit_cell=uc, pre_determined_n_real=pandda.get_reference_grid().get_grid_size())
###
#### Get the data from the pandda
###z_map = pandda.get_z_maps()[7].as_dense_vector()
###flex_object = flex.int([1 if v>3 else 0 for v in z_map])
#### Reshape the data to a grid
###flex_object.reshape(flex.grid(gridding.n_real()))
###
###ff = masks.flood_fill(flex_object, uc)

# # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # #

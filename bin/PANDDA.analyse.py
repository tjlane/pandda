#!/home/npearce/bin/Developer/cctbx/cctbx_build/bin/cctbx.python

import os, sys, glob, time
import numpy

from scitbx.array_family import flex

from PANDDAs.Parser import build_pandda_parser
from PANDDAs.Main import multi_dataset_analyser, spherical_mask, atomic_mask

# ============================================================================>
#####
# MANUAL SETTINGS
#####
# ============================================================================>

testing = False

# ============================================================================>
#####
# Parse Input Args
#####
# ============================================================================>

pandda_arg_parser = build_pandda_parser()
if not testing:
    pandda_args = pandda_arg_parser.parse_args()
else:
    pandda_args = pandda_arg_parser.parse_args(['--data-dirs','/home/npearce/Analyses/FragmentSoaks/BAZ2BA-Feb13-Soak/ProcessedFragmentSoak/BAZ2BA-*/1-apo',
                                                '--outdir','/home/npearce/Analyses/FragmentSoaks/BAZ2BA-Feb13-Soak/pandda',
                                                '--pdb-style','apo-BAZ2BA-*-refmac.pdb',
                                                '--ref-pdb','/home/npearce/Analyses/FragmentSoaks/BAZ2BA-Feb13-Soak/reference.pdb',
                                                '--ref-mtz','/home/npearce/Analyses/FragmentSoaks/BAZ2BA-Feb13-Soak/reference.mtz',
                                                '--cpus', '6', '--verbose'])

for key in pandda_args.__dict__.keys():
    print '{!s:>30}\t\t{!s}'.format(key, pandda_args.__dict__[key])

# ============================================================================>
#####
# Manual Settings
#####
# ============================================================================>

# Minimum size of cluster
min_cluster_volume = 10
# Z cutoff for maps
z_cutoff = 3
# Cutoff for separation of clusters (sqrt((2x)**2 + (2y)**2 + (2z)**2))
clustering_cutoff = 5
# Clustering Methods
clustering_criterion='distance'
clustering_metric='euclidean'
clustering_method='average'

# Cutoff for how close the "discovered" points have to be to the centroid of the ligand to be correct
max_peak_dist = 5

# ============================================================================>
#####
# Initialise and Settings
#####
# ============================================================================>
# TODO Change the name of the main processor object
# TODO Change the way the datasets are scaled
# TODO Choose the resolution limit better

pandda = multi_dataset_analyser(pandda_args)

pandda.set_map_type(map_type='2mFo-DFc')
pandda.set_map_scaling(scaling='none')
pandda.set_cut_resolution(d_min=2)
pandda.set_border_padding(6)

pandda.run_pandda_init()

# Build and load list of fragment screen datasets
input_files = pandda.build_input_list()

if (not pandda.get_all_datasets()) and (not input_files):
    # No datasets loaded - exit
    raise SystemExit('NO DATASETS LOADED')

elif input_files:

    # ============================================================================>
    #####
    # Add new files and load datasets
    #####
    # ============================================================================>

    pandda.add_new_files(input_files)
    pandda.load_new_datasets()
    pandda.initialise_masks()

    # ============================================================================>
    #####
    # Set Reference Dataset
    #####
    # Select the reference dataset
    # ============================================================================>

    if not pandda.reference_dataset():
        # Select the reference dataset
        ref_pdb, ref_mtz = pandda.select_reference_dataset()
        # Load the reference dataset
        pandda.load_reference_dataset(ref_pdb=ref_pdb, ref_mtz=ref_mtz)

    # ============================================================================>
    #####
    # Scale, Align and Filter All Data
    #####
    # TODO Revisit Scaling
    # ============================================================================>

    # Filter out datasets with different protein structures
    pandda.filter_datasets_1()

    # Scale and align the datasets to the reference
    pandda.scale_datasets()
    pandda.align_datasets()

else:
    # Rebuild the masks of the rejected datasets (quick)
    pandda.initialise_masks()
    pandda.filter_datasets_1()

# Collate many variables across the datasets to be used for filtering
pandda.collect_dataset_variation_data()

# Filter out the datasets that are not isomorphous and therefore incomparable
pandda.filter_datasets_2()

# Analyses the crystallographic and structural variability of the datasets
pandda.calculate_mean_structure_and_protein_masks(deviation_cutoff=0.5)

# Analyse the structural variation in the datasets
pandda.collect_structure_variation_data()

# ============================================================================>
#####
# Create Sample Grid
#####
# Choose the resolution limit for the maps
# Create reference grid based on the reference structure
#
# ============================================================================>

if pandda.reference_grid() is None:
    pandda.create_reference_grid(res_factor=0.25, expand_to_origin=False, buffer=pandda.get_border_padding())

# ============================================================================>
#####
# Create Local and Global Masks
#####
# Local mask used for forming groups of points around a grid point
# Global mask used for removing points in the bulk solvent regions
# ============================================================================>

if pandda.reference_grid().local_mask() is None:
    print('===================================>>>')
    print('Generating Local Mask')
    local_mask = spherical_mask(grid_spacing    = pandda.reference_grid().grid_spacing(),
                                distance_cutoff = 1.2,
                                grid_jump       = 2 )
    pandda.reference_grid().set_local_mask(local_mask)

if pandda.reference_grid().global_mask() is None:
    print('===================================>>>')
    print('Generating Global Mask')
    global_mask = atomic_mask(  cart_sites = pandda.reference_dataset().get_heavy_atom_sites(),
                                grid_size  = pandda.reference_grid().grid_size(),
                                unit_cell  = pandda.reference_grid().fake_unit_cell(),
                                max_dist   = pandda.get_border_padding(),
                                min_dist   = 1.5 )
    pandda.reference_grid().set_global_mask(global_mask)

if pandda.reference_grid().masked_grid_points() is None:
    pandda.mask_resampled_reference_grid()

# ============================================================================>
#####
# Local Map Handler Objects and Extract Map Values
#####
# ============================================================================>

pandda.load_reference_map_handler()
pandda.extract_reference_map_values()

# This always needs to be done, even if the objects are retrieved from a pickle
# TODO Make these load dynamically...
pandda.load_all_map_handlers()

# ============================================================================>
#####
# Scale the Maps, Calculate Statistical Maps
#####
# ============================================================================>

pandda.analyse_raw_maps()

if not (pandda.get_mean_map() and pandda.get_stds_map() and pandda.get_skew_map() and pandda.get_kurt_map()):
    pandda.calculate_map_statistics()

# ============================================================================>
#####
# Write Grid Point Distributions
#####
# ============================================================================>

try:
    pandda.write_grid_point_distributions(  grid_points = [ (100,  10, 32),
                                                            ( 90,  20, 36),
                                                            ( 80,  30, 40),
                                                            ( 70,  40, 44),
                                                            ( 60,  50, 48),
                                                            ( 50,  60, 52),
                                                            ( 40,  70, 56),
                                                            ( 30,  80, 60),
                                                            ( 20,  90, 64),
                                                            ( 10, 100, 68)],
                                            output_filename='point_line_distributions.csv' )
except:
    pass

# ============================================================================>
#####
# STORING FOR REUSE
#####
# Pickle all of the large arrays so they can be reloaded
# ============================================================================>

pandda.pickle_the_pandda()
pandda.update_pandda_size(tag='After Pre-processing')

# ============================================================================>
#####
# DATA PROCESSING
#####
# Calculate the moments of the distributions at the grid points
# Use the means and the stds to convert the maps to z-maps
# Use the local mask to look for groups of significant z-values
# ============================================================================>

pandda.print_clustering_settings(   z_cutoff             = z_cutoff,
                                    min_cluster_volume   = min_cluster_volume,
                                    clustering_cutoff    = clustering_cutoff,
                                    clustering_criterion = clustering_criterion,
                                    clustering_metric    = clustering_metric,
                                    clustering_method    = clustering_method  )

for d_handler in pandda.get_masked_datasets(mask_name='rejected - total', invert=True):

    if d_handler.raw_cluster_hits is not None:
        print '{!s} - ALREADY PROCESSED!'.format(d_handler.d_tag)
        continue

    pandda.log('===================================>>>', True)
    pandda.log('Calculating Z-MAPs for Dataset {!s} ({!s}/{!s})'.format(d_handler.d_tag, d_handler.d_num+1, self.get_number_of_datasets()), True)

    ##################################
    # EXTRACT MAP VALUES             #
    ##################################

    s_map = pandda.get_map(d_handler=d_handler)

    # Write the sampled map
    pandda.write_array_to_map(  output_file = d_handler.get_mtz_filename().replace('.mtz','.sampled.ccp4'),
                                map_data    = s_map)

    ##################################
    # CALCULATE Z-MAPS               #
    ##################################

    z_map = pandda.calculate_z_map(map_vals=s_map)

    # Write z map
    pandda.write_array_to_map(  output_file = d_handler.get_mtz_filename().replace('.mtz','.zmap.ccp4'),
                                map_data    = z_map)

    ##################################
    # LOOK FOR CLUSTERS OF Z-SCORES  #
    ##################################

    pandda.cluster_high_z_values(  d_handler            = d_handler,
                                   z_map                = z_map,
                                   z_cutoff             = z_cutoff,
                                   min_cluster_volume   = min_cluster_volume,
                                   clustering_cutoff    = clustering_cutoff,
                                   clustering_criterion = clustering_criterion,
                                   clustering_metric    = clustering_metric,
                                   clustering_method    = clustering_method)

#    ##################################
#    # POST-PROCESS Z-MAPS            #
#    ##################################
#
#    mod_z_map, resamp_mod_z_map = pandda.process_z_map(z_map=z_map)
#
#    # Write map
#    pandda.write_array_to_map(  output_file  = d_handler.get_mtz_filename().replace('.mtz','.processed.zvalues.ccp4'),
#                                map_data     = flex.double(mod_z_map))
#    # Write down-sampled map
#    pandda.write_array_to_map(  output_file  = d_handler.get_mtz_filename().replace('.mtz','.resamp.processed.zvalues.ccp4'),
#                                map_data     = flex.double(resamp_mod_z_map),
#                                grid_size    = pandda.reference_grid().resampled_grid_size(),
#                                grid_spacing = pandda.reference_grid().resampled_grid_spacing())


cluster_total, cluster_num, all_dataset_clusters = pandda.collate_all_clusters()

pandda.pickle_the_pandda()
pandda.update_pandda_size(tag='After Processing')

# ============================================================================>
#####
# ANALYSIS
#####
# Analyse the processed data
# TODO Create a `hit processor` class
# ============================================================================>

combined_cluster = pandda.process_z_value_clusters()

# ============================================================================>
#####
# SUMMARIES ------------------------------>>>
#####
# ============================================================================>

# TODO Write Summaries to be written the log/output file...

pandda.write_analysis_summary(output_file=None)

# ============================================================================>
#####
# XXX XXX XXX Setting Analysis Variables XXX XXX XXX
#####
# ============================================================================>

if testing:

    solution_file = os.path.join(pandda.output_handler.get_dir('root'), 'correct_hits.summary')
    if not os.path.exists(solution_file):
        print('CORRECT LIGANDS FILE DOES NOT EXIST!!!')
        testing = False
    else:
        solutions = [line.strip().split(', ') for line in open(solution_file, 'r').readlines()]

        correct_ligands = {}
        [correct_ligands.setdefault(xtal,[]) for xtal in [d.d_tag for d in pandda.get_all_datasets()]]
        [correct_ligands[xtal].append([(int(x), int(y), int(z)), smile]) for xtal, model, smile, x, y, z in solutions if xtal in correct_ligands.keys()]

# XXX XXX XXX -------------------------- XXX XXX XXX

# ============================================================================>

if testing:

    # Iterate through the correct ligands and check that there is a cluster near them
    print('===================================>>>')
    print('TRAINING - Checking Hits have been identified')

    missed_datasets = 0
    missed_ligands = 0
    missed_glycols = 0
    minimum_dists = []

    for d_tag in sorted(correct_ligands.keys()):
        print 'Checking for the ligands in Dataset {!s}'.format(d_tag)

        # Get the list of ligands that we should be detecting
        dataset_ligands = correct_ligands[d_tag]

        # Get the clusters of identified points for this dataset
        clust = pandda.get_dataset(d_tag=d_tag).clustered_hits

        for loc_cart_tuple, smile in dataset_ligands:

            loc_grid_tuple = tuple(flex.double(loc_cart_tuple)/pandda.reference_grid().grid_spacing())
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
                if dists[-1] <  max_peak_dist/pandda.reference_grid().grid_spacing():
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
                minimum_dists.append(min(dists)*pandda.reference_grid().grid_spacing())

    print('----------------------------------->>>')
    print('MISSED DATASETS: {!s}'.format(missed_datasets))
    print('MISSED LIGANDS: {!s}'.format(missed_ligands))
    print('MISSED GLYCOLS: {!s}'.format(missed_glycols))
    print('----------------------------------->>>')
    print('MINIMUM DISTANCES TO LIGAND (A): \n{!s}'.format(map(int,minimum_dists)))
    print('----------------------------------->>>')

if testing:

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
            d_tag, c_num = combined_cluster.get_keys([c_index])[0]
            c_centroid = combined_cluster.get_centroids([c_index])[0]
            print 'Checking Dataset {:4}, Cluster {:4}, Rank: {:4}, Val: {:6.2f}'.format(d_tag, c_num, c_rank, rank_vals[c_rank])

            is_a_ligand = 0
            skip = False

            # Check if there is a ligand in the dataset
            if d_tag not in correct_ligands.keys():
                pass
            else:
                # Get the ligands present in the dataset
                dataset_ligands = correct_ligands[d_tag]

                dists = []

                for l_num, (loc_cart_tuple, smile) in enumerate(dataset_ligands):
                    # Don't bother for Glycols
                    if smile == 'OCCO':
                        continue

                    # Convert to grid index for distance comparison
                    loc_grid_tuple = tuple(flex.double(loc_cart_tuple)/pandda.reference_grid().grid_spacing())
                    # Calculate distance to ligand
                    dists.append((flex.double(loc_grid_tuple) - flex.double(c_centroid)).norm())

                    # Check if correct
                    if dists[-1] < max_peak_dist/pandda.reference_grid().grid_spacing():
                        is_a_ligand = 1

                    # If it's a ligand, check to see if it's been previously identified
                    if is_a_ligand:
                        if (d_tag, l_num) in identified_ligands:
                            print '\t. - Ligand already found: Dataset {!s}'.format(d_tag)
                            skip = True
                        else:
                            print "\t1 - Ligand found: Dataset {!s}".format(d_tag)
                            identified_ligands.append((d_tag, l_num))

            if not skip:
                roc_testing.append([test_type, d_tag, d_tag, c_num, is_a_ligand, c_rank, rank_vals[c_rank]])

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
analyses_dir = pandda.output_handler.get_dir('manual_analyses')
# ============================================================================>

if 0:
    print('===================================>>>')
    print('Getting Grid Points distributions')

    grid_indexer = pandda.reference_grid().grid_indexer()

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

    data_summary = pandda.get_dataset_observations()

    with open(os.path.join(analyses_dir,'dataset_xtal_summary.csv'), 'w') as fh:
        fh.write('res_low, res_high, a, b, c, alpha, beta, gamma, cell volume, rmsd to mean, rfree, rwork\n')
        # Write out parameters for each dataset
        for d_num in range(203):
            out_list = list(data_summary.data['resolution'][d_num]) + list(data_summary.data['cell_params'][d_num]) + \
                                [data_summary.data['cell_volume'][d_num], data_summary.data['rmsd_to_mean'][d_num]] + \
                                list(data_summary.data['rfree_rwork'][d_num])
            out_line = ', '.join(map(str,out_list)) + '\n'
            fh.write(out_line)

# ============================================================================>

if 0:
    print('===================================>>>')
    print('Writing Dataset Map Summaries')

    data_summary = pandda.get_dataset_observations()

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

pandda.exit()

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
###uc = uctbx.unit_cell('{!s} {!s} {!s} 90 90 90'.format(*pandda.reference_grid().get_cart_max()))
###
#### Gridding object
#### Get the size of the cartesian grid in grid points
###gridding = maptbx.crystal_gridding(unit_cell=uc, pre_determined_n_real=pandda.reference_grid().get_grid_size())
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

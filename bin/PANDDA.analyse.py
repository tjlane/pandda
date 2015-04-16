#!/home/npearce/bin/Developer/cctbx/cctbx_build/bin/cctbx.python

import os, sys, glob, time
import numpy

from scitbx.array_family import flex
from scitbx.math import basic_statistics

from Giant.Xray.Symmetry import combine_hierarchies, generate_adjacent_symmetry_copies

from Giant.Grid.Masks import spherical_mask, atomic_mask, non_symmetrical_atomic_mask

from Giant.Stats.Cluster import cluster_data

from PANDDAs.Main import multi_dataset_analyser

def pandda_main(args):
    """Run the PANDDA algorithm, using supplied args object"""

    try:

        # ============================================================================>
        #####
        # MANUAL SETTINGS
        #####
        # ============================================================================>

        # None!

        # ============================================================================>
        #####
        # Initialise
        #####
        # ============================================================================>

        pandda = multi_dataset_analyser(args)

        # ============================================================================>
        #####
        # Initialise Settings
        #####
        # ============================================================================>

        # XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX
        # MOVE THIS INTO INIT FUNCTION AND POPULATE FROM THE ARG PARSER
        pandda.set_low_resolution(pandda.params.analysis.high_res_lower_limit)
        pandda.set_high_resolution(pandda.params.analysis.high_res_upper_limit)
        pandda.set_cut_resolution(2)
        # XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX

        pandda.run_pandda_init()

        # ============================================================================>
        #####
        # Build list of files in data directories
        #####
        # ============================================================================>

        input_files = pandda.build_input_list()

        if (not pandda.datasets.all()) and (not input_files):
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
            pandda.initialise_analysis()

            # ============================================================================>
            #####
            # Set Reference Dataset
            #####
            # Select the reference dataset
            # ============================================================================>

            if not pandda.reference_dataset():
                # Select the reference dataset
                ref_pdb, ref_mtz = pandda.select_reference_dataset(method='resolution')
                # Load the reference dataset
                pandda.load_reference_dataset(ref_pdb=ref_pdb, ref_mtz=ref_mtz)

            # ============================================================================>
            #####
            # Scale, Align and Initial-Filter All Data
            #####
            # TODO Revisit Scaling
            # ============================================================================>

            # Filter out datasets with different protein structures
            pandda.filter_datasets_1()

            # Scale and align the datasets to the reference
            pandda.scale_datasets(  ampl_label=pandda.params.maps.ampl_label,
                                    phas_label=pandda.params.maps.phas_label    )
            pandda.align_datasets(method=pandda.params.alignment.method)
    #        pandda.generate_crystal_contacts()

        else:
            # Rebuild the masks of the rejected datasets (quick)
            pandda.initialise_analysis()
            pandda.filter_datasets_1()

        # XXX XXX XXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXX XXX XXX
        # Use symmetry operations to create the symmetry mates of the reference structure
        sym_ops, sym_op_contacts, sym_hierarchies, chain_mappings = generate_adjacent_symmetry_copies(ref_hierarchy    = pandda.reference_dataset().new_structure().hierarchy,
                                                                                                      crystal_symmetry = pandda.reference_dataset().get_input().crystal_symmetry(),
                                                                                                      buffer_thickness = pandda.params.maps.border_padding+5,
                                                                                                      method=2)
        # Record the symmetry operations that generate the crystal contacts
        pandda.crystal_contact_generators = sym_ops

        # Create a combined hierarchy of the crystal contacts
        symmetry_root = combine_hierarchies(sym_hierarchies)
        symmetry_root.atoms().set_xyz(symmetry_root.atoms().extract_xyz() + pandda.get_reference_origin_shift())

        # Write out the symmetry sites
        if not os.path.exists(pandda.output_handler.get_file('reference_symmetry')):
            symmetry_root.write_pdb_file(pandda.output_handler.get_file('reference_symmetry'))
        # XXX XXX XXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXX XXX XXX

        # ============================================================================>
        #####
        # Filter and Analyse the Datasets
        #####
        # ============================================================================>

        # Collate many variables across the datasets to be used for filtering
        # Rename? - This should be the only thing that populates pandda.get_dataset_observations object - simple class?
        pandda.analyse_dataset_variability_1()
        # Filter out the datasets that are not isomorphous and therefore incomparable
        pandda.filter_datasets_2()

        # Analyses the crystallographic and structural variability of the datasets
        pandda.calculate_mean_structure_and_protein_masks(deviation_cutoff=0.5)
        pandda.analyse_dataset_variability_2()

        # Analyse the structural variation in the datasets
        pandda.analyse_structure_variability_1()

        # ============================================================================>
        #####
        # Update Settings
        #####
        # ============================================================================>

        # Update the resolution limits using the resolution limits from the datasets supplied
        pandda.set_low_resolution(   min(   pandda.get_low_resolution(),
                                            max(pandda.datasets_summary.get_data('high_res_limit'))   )   )
        pandda.set_high_resolution(  max(   pandda.get_high_resolution(),
                                            min(pandda.datasets_summary.get_data('high_res_limit'))   )   )

        # TODO POTENTIALLY MOVE THIS ...
        pandda.set_cut_resolution(   min(   pandda.get_low_resolution(),
                                            pandda.get_cut_resolution()   )   )

        pandda.log('===================================>>>')
        pandda.log('UPDATED RESOLUTION LIMITS')
        pandda.log('LOW RESOLUTION:  {!s}'.format(pandda.get_low_resolution()))
        pandda.log('HIGH RESOLUTION: {!s}'.format(pandda.get_high_resolution()))
        pandda.log('CUT RESOLUTION:  {!s}'.format(pandda.get_cut_resolution()))

        # ============================================================================>
        #####
        # Create Sample Grid
        #####
        # Create reference grid based on the reference structure
        # ============================================================================>

        if pandda.reference_grid() is None:
            pandda.create_reference_grid(   grid_spacing     = pandda.params.maps.resolution_factor*pandda.get_cut_resolution(),
                                            expand_to_origin = False,
                                            buffer           = pandda.params.maps.border_padding    )

        # ============================================================================>
        #####
        # Create Local and Global Masks
        #####
        # Local mask used for forming groups of points around a grid point
        # Global mask used for removing points in the bulk solvent regions
        # ============================================================================>

        if pandda.reference_grid().local_mask() is None:
            pandda.log('===================================>>>')
            pandda.log('Generating Local Mask')
            local_mask = spherical_mask(grid_spacing    = pandda.reference_grid().grid_spacing(),
                                        distance_cutoff = 1.2,
                                        grid_jump       = 1 )
            pandda.reference_grid().set_local_mask(local_mask)

        if pandda.reference_grid().global_mask() is None:
            pandda.log('===================================>>>')
            pandda.log('Generating Protein Mask')
            # Select the masking atoms from the reference structure
            cache = pandda.reference_dataset().get_hierarchy().atom_selection_cache()
            pro_sites_cart = pandda.reference_dataset().get_hierarchy().select(cache.selection('pepnames and not element H')).atoms().extract_xyz()
            # Generate the main protein mask
            global_mask = atomic_mask(  cart_sites   = pro_sites_cart,
                                        grid_size    = pandda.reference_grid().grid_size(),
                                        unit_cell    = pandda.reference_grid().fake_unit_cell(),
                                        max_dist     = pandda.params.maps.border_padding,
                                        min_dist     = 1.8 )
            pandda.reference_grid().set_global_mask(global_mask)

        if pandda.reference_grid().symmetry_mask() is None:
            pandda.log('===================================>>>')
            pandda.log('Generating Symmetry Mask')
            # Pull out the cartesian sites of the symmetry mates
            cache = symmetry_root.atom_selection_cache()
            sym_sites_cart = symmetry_root.select(cache.selection('pepnames and not element H')).atoms().extract_xyz()
            # Generate the symmetry mask
            symmetry_mask = non_symmetrical_atomic_mask(  cart_sites = sym_sites_cart,
                                                          grid_spacing = pandda.reference_grid().grid_spacing(),
                                                          grid_size  = pandda.reference_grid().grid_size(),
                                                          unit_cell  = pandda.reference_grid().fake_unit_cell(),
                                                          max_dist   = pandda.params.maps.border_padding,
                                                          min_dist   = 1.8 )
            pandda.reference_grid().set_symmetry_mask(symmetry_mask)

        # Print the summaries
        pandda.log('===================================>>>', True)
        pandda.log('Grid Summary: ', True)
        pandda.log(pandda.reference_grid().summary(), True)
        pandda.log(pandda.reference_grid().local_mask().summary(), True)
        pandda.log(pandda.reference_grid().global_mask().summary(), True)

        # TODO TODO TODO INCORPORATE THE COMBINATION OF MASKS INTO THIS FUNCTION
        # Create various masks to define regions of the grid by distance to the protein and symmetry copies
        pandda.mask_reference_grid(d_handler=pandda.reference_dataset())

    #    if pandda.reference_grid().masked_grid_points() is None:
    #        pandda.mask_resampled_reference_grid()
    #    # TODO TODO TODO

        # ============================================================================>
        #####
        # Pickle after creating the grid as this is time-consuming
        #####
        # ============================================================================>

        pandda.pickle_the_pandda(components=['grid'])
        pandda.update_pandda_size(tag='After Creating Grid')

        # ============================================================================>
        #####
        # Local Map Handler Objects and Extract Map Values
        #####
        # ============================================================================>

#        # DON'T ACTUALLY NEED THESE...
#        pandda.load_reference_map_handler()
#        pandda.extract_reference_map_values()

# XXX        # ============================================================================>
# XXX        #####
# XXX        # PUT IN THE RESOLUTION LOOP HERE!
# XXX        #####
# XXX        # ============================================================================>
# XXX        # BUT DO WE STILL USE THE SAME GRID? - YES SO THAT WE CAN COMPARE THE GRIDS
# XXX
# XXX
# XXX   for blah in blah:
# XXX       loop_over_resolution()

        # Select datasets based on resolution - this builds the 'selected for analysis' mask
        pandda.filter_datasets_3(resolution=pandda.get_cut_resolution())

        # ============================================================================>
        #####
        # Scale the Maps, Calculate Statistical Maps
        #####
        # ============================================================================>

        # Truncate the data to a particular resolution
        pandda.truncate_scaled_data()

        # Need to load all to calculate statistical maps
        pandda.load_and_morph_maps()

        # Only need to calculate mean if requested
        if (not pandda.stat_maps.mean_map) or pandda.args.method.recalculate_statistical_maps:
            pandda.calculate_mean_map()

        # Calculate the uncertainty of all loaded maps (needs the mean map to have been calculated)
        pandda.calculate_map_uncertainties()

        # Only need to calculate map statistics if requested
        if (not pandda.stat_maps.adj_stds_map) or pandda.args.method.recalculate_statistical_maps:
            pandda.calculate_map_statistics()

        # TODO TODO TODO
        # TODO IMPLEMENT - Collect all of the map statistics from the different dataset_handlers
        # So that the data_collections can be populated independently of when the uncertainties etc were calculated
        #pandda.analyse_dataset_variability_3()
        # TODO TODO TODO

        # ============================================================================>
        #####
        # STORING FOR REUSE
        #####
        # Pickle all of the large arrays so they can be reloaded
        # ============================================================================>

        pandda.pickle_the_pandda(all=True)
        pandda.update_pandda_size(tag='After Pre-processing')

        # ============================================================================>
        #####
        # Write out Grid Point Distributions for interesting grid points (high modality, etc...)
        #####
        # ============================================================================>

        # TODO TODO TODO
        # TODO TODO TODO

        # ============================================================================>
        #####
        # Write Grid Point Distributions - Standard Function, just to provide an output
        #####
        # ============================================================================>

        try:
            from libtbx.math_utils import iceil
            grid_size = pandda.reference_grid().grid_size()
            num_points = min(10, min(grid_size))
            assert num_points != 0
            grid_points = zip(*[range(0, s, iceil(s/num_points)) for s in grid_size])
            pandda.write_grid_point_distributions(  grid_points     = grid_points,
                                                    output_filename = None          )
        except:
            print('UNIMPORTANT: FAILED TO WRITE AUTOMATIC DISTRIBUTION OF GRID POINTS')
            raise

        # ============================================================================>
        #####
        # Manual Settings
        #####
        # ============================================================================>

        # Minimum size of cluster
        min_cluster_volume   = pandda.params.blob_search.min_blob_volume
        # Z cutoff for maps
        z_cutoff             = pandda.params.blob_search.contour_level
        # Clustering Methods
        clustering_criterion = pandda.params.blob_search.clustering.criterion
        clustering_metric    = pandda.params.blob_search.clustering.metric
        clustering_method    = pandda.params.blob_search.clustering.linkage
        # Cutoff for separation of clusters (sqrt((2x)**2 + (2y)**2 + (2z)**2)) -- allows diagonal grid points to connect
        clustering_cutoff = 1.1 * numpy.math.sqrt(3) * pandda.reference_grid().grid_spacing()

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

        # TODO TODO TODO - MOVE?
        # Combine the symmetry mask and the atomic mask to create the mask for clustering
        pandda.log('===================================>>>', True)
        pandda.log('Combining Atomic Mask and Symmetry Mask', True)
        grid_idxr = pandda.reference_grid().grid_indexer()
        combined_point_mask = [gp for gp in pandda.reference_grid().global_mask().total_mask() if pandda.reference_grid().symmetry_mask().inner_mask_binary()[grid_idxr(gp)] == 0]

        pandda.log('GLOBAL TOTAL MASK: {!s}'.format(len(pandda.reference_grid().global_mask().total_mask())))
        pandda.log('SYMMETRY INNER MASK: {!s}'.format(len(pandda.reference_grid().symmetry_mask().inner_mask())))
        pandda.log('COMBINED MASK: {!s}'.format(len(combined_point_mask)))
        # TODO TODO TODO

        # Time the processing of the dataset maps
        t_start = time.time()

        for d_handler in pandda.datasets.mask(mask_name='selected for analysis'):

            if d_handler.raw_cluster_hits is None:

                pandda.log('===================================>>>', True)
                pandda.log('Calculating Z-MAPs for Dataset {!s} ({!s}/{!s})'.format(d_handler.d_tag, d_handler.d_num+1, pandda.datasets.size()), True)

                ##################################
                # EXTRACT MAP VALUES             #
                ##################################

                s_map = d_handler.morphed_map
                assert s_map is not None, 'NO MAP FOUND'

                # Write the sampled map
                if not os.path.exists(d_handler.output_handler.get_file('sampled_map')):
                    pandda.write_array_to_map(  output_file = d_handler.output_handler.get_file('sampled_map'),
                                                map_data    = s_map)

                # Write distribution of the map values
                if not os.path.exists(d_handler.output_handler.get_file('s_map_png')):
                    pandda.write_map_value_distribution(map_vals=s_map,
                                                        grid_indices=pandda.reference_grid().global_mask().outer_mask_indices(),
                                                        output_file = d_handler.output_handler.get_file('s_map_png'))

                ##################################
                # CALCULATE MEAN-DIFF MAPS       #
                ##################################

                d_map = d_handler.morphed_map - pandda.stat_maps.mean_map

                # Write the map
                if not os.path.exists(d_handler.output_handler.get_file('mean_diff_map')):
                    pandda.write_array_to_map(  output_file = d_handler.output_handler.get_file('mean_diff_map'),
                                                map_data    = d_map)

                # Write distribution of the map values
                if not os.path.exists(d_handler.output_handler.get_file('d_mean_map_png')):
                    pandda.write_map_value_distribution(map_vals=d_map,
                                                        grid_indices=pandda.reference_grid().global_mask().outer_mask_indices(),
                                                        output_file = d_handler.output_handler.get_file('d_mean_map_png'))

                ##################################
                # CALCULATE Z-MAPS               #
                ##################################

                ###################################################################
                # NAIVE Z-MAP - NOT USING UNCERTAINTY ESTIMATION OR ADJUSTED STDS #
                ###################################################################

                z_map_naive = pandda.calculate_z_map(map_vals=s_map, method='naive')

                # Write z map
                if not os.path.exists(d_handler.output_handler.get_file('z_map_naive')):
                    pandda.write_array_to_map(  output_file = d_handler.output_handler.get_file('z_map_naive'),
                                                map_data    = z_map_naive)

                # Write distribution of the map values
                if not os.path.exists(d_handler.output_handler.get_file('z_map_naive_png')):
                    pandda.write_map_value_distribution(map_vals = z_map_naive,
                                                        grid_indices = pandda.reference_grid().global_mask().outer_mask_indices(),
                                                        output_file = d_handler.output_handler.get_file('z_map_naive_png'),
                                                        plot_normal = True)

                ##################################

                # Normalise this map to N(0,1)
                z_map_naive_masked = [z_map_naive[i] for i in pandda.reference_grid().global_mask().outer_mask_indices()]
                z_map_naive_normalised = (z_map_naive - numpy.mean(z_map_naive_masked)) / numpy.std(z_map_naive_masked)

                # Write z map
                if not os.path.exists(d_handler.output_handler.get_file('z_map_naive_normalised')):
                    pandda.write_array_to_map(  output_file = d_handler.output_handler.get_file('z_map_naive_normalised'),
                                                map_data    = z_map_naive_normalised)

                # Write distribution of the map values
                if not os.path.exists(d_handler.output_handler.get_file('z_map_naive_normalised_png')):
                    pandda.write_map_value_distribution(map_vals = z_map_naive_normalised,
                                                        grid_indices = pandda.reference_grid().global_mask().outer_mask_indices(),
                                                        output_file = d_handler.output_handler.get_file('z_map_naive_normalised_png'),
                                                        plot_normal = True)

                ##################################
                # ADJUSTED+UNCERTAINTY Z-MAP     #
                ##################################

                z_map_corrected = pandda.calculate_z_map(map_vals=s_map, method='adjusted+uncertainty', map_uncertainty=d_handler.get_map_uncertainty())

                # Write z map
                if not os.path.exists(d_handler.output_handler.get_file('z_map_corrected')):
                    pandda.write_array_to_map(  output_file = d_handler.output_handler.get_file('z_map_corrected'),
                                                map_data    = z_map_corrected)

                # Write distribution of the map values
                if not os.path.exists(d_handler.output_handler.get_file('z_map_corrected_png')):
                    pandda.write_map_value_distribution(map_vals = z_map_corrected,
                                                        grid_indices = pandda.reference_grid().global_mask().outer_mask_indices(),
                                                        output_file = d_handler.output_handler.get_file('z_map_corrected_png'),
                                                        plot_normal = True)

                ##################################

                # Normalise this map to N(0,1)
                z_map_corrected_masked = [z_map_corrected[i] for i in pandda.reference_grid().global_mask().outer_mask_indices()]
                z_map_corrected_normalised = (z_map_corrected - numpy.mean(z_map_corrected_masked)) / numpy.std(z_map_corrected_masked)

                # Write z map
                if not os.path.exists(d_handler.output_handler.get_file('z_map_corrected_normalised')):
                    pandda.write_array_to_map(  output_file = d_handler.output_handler.get_file('z_map_corrected_normalised'),
                                                map_data    = z_map_corrected_normalised)

                # Write distribution of the map values
                if not os.path.exists(d_handler.output_handler.get_file('z_map_corrected_normalised_png')):
                    pandda.write_map_value_distribution(map_vals = z_map_corrected_normalised,
                                                        grid_indices = pandda.reference_grid().global_mask().outer_mask_indices(),
                                                        output_file = d_handler.output_handler.get_file('z_map_corrected_normalised_png'),
                                                        plot_normal = True)

                ##########################################
                # ANALYSE Z-MAP FOR STATISTICAL VALIDITY #
                ##########################################

                z_map_stats = basic_statistics(flex.double([z_map_corrected[i] for i in pandda.reference_grid().global_mask().outer_mask_indices()]))

                d_handler.z_map_stats['z_map_mean']     = z_map_stats.mean
                d_handler.z_map_stats['z_map_std']      = z_map_stats.bias_corrected_standard_deviation
                d_handler.z_map_stats['z_map_skew']     = z_map_stats.skew
                d_handler.z_map_stats['z_map_kurtosis'] = z_map_stats.kurtosis

                #################################################
                # XXX WHICH MAP TO DO THE BLOB SEARCHING ON XXX #
                #################################################

                z_map = z_map_corrected_normalised

                ##################################
                # LOOK FOR CLUSTERS OF Z-SCORES  #
                ##################################

                d_handler.raw_cluster_hits = pandda.cluster_high_z_values(
                                                d_handler            = d_handler,
                                                z_map                = z_map,
                                                z_cutoff             = z_cutoff,
                                                point_mask           = combined_point_mask,
                                                min_cluster_volume   = min_cluster_volume,
                                                clustering_cutoff    = clustering_cutoff,
                                                clustering_criterion = clustering_criterion,
                                                clustering_metric    = clustering_metric,
                                                clustering_method    = clustering_method
                                            )

            # Create cluster object from the clustered points
            if d_handler.raw_cluster_hits:
                d_handler.clustered_hits = cluster_data(d_handler.raw_cluster_hits)

            ################################################
            # AUTOGENERATE SCRIPTS FOR VIEWING THE DATASET #
            ################################################

            pandda.write_pymol_scripts(d_handler=d_handler)

            if d_handler.raw_cluster_hits:

                pandda.log('===================================>>>', True)
                pandda.log('Imaging blobs in Dataset {!s}'.format(d_handler.d_tag), True)

                sorted_blob_indices = d_handler.clustered_hits.sort(sorting_function=max)

                for blob_rank, blob_grid_peak in enumerate(d_handler.clustered_hits.get_centroids(indices=sorted_blob_indices)):

                    blob_cart_peak = [g*pandda.reference_grid().grid_spacing() for g in blob_grid_peak]

                    # Make images of the blob
                    pandda.image_blob(  script    = d_handler.output_handler.get_file('ccp4mg_script'),
                                        image     = d_handler.output_handler.get_file('ccp4mg_png'),
                                        d_handler = d_handler,
                                        point_no  = blob_rank+1,
                                        point     = blob_cart_peak,
                                        towards   = [10,10,10]
                                    )

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

        # ============================================================================>

        t_end = time.time()
        pandda.log('===================================>>>', True)
        pandda.log('Map Processing Time: {!s}'.format(time.strftime("%H hours:%M minutes:%S seconds", time.gmtime(t_end - t_start))), True)

        pandda.pickle_the_pandda(components=['datasets'])
        pandda.update_pandda_size(tag='After Processing')

        # ============================================================================>
        #####
        # ANALYSIS
        #####
        # Analyse the processed data
        # TODO Create a `hit processor` class
        # ============================================================================>

        cluster_total, cluster_num, all_dataset_clusters = pandda.collate_all_clusters()

        # Print a summary of the number of identified clusters
        pandda.log('===================================>>>', True)
        pandda.log('Total Datasets with Clusters: {!s}'.format(len(cluster_num)), True)
        pandda.log('Total Clusters: {!s}'.format(cluster_total), True)
        # Print some slightly less important information
        pandda.log('===================================>>>', False)
        for d_tag, cluster_count in cluster_num:
            pandda.log('Dataset {!s}: {!s} Clusters'.format(d_tag, cluster_count), False)

        combined_cluster = pandda.process_z_value_clusters()

        # Write a combined file of all of the result hits
        pandda.write_ranked_cluster_csv(
                                cluster=combined_cluster,
                                sorted_indices=combined_cluster.sort(sorting_data='values', sorting_function=max, decreasing=True),
                                outfile=pandda.output_handler.get_file('blob_summaries')
                                )

        # ============================================================================>
        #####
        # SUMMARIES ------------------------------>>>
        #####
        # ============================================================================>

        pandda.collect_map_statistics()

        pandda.write_analysis_summary()
        pandda.write_summary_csvs()
        pandda.write_summary_graphs()

        if pandda.args.settings.testing:
            test_pandda_results(pandda)

        pandda.exit()

    except KeyboardInterrupt:
        raise
    except AssertionError:
        raise
#    except Exception as err:
#        return (pandda, err)

    return pandda, None

def test_pandda_results(pandda):
    """If a list of hits is available, test to see whether the pandda identified them"""

    from Giant.Stats.Utils import calculate_roc_vals

    # Cutoff for how close the "discovered" points have to be to the centroid of the ligand to be correct
    max_peak_dist = 5

    try:

        # ============================================================================>
        #####
        # Extract list of correct hits
        #####
        # ============================================================================>

        solution_file = os.path.join(pandda.output_handler.get_dir('root'), '../correct_hits.summary')
        if not os.path.exists(solution_file):
            print('!!!        NOT PERFORMING TESTING       !!!')
            print('!!! CORRECT LIGANDS FILE DOES NOT EXIST !!!')
            print('!!!        NOT PERFORMING TESTING       !!!')
            return

        solutions = [line.strip().split(', ') for line in open(solution_file, 'r').readlines()]
        correct_ligands = {}
        [correct_ligands.setdefault(xtal,[]) for xtal in [d.d_tag for d in pandda.datasets.all()]]
        [correct_ligands[xtal].append([(float(x), float(y), float(z)), smile]) for xtal, model, smile, x, y, z in solutions if xtal in correct_ligands.keys()]

        # ============================================================================>

        combined_cluster = pandda.process_z_value_clusters()

        # Iterate through the correct ligands and check that there is a cluster near them
        pandda.log('===================================>>>')
        pandda.log('TRAINING - Checking Hits have been identified')

        missed_datasets = 0
        missed_ligands = 0
        missed_glycols = 0
        minimum_dists = []

        # Go through the datasets with ligands in them
        for d_tag in sorted(correct_ligands.keys()):
            pandda.log('Checking for the ligands in Dataset {!s}'.format(d_tag))

            # Get the list of ligands that we should be detecting
            dataset_ligands = correct_ligands[d_tag]

            # Get the clusters of identified points for this dataset
            clust = pandda.datasets.get(d_tag=d_tag).clustered_hits

            # Go through the ligands in the dataset
            for loc_cart_tuple, smile in dataset_ligands:

                loc_grid_tuple = tuple(flex.double(loc_cart_tuple)/pandda.reference_grid().grid_spacing())
                pandda.log('\tLooking for {!s} at {!s}'.format(smile, loc_grid_tuple))

                # Check to see if there are any identified points in this dataset
                if clust==None:
                    pandda.log("\t\tTHIS DATASET WASN'T IDENTIFIED AS INTERESTING - MISSED IT")
                    missed_datasets += 1
                    continue

                # Calculate the distances between ligands and the identified clusters
                dists = []
                found_it = False

                # Iterate through the clusters in this dataset
                for clust_points in clust.get_points():
                    # Iterate through the points in this cluster
                    for grid_point in clust_points:
                        dists.append((flex.double(grid_point) - flex.double(loc_grid_tuple)).norm())
                        #print '\t\tHow about {!s}?'.format(grid_point)
                        if dists[-1] <  max_peak_dist/pandda.reference_grid().grid_spacing():
                            pandda.log('\t\t\tFOUND IT!')
                            found_it = True
                            break
                    if found_it == True:
                        break

                # Check to see if the ligands has been found
                if not found_it:
                    pandda.log('\t\t\tMISSED IT...')
                    if smile == 'OCCO':
                        missed_glycols += 1
                    else:
                        missed_ligands += 1

                if smile == 'OCCO':
                    pass
                else:
                    minimum_dists.append(min(dists)*pandda.reference_grid().grid_spacing())

        pandda.log('----------------------------------->>>')
        pandda.log('MISSED DATASETS: {!s}'.format(missed_datasets))
        pandda.log('MISSED LIGANDS: {!s}'.format(missed_ligands))
        pandda.log('MISSED GLYCOLS: {!s}'.format(missed_glycols))
        pandda.log('----------------------------------->>>')
        pandda.log('MINIMUM DISTANCES TO LIGAND (A): \n{!s}'.format(map(int,minimum_dists)))
        pandda.log('----------------------------------->>>')

        pandda.log('===================================>>>')
        pandda.log('Sorting Clusters')

        # Sort the clusters by size
        size_sorted_indices = combined_cluster.sort(sorting_data='sizes', sorting_function=None, decreasing=True)
        # Sort the clusters by mean
        mean_sorted_indices = combined_cluster.sort(sorting_data='values', sorting_function=numpy.mean, decreasing=True)
        # Sort the clusters by max
        max_sorted_indices = combined_cluster.sort(sorting_data='values', sorting_function=max, decreasing=True)

        pandda.log('===================================>>>')
        pandda.log('TRAINING - Calculating ROC Curves')

        COMBINED_ROC_RESULTS = {}

        for test_num, (c_indices, rank_vals) in enumerate([(size_sorted_indices, combined_cluster.get_sizes(size_sorted_indices)),
                                                           (mean_sorted_indices, combined_cluster.get_means(mean_sorted_indices)),
                                                           (max_sorted_indices,  combined_cluster.get_maxima(max_sorted_indices))]):

            test_type = ['c_size','z_mean','z_peak'][test_num]

#           if test_type != 'z_peak':
#               print('Skipping test: {!s}'.format(test_type))
#               continue

            # Output ROC information
            ROC_OUTPUT = []

            pandda.log('===================================>>>')
            pandda.log('Ranking by: {!s}'.format(test_type))
            pandda.log('===================================>>>')

            # List to see when a ligand is identified - avoids multiple identifications
            identified_ligands = []
            identified_idx = 1

            # Calculate ROC Curves etc - Iterate through clusters and check if near ligand
            for c_rank, c_index in enumerate(c_indices):
                # Pull out information to identify the cluster
                d_tag, c_num = combined_cluster.get_keys([c_index])[0]
                #c_centroid = combined_cluster.get_centroids([c_index])[0]
                pandda.log('Checking Dataset {:4}, Cluster {:4}, Rank: {:4}, Val: {:6.2f}'.format(d_tag, c_num, c_rank, rank_vals[c_rank]))

                # Reset
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

                        # Iterate through and see if this overlaps with the ligand
                        for c_point in combined_cluster.get_points([c_index])[0]:
                            # Calculate distance to ligand
                            dists.append((flex.double(loc_grid_tuple) - flex.double(c_point)).norm())
                            # Check if correct
                            if dists[-1] < max_peak_dist/pandda.reference_grid().grid_spacing():
                                is_a_ligand = 1
                                break

                        # If it's a ligand, check to see if it's been previously identified
                        if is_a_ligand:
                            if (d_tag, l_num) in identified_ligands:
                                pandda.log('\t. - Ligand already found: Dataset {!s}'.format(d_tag))
                                skip = True
                            else:
                                pandda.log("\t{!s} - Ligand found: Dataset {!s}".format(identified_idx, d_tag))
                                identified_idx += 1
                                identified_ligands.append((d_tag, l_num))

                if not skip:
                    ROC_OUTPUT.append({'d_tag': d_tag, 'c_num': c_num, 'is_a_ligand': is_a_ligand, 'rank': c_rank, 'val': rank_vals[c_rank]})

            # Calculate the ROC Curve
            correct_class = [t['is_a_ligand'] for t in ROC_OUTPUT]
            rank_vals     = [t['val'] for t in ROC_OUTPUT]
            ROC_RESULTS = calculate_roc_vals(correct_class=correct_class, rank_vals=rank_vals)
            COMBINED_ROC_RESULTS[test_type] = ROC_RESULTS

            pandda.log('Ligands Identified: {!s}'.format([l[0] for l in identified_ligands]))

        # Load plotting
        try:
            import matplotlib
            # Setup so that we can write without a display connected
            matplotlib.interactive(0)
            default_backend, validate_function = matplotlib.defaultParams['backend']
            from matplotlib import pyplot
            output_graphs = True
        except:
            output_graphs = False

        # OUTPUT ROC CURVE
        if output_graphs:

            for test_type in COMBINED_ROC_RESULTS.keys():

                roc = COMBINED_ROC_RESULTS[test_type]

                fig = pyplot.figure()
                pyplot.title('ROC CURVE FOR {!s}'.format(test_type))
                pyplot.plot(roc['VALS'], roc['SENS'], 'ko-')
                pyplot.plot(roc['VALS'], roc['SPEC'], 'ro-')
                pyplot.plot(roc['VALS'], roc['PREC'], 'go-')
                pyplot.plot(roc['VALS'], roc['ACC'],  'bo-')
                pyplot.xlabel('SCORE CUTOFF ({!s})'.format(test_type))
                pyplot.ylabel('ROC RESPONSE')
                # Apply tight layout to prevent overlaps
                pyplot.tight_layout()
                # Save
                pyplot.savefig(os.path.join(pandda.output_handler.get_dir('analyses'), 'roc-cutoff-{!s}.png'.format(test_type)))
                pyplot.close(fig)

                fig = pyplot.figure()
                pyplot.title('ROC CURVE FOR {!s}'.format(test_type))
                pyplot.plot(roc['FPR'], roc['TPR'],  'go-')
                pyplot.xlabel('FPR')
                pyplot.ylabel('TPR')
                # Apply tight layout to prevent overlaps
                pyplot.tight_layout()
                # Save
                pyplot.savefig(os.path.join(pandda.output_handler.get_dir('analyses'), 'roc-tpr-fpr-{!s}.png'.format(test_type)))
                pyplot.close(fig)

        # ============================================================================>
        # ======================================>
        # Manual Analyses - These only need processing once
        # ======================================>
        # ============================================================================>
        analyses_dir = pandda.output_handler.get_dir('analyses')
        # ============================================================================>

        if 0:
            pandda.log('===================================>>>')
            pandda.log('Calculating Deviations of C-alphas between structures')

            rms = lambda vals: numpy.sqrt(numpy.mean(numpy.abs(vals)**2))
            norm = lambda vals: numpy.sqrt(numpy.sum(numpy.abs(vals)**2))

            # Pull all c-alpha sites for each structure
            all_sites = numpy.array([d.transform_points_to_reference(d.get_calpha_sites()) for d in pandda.datasets.all()])
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

            pandda.log('Largest deviation from the mean site: {!s}'.format(diff_norms.max()))
            pandda.log('Average deviation from the mean site: {!s}'.format(diff_norms.mean()))

        # ============================================================================>

        if 0:
            pandda.log('===================================>>>')
            pandda.log('Clustering the Refined Structures')

            distance_matrix = []
            for d1 in pandda.datasets.all():
               distance_matrix.append([d1.transform_points_to_reference(d1.get_calpha_sites()).rms_difference(d2.transform_points_to_reference(d2.get_calpha_sites())) for d2 in pandda.datasets.all()])

            distance_matrix = numpy.array(distance_matrix)

            with open(os.path.join(analyses_dir,'calpha_distance_matrix.csv'), 'w') as fh:
                for row in distance_matrix:
                    out_list = row.round(3).tolist()
                    out_line = ', '.join(map(str,out_list)) + '\n'
                    fh.write(out_line)

        # ============================================================================>

    except:
        pandda.log('FAILURE DURING TESTING', True)
        raise

    return 0

# ============================================================================>
#
#   GENERATE ARGS
#
# ============================================================================>

def process_args_2(args=None):

    from PANDDAs.Parser import build_pandda_parser

    pandda_arg_parser = build_pandda_parser()

    if args:
        pandda_args = pandda_arg_parser.parse_args(args)
    else:
        pandda_args = pandda_arg_parser.parse_args()

    return pandda_args

# ============================================================================>
#
#   COMMAND LINE RUN
#
# ============================================================================>

if __name__ == '__main__':

    pandda, err = pandda_main(args=None)

    if err:
        raise err

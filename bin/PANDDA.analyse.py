#!/home/npearce/bin/Developer/cctbx/cctbx_build/bin/cctbx.python

import os, sys, glob, time
import numpy

from scitbx.array_family import flex
from scitbx.math import basic_statistics

from Giant.Xray.Symmetry import combine_hierarchies, generate_adjacent_symmetry_copies

from PANDDAs.Parser import build_pandda_parser
from PANDDAs.Main import multi_dataset_analyser, spherical_mask, atomic_mask

def pandda_main(args=None):

    # ============================================================================>
    #####
    # MANUAL SETTINGS
    #####
    # ============================================================================>

    testing = True

    # ============================================================================>
    #####
    # Parse Input Args
    #####
    # ============================================================================>

    pandda_arg_parser = build_pandda_parser()
    if testing:
        pandda_args = pandda_arg_parser.parse_args(['--data-dirs','/home/npearce/Analyses/FragmentSoaks/BAZ2BA-Feb13-Soak/ProcessedFragmentSoak/BAZ2BA-*/1-apo',
                                                    '--outdir','/home/npearce/Analyses/FragmentSoaks/BAZ2BA-Feb13-Soak/pandda',
                                                    '--pdb-style','apo-BAZ2BA-*-refmac.pdb',
                                                    '--ref-pdb','/home/npearce/Analyses/FragmentSoaks/BAZ2BA-Feb13-Soak/reference.pdb',
                                                    '--ref-mtz','/home/npearce/Analyses/FragmentSoaks/BAZ2BA-Feb13-Soak/reference.mtz',
                                                    '--cpus', '6', '--verbose'])
    elif args:
        pandda_args = pandda_arg_parser.parse_args(args)
    else:
        pandda_args = pandda_arg_parser.parse_args()

    print('===================================>>>')
    print('INPUT SETTINGS')
    print('===================================>>>')

    for key in pandda_args.__dict__.keys():
        print '{!s:>30}\t\t{!s}'.format(key, pandda_args.__dict__[key])

    # ============================================================================>
    #####
    # Initialise
    #####
    # ============================================================================>

    pandda = multi_dataset_analyser(pandda_args)

    # ============================================================================>
    #####
    # Initialise Settings
    #####
    # ============================================================================>

    pandda.set_alignment_method(method='global')

    #pandda.set_obs_map_type(map_type='2mFo-DFc')
    pandda.set_obs_map_type(map_type='mFo')
    pandda.set_map_scaling(scaling='volume')
    pandda.set_res_factor(0.4)
    pandda.set_border_padding(5)

    # XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX
    # MOVE THIS INTO INIT FUNCTION AND POPULATE FROM THE ARG PARSER
    pandda.set_low_resolution(2.5)
    pandda.set_high_resolution(0)
    pandda.set_cut_resolution(2)
    # XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX

    pandda.run_pandda_init()

    # ============================================================================>
    #####
    # Build list of files in data directories
    #####
    # ============================================================================>

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
        # TODO XXX ADD `data_scaling` METHOD HERE (ALLOW NONE) XXX TODO
        pandda.scale_datasets()
        pandda.align_datasets(method=pandda.get_alignment_method())
        pandda.generate_crystal_contacts()

    else:
        # Rebuild the masks of the rejected datasets (quick)
        pandda.initialise_analysis()
        pandda.filter_datasets_1()

    # XXX XXX XXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXX XXX XXX
    # Use symmetry operations to create the symmetry mates of the reference structure
    sym_ops, sym_op_contacts, sym_hierarchies, chain_mappings = generate_adjacent_symmetry_copies(ref_hierarchy=pandda.reference_dataset().new_structure().hierarchy,
                                                                                                  crystal_symmetry=pandda.reference_dataset().get_input().crystal_symmetry(),
                                                                                                  buffer_thickness=pandda.get_border_padding())
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
    #pandda.collate_dataset_parameters()
    pandda.collect_dataset_variation_data()

    # Filter out the datasets that are not isomorphous and therefore incomparable
    pandda.filter_datasets_2()

    # Analyses the crystallographic and structural variability of the datasets
    pandda.calculate_mean_structure_and_protein_masks(deviation_cutoff=0.5)

    # Analyse the structural variation in the datasets
    # Rename? - This should be the only thing that populates pandda.get_residue_observations object - simple class?
    #pandda.collate_residue_parameters()
    pandda.collect_structure_variation_data()

    # ============================================================================>
    #####
    # Update Settings
    #####
    # ============================================================================>

    # Update the resolution limits using the resolution limits from the datasets supplied
    pandda.set_low_resolution(   min(   pandda.get_low_resolution(),
                                        max(pandda.get_dataset_observations().get_data('high_res_limit'))   )   )
    pandda.set_high_resolution(  max(   pandda.get_high_resolution(),
                                        min(pandda.get_dataset_observations().get_data('high_res_limit'))   )   )

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
        pandda.create_reference_grid(grid_spacing=pandda.get_res_factor()*pandda.get_cut_resolution(), expand_to_origin=False, buffer=pandda.get_border_padding()+0.1)

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
        global_mask = atomic_mask(  cart_sites = pro_sites_cart,
                                    grid_size  = pandda.reference_grid().grid_size(),
                                    unit_cell  = pandda.reference_grid().fake_unit_cell(),
                                    max_dist   = pandda.get_border_padding(),
                                    min_dist   = 1.8 )
        pandda.reference_grid().set_global_mask(global_mask)

    if pandda.reference_grid().symmetry_mask() is None:
        pandda.log('===================================>>>')
        pandda.log('Generating Symmetry Mask')
        # Pull out the cartesian sites of the symmetry mates
        cache = symmetry_root.atom_selection_cache()
        sym_sites_cart = symmetry_root.select(cache.selection('pepnames and not element H')).atoms().extract_xyz()
        # Generate the symmetry mask
        symmetry_mask = atomic_mask(cart_sites = sym_sites_cart,
                                    grid_size  = pandda.reference_grid().grid_size(),
                                    unit_cell  = pandda.reference_grid().fake_unit_cell(),
                                    max_dist   = pandda.get_border_padding(),
                                    min_dist   = 1.8 )
        pandda.reference_grid().set_symmetry_mask(symmetry_mask)

    # Print the summaries
    pandda.log('===================================>>>')
    pandda.log('Grid Summary: ')
    pandda.log(pandda.reference_grid().summary())
    pandda.log(pandda.reference_grid().local_mask().summary())
    pandda.log(pandda.reference_grid().global_mask().summary())

    # TODO TODO TODO INCORPORATE THE COMBINATION OF MASKS INTO THIS FUNCTION
    # Create various masks to define regions of the grid by distance to the protein and symmetry copies
    pandda.mask_reference_grid(d_handler=pandda.reference_dataset())

#    if pandda.reference_grid().masked_grid_points() is None:
#        pandda.mask_resampled_reference_grid()
#    # TODO TODO TODO

    # ============================================================================>
    #####
    # Local Map Handler Objects and Extract Map Values
    #####
    # ============================================================================>

    pandda.load_reference_map_handler()
    pandda.extract_reference_map_values()

    # TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
    # TODO  SAVE AND LOAD MAPS FROM FILES (WOULD SAVE SOOOOO MUCH TIME!!!)  TODO
    # TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO

    pandda.load_all_map_handlers()

    # ============================================================================>
    #####
    # Scale the Maps, Calculate Statistical Maps
    #####
    # ============================================================================>

    if pandda.is_new_pandda():
        pandda.scale_raw_maps()
        pandda.calculate_mean_map()
        pandda.calculate_map_uncertainties()
        pandda.calculate_map_statistics()

    # TODO TODO TODO
    # TODO IMPLEMENT - Collect all of the map statistics from the different dataset_handlers
    # So that the data_collections can be populated independently of when the uncertainties etc were calculated
    #pandda.collect_map_data()
    # TODO TODO TODO

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
                                                output_filename=None )
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
    # Manual Settings
    #####
    # ============================================================================>

    # Minimum size of cluster
    min_cluster_volume = 10
    # Z cutoff for maps
    z_cutoff = 2
    # Clustering Methods
    clustering_criterion='distance'
    clustering_metric='euclidean'
    clustering_method='single'
    # Cutoff for separation of clusters (sqrt((2x)**2 + (2y)**2 + (2z)**2))
    clustering_cutoff = 1.1 * numpy.math.sqrt(3) * pandda.reference_grid().grid_spacing()

    # Cutoff for how close the "discovered" points have to be to the centroid of the ligand to be correct
    max_peak_dist = 5

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
        pandda.log('Calculating Z-MAPs for Dataset {!s} ({!s}/{!s})'.format(d_handler.d_tag, d_handler.d_num+1, pandda.get_number_of_datasets()), True)

        ##################################
        # EXTRACT MAP VALUES             #
        ##################################

        s_map = pandda.get_map(d_handler=d_handler, map_type=pandda.get_obs_map_type())

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
        # EXTRACT DIFFERENCE MAP         #
        ##################################

        #diff_map = pandda.get_map(d_handler=d_handler, map_type=pandda.get_diff_map_type())
        #
        ## Write the sampled map
        #pandda.write_array_to_map(  output_file = d_handler.output_handler.get_file('difference_map'),
        #                            map_data    = diff_map)
        #
        ## Write distribution of the map values
        #pandda.write_map_value_distribution(map_vals=diff_map,
        #                                    grid_indices=pandda.reference_grid().global_mask().inner_mask_indices(),
        #                                    output_file = d_handler.output_handler.get_file('d_map_png'))

        ##################################
        # CALCULATE MEAN-DIFF MAPS       #
        ##################################

        d_map = pandda.get_map(d_handler=d_handler, map_type=pandda.get_obs_map_type()) - pandda.get_mean_map()

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

        pandda.cluster_high_z_values(  d_handler            = d_handler,
                                       z_map                = z_map,
                                       z_cutoff             = z_cutoff,
                                       point_mask           = pandda.reference_grid().global_mask().total_mask(),
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

    # ============================================================================>

    pandda.pickle_the_pandda()
    pandda.update_pandda_size(tag='After Processing')

    # ============================================================================>
    #####
    # ANALYSIS
    #####
    # Analyse the processed data
    # TODO Create a `hit processor` class
    # ============================================================================>

    cluster_total, cluster_num, all_dataset_clusters = pandda.collate_all_clusters()

    combined_cluster = pandda.process_z_value_clusters()

    # ============================================================================>
    #####
    # SUMMARIES ------------------------------>>>
    #####
    # ============================================================================>

    pandda.collect_map_statistics()

    pandda.write_analysis_summary()
    pandda.write_pymol_scripts()
    pandda.write_summary_csvs()
    pandda.write_summary_graphs()

    # ============================================================================>
    #####
    # XXX XXX XXX Setting Analysis Variables XXX XXX XXX
    #####
    # ============================================================================>


    solution_file = os.path.join(pandda.output_handler.get_dir('root'), 'correct_hits.summary')
    if not os.path.exists(solution_file):
        testing = False
        print('CORRECT LIGANDS FILE DOES NOT EXIST!!!')
    else:
        testing = True
        solutions = [line.strip().split(', ') for line in open(solution_file, 'r').readlines()]

        correct_ligands = {}
        [correct_ligands.setdefault(xtal,[]) for xtal in [d.d_tag for d in pandda.get_all_datasets()]]
        [correct_ligands[xtal].append([(int(x), int(y), int(z)), smile]) for xtal, model, smile, x, y, z in solutions if xtal in correct_ligands.keys()]

    # XXX XXX XXX -------------------------- XXX XXX XXX

    # ============================================================================>

    if testing:

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
            clust = pandda.get_dataset(d_tag=d_tag).clustered_hits

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

    if testing:

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

        # Output ROC information
        roc_testing = []

        for test_num, (c_indices, rank_vals) in enumerate([(size_sorted_indices, combined_cluster.get_sizes(size_sorted_indices)),
                                                           (mean_sorted_indices, combined_cluster.get_means(mean_sorted_indices)),
                                                           (max_sorted_indices,  combined_cluster.get_maxima(max_sorted_indices))]):

            test_type = ['c_size','z_mean','z_peak'][test_num]

            pandda.log('===================================>>>')
            pandda.log('Ranking by: {!s}'.format(test_type))
            pandda.log('===================================>>>')

            # List to see when a ligand is identified - avoids multiple identifications
            identified_ligands = []

            # Calculate ROC Curves etc - Iterate through clusters and check if near ligand
            for c_rank, c_index in enumerate(c_indices):
                # Pull out information to identify the cluster
                d_tag, c_num = combined_cluster.get_keys([c_index])[0]
                #c_centroid = combined_cluster.get_centroids([c_index])[0]
                pandda.log('Checking Dataset {:4}, Cluster {:4}, Rank: {:4}, Val: {:6.2f}'.format(d_tag, c_num, c_rank, rank_vals[c_rank]))

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
                                pandda.log("\t1 - Ligand found: Dataset {!s}".format(d_tag))
                                identified_ligands.append((d_tag, l_num))

                if not skip:
                    # TODO WRITE OUT THIS INFORMATION TODO
                    roc_testing.append([test_type, d_tag, d_tag, c_num, is_a_ligand, c_rank, rank_vals[c_rank]])

            pandda.log('Ligands Identified: {!s}'.format([l[0] for l in identified_ligands]))

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

        pandda.log('Largest deviation from the mean site: {!s}'.format(diff_norms.max()))
        pandda.log('Average deviation from the mean site: {!s}'.format(diff_norms.mean()))

    # ============================================================================>

    if 0:
        pandda.log('===================================>>>')
        pandda.log('Clustering the Refined Structures')

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

    return pandda




# ============================================================================>
#
#   COMMAND LINE RUN
#
# ============================================================================>

if __name__ == '__main__':

    pandda = pandda_main()


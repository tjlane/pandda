import os, sys, glob, time, gc

#################################
try:
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.interactive(0)
    from matplotlib import pyplot
    pyplot.style.use('ggplot')
except:
    pass
#################################

import numpy

from scitbx.array_family import flex
from scitbx.math import basic_statistics

from libtbx import easy_mp
from libtbx.utils import Sorry, Failure

from bamboo.common import Meta, Info
from bamboo.common.logs import Log
from bamboo.common.path import rel_symlink

from giant.jiffies import extract_params_default
from giant.grid.masks import grid_mask
from giant.xray.maps.bdc import calculate_bdc_subtracted_map, estimate_event_background_correction

from pandda import welcome
from pandda import analyse_html, analyse_graphs
from pandda.phil import pandda_phil
from pandda.analyse_main import PanddaMultiDatasetAnalyser, PanddaMapAnalyser, PanddaZMapAnalyser, MapHolder
from pandda.handlers import DatasetHandler
from pandda.events import PointCluster, Event
from pandda.misc import *

# ============================================================================>
#
###                 Dataset Processing Function
#
# ============================================================================>

def process_dataset_map_func(proc_args):
    """The main function for comparing each dataset-map to the ensemble-maps"""

    # ============================================================================>
    # Unpack input args
    # ============================================================================>
    d_handler              = proc_args[0]
    m_handler              = proc_args[1]
    reference_grid         = proc_args[2]
    reference_dataset      = proc_args[3]
    map_analyser           = proc_args[4]
    args, params, settings = proc_args[5]

    log_strs = []

    # ============================================================================>
    # Build new blob search object
    # ============================================================================>
    print('Analysing Dataset {!s}'.format(d_handler.tag))

    blob_finder = PanddaZMapAnalyser( params       = params.blob_search,
                                      grid_spacing = reference_grid.grid_spacing(),
                                      log          = Log(log_file=d_handler.output_handler.get_file('dataset_log'), verbose=False) )

    # ============================================================================>
    # Generate masks for this dataset
    # ============================================================================>
    print('Generating Symmetry Contacts for Dataset {!s}'.format(d_handler.tag))

    # Generate symmetry contacts for this dataset
    dataset_sym_copies = d_handler.generate_symmetry_copies(rt_method      = params.alignment.method,
                                                            save_operators = False,
                                                            buffer         = params.masks.outer_mask+5   )
    # Extract protein atoms from the symmetry copies
    cache = dataset_sym_copies.atom_selection_cache()
    dataset_sym_sites_cart = dataset_sym_copies.select(cache.selection('pepnames and not element H')).atoms().extract_xyz()
    # Generate custom grid mask for this dataset
    dataset_contact_mask = grid_mask(   cart_sites = dataset_sym_sites_cart,
                                        grid_size  = reference_grid.grid_size(),
                                        unit_cell  = reference_grid.unit_cell(),
                                        max_dist   = params.masks.outer_mask,
                                        min_dist   = params.masks.inner_mask_symmetry )
    # Combine the standard mask with the custom mask
    grid_idxr = reference_grid.grid_indexer()
    dataset_total_mask = [gp for gp in reference_grid.global_mask().total_mask() if dataset_contact_mask.inner_mask_binary()[grid_idxr(gp)] < 0.5]
    dataset_outer_mask = [gp for gp in reference_grid.global_mask().outer_mask() if dataset_contact_mask.inner_mask_binary()[grid_idxr(gp)] < 0.5]
    # Write coordinates of symmetry copies
    dataset_sym_copies.write_pdb_file(d_handler.output_handler.get_file('symmetry_copies'))

    # ============================================================================>
    #####
    # CALCULATE Z-MAPS AND LOOK FOR LARGE BLOBS
    #####
    # ============================================================================>
    print('Calculating Z-MAPs for Dataset {!s}'.format(d_handler.tag))

    # ============================================================================>
    # CALCULATE MEAN-DIFF MAPS
    # ============================================================================>
    assert m_handler.map is not None, 'NO MAP FOUND'
    d_map = m_handler.map - map_analyser.statistical_maps.mean_map

    # ============================================================================>
    # NAIVE Z-MAP - NOT USING UNCERTAINTY ESTIMATION OR ADJUSTED STDS
    # ============================================================================>
    z_map_naive = map_analyser.calculate_z_map(map=m_handler.map, method='naive')
    # Normalise this map to N(0,1)
    z_map_naive_masked     = [z_map_naive[i] for i in reference_grid.global_mask().outer_mask_indices()]
    z_map_naive_normalised = (z_map_naive - numpy.mean(z_map_naive_masked)) / numpy.std(z_map_naive_masked)

    # ============================================================================>
    # UNCERTAINTY Z-MAP - NOT USING ADJUSTED STDS
    # ============================================================================>
    z_map_uncty = map_analyser.calculate_z_map(map=m_handler.map, uncertainty=m_handler.meta.map_uncertainty, method='uncertainty')
    # Normalise this map to N(0,1)
    z_map_uncty_masked     = [z_map_uncty[i] for i in reference_grid.global_mask().outer_mask_indices()]
    z_map_uncty_normalised = (z_map_uncty - numpy.mean(z_map_uncty_masked)) / numpy.std(z_map_uncty_masked)

    # ============================================================================>
    # ADJUSTED+UNCERTAINTY Z-MAP
    # ============================================================================>
    z_map_compl = map_analyser.calculate_z_map(map=m_handler.map, uncertainty=m_handler.meta.map_uncertainty, method='adjusted+uncertainty')
    # Normalise this map to N(0,1)
    z_map_compl_masked     = [z_map_compl[i] for i in reference_grid.global_mask().outer_mask_indices()]
    z_map_compl_normalised = (z_map_compl - numpy.mean(z_map_compl_masked)) / numpy.std(z_map_compl_masked)

    # ============================================================================>
    # ANALYSE Z-MAP FOR STATISTICAL VALIDITY
    # ============================================================================>
    # Calculate statistics of z-maps
    z_map_stats = basic_statistics(flex.double(z_map_compl_masked))
    m_handler.meta.z_mean = z_map_stats.mean
    m_handler.meta.z_stdv = z_map_stats.bias_corrected_standard_deviation
    m_handler.meta.z_skew = z_map_stats.skew
    m_handler.meta.z_kurt = z_map_stats.kurtosis

    # ============================================================================>
    # SELECT WHICH MAP TO DO THE BLOB SEARCHING ON
    # ============================================================================>
    if params.z_map.map_type   == 'naive':
        z_map = z_map_naive_normalised
    elif params.z_map.map_type == 'adjusted':
        z_map = z_map_adjst_normalised
    elif params.z_map.map_type == 'uncertainty':
        z_map = z_map_uncty_normalised
    elif params.z_map.map_type == 'adjusted+uncertainty':
        z_map = z_map_compl_normalised

    # ============================================================================>
    # STORE THE Z MAP
    # ============================================================================>
    z_map_holder = MapHolder(   num         = m_handler.num,
                                tag         = m_handler.tag,
                                map         = z_map,
                                # Change these for the 'fake' grid unit_cell and a P1 space_group
                                unit_cell   = reference_grid.unit_cell(),
                                space_group = reference_grid.space_group(),
                                meta        = Meta({'type'          : 'z-map',
                                                    'resolution'    : m_handler.meta.resolution}),
                                parent      = m_handler  )

    # ============================================================================>
    #####
    # WRITE ALL MAP DISTRIBUTIONS (THESE DON'T USE MUCH SPACE)
    #####
    # ============================================================================>
    # Sampled Map
    write_map_value_distribution(map_vals     = m_handler.map.select(reference_grid.global_mask().outer_mask_indices()),
                                 output_file  = d_handler.output_handler.get_file('s_map_png'))
    # Mean-Difference
    write_map_value_distribution(map_vals     = d_map.select(reference_grid.global_mask().outer_mask_indices()),
                                 output_file  = d_handler.output_handler.get_file('d_mean_map_png'))
    # Naive Z-Map
    write_map_value_distribution(map_vals     = z_map_naive.select(reference_grid.global_mask().outer_mask_indices()),
                                 output_file  = d_handler.output_handler.get_file('z_map_naive_png'),
                                 plot_normal  = True)
    # Normalised Naive Z-Map
    write_map_value_distribution(map_vals     = z_map_naive_normalised.select(reference_grid.global_mask().outer_mask_indices()),
                                 output_file  = d_handler.output_handler.get_file('z_map_naive_normalised_png'),
                                 plot_normal  = True)
    # Uncertainty Z-Map
    write_map_value_distribution(map_vals     = z_map_uncty.select(reference_grid.global_mask().outer_mask_indices()),
                                 output_file  = d_handler.output_handler.get_file('z_map_uncertainty_png'),
                                 plot_normal  = True)
    # Normalised Uncertainty Z-Map
    write_map_value_distribution(map_vals     = z_map_uncty_normalised.select(reference_grid.global_mask().outer_mask_indices()),
                                 output_file  = d_handler.output_handler.get_file('z_map_uncertainty_normalised_png'),
                                 plot_normal  = True)
    # Corrected Z-Map
    write_map_value_distribution(map_vals     = z_map_compl.select(reference_grid.global_mask().outer_mask_indices()),
                                 output_file  = d_handler.output_handler.get_file('z_map_corrected_png'),
                                 plot_normal  = True)
    # Normalised Corrected Z-Map
    write_map_value_distribution(map_vals     = z_map_compl_normalised.select(reference_grid.global_mask().outer_mask_indices()),
                                 output_file  = d_handler.output_handler.get_file('z_map_corrected_normalised_png'),
                                 plot_normal  = True)
    # Plot Q-Q Plot of Corrected Z-Map to see how normal it is
    write_qq_plot_against_normal(map_vals     = z_map_compl_normalised.select(reference_grid.global_mask().outer_mask_indices()),
                                 output_file  = d_handler.output_handler.get_file('z_map_qq_plot_png'))

    # ============================================================================>
    #####
    # LOOK FOR CLUSTERS OF LARGE Z-SCORES
    #####
    # ============================================================================>
    # Contour the grid at a particular Z-Value
    # ============================================================================>
    num_clusters, z_clusters = blob_finder.cluster_high_z_values(z_map=z_map, point_mask=dataset_total_mask)
    # ============================================================================>
    # Too many points to cluster -- probably a bad dataset
    # ============================================================================>
    if num_clusters == -1:
        # This dataset is too noisy to analyse - flag!
        #pandda.datasets.all_masks().set_mask_value(mask_name='noisy zmap', entry_id=d_handler.tag, value=True)
        log_strs.append('Z-Map too noisy to analyse')
        return d_handler, m_handler, log_strs

    # ============================================================================>
    #####
    # FILTER/SELECT CLUSTERS OF Z-SCORES
    #####
    # ============================================================================>
    # Filter the clusters by size and peak height
    # ============================================================================>
    if num_clusters > 0:
        num_clusters, z_clusters = blob_finder.filter_z_clusters_1(z_clusters=z_clusters)
        blob_finder.validate_clusters(z_clusters)
        if num_clusters == 0: log_strs.append('===> Minimum cluster peak/size not reached.')
    # ============================================================================>
    # Filter the clusters by distance from protein
    # ============================================================================>
    if num_clusters > 0:
        num_clusters, z_clusters = blob_finder.filter_z_clusters_2(z_clusters       = z_clusters,
                                                                   grid_origin_cart = reference_dataset.origin_shift(),
                                                                   ref_structure    = reference_dataset.new_structure().hierarchy)
        blob_finder.validate_clusters(z_clusters)
        if num_clusters == 0: log_strs.append('===> Clusters too far from protein.')
    # ============================================================================>
    # Group Nearby Clusters Together
    # ============================================================================>
    if num_clusters > 0:
        num_clusters, z_clusters = blob_finder.group_clusters(z_clusters=z_clusters)
        blob_finder.validate_clusters(z_clusters)
    # ============================================================================>
    # Filter the clusters by symmetry equivalence
    # ============================================================================>
    if num_clusters > 0:
        num_clusters, z_clusters = blob_finder.filter_z_clusters_3(z_clusters       = z_clusters,
                                                                   grid_origin_cart = reference_dataset.origin_shift(),
                                                                   ref_unit_cell    = reference_dataset.unit_cell,
                                                                   ref_sym_ops      = reference_dataset.crystal_contact_generators,
                                                                   ref_structure    = reference_dataset.new_structure().hierarchy     )
        blob_finder.validate_clusters(z_clusters)

    # ============================================================================>
    #####
    # WRITE MAPS
    #####
    # ============================================================================>
    if (num_clusters != 0) or args.output.developer.write_maps_for_empty_datasets:
        # ============================================================================>
        # NATIVE MAPS (ROTATED)
        # ============================================================================>
        # Z-map
        write_array_to_map( output_file = d_handler.output_handler.get_file('native_z_map'),
                            map_data    = rotate_map(grid=reference_grid, d_handler=d_handler, map_data=z_map),
                            grid        = reference_grid   )
        # ============================================================================>
        # REFERENCE FRAME MAPS (NOT ROTATED)
        # ============================================================================>
        # Write the sampled map
        write_array_to_map( output_file = d_handler.output_handler.get_file('sampled_map'),
                            map_data    = m_handler.map,
                            grid        = reference_grid   )
        # Write the mean-difference map
        write_array_to_map( output_file = d_handler.output_handler.get_file('mean_diff_map'),
                            map_data    = d_map,
                            grid        = reference_grid   )
        # Write Chosen Z-Map (ALWAYS WRITE THIS MAP)
        write_array_to_map( output_file = d_handler.output_handler.get_file('z_map'),
                            map_data    = z_map,
                            grid        = reference_grid   )
        # ============================================================================>
        # MASKS/BLOBS/GENERIC MAPS
        # ============================================================================>
        if args.output.developer.write_grid_masks:
            # Write map of grid + symmetry mask
            map_mask = numpy.zeros(reference_grid.grid_size_1d(), dtype=int)
            map_mask.put(map(reference_grid.grid_indexer(), dataset_outer_mask), 1)
            map_mask = flex.double(map_mask.tolist()); map_mask.reshape(flex.grid(reference_grid.grid_size()))
            write_array_to_map( output_file = d_handler.output_handler.get_file('grid_mask'),
                                map_data    = map_mask,
                                grid        = reference_grid    )
            # Write map of where the blobs are (high-Z mask)
            highz_points = []; [highz_points.extend(list(x[0])) for x in z_clusters]
            highz_points = [map(int, v) for v in highz_points]
            highz_map_array = numpy.zeros(reference_grid.grid_size_1d(), dtype=int)
            highz_map_array.put(map(reference_grid.grid_indexer(), list(highz_points)), 1)
            map_mask = flex.double(highz_map_array.tolist()); map_mask.reshape(flex.grid(reference_grid.grid_size()))
            write_array_to_map( output_file = d_handler.output_handler.get_file('high_z_mask'),
                                map_data    = map_mask,
                                grid        = reference_grid    )
        # ============================================================================>
        # Write different Z-Maps? (Probably only needed for testing)
        # ============================================================================>
        if args.output.developer.write_all_z_map_types:
            # Write Naive Z-Map
            write_array_to_map( output_file = d_handler.output_handler.get_file('z_map_naive'),
                                map_data    = z_map_naive,
                                grid        = reference_grid   )
            # Write Normalised Naive Z-Map
            write_array_to_map( output_file = d_handler.output_handler.get_file('z_map_naive_normalised'),
                                map_data    = z_map_naive_normalised,
                                grid        = reference_grid   )
            # Write Uncertainty Z-Map
            write_array_to_map( output_file = d_handler.output_handler.get_file('z_map_uncertainty'),
                                map_data    = z_map_uncertainty,
                                grid        = reference_grid   )
            # Write Normalised Uncertainty Z-Map
            write_array_to_map( output_file = d_handler.output_handler.get_file('z_map_uncertainty_normalised'),
                                map_data    = z_map_uncertainty_normalised,
                                grid        = reference_grid   )
            # Write Corrected Z-Map
            write_array_to_map( output_file = d_handler.output_handler.get_file('z_map_corrected'),
                                map_data    = z_map_compl,
                                grid        = reference_grid   )
            # Write Normalised Corrected Z-Map
            write_array_to_map( output_file = d_handler.output_handler.get_file('z_map_corrected_normalised'),
                                map_data    = z_map_compl_normalised,
                                grid        = reference_grid   )

    # ============================================================================>
    # Skip to next dataset if no clusters found
    # ============================================================================>
    if num_clusters > 0:
        log_strs.append('===> {!s} Cluster(s) found.'.format(num_clusters))
    else:
        log_strs.append('===> No Clusters found.')
        return (d_handler, m_handler, log_strs)
    assert num_clusters > 0, 'NUMBER OF CLUSTERS AFTER FILTERING == 0!'

    # ============================================================================>
    # Process the identified features
    # ============================================================================>
    for event_idx, (event_points, event_values) in enumerate(z_clusters):
        # Number events from 1
        event_num = event_idx + 1
        # Create a unique identifier for this event
        event_key = (d_handler.tag, event_num)
        # ============================================================================>
        # Create a point cluster object
        # ============================================================================>
        point_cluster = PointCluster(id=event_key, points=event_points, values=event_values)
        # ============================================================================>
        # Estimate the background correction of the detected feature
        # ============================================================================>
        # Extract sites for this cluster and estimate the background correction for the event
        log_strs.append('----------------------------------->>>')
        log_strs.append('Estimating Event {!s} Background Correction'.format(event_num))
        # Generate custom grid mask for this dataset
        event_mask = grid_mask( cart_sites = flex.vec3_double(point_cluster.points)*reference_grid.grid_spacing(),
                                grid_size  = reference_grid.grid_size(),
                                unit_cell  = reference_grid.unit_cell(),
                                max_dist   = 1.0,
                                min_dist   = 0.0 )
        expanded_points = list(event_mask.outer_mask())
        log_strs.append('=> Event sites ({!s} points) expanded to {!s} points'.format(len(point_cluster.points), len(expanded_points)))
        # Use the outer mask to provide a reference point for the correlation to the mean map
        reference_points = list(reference_grid.global_mask().inner_mask())
        # Calculate the correlation with the mean map as different amounts of the mean map are subtracted
        event_remain_est = estimate_event_background_correction(
            ref_map          = map_analyser.statistical_maps.mean_map,
            query_map        = m_handler.map,
            feature_region   = expanded_points,
            reference_region = reference_points,
            min_remain       = 1.0-params.background_correction.max_bdc,
            max_remain       = 1.0-params.background_correction.min_bdc,
            bdc_increment    = params.background_correction.increment,
            method           = 'value',
            verbose          = settings.verbose)
        log_strs.append('=> Event Background Correction estimated as {!s}'.format(1-event_remain_est))
        # ============================================================================>
        # Calculate background-corrected map
        # ============================================================================>
        event_map = calculate_bdc_subtracted_map(
            ref_map   = map_analyser.statistical_maps.mean_map,
            query_map = m_handler.map,
            bdc       = 1.0 - event_remain_est)
        # Write out these array (reference and native frames)
        write_array_to_map( output_file = d_handler.output_handler.get_file('event_map').format(event_num, event_remain_est),
                            map_data    = event_map,
                            grid        = reference_grid     )
        write_array_to_map( output_file = d_handler.output_handler.get_file('native_event_map').format(event_num, event_remain_est),
                            map_data    = rotate_map(grid=reference_grid, d_handler=d_handler, map_data=event_map),
                            grid        = reference_grid     )
#                            map_data    = rotate_map(grid=reference_grid, d_handler=d_handler, map_data=event_map, align_on_grid_point=map(int,point_cluster.centroid)))
        # ============================================================================>
        # Find the nearest calpha to the event
        # ============================================================================>
        rt_lab = reference_grid.partition().query_by_grid_points([map(int,point_cluster.centroid)])[0]
        log_strs.append('=> Nearest C-alpha to event: Chain {}, Residue {}'.format(rt_lab[0], rt_lab[1]))
        # ============================================================================>
        # Create an event object
        # ============================================================================>
        event_obj = Event(id=point_cluster.id, cluster=point_cluster)
        event_obj.info.estimated_pseudo_occupancy = event_remain_est
        event_obj.info.estimated_bdc              = 1.0 - event_remain_est
        # ============================================================================>
        # Append to dataset handler
        # ============================================================================>
        d_handler.events.append(event_obj)

        # ============================================================================>
        # Write out the graph of the background correction correlations
        # ============================================================================>
        #if pandda.settings.plot_graphs:
        #    pandda.something()

    return (d_handler, m_handler, log_strs)

# ============================================================================>
#
###                 PanDDA Initialisation Functions
#
# ============================================================================>

def pandda_dataset_setup(pandda):
    """Initialise the pandda object and load the input datasets"""

    pandda.log('',True)
    pandda.log('=====================================================>>>', True)
    pandda.log('Running PanDDA Setup', True)
    pandda.log('=====================================================>>>', True)
    pandda.log('',True)

    # ============================================================================>
    #####
    # Build list of files in data directories
    #####
    # ============================================================================>
    input_files = pandda.build_input_list()
    # Check that some datasets have been found or already loaded
    if (not pandda.datasets.all()) and (not input_files):
        raise Sorry('NO DATASETS HAVE BEEN SELECTED FOR ANALYSIS OR LOADED FROM PREVIOUS RUNS')
    # Check to see if we're reusing statistical maps
    if (not pandda.args.method.recalculate_statistical_maps) and pandda.stat_maps.get_resolutions():
        pandda.log('----------------------------------->>>', True)
        pandda.log('Pre-existing statistical maps (from previous runs) have been found and will be reused:', True)
        pandda.log('Resolutions of reloaded maps: {!s}'.format(', '.join(map(str,pandda.stat_maps.get_resolutions()))), True)
    # Check that enough datasets have been found
    elif pandda.datasets.size() + len(input_files) < pandda.params.analysis.min_build_datasets:
        pandda.log('----------------------------------->>>', True)
        pandda.log('NOT ENOUGH DATASETS HAVE BEEN LOADED FOR ANALYSIS', True)
        pandda.log('Number loaded ({!s}) is less than the {!s} needed.'.format(pandda.datasets.size()+len(input_files), pandda.params.analysis.min_build_datasets), True)
        pandda.log('This value is controlled by changing pandda.params.analysis.min_build_datasets', True)
        raise Sorry('NOT ENOUGH DATASETS HAVE BEEN LOADED FOR ANALYSIS')
    # If dry_run, exit after initial search
    if pandda.args.exit_flags.dry_run:
        raise SystemExit('Dry Run Only: Exiting')
    # Load and process input files
    if input_files and (pandda.args.input.max_new_datasets > 0):
        # ============================================================================>
        # Limit the number of files that can be added at once to a pandda
        # ============================================================================>
        if len(input_files) > pandda.args.input.max_new_datasets:
            pandda.log('----------------------------------->>>', True)
            pandda.log('Limiting the number of new datasets that are added to the pandda analysis (controlled by input.max_new_datasets)')
            pandda.log('Limiting analysis to the first {} of {} datasets'.format(pandda.args.input.max_new_datasets, len(input_files)))
            input_files = input_files[:pandda.args.input.max_new_datasets]
            assert len(input_files) == pandda.args.input.max_new_datasets, 'Something has gone wrong'
        # ============================================================================>
        #####
        # Add new files and load datasets
        #####
        # ============================================================================>
        pandda.add_new_files(input_files)
        pandda.load_new_datasets()
        pandda.initialise_dataset_masks_and_tables()
        pandda.check_loaded_datasets(datasets=pandda.datasets.all())
        # ============================================================================>
        #####
        # Set Reference Dataset
        #####
        # Select the reference dataset
        # ============================================================================>
        if not pandda.reference_dataset():
            # Filter datasets against the provided filter pdb if given
            if pandda.args.input.filter.pdb is not None:
                pandda.log('----------------------------------->>>', True)
                pandda.log('Filtering datasets against the provided pdb structure (defined by pandda.input.filter.pdb)', True)
                pandda.filter_datasets_1(filter_dataset=DatasetHandler(dataset_number=0, pdb_filename=pandda.args.input.filter.pdb))
            # Use given reference dataset, or select reference dataset
            if pandda.args.input.reference.pdb and pandda.args.input.reference.mtz:
                pandda.log('----------------------------------->>>', True)
                pandda.log('Reference Provided by User', True)
                ref_pdb, ref_mtz = pandda.args.input.reference.pdb, pandda.args.input.reference.mtz
            else:
                pandda.log('----------------------------------->>>', True)
                pandda.log('Selecting reference dataset from loaded datasets', True)
                ref_pdb, ref_mtz = pandda.select_reference_dataset(method='resolution')
            # Load the reference dataset
            pandda.load_reference_dataset(ref_pdb=ref_pdb, ref_mtz=ref_mtz)
        # ============================================================================>
        #####
        # Scale, Align and Initial-Filter All Data
        #####
        # ============================================================================>
        # Filter out datasets with different protein structures
        pandda.filter_datasets_1()
        # Load reflection data for each dataset
        pandda.load_reflection_data(ampl_label=pandda.params.maps.ampl_label,
                                    phas_label=pandda.params.maps.phas_label)
        # Align structures to the reference
        pandda.align_datasets(method=pandda.params.alignment.method)
        # Pickle the new loaded datasets
        pandda.pickle_the_pandda(components=['datasets'], datasets=pandda.datasets.mask(mask_name='old datasets', invert=True))
    else:
        # ============================================================================>
        # Rebuild the masks of the rejected datasets (quick)
        # ============================================================================>
        pandda.initialise_dataset_masks_and_tables()
        pandda.check_loaded_datasets(datasets=pandda.datasets.all())
        pandda.filter_datasets_1()

    # ============================================================================>
    # Check that enough VALID datasets have been found
    # ============================================================================>
    if (not pandda.args.method.recalculate_statistical_maps) and pandda.stat_maps.get_resolutions():
        # Using existing maps - don't need to check
        pass
    elif pandda.datasets.size(mask_name='rejected - total', invert=True) < pandda.params.analysis.min_build_datasets:
        pandda.log('----------------------------------->>>', True)
        pandda.log('NOT ENOUGH (NON-REJECTED) DATASETS TO BUILD DISTRIBUTIONS!', True)
        pandda.log('Number loaded ({!s}) is less than the {!s} needed.'.format(pandda.datasets.size(mask_name='rejected - total', invert=True), pandda.params.analysis.min_build_datasets), True)
        pandda.log('This value is defined by pandda.params.analysis.min_build_datasets', True)
        raise Sorry('NOT ENOUGH DATASETS LOADED')

    # ============================================================================>
    #####
    # Filter and Analyse the Datasets
    #####
    # ============================================================================>
    # Collate dataset variables and parameters
    pandda.collate_dataset_variables()
    # Filter out the datasets that are not isomorphous and therefore incomparable
    pandda.filter_datasets_2()
    # Calculate the calpha rmsd to the reference dataset
    pandda.calculate_dataset_rmsds_to_reference()
    # Calculate the variation between the alignments of neighbouring C-alphas
    pandda.analyse_alignment_variation()

    return

def pandda_grid_setup(pandda):
    """Generate the grid objects for the analysis"""

    pandda.log('',True)
    pandda.log('=====================================================>>>', True)
    pandda.log('Running Grid Setup',True)
    pandda.log('=====================================================>>>', True)
    pandda.log('',True)

    # ============================================================================>
    # Re-select the reference dataset - NOT YET IMPLEMENTED (TODO)
    # ============================================================================>
#    highest_res_dh = sorted(pandda.datasets.mask(mask_name='rejected - total', invert=True), key=lambda d: d.reflection_data().max_min_resolution()[1])[0]
#    pandda.log('Reselected Reference Dataset (used for map scaling): {}'.format(highest_res_dh.tag))

    # ============================================================================>
    #####
    # Create Sampling Grid (for generated maps)
    #####
    # Create reference grid based on the reference structure
    # ============================================================================>
    if pandda.reference_grid() is None:
        # Create parameter for setting grid spacing (multiple grids?)
        pandda.create_reference_grid(   grid_spacing     = pandda.params.maps.grid_spacing,
                                        expand_to_origin = False,
                                        buffer           = pandda.params.masks.outer_mask+pandda.params.maps.padding    )
        # Create various masks to define regions of the grid by distance to the protein and symmetry copies
        pandda.mask_reference_grid()
    # ============================================================================>
    # Store for reuse
    # ============================================================================>
    # Pickle all of the large arrays so they can be reloaded
    # ============================================================================>
    pandda.pickle_the_pandda(components=['grid'])

    # ============================================================================>
    # If setup_only, exit after initial search
    # ============================================================================>
    if pandda.args.exit_flags.setup_only:
        pandda.exit(error=False)
        raise SystemExit('Setup Only: Exiting')

    return

# ============================================================================>
#
###                 PanDDA Processing Functions
#
# ============================================================================>

def pandda_main_loop(pandda):
    """Calculate the statistical maps, and then look for events in each dataset"""

    pandda.log('',True)
    pandda.log('=====================================================>>>', True)
    pandda.log('Running Main PanDDA Analysis', True)
    pandda.log('=====================================================>>>', True)
    pandda.log('',True)

    # ============================================================================>
    # Validate/Reset the loaded datasets (those to be analysed)
    # ============================================================================>
    pandda.reset_loaded_datasets()

    # ============================================================================>
    #####
    # PRE-ANALYSIS ANALYSIS (DUMP OF DATASET PARAMETERS)
    #####
    # ============================================================================>
    pandda.write_output_csvs()
    analyse_graphs.write_dataset_summary_graphs(pandda)
    analyse_html.write_initial_html(pandda)

    # ============================================================================>
    #####
    # Update Settings
    #####
    # ============================================================================>
    # Update the resolution limits using the resolution limits from the datasets supplied
    if pandda.params.analysis.dynamic_res_limits:
        pandda.log('----------------------------------->>>')
        pandda.log('UPDATING RESOLUTION LIMITS -')
        pandda.set_low_resolution(  min(pandda.params.analysis.high_res_lower_limit,
                                    max(pandda.tables.dataset_info['high_resolution'])))
        pandda.set_high_resolution( max(pandda.params.analysis.high_res_upper_limit,
                                    pandda.reference_dataset().mtz_summary.high_res,
                                    min(pandda.tables.dataset_info['high_resolution'])))
    else:
        pandda.log('----------------------------------->>>')
        pandda.log('**NOT** UPDATING RESOLUTION LIMITS -')
        pandda.set_low_resolution(  pandda.params.analysis.high_res_lower_limit)
        pandda.set_high_resolution( pandda.params.analysis.high_res_upper_limit)
    # Print new resolution limits
    pandda.log('LOW RESOLUTION:  {!s}'.format(pandda.get_low_resolution()))
    pandda.log('HIGH RESOLUTION: {!s}'.format(pandda.get_high_resolution()))

    # ============================================================================>
    #####
    # PREPARE VARIABLES TO LOOP OVER RESOLUTION SHELLS
    #####
    # ============================================================================>
    # Calculate cutoffs for resolution shells
    # ============================================================================>
    if (not pandda.args.method.recalculate_statistical_maps) and pandda.stat_maps.get_resolutions():
        # Use the resolution limits of previously run campaigns
        res_limits = pandda.stat_maps.get_resolutions()
        assert res_limits, 'No Resolution Limits found from statistical maps: {!s}'.format(res_limits)
        min_limit = min(res_limits)
        max_limit = max(res_limits)
    else:
        # Set pandda.args.method.recalculate_statistical_maps to True as no statistical maps have been found
        if not pandda.args.method.recalculate_statistical_maps:
            pandda.log('No Statistical Maps Found: Setting pandda.args.method.recalculate_statistical_maps to True', True)
            pandda.args.method.recalculate_statistical_maps = True
        # ============================================================================>
        # Select resolution limits based on dataset resolutions and given arguments
        # ============================================================================>
        if pandda.params.analysis.dynamic_res_limits:
            # Round the high limit DOWN to the nearest 0.01
            min_limit = round(pandda.get_high_resolution()- 0.005, 2)    # i.e. 1.344 -> 1.34
            # Round the low limit UP to the nearest 0.01
            max_limit = round(pandda.get_low_resolution() + 0.005, 2)    # i.e. 3.423 -> 3.43
        else:
            # Take the arguments as given by the user
            min_limit = pandda.get_high_resolution()
            max_limit = pandda.get_low_resolution()
        # ============================================================================>
        # Create resolution shells (or single limit)
        # ============================================================================>
        if not pandda.params.analysis.high_res_increment:
            # No variable cutoff - select all
            res_limits = [max_limit]  # i.e. [2]
        else:
            # Calculate a range of resolution limits
            shell_width = pandda.params.analysis.high_res_increment
            res_limits = [round(x, 4) for x in numpy.arange(min_limit, max_limit, shell_width).tolist()]
            # Append the rounded max_limit to the end to ensure that the last dataset is processed
            if (not res_limits) or (res_limits[-1] != max_limit):
                res_limits.append(max_limit)
    # ==================================================>
    #####
    # Initialise for iterations over shells
    #####
    # ============================================================================>
    # Analyse all datasets from high_shell_limit -> cut_resolution (initialise high_shell_limit to 0)
    high_shell_limit = 0
    # Record how many datasets are processed at each resolution
    resolution_count = {}
    # ============================================================================>
    # Report
    # ============================================================================>
    pandda.log('----------------------------------->>>', True)
    if len(res_limits)==1:
        pandda.log('Analysing All Maps at {!s}A'.format(max_limit), True)
    else:
        pandda.log('Analysing Resolution Shells from {!s} -> {!s}A'.format(min_limit, max_limit), True)
        pandda.log('Limits: {!s}'.format(', '.join(map(str,res_limits))), True)

    # ============================================================================>
    #####
    # ANALYSE DATASETS - ITERATE THROUGH RESOLUTION SHELLS
    #####
    # ============================================================================>
    t_analysis_start = time.time()
    # ==================================================>
    pandda.log('----------------------------------->>>', True)
    pandda.log('Dataset Analysis Started: {!s}'.format(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime(t_analysis_start))), True)

    for cut_resolution in res_limits:

        # Which resolutions will be processed in this shell
        pandda.log('', True)
        pandda.log('=====================================================>>>', True)
        pandda.log('Processing Datasets from {!s}A -> {!s}A'.format(high_shell_limit, cut_resolution), True)
        pandda.log('=====================================================>>>', True)
        pandda.log('',True)

        # ============================================================================>
        # Select datasets to build the distributions
        # ============================================================================>
        building_mask_name = 'selected for building @ {!s}A'.format(cut_resolution)
        if pandda.args.method.recalculate_statistical_maps:
            pandda.log('----------------------------------->>>', True)
            pandda.log('Selecting Building Mask', True)
            building_mask = pandda.select_for_building_distributions(high_res_cutoff = cut_resolution,
                                                                     building_mask_name = building_mask_name)
            # Check that we have enough datasets to build distributions
            if sum(building_mask) < pandda.params.analysis.min_build_datasets:
                # Don't have enough to construct robust distributions
                pandda.log('NOT ENOUGH DATASETS TO CALCULATE DISTRIBUTIONS ({!s}<{!s})'.format(sum(building_mask),pandda.params.analysis.min_build_datasets), True)
                continue
            else:
                # Have enough to generate robust distributions
                pandda.log('ENOUGH DATASETS -> PROCESSING THIS RESOLUTION', True)
                pandda.log('Building Distributions using {!s} Datasets'.format(sum(building_mask)), True)
        else:
            pandda.log('----------------------------------->>>', True)
            pandda.log('**NOT** Selecting Building Mask (Using Existing Statistical Maps)', True)
            # Create a dummy mask as we won't be using any datasets for building
            building_mask = [False]*pandda.datasets.size()
            pandda.datasets.all_masks().add_mask(mask_name=building_mask_name, mask=building_mask)

        # ============================================================================>
        # Select the datasets to analyse
        # ============================================================================>
        analysis_mask_name = 'selected for analysis @ {!s}A'.format(cut_resolution)
        pandda.log('----------------------------------->>>', True)
        pandda.log('Selecting Analysis Mask', True)
        analysis_mask = pandda.select_for_analysis(high_res_large_cutoff = cut_resolution,
                                                   high_res_small_cutoff = high_shell_limit,
                                                   analysis_mask_name = analysis_mask_name)
        # Check that there're some datasets to analyse
        if sum(analysis_mask) == 0:
            pandda.log('NO DATASETS TO ANALYSE @ {!s}'.format(cut_resolution))
            continue
        else:
            pandda.log('Calculating Z-Maps for {!s} Datasets'.format(sum(analysis_mask)), True)

        # ============================================================================>
        # Combine the masks as we will need to load maps for all datasets
        # ============================================================================>
        pandda.log('----------------------------------->>>', True)
        pandda.log('Combining (Analysis and Building) Masks', True)
        map_load_mask = pandda.datasets.all_masks().combine_masks([analysis_mask_name, building_mask_name])
        map_load_mask_name = 'selected for loading maps @ {!s}A'.format(cut_resolution)
        pandda.datasets.all_masks().add_mask(mask_name=map_load_mask_name, mask=map_load_mask)

        # ============================================================================>
        # Report
        # ============================================================================>
        pandda.log('----------------------------------->>>')
        pandda.log('Mask Names for Building, Loading, and Analysis')
        pandda.log('Building ({!s} datasets): {!s}'.format(sum(building_mask), building_mask_name))
        pandda.log('Load Map ({!s} datasets): {!s}'.format(sum(map_load_mask), map_load_mask_name))
        pandda.log('Analysis ({!s} datasets): {!s}'.format(sum(analysis_mask), analysis_mask_name))
        pandda.log('----------------------------------->>>', True)
        pandda.log('Loading Maps for {!s} Datasets at {!s}A'.format(pandda.datasets.size(mask_name=map_load_mask_name), cut_resolution), True)

        # ============================================================================>
        #####
        # LOAD AND ANALYSE MAPS
        #####
        # ============================================================================>
        t_loop_start = time.time()
        # ==================================================>
        pandda.log('{!s}A Analysis Started: {!s}'.format(cut_resolution, time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime(t_loop_start))), True)

        # ============================================================================>
        # Update limit (cut resolution becomes top limit in next shell)
        # ============================================================================>
        high_shell_limit = cut_resolution

        # ============================================================================>
        #####
        # Truncate data, load and scale maps to reference dataset
        #####
        # ============================================================================>
        # Truncate the data to a particular resolution
        # ============================================================================>
        pandda.truncate_scaled_data(dataset_handlers = pandda.datasets.mask(mask_name=map_load_mask_name))
        # ============================================================================>
        # Load the reference map so that we can scale the individual maps to this
        # ============================================================================>
        ref_map_holder = pandda.load_reference_map( map_resolution = cut_resolution )
        # ============================================================================>
        # Load the required maps
        # ============================================================================>
        map_holder_list = pandda.load_and_morph_maps( dataset_handlers = pandda.datasets.mask(mask_name=map_load_mask_name),
                                                      ref_map_holder   = ref_map_holder,
                                                      map_resolution   = cut_resolution     )
        # ============================================================================>
        # Extract the statistical maps at this resolution (if requested)
        # ============================================================================>
        if pandda.args.method.recalculate_statistical_maps:
            # No Statistical maps - set to None so that they are created
            statistical_maps = None
        else:
            # Try to find a statistical map at this resolution
            statistical_maps = pandda.stat_maps.get(cut_resolution)
        # ============================================================================>
        # Create an object to hold all of the maps, and can be used to calculate the mean maps, etc...
        # ============================================================================>
        map_analyser = PanddaMapAnalyser(   dataset_maps     = map_holder_list,
                                            meta             = Meta({'resolution'   : cut_resolution,
                                                                     'grid_size'    : pandda.reference_grid().grid_size(),
                                                                     'grid_size_1d' : pandda.reference_grid().grid_size_1d()}),
                                            statistical_maps = statistical_maps,
                                            parent           = pandda,
                                            log              = pandda.log   )
        # ============================================================================>
        # Add the analysis mask to the map analyser
        # ============================================================================>
        map_analyser.dataset_maps.all_masks().add_mask(mask_name=analysis_mask_name, mask=[False]*map_analyser.dataset_maps.size())
        for dh in pandda.datasets.mask(mask_name=analysis_mask_name):
            map_analyser.dataset_maps.all_masks().set_mask_value(mask_name=analysis_mask_name, entry_id=dh.tag, value=True)
        # ============================================================================>
        # Add the building mask to the map analyser
        # ============================================================================>
        map_analyser.dataset_maps.all_masks().add_mask(mask_name=building_mask_name, mask=[False]*map_analyser.dataset_maps.size())
        for dh in pandda.datasets.mask(mask_name=building_mask_name):
            map_analyser.dataset_maps.all_masks().set_mask_value(mask_name=building_mask_name, entry_id=dh.tag, value=True)

        # ============================================================================>
        #####
        # Calculate Statistical Maps (if required) and Calculate Dataset Uncertainties
        #####
        # ============================================================================>
        if pandda.args.method.recalculate_statistical_maps:
            pandda.log('----------------------------------->>>', True)
            pandda.log('Building Map Distributions for {!s} Datasets at {!s}A'.format(pandda.datasets.size(mask_name=building_mask_name), cut_resolution), True)
        else:
            pandda.log('----------------------------------->>>', True)
            pandda.log('Using Existing Map Distributions at {!s}A'.format(cut_resolution), True)
            assert pandda.datasets.size(mask_name=building_mask_name) == 0, 'BUILDING MASKS HAVE BEEN SELECTED WHEN MAPS ALREADY EXIST'
            assert map_analyser.dataset_maps.size(mask_name=building_mask_name) == 0, 'BUILDING MASKS HAVE BEEN SELECTED WHEN MAPS ALREADY EXIST'
        # ============================================================================>
        # Calculate the mean map
        # ============================================================================>
        if (not map_analyser.statistical_maps.mean_map) or pandda.args.method.recalculate_statistical_maps:
            map_analyser.calculate_mean_map(masked_idxs = pandda.reference_grid().global_mask().outer_mask_indices(),
                                            mask_name   = building_mask_name)
        # ============================================================================>
        # If only Mean Map Requested -- exit
        # ============================================================================>
        if pandda.args.exit_flags.calculate_first_mean_map_only:
            write_array_to_map( output_file = pandda.output_handler.get_file('mean_map').format(cut_resolution),
                                map_data    = map_analyser.statistical_maps.mean_map,
                                grid        = pandda.reference_grid()     )
            raise SystemExit('Calculating First Mean Map Only: Exiting')
        # ============================================================================>
        # Plot Mean Map against Reference Map - should be fairly similar...
        # ============================================================================>
        # Plot the mean map against the reference map (unsorted)
        analyse_graphs.mean_obs_scatter(        f_name    = os.path.join(pandda.output_handler.get_dir('reference'), 'reference_against_mean_unsorted.png'),
                                                mean_vals = map_analyser.statistical_maps.mean_map.select(pandda.reference_grid().global_mask().outer_mask_indices()),
                                                obs_vals  = ref_map_holder.map.select(pandda.reference_grid().global_mask().outer_mask_indices())           )
        # Plot the mean map against the reference map (sorted)
        analyse_graphs.sorted_mean_obs_scatter( f_name    = os.path.join(pandda.output_handler.get_dir('reference'), 'reference_against_mean_sorted.png'),
                                                mean_vals = sorted(map_analyser.statistical_maps.mean_map.select(pandda.reference_grid().global_mask().outer_mask_indices())),
                                                obs_vals  = sorted(ref_map_holder.map.select(pandda.reference_grid().global_mask().outer_mask_indices()))   )
        # Plot the reference map distribution
        write_map_value_distribution( map_vals    = ref_map_holder.map.select(pandda.reference_grid().global_mask().outer_mask_indices()),
                                      output_file = os.path.join(pandda.output_handler.get_dir('reference'), 'reference_map_distribution.png')     )
        # ============================================================================>
        # Calculate the uncertainty of all loaded maps (needs the mean map to have been calculated)
        # ============================================================================>
        # If not pandda.args.method.recalculate_statistical_maps, then no building_mask datasets will have been added, so don't need to be selective in this function
        map_analyser.calculate_map_uncertainties(masked_idxs=pandda.reference_grid().global_mask().inner_mask_indices())
        # ============================================================================>
        # Plot uncertainties of maps
        # ============================================================================>
        try:
            from ascii_graph import Pyasciigraph
            g=Pyasciigraph()
            graph_data = [(mh.tag, round(mh.meta.map_uncertainty,3)) for mh in map_analyser.dataset_maps.all()]
            pandda.log('----------------------------------->>>', True)
            for l in g.graph(label='Sorted Map Uncertainties (Ascending Order)', data=graph_data, sort=1):
                if l.startswith('#######'): continue
                pandda.log(l.replace(u"\u2588", '=').replace('= ','> '), True)
            pandda.log('----------------------------------->>>', True)
            pandda.log('')
        except ImportError:
            print('IMPORT ERROR (ascii_graph) - CANNOT GENERATE UNCERTAINTY GRAPH')
        except:
            pass
        # ============================================================================>
        # Calculate the statistics of the maps
        # ============================================================================>
        if (not map_analyser.statistical_maps.sadj_map) or pandda.args.method.recalculate_statistical_maps:
            map_analyser.calculate_statistical_maps(masked_idxs = pandda.reference_grid().global_mask().outer_mask_indices(),
                                                    mask_name   = building_mask_name,
                                                    cpus        = pandda.settings.cpus)

#        # ============================================================================>
#        # PICKLE THE DATASETS THAT HAVE JUST BEEN PROCESSED
#        # ============================================================================>
#        pandda.sync_datasets(datasets=pandda.datasets.mask(mask_name=analysis_mask_name), overwrite_dataset_meta=True)
#        pandda.pickle_the_pandda(components=['datasets'], datasets=pandda.datasets.mask(mask_name=analysis_mask_name))
        # ============================================================================>
        # Extract and store the statistical map objects
        # ============================================================================>
        # Only need to add if we're calculating new statistical maps
        if pandda.args.method.recalculate_statistical_maps:
            # Check to see if already in there
            if cut_resolution in pandda.stat_maps.get_resolutions():
                pandda.log('Overwriting existing statistical maps @ {!s}A'.format(cut_resolution))
            pandda.stat_maps.add(stat_map_list=map_analyser.statistical_maps, resolution=cut_resolution, overwrite=True)
            pandda.pickle_the_pandda(components=['stat_maps'])
        # ============================================================================>
        # Pickle map analyser if required
        # ============================================================================>
        if pandda.args.output.pickling.pickle_map_analysers:
            pandda.pickle(pickle_file=pandda.pickle_handler.get_file('map_analyser').format(cut_resolution), pickle_object=map_analyser, overwrite=True)
        # ============================================================================>
        # Write Grid Point Distributions - Standard Function, just to provide an output
        # ============================================================================>
        try:
            from libtbx.math_utils import iceil
            grid_size = pandda.reference_grid().grid_size()
            num_points = min(10, min(grid_size))
            grid_points = zip(*[range(0, s, iceil(s/num_points)) for s in grid_size])
            pandda.write_grid_point_distributions(  grid_points     = grid_points,
                                                    map_analyser    = map_analyser,
                                                    output_filename = None          )
        except:
            print('UNIMPORTANT: FAILED TO WRITE AUTOMATIC DISTRIBUTION OF GRID POINTS')
            raise
        # ============================================================================>
        # Write out Grid Point Distributions for interesting grid points (high modality, etc...)
        # ============================================================================>
        # TODO TODO TODO
        # TODO TODO TODO

        # ============================================================================>
        #####
        # DATA PROCESSING
        #####
        # ============================================================================>
        assert cut_resolution not in resolution_count
        resolution_count[cut_resolution] = []
        # ============================================================================>
        # Create a new "dummy" map analyser for the parallel steps
        # ============================================================================>
        dummy_map_analyser = PanddaMapAnalyser( dataset_maps = None,
                                                meta         = Meta({'resolution'   : cut_resolution,
                                                                     'grid_size'    : pandda.reference_grid().grid_size(),
                                                                     'grid_size_1d' : pandda.reference_grid().grid_size_1d()}),
                                                statistical_maps = map_analyser.statistical_maps,
                                                parent           = None,
                                                log              = None   )
        # ============================================================================>
        # Blob Search Object
        # ============================================================================>
        blob_finder = PanddaZMapAnalyser( params = pandda.params.blob_search,
                                          grid_spacing = pandda.reference_grid().grid_spacing(),
                                          log = pandda.log )
        blob_finder.print_settings()

        # ============================================================================>
        #####
        # Calculate Z-Maps
        #####
        # ============================================================================>
        t_anal_start = time.time()
        # ==================================================>
        pandda.log('----------------------------------->>>', True)
        pandda.log('Preparing to Analyse {!s} Datasets at {!s}A'.format(pandda.datasets.size(mask_name=analysis_mask_name), cut_resolution), True)

        proc_dataset_inp_list = []

        # Iterate through and prepare to calculate z-maps
        for i_dh, d_handler in enumerate(pandda.datasets.mask(mask_name=analysis_mask_name)):

            # Record the which resolution this dataset was analysed at
            resolution_count[cut_resolution].append(d_handler.tag)

            # ============================================================================>
            # Validate/Check/Zero the dataset
            # ============================================================================>
            # Dataset should not have any events
            assert d_handler.events == []

            # ============================================================================>
            # Update datasets masks flag - for this analysis
            # ============================================================================>
            pandda.datasets.all_masks().set_mask_value(mask_name='analysed', entry_id=d_handler.tag, value=True)

            # ============================================================================>
            # Update the dataset meta object -- this is persistent
            # ============================================================================>
            d_handler.meta.analysed = True

            # ============================================================================>
            # Extract the map holder for this dataset
            # ============================================================================>
            m_handler = map_analyser.dataset_maps.get(tag=d_handler.tag)

            # ============================================================================>
            # Compile arguments for this datasets
            # ============================================================================>
            dataset_args = (    d_handler,
                                m_handler,
                                pandda.reference_grid(),
                                pandda.reference_dataset(),
                                dummy_map_analyser,
                                (pandda.args, pandda.params, pandda.settings)   )

            proc_dataset_inp_list.append(dataset_args)

#                # Input for process function
#                # Need to create map_analyser for each dataset (can fill with only stat maps to calculate Z-map)
#                d_handler = args[0]
#                m_handler = args[1]
#                ref_grid  = args[2]
#                map_alysr = args[3]
#                params    = args[4]

            ### XXX ------- preprocessing for the map func above here ------- XXX ###

        pandda.log('----------------------------------->>>', True)
        pandda.log('Calculating and Analysing Z-Maps for {!s} Dataset(s) at {!s}A'.format(pandda.datasets.size(mask_name=analysis_mask_name), cut_resolution), True)

        proc_results = easy_mp.pool_map(func=process_dataset_map_func, args=proc_dataset_inp_list,
                                        processes=pandda.settings.cpus, maxtasksperchild=1)

        pandda.log('----------------------------------->>>', True)
        pandda.log('Updating with results from analysing {!s} Dataset(s) at {!s}A'.format(pandda.datasets.size(mask_name=analysis_mask_name), cut_resolution), True)

        for results in proc_results:

            # Unpack results
            d_handler, m_handler, log_strs = results

            pandda.log('', True)
            pandda.log('======================================>>>', True)
            pandda.log('Z-Map Analysis Results for {}'.format(d_handler.tag), True)
            pandda.log('======================================>>>', True)
            pandda.log('\n'.join(log_strs), True)

            # ============================================================================>
            # STORE ANALYSIS DATA IN DATASET MAP TABLE
            # ============================================================================>
            # Add to the dataset map summary table
            pandda.tables.dataset_map_info.set_value(d_handler.tag, 'analysed_resolution', m_handler.meta.resolution)
            pandda.tables.dataset_map_info.set_value(d_handler.tag, 'map_uncertainty',     round(m_handler.meta.map_uncertainty,3))
            pandda.tables.dataset_map_info.set_value(d_handler.tag, 'obs_map_mean',        round(m_handler.meta.obs_map_mean,3))
            pandda.tables.dataset_map_info.set_value(d_handler.tag, 'obs_map_rms',         round(m_handler.meta.obs_map_rms,3))
            pandda.tables.dataset_map_info.set_value(d_handler.tag, 'z_map_mean',          round(m_handler.meta.z_mean,3))
            pandda.tables.dataset_map_info.set_value(d_handler.tag, 'z_map_std',           round(m_handler.meta.z_stdv,3))
            pandda.tables.dataset_map_info.set_value(d_handler.tag, 'z_map_skew',          round(m_handler.meta.z_skew,3))
            pandda.tables.dataset_map_info.set_value(d_handler.tag, 'z_map_kurt',          round(m_handler.meta.z_kurt,3))

            # ============================================================================>
            # WRITE OUT DATASET INFORMATION TO CSV FILE
            # ============================================================================>
            out_list = pandda.tables.dataset_info.loc[d_handler.tag].append(pandda.tables.dataset_map_info.loc[d_handler.tag])
            out_list.to_csv(path=d_handler.output_handler.get_file('dataset_info'), header=True, index_label='dtag')

            # ============================================================================>
            # Link outputs for datasets with hits
            # ============================================================================>
            if d_handler.events:
                # Create a link to the interesting directories in the initial results directory
                hit_dir = os.path.join(pandda.output_handler.get_dir('interesting_datasets'), d_handler.tag)
                if not os.path.exists(hit_dir): rel_symlink(orig=d_handler.output_handler.get_dir('root'), link=hit_dir)
                # ============================================================================>
                # Add event to the event table
                # ============================================================================>
                for e in d_handler.events:
                    pandda.add_event_to_event_table(d_handler=d_handler, event=e)

            # ============================================================================>
            # Update the master copy of the dataset object
            # ============================================================================>
            master_d_handler = pandda.datasets.get(tag=d_handler.tag)
            master_d_handler.events = d_handler.events
            master_d_handler.child  = m_handler

        # ============================================================================>
        #####
        # Write the summaries
        #####
        # ============================================================================>
        # Write statistical maps
        pandda.write_map_analyser_maps(map_analyser=map_analyser, analysis_mask_name=analysis_mask_name)
        # Write output graphs for map analysis
        if pandda.settings.plot_graphs:
            analyse_graphs.write_map_analyser_graphs(pandda=pandda, map_analyser=map_analyser, analysis_mask_name=analysis_mask_name)

        # ============================================================================>
        #####
        # Print summaries
        #####
        # ============================================================================>
        t_loop_end = time.time()
        # ==================================================>
        pandda.log('')
        pandda.log('=====================================================>>>', True)
        pandda.log('{!s}A Z-Map Processing Time: {!s}'.format(cut_resolution, time.strftime("%H hours:%M minutes:%S seconds", time.gmtime(t_loop_end - t_anal_start))), True)
        pandda.log('{!s}A Total Processing Time: {!s}'.format(cut_resolution, time.strftime("%H hours:%M minutes:%S seconds", time.gmtime(t_loop_end - t_loop_start))), True)
        pandda.log('=====================================================>>>', True)
        pandda.log('@{!s}A:\t {!s}/{!s} New Datasets Analysed'.format(cut_resolution, pandda.datasets.size(mask_name=analysis_mask_name), pandda.datasets.size(mask_name='old datasets', invert=True)))
        pandda.log('Total:\t {!s}/{!s} New Datasets Analysed'.format(pandda.datasets.size(mask_name='analysed'), pandda.datasets.size(mask_name='old datasets', invert=True)))
        pandda.log('=====================================================>>>', True)
        pandda.log('')

        # ============================================================================>
        #####
        # PICKLE THE DATASETS THAT HAVE JUST BEEN PROCESSED
        #####
        # ============================================================================>
        # Clear the linked maps if requested
        if not pandda.args.output.pickling.pickle_dataset_maps:
            for d_handler in pandda.datasets.mask(mask_name=analysis_mask_name):
                d_handler.child = None
        # Pickle the datasets
        pandda.sync_datasets(datasets=pandda.datasets.mask(mask_name=analysis_mask_name), overwrite_dataset_meta=True)
        pandda.pickle_the_pandda(components=['datasets'], datasets=pandda.datasets.mask(mask_name=analysis_mask_name))

        # ============================================================================>
        # DELETE THE MAP ANALYSER TO FREE UP MEMORY
        # ============================================================================>
        map_analyser.parent = None
        del map_analyser
        map_analyser = None
        # ============================================================================>
        # DELETE THE PROCESSED MAPS + Z-MAPS TO SAVE MEMORY (JUST BE CAREFUL NOT TO OVERWRITE THEIR PICKLED FILES!)
        # ============================================================================>
        for d_handler in pandda.datasets.mask(mask_name=analysis_mask_name):
            d_handler.child = None
        # ============================================================================>
        # Launch Garbage-Collection Manually just to be sure
        # ============================================================================>
        gc.collect()

        # ============================================================================>
        #####
        # LIVE ANALYSIS - RUN AT THE END OF EACH LOOP
        #####
        # Analyse the processed data
        # ============================================================================>
        # Extract all events from datasets
        all_events=[]; [all_events.extend(d.events) for d in pandda.datasets.all()]
        # Process any identified events and update output files
        if all_events:
            # ==================================================>
            # Cluster events and add to pandda tables
            # ==================================================>
            site_list = pandda.cluster_events_and_update(events=all_events)
            # ==================================================>
            # Update the output analysis files
            # ==================================================>
            pandda.write_output_csvs()
            analyse_html.write_analyse_html(pandda)
            # ==================================================>

    # ============================================================================>
    # Ensure that the collation step happens even if no resolutions are processed
    # ============================================================================>
    if not resolution_count:
        # Extract all events from datasets
        all_events=[]; [all_events.extend(d.events) for d in pandda.datasets.all()]
        # Process any identified events and update output files
        if all_events:
            # ==================================================>
            # Cluster events and add to pandda tables
            # ==================================================>
            site_list = pandda.cluster_events_and_update(events=all_events)
            # ==================================================>
            # Update the output analysis files
            # ==================================================>
            pandda.write_output_csvs()
            analyse_html.write_analyse_html(pandda)
            # ==================================================>

    # Pickle statistical maps again just to be sure
    pandda.pickle_the_pandda(components=['stat_maps'])

    # ============================================================================>
    #####
    # END OF MAIN LOOP
    #####
    # ============================================================================>
    t_analysis_end = time.time()
    # ==================================================>
    pandda.log('=====================================================>>>', True)
    pandda.log('Total Analysis Time: {!s}'.format(time.strftime("%H hours:%M minutes:%S seconds", time.gmtime(t_analysis_end - t_analysis_start))), True)
    pandda.log('=====================================================>>>', True)

    return

# ============================================================================>
#
###                 PanDDA Output + Wrap-Up Functions
#
# ============================================================================>

def pandda_end(pandda):
    """Write the output summaries and close"""

    pandda.log('',True)
    pandda.log('=====================================================>>>', True)
    pandda.log('Writing PanDDA End-of-Analysis Summary', True)
    pandda.log('=====================================================>>>', True)
    pandda.log('',True)

    # ============================================================================>
    #####
    # FINAL ANALYSIS - RUN ONCE AT THE END OF THE PROGRAM
    #####
    # Analyse all of the processed data
    # ============================================================================>

    # Collate events
    event_total, event_num, all_dataset_events = pandda.collate_event_counts()
    # Print some slightly less important information
    pandda.log('----------------------------------->>>', False)
    for d_tag, event_count in event_num:
        pandda.log('Dataset {!s}: {!s} Events'.format(d_tag, event_count), False)
    # Print a summary of the number of identified events
    pandda.log('----------------------------------->>>', True)
    pandda.log('Total Datasets with Events: {!s}'.format(len(event_num)), True)
    pandda.log('Total Events: {!s}'.format(event_total), True)

    pandda.log('----------------------------------->>>', True)
    pandda.log('Potentially Useful Shortcuts for Future Runs')
    if event_num:
        pandda.log('no_build={!s}'.format(','.join(zip(*event_num)[0])))

    # ============================================================================>
    #####
    # SUMMARIES ------------------------------>>>
    #####
    # ============================================================================>

    pandda.log('----------------------------------->>>', True)
    pandda.log('Writing final output files', True)
    pandda.write_output_csvs()
    analyse_html.write_analyse_html(pandda)

    # ============================================================================>
    # SCREEN GRAPHS -------------------------->>>
    # ============================================================================>

    resolution_counts = pandda.tables.dataset_map_info['analysed_resolution'].value_counts().sort_index()
    graph_data = [(str(r), c) for r,c in resolution_counts.iteritems()]
    if graph_data:
        try:
            from ascii_graph import Pyasciigraph
            g=Pyasciigraph()
            pandda.log('----------------------------------->>>', True)
            for l in g.graph(label='Datasets analysed at each resolution:', data=graph_data, sort=0):
                if l.startswith('#######'): continue
                pandda.log(l.replace(u"\u2588", '=').replace('= ','> '), True)
            pandda.log('----------------------------------->>>', True)
            pandda.log('')
        except ImportError: print('IMPORT ERROR (ascii_graph) - CANNOT GENERATE MAP ANALYSIS GRAPH')
        except:             pass
    else:
        pandda.log('> No Resolutions Analysed')

    pandda.log('=====================================================>>>', True)
    pandda.log('Datasets Processed: {!s}'.format(sum(resolution_counts)))
    pandda.log('Datasets Loaded {!s}'.format(pandda.datasets.size(mask_name='rejected - total', invert=True)))
    pandda.log('=====================================================>>>', True)

    return

# ============================================================================>
#
###                 PanDDA Command-Line Function
#
# ============================================================================>

def pandda_analyse_main(args):
    """Run the PANDDA algorithm, using supplied args"""

    working_phil = extract_params_default(master_phil=pandda_phil, args=args)

    try:
        # ============================================================================>
        #####
        # Initialise
        #####
        # ============================================================================
        pandda = PanddaMultiDatasetAnalyser(params=working_phil.extract())
        pandda.run_analysis_init()
        # ============================================================================>
        #####
        # Load and pre-process datasets. Generate grid and grid masks.
        #####
        # ============================================================================>
        pandda_dataset_setup(pandda=pandda)
        pandda_grid_setup(pandda=pandda)
        # ============================================================================>
        #####
        # Run the main analysis loop
        #####
        # ============================================================================>
        pandda_main_loop(pandda=pandda)
        # ============================================================================>
        #####
        # Write summaries and post-process
        #####
        # ============================================================================>
        pandda_end(pandda=pandda)
        # ============================================================================>
        #####
        # End
        #####
        # ============================================================================>
    except KeyboardInterrupt:
        raise
    except SystemExit:
        try:    pandda.log('Exited Normally')
        except: print '<<< Pandda exited before being initialised >>>'
    except:
        pandda.exit(error=True)
        raise
    else:
        pandda.exit(error=False)

    return pandda

# ============================================================================>
#
#   COMMAND LINE RUN
#
# ============================================================================>

if __name__ == '__main__':

    welcome()
    pandda = pandda_analyse_main(args=sys.argv[1:])




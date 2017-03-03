import os, sys, glob, time, gc

#################################
try:
    import matplotlib
    matplotlib.interactive(False)
    from matplotlib import pyplot
    pyplot.style.use('ggplot')
except Exception as e:
    print e
#################################

import numpy

from libtbx import easy_mp
from libtbx.utils import Sorry, Failure

from bamboo.common import Meta, Info
from bamboo.common.path import rel_symlink

from giant.jiffies import extract_params_default

from pandda import welcome
from pandda.phil import pandda_phil
from pandda.analyse.classes import PanddaMultiDatasetAnalyser, PanddaMapAnalyser, PanddaDataset
from pandda.analyse.functions import DatasetProcessor, NativeMapMaker, wrapper_run
from pandda.analyse.z_maps import PanddaZMapAnalyser
from pandda.analyse import graphs as analyse_graphs
from pandda.analyse import html as analyse_html
from pandda.misc import *


# ============================================================================>
#
###                 PanDDA Initialisation Functions
#
# ============================================================================>

def pandda_dataset_setup(pandda):
    """Initialise the pandda object and load the input datasets"""

    pandda.log('',True)
    pandda.log('################################### <~> ###################################', True)
    pandda.log('                           Running PanDDA Setup', True)
    pandda.log('################################### <~> ###################################', True)
    pandda.log('',True)

    # ============================================================================>
    # Build list of files in data directories
    # ============================================================================>
    input_files = pandda.build_input_list()
    # Check that some datasets have been found or already loaded
    if (not pandda.datasets.all()) and (not input_files):
        raise Sorry('NO DATASETS HAVE BEEN FOUND FOR ANALYSIS OR LOADED FROM PREVIOUS RUNS')
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
            if len(input_files) != pandda.args.input.max_new_datasets:
                raise Sorry('Something has gone wrong, number of selected datasets ({}) is not equal to the maximum ({})'.format(len(input_files), pandda.args.input.max_new_datasets))
        # ============================================================================>
        #####
        # Add new files and load datasets
        #####
        # ============================================================================>
        pandda.add_new_files(input_files)
        pandda.load_new_datasets()
        pandda.initialise_dataset_masks_and_tables()
        if pandda.args.method.reprocess_existing_datasets or pandda.args.method.reprocess_selected_datasets:
            pandda.check_loaded_datasets(datasets=pandda.datasets.all())
        else:
            pandda.check_loaded_datasets(datasets=pandda.datasets.mask(mask_name='old datasets', invert=True))
        # ============================================================================>
        #####
        # Set Reference Dataset
        #####
        # Select the reference dataset
        # ============================================================================>
        if not pandda.datasets.reference():
            # Filter datasets against the provided filter pdb if given
            if pandda.args.input.filter.pdb is not None:
                pandda.log('----------------------------------->>>', True)
                pandda.log('Filtering datasets against the provided pdb structure (defined by pandda.input.filter.pdb)', True)
                pandda.filter_datasets_1(filter_dataset=PanddaDataset.from_file(model_filename=pandda.args.input.filter.pdb))
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
            pandda.log('----------------------------------->>>', True)
        # ============================================================================>
        #####
        # Scale, Align and Initial-Filter All Data
        #####
        # ============================================================================>
        # Filter out datasets with different protein structures
        pandda.filter_datasets_1()
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
    # Write output csvs
    pandda.write_output_csvs()

    # Update masks
    pandda.datasets.all_masks().add_mask(name='valid - new', values=pandda.datasets.all_masks().combine_masks(
                                                                                                    names=['rejected - total', 'old datasets'],
                                                                                                    invert=[True, True],
                                                                                                    operation='and',
                                                                                                    invert_output=False))
    pandda.datasets.all_masks().add_mask(name='valid - old', values=pandda.datasets.all_masks().combine_masks(
                                                                                                    names=['rejected - total', 'old datasets'],
                                                                                                    invert=[True, False],
                                                                                                    operation='and',
                                                                                                    invert_output=False))

    return

def pandda_variation_analysis(pandda):
    """Analyse the variation in the loaded datasets"""

    pandda.log('',True)
    pandda.log('################################### <~> ###################################', True)
    pandda.log('                        Analysing dataset variation',True)
    pandda.log('################################### <~> ###################################', True)
    pandda.log('',True)

    # ============================================================================>
    #####
    # Perform variation analysis on the loaded datasets
    #####
    # ============================================================================>
    # Calculate the variation between the alignments of neighbouring C-alphas
    pandda.analyse_alignment_variation()
    pandda.analyse_overall_bfactors()

def pandda_grid_setup(pandda):
    """Generate the grid objects for the analysis"""

    pandda.log('',True)
    pandda.log('################################### <~> ###################################', True)
    pandda.log('                            Running Grid Setup',True)
    pandda.log('################################### <~> ###################################', True)
    pandda.log('',True)

    # ============================================================================>
    #####
    # Create Sampling Grid (for generated maps)
    #####
    # Create reference grid based on the reference structure
    # ============================================================================>
    if pandda.grid is None:
        # Create grid object - TODO multiple grids TODO
        pandda.create_reference_grid(
            dataset=pandda.datasets.reference(),
            grid_spacing=pandda.params.maps.grid_spacing)
        pandda.mask_reference_grid(
            dataset=pandda.datasets.reference())
        pandda.partition_reference_grid(
            dataset=pandda.datasets.reference())
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
    pandda.log('################################### <~> ###################################', True)
    pandda.log('                       Running Main PanDDA Analysis', True)
    pandda.log('################################### <~> ###################################', True)
    pandda.log('',True)

    # ============================================================================>
    # Validate/Reset the loaded datasets (those to be analysed)
    # ============================================================================>
    pandda.reset_loaded_datasets()
    # Load reflection data for each dataset
    if pandda.args.method.reprocess_existing_datasets or pandda.args.method.reprocess_selected_datasets:
        pandda.load_diffraction_data(datasets=pandda.datasets.all())
    else:
        pandda.load_diffraction_data(datasets=pandda.datasets.mask(mask_name='old datasets', invert=True))

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
                                        pandda.datasets.reference().data.summary.high_res,
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
        pandda.log('################################### <~> ###################################', True)
        pandda.log('#####                                                                 #####', True)
        pandda.log('#####{:^65}#####'.format('Processing Datasets from {!s}A -> {!s}A'.format(high_shell_limit, cut_resolution)), True)
        pandda.log('#####                                                                 #####', True)
        pandda.log('################################### <~> ###################################', True)
        pandda.log('',True)

        # ============================================================================>
        # Select datasets to build the distributions
        # ============================================================================>
        if pandda.args.method.recalculate_statistical_maps:
            pandda.log('----------------------------------->>>', True)
            pandda.log('Selecting Building Mask', True)
            pandda.log('----------------------------------->>>', True)
            building_mask_name, building_mask = pandda.select_for_building_distributions(high_res_cutoff=cut_resolution)
            # Check that we have enough datasets to build distributions
            if sum(building_mask) < pandda.params.analysis.min_build_datasets:
                # Don't have enough to construct robust distributions
                pandda.log('NOT ENOUGH DATASETS TO CALCULATE DISTRIBUTIONS ({!s}<{!s})'.format(sum(building_mask),pandda.params.analysis.min_build_datasets), True)
                pandda.log('----------------------------------->>>', True)
                continue
            else:
                # Have enough to generate robust distributions
                pandda.log('ENOUGH DATASETS -> PROCESSING THIS RESOLUTION', True)
                pandda.log('----------------------------------->>>', True)
                pandda.log('Using {} datasets for statistical electron density characterisation'.format(sum(building_mask)), True)
                pandda.log('BUILDING Datasets: {}'.format(','.join(['\n\t'*(not i%5)+d.tag for i,d in enumerate(pandda.datasets.mask(mask_name=building_mask_name))])))
        else:
            pandda.log('----------------------------------->>>', True)
            pandda.log('**NOT** Selecting Building Mask (Using Existing Statistical Maps)', True)
            # Create a dummy mask as we won't be using any datasets for building
            building_mask_name = 'dummy mask @ {}A'.format(cut_resolution)
            building_mask = [False]*pandda.datasets.size()
            pandda.datasets.all_masks().add_mask(name=building_mask_name, values=building_mask)

        # ============================================================================>
        # Select the datasets to analyse
        # ============================================================================>
        pandda.log('----------------------------------->>>', True)
        pandda.log('Selecting Analysis Mask', True)
        pandda.log('----------------------------------->>>', True)
        analysis_mask_name, analysis_mask = pandda.select_for_analysis(high_res_large_cutoff=cut_resolution, high_res_small_cutoff=high_shell_limit)
        # Check that there're some datasets to analyse
        if sum(analysis_mask) == 0:
            pandda.log('NO DATASETS TO ANALYSE @ {!s}'.format(cut_resolution))
            pandda.log('----------------------------------->>>', True)
            continue
        else:
            pandda.log('Comparing {!s} datasets against the characterised electron density'.format(sum(analysis_mask)), True)
            pandda.log('ANALYSIS Datasets: {}'.format(','.join(['\n\t'*(not i%5)+d.tag for i,d in enumerate(pandda.datasets.mask(mask_name=analysis_mask_name))])))

        # ============================================================================>
        # Combine the masks as we will need to load maps for all datasets
        # ============================================================================>
        map_load_mask = pandda.datasets.all_masks().combine_masks(names=[analysis_mask_name, building_mask_name], invert=False, operation='or', invert_output=False)
        map_load_mask_name = 'Loading @ {!s}A'.format(cut_resolution)
        pandda.datasets.all_masks().add_mask(name=map_load_mask_name, values=map_load_mask.values)

        # ============================================================================>
        #####
        # LOAD AND ANALYSE MAPS
        #####
        # ============================================================================>
        t_loop_start = time.time()
        # ==================================================>
        pandda.log('', True)
        pandda.log('----------------------------------->>>')
        pandda.log('{!s}A Analysis Started: {!s}'.format(cut_resolution, time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime(t_loop_start))), True)
        pandda.log('----------------------------------->>>')

        # ============================================================================>
        # Report
        # ============================================================================>
        pandda.log('', True)
        pandda.log('################################### <~> ###################################', True)
        pandda.log('                   Loading and aligning maps at {!s}A'.format(cut_resolution), True)
        pandda.log('################################### <~> ###################################', True)
        pandda.log('', True)
        pandda.log('----------------------------------->>>')
        pandda.log('Datasets for Building, Loading, and Analysis')
        pandda.log('Building: {!s} datasets'.format(sum(building_mask)))
        pandda.log('Analysis: {!s} datasets'.format(sum(analysis_mask)))
        pandda.log('Load Map: {!s} datasets'.format(sum(map_load_mask)))
        pandda.log('----------------------------------->>>', True)
        pandda.log('', True)

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
        pandda.truncate_scaled_data(datasets=pandda.datasets.mask(mask_name=map_load_mask_name), res_truncate=cut_resolution)
        # ============================================================================>
        # Load the reference map so that we can scale the individual maps to this
        # ============================================================================>
        ref_map = pandda.load_reference_map(map_resolution=cut_resolution)
        # ============================================================================>
        # Load the required maps
        # ============================================================================>
        map_holder_list = pandda.load_and_morph_maps(datasets=pandda.datasets.mask(mask_name=map_load_mask_name),
                                                     ref_map=ref_map, map_resolution=cut_resolution)
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
                                            meta             = Meta({'resolution'    : cut_resolution,
                                                                     'map_data_size' : pandda.grid.global_mask().outer_mask_indices().size()}),
                                            statistical_maps = statistical_maps,
                                            parent           = pandda,
                                            log              = pandda.log   )
        # ============================================================================>
        # If new statistical maps, add the meta data about which datasets are used
        # ============================================================================>
        if statistical_maps is None:
            map_analyser.statistical_maps.meta.characterisation_datasets = [d.tag for d in pandda.datasets.mask(mask_name=building_mask_name)]
        # ============================================================================>
        # Add the analysis mask to the map analyser
        # ============================================================================>
        map_analyser.dataset_maps.all_masks().add_mask(name=analysis_mask_name, values=False)
        for d in pandda.datasets.mask(mask_name=analysis_mask_name):
            map_analyser.dataset_maps.all_masks().set_value(name=analysis_mask_name, id=d.tag, value=True)
        # ============================================================================>
        # Add the building mask to the map analyser
        # ============================================================================>
        map_analyser.dataset_maps.all_masks().add_mask(name=building_mask_name, values=False)
        for d in pandda.datasets.mask(mask_name=building_mask_name):
            map_analyser.dataset_maps.all_masks().set_value(name=building_mask_name, id=d.tag, value=True)

        # ============================================================================>
        # Report
        # ============================================================================>
        pandda.log('', True)
        pandda.log('################################### <~> ###################################', True)
        pandda.log('               Performing statistical map analysis at {!s}A'.format(cut_resolution), True)
        pandda.log('################################### <~> ###################################', True)
        pandda.log('', True)

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
            map_analyser.calculate_mean_map(mask_name=building_mask_name)
        # ============================================================================>
        # If only Mean Map Requested -- exit
        # ============================================================================>
        if pandda.args.exit_flags.calculate_first_mean_map_only:
            NativeMapMaker.process(dataset  = pandda.datasets.reference(),
                                   map      = map_analyser.statistical_maps['mean_map'],
                                   filename = pandda.file_manager.get_file('mean_map').format(cut_resolution),
                                   args     = pandda.args,
                                   verbose  = pandda.settings.verbose).run()
#            write_array_as_map(grid=pandda.grid, array=map_analyser.statistical_maps.mean_map.copy().as_dense().data,
#                               f_name=pandda.file_manager.get_file('mean_map').format(cut_resolution))
            raise SystemExit('Calculating First Mean Map Only: Exiting')
        # ============================================================================>
        # Plot Mean Map against Reference Map - should be fairly similar...
        # ============================================================================>
        # Plot the mean map against the reference map (unsorted)
        analyse_graphs.mean_obs_scatter(f_name      = pandda.file_manager.get_file('ref_v_mean_map_unsort').format(cut_resolution),
                                        mean_vals   = map_analyser.statistical_maps.mean_map.as_sparse().data, obs_vals=ref_map.as_sparse().data)
        # Plot the mean map against the reference map (sorted)
        analyse_graphs.sorted_mean_obs_scatter(f_name    = pandda.file_manager.get_file('ref_v_mean_map_sort').format(cut_resolution),
                                               mean_vals = sorted(map_analyser.statistical_maps.mean_map.as_sparse().data), obs_vals=sorted(ref_map.as_sparse().data))
        # Plot the reference map distribution
        write_map_value_distribution(map_vals       = ref_map.as_sparse().data,
                                     output_file    = pandda.file_manager.get_file('ref_map_dist').format(cut_resolution))
        # ============================================================================>
        # Calculate the uncertainty of all loaded maps (needs the mean map to have been calculated)
        # ============================================================================>
        inner_mask_idxs_reindx = pandda.grid.index_on_other(query=pandda.grid.global_mask().inner_mask_indices(), other=pandda.grid.global_mask().outer_mask_indices())
        assert len(inner_mask_idxs_reindx) == len(pandda.grid.global_mask().inner_mask_indices())
        map_analyser.calculate_map_uncertainties(masked_idxs=inner_mask_idxs_reindx, cpus=pandda.settings.cpus)
        # ============================================================================>
        # Plot uncertainties of maps
        # ============================================================================>
        try:
            from ascii_graph import Pyasciigraph
            g=Pyasciigraph(float_format='{0:.3f}')
            graph_data = [(m.meta.tag, round(m.meta.map_uncertainty,3)) for m in map_analyser.dataset_maps.all()]
            pandda.log('----------------------------------->>>', True)
            for l in g.graph(label='Sorted Map Uncertainties (Ascending Order)', data=graph_data):
                if l.startswith('#######'): continue
                pandda.log(l.replace(u"\u2588", '=').replace('= ','> '), True)
            pandda.log('----------------------------------->>>', True)
            pandda.log('')
        except ImportError:
            print('IMPORT ERROR (ascii_graph) - CANNOT GENERATE UNCERTAINTY GRAPH')
        except:
            raise
        # ============================================================================>
        # Calculate the statistics of the maps
        # ============================================================================>
        if (not map_analyser.statistical_maps.sadj_map) or pandda.args.method.recalculate_statistical_maps:
            map_analyser.calculate_statistical_maps(mask_name=building_mask_name, cpus=pandda.settings.cpus)

        # ============================================================================>
        # PICKLE THE DATASETS THAT HAVE JUST BEEN PROCESSED
        # ============================================================================>
#        pandda.sync_datasets(datasets=pandda.datasets.mask(mask_name=analysis_mask_name), overwrite_dataset_meta=True)
#        pandda.pickle_the_pandda(components=['datasets'], datasets=pandda.datasets.mask(mask_name=analysis_mask_name))

        # ============================================================================>
        # Write the statistical maps
        # ============================================================================>
        pandda.write_map_analyser_maps(map_analyser=map_analyser)
        # ============================================================================>
        # Pickle the statistical map objects for re-use
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

        # TODO TODO TODO
        # ============================================================================>
        # Write out Grid Point Distributions for interesting grid points (high modality, etc...)
        # ============================================================================>
        # map_analyser.analyse_point_distributions()
        # TODO TODO TODO

        # ============================================================================>
        # Report
        # ============================================================================>
        pandda.log('', True)
        pandda.log('################################### <~> ###################################', True)
        pandda.log('               Preparing to analyse {!s} datasets at {!s}A'.format(pandda.datasets.size(mask_name=analysis_mask_name), cut_resolution), True)
        pandda.log('################################### <~> ###################################', True)
        pandda.log('', True)

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
        dummy_map_analyser = PanddaMapAnalyser(dataset_maps=None,
                                               meta=Meta({'resolution':cut_resolution,
                                                          'map_data_size':pandda.grid.global_mask().outer_mask_indices().size()}),
                                               statistical_maps=map_analyser.statistical_maps,
                                               parent=None, log=None)
        # ============================================================================>
        # Blob Search Object
        # ============================================================================>
        dummy_blob_finder = PanddaZMapAnalyser(params=pandda.params.blob_search, grid=pandda.grid, log=pandda.log)
        dummy_blob_finder.print_settings()

        # ============================================================================>
        #####
        # Calculate Z-Maps
        #####
        # ============================================================================>
        t_anal_start = time.time()
        # ==================================================>

        dataset_processor_list = []

        # Iterate through and prepare to calculate z-maps
        for i_d, dataset in enumerate(pandda.datasets.mask(mask_name=analysis_mask_name)):

            # ============================================================================>
            # Record which resolution this dataset was analysed at
            # ============================================================================>
            resolution_count[cut_resolution].append(dataset.tag)

            # ============================================================================>
            # Validate/Check/Zero the dataset
            # ============================================================================>
            # Dataset should not have any events
            assert dataset.events == []
            # Update datasets masks flag - for this analysis
            pandda.datasets.all_masks().set_value(name='analysed', id=dataset.tag, value=True)
            # Update the dataset meta object -- this is persistent
            dataset.meta.analysed = True
            # ============================================================================>
            # Extract the map for this dataset
            # ============================================================================>
            dataset_map = map_analyser.dataset_maps.get(tag=dataset.tag)
            # ============================================================================>
            # Compile arguments for this datasets
            # ============================================================================>
            dp = DatasetProcessor(dataset=dataset, dataset_map=dataset_map.as_sparse(),
                                  grid=pandda.grid, map_analyser=dummy_map_analyser,
                                  args=pandda.args, verbose=pandda.settings.verbose)
            dataset_processor_list.append(dp)

        pandda.log('', True)
        pandda.log('################################### <~> ###################################', True)
        pandda.log('{:^75}'.format('Calculating Z-maps and event maps for {!s} dataset(s) at {!s}A'.format(pandda.datasets.size(mask_name=analysis_mask_name), cut_resolution)), True)
        pandda.log('################################### <~> ###################################', True)
        pandda.log('', True)

        proc_results = easy_mp.pool_map(func=wrapper_run, args=dataset_processor_list, processes=pandda.settings.cpus, chunksize=1)

        pandda.log('----------------------------------->>>', True)
        pandda.log('Updating with results from analysing {!s} Dataset(s) at {!s}A'.format(pandda.datasets.size(mask_name=analysis_mask_name), cut_resolution), True)

        for results in proc_results:

            # Unpack results
            dataset, dataset_meta, log_strs = results

            pandda.log('', True)
            pandda.log('======================================>>>', True)
            pandda.log('Z-Map Analysis Results for {}'.format(dataset.tag), True)
            pandda.log('======================================>>>', True)
            pandda.log('\n'.join(log_strs), True)

            # ============================================================================>
            # STORE ANALYSIS DATA IN DATASET MAP TABLE
            # ============================================================================>
            # Add to the dataset map summary table
            pandda.tables.dataset_map_info.set_value(dataset.tag, 'analysed_resolution', dataset_meta.resolution)
            pandda.tables.dataset_map_info.set_value(dataset.tag, 'map_uncertainty',     round(dataset_meta.map_uncertainty,3))
            pandda.tables.dataset_map_info.set_value(dataset.tag, 'obs_map_mean',        round(dataset_meta.obs_map_mean,3))
            pandda.tables.dataset_map_info.set_value(dataset.tag, 'obs_map_rms',         round(dataset_meta.obs_map_rms,3))
            pandda.tables.dataset_map_info.set_value(dataset.tag, 'scl_map_mean',        round(dataset_meta.scl_map_mean,3))
            pandda.tables.dataset_map_info.set_value(dataset.tag, 'scl_map_rms',         round(dataset_meta.scl_map_rms,3))
            pandda.tables.dataset_map_info.set_value(dataset.tag, 'z_map_mean',          round(dataset_meta.z_mean,3))
            pandda.tables.dataset_map_info.set_value(dataset.tag, 'z_map_std',           round(dataset_meta.z_stdv,3))
            pandda.tables.dataset_map_info.set_value(dataset.tag, 'z_map_skew',          round(dataset_meta.z_skew,3))
            pandda.tables.dataset_map_info.set_value(dataset.tag, 'z_map_kurt',          round(dataset_meta.z_kurt,3))

            # ============================================================================>
            # WRITE OUT DATASET INFORMATION TO CSV FILE
            # ============================================================================>
            out_list = pandda.tables.dataset_info.loc[dataset.tag].append(pandda.tables.dataset_map_info.loc[dataset.tag])
            out_list.to_csv(path=dataset.file_manager.get_file('dataset_info'), header=True, index_label='dtag')

            # ============================================================================>
            # Mark as interesting and add events to the event table
            # ============================================================================>
            if dataset.events:
#                # Create a link to the interesting directories in the initial results directory
#                hit_dir = os.path.join(pandda.file_manager.get_dir('interesting_datasets'), dataset.tag)
#                if not os.path.exists(hit_dir): rel_symlink(orig=dataset.file_manager.get_dir('root'), link=hit_dir)
                pandda.datasets.all_masks().set_value(name='interesting', id=dataset.tag, value=True)
                for e in dataset.events:
                    pandda.add_event_to_event_table(dataset=dataset, event=e)

            # ============================================================================>
            # Update the master copy of the dataset object
            # ============================================================================>
            master_dataset = pandda.datasets.get(tag=dataset.tag)
            master_dataset.events = dataset.events
            master_dataset.child  = map_analyser.dataset_maps.get(tag=dataset.tag)

        # ============================================================================>
        #####
        # Generate native-aligned maps (in the crystallographic unit cell)
        #####
        # ============================================================================>
        native_map_maker_list = []
        for i_d, dataset in enumerate(pandda.datasets.mask(mask_name=analysis_mask_name)):
            # ============================================================================>
            # Make Z-map for each dataset (if events or write_z_maps_for_all_datasets)
            # ============================================================================>
            if (pandda.args.output.maps.write_z_maps=='interesting' and dataset.events) or (pandda.args.output.maps.write_z_maps=='all'):
                ref_z_map = map_analyser.calculate_z_map(map         = dataset.child,
                                                         uncertainty = dataset.child.meta.map_uncertainty,
                                                         method      = pandda.args.params.z_map.map_type)
                ref_z_map = (ref_z_map - numpy.mean(ref_z_map.data)) * (1.0/numpy.std(ref_z_map.data))
                map_maker = NativeMapMaker(dataset  = dataset,
                                           map      = ref_z_map,
                                           filename = dataset.file_manager.get_file('native_z_map'),
                                           args     = pandda.args,
                                           verbose  = pandda.settings.verbose)
                native_map_maker_list.append(map_maker)
            # ============================================================================>
            # Make Event-map for each event
            # ============================================================================>
            for i,e in enumerate(dataset.events):
                ref_event_map = ( dataset.child - map_analyser.statistical_maps.mean_map * e.info.estimated_bdc )
                map_maker = NativeMapMaker(dataset  = dataset,
                                           map      = ref_event_map,
                                           filename = dataset.file_manager.get_file('native_event_map').format(e.id[1], 1-e.info.estimated_bdc),
                                           args     = pandda.args,
                                           verbose  = pandda.settings.verbose)
                native_map_maker_list.append(map_maker)
            # ============================================================================>
            # Mean (ground-state) map for this resolution
            # ============================================================================>
            if (pandda.args.output.maps.write_mean_map=='interesting' and dataset.events) or (pandda.args.output.maps.write_mean_map=='all'):
                map_maker = NativeMapMaker(dataset  = dataset,
                                           map      = map_analyser.statistical_maps.mean_map,
                                           filename = dataset.file_manager.get_file('native_mean_map'),
                                           args     = pandda.args,
                                           verbose  = pandda.settings.verbose)
                native_map_maker_list.append(map_maker)
        # ============================================================================>
        # Statistical maps in the native frame of the reference dataset
        # ============================================================================>
        if pandda.args.output.maps.write_statistical_maps and pandda.args.method.recalculate_statistical_maps:
            for m in ['mean_map', 'medn_map','stds_map','sadj_map','skew_map','kurt_map','bimo_map']:
                map_maker = NativeMapMaker(dataset  = pandda.datasets.reference(),
                                           map      = map_analyser.statistical_maps[m],
                                           filename = pandda.file_manager.get_file(m).format(cut_resolution),
                                           args     = pandda.args,
                                           verbose  = pandda.settings.verbose)
                native_map_maker_list.append(map_maker)
        # ============================================================================>
        # Write the compiled list of maps
        # ============================================================================>
        if native_map_maker_list:
            pandda.log('', True)
            pandda.log('################################### <~> ###################################', True)
            pandda.log('{:^75}'.format('Outputting {!s} native maps at {!s}A'.format(len(native_map_maker_list), cut_resolution)), True)
            pandda.log('################################### <~> ###################################', True)
            pandda.log('', True)
            for n in native_map_maker_list:
                pandda.log('{:<30} -> {}'.format(n.data[0].tag, os.path.split(n.data[2])[-1]), True)
            proc_results = easy_mp.pool_map(func=wrapper_run, args=native_map_maker_list, processes=pandda.settings.cpus, chunksize=1)

        # ============================================================================>
        #####
        # Write the summaries
        #####
        # ============================================================================>
        pandda.sync_datasets(datasets=pandda.datasets.mask(mask_name=analysis_mask_name), overwrite_dataset_meta=True)
        # Write output graphs for map analysis
        if pandda.settings.plot_graphs:

            pandda.log('', True)
            pandda.log('################################### <~> ###################################', True)
            pandda.log('                     Writing Statistical Maps Summary'.format(cut_resolution), True)
            pandda.log('################################### <~> ###################################', True)
            pandda.log('', True)

            analyse_graphs.write_map_analyser_graphs(pandda=pandda, map_analyser=map_analyser, analysis_mask_name=analysis_mask_name)

        # ============================================================================>
        #####
        # Print summaries
        #####
        # ============================================================================>
        t_loop_end = time.time()
        # ==================================================>
        pandda.log('', True)
        pandda.log('################################### <~> ###################################', True)
        pandda.log('                     Finished Analysis at {!s}A'.format(cut_resolution), True)
        pandda.log('################################### <~> ###################################', True)
        pandda.log('', True)
        pandda.log('----------------------------------->>>', True)
        pandda.log('{!s}A Z-Map Processing Time: {!s}'.format(cut_resolution, time.strftime("%H hours:%M minutes:%S seconds", time.gmtime(t_loop_end - t_anal_start))), True)
        pandda.log('{!s}A Total Processing Time: {!s}'.format(cut_resolution, time.strftime("%H hours:%M minutes:%S seconds", time.gmtime(t_loop_end - t_loop_start))), True)
        pandda.log('----------------------------------->>>', True)
        pandda.log('@{!s}A:\t {!s}/{!s} New Datasets Analysed'.format(cut_resolution,
                                                                      pandda.datasets.all_masks().combine_masks(
                                                                                                    names=[analysis_mask_name, 'old datasets'],
                                                                                                    invert=[False,True],
                                                                                                    operation='and',
                                                                                                    invert_output=False).sum(),
                                                                      pandda.datasets.size(mask_name='valid - new')))
        pandda.log('@{!s}A:\t {!s}/{!s} Old Datasets Analysed'.format(cut_resolution,
                                                                      pandda.datasets.all_masks().combine_masks(
                                                                                                    names=[analysis_mask_name, 'old datasets'],
                                                                                                    invert=[False,False],
                                                                                                    operation='and',
                                                                                                    invert_output=False).sum(),
                                                                      pandda.datasets.size(mask_name='valid - old')))
        pandda.log('Total:\t {!s}/{!s} Datasets Analysed'.format(pandda.datasets.size(mask_name='analysed'), pandda.datasets.size(mask_name='rejected - total', invert=True)))
        pandda.log('----------------------------------->>>', True)

        # ============================================================================>
        #####
        # PICKLE THE DATASETS THAT HAVE JUST BEEN PROCESSED
        #####
        # ============================================================================>
        # Clear the linked maps if requested
        if not pandda.args.output.pickling.pickle_dataset_maps:
            for dataset in pandda.datasets.mask(mask_name=analysis_mask_name):
                dataset.child = None
        # Pickle the datasets
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
        for dataset in pandda.datasets.mask(mask_name=map_load_mask_name):
            dataset.child = None
        # No maps should still be loaded for memory consumption reasons
        for dataset in pandda.datasets.all():
            assert dataset.child is None
        # ============================================================================>
        # Launch Garbage-Collection Manually just to be sure
        # ============================================================================>
        gc.collect(); gc.collect(); gc.collect()

        # ============================================================================>
        #####
        # LIVE ANALYSIS - RUN AT THE END OF EACH LOOP
        #####
        # Analyse the processed data
        # ============================================================================>
        pandda.log('', True)
        pandda.log('################################### <~> ###################################', True)
        pandda.log('       Clustering identified events and updating central data tables', True)
        pandda.log('################################### <~> ###################################', True)
        pandda.log('', True)

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
    pandda.log('', True)
    pandda.log('################################### <~> ###################################', True)
    pandda.log('                             ANALYSIS COMPLETE', True)
    pandda.log('################################### <~> ###################################', True)
    pandda.log('', True)

    pandda.log('Total Analysis Time: {!s}'.format(time.strftime("%H hours:%M minutes:%S seconds", time.gmtime(t_analysis_end - t_analysis_start))), True)

    return

# ============================================================================>
#
###                 PanDDA Output + Wrap-Up Functions
#
# ============================================================================>

def pandda_end(pandda):
    """Write the output summaries and close"""

    pandda.log('',True)
    pandda.log('################################### <~> ###################################', True)
    pandda.log('                  Writing PanDDA End-of-Analysis Summary', True)
    pandda.log('################################### <~> ###################################', True)
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
    pandda.log('----------------------------------->>>', False)
    if event_num:
        pandda.log('')
        pandda.log('----------------------------------->>>', True)
        pandda.log('Potentially Useful Shortcuts for Future Runs', True)
        pandda.log('---------->>>', True)
        pandda.log('All datasets with events:')
        pandda.log('no_build={!s}'.format(','.join(zip(*event_num)[0])), True)
        for site_num, event_list in pandda.tables.event_info.groupby('site_idx').groups.items():
            pandda.log('---------->>>', True)
            pandda.log('Datasets with events at Site {}'.format(site_num))
            pandda.log('no_build={!s}'.format(','.join(zip(*event_list)[0])))
        pandda.log('----------------------------------->>>', False)

    pandda.log('')
    pandda.log('----------------------------------->>>', True)
    pandda.log('Lists of datasets used to generate the statistical maps', True)
    for res in sorted(pandda.stat_maps.get_resolutions()):
        sm = pandda.stat_maps.get(res)
        pandda.log('---------->>>', True)
        pandda.log('Statistical Electron Density Characterisation at {}A'.format(res))
        pandda.log('> Density characterised using {} datasets'.format(len(sm.meta.characterisation_datasets)))
        pandda.log('> Dataset IDs: {}'.format(','.join(['\n\t'*(not i%5)+d for i,d in enumerate(sm.meta.characterisation_datasets)])))

    # ============================================================================>
    #####
    # SUMMARIES ------------------------------>>>
    #####
    # ============================================================================>

    pandda.log('----------------------------------->>>', True)
    pandda.log('Writing final output files', True)
    pandda.write_output_csvs()
    analyse_html.write_analyse_html(pandda)
    pandda.log('----------------------------------->>>', True)

    # ============================================================================>
    # SCREEN GRAPHS -------------------------->>>
    # ============================================================================>

    resolution_counts = pandda.tables.dataset_map_info['analysed_resolution'].value_counts().sort_index()
    graph_data = [(str(r), c) for r,c in resolution_counts.iteritems()]
    if graph_data:
        try:
            from ascii_graph import Pyasciigraph
            g=Pyasciigraph()
            for l in g.graph(label='Datasets analysed at each resolution:', data=graph_data):
                if l.startswith('#######'): continue
                pandda.log(l.replace(u"\u2588", '=').replace('= ','> '), True)
        except ImportError: print('IMPORT ERROR (ascii_graph) - CANNOT GENERATE MAP ANALYSIS GRAPH')
        except:             raise
        pandda.log('----------------------------------->>>', True)
        pandda.log('')
    else:
        pandda.log('> No Resolutions Analysed')

    pandda.log('')
    pandda.log('----------------------------------->>>', True)
    pandda.log('Datasets Processed: {!s}'.format(sum(resolution_counts)))
    pandda.log('Datasets Loaded {!s}'.format(pandda.datasets.size(mask_name='rejected - total', invert=True)))
    pandda.log('----------------------------------->>>', True)
    pandda.log('')

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




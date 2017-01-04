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

from libtbx import easy_mp
from libtbx.utils import Sorry, Failure

from bamboo.common import Meta, Info
from bamboo.common.path import rel_symlink

from giant.jiffies import extract_params_default

from pandda import welcome
from pandda.phil import pandda_phil
from pandda.analyse.classes import PanddaMultiDatasetAnalyser, PanddaMapAnalyser, PanddaDataset
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
    pandda.log('=====================================================>>>', True)
    pandda.log('Running PanDDA Setup', True)
    pandda.log('=====================================================>>>', True)
    pandda.log('',True)

    # ============================================================================>
    # Build list of files in data directories
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
        # ============================================================================>
        #####
        # Scale, Align and Initial-Filter All Data
        #####
        # ============================================================================>
        # Filter out datasets with different protein structures
        pandda.filter_datasets_1()
        # Load reflection data for each dataset
        pandda.load_diffraction_data()
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
    #####
    # Create Sampling Grid (for generated maps)
    #####
    # Create reference grid based on the reference structure
    # ============================================================================>
    if pandda.grid is None:
        # Create grid object - TODO multiple grids TODO
        pandda.create_reference_grid(dataset=pandda.datasets.reference(), grid_spacing=pandda.params.maps.grid_spacing)
        # Create various masks to define regions of the grid by distance to the protein and symmetry copies
        pandda.mask_reference_grid(dataset=pandda.datasets.reference())
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
            building_mask = pandda.select_for_building_distributions(high_res_cutoff=cut_resolution, building_mask_name=building_mask_name)
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
        analysis_mask = pandda.select_for_analysis(high_res_large_cutoff=cut_resolution, high_res_small_cutoff=high_shell_limit, analysis_mask_name=analysis_mask_name)
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
        pandda.truncate_scaled_data(datasets=pandda.datasets.mask(mask_name=map_load_mask_name))
        # ============================================================================>
        # Load the reference map so that we can scale the individual maps to this
        # ============================================================================>
        ref_map = pandda.load_reference_map(map_resolution=cut_resolution)
        # ============================================================================>
        # Load the required maps
        # ============================================================================>
        map_holder_list = pandda.load_and_morph_maps(datasets=pandda.datasets.mask(mask_name=map_load_mask_name), ref_map=ref_map, map_resolution=cut_resolution)
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
        # Add the analysis mask to the map analyser
        # ============================================================================>
        map_analyser.dataset_maps.all_masks().add_mask(mask_name=analysis_mask_name, mask=[False]*map_analyser.dataset_maps.size())
        for d in pandda.datasets.mask(mask_name=analysis_mask_name):
            map_analyser.dataset_maps.all_masks().set_mask_value(mask_name=analysis_mask_name, entry_id=d.tag, value=True)
        # ============================================================================>
        # Add the building mask to the map analyser
        # ============================================================================>
        map_analyser.dataset_maps.all_masks().add_mask(mask_name=building_mask_name, mask=[False]*map_analyser.dataset_maps.size())
        for d in pandda.datasets.mask(mask_name=building_mask_name):
            map_analyser.dataset_maps.all_masks().set_mask_value(mask_name=building_mask_name, entry_id=d.tag, value=True)

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
            write_array_to_map(output_file=pandda.file_manager.get_file('mean_map').format(cut_resolution),
                            map_data=map_analyser.statistical_maps.mean_map, grid=pandda.grid)
            raise SystemExit('Calculating First Mean Map Only: Exiting')
        # ============================================================================>
        # Plot Mean Map against Reference Map - should be fairly similar...
        # ============================================================================>
        # Plot the mean map against the reference map (unsorted)
        analyse_graphs.mean_obs_scatter(f_name=os.path.join(pandda.file_manager.get_dir('reference'), 'reference_against_mean_unsorted.png'),
                                        mean_vals=map_analyser.statistical_maps.mean_map.data, obs_vals=ref_map.data)
        # Plot the mean map against the reference map (sorted)
        analyse_graphs.sorted_mean_obs_scatter(f_name=os.path.join(pandda.file_manager.get_dir('reference'), 'reference_against_mean_sorted.png'),
                                        mean_vals=sorted(map_analyser.statistical_maps.mean_map.data), obs_vals=sorted(ref_map.data))
        # Plot the reference map distribution
        write_map_value_distribution(map_vals=ref_map.data, output_file=os.path.join(pandda.file_manager.get_dir('reference'), 'reference_map_distribution.png'))
        # ============================================================================>
        # Calculate the uncertainty of all loaded maps (needs the mean map to have been calculated)
        # ============================================================================>
        # Reindex the inner mask indices on the outer mask indices
        outer_mask_idxs = list(pandda.grid.global_mask().outer_mask_indices())
        inner_mask_idxs = list(pandda.grid.global_mask().inner_mask_indices())
        tmp_inner_mask_idxs = [outer_mask_idxs.index(i) for i in inner_mask_idxs]
        # If not pandda.args.method.recalculate_statistical_maps, then no building_mask datasets will have been added, so don't need to be selective in this function
        map_analyser.calculate_map_uncertainties(masked_idxs=tmp_inner_mask_idxs, cpus=pandda.settings.cpus)
        # ============================================================================>
        # Plot uncertainties of maps
        # ============================================================================>
        try:
            from ascii_graph import Pyasciigraph
            g=Pyasciigraph()
            graph_data = [(m.tag, round(m.meta.map_uncertainty,3)) for m in map_analyser.dataset_maps.all()]
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
            map_analyser.calculate_statistical_maps(mask_name=building_mask_name, cpus=pandda.settings.cpus)

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
        dummy_map_analyser = PanddaMapAnalyser(dataset_maps=None, meta=Meta({'resolution':cut_resolution,'map_data_size':pandda.grid.global_mask().outer_mask_indices().size()}),
                                                statistical_maps=map_analyser.statistical_maps, parent=None, log=None)
        # ============================================================================>
        # Blob Search Object
        # ============================================================================>
        blob_finder = PanddaZMapAnalyser(params=pandda.params.blob_search, grid_spacing=pandda.reference_grid().grid_spacing(), log=pandda.log)
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

        dataset_processor_list = []

        # Iterate through and prepare to calculate z-maps
        for i_d, dataset in enumerate(pandda.datasets.mask(mask_name=analysis_mask_name)):

            # ============================================================================>
            # Record the which resolution this dataset was analysed at
            # ============================================================================>
            resolution_count[cut_resolution].append(dataset.tag)

            # ============================================================================>
            # Validate/Check/Zero the dataset
            # ============================================================================>
            # Dataset should not have any events
            assert dataset.events == []
            # Update datasets masks flag - for this analysis
            pandda.datasets.all_masks().set_mask_value(mask_name='analysed', entry_id=dataset.tag, value=True)
            # Update the dataset meta object -- this is persistent
            dataset.meta.analysed = True
            # ============================================================================>
            # Extract the map for this dataset
            # ============================================================================>
            dataset_map = map_analyser.dataset_maps.get(tag=dataset.tag)
            # ============================================================================>
            # Compile arguments for this datasets
            # ============================================================================>
            dataset_processor = DatasetProcessor(dataset=dataset, dataset_map=dataset_map, grid=pandda.grid, reference_dataset=pandda.datasets.reference(),
                                            map_analyser=dummy_map_analyser, args=pandda.args, verbose=pandda.settings.verbose)

            dataset_processor_list.append(dataset_processor)

            ### XXX ------- preprocessing for the map func above here ------- XXX ###

        pandda.log('----------------------------------->>>', True)
        pandda.log('Calculating and Analysing Z-Maps for {!s} Dataset(s) at {!s}A'.format(pandda.datasets.size(mask_name=analysis_mask_name), cut_resolution), True)

        proc_results = easy_mp.pool_map(func=wrapper_run, args=proc_dataset_inp_list, processes=pandda.settings.cpus, maxtasksperchild=1)

        pandda.log('----------------------------------->>>', True)
        pandda.log('Updating with results from analysing {!s} Dataset(s) at {!s}A'.format(pandda.datasets.size(mask_name=analysis_mask_name), cut_resolution), True)

        for results in proc_results:

            # Unpack results
            dataset, dataset_map, log_strs = results

            pandda.log('', True)
            pandda.log('======================================>>>', True)
            pandda.log('Z-Map Analysis Results for {}'.format(dataset.tag), True)
            pandda.log('======================================>>>', True)
            pandda.log('\n'.join(log_strs), True)

            # ============================================================================>
            # STORE ANALYSIS DATA IN DATASET MAP TABLE
            # ============================================================================>
            # Add to the dataset map summary table
            pandda.tables.dataset_map_info.set_value(dataset.tag, 'analysed_resolution', dataset_map.meta.resolution)
            pandda.tables.dataset_map_info.set_value(dataset.tag, 'map_uncertainty',     round(dataset_map.meta.map_uncertainty,3))
            pandda.tables.dataset_map_info.set_value(dataset.tag, 'obs_map_mean',        round(dataset_map.meta.obs_map_mean,3))
            pandda.tables.dataset_map_info.set_value(dataset.tag, 'obs_map_rms',         round(dataset_map.meta.obs_map_rms,3))
            pandda.tables.dataset_map_info.set_value(dataset.tag, 'z_map_mean',          round(dataset_map.meta.z_mean,3))
            pandda.tables.dataset_map_info.set_value(dataset.tag, 'z_map_std',           round(dataset_map.meta.z_stdv,3))
            pandda.tables.dataset_map_info.set_value(dataset.tag, 'z_map_skew',          round(dataset_map.meta.z_skew,3))
            pandda.tables.dataset_map_info.set_value(dataset.tag, 'z_map_kurt',          round(dataset_map.meta.z_kurt,3))

            # ============================================================================>
            # WRITE OUT DATASET INFORMATION TO CSV FILE
            # ============================================================================>
            out_list = pandda.tables.dataset_info.loc[dataset.tag].append(pandda.tables.dataset_map_info.loc[dataset.tag])
            out_list.to_csv(path=dataset.file_manager.get_file('dataset_info'), header=True, index_label='dtag')

            # ============================================================================>
            # Link outputs for datasets with hits
            # ============================================================================>
            if dataset.events:
                # Create a link to the interesting directories in the initial results directory
                hit_dir = os.path.join(pandda.file_manager.get_dir('interesting_datasets'), dataset.tag)
                if not os.path.exists(hit_dir): rel_symlink(orig=dataset.file_manager.get_dir('root'), link=hit_dir)
                # ============================================================================>
                # Add event to the event table
                # ============================================================================>
                for e in dataset.events:
                    pandda.add_event_to_event_table(dataset=dataset, event=e)

            # ============================================================================>
            # Update the master copy of the dataset object
            # ============================================================================>
            master_dataset = pandda.datasets.get(tag=dataset.tag)
            master_dataset.events = dataset.events
            master_dataset.child  = dataset_map

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
            for dataset in pandda.datasets.mask(mask_name=analysis_mask_name):
                dataset.child = None
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
        for dataset in pandda.datasets.mask(mask_name=analysis_mask_name):
            dataset.child = None
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




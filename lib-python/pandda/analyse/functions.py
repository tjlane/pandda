from __future__ import print_function

import os, sys, copy

import numpy

from scitbx.array_family import flex
from scitbx.math import basic_statistics

from libtbx.utils import Sorry, Failure

from bamboo.common import Meta, Info
from bamboo.common.logs import Log
from bamboo.common.path import rel_symlink

from giant.grid.masks import GridMask
from giant.stats import quantile_quantile_scaling
from giant.stats.ospina import estimate_true_underlying_sd
from giant.structure.select import protein, find_nearest_atoms
from giant.xray.maps import scale_map_to_reference
from giant.xray.maps.bdc import calculate_bdc_subtracted_map, estimate_event_background_correction
from giant.xray.local_align_maps import create_native_map

from pandda import welcome
from pandda.phil import pandda_phil
from pandda.analyse.z_maps import PanddaZMapAnalyser
from pandda.analyse.events import PointCluster, Event
from pandda.analyse import graphs as analyse_graphs
from pandda.misc import *

def wrapper_run(c):
    if c is not None:
        return c.run()
    else:
        return c

# ============================================================================>


class DatasetProcesser(object):


    def __init__(self, dataset, dataset_map, grid, reference_dataset, map_analyser, args, verbose):
        """
        The main object for comparing each dataset-map to the ensemble-maps.
        Constructed so that the class can be initialised and then called within a multiprocessing function.

        Usage:
            new_proc = DatasetProcessor(...)
            output   = new_proc.run()
        or:
            output   = DatasetProcessor.process(...)
        """

        self.data = (dataset,dataset_map,reference_dataset,grid,map_analyser,args,verbose)

    @classmethod
    def process(cls, dataset, dataset_map, grid, reference_dataset, map_analyser, args, verbose):
        """Process the dataset immediately and return output"""
        return cls(dataset=dataset, dataset_map=dataset_map, reference_dataset=reference_dataset,
                    grid=grid, map_analyser=map_analyser, args=args, verbose=verbose).run()

    def run(self):
        """Process the dataset"""

        dataset,dataset_map,reference_dataset,grid,map_analyser,args,verbose = self.data

        # ============================================================================>
        # Prepare output objects
        # ============================================================================>
        log_strs = []

        # ============================================================================>
        # Build new blob search object
        # ============================================================================>
        print('Creating ZMapAnalyser for Dataset {!s}'.format(dataset.tag))

        blob_finder = PanddaZMapAnalyser(params=args.params.blob_search, grid=grid,
                            log=Log(log_file=dataset.file_manager.get_file('dataset_log'), verbose=False))

        print('Writing Log for Dataset {!s} to '.format(dataset.tag, dataset.file_manager.get_file('dataset_log')))

        # ============================================================================>
        # Generate masks for this dataset
        # ============================================================================>
        print('Generating Symmetry Contacts for Dataset {!s}'.format(dataset.tag))

        # Generate symmetry contacts for this dataset and align to reference frame
        dataset_sym_copies = dataset.model.crystal_contacts(distance_cutoff=args.params.masks.outer_mask+5, combine_copies=True)
        dataset_sym_copies.atoms().set_xyz(dataset.model.alignment.nat2ref(dataset_sym_copies.atoms().extract_xyz()))
        dataset_sym_copies.write_pdb_file(dataset.file_manager.get_file('symmetry_copies'))

        # Extract protein atoms from the symmetry copies
        dataset_sym_sites_cart = protein(dataset_sym_copies).atoms().extract_xyz()

        # Generate custom grid mask for this dataset
        dataset_mask = GridMask(
                            parent=grid, sites_cart=dataset_sym_sites_cart,
                            max_dist=args.params.masks.outer_mask,
                            min_dist=args.params.masks.inner_mask_symmetry)
        # TODO REIMPLEMENT TODO #
        # Combine the standard mask with the custom mask
        grid_idxr = grid.indexer()
        dataset_total_mask = [gp for gp in grid.global_mask().total_mask() if dataset_mask.inner_mask_binary()[grid_idxr(gp)] < 0.5]
        dataset_outer_mask = [gp for gp in grid.global_mask().outer_mask() if dataset_mask.inner_mask_binary()[grid_idxr(gp)] < 0.5]
        # TODO REIMPLEMENT TODO #

        # ============================================================================>
        #####
        # CALCULATE Z-MAPS AND LOOK FOR LARGE BLOBS
        #####
        # ============================================================================>
        print('Calculating Z-MAPs for Dataset {!s}'.format(dataset.tag))

        # ============================================================================>
        # CALCULATE MEAN-DIFF MAPS
        # ============================================================================>
        assert dataset_map.data is not None, 'Something has gone wrong - this dataset has no loaded map'
        mean_diff_map = dataset_map - map_analyser.statistical_maps.mean_map.data
        # ============================================================================>
        # NAIVE Z-MAP - NOT USING UNCERTAINTY ESTIMATION OR ADJUSTED STDS
        # ============================================================================>
        z_map_naive = map_analyser.calculate_z_map(map=dataset_map, method='naive')
        z_map_naive_normalised = (z_map_naive - numpy.mean(z_map_naive.data)) / numpy.std(z_map_naive.data)
        # ============================================================================>
        # UNCERTAINTY Z-MAP - NOT USING ADJUSTED STDS
        # ============================================================================>
        z_map_uncty = map_analyser.calculate_z_map(map=dataset_map, uncertainty=dataset_map.meta.map_uncertainty, method='uncertainty')
        z_map_uncty_normalised = (z_map_uncty - numpy.mean(z_map_uncty.data)) / numpy.std(z_map_uncty.data)
        # ============================================================================>
        # ADJUSTED+UNCERTAINTY Z-MAP
        # ============================================================================>
        z_map_compl = map_analyser.calculate_z_map(map=dataset_map, uncertainty=dataset_map.meta.map_uncertainty, method='adjusted+uncertainty')
        z_map_compl_normalised = (z_map_compl - numpy.mean(z_map_compl.data)) / numpy.std(z_map_compl.data)

        # ============================================================================>
        # SELECT WHICH MAP TO DO THE BLOB SEARCHING ON
        # ============================================================================>
        if   args.params.z_map.map_type == 'naive':                z_map = z_map_naive_normalised
        elif args.params.z_map.map_type == 'adjusted':             z_map = z_map_adjst_normalised
        elif args.params.z_map.map_type == 'uncertainty':          z_map = z_map_uncty_normalised
        elif args.params.z_map.map_type == 'adjusted+uncertainty': z_map = z_map_compl_normalised
        else: raise Exception('Invalid Z-map type')

        # ============================================================================>
        # RECORD Z-MAP FOR STATISTICS
        # ============================================================================>
        # Calculate statistics of z-maps
        z_map_stats = basic_statistics(flex.double(z_map.data))
        dataset_map.meta.z_mean = z_map_stats.mean
        dataset_map.meta.z_stdv = z_map_stats.bias_corrected_standard_deviation
        dataset_map.meta.z_skew = z_map_stats.skew
        dataset_map.meta.z_kurt = z_map_stats.kurtosis
        # ============================================================================>
        # STORE THE Z MAP
        # ============================================================================>
        z_map_holder.meta.type = 'z-map'
        z_map_holder.parent = dataset_map

        # ============================================================================>
        #####
        # WRITE ALL MAP DISTRIBUTIONS (THESE DON'T USE MUCH SPACE)
        #####
        # ============================================================================>
        # Sampled Map
        write_map_value_distribution(map_vals=dataset_map.data,
                            output_file=dataset.file_manager.get_file('s_map_png'))
        # Mean-Difference
        write_map_value_distribution(map_vals=mean_diff_map.data,
                            output_file=dataset.file_manager.get_file('d_mean_map_png'))
        # Naive Z-Map
        write_map_value_distribution(map_vals=z_map_naive.data,
                            output_file=dataset.file_manager.get_file('z_map_naive_png'),
                            plot_normal=True)
        # Normalised Naive Z-Map
        write_map_value_distribution(map_vals=z_map_naive_normalised.data,
                            output_file=dataset.file_manager.get_file('z_map_naive_normalised_png'),
                            plot_normal=True)
        # Uncertainty Z-Map
        write_map_value_distribution(map_vals=z_map_uncty.data,
                            output_file=dataset.file_manager.get_file('z_map_uncertainty_png'),
                            plot_normal=True)
        # Normalised Uncertainty Z-Map
        write_map_value_distribution(map_vals=z_map_uncty_normalised.data,
                            output_file=dataset.file_manager.get_file('z_map_uncertainty_normalised_png'),
                            plot_normal=True)
        # Corrected Z-Map
        write_map_value_distribution(map_vals=z_map_compl.data,
                            output_file=dataset.file_manager.get_file('z_map_corrected_png'),
                            plot_normal=True)
        # Normalised Corrected Z-Map
        write_map_value_distribution(map_vals=z_map_compl_normalised.data,
                            output_file=dataset.file_manager.get_file('z_map_corrected_normalised_png'),
                            plot_normal=True)
        # Plot Q-Q Plot of Corrected Z-Map to see how normal it is
        write_qq_plot_against_normal(map_vals=z_map_compl_normalised.data,
                            output_file=dataset.file_manager.get_file('z_map_qq_plot_png'))

        # ============================================================================>
        #####
        # LOOK FOR CLUSTERS OF LARGE Z-SCORES
        #####
        # ============================================================================>
        # Contour the grid at a particular Z-Value
        # ============================================================================>
        num_clusters, z_clusters = blob_finder.cluster_high_z_values(z_map=z_map.as_dense().data, point_mask=dataset_total_mask)
        # ============================================================================>
        # Too many points to cluster -- probably a bad dataset
        # ============================================================================>
        if num_clusters == -1:
            # This dataset is too noisy to analyse - flag!
            #pandda.datasets.all_masks().set_mask_value(mask_name='noisy zmap', entry_id=d_handler.tag, value=True)
            log_strs.append('Z-Map too noisy to analyse')
            return dataset, dataset_map.as_sparse(), log_strs

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
            num_clusters, z_clusters = blob_finder.filter_z_clusters_2(z_clusters=z_clusters, dataset=dataset)
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
            num_clusters, z_clusters = blob_finder.filter_z_clusters_3(z_clusters=z_clusters, dataset=dataset)
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
            native_z_map_data = create_native_map(
                                    native_crystal_symmetry=dataset.model.unit_cell,
                                    native_sites=dataset.model.hierarchy.atoms().extract_xyz(),
                                    native_hierarchy=dataset.model.hierarchy,
                                    reference_map=z_map.as_map(), step=0.75,
                                    filename=dataset.file_manager.get_file('native_z_map'))
            # ============================================================================>
            # WRITE REFERENCE FRAME MAPS (NOT ROTATED)
            # ============================================================================>
            dataset_map.to_file(filename=dataset.file_manager.get_file('sampled_map'))
            mean_diff_map.to_file(filename=dataset.file_manager.get_file('mean_diff_map'))
            z_map.to_file(filename=dataset.file_manager.get_file('z_map'))
            # ============================================================================>
            # MASKS/BLOBS/GENERIC MAPS
            # ============================================================================>
            if args.output.developer.write_grid_masks:
                # Write map of grid + symmetry mask
                map_mask = numpy.zeros(grid.grid_size_1d(), dtype=int)
                map_mask.put(map(grid.indexer(), dataset_outer_mask), 1)
                map_mask = flex.double(map_mask.tolist()); map_mask.reshape(flex.grid(grid.grid_size()))
                write_array_to_map(output_file=dataset.file_manager.get_file('grid_mask'), map_data=map_mask, grid=grid)
                # Write map of where the blobs are (high-Z mask)
                highz_points = []; [highz_points.extend(list(x[0])) for x in z_clusters]
                highz_points = [map(int, v) for v in highz_points]
                highz_map_array = numpy.zeros(grid.grid_size_1d(), dtype=int)
                highz_map_array.put(map(grid.indexer(), list(highz_points)), 1)
                map_mask = flex.double(highz_map_array.tolist()); map_mask.reshape(flex.grid(grid.grid_size()))
                write_array_to_map(output_file = d_handler.file_manager.get_file('high_z_mask'), map_data=map_mask, grid=grid)
            # ============================================================================>
            # Write different Z-Maps? (Probably only needed for testing)
            # ============================================================================>
            if args.output.developer.write_all_z_map_types:
                z_map_naive.to_file(filename=dataset.file_manager.get_file('z_map_naive'))
                z_map_naive_normalised.to_file(filename=dataset.file_manager.get_file('z_map_naive_normalised'))
                z_map_uncty.to_file(filename=dataset.file_manager.get_file('z_map_uncertainty'))
                z_map_uncty_normalised.to_file(filename=dataset.file_manager.get_file('z_map_uncertainty_normalised'))
                z_map_compl.to_file(filename=dataset.file_manager.get_file('z_map_corrected'))
                z_map_compl_normalised.to_file(filename=dataset.file_manager.get_file('z_map_corrected_normalised'))

        # ============================================================================>
        # Skip to next dataset if no clusters found
        # ============================================================================>
        if num_clusters > 0:
            log_strs.append('===> {!s} Cluster(s) found.'.format(num_clusters))
        else:
            log_strs.append('===> No Clusters found.')
            return (d_handler, dataset_map.as_sparse(), log_strs)
        assert num_clusters > 0, 'NUMBER OF CLUSTERS AFTER FILTERING == 0!'

        # ============================================================================>
        # Process the identified features
        # ============================================================================>
        for event_idx, (event_points, event_values) in enumerate(z_clusters):
            # Number events from 1
            event_num = event_idx + 1
            # Create a unique identifier for this event
            event_key = (dataset.tag, event_num)
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
            event_mask = GridMask(parent=grid, sites_cart=grid.grid2cart(point_cluster.points, origin=False), max_dist=1.0, min_dist=0.0)
            log_strs.append('=> Event sites ({!s} points) expanded to {!s} points'.format(len(point_cluster.points), len(event_mask.outer_mask_indices())))
            # Select masks to define regions for bdc calculation
            exp_event_idxs = event_mask.outer_mask_indices()
            reference_idxs = grid.global_mask().inner_mask_indices()
            # Calculate the correlation with the mean map as different amounts of the mean map are subtracted
            event_remain_est = estimate_event_background_correction(
                ref_map_data   = map_analyser.statistical_maps.mean_map.as_dense().data,
                query_map_data = dataset_map.as_dense().data,
                feature_idxs   = exp_event_idxs,
                reference_idxs = reference_idxs,
                min_remain     = 1.0-params.background_correction.max_bdc,
                max_remain     = 1.0-params.background_correction.min_bdc,
                bdc_increment  = params.background_correction.increment,
                method         = 'value',
                verbose        = verbose)
            log_strs.append('=> Event Background Correction estimated as {!s}'.format(1-event_remain_est))
            # ============================================================================>
            # Calculate background-corrected map
            # ============================================================================>
            event_map_data = calculate_bdc_subtracted_map(
                                    ref_map_data   = map_analyser.statistical_maps.mean_map.as_dense().as_data,
                                    query_map_data = dataset_map.as_dense().data,
                                    bdc            = 1.0 - event_remain_est)
            event_map = dataset_map.new_from_template(event_map_data, sparse=False)
            # Write out these array (reference and native frames)
            event_map.to_file(filename=dataset.file_manager.get_file('event_map').format(event_num, event_remain_est))
            native_event_map_data = create_native_map(
                                        native_crystal_symmetry=dataset.model.crystal_symmetry,
                                        native_sites=dataset.model.hierarchy.atoms().extract_xyz(),
                                        native_hierarchy=dataset.model.hierarchy,
                                        reference_map=event_map.as_map(), step=0.75,
                                        filename=dataset.file_manager.get_file('native_event_map').format(event_num, event_remain_est))

            # ============================================================================>
            # Find the nearest atom to the event
            # ============================================================================>
            atm = find_nearest_atoms(atoms=protein(dataset.model.hierarchy).atoms(),
                                points=dataset.model.alignment.ref2nat(grid.grid2cart(sites_grid=[map(int,point_cluster.centroid)], origin=False)))[0]
            log_strs.append('=> Nearest Residue to event: Chain {}, Residue {} {}'.format(atm.parent().parent().parent().id,
                                                    atm.parent().parent().unique_resnames()[0], atm.parent().parent().resid()))
            # ============================================================================>
            # Create an event object
            # ============================================================================>
            event_obj = Event(id=point_cluster.id, cluster=point_cluster)
            event_obj.info.estimated_pseudo_occupancy = event_remain_est
            event_obj.info.estimated_bdc              = 1.0 - event_remain_est
            # ============================================================================>
            # Append to dataset handler
            # ============================================================================>
            dataset.events.append(event_obj)

        # ============================================================================>
        # Write out the graph of the background correction correlations
        # ============================================================================>
        #if pandda.settings.plot_graphs:
        #    pandda.something()

        return (dataset, dataset_map.as_sparse(), log_strs)


class DatasetAligner(object):


    def __init__(self, model, other, method, id):
        """
        Shortcut object for aligning datasets.
        Constructed so that the class can be initialised and then called within a multiprocessing function.

        Usage:
            new    = DatasetAligner(...)
            output = new.run()
        or:
            output = DatasetAligner.process(...)
        """

        self.data = (model, other, method, id)

    @classmethod
    def process(cls, model, other, method, id):
        """Process the dataset immediately and return output"""
        return cls(model=model, other=other, method=method, id=id).run()

    def run(self):
        model, other, method, id = self.data
        alignment = model.align_to(other_hierarchy=other.hierarchy, method=method)
        assert alignment.id is None
        alignment.id = id
        return alignment


class MapLoader(object):


    def __init__(self, dataset, grid, reference_map, args, verbose):
        """
        The main object for loading the maps for PanDDA.
        Constructed so that the class can be initialised and then called within a multiprocessing function.

        Usage:
            new    = MapLoader(...)
            output = new.run()
        or:
            output = MapLoader.process(...)
        """

        self.data = (dataset, grid, reference_map, args, verbose)

    @classmethod
    def process(cls, dataset, grid, reference_map, args, verbose):
        """Process the dataset immediately and return output"""
        return cls(dataset=dataset, grid=grid, reference_map=reference_map, args=args, verbose=verbose)

    def run(self):

        dataset, grid, reference_map, args, verbose = self.data

        assert reference_map.is_sparse(), 'Reference map is not in sparse form'

        print('\rLoading Maps for Dataset {!s}'.format(dataset.tag), end=''); sys.stdout.flush()

        # Extract the map data
        fft_map = dataset.data.get_fft_map(structure_factors='truncated')

        # Scale the map
        if   args.params.maps.scaling == 'none':   pass
        elif args.params.maps.scaling == 'sigma':  fft_map.apply_sigma_scaling()
        elif args.params.maps.scaling == 'volume': fft_map.apply_volume_scaling()

        # ============================================================================>
        # MORPH MAPS TO REFERENCE FRAME
        # ============================================================================>

        # Create map handler and map to the reference frame
        native_map_true = dataset.data.get_electron_density_map(fft_map_name='truncated')

        # Extract the map sites from the grid partition
        point_mappings_grid = grid.partition.nn_groups[grid.partition.nn_groups!=-1]
        sites_cart_map = grid.grid2cart(grid.global_mask().outer_mask(), origin=False)
        # Translate the grid partition mappings to the dataset alignment mappings
        mappings_grid2dataset, mappings_dataset2grid = get_closest_mapping_between_lists(list_1=list(grid.partition.sites_cart), list_2=list(dataset.model.alignment.reference_sites))
        point_mappings_dataset = numpy.array([mappings_grid2dataset[i] for i in point_mappings_grid])
        assert point_mappings_dataset.count(-1) == 0
        sites_cart_map_d = dataset.model.alignment.nat2ref(coordinates=sites_cart_map, mappings=point_mappings_dataset)
        map_data = native_map_true.as_map().get_cart_values(sites_cart_map_d)

        # Calculate mean and rms of map values
        map_mean = map_data.min_max_mean().mean
        map_rms  = map_data.standard_deviation_of_the_sample()

        # Scale map to reference
        morphed_map.data = scale_map_to_reference(ref_vals=reference_map.map_data, vals=morphed_map.data)
        # Create map holder
        morphed_map = reference_map.new_from_template(map_data=morphed_map, sparse=True)
        morphed_map.meta.num = dataset.num
        morphed_map.meta.tag = dataset.tag
        morphed_map.meta.type = 'observed map'
        morphed_map.meta.resolution = reference_map.meta.resolution
        morphed_map.meta.map_uncertainty = None
        morphed_map.meta.obs_map_mean = d_map_mean
        morphed_map.meta.obs_map_rms = d_map_rms

        return morphed_map.as_sparse()


class UncertaintyCalculator(object):


    def __init__(self, query_values, ref_values, q_cut=1.5, file_manager=None):
        """
        Object for calculating the uncertinty of query_map relative to ref_map.
        Constructed so that the class can be initialised and then called within a multiprocessing function.

        Usage:
            new    = UncertaintyCalculator(...)
            output = new.run()
        or:
            output = UncertaintyCalculator.process(...)
        """

        self.data = (query_values, ref_values, q_cut, file_manager)

    @classmethod
    def process(cls, query_values, ref_values, q_cut, file_manager):
        return cls(query_values=query_values, ref_values=ref_values, q_cut=q_cut, file_manager=file_manager).run()

    def run(self):
        """Calculate statistics of array of observed values with uncertainties"""

        query_values, ref_values, q_cut, file_manager = self.data

        # Extract the theoretical quantiles that we would expect if these values were from a normal distribution
        the_diff_vals = normal_distribution().quantiles(len(query_values))

        # Select the points in the middle of the distribution
        mid_idxs = (the_diff_vals < q_cut).iselection().intersection((the_diff_vals > -1*q_cut).iselection())
        mid_the_diff_vals = the_diff_vals.select(mid_idxs)

        # Calculate the difference from the reference values
        act_diff_vals = query_values - ref_values
        srt_act_diff_vals = flex.double(sorted(act_diff_vals))
        mid_act_diff_vals = srt_act_diff_vals.select(mid_idxs)

        # Calculate the slope of the centre of the graph
        map_unc, map_off = numpy.polyfit(x=mid_the_diff_vals, y=mid_act_diff_vals, deg=1)

        try:
            import matplotlib
            matplotlib.use('Agg')
            matplotlib.interactive(0)
            from matplotlib import pyplot
            pyplot.style.use('ggplot')
            output_graphs = True
        except:
            output_graphs = False

        if output_graphs and file_manager:
            # Sort query and ref values for plotting
            srt_query_vals = sorted(query_values)
            srt_ref_vals = sorted(ref_values)

            analyse_graphs.mean_obs_scatter(f_name=file_manager.get_file('obs_qqplot_unsorted_png'),
                                            mean_vals=ref_vals, obs_vals=query_values)

            analyse_graphs.sorted_mean_obs_scatter(f_name=file_manager.get_file('obs_qqplot_sorted_png'),
                                            mean_vals=srt_ref_vals, obs_vals=srt_query_vals)

            analyse_graphs.diff_mean_qqplot(f_name=file_manager.get_file('unc_qqplot_png'),
                                            map_off=map_off, map_unc=map_unc,
                                            q_cut=q_cut, obs_diff=srt_act_diff_vals,
                                            quantile=the_diff_vals)

        return map_unc


class DensityStatistics(object):


    def __init__(self, observations_array, uncertainties):
        """
        The main object for loading the maps for PanDDA.
        Constructed so that the class can be initialised and then called within a multiprocessing function.

        Usage:
            new    = MapLoader(...)
            output = new.run()
        or:
            output = MapLoader.process(...)
        """

        assert observations_array.shape[1] == len(uncertainties), 'Arrays are not compatible sizes'

        self.data = (observations_array, uncertainties)

    @classmethod
    def process(cls, observations_array, uncertainties):
        return cls(observations_array=observations_array, uncertainties=uncertainties).run()

    def run(self):
        """Calculate statistics of array of observed values with uncertainties"""
        observations_array, uncertainties = self.data
        out = numpy.empty((observations_array.shape[0], 5))
        for i, observations in enumerate(observations_array):
            assert len(uncertainties) == len(observations)
            out[i,:] = self._calculate_statistics(observations=observations, uncertainties=uncertainties)
        return out

    def _calculate_statistics(self, observations, uncertainties):
        """Calculate statistics for one set of observations and uncertainties"""
        guess_factor = 0.001
        stats_obj = basic_statistics(flex.double(observations))
        stdv = stats_obj.bias_corrected_standard_deviation
        sadj = estimate_true_underlying_sd(obs_vals=observations, obs_error=uncertainties, est_sigma=stdv*guess_factor)
        skew = stats_obj.skew
        kurt = stats_obj.kurtosis
        bimo = (skew**2 + 1)/kurt
        return (stdv, sadj, skew, kurt, bimo)


# ============================================================================>

def get_closest_mapping_between_lists(list_1, list_2):
    """
    Take one set of list_1 and find the identical sites in list_2 (and vice versa).
    Missing sites will be mapped to the closest neighbouring site.
    Return list of indices mapping the site in one to the closest site in the other.
    """

    tmp_idxs_1_to_2 = [list_2.index(s1) if s1 in list_2 else -1 for s1 in list_1]
    out_idxs_1_to_2 = copy.copy(tmp_idxs_1_to_2)
    tmp_idxs_2_to_1 = [list_1.index(s2) if s2 in list_1 else -1 for s2 in list_2]
    out_idxs_2_to_1 = copy.copy(tmp_idxs_2_to_1)

    for in_idxs, out_idxs in ((tmp_idxs_1_to_2,out_idxs_1_to_2),
                              (tmp_idxs_2_to_1,out_idxs_2_to_1)):
        assert in_idxs.count(-1) != len(in_idxs), 'no matching sites found between mappings'

        l = len(in_idxs)
        for i, idx in enumerate(in_idxs):
            d_i = 0
            while out_idxs[i] == -1:
                d_i += 1
                p_i = i+d_i
                n_i = i-d_i
                if   (p_i<l)  and (in_idxs[p_i] != -1):
                    out_idxs[i] = out_idxs[p_i]
                elif (n_i>=0) and (in_idxs[n_i] != -1):
                    out_idxs[i] = out_idxs[n_i]

    return (out_idxs_1_to_2, out_idxs_2_to_1)





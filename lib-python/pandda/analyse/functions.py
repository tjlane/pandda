from __future__ import print_function

import os, sys, copy, time, glob
import traceback

import numpy

from scitbx.array_family import flex
from scitbx.math import basic_statistics
from scitbx.math.distributions import normal_distribution

from libtbx.utils import Sorry, Failure

from bamboo.common import Meta, Info
from bamboo.common.logs import Log
from bamboo.common.path import rel_symlink
from bamboo.pymol_utils import PythonScript
from bamboo.stats.ospina import estimate_true_underlying_sd

from giant.grid.masks import GridMask
from giant.structure.select import protein, non_water, find_nearest_atoms
from giant.xray.maps import ElectronDensityMap
from giant.xray.maps.scale import scale_map_to_reference
from giant.xray.maps.bdc import calculate_varying_bdc_correlations, calculate_maximum_series_discrepancy, calculate_bdc_subtracted_map
from giant.xray.maps.local_align import create_native_map

from pandda import welcome
from pandda.phil import pandda_phil
from pandda.analyse.z_maps import PanddaZMapAnalyser
from pandda.analyse.events import PointCluster, Event
from pandda.analyse import graphs as analyse_graphs

def wrapper_run(c):
    if c is not None:
        try:
            return c.run()
        except:
            return traceback.format_exc()
    else:
        return c

# ============================================================================>


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

        try:
            alignment = model.align_to(other_hierarchy=other.hierarchy, method=method, require_hierarchies_identical=False)
            assert alignment.id is None
            alignment.id = id
            return alignment
        except:
            return traceback.format_exc()

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

        log_file = dataset.file_manager.get_file('dataset_log')

        assert reference_map.is_sparse(), 'Reference map is not in sparse form'

        # ============================================================================>
        # Create map handler in the native coordinate frame
        # ============================================================================>
        # Extract the map data
        fft_map = dataset.data.fft_maps['truncated']
        # Scale the map
        if   args.params.maps.density_scaling == 'none':   pass
        elif args.params.maps.density_scaling == 'sigma':  fft_map.apply_sigma_scaling()
        elif args.params.maps.density_scaling == 'volume': fft_map.apply_volume_scaling()

        # ============================================================================>
        # Morph the map to the reference frame
        # ============================================================================>
        # Extract the map sites from the grid partition
        point_mappings_grid = grid.partition.nn_groups[grid.global_mask().outer_mask_indices()]
        assert sum(point_mappings_grid == -1) == 0
        sites_cart_map = grid.grid2cart(grid.global_mask().outer_mask(), origin_shift=True)
        # Translate the grid partition mappings to the dataset alignment mappings
        mappings_grid2dataset = get_interpolated_mapping_between_coordinates(query_list=grid.partition.sites_cart,
                                                                             ref_list=dataset.model.alignment.reference_sites,
                                                                             tol=0.01)
        point_mappings_dataset = numpy.array([mappings_grid2dataset[i] for i in point_mappings_grid])
        assert sum(point_mappings_dataset == -1) == 0
        sites_cart_map_d = dataset.model.alignment.ref2nat(coordinates=sites_cart_map, mappings=point_mappings_dataset)
        # Create and sample from map object
        native_map_true = ElectronDensityMap.from_fft_map(fft_map)
        morphed_map_data = native_map_true.get_cart_values(sites_cart_map_d)

        # Scale map to reference
        scale_mask = grid.index_on_other(query=grid.global_mask().inner_mask_indices(),
                                         other=grid.global_mask().outer_mask_indices())
        scaled_map_data = scale_map_to_reference(ref_vals   = reference_map.data,
                                                 vals       = morphed_map_data,
                                                 mask_idxs  = flex.size_t(scale_mask))
        # Create map holder
        morphed_map = reference_map.new_from_template(map_data=scaled_map_data, sparse=reference_map.is_sparse())
        morphed_map.meta.num = dataset.num
        morphed_map.meta.tag = dataset.tag
        morphed_map.meta.type = 'observed map'
        morphed_map.meta.resolution = reference_map.meta.resolution
        morphed_map.meta.map_uncertainty = None
        morphed_map.meta.obs_map_mean = morphed_map_data.min_max_mean().mean
        morphed_map.meta.obs_map_rms  = morphed_map_data.standard_deviation_of_the_sample()
        morphed_map.meta.scl_map_mean = scaled_map_data.min_max_mean().mean
        morphed_map.meta.scl_map_rms  = scaled_map_data.standard_deviation_of_the_sample()

        #with open('{}-map-load.log'.format(dataset.tag), 'w') as fh:
        #    fh.write('Dx:'+str(numpy.round(numpy.array(dataset.model.alignment.reference_sites),2))+'\n')
        #    fh.write('Gx:'+str(numpy.round(numpy.array(grid.partition.sites_cart),2))+'\n')
        #    fh.write('G2D'+str(mappings_grid2dataset)+'\n')
        #    fh.write('D2G'+str(mappings_dataset2grid)+'\n')

        # Print a running row of dots
        print('>', end=''); sys.stdout.flush()

        return morphed_map.make_sparse()


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

        assert len(query_values) == len(ref_values)

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
            matplotlib.interactive(False)
            from matplotlib import pyplot
            pyplot.style.use('ggplot')
            output_graphs = True
        except:
            output_graphs = False

        if output_graphs and file_manager:
            # Sort query and ref values for plotting
            srt_query_vals = sorted(query_values)
            srt_ref_vals = sorted(ref_values)

            analyse_graphs.mean_obs_scatter(f_name    = file_manager.get_file('obs_qqplot_unsorted_png'),
                                            mean_vals = ref_values,
                                            obs_vals  = query_values)

            analyse_graphs.sorted_mean_obs_scatter(f_name    = file_manager.get_file('obs_qqplot_sorted_png'),
                                                   mean_vals = srt_ref_vals,
                                                   obs_vals  = srt_query_vals)

            analyse_graphs.uncertainty_qqplot(f_name   = file_manager.get_file('unc_qqplot_png'),
                                              map_off  = map_off,
                                              map_unc  = map_unc,
                                              q_cut    = q_cut,
                                              obs_diff = srt_act_diff_vals,
                                              quantile = the_diff_vals)

        # Print a running row of dots
        print('>', end=''); sys.stdout.flush()

        return map_unc


class DensityStatistics(object):


    def __init__(self, observations_array, uncertainties):
        """
        The main object for loading the maps for PanDDA.
        Constructed so that the class can be initialised and then called within a multiprocessing function.

        Usage:
            new    = DensityStatistics(...)
            output = new.run()
        or:
            output = DensityStatistics.process(...)
        """
        assert observations_array.shape[1] == len(uncertainties), 'Arrays are not compatible sizes'
        assert len(observations_array[0]) == len(uncertainties), 'Sorry. Something has gone wrong...'
        assert len(observations_array[0,:]) == len(uncertainties), 'Sorry. Something has gone wrong...'

        self.data = (observations_array, uncertainties)

    @classmethod
    def process(cls, observations_array, uncertainties):
        return cls(observations_array=observations_array, uncertainties=uncertainties).run()

    def run(self):
        """Calculate statistics of array of observed values with uncertainties"""
        observations_array, uncertainties = self.data
        out = numpy.zeros((observations_array.shape[0], 5), dtype=numpy.float64)
        for i, observations in enumerate(observations_array):
            assert len(uncertainties) == len(observations)
            assert (out[i,:]==0.0).all()
            out[i,:] = self._calculate_statistics(observations=observations, uncertainties=uncertainties)
        return out

    def _calculate_statistics(self, observations, uncertainties):
        """Calculate statistics for one set of observations and uncertainties"""
        guess_factor = 0.001
        stats_obj = basic_statistics(flex.double(numpy.ascontiguousarray(observations)))
        stdv = stats_obj.bias_corrected_standard_deviation
        sadj = estimate_true_underlying_sd(obs_vals=observations, obs_error=uncertainties, est_sigma=stdv*guess_factor)
        skew = stats_obj.skew
        kurt = stats_obj.kurtosis
        bimo = (skew**2 + 1)/kurt
        return (stdv, sadj, skew, kurt, bimo)


class DatasetProcessor(object):


    def __init__(self, dataset, dataset_map, grid, map_analyser, args, verbose):
        """
        The main object for comparing each dataset-map to the ensemble-maps.
        Constructed so that the class can be initialised and then called within a multiprocessing function.

        Usage:
            new_proc = DatasetProcessor(...)
            output   = new_proc.run()
        or:
            output   = DatasetProcessor.process(...)
        """

        self.data = (dataset, dataset_map, grid, map_analyser, args, verbose)

    @classmethod
    def process(cls, dataset, dataset_map, grid, map_analyser, args, verbose):
        """Process the dataset immediately and return output"""
        return cls(dataset=dataset, dataset_map=dataset_map, grid=grid, map_analyser=map_analyser, args=args, verbose=verbose).run()

    def run(self):
        """Process the dataset"""

        dataset, dataset_map, grid, map_analyser, args, verbose = self.data

        # TODO Hardcoded check - to be removed? TODO
        assert dataset_map.is_sparse()

        # ============================================================================>
        # Prepare output objects
        # ============================================================================>
        log_strs = []
        log_file = dataset.file_manager.get_file('dataset_log')
        log = Log(log_file=log_file, verbose=False, silent=True)

        # ============================================================================>
        # Build new blob search object
        # ============================================================================>
        blob_finder = PanddaZMapAnalyser(params = args.params.z_map_analysis,
                                         grid   = grid,
                                         log    = log)
        print('Writing log for dataset {!s} to ...{}'.format(dataset.tag, log_file[log_file.index('processed'):]))

        # ============================================================================>
        # Extract the global mask object from the grid
        # ============================================================================>
        dset_total_temp = grid.global_mask().total_mask_binary().copy()

        # ============================================================================>
        # Generate symmetry masks for this dataset
        # ============================================================================>
        log.bar()
        log('Masking symetry contacts from Z-map.')
        # Generate symmetry contacts for this dataset and align to reference frame
        dataset_sym_copies = dataset.model.crystal_contacts(distance_cutoff=args.params.masks.outer_mask+5, combine_copies=True)
        dataset_sym_copies.atoms().set_xyz(dataset.model.alignment.nat2ref(dataset_sym_copies.atoms().extract_xyz()))
        # Only need to write if writing reference frame maps
        if args.output.developer.write_reference_frame_maps:
            dataset_sym_copies.write_pdb_file(dataset.file_manager.get_file('symmetry_copies'))
        # Extract protein atoms from the symmetry copies
        dataset_sym_sites_cart = non_water(dataset_sym_copies).atoms().extract_xyz()
        # Generate symmetry contacts grid mask
        dataset_mask = GridMask(parent     = grid,
                                sites_cart = dataset_sym_sites_cart,
                                max_dist   = args.params.masks.outer_mask,
                                min_dist   = args.params.masks.inner_mask_symmetry)
        # Combine with the total mask to generate custom mask for this dataset
        dset_total_temp.put(dataset_mask.inner_mask_indices(), 0)
        dset_total_idxs = numpy.where(dset_total_temp)[0]
        log('After masking with symmetry contacts: {} points for Z-map analysis'.format(len(dset_total_idxs)))
        # Write map of grid + symmetry mask
        if args.output.developer.write_reference_frame_grid_masks:
            grid.write_indices_as_map(indices = dset_total_idxs,
                                      f_name  = dataset.file_manager.get_file('grid_mask'),
                                      origin_shift = True)

        # ============================================================================>
        # Generate custom masks for this dataset
        # ============================================================================>
        if args.params.z_map_analysis.masks.selection_string is not None:
            log.bar()
            log('Applying custom mask to the Z-map: "{}"'.format(args.params.z_map_analysis.masks.selection_string))
            cache = dataset.model.hierarchy.atom_selection_cache()
            custom_mask_selection = cache.selection(args.params.z_map_analysis.masks.selection_string)
            custom_mask_sites = dataset.model.hierarchy.select(custom_mask_selection).atoms().extract_xyz()
            log('Masking with {} atoms'.format(len(custom_mask_sites)))
            # Generate custom grid mask
            dataset_mask = GridMask(parent     = grid,
                                    sites_cart = custom_mask_sites,
                                    max_dist   = args.params.z_map_analysis.masks.outer_mask,
                                    min_dist   = args.params.z_map_analysis.masks.inner_mask)
            # Combine with the total mask to generate custom mask for this dataset
            dset_total_temp *= dataset_mask.total_mask_binary()
            dset_total_idxs = numpy.where(dset_total_temp)[0]
            log('After masking with custom mask: {} points for Z-map analysis'.format(len(dset_total_idxs)))
            # Write out mask
            grid.write_indices_as_map(indices = dset_total_idxs,
                                      f_name  = dataset.file_manager.get_file('z_map_mask'),
                                      origin_shift = True)

        # ============================================================================>
        #####
        # CALCULATE Z-MAPS AND LOOK FOR LARGE BLOBS
        #####
        # ============================================================================>
        # Check maps and that all maps are sparse
        # ============================================================================>
        assert dataset_map.data is not None, 'Something has gone wrong - this dataset has no loaded map'
        assert dataset_map.is_sparse() is map_analyser.statistical_maps.mean_map.is_sparse()
        assert dataset_map.is_sparse() is map_analyser.statistical_maps.medn_map.is_sparse()
        assert dataset_map.is_sparse() is map_analyser.statistical_maps.stds_map.is_sparse()
        assert dataset_map.is_sparse() is map_analyser.statistical_maps.sadj_map.is_sparse()
        # ============================================================================>
        # CALCULATE MEAN-DIFF MAPS
        # ============================================================================>
        mean_diff_map = map_analyser.calculate_z_map(map=dataset_map, method='none')
#        # ============================================================================>
#        # NAIVE Z-MAP - NOT USING UNCERTAINTY ESTIMATION OR ADJUSTED STDS
#        # ============================================================================>
#        z_map_naive = map_analyser.calculate_z_map(map=dataset_map, method='naive')
#        z_map_naive_normalised = z_map_naive.normalised_copy()
        # ============================================================================>
        # UNCERTAINTY Z-MAP - NOT USING ADJUSTED STDS
        # ============================================================================>
        z_map_uncty = map_analyser.calculate_z_map(map=dataset_map, uncertainty=dataset_map.meta.map_uncertainty, method='uncertainty')
        z_map_uncty_normalised = z_map_uncty.normalised_copy()
        # ============================================================================>
        # ADJUSTED+UNCERTAINTY Z-MAP
        # ============================================================================>
        z_map_compl = map_analyser.calculate_z_map(map=dataset_map, uncertainty=dataset_map.meta.map_uncertainty, method='adjusted+uncertainty')
        z_map_compl_normalised = z_map_compl.normalised_copy()

        # ============================================================================>
        # SELECT WHICH MAP TO DO THE BLOB SEARCHING ON
        # ============================================================================>
#        if args.params.statistical_maps.z_map_type == 'naive':
#            z_map = z_map_naive_normalised
#            z_map_stats = basic_statistics(flex.double(z_map_naive.data))
        if args.params.statistical_maps.z_map_type == 'uncertainty':
            z_map = z_map_uncty_normalised
            z_map_stats = basic_statistics(flex.double(z_map_uncty.data))
        elif args.params.statistical_maps.z_map_type == 'adjusted+uncertainty':
            z_map = z_map_compl_normalised
            z_map_stats = basic_statistics(flex.double(z_map_compl.data))
        else: raise Exception('Invalid Z-map type')

        # ============================================================================>
        # RECORD Z-MAP FOR STATISTICS
        # ============================================================================>
        # Calculate statistics of z-maps
        dataset_map.meta.z_mean = z_map_stats.mean
        dataset_map.meta.z_stdv = z_map_stats.bias_corrected_standard_deviation
        dataset_map.meta.z_skew = z_map_stats.skew
        dataset_map.meta.z_kurt = z_map_stats.kurtosis
        # ============================================================================>
        z_map.meta.type = 'z-map'
        # ============================================================================>

        # ============================================================================>
        #####
        # WRITE ALL MAP DISTRIBUTIONS (THESE DON'T USE MUCH SPACE)
        #####
        # ============================================================================>
        # Sampled Map
        analyse_graphs.map_value_distribution(f_name    = dataset.file_manager.get_file('s_map_png'),
                                              plot_vals = dataset_map.get_map_data(sparse=True))
        # Mean-Difference
        analyse_graphs.map_value_distribution(f_name    = dataset.file_manager.get_file('d_mean_map_png'),
                                              plot_vals = mean_diff_map.get_map_data(sparse=True))
#        # Naive Z-Map
#        analyse_graphs.map_value_distribution(f_name      = dataset.file_manager.get_file('z_map_naive_png'),
#                                              plot_vals   = z_map_naive.get_map_data(sparse=True),
#                                              plot_normal = True)
#        # Normalised Naive Z-Map
#        analyse_graphs.map_value_distribution(f_name      = dataset.file_manager.get_file('z_map_naive_normalised_png'),
#                                              plot_vals   = z_map_naive_normalised.get_map_data(sparse=True),
#                                              plot_normal = True)
        # Uncertainty Z-Map
        analyse_graphs.map_value_distribution(f_name      = dataset.file_manager.get_file('z_map_uncertainty_png'),
                                              plot_vals   = z_map_uncty.get_map_data(sparse=True),
                                              plot_normal = True)
        # Normalised Uncertainty Z-Map
        analyse_graphs.map_value_distribution(f_name      = dataset.file_manager.get_file('z_map_uncertainty_normalised_png'),
                                              plot_vals   = z_map_uncty_normalised.get_map_data(sparse=True),
                                              plot_normal = True)
        # Corrected Z-Map
        analyse_graphs.map_value_distribution(f_name      = dataset.file_manager.get_file('z_map_corrected_png'),
                                              plot_vals   = z_map_compl.get_map_data(sparse=True),
                                              plot_normal = True)
        # Normalised Corrected Z-Map
        analyse_graphs.map_value_distribution(f_name      = dataset.file_manager.get_file('z_map_corrected_normalised_png'),
                                              plot_vals   = z_map_compl_normalised.get_map_data(sparse=True),
                                              plot_normal = True)
        # Plot Q-Q Plot of Corrected Z-Map to see how normal it is
        analyse_graphs.qq_plot_against_normal(f_name    = dataset.file_manager.get_file('z_map_qq_plot_png'),
                                              plot_vals = z_map_compl_normalised.get_map_data(sparse=True))

        # ============================================================================>
        #####
        # LOOK FOR CLUSTERS OF LARGE Z-SCORES
        #####
        # ============================================================================>
        # Contour the grid at a particular Z-Value
        # ============================================================================>
        num_clusters, z_clusters = blob_finder.cluster_high_z_values( z_map_data     = z_map.get_map_data(sparse=False),
                                                                      point_mask_idx = dset_total_idxs)
        # ============================================================================>
        # Too many points to cluster -- probably a bad dataset
        # ============================================================================>
        if num_clusters == -1:
            # This dataset is too noisy to analyse - flag!
            log_strs.append('Z-Map too noisy to analyse -- not sure what has gone wrong here...')
            return dataset, dataset_map.meta, log_strs

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
        # write dataset maps in the reference frame
        # ============================================================================>
        if args.output.developer.write_reference_frame_maps:
            dataset_map.to_file(filename=dataset.file_manager.get_file('sampled_map'), space_group=grid.space_group())
            mean_diff_map.to_file(filename=dataset.file_manager.get_file('mean_diff_map'), space_group=grid.space_group())
            z_map.to_file(filename=dataset.file_manager.get_file('z_map'), space_group=grid.space_group())
        # ============================================================================>
        # Write out mask of the high z-values
        # ============================================================================>
        if args.output.developer.write_reference_frame_grid_masks:
            # Write map of where the blobs are (high-Z mask)
            highz_points = []; [highz_points.extend(list(x[0])) for x in z_clusters]
            highz_points = [map(int, v) for v in highz_points]
            highz_indices = map(grid.indexer(), list(highz_points))
            grid.write_indices_as_map(indices = highz_indices,
                                      f_name  = dataset.file_manager.get_file('high_z_mask'),
                                      origin_shift=True)
        # ============================================================================>
        # Write different Z-Maps? (Probably only needed for testing)
        # ============================================================================>
        if args.output.developer.write_reference_frame_all_z_map_types:
#            z_map_naive.to_file(filename=dataset.file_manager.get_file('z_map_naive'), space_group=grid.space_group())
#            z_map_naive_normalised.to_file(filename=dataset.file_manager.get_file('z_map_naive_normalised'), space_group=grid.space_group())
            z_map_uncty.to_file(filename=dataset.file_manager.get_file('z_map_uncertainty'), space_group=grid.space_group())
            z_map_uncty_normalised.to_file(filename=dataset.file_manager.get_file('z_map_uncertainty_normalised'), space_group=grid.space_group())
            z_map_compl.to_file(filename=dataset.file_manager.get_file('z_map_corrected'), space_group=grid.space_group())
            z_map_compl_normalised.to_file(filename=dataset.file_manager.get_file('z_map_corrected_normalised'), space_group=grid.space_group())

        # ============================================================================>
        # Skip to next dataset if no clusters found
        # ============================================================================>
        if num_clusters > 0:
            log_strs.append('===> {!s} Cluster(s) found.'.format(num_clusters))
        else:
            log_strs.append('===> No Clusters found.')
            return (dataset, dataset_map.meta, log_strs)
        assert num_clusters > 0, 'NUMBER OF CLUSTERS AFTER FILTERING == 0!'

        # ============================================================================>
        # Extract the map data in non-sparse format
        # ============================================================================>
        dset_map_data = dataset_map.get_map_data(sparse=False)
        avrg_map_data = map_analyser.average_map().get_map_data(sparse=False)
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
            event_mask = GridMask(parent=grid, sites_cart=grid.grid2cart(point_cluster.points, origin_shift=True), max_dist=2.0, min_dist=0.0)
            log_strs.append('=> Event sites ({!s} points) expanded to {!s} points'.format(len(point_cluster.points), len(event_mask.outer_mask_indices())))
            # Select masks to define regions for bdc calculation
            exp_event_idxs = flex.size_t(event_mask.outer_mask_indices())
            reference_idxs = flex.size_t(grid.global_mask().inner_mask_indices())
            # ============================================================================>
            # Generate BDC-estimation curve and estimate BDC
            # ============================================================================>
            event_remains, event_corrs, global_corrs = calculate_varying_bdc_correlations(
                ref_map_data   = avrg_map_data,
                query_map_data = dset_map_data,
                feature_idxs   = exp_event_idxs,
                reference_idxs = reference_idxs,
                min_remain     = 1.0-args.params.background_correction.max_bdc,
                max_remain     = 1.0-args.params.background_correction.min_bdc,
                bdc_increment  = args.params.background_correction.increment,
                verbose        = verbose)
            event_remain_est = calculate_maximum_series_discrepancy(
                labels   = event_remains,
                series_1 = global_corrs,
                series_2 = event_corrs)
            analyse_graphs.write_occupancy_graph(
                f_name=dataset.file_manager.get_file('bdc_est_png').format(event_num),
                x_values=event_remains,
                global_values=global_corrs,
                local_values=event_corrs)
            log_strs.append('=> Event Background Correction estimated as {!s}'.format(1-event_remain_est))
            # Reporting (log is normally silenced)
            blob_finder.log('Min-Max: {} {}'.format(1.0-args.params.background_correction.max_bdc, 1.0-args.params.background_correction.min_bdc))
            blob_finder.log('Event number: {}'.format(event_num))
            blob_finder.log('Event Remains: {}'.format(','.join(map(str,event_remains))))
            blob_finder.log('Event Corrs:  {}'.format(','.join(map(str,event_corrs))))
            blob_finder.log('Global Corrs: {}'.format(','.join(map(str,global_corrs))))
            # Apply multiplier if provided
            blob_finder.log('Applying multiplier to output 1-BDC: {}'.format(args.params.background_correction.output_multiplier))
            event_remain_est = min(event_remain_est*args.params.background_correction.output_multiplier, 1.0-args.params.background_correction.min_bdc)
            # ============================================================================>
            # Calculate the map correlations at the selected BDC
            # ============================================================================>
            event_map_data = calculate_bdc_subtracted_map(
                                    ref_map_data   = avrg_map_data,
                                    query_map_data = dset_map_data,
                                    bdc            = 1.0 - event_remain_est)
            global_corr = numpy.corrcoef(event_map_data.select(reference_idxs), avrg_map_data.select(reference_idxs))[0,1]
            local_corr  = numpy.corrcoef(event_map_data.select(exp_event_idxs), avrg_map_data.select(exp_event_idxs))[0,1]
            # ============================================================================>
            # Write out EVENT map (in the reference frame) and grid masks
            # ============================================================================>
            if args.output.developer.write_reference_frame_maps:
                event_map = dataset_map.new_from_template(event_map_data, sparse=False)
                event_map.to_file(filename=dataset.file_manager.get_file('event_map').format(event_num, event_remain_est), space_group=grid.space_group())
            if args.output.developer.write_reference_frame_grid_masks:
                grid.write_indices_as_map(indices=event_mask.outer_mask_indices(), f_name=dataset.file_manager.get_file('grid_mask').replace('.ccp4','')+'-event-mask-{}.ccp4'.format(event_num))

            # ============================================================================>
            # Find the nearest atom to the event
            # ============================================================================>
            atm = find_nearest_atoms(atoms=list(protein(dataset.model.hierarchy).atoms_with_labels()),
                                     query=dataset.model.alignment.ref2nat(grid.grid2cart(sites_grid=[map(int,point_cluster.centroid)], origin_shift=True)))[0]
            log_strs.append('=> Nearest Residue to event: Chain {}, Residue {} {}'.format(atm.chain_id, atm.resname, atm.resid()))
            # ============================================================================>
            # Create an event object
            # ============================================================================>
            event_obj = Event(id=point_cluster.id, cluster=point_cluster)
            event_obj.info.estimated_pseudo_occupancy = event_remain_est
            event_obj.info.estimated_bdc              = 1.0 - event_remain_est
            event_obj.info.global_correlation = global_corr
            event_obj.info.local_correlation  = local_corr
            # ============================================================================>
            # Append to dataset handler
            # ============================================================================>
            dataset.events.append(event_obj)

        # ============================================================================>
        # Write out pymol script to load all of the maps easily
        # ============================================================================>
        pml = PythonScript()
        pml.set_normalise_maps(False)
        # Load Structures
        name = pml.load_pdb(f_name=dataset.file_manager.get_file('aligned_model'))
        pml.repr_as(obj=name, style='sticks')
        name = pml.load_pdb(f_name=dataset.file_manager.get_file('symmetry_copies'))
        pml.repr_hide(obj=name)
        # Load Sampled Map
        name = pml.load_map(f_name=dataset.file_manager.get_file('sampled_map'))
        mesh = pml.make_mesh(obj=name, contour_level=1.0, colour='blue')
        # Load Z-maps
        name = pml.load_map(f_name=dataset.file_manager.get_file('z_map'))
        mesh = pml.make_mesh(obj=name, mesh_suffix='.plus', contour_level=3.0, colour='green')
        mesh = pml.make_mesh(obj=name, mesh_suffix='.mins', contour_level=-3.0, colour='red')
        # Load Event maps
        for f in sorted(glob.glob(dataset.file_manager.get_file('event_map').format('*','*'))):
            name = pml.load_map(f_name=f)
            mesh = pml.make_mesh(obj=name, contour_level=float(f.split('_')[-2]), colour='hotpink')
        # Load Miscellaneous maps (e.g. masks)
        for f in sorted(glob.glob(os.path.join(dataset.file_manager.get_dir('root'),'*mask*.ccp4'))):
            name = pml.load_map(f_name=f)
            mesh = pml.make_mesh(obj=name, contour_level=0.0, colour='grey')

        pml.write_script(f_name=dataset.file_manager.get_file('pymol_script'), overwrite=True)

        return (dataset, dataset_map.meta, log_strs)


class NativeMapMaker(object):


    def __init__(self, dataset, map_obj, sites_mask, filename, args, verbose):
        """
        The main object for comparing each dataset-map to the ensemble-maps.
        Constructed so that the class can be initialised and then called within a multiprocessing function.

        Usage:
            new_proc = DatasetProcessor(...)
            output   = new_proc.run()
        or:
            output   = DatasetProcessor.process(...)
        """

        self.data = (dataset, map_obj, sites_mask, filename, args, verbose)

    @classmethod
    def process(cls, dataset, map_obj, sites_mask, filename, args, verbose):
        """Process the dataset immediately and return output"""
        return cls(dataset=dataset, map_obj=map_obj, sites_mask=sites_mask, filename=filename, args=args, verbose=verbose).run()

    def run(self):
        """Process the dataset"""

        t1 = time.time()
        dataset, map_obj, sites_mask, filename, args, verbose = self.data

        native_map_data = create_native_map(
                            native_crystal_symmetry = dataset.model.crystal_symmetry,
                            native_sites            = dataset.model.alignment.ref2nat(sites_mask),
                            alignment               = dataset.model.alignment,
                            reference_map           = map_obj.make_dense(),
                            site_mask_radius        = args.params.masks.outer_mask,
                            step                    = args.params.maps.grid_spacing,
                            filename                = filename
                        )
#        print('Time: '+str(time.time()-t1))

        return None

# ============================================================================>

# TODO move to more sensible permanent place

def get_interpolated_mapping_between_coordinates(query_list, ref_list, tol=0.01):
    """
    Take each of query_list and find the sites in ref_list (within tolerance).
    Missing sites will be interpolated to the closest neighbouring site.
    Return list of indices mapping the site in one to the closest site in the other.
    """
    ref_list = flex.vec3_double(ref_list)
    tmp_idxs_q_to_r = [closest_point_within_tolerance(query=q, ref_list=ref_list, tol=tol) for q in query_list]
    assert tmp_idxs_q_to_r.count(-1) != len(tmp_idxs_q_to_r), 'no matching sites found between mappings'
    out_idxs_q_to_r = copy.copy(tmp_idxs_q_to_r)
    l = len(tmp_idxs_q_to_r)
    # Populate the missing values with the nearest value
    for i in range(l):
        d_i = 0
        while out_idxs_q_to_r[i] == -1:
            d_i += 1; p_i = i+d_i; n_i = i-d_i
            if   (p_i<l)  and (tmp_idxs_q_to_r[p_i] != -1):
                out_idxs_q_to_r[i] = out_idxs_q_to_r[p_i]
            elif (n_i>=0) and (tmp_idxs_q_to_r[n_i] != -1):
                out_idxs_q_to_r[i] = out_idxs_q_to_r[n_i]

    return out_idxs_q_to_r

def closest_point_within_tolerance(query, ref_list, tol):
    dist_sq = list((ref_list-query).dot())
    dmin_sq = min(dist_sq)
    if dmin_sq > tol**2:
        return -1
    return dist_sq.index(dmin_sq)


from __future__ import print_function

import os, sys, glob, time, re
import copy, warnings, traceback

import scipy.cluster
import scipy.linalg

import numpy, pandas

from libtbx import easy_mp
from libtbx.utils import Sorry, Failure

import cctbx.maptbx
import scitbx.sparse
import scitbx.matrix

from scitbx.array_family import flex
from scitbx.math import basic_statistics
from scitbx.math.distributions import normal_distribution

from bamboo.common import Meta, Info
from bamboo.common.logs import Log
from bamboo.common.file import FileManager
from bamboo.common.path import easy_directory, rel_symlink, delete_with_glob
from bamboo.common.status import status_bar, status_bar_2
from bamboo.common.command import CommandManager

from bamboo.plot import bar

from giant.manager import Program

from giant.grid import grid_handler
from giant.grid.masks import spherical_mask, atomic_mask, grid_mask

from giant.xray.data.utils import extract_structure_factors
from giant.xray.maps import scale_map_to_reference
from giant.xray.symmetry import find_symmetry_equivalent_groups
from giant.stats.ospina import estimate_true_underlying_sd
from giant.stats.cluster import find_connected_groups, generate_group_idxs

from pandda import analyse_graphs
from pandda.phil import pandda_phil
from pandda.misc import write_map_value_distribution, write_qq_plot_against_normal, write_array_to_map, rotate_map
from pandda.lists import MapList, PanddaStatMapList, PanddaMultipleStatMapList
from pandda.events import cluster_events
from pandda.holders import MapHolder, MapHolderList, DatasetHandlerList
from pandda.handlers import DatasetHandler, ReferenceDatasetHandler, map_handler, align_dataset_to_reference

from pandda import PANDDA_TEXT, PANDDA_VERSION
from pandda.constants import *

def round_no_fail(a, decimals=0):
    try:    return numpy.round(a, decimals)
    except: return None

def map_statistics_map_func(arg_dict):
    map_vals          = arg_dict['map_vals']
    map_uncertainties = arg_dict['map_uncertainties']
    # Starting guess for the underlying sigma is raw_std/guess_factor
    guess_factor = 0.001
    # Calculate statistics from extracted map values
    stats_obj = basic_statistics(flex.double(map_vals.tolist()))
    # Calculate Standard Deviation of the map values
    stdv  = stats_obj.bias_corrected_standard_deviation
    # Calculate the adjusted standard deviation
    sadj = estimate_true_underlying_sd(obs_vals=map_vals.tolist(), obs_error=map_uncertainties, est_sigma=stdv*guess_factor)
    # Calculate other statistics
    skew = stats_obj.skew
    kurt = stats_obj.kurtosis
    bimo = (skew**2 + 1)/kurt
    return (stdv, sadj, skew, kurt, bimo)

def load_dataset_map_func(arg_dict):
    return DatasetHandler(**arg_dict)

def align_map_func(arg_dict):
    global_mx, local_mxs = align_dataset_to_reference(
                                    d_handler   = arg_dict['d_handler'],
                                    ref_handler = arg_dict['r_handler'],
                                    method      = arg_dict['method']   )
    return (arg_dict['d_handler'].tag, (global_mx, local_mxs))

def load_maps_map_func(arg_dict):
    """Function to allow maps to be loaded in parallel"""

    d_handler       = arg_dict['d_handler']
    grid            = arg_dict['grid']
    params          = arg_dict['params']
    map_resolution  = arg_dict['map_resolution']
    ref_map_holder  = arg_dict['ref_map_holder']
    masked_idxs     = arg_dict['masked_idxs']
    masked_cart_ref = arg_dict['masked_cart_ref']
    masked_cart_map = arg_dict['masked_cart_map']
    scaling_mask    = arg_dict['scaling_mask']

    print('\rLoading Maps for Dataset {!s}'.format(d_handler.tag), end=''); sys.stdout.flush()

    fft_map = d_handler.fft_map

    # Scale the map
    if   params.maps.scaling == 'none':   pass
    elif params.maps.scaling == 'sigma':  fft_map.apply_sigma_scaling()
    elif params.maps.scaling == 'volume': fft_map.apply_volume_scaling()

    #################################
    # MORPH MAPS TO REFERENCE FRAME #
    #################################
    # Create map handler and map to the reference frame
    native_map_handler = map_handler(   map_data  = fft_map.real_map(),
                                        unit_cell = fft_map.unit_cell() )
    # Transform Coordinates
    masked_cart_d = d_handler.transform_from_reference( points = masked_cart_ref,
                                                        method = params.alignment.method,
                                                        point_mappings = masked_cart_map )
    # Sample the map at these points
    masked_vals_d = native_map_handler.get_cart_values(masked_cart_d)

    # Calculate mean and rms of map values
    d_map_mean = masked_vals_d.min_max_mean().mean
    d_map_rms  = masked_vals_d.standard_deviation_of_the_sample()

    # Calculate the sparse vector of the masked map values
    morphed_map_sparse = scitbx.sparse.vector(grid.size_1d(), dict(zip(masked_idxs, masked_vals_d)))
    morphed_map        = morphed_map_sparse.as_dense_vector()
    # Reshape into right shape of the grid
    morphed_map.reshape(grid)
    # Scale map to reference
    morphed_map = scale_map_to_reference(ref_map=ref_map_holder.map, map=morphed_map, mask_idxs=scaling_mask)
    # Create map holder
    map_holder = MapHolder( num = d_handler.num,
                            tag = d_handler.tag,
                            map = morphed_map,
                            # Change these for the 'fake' grid unit_cell and a P1 space_group
                            unit_cell   = fft_map.unit_cell(),
                            space_group = fft_map.space_group(),
                            meta        = Meta({'type'            : 'observed map',
                                                'resolution'      : map_resolution,
                                                'map_uncertainty' : None,
                                                'obs_map_mean'    : d_map_mean,
                                                'obs_map_rms'     : d_map_rms }),
                            parent      = None  )

    return map_holder

class PanddaMapAnalyser(object):
    """Class to hold dataset maps, statistical maps and meta data for a set of related maps. Also holds functions for analysing the maps."""
    def __init__(self, dataset_maps, meta, statistical_maps=None, parent=None, log=None):
        # Validate the meta
        assert isinstance(meta, Meta), 'meta must be of type Meta. Type given: {!s}'.format(type(meta))
        assert hasattr(meta, 'resolution')
        assert hasattr(meta, 'grid_size')
        assert hasattr(meta, 'grid_size_1d')
        self.meta = meta
        # Validate the dataset maps
        if dataset_maps:
            assert isinstance(dataset_maps, MapHolderList), 'dataset_maps must be stored in a MapHolderList. Type given: {!s}'.format(type(dataset_maps))
            self.dataset_maps = dataset_maps
        else:
            self.dataset_maps = None
        if statistical_maps:
            assert isinstance(statistical_maps, PanddaStatMapList), 'statistical_maps must be of type MapList. Type given: {!s}'.format(type(statistical_maps))
            self.statistical_maps = statistical_maps
        else:
            self.statistical_maps = PanddaStatMapList()
            self.statistical_maps.meta = self.meta
        # Validate the parent object (main pandda object)
        if parent:
            assert isinstance(parent, PanddaMultiDatasetAnalyser), 'parent must be of type PanddaMultiDatasetAnalyser. Type given: {!s}'.format(type(parent))
            self.parent = parent
        else:
            self.parent = None

        if log:
            self.log = log
        elif (self.parent and self.parent.log):
            self.log = self.parent.log
        else:
            self.log = Log(verbose=False)

    def validate_maps(self):
        """Check that all of the added maps are the same size etc..."""

        for mh in self.dataset_maps.all():
            print('Checking Map {!s}'.format(mh.tag))

    def calculate_mean_map(self, masked_idxs=None, mask_name=None):
        """Calculate the mean map from all of the different observations"""

        if not masked_idxs: masked_idxs = range(0, self.meta.grid_size_1d)
        else:               assert max(masked_idxs) < self.meta.grid_size_1d, 'masked_idxs out of range of map'

        # Create statistics objects for each grid point
        self.log('----------------------------------->>>', True)
        self.log('Calculating Mean Map', True)

        # Chunk the points into groups - Compromise between cpu time and memory usage - ~200 dataset -> chunksize of 5000
        chunk_size = 1000*int(2000/self.dataset_maps.size(mask_name=mask_name))
        chunked_indices = [masked_idxs[i:i + chunk_size] for i in range(0, len(masked_idxs), chunk_size)]
        num_chunks = len(chunked_indices)

        self.log('Iterating through {!s} points in {!s} chunks'.format(len(masked_idxs), num_chunks), True)
        t1 = time.time()

        # Calculate Mean Maps
        mean_map_vals = numpy.zeros(self.meta.grid_size_1d)
        medn_map_vals = numpy.zeros(self.meta.grid_size_1d)
        # Calculate the mean map across the datasets
        for i_chunk, chunk_idxs in enumerate(chunked_indices):
            status_bar_2(n=i_chunk, n_max=num_chunks)
            # Pull out the values for this chunk
            p_map_vals = numpy.array([mh.map.select(chunk_idxs) for mh in self.dataset_maps.mask(mask_name=mask_name)])
            # Check that the output values are the expected dimensions
            if i_chunk+1 < num_chunks:
                assert len(p_map_vals) == self.dataset_maps.size(mask_name=mask_name)
                assert len(p_map_vals.T) == chunk_size
            # Calculate the mean of this chunk
            p_map_means = numpy.mean(p_map_vals, axis=0)
            mean_map_vals.put(chunk_idxs, p_map_means)
            # Calculate the median of this chunk
            p_map_medns = numpy.median(p_map_vals, axis=0)
            medn_map_vals.put(chunk_idxs, p_map_medns)
        # Print full bar
        status_bar_2(n=num_chunks, n_max=num_chunks)

        t2 = time.time()
        self.log('> MEAN MAP CALCULATION > Time Taken: {!s} seconds'.format(int(t2-t1)))

        # Save and reshape the array to be the shape of the map
        self.statistical_maps.mean_map = flex.double(mean_map_vals.tolist())
        self.statistical_maps.mean_map.reshape(flex.grid(self.meta.grid_size))
        # Save and reshape the array to be the shape of the map
        self.statistical_maps.medn_map = flex.double(medn_map_vals.tolist())
        self.statistical_maps.medn_map.reshape(flex.grid(self.meta.grid_size))

        return self.statistical_maps.mean_map

    def calculate_map_uncertainties(self, masked_idxs=None, mask_name=None, output_graphs=True):
        """Calculate the uncertainty in each of the different maps"""

        try:
            import matplotlib
            matplotlib.use('Agg')
            matplotlib.interactive(0)
            from matplotlib import pyplot
            pyplot.style.use('ggplot')
        except:
            output_graphs = False

        if not masked_idxs: masked_idxs = range(0, self.meta.grid_size_1d)
        else:               assert max(masked_idxs) < self.meta.grid_size_1d, 'masked_idxs out of range of map'

        # Section of q-q plot to calculate slope over (-q_cut -> +q_cut)
        q_cut = 1.5
        # Extract the mean map values
        masked_mean_vals = self.statistical_maps.mean_map.select(masked_idxs)
        # ...and sort them
        sorted_mean_vals = sorted(masked_mean_vals)
        # Extract the theoretical quantiles that we would expect if these values were from a normal distribution
        expected_diff_vals = normal_distribution().quantiles(len(masked_idxs))
        # Select the points in the middle of the distribution
        middle_indices = (expected_diff_vals < q_cut).iselection().intersection((expected_diff_vals > -1*q_cut).iselection())
        middle_expected_diff_vals = expected_diff_vals.select(middle_indices)

        self.log('----------------------------------->>>')
        self.log('Calculating Map Uncertainties')
        t1 = time.time()

        for i_mh, mh in enumerate(self.dataset_maps.mask(mask_name=mask_name)):

            if mh.meta.map_uncertainty is not None:
                print('SKIPPING Dataset {!s} ({!s}/{!s})                                 '.format(mh.tag, i_mh+1, self.dataset_maps.size(mask_name=mask_name)))
                continue
            else:
                status_bar_2(n=i_mh+1, n_max=self.dataset_maps.size(mask_name=mask_name))

            # Extract the map values as the masked indices
            masked_map_vals = mh.map.select(masked_idxs)
            # Subtract the mean map from the observed map
            diff_mean_map_vals = masked_map_vals - masked_mean_vals
            # Sort the map values for both the masked map and the deviations from the mean map
            sorted_observed_vals = sorted(masked_map_vals)
            sorted_observed_diff_vals = sorted(diff_mean_map_vals)
            # Select the middle indices
            middle_observed_diff_vals = flex.double(sorted_observed_diff_vals).select(middle_indices)
            # Calculate the slope of the centre of the graph
            fit_coeffs = numpy.polyfit(x=middle_expected_diff_vals, y=middle_observed_diff_vals, deg=1)
            map_unc, map_off = fit_coeffs[0:2]

            # Save the uncertainty
            mh.meta.map_uncertainty = map_unc

            if output_graphs and mh.parent:

                analyse_graphs.mean_obs_scatter( f_name    = mh.parent.output_handler.get_file('obs_qqplot_unsorted_png'),
                                                 mean_vals = masked_mean_vals,
                                                 obs_vals  = masked_map_vals     )

                analyse_graphs.sorted_mean_obs_scatter( f_name    = mh.parent.output_handler.get_file('obs_qqplot_sorted_png'),
                                                        mean_vals = sorted_mean_vals,
                                                        obs_vals  = sorted_observed_vals    )

                analyse_graphs.diff_mean_qqplot( f_name   = mh.parent.output_handler.get_file('unc_qqplot_png'),
                                                 map_off  = map_off,
                                                 map_unc  = map_unc,
                                                 q_cut    = q_cut,
                                                 obs_diff = sorted_observed_diff_vals,
                                                 quantile = expected_diff_vals  )

        t2 = time.time()
        self.log('> MAP UNCERTAINTY CALCULATION > Time Taken: {!s} seconds'.format(int(t2-t1)))

        return [mh.meta.map_uncertainty for mh in self.dataset_maps.mask(mask_name=mask_name)]

    def calculate_statistical_maps(self, masked_idxs=None, mask_name=None, cpus=1, ignore_warnings=True):
        """Take the sampled maps and calculate statistics for each grid point across the datasets"""

        # Create statistics objects for each grid point
        self.log('----------------------------------->>>', True)
        self.log('Calculating Statistics of Grid Points (using {!s} cores)'.format(cpus), True)

        if ignore_warnings:
            self.log('Suppressing Numpy Iteration Warnings... (Nothing to worry about)')
            warnings.simplefilter('ignore', category=RuntimeWarning)

        if not masked_idxs: masked_idxs = range(0, self.meta.grid_size_1d)
        else:               assert max(masked_idxs) < self.meta.grid_size_1d, 'masked_idxs out of range of map'

        # Extract the map uncertainties
        all_uncertainties = [mh.meta.map_uncertainty for mh in self.dataset_maps.mask(mask_name=mask_name)]

        # Chunk the points into groups - Compromise between cpu time and memory usage - ~200 dataset -> chunksize of 5000
        chunk_size = 1000*int(2000/self.dataset_maps.size(mask_name=mask_name))
        chunked_indices = [masked_idxs[i:i + chunk_size] for i in range(0, len(masked_idxs), chunk_size)]
        num_chunks = len(chunked_indices)

        self.log('Iterating through {!s} points in {!s} chunks'.format(len(masked_idxs), num_chunks), True)
        t1 = time.time()

        # Statistics objects for each of the points we're interested in
        masked_point_statistics = []

        # Calculate the statistics of the map values across the datasets
        for i_chunk, chunk_idxs in enumerate(chunked_indices):
            status_bar_2(n=i_chunk, n_max=num_chunks)
            # Select map values from each map
            p_map_vals = numpy.array([mh.map.select(chunk_idxs) for mh in self.dataset_maps.mask(mask_name=mask_name)])
            if i_chunk+1 < num_chunks:
                assert len(p_map_vals) == self.dataset_maps.size(mask_name=mask_name)
                assert len(p_map_vals.T) == chunk_size

            # Generate arg_list for analysis
            arg_list = [{'map_vals':map_vals, 'map_uncertainties':all_uncertainties} for map_vals in p_map_vals.T]
            point_statistics = easy_mp.pool_map(func=map_statistics_map_func, args=arg_list, processes=cpus)

            masked_point_statistics.extend(point_statistics)

        status_bar_2(n=num_chunks, n_max=num_chunks)
        assert len(masked_point_statistics) == len(masked_idxs)

        t2 = time.time()
        self.log('> MAP POINT ANALYSIS > Time Taken: {!s} seconds'.format(int(t2-t1)))

        # Calculate Stds Maps - Set the background to be tiny but non-zero so we can still divide by it
        stds_map_vals = 1e-18*numpy.ones(self.meta.grid_size_1d)
        stds_map_vals.put(masked_idxs, [ps[0] for ps in masked_point_statistics])
        stds_map_vals = flex.double(stds_map_vals.tolist())
        self.statistical_maps.stds_map = stds_map_vals
        # Reshape the array to be the shape of the map
        self.statistical_maps.stds_map.reshape(flex.grid(self.meta.grid_size))

        # Calculate ADJUSTED Stds Maps - Set the background to be tiny but non-zero so we can still divide by it
        sadj_map_vals = 1e-18*numpy.ones(self.meta.grid_size_1d)
        sadj_map_vals.put(masked_idxs, [ps[1] for ps in masked_point_statistics])
        sadj_map_vals = flex.double(sadj_map_vals.tolist())
        self.statistical_maps.sadj_map = sadj_map_vals
        # Reshape the array to be the shape of the map
        self.statistical_maps.sadj_map.reshape(flex.grid(self.meta.grid_size))

        # Calculate Skew Maps
        skew_map_vals = numpy.zeros(self.meta.grid_size_1d)
        skew_map_vals.put(masked_idxs, [ps[2] for ps in masked_point_statistics])
        skew_map_vals = flex.double(skew_map_vals.tolist())
        self.statistical_maps.skew_map = skew_map_vals
        # Reshape the array to be the shape of the map
        self.statistical_maps.skew_map.reshape(flex.grid(self.meta.grid_size))

        # Calculate Kurtosis Maps
        kurt_map_vals = numpy.zeros(self.meta.grid_size_1d)
        kurt_map_vals.put(masked_idxs, [ps[3] for ps in masked_point_statistics])
        kurt_map_vals = flex.double(kurt_map_vals.tolist())
        self.statistical_maps.kurt_map = kurt_map_vals
        # Reshape the array to be the shape of the map
        self.statistical_maps.kurt_map.reshape(flex.grid(self.meta.grid_size))

        # Calculate Bimodality Maps
        bimo_map_vals = numpy.zeros(self.meta.grid_size_1d)
        bimo_map_vals.put(masked_idxs, [ps[4] for ps in masked_point_statistics])
        bimo_map_vals = flex.double(bimo_map_vals.tolist())
        self.statistical_maps.bimo_map = bimo_map_vals
        # Reshape the array to be the shape of the map
        self.statistical_maps.bimo_map.reshape(flex.grid(self.meta.grid_size))

        return self.statistical_maps

    def calculate_z_map(self, method, map=None, tag=None, uncertainty=None):
        """Calculate the z-map relative to the mean and std map"""

        assert method in ['naive','adjusted','uncertainty','adjusted+uncertainty']
        assert [tag, map].count(None) == 1, 'Must provide tag OR map'

        if tag:
            mh = self.dataset_maps.get(tag=tag)
            if not map:
                map = mh.map
            if not uncertainty:
                uncertainty = mh.meta.map_uncertainty

        if 'uncertainty' in method:
            assert uncertainty is not None

        # Calculate Z-values
        if method == 'naive':
            z_map_vals = (map - self.statistical_maps.mean_map)/self.statistical_maps.stds_map
        elif method == 'adjusted':
            z_map_vals = (map - self.statistical_maps.mean_map)/self.statistical_maps.sadj_map
        elif method == 'uncertainty':
            z_map_vals = (map - self.statistical_maps.mean_map)/uncertainty
        elif method == 'adjusted+uncertainty':
            z_map_vals = (map - self.statistical_maps.mean_map)/flex.sqrt(self.statistical_maps.sadj_map**2 + uncertainty**2)
        else:
            raise Exception('method not found: {!s}'.format(method))

        return z_map_vals

class PanddaMultiDatasetAnalyser(Program):
    """Class for the identification of unique "events" in a series of isomorphous-ish crystals"""

    _NAME    = 'pandda.analyse'
    _VERSION = PANDDA_VERSION
    _TEXT    = PANDDA_TEXT

    def __init__(self, params):
        """Class for the identification of unique "events" in a series of isomorphous-ish crystals"""

        # Log init time
        self._init_time = time.time()

        # ===============================================================================>
        # PROCESS INPUT ARGUMENTS
        # ===============================================================================>

        # Read in the master phil
        self.master_phil = pandda_phil
        # Store the whole params object
        self._input_params = params
        # Pull out the python object of the arguments (at the `pandda` level)
        self.args = self._input_params.pandda
        # Most of the variables are contained within the params object - create a shortcut to save typing
        self.params = self._input_params.pandda.params
        # Program settings (num cpus, etc...)
        self.settings = self._input_params.settings

        # ===============================================================================>
        # OUTPUT FILES STUFF
        # ===============================================================================>

        assert self.args.output.out_dir, 'pandda.output.out_dir IS NOT DEFINED'
        self.out_dir = easy_directory(os.path.abspath(self.args.output.out_dir))

        # Create a log for the object
        self.log = Log(log_file=os.path.join(self.out_dir, 'pandda-{}.log'.format(time.strftime("%Y-%m-%d-%H%M", time.gmtime()))), verbose=self.settings.verbose)

        # ===============================================================================>
        # SETTINGS STUFF
        # ===============================================================================>
        self._high_resolution = None
        self._low_resolution = None

        # ===============================================================================>
        # DATA AND MAPS STUFF
        # ===============================================================================>

        # File names
        self._input_files = []
        # Dataset Objects
        self.datasets = DatasetHandlerList()

        # Reference Objects
        self._ref_dataset_index = None
        self._ref_dataset = None
        self._ref_grid = None

        # Variables for storing loaded pickle data from previous runs
        self.pickled_dataset_meta = None

        # Map Statistics
        self.stat_maps = PanddaMultipleStatMapList()

        # ===============================================================================>
        # ANALYSIS OBJECTS
        # ===============================================================================>

        # Create tables object (to hold pandas dataframe objects)
        self.tables = Meta()
        # Record global information about the datasets
        self.tables.dataset_info        = pandas.DataFrame( data    = None,
                                                            index   = pandas.Index(data=[], name='dtag'),
                                                            columns = PanddaTableFields.all_dataset_fields      )
        # Record information about the created maps for each dataset
        self.tables.dataset_map_info    = pandas.DataFrame( data    = None,
                                                            index   = pandas.Index(data=[], name='dtag'),
                                                            columns = PanddaTableFields.all_dataset_map_fields  )
        # Record the events detected in each dataset
        self.tables.event_info          = pandas.DataFrame( data    = None,
                                                            index   = pandas.MultiIndex(levels=[[],[]], labels=[[],[]], names=['dtag','event_idx']),
                                                            columns = PanddaTableFields.all_event_fields        )
        # Record information about the clustered events (cluster of events = site)
        self.tables.site_info           = pandas.DataFrame( data    = None,
                                                            index   = pandas.Index(data=[], name='site_idx'),
                                                            columns = PanddaTableFields.all_site_fields         )

        # ===============================================================================>
        # DIRECTORY STUFF
        # ===============================================================================>

        # Set up the output folders and files
        self._directory_setup()
        # Set up the pickled object filenames
        self._pickle_setup()

        if os.path.exists(self.pickle_handler.get_file('dataset_meta')):
            self._new_pandda = False
        else:
            self._new_pandda = True

    def _directory_setup(self):
        """Initialise the pandda directory system"""

        # Create a file and directory organiser
        self.output_handler = FileManager(rootdir=self.out_dir)

        # Filename templates
        f = PanddaAnalyserFilenames

        # ===============================================================================>
        # Global Directories that do not change from run to run
        # ===============================================================================>

        self.output_handler.add_dir(dir_name='empty_directories',    dir_tag='empty_directories',    top_dir_tag='root', create=False, exists=False)
        self.output_handler.add_dir(dir_name='processed_datasets',   dir_tag='processed_datasets',   top_dir_tag='root', create=False, exists=False)
        self.output_handler.add_dir(dir_name='rejected_datasets',    dir_tag='rejected_datasets',    top_dir_tag='root', create=False, exists=False)
        self.output_handler.add_dir(dir_name='interesting_datasets', dir_tag='interesting_datasets', top_dir_tag='root', create=False, exists=False)
        self.output_handler.add_dir(dir_name='aligned_structures',   dir_tag='aligned_structures',   top_dir_tag='root', create=False, exists=False)
        self.output_handler.add_dir(dir_name='analyses',             dir_tag='analyses',             top_dir_tag='root', create=False, exists=False)

        # ================================================>
        # Input + Status parameters
        # ================================================>
        self.output_handler.add_file(file_name='pandda.eff',        file_tag='settings',         dir_tag='root')
        self.output_handler.add_file(file_name='pandda.{}',         file_tag='status',           dir_tag='root')

        # Store dataset summary graphs
        self.output_handler.add_dir(dir_name='dataset_graph_summaries', dir_tag='d_graphs', top_dir_tag='analyses', create=False, exists=False)
        self.output_handler.add_file(file_name='dataset_resolutions.png',           file_tag='d_resolutions',           dir_tag='d_graphs')
        self.output_handler.add_file(file_name='dataset_rfactors.png',              file_tag='d_rfactors',              dir_tag='d_graphs')
        self.output_handler.add_file(file_name='dataset_global_rmsd_to_ref.png',    file_tag='d_global_rmsd_to_ref',    dir_tag='d_graphs')
        self.output_handler.add_file(file_name='dataset_cell_axes.png',             file_tag='d_cell_axes',             dir_tag='d_graphs')
        self.output_handler.add_file(file_name='dataset_cell_angles.png',           file_tag='d_cell_angles',           dir_tag='d_graphs')
        self.output_handler.add_file(file_name='dataset_cell_volumes.png',          file_tag='d_cell_volumes',          dir_tag='d_graphs')
        # Somewhere to store the dataset information (general values)
        self.output_handler.add_file(file_name=f.dataset_info,                      file_tag='dataset_info',            dir_tag='analyses')
        self.output_handler.add_file(file_name=f.dataset_map_info,                  file_tag='dataset_map_info',        dir_tag='analyses')
        self.output_handler.add_file(file_name=f.dataset_combined_info,             file_tag='dataset_combined_info',   dir_tag='analyses')
        # Somewhere to store the dataset information (identified events + sites)
        self.output_handler.add_file(file_name=f.event_info,                        file_tag='event_info',              dir_tag='analyses')
        self.output_handler.add_file(file_name=f.site_info,                         file_tag='site_info',               dir_tag='analyses')
        self.output_handler.add_file(file_name='_point_distributions.csv',          file_tag='point_distributions',     dir_tag='analyses')
        # Somewhere to store the analysis summaries - for the user
        self.output_handler.add_dir(dir_name='results_summaries', dir_tag='output_summaries', top_dir_tag='root', create=False, exists=False)
        self.output_handler.add_file(file_name=f.initial_html,                      file_tag='initial_html',            dir_tag='output_summaries')
        self.output_handler.add_file(file_name=f.analyse_html,                      file_tag='analyse_html',            dir_tag='output_summaries')
        self.output_handler.add_file(file_name=f.analyse_site_graph,                file_tag='analyse_site_graph',      dir_tag='output_summaries')
        self.output_handler.add_file(file_name=f.analyse_site_graph_mult,           file_tag='analyse_site_graph_mult', dir_tag='output_summaries')
        self.output_handler.add_file(file_name=f.pymol_sites_py,                    file_tag='pymol_sites_py',          dir_tag='output_summaries')
        self.output_handler.add_file(file_name=f.pymol_sites_pml,                   file_tag='pymol_sites_pml',         dir_tag='output_summaries')
        self.output_handler.add_file(file_name=f.pymol_sites_png_1,                 file_tag='pymol_sites_png_1',       dir_tag='output_summaries')
        self.output_handler.add_file(file_name=f.pymol_sites_png_2,                 file_tag='pymol_sites_png_2',       dir_tag='output_summaries')

        # Somewhere to store the pickled objects
        self.output_handler.add_dir(dir_name='pickled_panddas', dir_tag='pickles', top_dir_tag='root', create=False, exists=False)

        # ===============================================================================>
        # Reference Structure Files (should only be needed once for writing and then only for reloading)
        # ===============================================================================>

        # Reference Structure and Dataset
        self.output_handler.add_dir(dir_name='reference', dir_tag='reference', top_dir_tag='root', create=False, exists=False)
        self.output_handler.add_file(file_name=f.reference_structure,               file_tag='reference_structure', dir_tag='reference')
        self.output_handler.add_file(file_name=f.reference_dataset,                 file_tag='reference_dataset',   dir_tag='reference')
        self.output_handler.add_file(file_name=f.reference_on_origin,               file_tag='reference_on_origin', dir_tag='reference')
        self.output_handler.add_file(file_name=f.reference_symmetry,                file_tag='reference_symmetry',  dir_tag='reference')

        # ===============================================================================>
        # Standard template files that will be populated when needed
        # ===============================================================================>

        # Statistical Maps
        self.output_handler.add_dir(dir_name='statistical_maps', dir_tag='statistical_maps', top_dir_tag='root', create=False, exists=False)
        self.output_handler.add_file(file_name=f.mean_map,                          file_tag='mean_map',            dir_tag='statistical_maps')
        self.output_handler.add_file(file_name=f.medn_map,                          file_tag='medn_map',            dir_tag='statistical_maps')
        self.output_handler.add_file(file_name=f.stds_map,                          file_tag='stds_map',            dir_tag='statistical_maps')
        self.output_handler.add_file(file_name=f.sadj_map,                          file_tag='sadj_map',            dir_tag='statistical_maps')
        self.output_handler.add_file(file_name=f.skew_map,                          file_tag='skew_map',            dir_tag='statistical_maps')
        self.output_handler.add_file(file_name=f.kurt_map,                          file_tag='kurt_map',            dir_tag='statistical_maps')
        self.output_handler.add_file(file_name=f.bimo_map,                          file_tag='bimo_map',            dir_tag='statistical_maps')

    def _pickle_setup(self):
        """Initialise all of the pickle filenames"""

        # Pickle Handler
        self.pickle_handler = FileManager(rootdir=self.output_handler.get_dir('pickles'))
        # Pickled Reference Objects
        self.pickle_handler.add_file(file_name='reference_grid.pickle',    file_tag='reference_grid')
        self.pickle_handler.add_file(file_name='reference_dataset.pickle', file_tag='reference_dataset')
        # Pickled Datasets
        self.pickle_handler.add_file(file_name='dataset_masks.pickle',     file_tag='dataset_masks')
        self.pickle_handler.add_file(file_name='dataset_meta.pickle',      file_tag='dataset_meta')
        # Pickled Statistical Maps
        self.pickle_handler.add_file(file_name='statistical_maps.pickle',  file_tag='stat_maps')
        # Map Analysers
        self.pickle_handler.add_file(file_name='map_analyser_{}A.pickle',  file_tag='map_analyser')
        # Pickled Pandda -- (main object)
        self.pickle_handler.add_file(file_name='my_pandda.pickle',         file_tag='my_pandda')

    def run_analysis_init(self):
        """Set up the pandda for a new analysis (doing this will override links to analyses)"""

        # ================================================>
        # Validate the input parameters
        # ================================================>
        self._validate_parameters()

        # ================================================>
        # Create a new analysis directory for analyses/summaries
        # ================================================>
        # New directories will be created for each run (so that data is not overwritten) by the time of the run
        if self.args.output.new_analysis_dir or (not os.path.exists(self.output_handler.get_dir('analyses'))):
            analysis_time_name = 'analyses-{}'.format(time.strftime("%Y-%m-%d-%H%M", time.gmtime(self._init_time)))
            analysis_time_path = easy_directory(os.path.join(self.output_handler.get_dir('root'), analysis_time_name))
            analysis_link_path = self.output_handler.get_dir('analyses')
            # Remove old analysis link if it exists and link in the new analysis directory
            if os.path.exists(analysis_link_path) and os.path.islink(analysis_link_path): os.unlink(analysis_link_path)
            rel_symlink(orig=analysis_time_path, link=analysis_link_path)

        assert os.path.exists(self.output_handler.get_dir('analyses')), 'Output analysis directory does not exist'

        # ===============================================================================>
        # Update the FileManager to make sure all directories are now created
        # ===============================================================================>
        self.output_handler.check_and_create_directories()

        # ===============================================================================>
        # Write the header to the log file
        # ===============================================================================>
        # Print logo and parameters to log
        self.log(PANDDA_TEXT, True)
        self.write_running_parameters_to_log()
        # Write the used parameters to file
        with open(self.output_handler.get_file('settings'), 'w') as out_file:
            out_file.write( '\n'.join([ '# Command Line Args',
                                        '# ',
                                        # This line assures that strings are quoted for easy copy-pasting
                                        '# pandda.analyse '+' '.join(sys.argv[1:]).replace(' ','" ').replace('=','="')+'"',
                                        '',
                                        '# Used Settings:',
                                        '',
                                        self.master_phil.format(python_object=self._input_params).as_str() ]))

        # ===============================================================================>
        # PRINT SOME HELPFUL INFORMATION
        # ===============================================================================>

        self.log('===================================>>>', True)
        self.log('RUNNING FROM: {!s}'.format(sys.argv[0]), True)
        self.log('----------------------------------->>>', True)
        self.log('READING INPUT FROM : {!s}'.format(self.args.input.data_dirs), True)
        self.log('----------------------------------->>>', True)
        self.log('WRITING OUTPUT TO: {!s}'.format(self.out_dir), True)
        self.log('===================================>>>', True)
        self.log('', True)

        # ===============================================================================>
        # LOOK FOR MATPLOTLIB TO SEE IF WE CAN GENERATE GRAPHS
        # ===============================================================================>

        if self.settings.plot_graphs:
            if self.check_for_matplotlib(): pass
            else: self.settings.plot_graphs = False

        # ===============================================================================>
        # CHANGE INTO OUTPUT DIRECTORY
        # ===============================================================================>
        os.chdir(self.out_dir)

        # ===============================================================================>
        # REPOPULATE PANDDA FROM PREVIOUS RUNS
        # ===============================================================================>
        # Load any objects from previous runs
        self.load_pickled_objects()

        # Reload reference dataset
        if (not self.reference_dataset()) and os.path.exists(self.output_handler.get_file('reference_structure')) and os.path.exists(self.output_handler.get_file('reference_dataset')):
            self.log('----------------------------------->>>', True)
            self.log('Loading Reference Dataset', True)
            self.load_reference_dataset(ref_pdb=self.output_handler.get_file('reference_structure'), ref_mtz=self.output_handler.get_file('reference_dataset'))

        # ===============================================================================>
        # LOG THE START TIME
        # ===============================================================================>
        # Store start time and print to log
        self.log('----------------------------------->>>', True)
        self.log('Analysis Started: {!s}'.format(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime(self._init_time))), True)
        # ===============================================================================>
        # Update Status
        # ===============================================================================>
        self.update_status('running')

    def _validate_parameters(self):
        """Validate and preprocess the loaded parameters"""

        self.log('')
        self.log('===================================>>>', True)
        self.log('Validating Input Parameters:', True)
        self.log('===================================>>>', True)
        self.log('')

        p = self.args

        # Input
        assert p.input.data_dirs is not None, 'pandda.input.data_dirs IS NOT DEFINED'
        assert p.input.pdb_style, 'pandda.input.pdb_style IS NOT DEFINED'
        assert p.input.regex.pdb_regex or p.input.regex.mtz_regex or p.input.regex.dir_regex or \
               (p.input.pdb_style and ('*' in p.input.pdb_style)) or \
               (p.input.mtz_style and ('*' in p.input.mtz_style)) or \
               (p.input.data_dirs and ('*' in p.input.data_dirs))
        assert p.input.lig_style, 'pandda.input.lig_style IS NOT DEFINED'

        # Make fullpath so we can run on the eff file from anywhere and change directories without worrying about relative paths
        p.input.data_dirs = os.path.abspath(p.input.data_dirs)
        p.output.out_dir  = os.path.abspath(p.output.out_dir)
        if p.input.filter.pdb:
            p.input.filter.pdb = os.path.abspath(p.input.filter.pdb)

        # If any datasets are set to be reprocessed, reload all datasets (need to change this to allow for "reload_selected_datasets")
        if p.method.reprocess_existing_datasets or self.args.method.reprocess_selected_datasets:
            self.log('----------------------------------->>>', True)
            self.log('Setting method.reload_existing_datasets = True')
            self.log('Old (previously processed) datasets will be reloaded')
            p.method.reload_existing_datasets = True

        if self.is_new_pandda() or self.args.method.reprocess_existing_datasets:
            self.log('----------------------------------->>>', True)
            self.log('Setting output.new_analysis_dir = True')
            self.log('A new output analysis directory will be created')
            p.output.new_analysis_dir = True

    def load_pickled_objects(self):
        """Loads any pickled objects it finds"""

        self.log('----------------------------------->>>', True)
        self.log('Looking for Pickled Files in Input Directory: {!s}'.format(os.path.relpath(self.pickle_handler.get_dir('root'))), True)

        # Record whether any pickled objects are loaded
        pickles_found = False

        # Load Reference Grid
        if os.path.exists(self.pickle_handler.get_file('reference_grid')):
            pickles_found = True
            self.log('-> Loading reference grid')
            self.set_reference_grid(self.unpickle(self.pickle_handler.get_file('reference_grid')))

        # Load Reference Dataset
        if os.path.exists(self.pickle_handler.get_file('reference_dataset')):
            pickles_found = True
            self.log('-> Loading reference dataset')
            self.set_reference_dataset(self.unpickle(self.pickle_handler.get_file('reference_dataset')))

        # Load the datasets
        if os.path.exists(self.pickle_handler.get_file('dataset_meta')):
            pickles_found = True
            # Unpickle the list of the pickled datasets from the directory structure
            self.log('-> Loading old dataset information (existing datasets)')
            self.pickled_dataset_meta = self.unpickle(self.pickle_handler.get_file('dataset_meta'))

            if self.args.method.reload_existing_datasets:
                # Extract the paths of the pickled dataset objects
                pickled_dataset_list = self.pickled_dataset_meta.dataset_pickle_list
                # Check they all exist - should be relative to the outdirectory
                for filename in pickled_dataset_list:
                    assert os.path.isfile(os.path.join(self.out_dir, filename)), 'File does not exist: {!s}'.format(filename)
                # Unpickle the datasets and add them to the dataset handler list
                self.log('-> Reloading old datasets')
                self.datasets.add([self.unpickle(os.path.join(self.out_dir,f)) for f in pickled_dataset_list])
            else:
                self.log('-> Not reloading old datasets')
        else:
            # No datasets to load - this must be False
            self.args.method.reload_existing_datasets = False
            self.log('-> No old datasets found')

        # Load Statistical Maps
        if os.path.exists(self.pickle_handler.get_file('stat_maps')):
            pickles_found = True
            self.log('-> Loading old statistical maps')
            self.stat_maps = self.unpickle(self.pickle_handler.get_file('stat_maps'))

        if not pickles_found:
            self.log('-> No Pickles Found', True)

    def pickle_the_pandda(self, components=None, all=False, datasets=None):
        """Pickles it's major components for quick loading..."""

        if all == True:
            self.log('----------------------------------->>>', True)
            self.log('Pickling the Pandda', True)
        elif not components:
            self.log('----------------------------------->>>', True)
            self.log('Pickling NOTHING', True)
            return
        else:
            self.log('----------------------------------->>>', True)
            self.log('Selective Pickling: {!s}'.format(', '.join(components).upper()), True)

        if all or ('grid' in components):
            self.log('----------------------------------->>>')
            if self.reference_grid() is not None:
                self.log('Pickling Reference Grid')
                self.pickle(pickle_file   = self.pickle_handler.get_file('reference_grid'),
                            pickle_object = self.reference_grid(),
                            overwrite = False)
            else:
                self.log('No Reference Grid to Pickle')

        if all or ('datasets' in components):
            self.log('----------------------------------->>>')

            if self.reference_dataset():
                self.log('Pickling Reference Dataset')

                self.pickle(pickle_file   = self.pickle_handler.get_file('reference_dataset'),
                            pickle_object = self.reference_dataset().get_pickle_copy(),
                            overwrite = True)

            # If no datasets given, pickle them all
            if not datasets:
                datasets = self.datasets.all()
            # Pickle the datasets (individual pickle files)
            if datasets:
                self.log('Pickling Datasets')
                # Need to get a 'pickle_copy' as some of the dataset handler may not be pickleable
                for d_handler in datasets:
                    self.pickle(pickle_file   = d_handler.output_handler.get_file('dataset_pickle'),
                                pickle_object = d_handler.get_pickle_copy(),
                                overwrite = True)
            else:
                self.log('No Datasets to Pickle')

        if all or ('stat_maps' in components):
            self.log('----------------------------------->>>')
            if self.stat_maps is not None:
                self.log('Pickling Statistical Maps')
                self.pickle(pickle_file   = self.pickle_handler.get_file('stat_maps'),
                            pickle_object = self.stat_maps,
                            overwrite = True)
            else:
                self.log('No Statistical Maps to Pickle')

    def exit(self, error=False):
        """Exit the PANDDA, record runtime etc..."""

        self._finish_time = time.time()
        self.log('Runtime: {!s}'.format(time.strftime("%H hours:%M minutes:%S seconds", time.gmtime(self._finish_time - self._init_time))))

        # If error, don't make meta or self pickle
        if error:
            self.update_status('errored')
            self.log('', True)
            self.log('===================================>>>', True)
            self.log('PANDDA exited with an error')
            self.log('----------------------------------->>>', True)
            self.log('Error Traceback: ')
            self.log('----------------------------------->>>', True)
            self.log(traceback.format_exc())
            self.log('===================================>>>', True)
            return
        else:
            self.update_status('done')
            self.log('', True)
            self.log('===================================>>>', True)
            self.log('.. FINISHED! PANDDA EXITED NORMALLY ..', True)
            self.log('===================================>>>', True)

        try:
            # Extract meta about the datasets
            if self.pickled_dataset_meta and (not self.args.method.reload_existing_datasets):
                # Combine the old meta with the new
                self.log('Combining old dataset meta with new meta for pickle')
                number_of_datasets  = self.pickled_dataset_meta.number_of_datasets  + self.datasets.size()
                dataset_labels      = self.pickled_dataset_meta.dataset_labels      + [d.tag for d in self.datasets.all()]
                dataset_pickle_list = self.pickled_dataset_meta.dataset_pickle_list + [os.path.relpath(d.output_handler.get_file('dataset_pickle'), start=self.out_dir) for d in self.datasets.all()]
            else:
                # Don't need to consider the pickle as the datasets have been reloaded
                self.log('Creating new meta for pickle')
                number_of_datasets  = self.datasets.size()
                dataset_labels      = [d.tag for d in self.datasets.all()]
                dataset_pickle_list = [os.path.relpath(d.output_handler.get_file('dataset_pickle'), start=self.out_dir) for d in self.datasets.all()]
            # Create a dictionary to be stored
            dataset_meta = Meta({'number_of_datasets'    : number_of_datasets,
                                 'dataset_labels'        : dataset_labels,
                                 'dataset_pickle_list'   : dataset_pickle_list
                            })
            # Pickle the list of locations of the dataset pickles
            self.pickle(pickle_file=self.pickle_handler.get_file('dataset_meta'), pickle_object=dataset_meta, overwrite=True)
        except:
            self.log('FAILED TO PICKLE META')

        # Lastly, pickle myself (if required)
        try:
            if self.args.output.pickling.pickle_complete_pandda:
                # Pickle myself
                self.log('----------------------------------->>>', True)
                self.log('Pickling the PANDDA Results')
                self.pickle(pickle_file=self.pickle_handler.get_file('my_pandda'), pickle_object=self, overwrite=True)
        except:
            self.log('FAILED TO PICKLE MYSELF')

    def is_new_pandda(self):
        """Is this the first time the program has been run?"""
        return self._new_pandda

    def set_high_resolution(self, res):
        self._high_resolution = res
    def get_high_resolution(self):
        return self._high_resolution

    def set_low_resolution(self, res):
        self._low_resolution = res
    def get_low_resolution(self):
        return self._low_resolution

    def set_reference_dataset(self, dataset):
        """Set a reference dataset created externally"""
        self._ref_dataset = dataset
    def reference_dataset(self):
        """Get the dataset used as the reference for scaling and aligning the other datasets"""
        return self._ref_dataset
    def set_reference_grid(self, ref_grid):
        """Set a grid created externally"""
        self._ref_grid = ref_grid
    def reference_grid(self):
        """Get the reference grid used to generate the points for sampling the maps"""
        return self._ref_grid

    def new_files(self):
        """Get all of the files that were added on this run"""
        return self._input_files

    def initialise_dataset_masks_and_tables(self):
        """Add blank masks to the mask objects, based on how many datasets have been loaded"""

        self.log('----------------------------------->>>', True)
        self.log('Initialising Dataset Masks.', True)

        # Set the dataset mask lengths and ids (dataset tags)
        self.datasets.all_masks().set_mask_length(mask_length=self.datasets.size())
        self.datasets.all_masks().set_entry_ids(entry_ids=[d.tag for d in self.datasets.all()])
        # Initialise standard blank masks
        for mask_name in PanddaMaskNames.all_mask_names:
            self.datasets.all_masks().add_mask(mask_name=mask_name, mask=[False]*self.datasets.size())

        # Initialise masks for datasets that shouldn't be analysed
        if self.args.input.flags.no_analyse:
            no_analyse_tags = self.args.input.flags.no_analyse.split(',')
            self.log('Not analysing {!s} Datasets: {!s}'.format(len(no_analyse_tags), ', '.join(no_analyse_tags)))
            no_analyse_mask = [True if d.tag in no_analyse_tags else False for d in self.datasets.all()]
            self.datasets.all_masks().add_mask(mask_name='no_analyse', mask=no_analyse_mask)

        # Initialise mask for datasets that shouldn't be used for building
        if self.args.input.flags.no_build:
            no_build_tags = self.args.input.flags.no_build.split(',')
            self.log('Not building distributions from {!s} Datasets: {!s}'.format(len(no_build_tags), ', '.join(no_build_tags)))
            no_build_mask = [True if d.tag in no_build_tags else False for d in self.datasets.all()]
            self.datasets.all_masks().add_mask(mask_name='no_build', mask=no_build_mask)

        # Initialise mask for datasets that have been previously pickled
        self.datasets.all_masks().add_mask(mask_name='old datasets', mask=[False]*self.datasets.size())

        if self.pickled_dataset_meta and self.args.method.reload_existing_datasets:
            for tag in self.pickled_dataset_meta.dataset_labels:
                self.datasets.all_masks().set_mask_value(mask_name='old datasets', entry_id=tag, value=True)
            self.log('Considering {!s} datasets as "New Datasets"'.format(self.datasets.size(mask_name='old datasets', invert=True)))
            self.log('Considering {!s} datasets as "Old Datasets"'.format(self.datasets.size(mask_name='old datasets')))
        else:
            self.log('Considering all {!s} datasets as "New Datasets"'.format(self.datasets.size(mask_name='old datasets', invert=True)))
            assert self.datasets.size(mask_name='old datasets', invert=True) == self.datasets.size(), 'Total datasets should be same as total new datasets'

        self.log('----------------------------------->>>', True)
        self.log('Initialising dataset data-tables.', True)

        # Add dataset tags as rows in the tables
        self.tables.dataset_info     = self.tables.dataset_info.append(pandas.DataFrame(index=[d.tag for d in self.datasets.all()]), verify_integrity=True)
        self.tables.dataset_map_info = self.tables.dataset_map_info.append(pandas.DataFrame(index=[d.tag for d in self.datasets.all()]), verify_integrity=True)

        old_datasets = self.datasets.mask(mask_name='old datasets')
        if (not self.args.method.reprocess_existing_datasets) and old_datasets:
            self.log('Syncing old dataset information to dataset tables.', True)
            self.sync_datasets(datasets=old_datasets)
            self.log('Syncing old dataset events to output tables.', True)
            for d_handler in old_datasets:
                if d_handler.events:
                    for e in d_handler.events:
                        self.add_event_to_event_table(d_handler=d_handler, event=e)

    def select_reference_dataset(self, method='resolution', max_rfree=0.4, min_resolution=5):
        """Select dataset to act as the reference - scaling, aligning etc"""

        assert method in ['resolution','rfree'], 'METHOD FOR SELECTING THE REFERENCE DATASET NOT RECOGNISED: {!s}'.format(method)

        # Create a mask of the datasets that can be selected as the reference dataset
        no_build_mask = self.datasets.all_masks().get_mask(mask_name='no_build')
        rejected_mask = self.datasets.all_masks().get_mask(mask_name='rejected - total')
        potential_reference_mask = self.datasets.all_masks().combine_masks_custom(masks=[no_build_mask,rejected_mask], invert=True)
        self.datasets.all_masks().add_mask(mask_name='potential reference datasets', mask=potential_reference_mask)
        # Get the potential reference datasets
        filtered_datasets = self.datasets.mask(mask_name='potential reference datasets')
        if not filtered_datasets: raise Exception('NO NON-REJECTED DATASETS REMAINING')

        self.log('---------->>>', True)
        self.log('Selecting Reference Dataset by: {!s}'.format(method), True)
        if method == 'rfree':
            # Get RFrees of datasets (set to dummy value of 999 if resolution is too high so that it is not selected)
            r_frees = [d.input().get_r_rfree_sigma().r_free if (d.reflection_data().max_min_resolution()[1] < min_resolution) else 999 for d in filtered_datasets]
            if len(resolns) == 0: raise Exception('NO DATASETS BELOW RESOLUTION CUTOFF {!s}A - CANNOT SELECT REFERENCE DATASET'.format(min_resolution))
            ref_dataset_index = r_frees.index(min(r_frees))
        elif method == 'resolution':
            # Get Resolutions of datasets (set to dummy value of 999 if r-free is too high so that it is not selected)
            resolns = [d.reflection_data().max_min_resolution()[1] if (d.input().get_r_rfree_sigma().r_free < max_rfree) else 999 for d in filtered_datasets]
            if len(resolns) == 0: raise Exception('NO DATASETS BELOW RFREE CUTOFF {!s} - CANNOT SELECT REFERENCE DATASET'.format(max_rfree))
            ref_dataset_index = resolns.index(min(resolns))

        reference = filtered_datasets[ref_dataset_index]
        self._ref_dataset_index = reference.num
        self.log('Reference Selected: {!s}'.format(reference.tag), True)
        self.log('Resolution: {!s}, RFree: {!s}'.format(reference.reflection_data().max_min_resolution()[1], reference.input().get_r_rfree_sigma().r_free), True)

        return reference.pdb_filename(), reference.mtz_filename()

    def load_reference_dataset(self, ref_pdb, ref_mtz):
        """Set the reference dataset, to which all other datasets will be aligned and scaled"""

        self.log('---------->>>', True)
        self.log('Loading Reference Dataset: {!s}'.format(ref_mtz), True)

        link_ref_pdb = self.output_handler.get_file('reference_structure')
        link_ref_mtz = self.output_handler.get_file('reference_dataset')

        # Remove old links?
        if os.path.abspath(ref_pdb) != os.path.abspath(link_ref_pdb):
            if os.path.exists(link_ref_pdb): os.unlink(link_ref_pdb)
            if os.path.exists(link_ref_mtz): os.unlink(link_ref_mtz)
        # Create links to dataset
        if not os.path.exists(link_ref_pdb): rel_symlink(orig=ref_pdb, link=link_ref_pdb)
        if not os.path.exists(link_ref_mtz): rel_symlink(orig=ref_mtz, link=link_ref_mtz)

        # Create and set reference dataset
        ref_handler = ReferenceDatasetHandler(
                        dataset_number = -1,
                        pdb_filename   = os.path.relpath(link_ref_pdb, start=self.out_dir),
                        mtz_filename   = os.path.relpath(link_ref_mtz, start=self.out_dir),
                        dataset_tag    = 'reference')
        self.set_reference_dataset(ref_handler)

        # Calculate the shift required to move the reference structure into the positive quadrant
        total_border_padding = self.params.maps.padding + self.params.masks.outer_mask
        self.reference_dataset().set_origin_shift(tuple(flex.double(3, total_border_padding) - flex.double(self.reference_dataset().input().atoms().extract_xyz().min())))
        self.log('Origin Shift for reference structure: {!s}'.format(tuple([round(s,3) for s in self.reference_dataset().origin_shift()])))
        # Shift the reference structure by this amount so that it is aligned with the reference grid
        ref_hierarchy = self.reference_dataset().hierarchy()
        ref_hierarchy.atoms().set_xyz(ref_hierarchy.atoms().extract_xyz() + self.reference_dataset().origin_shift())

        # Save this output -- this is essentially defines the reference frame
        if not os.path.exists(self.output_handler.get_file('reference_on_origin')):
            ref_hierarchy.write_pdb_file(self.output_handler.get_file('reference_on_origin'))

        # Align reference dataset to it's shifted self
        r = scitbx.matrix.rec([1,0,0,0,1,0,0,0,1], (3,3))
        t = scitbx.matrix.rec(self.reference_dataset().origin_shift(), (3,1))
        rt = scitbx.matrix.rt((r,t))
        self.reference_dataset().set_global_alignment(alignment=rt)

        # Create neighbouring symmetry copies of the reference structures
        ref_sym_copies = self.reference_dataset().generate_symmetry_copies(rt_method='global', save_operators=True, buffer=self.params.masks.outer_mask+5)
        # Write out the symmetry sites
        if not os.path.exists(self.output_handler.get_file('reference_symmetry')):
            ref_sym_copies.write_pdb_file(self.output_handler.get_file('reference_symmetry'))

        return self.reference_dataset()

    def create_reference_grid(self, grid_spacing, expand_to_origin, buffer=0):
        """Create a grid over the reference protein"""

        self.log('----------------------------------->>>', True)
        self.log('Creating Reference Grid', True)

        sites_cart = self.reference_dataset().input().atoms().extract_xyz()

        assert (flex.vec3_double(sorted(self.reference_dataset().input().atoms().extract_xyz())) -
                flex.vec3_double(sorted(self.reference_dataset().hierarchy().atoms().extract_xyz()))).dot().norm() == 0.0, 'EH? Coordinates should be the same?'

        self._ref_grid = grid_handler(verbose=self.settings.verbose)
        self._ref_grid.set_grid_spacing(spacing=grid_spacing)
        self._ref_grid.set_cart_extent(cart_min=tuple([s-buffer for s in sites_cart.min()]), cart_max=tuple([s+buffer for s in sites_cart.max()]))
        self._ref_grid.create_cartesian_grid(expand_to_origin=expand_to_origin)
        self.log(self._ref_grid.summary())

        if self.params.alignment.method == 'local':
            self.log('---------->>>', True)
            self.log('Partitioning Reference Grid', True)

            # Pull out the calphas
            calpha_hierarchy = self.reference_dataset().calphas()

            t1 = time.time()
            # Calculate the nearest residue for each point on the grid
            self.reference_grid().create_grid_partition(atomic_hierarchy=calpha_hierarchy)
            # Partition Grid
            self.reference_grid().partition().partition_grid(cpus=self.settings.cpus)
            t2 = time.time()
            self.log('> GRID PARTITIONING > Time Taken: {!s} seconds'.format(int(t2-t1)))

    def mask_reference_grid(self):
        """Create masks for the reference grid based on distances from atoms in the reference structure"""

        self.log('----------------------------------->>>', True)
        self.log('Masking Reference Grid', True)

        # ============================================================================>
        # Get neighbouring symmetry copies of the reference structures
        # ============================================================================>
        ref_sym_copies = self.reference_dataset().symmetry_copies
        # ============================================================================>
        # Local mask used for forming groups of points around a grid point
        # ============================================================================>
        if self.reference_grid().local_mask() is None:
            self.log('---------->>>', True)
            self.log('Generating Local Mask')
            local_mask = spherical_mask(grid_spacing    = self.reference_grid().grid_spacing(),
                                        distance_cutoff = 1.2,
                                        grid_jump       = 1 )
            self.reference_grid().set_local_mask(local_mask)
        # ============================================================================>
        # Global mask used for removing points in the bulk solvent regions
        # ============================================================================>
        if self.reference_grid().global_mask() is None:
            self.log('---------->>>', True)
            self.log('Generating Protein Mask')
            # Select the masking atoms from the reference structure
            cache = self.reference_dataset().hierarchy().atom_selection_cache()
            pro_sites_cart = self.reference_dataset().hierarchy().select(cache.selection('pepnames and not element H')).atoms().extract_xyz()
            # Generate the main protein mask
            global_mask = atomic_mask(  cart_sites   = pro_sites_cart,
                                        grid_size    = self.reference_grid().grid_size(),
                                        unit_cell    = self.reference_grid().unit_cell(),
                                        max_dist     = self.params.masks.outer_mask,
                                        min_dist     = self.params.masks.inner_mask )
            self.reference_grid().set_global_mask(global_mask)
        # ============================================================================>
        # Global mask used for removing points close to symmetry copies of the protein
        # ============================================================================>
        if self.reference_grid().symmetry_mask() is None:
            self.log('---------->>>', True)
            self.log('Generating Symmetry Mask')
            # Pull out the cartesian sites of the symmetry mates
            cache = ref_sym_copies.atom_selection_cache()
            sym_sites_cart = ref_sym_copies.select(cache.selection('pepnames and not element H')).atoms().extract_xyz()
            # Generate the symmetry mask
            symmetry_mask = grid_mask(cart_sites = sym_sites_cart,
                                      grid_size  = self.reference_grid().grid_size(),
                                      unit_cell  = self.reference_grid().unit_cell(),
                                      max_dist   = self.params.masks.outer_mask,
                                      min_dist   = self.params.masks.inner_mask_symmetry )
            self.reference_grid().set_symmetry_mask(symmetry_mask)

        # ============================================================================>
        # Write masked maps
        # ============================================================================>
        # Write protein masked map
        mask_map_file = self.output_handler.get_file('reference_dataset').replace('.mtz','.totalmask.ccp4')
        map_mask = flex.double(self.reference_grid().global_mask().total_mask_binary().astype(int))
        map_mask.reshape(flex.grid(self.reference_grid().grid_size()))
        write_array_to_map(output_file=mask_map_file, map_data=map_mask, grid=self.reference_grid())
        # Write symmetry masked map
        mask_map_file = self.output_handler.get_file('reference_dataset').replace('.mtz','.symmask.ccp4')
        map_mask = flex.double(self.reference_grid().symmetry_mask().total_mask_binary().astype(int))
        map_mask.reshape(flex.grid(self.reference_grid().grid_size()))
        write_array_to_map(output_file=mask_map_file, map_data=map_mask, grid=self.reference_grid())

        return

    def build_input_list(self):
        """Builds a list of input files from the command line arguments passed"""

        self.log('----------------------------------->>>', True)
        self.log('Building List of Datasets')

        dir_style = self.args.input.data_dirs.strip('./')
        pdb_style = self.args.input.pdb_style.lstrip('/')
        if self.args.input.mtz_style:
            mtz_style = self.args.input.mtz_style.lstrip('/')
        else:
            assert pdb_style.endswith('.pdb'), 'pdb_style does not end in .pdb'
            mtz_style = pdb_style.replace('.pdb','.mtz')

        # Datasets that are already added
        new_files = []
        empty_directories = []

        for dir in sorted(glob.glob(self.args.input.data_dirs)):
            pdb_files = [f for f in glob.glob(os.path.join(dir, pdb_style)) if os.path.exists(f)]
            mtz_files = [f for f in glob.glob(os.path.join(dir, mtz_style)) if os.path.exists(f)]
            if not (pdb_files and mtz_files):
                print('EMPTY DIRECTORY: {!s}'.format(dir))
                empty_directories.append(dir)
            elif not pdb_files:
                print('NO PDB IN DIRECTORY: {!s}'.format(dir))
                empty_directories.append(dir)
            elif not mtz_files:
                print('NO MTZ IN DIRECTORY: {!s}'.format(dir))
                empty_directories.append(dir)
            else:
                assert len(pdb_files) == 1, 'More than one matching PDB file found: {!s}'.format(os.path.join(dir, pdb_style))
                assert len(mtz_files) == 1, 'More than one matching MTZ file found: {!s}'.format(os.path.join(dir, mtz_style))

                new_pdb = pdb_files[0]
                new_mtz = mtz_files[0]
                dataset_tag = [None]

                # Do regex matching on the file pairs
                if '*' in pdb_style:
                    pdb_base = os.path.basename(new_pdb)
                    if self.args.input.regex.pdb_regex:
                        pdb_regex = self.args.input.regex.pdb_regex
                    else:
                        pdb_regex = pdb_style.replace('*', '(.*)')
                    pdb_tag = re.findall(pdb_regex, pdb_base)
                    assert pdb_tag, 'NO PDB TAG FOUND: {!s} -> {!s}'.format(pdb_regex, pdb_base)
                    if isinstance(pdb_tag[0], tuple):
                        self.log('More than one PDB TAG found - choosing the first one of {!s}'.format(pdb_tag[0]))
                        pdb_tag = list(pdb_tag[0])[0:1]
                else: pdb_regex = pdb_tag = None

                if '*' in mtz_style:
                    mtz_base = os.path.basename(new_mtz)
                    if self.args.input.regex.mtz_regex:
                        mtz_regex = self.args.input.regex.mtz_regex
                    else:
                        mtz_regex = mtz_style.replace('*', '(.*)')
                    mtz_tag = re.findall(mtz_regex, mtz_base)
                    assert mtz_tag, 'NO MTZ TAG FOUND: {!s} -> {!s}'.format(mtz_regex, mtz_base)
                    if isinstance(mtz_tag[0], tuple):
                        self.log('More than one MTZ TAG found - choosing the first one of {!s}'.format(mtz_tag[0]))
                        mtz_tag = list(mtz_tag[0])[0:1]
                else: mtz_regex = mtz_tag = None

                if '*' in dir_style:
                    dir_base = os.path.dirname(pdb_files[0])
                    if self.args.input.regex.dir_regex:
                        dir_regex = self.args.input.regex.dir_regex
                    else:
                        dir_regex = dir_style.replace('*', '(.*)')
                    dir_tag = re.findall(dir_regex, dir_base)
                    assert dir_tag, 'NO DIR TAG FOUND: {!s} -> {!s}'.format(dir_regex, dir_base)
                    if isinstance(dir_tag[0], tuple):
                        self.log('More than one DIR TAG found - choosing the first one of {!s}'.format(dir_tag[0]))
                        dir_tag = list(dir_tag[0])[0:1]
                else: dir_regex = dir_tag = None

                if pdb_tag and mtz_tag: assert pdb_tag == mtz_tag, 'PDB-MTZ TAGS ARE NOT IDENTICAL'
                if dir_tag and pdb_tag: assert dir_tag == pdb_tag, 'DIR-PDB TAGS ARE NOT IDENTICAL'
                if dir_tag and mtz_tag: assert dir_tag == mtz_tag, 'DIR-MTZ TAGS ARE NOT IDENTICAL'

                if   dir_tag: dataset_tag = dir_tag
                elif pdb_tag: dataset_tag = pdb_tag
                elif mtz_tag: dataset_tag = mtz_tag

                # Add prefix
                if isinstance(dataset_tag[0], str): dataset_tag = [self.args.output.dataset_prefix + dataset_tag[0]]
                else:                               assert dataset_tag[0] is None

                new_files.append(pdb_files+mtz_files+dataset_tag)

        # Filter out the already added files
        if self.pickled_dataset_meta:
            filtered_new_files = []
            for i, (pdb, mtz, tag) in enumerate(new_files):
                if tag in self.pickled_dataset_meta.dataset_labels:
                    self.log('Dataset with this tag has already been loaded: {!s} - Not loading'.format(tag))
                else:
                    filtered_new_files.append(new_files[i])
        else:
            filtered_new_files = new_files

        # Filter out manually labelled datasets to ignore
        if self.args.input.flags.ignore_datasets:
            ignore_tags = self.args.input.flags.ignore_datasets.split(',')
            self.log('Ignoring {!s} Datasets: {!s}'.format(len(ignore_tags), ', '.join(ignore_tags)))
            re_filtered_new_files = []
            for i, (pdb, mtz, tag) in enumerate(filtered_new_files):
                if tag in ignore_tags:
                    self.log('Ignoring Dataset: {!s}'.format(tag))
                else:
                    re_filtered_new_files.append(filtered_new_files[i])
            filtered_new_files = re_filtered_new_files

        # Get the list of already linked empty_directories
        empty_dir_prefix = 'Dir_'
        link_old_empty_dirs = glob.glob(os.path.join(self.output_handler.get_dir('empty_directories'),'*'))
        real_old_empty_dirs = [os.path.realpath(p) for p in link_old_empty_dirs]
        # Pull out the highest current idx
        emp_num_i = 0
        emp_num_offset = max([0]+[int(os.path.basename(v).strip(empty_dir_prefix)) for v in link_old_empty_dirs])
        assert emp_num_offset == len(link_old_empty_dirs), 'Numbering of empty directories is not consecutive'
        # Link the empty directories into the same directory
        for dir in empty_directories:
            # Already linked
            if os.path.realpath(dir) in real_old_empty_dirs:
                continue
            # Increment counter
            emp_num_i += 1
            # Create new dir
            empty_dir = os.path.join(self.output_handler.get_dir('empty_directories'), empty_dir_prefix+'{:05d}'.format(emp_num_i+emp_num_offset))
            if not os.path.exists(empty_dir):
                os.symlink(dir, empty_dir)
            else:
                raise Exception('THIS DIRECTORY SHOULD NOT EXIST: {!s}'.format(empty_dir))

        # Record number of empty datasets
        self.log('---------->>>', True)
        self.log('{!s} EMPTY DIRECTORIES FOUND (TOTAL)'.format(emp_num_i+emp_num_offset), True)
        self.log('{!s} EMPTY DIRECTORIES FOUND (NEW)'.format(emp_num_i), True)

        # Record total number of datasets, and total number of new datasets
        if self.pickled_dataset_meta: num_old = self.pickled_dataset_meta.number_of_datasets
        else:                         num_old = 0
        self.log('---------->>>', True)
        self.log('{!s} DATASETS FOUND (TOTAL)'.format(len(filtered_new_files)+num_old), True)
        self.log('{!s} DATASETS FOUND (NEW)'.format(len(filtered_new_files)), True)

        return filtered_new_files

    def add_new_files(self, input_files):
        """Add (pdb, mtz) file pairs to the datasets to be processed"""

        self._input_files = input_files

        self.log('----------------------------------->>>', True)
        self.log('{!s} Datasets Added'.format(len(input_files)), True)

    def load_new_datasets(self):
        """Read in maps for the input datasets"""

        if not self.datasets.all() and self.is_new_pandda():
            self.log('Adding First Datasets to Pandda')
        else:
            self.log('Adding more datasets')
            self.log('{!s} already loaded'.format(self.datasets.size()))
            self.log('{!s} not loaded'.format(self.pickled_dataset_meta.number_of_datasets - self.datasets.size()))

        # Counting offset for dataset index
        if self.pickled_dataset_meta: n_offset = self.pickled_dataset_meta.number_of_datasets
        else:                         n_offset = 0

        # Generate arg_list for loading
        arg_list = [{'dataset_number':dnum+n_offset, 'pdb_filename':pdb, 'mtz_filename':mtz, 'dataset_tag':dtag} for dnum, (pdb, mtz, dtag) in enumerate(self.new_files())]

        start = time.time()
        self.log('----------------------------------->>>', True)
        print('Adding Datasets... (using {!s} cores)'.format(self.settings.cpus))
        loaded_datasets = easy_mp.pool_map(func=load_dataset_map_func, args=arg_list, processes=self.settings.cpus)
        finish = time.time()
        self.log('> Adding Datasets > Time Taken: {!s} seconds'.format(int(finish-start)), True)

        lig_style = self.args.input.lig_style.strip('/')

        # Output Path Templates
        f = PanddaDatasetFilenames
        p = PanddaDatasetPNGFilenames

        for d_handler in loaded_datasets:

            # Intialise the meta for the dataset
            d_handler.meta.analysed = False
            d_handler.meta.dataset_info = None
            d_handler.meta.dataset_map_info = None

            # Create a file manager object
            d_handler.initialise_output_directory(outputdir=os.path.join(self.output_handler.get_dir('processed_datasets'), d_handler.tag))

            # Main input/output files
            d_handler.output_handler.add_file(file_name=f.input_structure.format(d_handler.tag),                    file_tag='input_structure'              )
            d_handler.output_handler.add_file(file_name=f.input_data.format(d_handler.tag),                         file_tag='input_data'                   )
            d_handler.output_handler.add_file(file_name=f.dataset_info.format(d_handler.tag),                       file_tag='dataset_info'                 )
            d_handler.output_handler.add_file(file_name=f.dataset_log.format(d_handler.tag),                        file_tag='dataset_log'                  )
            d_handler.output_handler.add_file(file_name=f.aligned_structure.format(d_handler.tag),                  file_tag='aligned_structure'            )
            d_handler.output_handler.add_file(file_name=f.symmetry_copies.format(d_handler.tag),                    file_tag='symmetry_copies'              )
            d_handler.output_handler.add_file(file_name=f.sampled_map.format(d_handler.tag),                        file_tag='sampled_map'                  )
            d_handler.output_handler.add_file(file_name=f.mean_diff_map.format(d_handler.tag),                      file_tag='mean_diff_map'                )
            d_handler.output_handler.add_file(file_name=f.z_map.format(d_handler.tag),                              file_tag='z_map'                        )
            d_handler.output_handler.add_file(file_name=f.z_map_naive.format(d_handler.tag),                        file_tag='z_map_naive'                  )
            d_handler.output_handler.add_file(file_name=f.z_map_naive_norm.format(d_handler.tag),                   file_tag='z_map_naive_normalised'       )
            d_handler.output_handler.add_file(file_name=f.z_map_uncertainty.format(d_handler.tag),                  file_tag='z_map_uncertainty'            )
            d_handler.output_handler.add_file(file_name=f.z_map_uncertainty_norm.format(d_handler.tag),             file_tag='z_map_uncertainty_normalised' )
            d_handler.output_handler.add_file(file_name=f.z_map_corrected.format(d_handler.tag),                    file_tag='z_map_corrected'              )
            d_handler.output_handler.add_file(file_name=f.z_map_corrected_norm.format(d_handler.tag),               file_tag='z_map_corrected_normalised'   )
            d_handler.output_handler.add_file(file_name=f.event_map.format(d_handler.tag, '{!s}', '{!s}'),          file_tag='event_map'                    )

            # Miscellaneous files
            d_handler.output_handler.add_file(file_name=f.high_z_mask.format(d_handler.tag), file_tag='high_z_mask')
            d_handler.output_handler.add_file(file_name=f.grid_mask.format(d_handler.tag),   file_tag='grid_mask')

            # Links to ligand files (if they've been found)
            d_handler.output_handler.add_dir(dir_name='ligand_files', dir_tag='ligand', top_dir_tag='root')
#            d_handler.output_handler.add_file(file_name=f.ligand_coordinates.format(d_handler.tag),                 file_tag='ligand_coordinates'   )
#            d_handler.output_handler.add_file(file_name=f.ligand_restraints.format(d_handler.tag),                  file_tag='ligand_restraints'    )
#            d_handler.output_handler.add_file(file_name=f.ligand_image.format(d_handler.tag),                       file_tag='ligand_image'         )

            # Native (back-rotated/transformed) maps
            d_handler.output_handler.add_file(file_name=f.native_obs_map.format(d_handler.tag),                     file_tag='native_obs_map'       )
            d_handler.output_handler.add_file(file_name=f.native_z_map.format(d_handler.tag),                       file_tag='native_z_map'         )
            d_handler.output_handler.add_file(file_name=f.native_event_map.format(d_handler.tag,'{!s}','{!s}'),     file_tag='native_event_map'     )

            # Fitted structures when modelled with pandda.inspect
            d_handler.output_handler.add_dir(dir_name='modelled_structures', dir_tag='models', top_dir_tag='root')

            # Output images
            d_handler.output_handler.add_dir(dir_name='output_images', dir_tag='images', top_dir_tag='root')
            # Smapled map
            d_handler.output_handler.add_file(file_name=p.s_map_png.format(d_handler.tag),                          file_tag='s_map_png',                        dir_tag='images')
            d_handler.output_handler.add_file(file_name=p.d_mean_map_png.format(d_handler.tag),                     file_tag='d_mean_map_png',                   dir_tag='images')
            d_handler.output_handler.add_file(file_name=p.z_map_naive_png.format(d_handler.tag),                    file_tag='z_map_naive_png',                  dir_tag='images')
            d_handler.output_handler.add_file(file_name=p.z_map_naive_norm_png.format(d_handler.tag),               file_tag='z_map_naive_normalised_png',       dir_tag='images')
            d_handler.output_handler.add_file(file_name=p.z_map_uncertainty_png.format(d_handler.tag),              file_tag='z_map_uncertainty_png',            dir_tag='images')
            d_handler.output_handler.add_file(file_name=p.z_map_uncertainty_norm_png.format(d_handler.tag),         file_tag='z_map_uncertainty_normalised_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name=p.z_map_corrected_png.format(d_handler.tag),                file_tag='z_map_corrected_png',              dir_tag='images')
            d_handler.output_handler.add_file(file_name=p.z_map_corrected_norm_png.format(d_handler.tag),           file_tag='z_map_corrected_normalised_png',   dir_tag='images')
            d_handler.output_handler.add_file(file_name=p.z_map_qq_plot_png.format(d_handler.tag),                  file_tag='z_map_qq_plot_png',                dir_tag='images')
            d_handler.output_handler.add_file(file_name=p.bdc_est_png.format(d_handler.tag, '{!s}'),                file_tag='bdc_est_png',                     dir_tag='images')
            d_handler.output_handler.add_file(file_name=p.unc_qqplot_png.format(d_handler.tag),                     file_tag='unc_qqplot_png',                   dir_tag='images')
            d_handler.output_handler.add_file(file_name=p.obs_qqplot_sorted_png.format(d_handler.tag),              file_tag='obs_qqplot_sorted_png',            dir_tag='images')
            d_handler.output_handler.add_file(file_name=p.obs_qqplot_unsorted_png.format(d_handler.tag),            file_tag='obs_qqplot_unsorted_png',          dir_tag='images')

            # Analysis files
            d_handler.output_handler.add_file(file_name=f.z_peaks_csv.format(d_handler.tag), file_tag='z_peaks_csv')

            # Scripts
            d_handler.output_handler.add_dir(dir_name='scripts', dir_tag='scripts', top_dir_tag='root')
            d_handler.output_handler.add_file(file_name=f.pymol_script,     file_tag='pymol_script',    dir_tag='scripts')
            d_handler.output_handler.add_file(file_name=f.ccp4mg_script,    file_tag='ccp4mg_script',   dir_tag='scripts')

            # Output blobs
            d_handler.output_handler.add_dir(dir_name='blobs', dir_tag='blobs', top_dir_tag='root')
            d_handler.output_handler.add_file(file_name=f.ccp4mg_png,       file_tag='ccp4mg_png',      dir_tag='blobs')

            # Pickled objects
            d_handler.output_handler.add_dir(dir_name='pickles', dir_tag='pickles', top_dir_tag='root')
            d_handler.output_handler.add_file(file_name=f.dataset_pickle,   file_tag='dataset_pickle',  dir_tag='pickles')

            ##############################################################################################################

            # Links for the dataset input files
            link_pdb = d_handler.output_handler.get_file('input_structure')
            link_mtz = d_handler.output_handler.get_file('input_data')
            # Link the input files to the output folder
            if not os.path.exists(link_pdb): rel_symlink(orig=d_handler.pdb_filename(), link=link_pdb)
            if not os.path.exists(link_mtz): rel_symlink(orig=d_handler.mtz_filename(), link=link_mtz)

            # Search for ligand files and link them to the output ligands folder
            lig_glob  = os.path.join(os.path.dirname(d_handler.pdb_filename()), lig_style)
            lig_files = glob.glob(lig_glob)
            for lig_file in lig_files:
                # Find all files with the same basename but allowing for different extensions. Then link to output folder.
                lig_base = os.path.splitext(lig_file)[0] + '.*'
                lig_matches = glob.glob(lig_base)
                for lig in lig_matches:
                    out_path = os.path.join(d_handler.output_handler.get_dir('ligand'), os.path.basename(lig))
                    if os.path.exists(lig) and (not os.path.exists(out_path)):
                        try: os.symlink(lig, out_path)
                        except: pass

            # Lastly: Update the pointer to the new path (relative to the pandda directory)
            d_handler._pdb_file = os.path.relpath(link_pdb, start=self.out_dir)
            d_handler._mtz_file = os.path.relpath(link_mtz, start=self.out_dir)

        self.datasets.add(loaded_datasets)
        self.log('{!s} Datasets Loaded (New).          '.format(len(loaded_datasets), True))
        self.log('{!s} Datasets Loaded (Total).        '.format(self.datasets.size(), True))

    #########################################################################################################
    #                                                                                                       #
    #                                     Dataset loading / processing                                      #
    #                                                                                                       #
    #########################################################################################################

    def load_reflection_data(self, ampl_label, phas_label):
        """Extract amplitudes and phases for creating map"""

        # Extract reflection data for the reference dataset
        ref_ampl_label = self.args.input.reference.amp_label if self.args.input.reference.amp_label else ampl_label
        ref_phas_label = self.args.input.reference.pha_label if self.args.input.reference.pha_label else phas_label
        self.reference_dataset().sfs = extract_structure_factors(self.reference_dataset().reflection_data(), ampl_label=ref_ampl_label, phas_label=ref_phas_label)

        t1 = time.time()
        self.log('----------------------------------->>>', True)
        for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True):
            # Extract miller array of structure factors
            if d_handler.sfs != None:
                print('\rAlready Loaded: Dataset {!s}          '.format(d_handler.tag), end=''); sys.stdout.flush()
            else:
                print('\rLoading Diffraction Data: Dataset {!s}                 '.format(d_handler.tag), end=''); sys.stdout.flush()
                d_handler.sfs = extract_structure_factors(mtz_object=d_handler.reflection_data(), ampl_label=ampl_label, phas_label=phas_label)
        t2 = time.time()
        self.log('\r> Structure Factors Extracted > Time Taken: {!s} seconds'.format(int(t2-t1)), True)

    def align_datasets(self, method):
        """Align each structure the reference structure"""

        assert method in ['local','global'], 'METHOD NOT DEFINED: {!s}'.format(method)

        # Select the datasets for alignment
        if method == 'local':
            # If local alignment has been chosen, also do a global alignment
            method = 'both'
            # Select the datasets that don't have local alignments
            align_datasets = [d_handler for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True) if not d_handler.local_alignment_transforms()]
        else:
            # Select the datasets that don't have global alignments
            align_datasets = [d_handler for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True) if not d_handler.global_alignment_transform()]

        if not align_datasets:
            self.log('----------------------------------->>>', True)
            self.log('No datasets to align')
            return

        # Align the datasets using multiple cores if possible
        self.log('----------------------------------->>>', True)
        print('Generating Alignments (using {!s} cores) for {} datasets'.format(self.settings.cpus, len(align_datasets)))
        start = time.time()
        arg_list = [{'r_handler':self.reference_dataset(), 'd_handler':d_handler, 'method':method} for d_handler in align_datasets]
        alignment_transforms = easy_mp.pool_map(func=align_map_func, args=arg_list, processes=self.settings.cpus)
        alignment_transforms = dict(alignment_transforms)
        finish = time.time()
        self.log('\r> Generating Alignments > Time Taken: {!s} seconds'.format(int(finish-start)), True)

        t1 = time.time()
        # Post-process the alignments (write out aligned structures etc) - only need to do it for "new" datasets
#        for d_handler in self.datasets.mask(mask_name='old datasets', invert=True):
        for d_handler in align_datasets:
            print('\rAligning Structures: Dataset {!s}          '.format(d_handler.tag), end=''); sys.stdout.flush()

            # Align to reference structure to get mapping transform
            #d_handler.align_to_reference(ref_handler=self.reference_dataset(), method=method)

            global_mx, local_mxs = alignment_transforms[d_handler.tag]

            # Store the alignments in the dataset handler
            d_handler.set_global_alignment(alignment=global_mx)
            d_handler.set_local_alignments(alignment=local_mxs)

            assert d_handler.global_alignment_transform()

            # Create a copy to transform
            aligned_struc = d_handler.hierarchy().deep_copy()

            # Always do a global alignment
            aligned_struc.atoms().set_xyz(d_handler.transform_to_reference(points=d_handler.hierarchy().atoms().extract_xyz(), method='global'))
            # Also write out into the aligned_structures folder
            aligned_struc.write_pdb_file(file_name=os.path.join(self.output_handler.get_dir('aligned_structures'), '{!s}-global-aligned.pdb'.format(d_handler.tag)))

            if method == 'both':
                # Write out the aligned structure from the local alignments
                aligned_struc.atoms().set_xyz(d_handler.transform_to_reference(points=d_handler.hierarchy().atoms().extract_xyz(), method='local'))
                # Also write out into the aligned_structures folder
                aligned_struc.write_pdb_file(file_name=os.path.join(self.output_handler.get_dir('aligned_structures'), '{!s}-local-aligned.pdb'.format(d_handler.tag)))

            # Write out next to the original files (if global, global, if local, local)
            aligned_struc.write_pdb_file(file_name=d_handler.output_handler.get_file('aligned_structure'))

            # Align to the reference to find the alignment rmsd
            mapped_ref_sites = d_handler.transform_from_reference(points=self.reference_dataset().hierarchy().atoms().extract_xyz(), method='global')
            alignment_rmsd = d_handler.hierarchy().atoms().extract_xyz().rms_difference(mapped_ref_sites)
            # Check to see if this is the reference dataset
            if alignment_rmsd == 0.0:
                # This must be the reference! Set the dataset number if it's not already set
                if self._ref_dataset_index is None:
                    self._ref_dataset_index = d_handler.num
                    print('REFERENCE FOUND! {!s}                          '.format(d_handler.tag))
                # Raise error if the reference has already been set and it's not this dataset
                elif self._ref_dataset_index != d_handler.num:
                    raise Exception('ALIGNED OBJECT EQUAL TO UNALIGNED OBJECT - THIS IS MOST UNLIKELY')

        t2 = time.time()
        self.log('\r> Aligning Structures > Time Taken: {!s} seconds'.format(int(t2-t1)), True)

    def collate_dataset_variables(self):
        """Go through all of the datasets and collect lots of different characteristics of the datasets for identifying odd datasets"""

        self.log('----------------------------------->>>', True)
        self.log('Collating Dataset Structure/Crystal Variables', True)

        for d in self.datasets.all():
            # Resolution info
            self.tables.dataset_info.set_value(d.tag, 'high_resolution', numpy.round(d.mtz_summary.high_res,3))
            self.tables.dataset_info.set_value(d.tag, 'low_resolution',  numpy.round(d.mtz_summary.low_res,3))
            # Unit cell info
            self.tables.dataset_info.set_value(d.tag, ['uc_a','uc_b','uc_c','uc_alpha','uc_beta','uc_gamma'],   numpy.round(d.mtz_summary.unit_cell.parameters(),3))
            self.tables.dataset_info.set_value(d.tag, 'uc_vol',                                                 numpy.round(d.mtz_summary.unit_cell.volume()),3)
            # Spacegroup info
            self.tables.dataset_info.set_value(d.tag, 'space_group', d.mtz_summary.space_group.info().type().lookup_symbol())
            # Quality info
            self.tables.dataset_info.set_value(d.tag, 'r_work', round_no_fail(d.input().get_r_rfree_sigma().r_work,3))
            self.tables.dataset_info.set_value(d.tag, 'r_free', round_no_fail(d.input().get_r_rfree_sigma().r_free,3))

    #########################################################################################################
    #                                                                                                       #
    #                                 Dataset variation analysis functions                                  #
    #                                                                                                       #
    #########################################################################################################

    def calculate_dataset_rmsds_to_reference(self):
        """Go through all of the datasets and collect lots of different characteristics of the datasets for identifying odd datasets"""

        self.log('----------------------------------->>>', True)
        self.log('Calculating Dataset CAlpha RMSDs to Reference', True)

        # Now calculate the variation in the structure, from the reference
        for d in self.datasets.mask(mask_name='rejected - total', invert=True):
            rmsd = d.calpha_sites().rms_difference(d.transform_from_reference(points=self.reference_dataset().calpha_sites(), method='global'))
            self.tables.dataset_info.set_value(d.tag, 'rmsd_to_reference', numpy.round(rmsd,3))

    def analyse_alignment_variation(self):
        """Look at all of the rotation matrices for the local alignments and calculate the rms between neighbours"""

        assert self.params.alignment.method in ['global', 'local']

        if self.params.alignment.method == 'global':
            self.log('GLOBAL ALIGNMENT SELECTED - NOT ANALYSING ROTATION MATRICES')
            return

        if self.settings.plot_graphs:
            import matplotlib
            matplotlib.interactive(0)
            from matplotlib import pyplot

        # Select datasets to analyse
        used_datasets = self.datasets.mask(mask_name='rejected - total', invert=True)

        # Reference c_alpha labels
        ref_c_alpha_labels = sorted(used_datasets[0].local_alignment_transforms().keys())

        # Array to hold the output data
        num_datasets = len(used_datasets)
        num_pairs =  len(ref_c_alpha_labels)-1
        output_diffs = numpy.zeros((num_datasets, num_pairs, 2))

        # Iterate through the datasets and pull out the alignment matrices
        for d_num, d_handler in enumerate(used_datasets):
            # Extract and sort dataset alignments
            alignments = d_handler.local_alignment_transforms()
            alignment_keys = sorted(alignments.keys())
            assert alignment_keys == ref_c_alpha_labels

            # Iterate through adjacent pairs of matrices
            for i in range(0, num_pairs):
                # Label and lsq fit for the current calpha
                calpha_1 = alignment_keys[i]
                rt_1 = alignments[calpha_1]
                # And for the next calpha
                calpha_2 = alignment_keys[i+1]
                rt_2 = alignments[calpha_2]

                assert calpha_1 == ref_c_alpha_labels[i]
                assert calpha_2 == ref_c_alpha_labels[i+1]

                # Calculate the difference in the angles of the alignment matrices
                theta_1 = scitbx.math.math.acos((rt_1.r.trace()-1)/2.0)
                theta_2 = scitbx.math.math.acos((rt_2.r.trace()-1)/2.0)
                # XXX Should we calculate the absolute of the difference?
                theta_rad = theta_2-theta_1
                theta_deg = theta_rad * 180.0/scitbx.math.math.pi
                # Calculate the difference in the translation
                t_shift = (rt_2.t-rt_1.t).norm_sq()**0.5

#                # Calculate the angles from the multiplication of one by the inverse of the other
#                rt_1_2 = rt_1 * rt_2.inverse()
#                # Calculate the angle of the rotation matrix
#                theta_rad = scitbx.math.math.acos((rt_1_2.r.trace()-1)/2.0)
#                theta_deg = theta_rad * 180.0/scitbx.math.math.pi
#                # Calculate the length of the shift
#                t_shift =  rt_1_2.t.norm_sq()**0.5

                # Append to the array
                output_diffs[d_num, i, :] = theta_deg, t_shift

        # Directory to write the output to
        var_out_dir = self.output_handler.get_dir('analyses')
        # Write out to file
        numpy.savetxt(  fname = os.path.join(var_out_dir, 'calpha_rt_r_variation.csv'), X=output_diffs[:,:,0], delimiter=',', newline='\n' )
        numpy.savetxt(  fname = os.path.join(var_out_dir, 'calpha_rt_t_variation.csv'), X=output_diffs[:,:,1], delimiter=',', newline='\n' )

        # Write out graphs
        if self.settings.plot_graphs:

            # Create labels
            labels = ['']*num_pairs
            for i in range(0, num_pairs, 5)+[num_pairs-1]:
                labels[i] = i+1
            # Clear the last n before the last one
            n = 4
            labels[-1-n:-1] = ['']*n

            # BOX PLOT OF ROTATION AND TRANSLATION SHIFTS
            fig = pyplot.figure()
            pyplot.title('Rotation-translation alignment matrix variation between adjacent C-alpha')
            # ADJACENT ANGLE VARIATION
            pyplot.subplot(2, 1, 1)
            pyplot.boxplot(x=output_diffs[:,:,0], notch=True, sym='.', widths=0.5, whis=[5,95], whiskerprops={'ls':'-'}, flierprops={'ms':1}, labels=labels) # whis='range'
            pyplot.xlabel('C-alpha index')
            pyplot.ylabel('Angle Difference\n(degrees)')
            # ADJACENT SHIFT VARIATION
            pyplot.subplot(2, 1, 2)
            pyplot.boxplot(x=output_diffs[:,:,1], notch=True, sym='.', widths=0.5, whis=[5,95], whiskerprops={'ls':'-'}, flierprops={'ms':1}, labels=labels) # whis='range'
            pyplot.xlabel('C-alpha Index')
            pyplot.ylabel('Translation Difference\n(angstroms)')
            # Apply tight layout to prevent overlaps
            pyplot.tight_layout()
            # Save both
            pyplot.savefig(os.path.join(var_out_dir, 'calpha_rt_variation.png'), format='png')
            pyplot.close(fig)

    #########################################################################################################
    #                                                                                                       #
    #                                             Dataset filtering                                         #
    #                                                                                                       #
    #########################################################################################################

    def filter_datasets_1(self, filter_dataset=None):
        """Filter out the datasets which contain different protein models (i.e. protein length, sequence, etc)"""

        self.log('----------------------------------->>>', True)
        self.log('Filtering Datasets (Datasets that are different to the reference dataset). Potential Classes:', True)
        for failure_class in PanddaMaskNames.structure_mask_names:
            self.log('\t{!s}'.format(failure_class), True)
        self.log('---------->>>', True)

        # If no filtering dataset given, filter against the reference dataset
        if not filter_dataset: filter_dataset = self.reference_dataset()

        # Check that the same protein structure is present in each dataset - THIS MASK SHOULD INCLUDE ALL DATASETS AT FIRST
        for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True):

            print('\rFiltering Dataset {!s}          '.format(d_handler.tag), end=''); sys.stdout.flush()
            # Check the space group of the dataset
            if d_handler.input().crystal_symmetry().space_group().info().symbol_and_number() != filter_dataset.input().crystal_symmetry().space_group().info().symbol_and_number():
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.tag))
                self.log('Different Space Group')
                self.log('Reference: {!s}, {!s}: {!s}'.format(filter_dataset.input().crystal_symmetry().space_group().info().symbol_and_number(),
                                                        d_handler.tag, d_handler.input().crystal_symmetry().space_group().info().symbol_and_number()))
                self.log('---------->>>', True)
                self.tables.dataset_info.set_value(d_handler.tag, 'rejection_reason', 'Different Space Group to Reference')
                self.datasets.all_masks().set_mask_value(mask_name='bad structure - different space group', entry_id=d_handler.tag, value=True)
            # Check that the hierarchies are identical
            if not d_handler.hierarchy().is_similar_hierarchy(filter_dataset.hierarchy()):
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.tag))
                self.log('Non-Identical Structure (Structures do not contain the same atoms)')
                self.log('---------->>>', True)
                self.tables.dataset_info.set_value(d_handler.tag, 'rejection_reason', 'Atoms present in the dataset are different to atoms present in the reference structure')
                self.datasets.all_masks().set_mask_value(mask_name='bad structure - non-identical structures', entry_id=d_handler.tag, value=True)

        # Combine structure_masks
        structure_reject_mask = self.datasets.all_masks().combine_masks(PanddaMaskNames.structure_mask_names)
        self.datasets.all_masks().add_mask(mask_name='rejected - structure', mask=structure_reject_mask)

        # Update the combined masks
        combined_reject_mask = self.datasets.all_masks().combine_masks(PanddaMaskNames.reject_mask_names)
        self.datasets.all_masks().add_mask(mask_name='rejected - total', mask=combined_reject_mask)

        self.log('\rDatasets Filtered.                          ', True)
        self.log('----------------------------------->>>')
        self.log('Rejected Datasets (Structure): {!s}'.format(sum(self.datasets.all_masks().get_mask(mask_name='rejected - structure'))), True)
        self.log('Rejected Datasets (Total):     {!s}'.format(sum(self.datasets.all_masks().get_mask(mask_name='rejected - total'))), True)
        self.log('----------------------------------->>>')

        reject_reasons = self.tables.dataset_info['rejection_reason'].value_counts().sort_index()
        if reject_reasons.any():
            self.log('Reasons for Rejection:')
            for reason, count in reject_reasons.iteritems():
                self.log('{} Dataset(s) - {}'.format(count, reason))

        # Link all rejected datasets into the rejected directory
        for d_handler in self.datasets.mask(mask_name='rejected - total'):
            reject_dir = os.path.join(self.output_handler.get_dir('rejected_datasets'), d_handler.tag)
            if not os.path.exists(reject_dir):
                rel_symlink(orig=d_handler.output_handler.get_dir('root'), link=reject_dir)

    def filter_datasets_2(self):
        """Filter out the non-isomorphous datasets"""

        self.log('----------------------------------->>>', True)
        self.log('Filtering Datasets (bad-quality or large rmsd structures). Potential Classes:', True)
        for failure_class in PanddaMaskNames.crystal_mask_names:
            self.log('\t{!s}'.format(failure_class), True)
        self.log('---------->>>', True)

        # Check that each dataset is similar enough to be compared
        for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True):

            print('\rFiltering Dataset {!s}          '.format(d_handler.tag), end=''); sys.stdout.flush()
            # Check that it correlates well with itself before and after scaling
            if d_handler.input().get_r_rfree_sigma().r_free > self.params.filtering.max_rfree:
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.tag))
                self.log('RFree is higher than cutoff: {!s}'.format(self.params.filtering.max_rfree))
                self.log('High RFree: {!s}'.format(d_handler.input().get_r_rfree_sigma().r_free))
                self.log('---------->>>', True)
                self.tables.dataset_info.set_value(d_handler.tag, 'rejection_reason', 'R-free is too high')
                self.datasets.all_masks().set_mask_value(mask_name='bad crystal - rfree', entry_id=d_handler.tag, value=True)
            # Check the deviation from the average sites
            elif d_handler.calpha_sites().rms_difference(d_handler.transform_from_reference(points=self.reference_dataset().calpha_sites(), method='global')) > self.params.filtering.max_rmsd_to_reference:
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.tag))
                self.log('C-alpha RMSD is too large')
                self.log('Aligned (Calpha) RMSD: {!s}'.format(d_handler.calpha_sites().rms_difference(d_handler.transform_from_reference(points=self.reference_dataset().calpha_sites(), method='global'))))
                self.log('---------->>>', True)
                self.tables.dataset_info.set_value(d_handler.tag, 'rejection_reason', 'High RMSD to aligned reference structure')
                self.datasets.all_masks().set_mask_value(mask_name='bad crystal - isomorphous structure', entry_id=d_handler.tag, value=True)
            else:
                pass

        # Combine crystal masks
        crystal_reject_mask = self.datasets.all_masks().combine_masks(PanddaMaskNames.crystal_mask_names)
        self.datasets.all_masks().add_mask(mask_name='rejected - crystal', mask=crystal_reject_mask)

        # Combine all of the masks
        combined_reject_mask = self.datasets.all_masks().combine_masks(PanddaMaskNames.reject_mask_names)
        self.datasets.all_masks().add_mask(mask_name='rejected - total', mask=combined_reject_mask)

        self.log('\rDatasets Filtered.               ', True)
        self.log('----------------------------------->>>')
        self.log('Rejected Datasets (Structure): {!s}'.format(sum(self.datasets.all_masks().get_mask(mask_name='rejected - structure'))), True)
        self.log('Rejected Datasets (Crystal):   {!s}'.format(sum(self.datasets.all_masks().get_mask(mask_name='rejected - crystal'))), True)
        self.log('Rejected Datasets (Total):     {!s}'.format(sum(self.datasets.all_masks().get_mask(mask_name='rejected - total'))), True)

        reject_reasons = self.tables.dataset_info['rejection_reason'].value_counts().sort_index()
        if reject_reasons.any():
            self.log('Reasons for Rejection:')
            for reason, count in reject_reasons.iteritems():
                self.log('{} Dataset(s) - {}'.format(count, reason))

        # Link all rejected datasets into the rejected directory
        for d_handler in self.datasets.mask(mask_name='rejected - total'):
            reject_dir = os.path.join(self.output_handler.get_dir('rejected_datasets'), d_handler.tag)
            if not os.path.exists(reject_dir):
                rel_symlink(orig=d_handler.output_handler.get_dir('root'), link=reject_dir)

    #########################################################################################################
    #                                                                                                       #
    #                              Dataset utility functions (e.g. reset/sync)                              #
    #                                                                                                       #
    #########################################################################################################

    def reset_loaded_datasets(self):
        """Check that pickled datasets are ready for reprocessing, etc, if required"""

        if self.args.method.reprocess_selected_datasets: datasets_for_reprocessing = self.args.method.reprocess_selected_datasets.split(',')
        else:                                            datasets_for_reprocessing = []

        for d_handler in self.datasets.all():
            if self.args.method.reprocess_existing_datasets or (d_handler.tag in datasets_for_reprocessing):

                # Reset the meta objects
                d_handler.meta.analysed = False
                d_handler.meta.dataset_map_info = None

                # Delete events from before
                d_handler.events = []

                # Reset the map information for the dataset
                self.tables.dataset_map_info.loc[d_handler.tag] = numpy.nan

    def check_loaded_datasets(self, datasets):
        """Check that the datasets are analysable (have the right mtz columns, etc)"""

        self.log('----------------------------------->>>', True)
        self.log('Performing checks on the loaded datasets', True)

        for d_handler in datasets:
            if not d_handler.reflection_data().has_column(self.params.maps.ampl_label):
                raise Sorry('Amplitude column {} was not found in the reflection data for dataset {}. You may need to change the pandda.params.maps.ampl_label option.'.format(self.params.maps.ampl_label, d_handler.tag))
            if not d_handler.reflection_data().has_column(self.params.maps.phas_label):
                raise Sorry('Phase column {} was not found in the reflection data for dataset {}. You may need to change the pandda.params.maps.phas_label option.'.format(self.params.maps.phas_label, d_handler.tag))

    def sync_datasets(self, datasets=None, overwrite_dataset_meta=False):
        """Sync the loaded datasets and the pandda dataset tables"""

        if not datasets: datasets = self.datasets.all()

        for d_handler in datasets:
            # Copy data from pandda dataset tables to dataset
            if (d_handler.meta.dataset_info is None) or overwrite_dataset_meta:
                d_handler.meta.dataset_info = self.tables.dataset_info.loc[d_handler.tag]
            # Copy data from dataset to pandda dataset tables
            else:
                for col,val in d_handler.meta.dataset_info.iteritems():
                    self.tables.dataset_info.set_value(index=d_handler.tag, col=col, value=val)

            # Copy data from pandda dataset tables to dataset
            if (d_handler.meta.dataset_map_info is None) or overwrite_dataset_meta:
                d_handler.meta.dataset_map_info = self.tables.dataset_map_info.loc[d_handler.tag]
            # Copy data from dataset to pandda dataset tables
            else:
                for col,val in d_handler.meta.dataset_map_info.iteritems():
                    self.tables.dataset_map_info.set_value(index=d_handler.tag, col=col, value=val)

    #########################################################################################################
    #                                                                                                       #
    #                                Analysis functions (Dataset selection)                                 #
    #                                                                                                       #
    #########################################################################################################

    def select_for_building_distributions(self, high_res_cutoff, building_mask_name):
        """Select all datasets with resolution better than high_res_cutoff"""

        # Create empty mask
        self.datasets.all_masks().add_mask(mask_name=building_mask_name, mask=[False]*self.datasets.size())
        # Counter for the number of datasets to select
        total_build = 0
        # Select from the datasets that haven't been rejected
        for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True):
            # Check the resolution of the dataset
            if self.tables.dataset_info.get_value(index=d_handler.tag, col='high_resolution') > high_res_cutoff:
                continue
            # Check to see if this has been excluded from building
            elif self.datasets.all_masks().get_mask_value(mask_name='no_build', entry_id=d_handler.tag) == True:
                self.log('Rejecting Dataset {!s}: Excluded from building'.format(d_handler.tag))
                continue
            else:
                self.datasets.all_masks().set_mask_value(mask_name=building_mask_name, entry_id=d_handler.tag, value=True)
                # Check to see if the number of datasets to use in building has been reached
                total_build += 1
                if total_build >= self.params.analysis.max_build_datasets:
                    self.log('Maximum number of datasets for building reached: {!s}={!s}'.format(total_build, self.params.analysis.max_build_datasets))
                    break
        return self.datasets.all_masks().get_mask(building_mask_name)

    def select_for_analysis(self, high_res_large_cutoff, high_res_small_cutoff, analysis_mask_name):
        """Select all datasets with resolution between high and low limits"""

        assert high_res_large_cutoff > high_res_small_cutoff, '{!s} must be larger than {!s}'.format(high_res_large_cutoff, high_res_small_cutoff)

        if self.args.method.reprocess_selected_datasets: datasets_for_reprocessing = self.args.method.reprocess_selected_datasets.split(',')
        else:                                            datasets_for_reprocessing = []

        # Create empty mask
        self.datasets.all_masks().add_mask(mask_name=analysis_mask_name, mask=[False]*self.datasets.size())
        # Select from the datasets that haven't been rejected
        for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True):
            # Check the resolution of the dataset (is not too low)
            if self.tables.dataset_info.get_value(index=d_handler.tag, col='high_resolution') > high_res_large_cutoff:
                continue
            # Check the resolution of the dataset (is not too high)
            elif self.tables.dataset_info.get_value(index=d_handler.tag, col='high_resolution') <= high_res_small_cutoff:
                continue
            # Check to see if this has been excluded from building
            elif self.datasets.all_masks().get_mask_value(mask_name='no_analyse', entry_id=d_handler.tag) == True:
                self.log('Rejecting Dataset {!s}: Excluded from analysis'.format(d_handler.tag))
                continue
            elif self.datasets.all_masks().get_mask_value(mask_name='old datasets', entry_id=d_handler.tag) and (not self.args.method.reprocess_existing_datasets) and (d_handler.tag not in datasets_for_reprocessing):
                self.log('Rejecting Dataset {!s}: Already Processed (Old Dataset)'.format(d_handler.tag))
                continue
            else:
                self.datasets.all_masks().set_mask_value(mask_name=analysis_mask_name, entry_id=d_handler.tag, value=True)
        return self.datasets.all_masks().get_mask(analysis_mask_name)

    def truncate_scaled_data(self, dataset_handlers, truncation_stuff=None):
        """Truncate data at the same indices across all the datasets"""

        self.log('----------------------------------->>>', True)
        self.log('Truncating Reflection Data', True)

        # Calculate which reflections are present in the reference dataset
        ref_size = self.reference_dataset().sfs.set().size()
        self.log('Number of Reflections in Reference Dataset: {!s}'.format(ref_size))

        # Truncate reflections to the common set (not including the reference dataset)
        common_set = dataset_handlers[0].sfs.set()
        for i_dh, d_handler in enumerate(dataset_handlers):
            common_set = common_set.common_set(d_handler.sfs, assert_is_similar_symmetry=False)

        self.log('----------------------------------->>>', True)
        self.log('Number of Common Reflections between Datasets: {!s} ({!s}% of reference)'.format(common_set.size(), int(100.0*common_set.size()/ref_size)))
        self.log('After Truncation - Reflections per dataset: {!s}'.format(common_set.size()))

        # Create maps for all of the datasets (including the reference dataset)
        for d_handler in [self.reference_dataset()]+dataset_handlers:
            d_handler.tr_sfs = d_handler.sfs.common_set(common_set, assert_is_similar_symmetry=False)

        reslns = [d.tr_sfs.d_min() for d in dataset_handlers]
        min_res = min(reslns)
        max_res = max(reslns)

        self.log('After Truncation - Resolution Range: {!s}-{!s}'.format(min_res, max_res))

        # TODO WRITE HISTOGRAMS OF RESOLUTIONS?!

    def load_reference_map(self, map_resolution=0):
        """Load the reference map, and calculate some map statistics"""

        # Get the reference dataset
        ref_handler = self.reference_dataset()

        # Take the scaled diffraction data for the dataset and create fft
        fft_map = ref_handler.tr_sfs.fft_map( resolution_factor = self.params.maps.resolution_factor,
                                              d_min = map_resolution,
                                              symmetry_flags = cctbx.maptbx.use_space_group_symmetry )
        # Scale the map
        if self.params.maps.scaling == 'none':     pass
        elif self.params.maps.scaling == 'sigma':  fft_map.apply_sigma_scaling()
        elif self.params.maps.scaling == 'volume': fft_map.apply_volume_scaling()
        # Store the unit cell and the space group for later
        ref_handler.unit_cell   = fft_map.unit_cell()
        ref_handler.space_group = fft_map.space_group()
        # Create map handler and map to the reference frame
        native_map_handler = map_handler(   map_data  = fft_map.real_map(),
                                            unit_cell = fft_map.unit_cell()   )
        # Extract the points for the morphed maps (in the reference frame)
        masked_gps      = list(self.reference_grid().global_mask().outer_mask())
        masked_idxs     = self.reference_grid().global_mask().outer_mask_indices()
        masked_cart_ref = flex.vec3_double(masked_gps)*self.reference_grid().grid_spacing()
        # Transform Coordinates
        masked_cart_nat = ref_handler.transform_from_reference( points = masked_cart_ref,
                                                                method = 'global' )
        # Sample the map at the masked points
        masked_vals_ref = native_map_handler.get_cart_values(masked_cart_nat)
        # Record the Mean and RMS values of these map values so that the other datasets can be scaled to the same
        map_mean = masked_vals_ref.min_max_mean().mean
        map_rms  = masked_vals_ref.standard_deviation_of_the_sample()
        # Calculate the sparse vector of the masked map values
        masked_map_sparse = scitbx.sparse.vector(self.reference_grid().grid_size_1d(), dict(zip(masked_idxs, masked_vals_ref)))
        masked_map        = masked_map_sparse.as_dense_vector()
        # Map Holder for the calculated map
        ref_map_holder = MapHolder( num         = ref_handler.num,
                                    tag         = ref_handler.tag,
                                    map         = masked_map,
                                    # Change these for the 'fake' grid unit_cell and a P1 space_group
                                    unit_cell   = fft_map.unit_cell(),
                                    space_group = fft_map.space_group(),
                                    meta        = Meta({'type'            : 'reference observed map',
                                                        'resolution'      : map_resolution,
#                                                        'masked_map_vals' : masked_vals_ref,
                                                        'map_mean'        : map_mean,
                                                        'map_rms'         : map_rms   }),
                                    parent      = ref_handler    )

        return ref_map_holder

    def load_and_morph_maps(self, dataset_handlers, ref_map_holder, map_resolution=0):
        """Create map from miller arrays. Transform map into the reference frame by sampling at the given points."""

        self.log('----------------------------------->>>', True)
        self.log('Loading Electron Density Maps @ {!s}A'.format(map_resolution), True)

        # Create holder for the output map objects
        map_holder_list = MapHolderList()

        # Extract the points for the morphed maps (in the reference frame)
        masked_gps      = list(self.reference_grid().global_mask().outer_mask())
        masked_idxs     = self.reference_grid().global_mask().outer_mask_indices()
        masked_cart_ref = flex.vec3_double(masked_gps)*self.reference_grid().grid_spacing()
        # Mapping of grid points to rotation matrix keys (residue CA labels)
        if self.params.alignment.method == 'local':
            masked_cart_map = self.reference_grid().partition().query_by_grid_indices(masked_idxs)
        else:
            masked_cart_map = None

        self.log('----------------------------------->>>', True)
        self.log('FFT-ing Maps')
        start = time.time()
        # FFT the maps from the truncated data
        for d_handler in dataset_handlers:
            # Take the truncated diffraction data for each dataset and create fft
            d_handler.fft_map = d_handler.tr_sfs.fft_map( resolution_factor = self.params.maps.resolution_factor,
                                                          d_min             = map_resolution,
                                                          symmetry_flags    = cctbx.maptbx.use_space_group_symmetry  )
        finish = time.time()
        self.log('\r> FFT-ing Maps ({!s} Datasets) > Time Taken: {!s} seconds'.format(len(dataset_handlers), int(finish-start)), True)

        # Load maps using multiple cores if possible
        self.log('----------------------------------->>>', True)
        print('Loading Maps (using {!s} cores)'.format(self.settings.cpus))
        start = time.time()
        arg_list = [{ 'd_handler' : d_handler,    'params' : self.params,     'map_resolution' : map_resolution, \
                      'grid'      : self.reference_grid().grid_indexer(),     'ref_map_holder' : ref_map_holder, \
                      'masked_cart_ref' : masked_cart_ref, \
                      'masked_cart_map' : masked_cart_map, \
                      'masked_idxs'     : masked_idxs, \
#                      'scaling_mask'    : self.reference_grid().global_mask().outer_mask_indices()      } for d_handler in dataset_handlers]
                      'scaling_mask'    : self.reference_grid().global_mask().inner_mask_indices()      } for d_handler in dataset_handlers]
        map_holders = easy_mp.pool_map(func=load_maps_map_func, args=arg_list, processes=self.settings.cpus, chunksize=1)
        # Append to the map holder list
        map_holder_list.add(map_holders)
        # Go through and assign parents
        for mh in map_holder_list.all():
            mh.parent = self.datasets.get(tag=mh.tag)
            # Clear the fft'd map
            mh.parent.fft_map = None

        finish = time.time()
        self.log('\r> Loading Maps ({!s} Datasets) > Time Taken: {!s} seconds'.format(map_holder_list.size(), int(finish-start)), True)

        # Set the mask lengths and entry ids
        map_holder_list.all_masks().set_mask_length(mask_length=map_holder_list.size())
        map_holder_list.all_masks().set_entry_ids(entry_ids=[mh.tag for mh in map_holder_list.all()])

        return map_holder_list

#    def process_z_map(self, z_map, local_mask_function=None):
#        """Process the z_maps, looking for groups of high z-values"""
#
#        # TODO DO IN PARALLEL!!!
#
#        # Translate between 3d grid point and 1d array index
#        grid_indexer = self.reference_grid().grid_indexer()
#
#        # Set local_mask function if not given - set to rms
#        if local_mask_function is None:
#            # Given arguments: (map values, grid points (relative to central point))
#            local_mask_function = lambda vals, gps: numpy.sqrt(numpy.mean(numpy.power(vals,2)))
#
#        # Initialise empty sparse vector
#        mod_z_map = sparse.vector(self.reference_grid().grid_size_1d())
#
#        self.log('Post Processing Z-Map.', True)
#
#        # Iterate through all grid points and calculate modified z-value at each point (0 at unsampled points or in buffer zone)
#        for p_num, gp in enumerate(self.reference_grid().masked_grid_points()):
#
#            status_bar(n=p_num, n_max=len(self.reference_grid().masked_grid_points()))
#
#            # Find nearby values within local mask
#            gp_masked = self.reference_grid().local_mask().apply_mask(gp)
#            # Find the map values at these points
#            local_map_vals = [z_map[grid_indexer(gp_m)] for gp_m in gp_masked]
#
#            # Use the local mask function to modify the extracted map values (use local mask vectors so relative to gp)
#            mod_z_val = local_mask_function(vals=local_map_vals, gps=self.reference_grid().local_mask().mask())
#
#            # Add to new map array
#            mod_z_map[grid_indexer(gp)] = mod_z_val
#
#        # Get down-sampled map
#        resamp_mod_z_map = [mod_z_map[grid_indexer(p)] for p in self.reference_grid().resampled_grid_points()]
#        assert flex.product(flex.int(self.reference_grid().resampled_grid_size())) == len(resamp_mod_z_map)
#
#        return mod_z_map.as_dense_vector(), flex.double(resamp_mod_z_map)

    def collate_event_counts(self):
        """Collate eventss from all of the datasets"""

        self.log('----------------------------------->>>', False)
        self.log('Collating Clusters', False)

        # List of points to be returned
        all_dataset_events = dict([(d.tag, d.events) for d in self.datasets.all()])

        # Print Cluster Summaries
        event_num = [(k, len(all_dataset_events[k])) for k in sorted(all_dataset_events.keys()) if all_dataset_events[k]]
        event_total = sum([a[1] for a in event_num])

        return event_total, event_num, all_dataset_events

    def cluster_events_and_update(self, events=[], update_tables=True, update_output=True):
        """Cluster events to sites and add information to the pandda tables"""

        if not events:
            print('No Events Found')
            return None

        self.log('----------------------------------->>>', True)
        self.log('Clustering identified events: {} Event(s)'.format(len(events)))
        # Cluster events to sites
        site_list = cluster_events(events=events, cutoff=15.0/self.reference_grid().grid_spacing(), linkage='average')
        # Sort sites by largest Z-Value
        site_list.sort(key=lambda s: (s.info.num_events, max([e.cluster.max for e in s.children])), reverse=True).renumber()
        # Add meta to the site list TODO implement this function -- blank at the moment TODO
        [s.find_protein_context(hierarchy=self.reference_dataset().hierarchy()) for s in site_list.children]
        # Update the pandda tables?
        if update_tables:
            self.update_site_table(site_list=site_list, clear_table=True)
            self.update_event_table_site_info(events=events)
        # Generate output images and graphs?
        if update_output:
            # Plot output graph of site list
            self.log('Deleting old images: ')
            delete_with_glob(glob_str=self.output_handler.get_file('analyse_site_graph_mult').format('*'))
            bar.multiple_bar_plot_over_several_images(
                                    f_template = self.output_handler.get_file('analyse_site_graph_mult'),
                                    plot_vals  = [sorted([e.cluster.max for e in s.children],reverse=True) for s in site_list.children]   )
            # Create pictures of the sites on the protein
            self.make_pymol_site_image_and_scripts(site_list=site_list, make_images=True)

        return site_list

#    def image_blob(self, script, image, d_handler, point, point_no, towards=[10,10,10]):
#        """Take pictures of the maps with ccp4mg"""
#
#        from giant.graphics import calculate_view_quaternion, multiply_quaternions
#
#        # Get the template to be filled in
#        template = PANDDA_HTML_ENV.get_template('ccp4mg-pic.py')
#
#        orientation = calculate_view_quaternion(towards, point)
#        rotate_1 = multiply_quaternions(orientation, (0.0, 0.5**0.5, 0.0, 0.5**0.5))
#        rotate_2 = multiply_quaternions(orientation, (0.5**0.5, 0.0, 0.0, 0.5**0.5))
#
#        for view_no, view in enumerate([orientation, rotate_1, rotate_2]):
#
#            view_script = script.format(point_no, view_no)
#            view_image  = image.format(point_no, view_no)
#
#            ccp4mg_script = template.render({
#                                                'view'  :{
#                                                                'camera_centre' : [-1*c for c in point],
#                                                                'orientation'   : list(view)
#                                                            },
#                                                'mol'   :{
#                                                                'path'  : d_handler.output_handler.get_file('aligned_structure'),
#                                                                'name'  : 'aligned_structure'
#                                                            },
#                                         #       'map'   :{
#                                         #                       'path'    : d_handler.output_handler.get_file('sampled_map'),
#                                         #                       'name'    : 'sampled_map',
#                                         #                       'contour' : [1]
#                                         #                   },
#                                                'diff_map' :{
#                                                                'path'    : d_handler.output_handler.get_file('z_map_corrected_normalised'),
#                                                                'name'    : 'diff_map',
#                                                            #    'neg-contour' : -3,
#                                                                'pos-contour' : [2,3,4,5]
#                                                            }
#                                            })
#
#            # Write out the ccp4mg script to the dataset's scripts folder
#            with open(view_script, 'w') as fh:
#                fh.write(ccp4mg_script)
#
#            # Make the images
#            c = CommandManager('ccp4mg')
#            c.SetArguments(['-norestore','-picture', view_script, '-R', view_image, '-RO', """'{"size":"1500x1500"}'""", '-quit'])
#            c.Run()
#
#            if not os.path.exists(view_image):
#                print('FAILED TO MAKE IMAGES')
#                print(c.err)

    def write_map_analyser_maps(self, map_analyser, analysis_mask_name):
        """Write statistical maps for a map_analyser object"""

        ########################################################

        map_res = map_analyser.meta.resolution

        self.log('----------------------------------->>>')
        self.log('=> Writing Statistical Maps of Analysis @ {!s}A'.format(map_res))

        ########################################################

        self.log('=> Writing Statistical Maps')

        write_array_to_map(output_file = self.output_handler.get_file('mean_map').format(map_res),
                           map_data    = map_analyser.statistical_maps.mean_map,
                           grid        = self.reference_grid()     )
        write_array_to_map(output_file = self.output_handler.get_file('stds_map').format(map_res),
                           map_data    = map_analyser.statistical_maps.stds_map,
                           grid        = self.reference_grid()     )
        write_array_to_map(output_file = self.output_handler.get_file('sadj_map').format(map_res),
                           map_data    = map_analyser.statistical_maps.sadj_map,
                           grid        = self.reference_grid()     )
        write_array_to_map(output_file = self.output_handler.get_file('skew_map').format(map_res),
                           map_data    = map_analyser.statistical_maps.skew_map,
                           grid        = self.reference_grid()     )
        write_array_to_map(output_file = self.output_handler.get_file('kurt_map').format(map_res),
                           map_data    = map_analyser.statistical_maps.kurt_map,
                           grid        = self.reference_grid()     )
        write_array_to_map(output_file = self.output_handler.get_file('bimo_map').format(map_res),
                           map_data    = map_analyser.statistical_maps.bimo_map,
                           grid        = self.reference_grid()     )

    def write_output_csvs(self):
        """Write CSV file of dataset variables"""

        self.log('----------------------------------->>>')
        self.log('Writing Dataset + Dataset Map Summary CSV')

        # Write the dataset information to csv file
        self.tables.dataset_info.to_csv(path_or_buf=self.output_handler.get_file('dataset_info'), index_label='dtag')
        self.tables.dataset_map_info.to_csv(path_or_buf=self.output_handler.get_file('dataset_map_info'), index_label='dtag')

        self.log('----------------------------------->>>')
        self.log('Writing COMBINED Dataset Summary CSV')

        # Join the tables on the index of the main table
        comb_tab = self.tables.dataset_info.join(self.tables.dataset_map_info, how='outer')
        comb_tab.to_csv(path_or_buf=self.output_handler.get_file('dataset_combined_info'), index_label='dtag')

        self.log('----------------------------------->>>')
        self.log('Writing Event+Site Summary CSVs')

        # Sort the event data by z-peak and write out
        sort_eve = self.tables.event_info.sort(columns=['site_idx',self.args.results.events.order_by], ascending=[1,0])
        sort_eve = sort_eve.join(comb_tab, how='right')
        sort_eve.to_csv(path_or_buf=self.output_handler.get_file('event_info'))
        # Sort the sites by number of events and write out
        sort_sit = self.tables.site_info.sort(columns=[self.args.results.sites.order_by],ascending=[0])
        sort_sit.to_csv( path_or_buf=self.output_handler.get_file('site_info'))

    def update_site_table(self, site_list, clear_table=True):
        """Add site entries to the site table"""

        # Clear an existing table
        if clear_table:
            self.tables.site_info = pandas.DataFrame(data    = None,
                                                     index   = self.tables.site_info.index.reindex([])[0],
                                                     columns = self.tables.site_info.columns)
        # Go through and update the site information
        for site in site_list.children:
            self.tables.site_info.loc[site.id,:] = None
            self.tables.site_info.set_value(site.id, 'centroid', tuple(flex.double(site.info.centroid)*self.reference_grid().grid_spacing()))
            self.tables.site_info.set_value(site.id, 'native_centroid', tuple(self.reference_dataset().transform_from_reference(
                                                                                                        points=flex.vec3_double([site.info.centroid])*self.reference_grid().grid_spacing(),
                                                                                                        method='global')[0]))

    def make_pymol_site_image_and_scripts(self, site_list, make_images=True):
        """Generate pymol script to mark the location of identified sites"""

        pymol_str =  '# Mark the identified sites on the protein\n'
        pymol_str += 'from pymol import cmd\n'
        pymol_str += 'from pymol.cgo import *\n'
        pymol_str += 'cmd.load("{}", "reference")\n'.format(os.path.relpath(self.output_handler.get_file('reference_on_origin'), start=self.output_handler.get_dir('output_summaries')))
        pymol_str += 'cmd.show_as("cartoon", "reference")\n'
        pymol_str += 'cmd.color("cyan", "reference")\n'
        # Add sphere at each of the sites
        for site in site_list.children:
            # Only print the site if it has more than one event
            if len(site.children) > 1:
                lab = 'site_{}'.format(site.id)
                com = tuple(flex.double(site.info.centroid)*self.reference_grid().grid_spacing())
                pymol_str += 'cmd.pseudoatom("{}", pos={}, vdw=2.5)\n'.format(lab, com)
                pymol_str += 'cmd.show("sphere", "{}")\n'.format(lab)
                pymol_str += 'cmd.label("{}", "{}")\n'.format(lab, site.id)
                pymol_str += 'cmd.color("deepteal", "{}")\n'.format(lab)
                pymol_str += 'cmd.set("label_color", "white", "{}")\n'.format(lab)
            # Label events as smaller spheres
            for event in site.children:
                lab = 'event'
                com = tuple(flex.double(event.cluster.centroid)*self.reference_grid().grid_spacing())
                pymol_str += 'cmd.pseudoatom("{}", pos={}, vdw=0.5)\n'.format(lab, com)
                pymol_str += 'cmd.show("sphere", "{}")\n'.format(lab)
                pymol_str += 'cmd.color("blue", "{}")\n'.format(lab)
        # Set label things...
        pymol_str += 'cmd.set("label_size", 25)\n'
        pymol_str += 'cmd.set("label_position", (0,0,4))\n'
        pymol_str += 'cmd.bg_color(color="white")\n'
        # Write as python script
        with open(self.output_handler.get_file(file_tag='pymol_sites_py'), 'w') as fh:
            fh.write(pymol_str)

        # Run Pymol to generate images and output to pngs
        if make_images:
            pymol_str =  '# Load the protein representation and output images of sites\n'
            pymol_str += 'run {}\n'.format(os.path.relpath(self.output_handler.get_file(file_tag='pymol_sites_py'), start=self.output_handler.get_dir('output_summaries')))
            pymol_str += 'set ray_opaque_background, off\n'
            pymol_str += 'set specular, off\n'
            pymol_str += 'orient\n'
            pymol_str += 'png {}, width=1200, dpi=300, ray=1\n'.format(os.path.relpath(self.output_handler.get_file(file_tag='pymol_sites_png_1'), start=self.output_handler.get_dir('output_summaries')))
            pymol_str += 'rotate y, 180\n'
            pymol_str += 'png {}, width=1200, dpi=300, ray=1\n'.format(os.path.relpath(self.output_handler.get_file(file_tag='pymol_sites_png_2'), start=self.output_handler.get_dir('output_summaries')))
            pymol_str += 'quit'

            with open(self.output_handler.get_file(file_tag='pymol_sites_pml'), 'w') as fh:
                fh.write(pymol_str)

            # Change into directory as script runs off of relative paths
            os.chdir(self.output_handler.get_dir('output_summaries'))
            c = CommandManager('pymol')
            c.add_command_line_arguments(['-c', self.output_handler.get_file(file_tag='pymol_sites_pml')])
            try:    c.run()
            except: print("Failed to start pymol - maybe it's not available?")
            # Change back to top directory
            os.chdir(self.out_dir)

        #os.remove(self.output_handler.get_file(file_tag='pymol_sites_pml'))
        #os.remove(self.output_handler.get_file(file_tag='pymol_sites_py'))

    def add_event_to_event_table(self, d_handler, event):
        """Add event entries to the event table"""

        # Check event has not been added previously
        assert event.id not in self.tables.event_info.index.values.tolist(), 'Event Already Added!: {!s}'.format(event.id)
        # Add values to a new row in the table
        self.tables.event_info.loc[event.id,:] = None
        # Default to site_idx of 0 if no site given
        if event.parent:    site_idx = event.parent.id
        else:               site_idx = 0
        self.tables.event_info.set_value(event.id, 'site_idx', site_idx)
        # Event and cluster information
        self.tables.event_info.set_value(event.id, '1-BDC',  round(1.0-event.info.estimated_bdc,2))
        self.tables.event_info.set_value(event.id, 'z_peak', round(event.cluster.max,2))
        self.tables.event_info.set_value(event.id, 'z_mean', round(event.cluster.mean,2))
        self.tables.event_info.set_value(event.id, 'cluster_size', event.cluster.size)
        self.tables.event_info.set_value(event.id, ['refx','refy','refz'], list(flex.double(event.cluster.peak)*self.reference_grid().grid_spacing()))
        if self.params.alignment.method=='local': mappings = self.reference_grid().partition().query_by_grid_points([map(int,event.cluster.peak)])
        else:                                     mappings = None
        self.tables.event_info.set_value(event.id, ['x','y','z'], list(d_handler.transform_from_reference(  points=flex.vec3_double([event.cluster.peak])*self.reference_grid().grid_spacing(),
                                                                                                            method=self.params.alignment.method,
                                                                                                            point_mappings=mappings       )[0]))

    def update_event_table_site_info(self, events):
        """Update the event table for pre-existing events"""
        for e in events:
            assert e.id, 'NO ID GIVEN: {!s}'.format(e.id)
            assert e.parent, 'EVENT HAS NO PARENT: {!s}'.format(e.parent)
            self.tables.event_info.set_value(e.id, 'site_idx', e.parent.id)

    def write_grid_point_distributions(self, grid_points, map_analyser, output_filename=None):
        """Write CSV file of grid points, dataset numbers and map values"""

        if not output_filename: output_filename = self.output_handler.get_file('point_distributions')

        self.log('----------------------------------->>>', True)
        self.log('Writing Grid Points distributions: {!s} points'.format(len(grid_points)), True)
        self.log('Writing to {!s}'.format(output_filename), True)

        # Add quotations around the grid points
        grid_points_str = ['\"'+str(g)+'\"' for g in grid_points]

        # Iterate through each dataset
        with open(output_filename, 'w') as fh:

            fh.write('grid_point, dataset_tag, dataset_num, map_val\n')
            for mh in map_analyser.dataset_maps.all():
                map_vals = [mh.map[g] for g in grid_points]
                # Create list of tags for zipping
                d_tags = [mh.tag]*len(grid_points)
                d_nums = [mh.num]*len(grid_points)
                # Format and write
                out_list = [', '.join(map(str,tup)) for tup in zip(grid_points_str, d_tags, d_nums, map_vals)]
                out_line = '\n'.join(map(str,out_list)) + '\n'
                fh.write(out_line)

class PanddaZMapAnalyser(object):
    def __init__(self, params, grid_spacing, log):
        self.log = log

        self.params = params
        self.grid_spacing = grid_spacing

        self.grid_clustering_cutoff = 1.1 * numpy.math.sqrt(3)
        self.real_clustering_cufoff = self.grid_clustering_cutoff * grid_spacing
        self.grid_minimum_volume = int(self.params.min_blob_volume/(grid_spacing**3))

    def print_settings(self):
        self.log('----------------------------------->>>', True)
        self.log('Z-Scores Clustering', True)
        self.log('----------------------------------->>>', True)
        self.log('Clustering Points with Z-Scores > {!s}'.format(self.params.contour_level), True)
        self.log('----------------------------------->>>', True)
        self.log('Clustering Cutoff (A):       {!s}'.format(self.real_clustering_cufoff), True)
        self.log('----------------------------------->>>', True)
        self.log('Minimum Cluster Z-Peak:      {!s}'.format(self.params.min_blob_z_peak), True)
        self.log('Minimum Cluster Volume (A):  {!s}'.format(self.params.min_blob_volume), True)
        self.log('Minimum Cluster Size:        {!s}'.format(self.grid_minimum_volume), True)

    def cluster_high_z_values(self, z_map, point_mask):
        """Finds all the points in the z-map above `z_cutoff`, points will then be clustered into groups of cutoff `clustering_cutoff` angstroms"""

        point_mask_gps = flex.vec3_double(point_mask)
        # Convert the gps to idxs
        point_mask_idx = flex.size_t(map(z_map.accessor(), point_mask))
        # Select these values from the map
        point_mask_val = z_map.select(point_mask_idx)
        # Find values above cutoff
        above_bool = (point_mask_val >= self.params.contour_level)
        above_val = point_mask_val.select(above_bool)
        above_gps = point_mask_gps.select(above_bool)
        above_len = len(above_val)

        # No Cluster points found
        if   above_len == 0:
            return 0, []
        # One Cluster point found
        elif above_len == 1:
            return 1, [(above_gps, above_val)]
        # Can't cluster if there are too many points
        elif above_len > 10000:
            return -1, [(above_gps, above_val)]
        # Cluster points if we have found them
        else:
            self.log('> Clustering {!s} Points.'.format(above_len))
            # Cluster the extracted points
            t1 = time.time()
            cluster_ids = scipy.cluster.hierarchy.fclusterdata( X = above_gps,
                                                                t = self.grid_clustering_cutoff,
                                                                criterion = 'distance',
                                                                metric    = 'euclidean',
                                                                method    = 'single' )
            cluster_ids = list(cluster_ids)
            t2 = time.time()
            self.log('> Clustering > Time Taken: {!s} seconds'.format(int(t2-t1)))

            # Get the number of clusters
            num_clusters = max(cluster_ids)
            # Group the values by cluster id
            z_clusters = []
            for c_id, c_idxs in generate_group_idxs(cluster_ids):
                c_idxs = flex.size_t(c_idxs)
                c_gps = above_gps.select(c_idxs)
                c_val = above_val.select(c_idxs)
                z_clusters.append((c_gps, c_val))
            assert num_clusters == len(z_clusters)
            return num_clusters, z_clusters

    def validate_clusters(self, z_clusters):
        for i, (gps, vals) in enumerate(z_clusters):
            #print('Cluster {}'.format(i))
            #print('Points: {}'.format(len(gps)))
            #print('Values: {}'.format(len(vals)))
            assert len(gps) == len(vals)

    def filter_z_clusters_1(self, z_clusters):
        """Filter the z-clusters on a variety of criteria (size, peak value)"""

        self.log('----------------------------------->>>')
        self.log('Filtering by blob size and peak value')

        filt_z_clusters = z_clusters
        # Filter out small clusters - get numbers of clusters satisfying the minimum cluster size
        large_clusters = (flex.int([x[1].size() for x in filt_z_clusters]) >= self.grid_minimum_volume).iselection()
        if large_clusters.size() == 0:  return 0, []
        filt_z_clusters = [filt_z_clusters[i] for i in large_clusters]
        # Filter out weak clusters - get numbers of clusters satisfying the minimum z_peak value
        strong_clusters = (flex.double([x[1].min_max_mean().max for x in filt_z_clusters]) >= self.params.min_blob_z_peak).iselection()
        if strong_clusters.size() == 0: return 0, []
        filt_z_clusters = [filt_z_clusters[i] for i in strong_clusters]

        self.log('Filtered {!s} Clusters to {!s} Clusters'.format(len(z_clusters), len(filt_z_clusters)))
        return len(filt_z_clusters), filt_z_clusters

    def filter_z_clusters_2(self, z_clusters, grid_origin_cart, ref_structure, min_contact_dist=6):
        """Find and remove clusters more than a minimum distance from the protein"""

        # min_contact_dist - blobs are rejected if they are more than this distance from the protein

        self.log('----------------------------------->>>')
        self.log('Filtering by minimum distance from protein')

        # Extract the protein sites
        ref_sites_cart = ref_structure.select(ref_structure.atom_selection_cache().selection('pepnames')).atoms().extract_xyz()
        # Save time - calculate the square of the contact distance
        min_contact_dist_sq = min_contact_dist**2

        # Remove any clusters that are more than min_contact_dist from the protein
        filtered_c_idxs = []
        for c_idx, (c_gps, c_val) in enumerate(z_clusters):
            # Extract points in cluster
            cluster_points_cart = (c_gps * self.grid_spacing) - grid_origin_cart
            # Calculate minimum distance to protein
            for r_site_cart in ref_sites_cart:
                diff_vecs_cart = cluster_points_cart - r_site_cart
                # Keep cluster if minimum distance is less than min_contact_dist
                if min(diff_vecs_cart.dot()) < min_contact_dist_sq:
                    filtered_c_idxs.append(c_idx)
                    break
            # Report
#            if self.log.verbose:
#                if filtered_c_idxs and (filtered_c_idxs[-1] == c_idx):
#                    print('KEEPING CLUSTER:', c_idx)
#                else:
#                    print('REJECTING CLUSTER:', c_idx)
        # Select filtered clusters
        filt_z_clusters = [z_clusters[i] for i in filtered_c_idxs]

        self.log('Filtered {!s} Clusters to {!s} Clusters'.format(len(z_clusters), len(filt_z_clusters)))
        return len(filt_z_clusters), filt_z_clusters

    def filter_z_clusters_3(self, z_clusters, grid_origin_cart, ref_unit_cell, ref_sym_ops, ref_structure, max_contact_dist=8):
        """Find and remove symmetry equivalent clusters"""

        if len(z_clusters) == 1:
            return 1, z_clusters
        else:
            self.log('----------------------------------->>>')
            self.log('Filtering symmetry equivalent clusters')

        # Extract the protein sites
        ref_sites_cart = ref_structure.select(ref_structure.atom_selection_cache().selection('pepnames')).atoms().extract_xyz()

        # Cartesianise and fractionalise the points in each of the clusters (in the crystallographic frame)
        points_cart = [None]*len(z_clusters)
        points_frac = [None]*len(z_clusters)
        for c_idx, (c_gps, c_val) in enumerate(z_clusters):
            # Extract points in cluster
            points_cart[c_idx] = (c_gps * self.grid_spacing) - grid_origin_cart
            # Fractionalise them to the unit cell of the reference structure
            points_frac[c_idx] = ref_unit_cell.fractionalize(points_cart[c_idx])
        # Find the sets of clusters that are symmetry related
        sym_equiv_groups = find_symmetry_equivalent_groups( points_frac = points_frac,
                                                            sym_ops     = ref_sym_ops,
                                                            unit_cell   = ref_unit_cell,
                                                            cutoff_cart = 1.05*1.7321*self.grid_spacing )
        # max_contact_dist - a point contacts an atom if the atoms is within this distance of it
        # Save time - calculate the square of the contact distance
        max_contact_dist_sq = max_contact_dist**2
        # Iterate through and chose one from each group to keep
        filt_z_clusters = []
        for g_id, g_idxs in generate_group_idxs(sym_equiv_groups):
            # Count the number of contact for each cluster in the group
            c_contacts = []
            # Iterate through cluster in the group
            for c_idx in g_idxs:
                # Initialise contact counter
                contacts = 0
                # Get the cartesian points for the cluster
                c_points_cart = points_cart[c_idx]
                # Again, use the brute force all-v-all method
                for rp in ref_sites_cart:
                    diffs_cart = c_points_cart - rp
                    # Check to see if site closer to cluster than minimum
                    if min(diffs_cart.dot()) < max_contact_dist_sq:
                        contacts += 1
                # Record the number of contacts (over size of cluster)
                c_contacts.append(1.0*contacts/len(c_points_cart))
#                if self.log.verbose:
#                    print('CLUSTER:', c_idx, ', CONTACTS PER POINT:', round(c_contacts[-1],3))

            # Find the cluster with the most contacts
            max_contacts = max(c_contacts)
            if max_contacts == 0:
                raise Exception('MAX CONTACTS IS 0!')
            else:
                cluster_to_keep = g_idxs[c_contacts.index(max_contacts)]
                filt_z_clusters.append(z_clusters[cluster_to_keep])
#                if self.log.verbose:
#                    print('KEEPING CLUSTER', cluster_to_keep)
        assert len(filt_z_clusters) == max(sym_equiv_groups), 'NUMBER OF UNIQUE GROUPS AND GROUPS TO BE RETURNED NOT THE SAME'

        self.log('Filtered {!s} Clusters to {!s} Clusters'.format(len(z_clusters), len(filt_z_clusters)))
        return len(filt_z_clusters), filt_z_clusters

    def group_clusters(self, z_clusters, separation_cutoff=5):
        """Join clusters that are separated by less than max_separation"""

        if len(z_clusters) == 1:
            return 1, z_clusters
        else:
            self.log('----------------------------------->>>')
            self.log('Grouping Nearby Clusters')

        # Minimum distance between grid points to be joined (squared)
        grid_cutoff_sq = (separation_cutoff/self.grid_spacing)**2

        # Record which clusters are to be joined
        connect_array = numpy.zeros((len(z_clusters),len(z_clusters)), dtype=int)
        for i_clust_1, (c_gps_1, c_val_1) in enumerate(z_clusters):
            for i_clust_2, (c_gps_2, c_val_2) in enumerate(z_clusters):
                # Skip if this is the same blob
                if i_clust_1 == i_clust_2:
                    connect_array[(i_clust_1, i_clust_2)] = 1
                    continue
                # Extract the minimum separation of the grid points
                min_dist_sq = min([min((c_gps_2 - gp).dot()) for gp in c_gps_1])
                # Check to see if they should be joined
                if min_dist_sq < grid_cutoff_sq:
                    connect_array[(i_clust_1, i_clust_2)] = 1
        # Cluster the connection array
        cluster_groupings = find_connected_groups(connection_matrix=connect_array)
        # Concatenate smaller clusters into larger clusters
        grouped_clusters = []
        for g_id, g_idxs in generate_group_idxs(cluster_groupings):
            g_gps = []; [g_gps.extend(z_clusters[i][0]) for i in g_idxs]
            g_gps = flex.vec3_double(g_gps)
            g_val = []; [g_val.extend(z_clusters[i][1]) for i in g_idxs]
            g_val = flex.double(g_val)
            grouped_clusters.append((g_gps, g_val))

        assert len(grouped_clusters) == max(cluster_groupings)

        self.log('Grouped {!s} Clusters together to form {!s} Clusters'.format(len(z_clusters), len(grouped_clusters)))
        return len(grouped_clusters), grouped_clusters


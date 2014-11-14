import os, sys, glob
import numpy

import iotbx.pdb as pdb_reader
import iotbx.map_tools as map_tools

from cctbx import maptbx
from scitbx.array_family import flex
from scitbx.math import superpose, basic_statistics
from libtbx.math_utils import ifloor, iceil
from iotbx.reflection_file_utils import extract_miller_array_from_file

from Giant.Xray.Miller.Utils import apply_simple_scaling
from Giant.Xray.Structure.Select import get_calpha_sites
from Giant.Xray.Maps.Utils import write_1d_array_as_p1_map
from Giant.Xray.Maps.Utils import get_fft_map_from_f_obs_and_structure
from Giant.Xray.Maps.Grid import get_bounding_box_for_structure, create_cartesian_grid, calculate_sampling_distance
from Giant.Stats.Normalise import normalise_array_to_z_scores

from Giant.Xray.Maps.Grid import get_grid_points_within_distance_cutoff_of_origin, combine_grid_point_and_grid_vectors
from Giant.Xray.Maps.Grid import get_grid_points_within_distance_cutoff_of_cart_sites
from Giant.Stats.Tests import test_significance_of_group_of_z_values, convert_pvalue_to_zscore
from Giant.Stats.Utils import resample_ordered_list_of_values

def status_bar(n, n_max):
    if (n+1==n_max):
        print '\r>>', 100, '%'
    elif (n%int(n_max/100)==0):
        print '\r>>', int(round(100.0*(n+1)/n_max,0)), '%',
    sys.stdout.flush()

class multi_dataset_analyser(object):
    def __init__(self, verbose=True):
        """Class for the processing of datasets from a fragment soaking campaign"""

        self.verbose = verbose

        self._map_type = '2mFo-DFc'
        self._res_factor = None
        self._cut_resolution = None

        self._raw_file_pairs = []
        self._file_pairs = []

        self._raw_datasets = []
        self._datasets = []

        self._ref_dataset = None
        self._ref_grid = None
        self._ref_map = None

        self._mean_map = None
        self._stds_map = None
        self._skew_map = None
        self._kurt_map = None

        self._maps = []
        self._z_maps = []
        self._mod_z_maps = []

        self._log = ''

    def log(self, message, show=False):
        """Log message to file, and mirror to stdout if verbose or force_print"""

        if self.verbose or show:
            print(message)
        self._log = self._log + message + '\n'

    def print_log(self):
        print self._log

    def set_map_type(self, map_type):
        self._map_type = map_type
    def get_map_type(self):
        return self._map_type
    def set_res_factor(self, res_factor):
        self._res_factor = res_factor
    def get_res_factor(self):
        return self._res_factor
    def set_cut_resolution(self, d_min):
        self._cut_resolution = d_min
    def get_cut_resolution(self):
        return self._cut_resolution

    def get_input_files(self):
        """Get all of the files that were added"""
        return self._raw_file_pairs
    def get_used_files(self):
        """Get the files that have been used to generate the distributions"""
        return self._file_pairs

    def set_reference_dataset(self, dataset):
        """Set a reference dataset created externally"""
        self._ref_dataset = dataset
    def get_reference_dataset(self):
        """Get the dataset used as the reference for scaling and aligning the other datasets"""
        return self._ref_dataset
    def set_reference_grid(self, ref_grid):
        """Set a grid created externally"""
        self._ref_grid = ref_grid
        self._ref_grid.show_summary()
    def get_reference_grid(self):
        """Get the reference grid used to generate the points for sampling the maps"""
        return self._ref_grid

    def get_all_datasets(self):
        """Return all datasets added"""
        return self._raw_datasets
    def get_used_datasets(self):
        """Return datasets used to build distributions"""
        return self._datasets

    def get_mean_map(self):
        """Returns the average map across the datasets"""
        return self._mean_map
    def get_stds_map(self):
        """Returns the std dev map across the datasets"""
        return self._stds_map
    def get_skew_map(self):
        """Returns the skewness map across the datasets"""
        return self._skew_map
    def get_kurt_map(self):
        """Returns the kurtosis map across the datasets"""
        return self._kurt_map
    # TODO BIMODALITY MAP

    def get_maps(self):
        """Returns maps for the used datasets"""
        return self._maps
    def get_z_maps(self):
        """Returns z-maps for the used datasets"""
        return self._z_maps
    def get_modified_z_maps(self):
        """Returns the post-processed, modified, z-maps for the used datasets"""
        return self._mod_z_maps

    def load_reference_dataset(self, ref_pdb, ref_mtz):
        """Set the reference dataset, to which all other datasets will be aligned and scaled"""

        self._ref_dataset = dataset_handler(ref_pdb, ref_mtz)
        self._ref_dataset.create_fft_map(map_type=self.get_map_type())
        self._ref_dataset.get_map().apply_volume_scaling()
        self._ref_dataset.create_map_handler()
    def extract_reference_map_values(self):
        """Read in map for the reference dataset"""

        ref_map_vals = self.get_reference_dataset().get_map_handler().get_cart_values(self.get_reference_grid().get_cart_points())
        ref_map_file = self.get_reference_dataset().get_mtz_filename().replace('.mtz','.sampled.ccp4')
        self.write_array_to_map(output_file=ref_map_file, map_data=ref_map_vals)
        self._ref_map = ref_map_vals
    def create_reference_grid(self, res_factor=None, include_origin=True, buffer=None):
        """Create a grid over the reference protein"""

        if res_factor:
            self.set_res_factor(res_factor)
        grid_spacing = calculate_sampling_distance(self.get_cut_resolution(), self.get_res_factor())
        min_max_sites = self.get_reference_dataset().get_structure_min_max_sites()

        self._ref_grid = grid_handler()
        self._ref_grid.set_grid_spacing(spacing=grid_spacing)
        self._ref_grid.set_cart_extent(cart_min=min_max_sites[0], cart_max=min_max_sites[1])
        self._ref_grid.create_cartesian_grid(include_origin=include_origin)

        self._ref_grid.show_summary()
    def mask_resampled_reference_grid(self):
        """Using the local and global masks, mask the resamples grid points"""

        # Create a grid mask around each point
        self.get_reference_grid().show_summary()
        self.get_reference_grid().get_local_mask().show_summary()

        # Get the grid points that are not masked, and points in the buffer
        sample_points = self.get_reference_grid().get_sample_grid_points()
        buffer_mask   = self.get_reference_grid().get_buffer_zone_points()
        global_mask   = self.get_reference_grid().get_global_mask_points()

        masked_sample_points_1 = sample_points
        print('Resampled Grid Size (3D): {!s}'.format(self.get_reference_grid().get_resampled_grid_size()))
        print('Resampled Grid Size (1D): {!s}'.format(len(masked_sample_points_1)))
        print('Filtering Buffer Mask... (Edge of Cell)')
        masked_sample_points_2 = [p for p in masked_sample_points_1 if p not in buffer_mask]
        print('Filtered Points: {!s}'.format(len(masked_sample_points_2)))
        print('Filtering Global Mask... (Protein)')
        masked_sample_points_3 = [p for p in masked_sample_points_2 if p in global_mask]
        print('Filtered Points: {!s}'.format(len(masked_sample_points_3)))

        self.get_reference_grid().set_masked_grid_points(resampled_points_3)

    def add_input_files(self, file_pairs):
        """Add (pdb, mtz) file pairs to the datasets to be processed"""

        print '===================================>>>'
        self._raw_file_pairs = file_pairs
        print '{!s} datasets added'.format(len(file_pairs))

    def load_input_datasets(self):
        """Read in maps for the input datasets"""

        for d_num, (pdb, mtz) in enumerate(self.get_input_files()):
            print '===================================>>>'
            print 'Loading Dataset {!s}'.format(d_num+1)
            print pdb
            print mtz
            # Create object to handle processing of the dataset
            new_dataset = dataset_handler(pdb, mtz)
            # Add dataset to the list of all datasets
            self._raw_datasets.append(new_dataset)

    def scale_datasets_and_load_maps(self):
        """Iterate through the datasets, scale them to the reference, and extract maps"""

        for d_num, new_dataset in enumerate(self.get_all_datasets()):
            print '===================================>>>'
            print 'Scaling Dataset {!s}'.format(d_num+1)
            print '===============>>>'
            # Scale new data to the reference dataset
            new_dataset.scale_fobs_to_reference(ref_miller=self.get_reference_dataset().get_fobs_miller_array())
            # Create maps
            new_dataset.create_fft_map(miller=new_dataset.get_scaled_fobs_miller_array(), map_type=self.get_map_type(), d_min=self.get_cut_resolution())
            new_dataset.get_map().apply_volume_scaling()
            new_dataset.create_map_handler()

    def align_and_filter_datasets(self):
        """Iterate through the datasets, and filter them by quality"""

        for d_num, new_dataset in enumerate(self.get_all_datasets()):
            # Retrieve associated files
            pdb, mtz = self.get_input_files()[d_num]
            print '===================================>>>'
            print 'Filtering Dataset {!s}'.format(d_num+1)
            print '===============>>>'
            # Check that it correlates well with itself before and after scaling
            scaled_correlation = new_dataset.get_scaled_fobs_miller_array().correlation(new_dataset.get_fobs_miller_array()).coefficient()
            print 'Scaled-Unscaled Correlation: {!s}'.format(scaled_correlation)
            if scaled_correlation < 0.9:
                print 'Low correlation between scaled and unscaled data - Removing Dataset!'
                raise Exception('Low correlation dataset?')
                continue
            # Check the resolution of the dataset
            print 'Ref Resolution: ', self.get_reference_dataset().get_map().d_min()
            print 'Map Resolution: ', new_dataset.get_map().d_min()
            if new_dataset.get_map().d_min() > self.get_cut_resolution():
                print 'Does not meet high-resolution cutoff - Removing Dataset!'
                raise Exception('Low resolution dataset?')
                continue

            # Align to reference structure to get mapping transform
            new_dataset.align_to_reference(reference_calphas=self.get_reference_dataset().get_calpha_sites())

            # TODO Check that the transformation is SMALL and filter for large movements?
#            alignment_fit = new_dataset.get_alignment_transform()
#            print 'R:\t', '\n\t'.join([' '.join(map(str,l)) for l in alignment_fit.r.as_list_of_lists()])
#            print 'T:\t', '\n\t'.join([' '.join(map(str,l)) for l in alignment_fit.t.as_list_of_lists()])

            print 'Initial RMSD: {!s}'.format(new_dataset.get_calpha_sites().rms_difference(self.get_reference_dataset().get_calpha_sites()))
            print 'Aligned RMSD: {!s}'.format(new_dataset.get_calpha_sites().rms_difference(new_dataset.transform_points_from_reference(self.get_reference_dataset().get_calpha_sites())))

            self._file_pairs.append((pdb, mtz))
            self._datasets.append(new_dataset)

    def extract_map_values(self):
        """Extract map values for the filtered datasets"""

        for d_num, new_dataset in enumerate(self.get_used_datasets()):
            # Retrieve associated files
            pdb, mtz = self.get_used_files()[d_num]
            print '===================================>>>'
            print 'Loading Dataset Maps {!s}'.format(d_num+1)
            print '===============>>>'
            # Extract the map values at the transformed grid points
            new_sample_points = new_dataset.transform_points_from_reference(self.get_reference_grid().get_cart_points())
            new_map_vals = new_dataset.get_map_handler().get_cart_values(new_sample_points)

            # Write the sampled map
            new_sampled_map_file = mtz.replace('.mtz','.sampled.ccp4')
            self.write_array_to_map(output_file = new_sampled_map_file,
                                    map_data = new_map_vals,
                                    grid_size = self.get_reference_grid().get_grid_size(),
                                    grid_spacing = self.get_reference_grid().get_grid_spacing())

            # Append the dataset, the map values to common lists
            self._maps.append(new_map_vals)

        print '===================================>>>'
        print 'After Filtering: {!s} Datasets'.format(len(self.get_used_files()))

    def calculate_map_statistics(self):
        """Take the sampled maps and calculate statistics for each grid point across the datasets"""

        print '===================================>>>'
        print 'Combining maps into array...',; sys.stdout.flush()
        map_array = numpy.array(self.get_maps())
        print 'finished.'

        print 'Calculating Statistics of Grid Points'

        # Calculate Mean Maps
#        mean_map_vals = map_array.mean(axis=0)
        mean_map_vals = numpy.array([basic_statistics(flex.double(map_array[:,i].tolist())).mean for i in xrange(self.get_reference_grid().get_grid_size_1d())])
        assert len(mean_map_vals) == self.get_reference_grid().get_grid_size_1d()
        mean_map_file = self.get_reference_dataset().get_mtz_filename().replace('.mtz','.mean.ccp4')
        self.write_array_to_map(output_file=mean_map_file, map_data=flex.double(mean_map_vals),
                                grid_size=self.get_reference_grid().get_grid_size(), grid_spacing=self.get_reference_grid().get_grid_spacing())

        # Calculate Stds Maps
#        stds_map_vals = map_array.std(axis=0)
        stds_map_vals = numpy.array([basic_statistics(flex.double(map_array[:,i].tolist())).biased_standard_deviation for i in xrange(self.get_reference_grid().get_grid_size_1d())])
        assert len(stds_map_vals) == self.get_reference_grid().get_grid_size_1d()
        stds_map_file = self.get_reference_dataset().get_mtz_filename().replace('.mtz','.stds.ccp4')
        self.write_array_to_map(output_file=stds_map_file, map_data=flex.double(stds_map_vals),
                                grid_size=self.get_reference_grid().get_grid_size(), grid_spacing=self.get_reference_grid().get_grid_spacing())

        # Calculate Skew Maps
        skew_map_vals = numpy.array([basic_statistics(flex.double(map_array[:,i].tolist())).skew for i in xrange(self.get_reference_grid().get_grid_size_1d())])
        assert len(skew_map_vals) == self.get_reference_grid().get_grid_size_1d()
        skew_map_file = self.get_reference_dataset().get_mtz_filename().replace('.mtz','.skew.ccp4')
        self.write_array_to_map(output_file=skew_map_file, map_data=flex.double(skew_map_vals),
                                grid_size=self.get_reference_grid().get_grid_size(), grid_spacing=self.get_reference_grid().get_grid_spacing())

        # Calculate Kurtosis Maps
        kurt_map_vals = numpy.array([basic_statistics(flex.double(map_array[:,i].tolist())).kurtosis for i in xrange(self.get_reference_grid().get_grid_size_1d())])
        assert len(kurt_map_vals) == self.get_reference_grid().get_grid_size_1d()
        kurt_map_file = self.get_reference_dataset().get_mtz_filename().replace('.mtz','.kurt.ccp4')
        self.write_array_to_map(output_file=kurt_map_file, map_data=flex.double(kurt_map_vals),
                                grid_size=self.get_reference_grid().get_grid_size(), grid_spacing=self.get_reference_grid().get_grid_spacing())

        # Store map vals
        self._mean_map = mean_map_vals.tolist()
        self._stds_map = stds_map_vals.tolist()
        self._skew_map = skew_map_vals.tolist()
        self._kurt_map = kurt_map_vals.tolist()

    def normalise_maps_to_z_maps(self):
        """Normalise the map values to z-values"""

        for d_num, map_data in enumerate(self.get_maps()):
            print '===================================>>>'
            print 'Normalising Dataset {!s}'.format(d_num+1)

            pdb, mtz = self.get_used_files()[d_num]

            # Create Z-scores
            z_array = normalise_array_to_z_scores(input_array=numpy.array(map_data), element_means=self.get_mean_map(), element_stds=self.get_stds_map())
            # Store data
            self._z_maps.append(z_array.tolist())

            # Calculate stats for the map
            z_stats = basic_statistics(flex.double(z_array.tolist()))

            # Write map
            z_map_file = mtz.replace('.mtz','.zvalues.ccp4')
            self.write_array_to_map(output_file=z_map_file, map_data=flex.double(z_array),
                                grid_size=self.get_reference_grid().get_grid_size(), grid_spacing=self.get_reference_grid().get_grid_spacing())

    def post_process_z_maps(self):
        """Process the z_maps, looking for groups of high z-values"""

        print '===================================>>>'
        print 'Combining maps into array...',; sys.stdout.flush()
        z_map_array = numpy.array(self.get_z_maps())
        print 'finished.'

        # Get the masked grid points
        masked_points = self.get_reference_grid().get_masked_grid_points()
        # Translate between 3d grid point and 1d array index
        grid_indexer = self.get_reference_grid().get_grid_indexer()

        # Get the size of the down-sampled grid
        resampled_grid_points = self.get_reference_grid().get_resampled_grid_points()
        resampled_grid_size = self.get_reference_grid().get_resampled_grid_size()

        for d_num, z_map_data in enumerate(self.get_z_maps()):
            print '===================================>>>'
            print 'Post-Processing Dataset {!s}'.format(d_num+1)

            mod_z_map_data = [0]*self.get_reference_grid().get_grid_size_1d()

            # Iterate through all grid points and calculate modified z-value at each point (0 at unsampled points or in buffer zone)
            for p_num, gp in enumerate(masked_points):

                status_bar(n=p_num, n_max=len(masked_points))

                # Find nearby values within local mask
                gp_masked = self.get_reference_grid().get_local_mask().apply_mask(gp)
                # Find the map values at these points
                local_map_vals = [z_map_data[grid_indexer(gp_m)] for gp_m in gp_masked]
                # Down-sample recorded values
                resamp_map_vals = resample_ordered_list_of_values(local_map_vals)
                # Calculate significance of values
                pval = test_significance_of_group_of_z_values(resamp_map_vals)
                # Convert to Z-score
                mod_z_val = convert_pvalue_to_zscore(pval)

                # Add to new map array
                mod_z_map_data[grid_indexer(gp)] = mod_z_val

            pdb, mtz = self.get_used_files()[d_num]

            # Write map
            mod_z_map_file = mtz.replace('.mtz','.modzvalues.ccp4')
            self.write_array_to_map(output_file=mod_z_map_file, map_data=flex.double(mod_z_map_data))

            # Get down-sampled map
            resamp_mod_z_map_data = [mod_z_map_data[grid_indexer(p)] for p in resampled_grid_points]
            assert resampled_grid_size[0]*resampled_grid_size[1]*resampled_grid_size[2] == len(samp_mod_z_map_data)

            # Calculate stats for the map
            # TODO Change this so that it's only over the masked points
#            resamp_mod_z_stats = basic_statistics(flex.double(samp_mod_z_map_data))

            # Write down-sampled map
            resamp_mod_z_map_file = mtz.replace('.mtz','.resamp.modzvalues.ccp4')
            self.write_array_to_map(output_file=resamp_mod_z_map_file, map_data=flex.double(resamp_mod_z_map_data),
                                    grid_size=resampled_grid_size, grid_spacing=self.get_reference_grid().get_resampled_grid_spacing())

    def write_array_to_map(self, output_file, map_data, grid_size=None, grid_spacing=None):
        """Takes a 1d array and writes it to a map"""
        if not grid_size:    grid_size    = self.get_reference_grid().get_grid_size()
        if not grid_spacing: grid_spacing = self.get_reference_grid().get_grid_spacing()
        print 'Writing Map: {!s}'.format(output_file)
        write_1d_array_as_p1_map(file_name=output_file, map_data=map_data, grid_size=grid_size, grid_spacing=grid_spacing)

    def write_array_to_file(self, output_file, array, nrow):
        """Takes an array and writes it to file"""

        print type(array)
        print 'Writing Array: {!s}'.format(output_file)
        with open(output_file, 'w') as output_fh:
            # Write (x,y,z) tuples on each line
            for i in xrange(nrow):
                status_bar(n=i, n_max=nrow)
                output_line = [str(round(x,3)) for x in array[i]]
                output_fh.write(', '.join(output_line)+'\n')

class dataset_handler(object):
    def __init__(self, pdb_filename, mtz_filename):
        """Create a dataset object to allow common functions to be applied easily to a number of datasets"""

        assert os.path.exists(pdb_filename), 'PDB file does not exist!'
        assert os.path.exists(mtz_filename), 'MTZ file does not exist!'

        # Store filenames
        self._pdb_file = os.path.abspath(pdb_filename)
        self._mtz_file = os.path.abspath(mtz_filename)
        # PDB Objects
        self._pdb_input = pdb_reader.input(source_info=None,lines=open(self._pdb_file,'r').read())
        self._pdb_struc = self._pdb_input.xray_structure_simple()
        self._min_max_sites = get_bounding_box_for_structure(self.get_pdb_structure())

        # Store C-alpha sites
        self._calpha_sites = get_calpha_sites(input_obj=self._pdb_input, structure_obj=self._pdb_struc)

        # Extract miller array
        self._fobs_miller = extract_miller_array_from_file(self._mtz_file, 'F,SIGF', log=open(os.devnull,'a'))
        self._scaled_fobs_miller = None

        # Initialise other variables
        self._fft_map = None
        self._unit_cell = None
        self._space_group = None
        self._basic_map = None
        self._rt_transform = None

    def get_pdb_filename(self):
        return self._pdb_file
    def get_mtz_filename(self):
        return self._mtz_file
    def get_pdb_input(self):
        return self._pdb_input
    def get_pdb_structure(self):
        return self._pdb_struc
    def get_calpha_sites(self):
        return self._calpha_sites
    def get_alignment_transform(self):
        return self._rt_transform
    def get_structure_min_max_sites(self):
        return self._min_max_sites
    def get_fobs_miller_array(self):
        return self._fobs_miller
    def get_scaled_fobs_miller_array(self):
        return self._scaled_fobs_miller
    def get_dataset_resolution(self):
        return self.get_fobs_miller_array().d_min()
    def get_map(self):
        return self._fft_map
    def get_map_handler(self):
        return self._basic_map

    def transform_points_from_reference(self, points):
        """Use alignment to map from reference frame to our frame"""
        return self.get_alignment_transform().rt().inverse() * points

    def transform_points_to_reference(self, points):
        """Use alignment to map to reference frame from our frame"""
        return self.get_alignment_transform().rt() * points

    def scale_fobs_to_reference(self, ref_miller):
        """Scale the observed data to another dataset"""
        self._scaled_fobs_miller = apply_simple_scaling(self.get_fobs_miller_array(), ref_miller)
        return self._scaled_fobs_miller

    def align_to_reference(self, reference_calphas):
        """Calculate the rotation and translation needed to align one structure to another"""
        self._rt_transform = superpose.least_squares_fit(reference_sites=reference_calphas, other_sites=self.get_calpha_sites())
        return self._rt_transform

    def create_fft_map(self, miller=None, map_type='2mFo-DFc', d_min=None):
        """Create an fft map to allow map values to be calculated"""
        if not miller:
            self._fft_map = get_fft_map_from_f_obs_and_structure(f_obs=self._fobs_miller, xray_structure=self._pdb_struc, map_type=map_type, d_min=d_min)
        else:
            self._fft_map = get_fft_map_from_f_obs_and_structure(f_obs=miller, xray_structure=self._pdb_struc, map_type=map_type, d_min=d_min)
        self._unit_cell = self._fft_map.unit_cell()
        self._space_group = self._fft_map.space_group()
        return self._fft_map

    def create_map_handler(self):
        """Create a map handler to allow for easy retrieval of map values"""
        if self._fft_map == None: self.create_fft_map()
        self._basic_map = maptbx.basic_map(maptbx.basic_map_unit_cell_flag(), self._fft_map.real_map(), self._fft_map.real_map().focus(),
                            self._unit_cell.orthogonalization_matrix(), maptbx.out_of_bounds_clamp(0).as_handle(), self._unit_cell)
        return self._basic_map

class grid_handler(object):
    def __init__(self):
        """Create and manage a grid object to be sampled across many aligned datasets"""

        self._grid_size = None
        self._grid_spacing = None

        self._cart_max = None
        self._cart_min = None
        self._cart_size = None
        self._cart_points = None

        self._local_mask = None
        self._global_mask = None

        self._resampled_grid_points = None
        self._resampled_grid_size = None
        self._masked_grid_points = None
        self._buffer_zone_points = None

    def set_grid_spacing(self, spacing):
        self._grid_spacing = spacing
    def get_grid_spacing(self):
        return self._grid_spacing
    def get_grid_size(self):
        return self._grid_size
    def get_grid_size_1d(self):
        return self._grid_size[0]*self._grid_size[1]*self._grid_size[2]
    def get_grid_point_volume(self):
        return self.get_grid_spacing()**3

    def set_cart_extent(self, cart_min, cart_max):
        self._cart_min = cart_min
        self._cart_max = cart_max
        self._cart_size = tuple([s1 - s2 for s1,s2 in zip(cart_max, cart_min)])
    def get_cart_max(self):
        return self._cart_max
    def get_cart_min(self):
        return self._cart_min
    def get_cart_extent(self):
        return (self._cart_min, self._cart_max)
    def get_cart_size(self):
        return self._cart_size

    def get_cart_points(self):
        return self._cart_points
    def get_grid_points(self):
        return flex.nested_loop(self.get_grid_size())

    def get_grid_indexer(self):
        """Translates between 3d grid coordinates and 1d array coordinates"""
        return flex.grid(self.get_grid_size())

    def get_local_mask(self):
        return self._local_mask
    def get_global_mask(self):
        return self._global_mask

    def set_masked_grid_points(self, masked_points):
        self._masked_grid_points = masked_points
    def get_masked_grid_points(self):
        return self._masked_grid_points

    def show_summary(self):
        print '===================================>>>'
        print('Grid Spacing:        {!s}'.format(round(self.get_grid_spacing(), 3)))
        print('Grid Point Volume:   {!s}'.format(round(self.get_grid_point_volume(),3)))
        print('Size of Grid (3D):   {!s}'.format(tuple([round(x,3) for x in self.get_grid_size()])))
        print('Size of Grid (1D):   {!s}'.format(self.get_grid_size_1d()))
        print('Size of Grid (Cart): {!s}'.format(tuple([round(x,3) for x in self.get_cart_size()])))

    def create_cartesian_grid(self, include_origin=True):
        if include_origin:
            box_size, self._grid_size, self._cart_points = create_cartesian_grid(min_carts=(0,0,0),
                                                                                 max_carts=self.get_cart_max(),
                                                                                 grid_spacing=self.get_grid_spacing())
        else:
            box_size, self._grid_size, self._cart_points = create_cartesian_grid(min_carts=self.get_cart_min(),
                                                                                 max_carts=self.get_cart_max(),
                                                                                 grid_spacing=self.get_grid_spacing())

        # Update max/min cart sizes as the grid will be slightly larger than the requested size
        self.set_cart_extent(cart_min=self.get_cart_points()[0], cart_max=self.get_cart_points()[-1])

        return self.get_grid_size()

    def create_local_mask(self, distance_cutoff, mask='sphere'):
        """Generate a sphere of grid points within distance_cutoff of the origin"""

        print('===================================>>>')
        print('Generating Local Mask')
        if mask == 'sphere':
            self._local_mask = masking_sphere(self.get_grid_spacing(), distance_cutoff=distance_cutoff)
        else:
            raise Exception('Mask not found')

    def create_global_mask(self, cart_sites, distance_cutoff):
        """Get the mask of all grid points within distance_cutoff of any of cart_sites"""

        print('===================================>>>')
        print('Generating Global Mask')
        self._global_mask = global_mask(cart_sites=cart_sites, grid_spacing=self.get_grid_spacing(), distance_cutoff=distance_cutoff)

    def get_resampled_grid_points(self):
        """Get a down-sampled list of grid points, based on the local mask used"""

        if not self._resampled_grid_points:
            # Local mask gives information about how frequently to sample
            grid_jump = self.get_local_mask().get_grid_jump()
            grid_size = self.get_grid_size()
            # Points to sample at
            self._resampled_grid_points = [gp for gp in flex.nested_loop(grid_size) if not [1 for coord in gp if (coord%grid_jump != 0)]]
            self._resampled_grid_size = tuple([int(1+(g-1)/grid_jump) for g in grid_size])
            self._resampled_grid_spacing = grid_jump*self.get_grid_spacing()

        return self._resampled_grid_points

    def get_resampled_grid_size(self):
        """Gets the size of the re-sampled grid"""
        return self._resampled_grid_size
    def get_resampled_grid_spacing(self):
        """Gets the grid spacing for the re-sampled grid"""
        return self._resampled_grid_spacing

    def get_buffer_zone_points(self):
        """Get a list of points in the buffer zone of the map where the values cannot be calculated"""

        if not self._buffer_zone_points:
            # Local mask gives information about how big the buffer needs to be
            buffer = self.get_local_mask().get_buffer_size()
            grid_size = self.get_grid_size()
            # Get which grid points are within the buffer zones
            self._buffer_zone_points = [gp for gp in flex.nested_loop(grid_size) if [1 for i_dim, coord in enumerate(gp) if (coord<buffer or coord>(grid_size[i_dim]-1-buffer))]]

        return self._buffer_zone_points

    def get_global_mask_points(self):
        """Mask the grid to a subset of points by distance cutoff from sites"""
        return self.get_global_mask().get_grid_points()

class global_mask(object):
    def __init__(self, cart_sites, grid_spacing, distance_cutoff):
        """Take a grid and calculate all grid points with a certain distance cutoff of any point in cart_sites"""
        self._mask_grid_points = get_grid_points_within_distance_cutoff_of_cart_sites(cart_sites=cart_sites, grid_spacing=grid_spacing, distance_cutoff=distance_cutoff)

    def get_grid_points(self):
        """Return the grid points allowed by the mask"""
        return self._mask_grid_points

    def get_masked_grid(self, grid_size):
        """Return a masked grid of size grid_size"""
        return [1 if p in self.get_grid_points() else 0 for p in flex.nested_loop(grid_size)]

class masking_sphere(object):
    def __init__(self, grid_spacing, distance_cutoff, grid_jump=None):
        """Sphere used to mask grid points within a certain distance of a point"""

        self._mask = get_grid_points_within_distance_cutoff_of_origin(grid_spacing=grid_spacing, distance_cutoff=distance_cutoff)
        self._radius = distance_cutoff
        self._buffer = max(max(self._mask))
        self._grid_spacing = grid_spacing
        if grid_jump:
            self._grid_jump = int(grid_jump)
        else:
            self._grid_jump = iceil(self._buffer)
        if self._grid_jump == 0:
            self._grid_jump = 1

    def get_mask(self):
        return self._mask
    def get_buffer_size(self):
        return self._buffer
    def get_grid_spacing(self):
        return self._grid_spacing
    def get_grid_jump(self):
        return self._grid_jump
    def get_radius(self):
        return self._radius
    def get_volume(self):
        return (4/3.0)*numpy.pi*(self.get_radius()**3)

    def apply_mask(self, grid_point):
        """Combine a grid point with all of the masking vectors"""
        return combine_grid_point_and_grid_vectors(start_point=grid_point, grid_vectors=self.get_mask())

    def show_summary(self):
        print('===================================>>>')
        print('Number of Mask Points:  {!s}'.format(len(self.get_mask())))
        print('Mask Radius (Cart):     {!s}'.format(self.get_radius()))
        print('Mask Volume (Cart):     {!s}'.format(self.get_volume()))
        print('Largest Mask Vector:    {!s}'.format(max(self.get_mask())))
        print('Req. Edge Buffer Zone:  {!s}'.format(self.get_buffer_size()))
        print('Sampling Fraction (1D): 1/{!s}'.format(self.get_grid_jump()))
        print('Sampling Fraction (3D): 1/{!s}'.format(self.get_grid_jump()**3))

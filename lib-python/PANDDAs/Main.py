import os, sys, glob, time, re
import copy, resource, gc
import multiprocessing

from scipy import stats as scipy_stats
from scipy import cluster as scipy_cluster
import numpy

import iotbx.pdb as pdb_reader
#import iotbx.map_tools as map_tools

from libtbx import easy_pickle
from libtbx import easy_mp

from cctbx import maptbx
from cctbx import crystal
from scitbx import sparse
from scitbx.array_family import flex
from scitbx.math import superpose, basic_statistics
from libtbx.math_utils import ifloor, iceil
from iotbx.reflection_file_utils import extract_miller_array_from_file

from Bamboo.Common import compare_dictionaries
from Bamboo.Common.File import output_file_object, easy_directory
from Bamboo.Common.Data import data_collection, multiple_data_collection
from Bamboo.Common.Masks import mask_collection

from Giant.Xray.Miller.Utils import apply_simple_scaling
from Giant.Xray.Structure.Select import get_calpha_sites, get_backbone_sites
from Giant.Xray.Maps.Utils import write_1d_array_as_p1_map
from Giant.Xray.Maps.Utils import get_fft_map_from_f_obs_and_structure
from Giant.Grid.Utils import create_cartesian_grid, calculate_sampling_distance
from Giant.Stats.Normalise import normalise_array_to_z_scores

from Giant.Grid.Utils import get_grid_points_within_distance_cutoff_of_origin, combine_grid_point_and_grid_vectors
from Giant.Grid.Utils import get_grid_points_within_distance_cutoff_of_cart_sites
from Giant.Grid.Utils import get_grid_points_within_distance_cutoff_of_cart_sites_2
from Giant.Stats.Tests import test_significance_of_group_of_z_values, convert_pvalue_to_zscore
from Giant.Stats.Utils import resample_ordered_list_of_values, calculate_minimum_redundancy

from Giant.Stats.Cluster import cluster_data, combine_clusters

from Giant.Utils import status_bar

from PANDDAs.HTML import PANDDA_HTML_ENV
from PANDDAs.Settings import PANDDA_TOP, PANDDA_TEXT

STRUCTURE_MASK_NAMES = [    'bad structure - chain counts',
                            'bad structure - chain ids',
                            'bad structure - chain sequences',
                            'bad structure - residue counts',
                            'bad structure - atom counts',
                            'bad structure - non-identical structures'  ]
CRYSTAL_MASK_NAMES   = [    'bad crystal - isomorphous crystal',
                            'bad crystal - isomorphous structure',
                            'bad crystal - data quality',
                            'bad crystal - data correlation',
                            'bad crystal - resolution',
                            'bad crystal - rfree'  ]
REJECT_MASK_NAMES    = [    'rejected - total',
                            'rejected - crystal',
                            'rejected - structure',
                            'rejected - unknown'  ]

FLAG_MASK_NAMES      = [    'noisy zmap',
                            'interesting'  ]

def structure_summary(hierarchy):
    """Return a summary of a structure"""

    # Dictionary to be returned
    s_dict = {}

    # Counts of different things in the structure
    o_counts = hierarchy.overall_counts()

    assert o_counts.chain_ids.values().count(1) == len(o_counts.chain_ids.values()), 'CHAINS MAY ONLY BE PRESENT ONCE IN STRUCTURE'

    # Chain Summaries
    s_dict['chain_sequences'] = {}
    s_dict['chain_classes'] = {}
    s_dict['chain_is_protein'] = {}

    for chain in hierarchy.chains():
        s_dict['chain_sequences'][chain.id] = ''.join(chain.as_sequence())
        s_dict['chain_classes'][chain.id] = chain.get_residue_names_and_classes()[1]
        s_dict['chain_is_protein'][chain.id] = chain.is_protein()

    return o_counts, s_dict

class structure_comparer(object):
    def __init__(self, hierarchy):
        """Object for summarising a structure and comparing it to others"""
        self.hierarchy = hierarchy

    def compare_all(self, other_hierarchy):
        """Compare lots of things about the structures to check they're identical"""

        sequence_diffs = self.compare_chain_sequences(other_hierarchy)
        count_diffs = self.compare_overall_counts(other_hierarchy)

        return sequence_diffs + count_diffs

    def compare_chain_sequences(self, other_hierarchy):
        """Extract the sequence for each chain and compare them"""

        # Extract sequences by chain
        my_seqs = dict([(c.id, c.as_sequence()) for c in self.hierarchy.chains()])
        other_seqs = dict([(c.id, c.as_sequence()) for c in other_hierarchy.chains()])

        differences = []

        # Compare the elements and counts thereof
        uniq_sequences, diff_sequences = compare_dictionaries(my_seqs, other_seqs)
        # Ignore the uniq chains - these are picked up by overall counts function
        for diff_seq in diff_sequences:
            print('DIFFERENT SEQUENCES: {!s}'.format(diff_seq))
            differences.append(('different sequence',)+diff_seq)

        return differences

    def compare_overall_counts(self, other_hierarchy):
        """Compare the counts of many variables in one structure to those in another structure"""

        # Extract counts
        my_counts = self.hierarchy.overall_counts()
        other_counts = other_hierarchy.overall_counts()

        differences = []

        # Counts of different variables
        count_vars = [  'n_models',
#                        'n_empty_models',
#                        'n_duplicate_model_ids',
                        'n_chains',
#                        'n_empty_chains',
#                        'n_duplicate_chain_ids',
#                        'n_chains_with_mix_of_proper_and_improper_alt_conf',
#                        'n_explicit_chain_breaks',
                        'n_residues',
                        'n_residue_groups',
#                        'n_empty_residue_groups',
                        'n_atoms',
#                        'n_empty_atom_groups',
#                        'n_duplicate_atom_labels',
#                        'n_anisou',
#                        'n_alt_conf',
#                        'n_alt_conf_none',
#                        'n_alt_conf_pure',
#                        'n_alt_conf_proper',
#                        'n_alt_conf_improper'
                                                ]

        for var_name in count_vars:
            if my_counts.__dict__[var_name] != other_counts.__dict__[var_name]:
                print('COUNT DIFFERENCE: `{!s}`'.format(var_name))
                differences.append(('count difference', var_name, my_counts.__dict__[var_name], other_counts.__dict__[var_name]))

        # Compare the chain ids and numbers of each chain
        uniq_chains, diff_chains = compare_dictionaries(my_counts.chain_ids, other_counts.chain_ids)
        for uniq_ch in uniq_chains:
            print('UNIQUE CHAIN: {!s}'.format(uniq_ch))
            differences.append(('unique chain', )+uniq_ch)
        for diff_ch in diff_chains:
            print('DIFFERENT CHAIN COUNT: {!s}'.format(diff_ch))
            differences.append(('different chain count',)+diff_ch)

        # Compare the resnames and counts thereof
        uniq_resnames, diff_resnames = compare_dictionaries(my_counts.resnames, other_counts.resnames)
        for uniq_res in uniq_resnames:
            print('UNIQUE RESIDUE: {!s}'.format(uniq_res))
            differences.append(('unique residue', )+uniq_res)
        for diff_res in diff_resnames:
            print('DIFFERENT RESIDUE COUNT: {!s}'.format(diff_res))
            differences.append(('different residue count',)+diff_res)

        # Compare the residues classes and counts thereof
        uniq_classes, diff_classes = compare_dictionaries(my_counts.resname_classes, other_counts.resname_classes)
        for uniq_cl in uniq_classes:
            print('UNIQUE RESIDUE CLASS: {!s}'.format(uniq_cl))
            differences.append(('unique residue class', )+uniq_cl)
        for diff_cl in diff_classes:
            print('DIFFERENT RESIDUE CLASS COUNT: {!s}'.format(diff_cl))
            differences.append(('different residue class count',)+diff_cl)

        # Compare the elements and counts thereof
        uniq_elems, diff_elems = compare_dictionaries(my_counts.element_charge_types, other_counts.element_charge_types)
        for uniq_el in uniq_elems:
            print('UNIQUE ELEMENT: {!s}'.format(uniq_el))
            differences.append(('unique element', )+uniq_el)
        for diff_el in diff_elems:
            print('DIFFERENT ELEMENT COUNT: {!s}'.format(diff_el))
            differences.append(('different element count',)+diff_el)

        return differences

class dataset_handler(object):
    def __init__(self, dataset_number, pdb_filename, mtz_filename, dataset_tag=None):
        """Create a dataset object to allow common functions to be applied easily to a number of datasets"""

        assert os.path.exists(pdb_filename), 'PDB file does not exist!'
        assert os.path.exists(mtz_filename), 'MTZ file does not exist!'

        # Store dataset number
        self.d_num = dataset_number
        # Store the name for the dataset
        if dataset_tag:
            self.d_tag = dataset_tag
        else:
            self.d_tag = 'D{:05d}'.format(self.d_num)

        # Store filenames
        self._pdb_file = os.path.abspath(pdb_filename)
        self._mtz_file = os.path.abspath(mtz_filename)
        # PDB Objects
        self._input = pdb_reader.input(source_info=None,lines=open(self._pdb_file,'r').read())
        self._hiera = self._input.construct_hierarchy()

        # Extract miller array
        self._fobs_miller = extract_miller_array_from_file(self._mtz_file, 'F,SIGF', log=open(os.devnull,'a'))
        self._scaled_fobs_miller = None

        # Initialise other variables
        self._fft_map = None
        self._unit_cell = None
        self._space_group = None
        self._basic_map = None
        self._rt_transform = None

        # Map Offsets for manual scaling
        self._raw_map_offset = None
        self._raw_map_scale = None

        # Map of the clusters in the dataset
        self.raw_cluster_hits = None
        self.clustered_hits = None

    def get_pdb_filename(self):
        return self._pdb_file
    def get_mtz_filename(self):
        return self._mtz_file
    def get_input(self):
        return self._input
    def get_hierarchy(self):
        return self._hiera
    def get_heavy_atom_sites(self):
        xray_structure = self.get_input().xray_structure_simple()
        return xray_structure.sites_cart().select(xray_structure.heavy_selection())
    def get_calpha_sites(self):
        xray_structure = self.get_input().xray_structure_simple()
        return xray_structure.sites_cart().select(xray_structure.backbone_selection(atom_names=['CA']))
    def get_backbone_sites(self):
        xray_structure = self.get_input().xray_structure_simple()
        return xray_structure.sites_cart().select(xray_structure.backbone_selection())
    def get_alignment_transform(self):
        return self._rt_transform
    def get_fobs_miller_array(self):
        return self._fobs_miller
    def get_scaled_fobs_miller_array(self):
        return self._scaled_fobs_miller
    def get_resolution(self):
        return self.get_fobs_miller_array().d_min()
    def get_map(self):
        return self._fft_map
    def get_map_handler(self):
        return self._basic_map

    def get_pickle_copy(self):
        """Get copy of self that can be pickled - some cctbx objects cannot be pickled..."""
        new = copy.copy(self)
        new._fft_map = None
        new._basic_map = None
        return new

    def get_structure_summary(self):
        return structure_summary(hierarchy=self._hiera)

    def compare_structure(self):
        """Return an object that can be used to compare structures"""
        return structure_comparer(hierarchy=self._hiera)

    def transform_from_reference(self, points):
        """Use alignment to map from reference frame to our frame"""
        return self.get_alignment_transform().rt().inverse() * points

    def transform_to_reference(self, points):
        """Use alignment to map to reference frame from our frame"""
        return self.get_alignment_transform().rt() * points

    def scale_fobs_to_reference(self, ref_miller):
        """Scale the observed data to another dataset"""
        self._scaled_fobs_miller = apply_simple_scaling(self.get_fobs_miller_array(), ref_miller)
        return self._scaled_fobs_miller

    def align_to_reference(self, ref_handler, sites='calpha', method='global'):
        """Calculate the rotation and translation needed to align one structure to another"""
        assert sites in ['calpha','backbone']
        assert method in ['global', 'local']
        # TODO ADD FLEXIBLE ALIGNMENT
        if method == 'global':
            if sites == 'calpha':
                my_sites = self.get_calpha_sites()
                ref_sites = ref_handler.get_calpha_sites()
            elif sites == 'backbone':
                my_sites = self.get_backbone_sites()
                ref_sites = ref_handler.get_backbone_sites()
            assert len(my_sites) == len(ref_sites)
            self._rt_transform = superpose.least_squares_fit(reference_sites=ref_sites, other_sites=my_sites)
        elif method == 'local':
            print 'USING FLEXIBLE ALIGNMENT'
            raise Exception('BUT NOT YET...')
        return self._rt_transform

    def create_fft_map(self, miller=None, map_type='2mFo-DFc', d_min=None):
        """Create an fft map to allow map values to be calculated"""
        if not miller:
            self._fft_map = get_fft_map_from_f_obs_and_structure(f_obs=self._fobs_miller, xray_structure=self.get_input().xray_structure_simple(), map_type=map_type, d_min=d_min)
        else:
            self._fft_map = get_fft_map_from_f_obs_and_structure(f_obs=miller, xray_structure=self.get_input().xray_structure_simple(), map_type=map_type, d_min=d_min)
        self._unit_cell = self._fft_map.unit_cell()
        self._space_group = self._fft_map.space_group()
        return self._fft_map

    def create_map_handler(self):
        """Create a map handler to allow for easy retrieval of map values"""
        if self._fft_map == None: self.create_fft_map()
        self._basic_map = maptbx.basic_map(maptbx.basic_map_unit_cell_flag(), self._fft_map.real_map(), self._fft_map.real_map().focus(),
                            self._unit_cell.orthogonalization_matrix(), maptbx.out_of_bounds_clamp(0).as_handle(), self._unit_cell)
        return self._basic_map

    def set_raw_map_offset_and_scale(self, offset, scale):
        """Set an offset and a scaling factor for all of the values in the map"""
        self._raw_map_offset = offset
        self._raw_map_scale = scale

    def get_cart_values(self, points, offset=None, scale=None):
        """Return a list of grid values shifted by `offset` and scaled by `scale` - (val+offset)*scale"""

        assert (self._raw_map_offset != None) and (self._raw_map_scale != None)

        # POTENTIAL TODO? - USE THIS AS THE METHOD TO EXTRACT MAP VALUES - THIS WOULD ALLOW THE MAP HANDLERS AND FFT TO BE DONE WHEN NEEDED!
        # Get unmodified map values
        orig_map_vals = self.get_map_handler().get_cart_values(points)
        # Shift and scale
        return (orig_map_vals+self._raw_map_offset)*self._raw_map_scale

class multi_dataset_analyser(object):
    def __init__(self, args):
#                    outdir='./pandda', datadir='./Processing', pdb_style='*/refine.pdb', mtz_style='*/refine.mtz',
#                    ref_pdb='./reference.pdb', ref_mtz='./reference.mtz', run_mol_subst=False,
#                    verbose=True, keep_maps_in_memory=False, maxmemory=25):
        """Class for the processing of datasets from a fragment soaking campaign"""

        self.args = args

        self.verbose = args.verbose
        self.keep_maps_in_memory = args.keep_maps_in_memory
        self.write_all_maps = False

        self._log = ''
        self._version = 0.1

        # ===============================================================================>
        # INPUT FILES STUFF
        # ===============================================================================>

        self.data_dirs = os.path.abspath(args.data_dirs)

        # ===============================================================================>
        # OUTPUT FILES STUFF
        # ===============================================================================>

        self.outdir = easy_directory(os.path.abspath(args.outdir))
        self.log_file = os.path.join(self.outdir, 'pandda.log')

        # Set up the output folders and files
        self._run_directory_setup()

        # ===============================================================================>
        # SETTINGS STUFF
        # ===============================================================================>

        self._map_type = '2mFo-DFc'
        self._res_factor = None
        self._cut_resolution = None
        self._border_padding = 10

        # Current and maximum sizes of the pandda in memory
        self._pandda_size = [('Zero',0)]
        self._max_pandda_size = args.max_memory*1024**3

        # ===============================================================================>
        # DATA AND MAPS STUFF
        # ===============================================================================>

        # File names
        self._input_files = []
        # Dataset Objects
        self._datasets = []
        # Reference Objects
        self._ref_dataset_index = None
        self._ref_dataset = None
        self._ref_grid = None
        self._ref_map = None
        # Distance to displace the reference model by so that its bounding box sits on the origin
        self._ref_origin_shift = None

        # Map Arrays
        self._maps = []
        self._z_maps = []
        self._mod_z_maps = []

        # Map Statistics
        self._mean_map = None
        self._stds_map = None
        self._skew_map = None
        self._kurt_map = None
        self._bimo_map = None

        # ===============================================================================>
        # ANALYSIS OBJECTS
        # ===============================================================================>

        # Masks for various stages of the pipeline
        self._dataset_masks = mask_collection()

        # Map and Dataset statistics
        self._map_summary = data_collection()
        self._dataset_observations = data_collection()

        # Summaries of detected points
        self._cluster_summary = data_collection()

        # Structure summary object
        self._structure_observations = multiple_data_collection()

        # Average Structure (and masks)
        self._average_calpha_sites = None
        self._residue_deviation_masks = None

        self.timings = data_collection()

        # ===============================================================================>
        # PICKLE STUFF
        # ===============================================================================>

        # Set up the pickled object filenames
        self._run_pickle_setup()

        if os.path.exists(self.pickle_handler.get_file('dataset_objs')):
            self._new_pandda = False
        else:
            self._new_pandda = True

    def run_pandda_init(self):
        """Set up the pandda"""

        # Print logo to log
        self.log(PANDDA_TEXT.format(self._version), True)

        # Store start time and print to log
        self._init_time = time.time()
        self.log('Analysis Started: {!s}'.format(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime(self._init_time))), True)

        # ===============================================================================>
        # REPOPULATE PANDDA FROM PREVIOUS RUNS
        # ===============================================================================>

        if os.path.exists(self.output_handler.get_file('reference_structure')) and os.path.exists(self.output_handler.get_file('reference_dataset')):
            self.log('===================================>>>', True)
            self.log('Loading Reference Dataset', True)
            self.load_reference_dataset(ref_pdb=self.output_handler.get_file('reference_structure'), ref_mtz=self.output_handler.get_file('reference_dataset'))

        # Load any objects from previous runs
        self.load_pickled_objects()

        # Get the size of the empty pandda
        self.update_pandda_size(tag='Initialised Pandda')

    def _run_pickle_setup(self):
        """Initialise all of the pickle filenames"""

        # Pickle Handler
        self.pickle_handler = output_file_object(rootdir=self.output_handler.get_dir('pickle'))
        # Pickled Reference Objects
        self.pickle_handler.add_file(file_name='reference_grid.pickle', file_tag='reference_grid')
        self.pickle_handler.add_file(file_name='global_mask.pickle', file_tag='global_mask')
        self.pickle_handler.add_file(file_name='local_mask.pickle', file_tag='local_mask')
        # Pickled Datasets
        self.pickle_handler.add_file(file_name='dataset_masks.pickle', file_tag='dataset_masks')
        self.pickle_handler.add_file(file_name='dataset_objs.pickle', file_tag='dataset_objs')
        # Pickled Information
        self.pickle_handler.add_file(file_name='dataset_observations.pickle', file_tag='dataset_observations')
        self.pickle_handler.add_file(file_name='map_summary.pickle', file_tag='map_summary')
        # Pickled Maps
        self.pickle_handler.add_file(file_name='map_values.pickle', file_tag='map_values')
        self.pickle_handler.add_file(file_name='z_map_values.pickle', file_tag='z_map_values')
        self.pickle_handler.add_file(file_name='mod_z_map_values.pickle', file_tag='mod_z_map_values')
        # Pickled Stats
        self.pickle_handler.add_file(file_name='mean_map.pickle', file_tag='mean_map')
        self.pickle_handler.add_file(file_name='stds_map.pickle', file_tag='stds_map')
        self.pickle_handler.add_file(file_name='skew_map.pickle', file_tag='skew_map')
        self.pickle_handler.add_file(file_name='kurt_map.pickle', file_tag='kurt_map')
        self.pickle_handler.add_file(file_name='bimo_map.pickle', file_tag='bimo_map')

    def _run_directory_setup(self):
        """Initialise the pandda directory system"""

        # Create a file and directory organiser
        self.output_handler = output_file_object(rootdir=self.outdir)
        # Somewhere to store all of the input files
        self.output_handler.add_dir(dir_name='input_files', dir_tag='input_files', top_dir_tag='root')
        # Somewhere to store all of the aligned structures
        self.output_handler.add_dir(dir_name='aligned_structures', dir_tag='aligned_structures', top_dir_tag='root')
        # Somewhere to store the interesting datasets
        self.output_handler.add_dir(dir_name='initial_hits', dir_tag='initial_hits', top_dir_tag='root')
        # Somewhere to store the manual analyses
        self.output_handler.add_dir(dir_name='manual_analyses', dir_tag='manual_analyses', top_dir_tag='root')
        # Somewhere to store the analysis summaries
        self.output_handler.add_dir(dir_name='summaries', dir_tag='output_summaries', top_dir_tag='root')
        self.output_handler.add_file(file_name='success_summary_table.html', file_tag='summary_table', dir_tag='output_summaries')
        # Somewhere to store the pickled objects
        self.output_handler.add_dir(dir_name='pickled_panddas', dir_tag='pickle', top_dir_tag='root')

        # Reference Structure and Dataset
        self.output_handler.add_dir(dir_name='reference', dir_tag='reference', top_dir_tag='root')
        self.output_handler.add_file(file_name='reference.pdb', file_tag='reference_structure', dir_tag='reference')
        self.output_handler.add_file(file_name='reference.mtz', file_tag='reference_dataset', dir_tag='reference')
        self.output_handler.add_file(file_name='reference.aligned.pdb', file_tag='reference_structure_on_origin', dir_tag='reference')

    def load_pickled_objects(self):
        """Loads any pickled objects it finds"""

        self.log('===================================>>>', True)
        self.log('Looking for Pickled Files in Input Directory: {!s}'.format(os.path.relpath(self.pickle_handler.get_dir('root'))), True)
        # Load Reference Grid
        if os.path.exists(self.pickle_handler.get_file('reference_grid')):
            self.set_reference_grid(self.unpickle(self.pickle_handler.get_file('reference_grid')))
        if self.reference_grid() is not None:
            if os.path.exists(self.pickle_handler.get_file('global_mask')) and (self.reference_grid().global_mask() is None):
                self.reference_grid().set_global_mask(self.unpickle(self.pickle_handler.get_file('global_mask')))
            if os.path.exists(self.pickle_handler.get_file('local_mask')) and (self.reference_grid().local_mask() is None):
                self.reference_grid().set_local_mask(self.unpickle(self.pickle_handler.get_file('local_mask')))

        # Load the datasets
        if os.path.exists(self.pickle_handler.get_file('dataset_objs')):
            self._datasets = self.unpickle(self.pickle_handler.get_file('dataset_objs'))
            self.update_pandda_size(tag='After Unpickling Dataset Objects')

        # Load the dataset summary object
#        if os.path.exists(self.pickle_handler.get_file('dataset_masks')):
#            self._dataset_masks = self.unpickle(self.pickle_handler.get_file('dataset_masks'))
#        if os.path.exists(self.pickle_handler.get_file('dataset_observations')):
#            self._dataset_observations = self.unpickle(self.pickle_handler.get_file('dataset_observations'))
#        if os.path.exists(self.pickle_handler.get_file('map_summary')):
#            self._map_summary = self.unpickle(self.pickle_handler.get_file('map_summary'))

        # Load Statistical Maps
        if os.path.exists(self.pickle_handler.get_file('mean_map')):
            self._mean_map = self.unpickle(self.pickle_handler.get_file('mean_map'))
        if os.path.exists(self.pickle_handler.get_file('stds_map')):
            self._stds_map = self.unpickle(self.pickle_handler.get_file('stds_map'))
        if os.path.exists(self.pickle_handler.get_file('skew_map')):
            self._skew_map = self.unpickle(self.pickle_handler.get_file('skew_map'))
        if os.path.exists(self.pickle_handler.get_file('kurt_map')):
            self._kurt_map = self.unpickle(self.pickle_handler.get_file('kurt_map'))
        if os.path.exists(self.pickle_handler.get_file('bimo_map')):
            self._bimo_map = self.unpickle(self.pickle_handler.get_file('bimo_map'))
            self.update_pandda_size(tag='After Unpickling Statistical Maps')

    def pickle_the_pandda(self):
        """Pickles it's major components for quick loading..."""

        self.log('===================================>>>', True)
        self.log('Pickling the PanDDA', True)

        self.log('===================================>>>')
        self.log('Pickling Reference Grid')
        if self._ref_grid is not None:
            self.pickle(pickle_file=self.pickle_handler.get_file('reference_grid'), pickle_object=self.reference_grid(), force=False)
            if self._ref_grid.global_mask() is not None:
                self.pickle(pickle_file=self.pickle_handler.get_file('global_mask'), pickle_object=self.reference_grid().global_mask(), force=False)
            if self._ref_grid.local_mask() is not None:
                self.pickle(pickle_file=self.pickle_handler.get_file('local_mask'), pickle_object=self.reference_grid().local_mask(), force=False)

        self.log('===================================>>>')
        self.log('Pickling Datasets')
        if self._datasets:
            self.pickle(pickle_file=self.pickle_handler.get_file('dataset_objs'), pickle_object=[d.get_pickle_copy() for d in self.get_all_datasets()])

#        self.log('===================================>>>')
#        self.log('Pickling Dataset Variation Summary')
#        if self._dataset_masks:
#            self.pickle(pickle_file=self.pickle_handler.get_file('dataset_masks'), pickle_object=self._dataset_masks)
#        if self._dataset_observations:
#            self.pickle(pickle_file=self.pickle_handler.get_file('dataset_observations'), pickle_object=self.get_dataset_observations())
#        if self._map_summary:
#            self.pickle(pickle_file=self.pickle_handler.get_file('map_summary'), pickle_object=self.get_map_summary())

        self.log('===================================>>>')
        self.log('Pickling Statistical Maps')
        if self._mean_map is not None:
            self.pickle(pickle_file=self.pickle_handler.get_file('mean_map'), pickle_object=self.get_mean_map())
        if self._stds_map is not None:
            self.pickle(pickle_file=self.pickle_handler.get_file('stds_map'), pickle_object=self.get_stds_map())
        if self._skew_map is not None:
            self.pickle(pickle_file=self.pickle_handler.get_file('skew_map'), pickle_object=self.get_skew_map())
        if self._kurt_map is not None:
            self.pickle(pickle_file=self.pickle_handler.get_file('kurt_map'), pickle_object=self.get_kurt_map())
        if self._bimo_map is not None:
            self.pickle(pickle_file=self.pickle_handler.get_file('bimo_map'), pickle_object=self.get_bimo_map())

    def exit(self):

        self.log('===================================>>>', True)
        self.log('...FINISHED!...', True)
        self.log('===================================>>>', True)
        self._finish_time = time.time()
        self.log('Runtime: {!s}'.format(time.strftime("%H hours:%M minutes:%S seconds", time.gmtime(self._finish_time - self._init_time))))
        raise SystemExit('PANDDA EXECUTED/TRANQUILISED.')

    def log(self, message, show=False):
        """Log message to file, and mirror to stdout if verbose or force_print"""
        if not isinstance(message, str):    message = str(message)
        # Print to stdout
        if self.verbose or show:            print(message)
        # Store in internal string
        self._log = self._log + message + '\n'
        # Write to file
        open(self.log_file, 'a').write(message+'\n')

    def print_log(self):
        print(self._log)

    def get_max_pandda_size(self):
        return self._max_pandda_size
    def set_max_pandda_size(self, max_bytes):
        self._max_pandda_size = max_bytes
    def get_pandda_size(self):
        """Returns the history of the memory consumption of the PANDDA"""
        return self._pandda_size
    def update_pandda_size(self, tag):
        pandda_size = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1000
        assert pandda_size < self.get_max_pandda_size(), 'PANDDA HAS EXCEEDED THE MAXIMUM AMOUNT OF ALLOWED MEMORY'
        self._pandda_size.append((tag, pandda_size))
        from humanize import naturalsize
        self.log(tag+': '+naturalsize(pandda_size, binary=True))

    def is_new_pandda(self):
        """Is this the first time the program has been run?"""
        return self._new_pandda

    def set_map_type(self, map_type):
        assert map_type in ['2mFo-DFc','mFo-DFc','mFo'], 'Map type not recognised: {!s}'.format(map_type)
        self._map_type = map_type
    def get_map_type(self):
        return self._map_type
    def set_map_scaling(self, scaling):
        assert scaling in ['none','sigma','volume','absolute'], 'Map scaling not recognised: {!s}'.format(scaling)
        self._map_scaling = scaling
    def get_map_scaling(self):
        return self._map_scaling

    def set_res_factor(self, res_factor):
        self._res_factor = res_factor
    def get_res_factor(self):
        return self._res_factor
    def set_cut_resolution(self, d_min):
        self._cut_resolution = d_min
    def get_cut_resolution(self):
        return self._cut_resolution

    def set_border_padding(self, border):
        self._border_padding = border
    def get_border_padding(self):
        return self._border_padding

    def set_reference_origin_shift(self, origin_shift):
        self._ref_origin_shift = origin_shift
    def get_reference_origin_shift(self):
        return self._ref_origin_shift

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

    def get_dataset_masks(self):
        """Get the mask object for the datasets"""
        return self._dataset_masks

    def get_added_files(self):
        """Get all of the files that were added"""
        return self._input_files

    def get_dataset(self, d_num=None, d_tag=None):
        """Get a dataset by dataset number or tag"""
        assert d_num or d_tag, 'NO NUM or TAG given!'
        assert not (d_num and d_tag), 'BOTH NUM and TAG given!'
        if d_num:
            matching = [d for d in self.get_all_datasets() if d.d_num == d_num]
        else:
            matching = [d for d in self.get_all_datasets() if d.d_tag == d_tag]
        assert len(matching) == 1, 'MORE THAN ONE MATCHING DATASET FOR DNUM'
        return matching[0]

    def get_all_datasets(self):
        """Return all datasets added"""
        return self._datasets
    def get_masked_datasets(self, mask_name, invert=False):
        """Use a custom mask to select datasets"""
        return self.get_dataset_masks().mask(mask_name=mask_name, input_list=self._datasets, invert=invert)

    def get_number_of_datasets(self, mask_name=None, invert=False):
        """Returns the number of datasets"""
        if mask_name:   return len(self.get_dataset_masks().mask(mask_name=mask_name, input_list=self._datasets, invert=invert))
        else:           return len(self._datasets)

    def get_dataset_observations(self):
        """Distributions of different variables across the datasets"""
        return self._dataset_observations
    def get_map_summary(self):
        """Distributions of different map variables across the datasets"""
        return self._map_summary
    def get_cluster_summary(self):
        """Summaries of difference clusters across the datasets"""
        return self._cluster_summary

    def get_structure_observations(self):
        """Distributions of different variables across the datasets"""
        return self._structure_observations
    def get_calpha_average_sites(self):
        return self._average_calpha_sites
    def get_calpha_deviation_masks(self):
        return self._residue_deviation_masks

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
    def get_bimo_map(self):
        """Returns the bimodality map across the datasets"""
        return self._bimo_map

    def initialise_masks(self):
        """Add blank masks to the mask objects, based on how many datasets have been loaded"""

        self.log('===================================>>>', True)
        self.log('Initialising Dataset Masks.', True)

        # Store lists of all of the possible masks
        self._structure_mask_names = STRUCTURE_MASK_NAMES
        self._crystal_mask_names   = CRYSTAL_MASK_NAMES
        self._reject_mask_names    = REJECT_MASK_NAMES
        self._flag_mask_names      = FLAG_MASK_NAMES
        self._all_mask_names = self._structure_mask_names + self._crystal_mask_names + self._reject_mask_names + self._flag_mask_names

        # Set the dataset mask lengths and ids (dataset tags)
        self.get_dataset_masks().set_mask_length(mask_length=self.get_number_of_datasets())
        self.get_dataset_masks().set_entry_ids(entry_ids=[d.d_tag for d in self.get_all_datasets()])
        # Initialise standard blank masks
        for mask_name in self._all_mask_names:
            self.get_dataset_masks().add_mask(mask_name=mask_name, mask=[False]*self.get_number_of_datasets())

        # Set the dataset observations lengths
        self.get_dataset_observations().set_data_length(self.get_number_of_datasets())
        self.get_dataset_observations().set_data_ids([d.d_tag for d in self.get_all_datasets()])

    def load_reference_dataset(self, ref_pdb, ref_mtz):
        """Set the reference dataset, to which all other datasets will be aligned and scaled"""

        self.log('===================================>>>', True)
        self.log('Loading Reference Dataset: {!s}'.format(ref_mtz), True)

        self._ref_dataset = dataset_handler(dataset_number=-1, pdb_filename=ref_pdb, mtz_filename=ref_mtz, dataset_tag='reference')

        if not os.path.exists(self.output_handler.get_file('reference_structure')):
            os.symlink(ref_pdb, self.output_handler.get_file('reference_structure'))
        if not os.path.exists(self.output_handler.get_file('reference_dataset')):
            os.symlink(ref_mtz, self.output_handler.get_file('reference_dataset'))

        # Calculate the shift required to move the reference structure into the positive quadrant
        self.set_reference_origin_shift(tuple(flex.double(3, self.get_border_padding()) - flex.double(self.reference_dataset().get_input().atoms().extract_xyz().min())))
        self.log('Origin Shift for reference structure: {!s}'.format(tuple([round(s,3) for s in self.get_reference_origin_shift()])))
        # Shift the reference structure by this amount so that it is aligned with the reference grid
        ref_hierarchy = self.reference_dataset().get_hierarchy()
        ref_hierarchy.atoms().set_xyz(ref_hierarchy.atoms().extract_xyz() + self.get_reference_origin_shift())

        if not os.path.exists(self.output_handler.get_file('reference_structure_on_origin')):
            ref_hierarchy.write_pdb_file(self.output_handler.get_file('reference_structure_on_origin'))

#             # Write out the translated structure
#            translated_struc = self.reference_dataset().get_hierarchy().deep_copy()
#            translated_struc.atoms().set_xyz(translated_struc.atoms().extract_xyz() + self.get_reference_origin_shift())
#            # Write out next to the original files
#            translated_struc.write_pdb_file(self.output_handler.get_file('reference_structure_on_origin'))

        self._initialise_structure_analysis()

    def _initialise_structure_analysis(self):
        """Analyse the reference structure and pull out structural information/summaries"""

        self.log('===================================>>>', True)
        self.log('Extracting Reference Structure Information', True)

        # Extract the reference hierarchy
        ref_hierarchy = self.reference_dataset().get_hierarchy()
        # Build a list of the residue labels
        self.get_structure_observations().set_collection_ids([(rg.parent().id, rg.resid(), rg.unique_resnames()[0]) for rg in ref_hierarchy.residue_groups()])

        # Build a list of residue types
        residue_types = [pdb_reader.common_residue_names_get_class(rg.unique_resnames()[0]) for rg in ref_hierarchy.residue_groups()]
        self.get_structure_observations().add_collection_info(info_name='residue types', info_values=residue_types)

        # TODO Add residue labels (chain, resnum, etc.)

    def load_reference_map_handler(self):
        """Load the reference map handlers"""

        self.log('===================================>>>', True)
        self.log('Loading Reference Map Handler', True)

        self._ref_dataset.create_fft_map(map_type=self.get_map_type())
        if self.get_map_scaling() == 'none':
            pass
        elif self.get_map_scaling() == 'sigma':
            self._ref_dataset.get_map().apply_sigma_scaling()
        elif self.get_map_scaling() == 'volume':
            self._ref_dataset.get_map().apply_volume_scaling()
        self._ref_dataset.create_map_handler()

    def extract_reference_map_values(self):
        """Read in map for the reference dataset"""

        self.log('===================================>>>', True)
        self.log('Loading Reference Map Values', True)

        # DON'T NEED TO SHIFT THESE COORDINATES BECAUSE THE MAPS ARE CREATED FROM THE MODEL, WHICH HAS BEEN TRANSLATED
        ref_map_vals = self.reference_dataset().get_map_handler().get_cart_values(self.reference_grid().cart_points())
        ref_map_file = self.output_handler.get_file('reference_dataset').replace('.mtz','.sampled.ccp4')
        self.write_array_to_map(output_file=ref_map_file, map_data=ref_map_vals)
        self._ref_map = ref_map_vals

    def create_reference_grid(self, res_factor, expand_to_origin, buffer=0):
        """Create a grid over the reference protein"""

        if res_factor:
            self.set_res_factor(res_factor)
        grid_spacing = calculate_sampling_distance(self.get_cut_resolution(), self.get_res_factor())
        sites_cart = self.reference_dataset().get_input().atoms().extract_xyz()

        assert (flex.vec3_double(sorted(self.reference_dataset().get_input().atoms().extract_xyz())) -
                flex.vec3_double(sorted(self.reference_dataset().get_hierarchy().atoms().extract_xyz()))).dot().norm() == 0.0, 'EH? Coordinates should be the same?'

        self._ref_grid = grid_handler(verbose=self.verbose)
        self._ref_grid.set_grid_spacing(spacing=grid_spacing)
        self._ref_grid.set_cart_extent(cart_min=tuple([s-buffer for s in sites_cart.min()]), cart_max=tuple([s+buffer for s in sites_cart.max()]))
        self._ref_grid.create_cartesian_grid(expand_to_origin=expand_to_origin)
        self._ref_grid.show_summary()

    def mask_resampled_reference_grid(self):
        """Using the local and global masks, mask the resampled grid points"""

        self.log('===================================>>>', True)
        self.log('Masking Resampled Grid (Global and Local Mask)', True)

        # Create a grid mask around each point
        self.log('===================================>>>', True)
        self.log('Grid Summary: ')
        self.reference_grid().show_summary()
        self.reference_grid().local_mask().show_summary()
        self.reference_grid().global_mask().show_summary()

        # Get the grid points that are not masked, and points in the buffer
        grid_indexer = self.reference_grid().grid_indexer()
        resampled_points = self.reference_grid().resampled_grid_points()
        buffer_mask_binary = self.reference_grid().buffer_mask_binary()
        global_mask_binary = self.reference_grid().global_mask().total_mask_binary()

        self.log('===================================>>>')
        self.log('Resampled Grid Size (3D): {!s}'.format(self.reference_grid().resampled_grid_size()))
        self.log('Resampled Grid Size (1D): {!s}'.format(len(resampled_points)))
        self.log('===================================>>>')

        # Remove points using the protein mask, and the mask around the edge of the grid
        masked_resampled_points_1 = resampled_points
        self.log('Filtering with Buffer Mask... (Edge of Cell)')
        masked_resampled_points_2 = [p for p in masked_resampled_points_1 if buffer_mask_binary[grid_indexer(p)] == 0]
        self.log('Filtered Points: {!s}'.format(len(masked_resampled_points_2)))
        self.log('Filtering with Global Mask... (Protein)')
        masked_resampled_points_3 = [p for p in masked_resampled_points_2 if global_mask_binary[grid_indexer(p)] == 1]
        self.log('Filtered Points: {!s}'.format(len(masked_resampled_points_3)))
        masked_resampled_points = masked_resampled_points_3

        # Store points in the reference grid object
        self.reference_grid().set_masked_grid_points(masked_resampled_points)

        # Calculate the resampled mask on the original grid
        full_binary_mask = numpy.zeros(self.reference_grid().grid_size_1d(), int)
        [full_binary_mask.put(grid_indexer(gp), 1) for gp in masked_resampled_points]
        # Downsample to the resampled grid
        resampled_binary_mask = [full_binary_mask[grid_indexer(gp)] for gp in resampled_points]
        # Write masked map
        mask_map_file = self.output_handler.get_file('reference_dataset').replace('.mtz','.totalmask.ccp4')
        self.write_array_to_map(output_file = mask_map_file,
                                map_data = flex.double(resampled_binary_mask),
                                grid_size = self.reference_grid().resampled_grid_size(),
                                grid_spacing = self.reference_grid().resampled_grid_spacing())

    def build_input_list(self):
        """Builds a list of input files from the command line arguments passed"""

        self.log('===================================>>>', True)
        self.log('Building List of Datasets')

        dir_style = self.args.data_dirs.strip('./')
        pdb_style = self.args.pdb_style.lstrip('/')
        if self.args.mtz_style:
            mtz_style = self.args.mtz_style.lstrip('/')
        else:
            mtz_style = pdb_style.replace('.pdb','.mtz')

        # Datasets that are already added
        already_added  = [(d.get_pdb_filename(), d.get_mtz_filename()) for d in self.get_all_datasets()]
        new_files = []

        for dir in sorted(glob.glob(self.data_dirs)):
            pdb_files = glob.glob(os.path.join(dir, pdb_style))
            mtz_files = glob.glob(os.path.join(dir, mtz_style))
            if not (pdb_files and mtz_files):
                print('EMPTY DIRECTORY: {!s}'.format(dir))
            elif not pdb_files:
                print('NO PDB IN DIRECTORY: {!s}'.format(dir))
            elif not mtz_files:
                print('NO MTZ IN DIRECTORY: {!s}'.format(dir))
            else:
                assert len(pdb_files) == 1, 'More than one matching PDB file found: {!s}'.format(os.path.join(dir, pdb_style))
                assert len(mtz_files) == 1, 'More than one matching MTZ file found: {!s}'.format(os.path.join(dir, mtz_style))

                new_pdb = pdb_files[0]
                new_mtz = mtz_files[0]
                dataset_tag = [None]

                if (os.path.abspath(new_pdb), os.path.abspath(new_mtz)) in already_added:
                    self.log('ALREADY ADDED (SKIPPING): {!s}'.format(dir))
                    continue

                # Do regex matching on the file pairs
                if '*' in pdb_style:
                    pdb_regex = pdb_style.replace('*', '(.*)')
                    pdb_tag = re.findall(pdb_regex, os.path.basename(pdb_files[0]))
                    assert len(pdb_tag) == 1, 'MORE THAN ONE PDB TAG FOUND'
                else: pdb_regex = pdb_tag = None

                if '*' in mtz_style:
                    mtz_regex = mtz_style.replace('*', '(.*)')
                    mtz_tag = re.findall(mtz_regex, os.path.basename(mtz_files[0]))
                    assert len(mtz_tag) == 1, 'MORE THAN ONE MTZ TAG FOUND'
                else: mtz_regex = mtz_tag = None

                if '*' in dir_style:
                    dir_regex = dir_style.replace('*', '(.*)')
                    dir_tag = re.findall(dir_regex, os.path.dirname(pdb_files[0]))
                    assert len(dir_tag) == 1, 'MORE THAN ONE DIR TAG FOUND'
                else: dir_regex = dir_tag = None

                if pdb_tag and mtz_tag: assert pdb_tag == mtz_tag, 'PDB-MTZ TAGS ARE NOT IDENTICAL'
                if dir_tag and pdb_tag: assert dir_tag == pdb_tag, 'DIR-PDB TAGS ARE NOT IDENTICAL'
                if dir_tag and mtz_tag: assert dir_tag == mtz_tag, 'DIR-MTZ TAGS ARE NOT IDENTICAL'

                if   dir_tag: dataset_tag = dir_tag
                elif pdb_tag: dataset_tag = pdb_tag
                elif mtz_tag: dataset_tag = mtz_tag

                new_files.append(pdb_files+mtz_files+dataset_tag)

        return new_files

    def add_new_files(self, input_files):
        """Add (pdb, mtz) file pairs to the datasets to be processed"""

        self.log('===================================>>>', True)
        self._input_files = input_files
        self.log('{!s} Datasets Added'.format(len(input_files)), True)

    def load_new_datasets(self):
        """Read in maps for the input datasets"""

        def map_func(arg_dict):
            return dataset_handler(**arg_dict)

        # TODO Change this so datasets can be added dynamically
        assert self._datasets == []

        # Generate arg_list for loading
        arg_list = [{'dataset_number':d_num, 'pdb_filename':pdb, 'mtz_filename':mtz, 'dataset_tag':tag} for d_num, (pdb, mtz, tag) in enumerate(self.get_added_files())]

        start = time.time()
        self.log('===================================>>>', True)
        print 'Loading Datasets... (using {!s} cores)'.format(self.args.cpus)
        loaded_datasets = easy_mp.pool_map(fixed_func=map_func, args=arg_list, processes=self.args.cpus)
        finish = time.time()
        print('-> Time Taken: {!s} seconds'.format(int(finish-start)))

        self._datasets.extend(loaded_datasets)
        self.log('{!s} Datasets Loaded.          '.format(len(loaded_datasets), True))

    def select_reference_dataset(self):
        """Select dataset to act as the reference - scaling, aligning etc"""

        r_frees = [d.get_input().get_r_rfree_sigma().r_free for d in self.get_all_datasets()]
        self._ref_dataset_index = r_frees.index(min(r_frees))
        reference = self.get_all_datasets()[self._ref_dataset_index]
        assert reference.d_num == self._ref_dataset_index, 'INDEX DOES NOT MATCH DNUM'
        self.log('===================================>>>', True)
        self.log('Reference Selected: Dataset {!s}'.format(self._ref_dataset_index+1), True)
        return reference.get_pdb_filename(), reference.get_mtz_filename()

    def scale_datasets(self):
        """Iterate through the datasets, and scale the reflections to the reference"""

        self.log('===================================>>>', True)
        for d_handler in self.get_masked_datasets(mask_name='rejected - total', invert=True):
            print '\rScaling Dataset {!s}'.format(d_handler.d_tag),; sys.stdout.flush()
            # Scale new data to the reference dataset
            d_handler.scale_fobs_to_reference(ref_miller=self.reference_dataset().get_fobs_miller_array())
        self.log('\rDatasets Scaled.          ', True)

    def align_datasets(self, sites='calpha', method='global'):
        """Align each structure the reference structure"""

        assert sites in ['calpha','backbone']

        self.log('===================================>>>', True)
        for d_handler in self.get_masked_datasets(mask_name='rejected - total', invert=True):
            print '\rAligning Dataset {!s}'.format(d_handler.d_tag),; sys.stdout.flush()
#            try:
            # Align to reference structure to get mapping transform
            d_handler.align_to_reference(ref_handler=self.reference_dataset(), sites=sites, method=method)
#            except:
#                self.log('> Alignment Error - skipping', True)
#                continue
            # Write out the aligned structure
            aligned_struc = d_handler.get_hierarchy().deep_copy()
            aligned_struc.atoms().set_xyz(d_handler.transform_to_reference(d_handler.get_hierarchy().atoms().extract_xyz()))
            aligned_file = d_handler.get_pdb_filename().replace('.pdb','.aligned.pdb')
            # Write out next to the original files
            aligned_struc.write_pdb_file(file_name=aligned_file)
            # Also write out into the aligned_structures folder
            aligned_struc.write_pdb_file(file_name=os.path.join(self.output_handler.get_dir('aligned_structures'), 'structure_{!s}.pdb'.format(d_handler.d_tag)))
            # Check that the objects have copied properly so that the original coordinates haven't been modified by .set_xyz()
            alignment_rmsd = d_handler.get_hierarchy().atoms().extract_xyz().rms_difference(aligned_struc.atoms().extract_xyz())
            # Check that the deep_copy has worked
            if alignment_rmsd == 0.0:
                # This must be the reference! Set the dataset number if it's not already set
                if self._ref_dataset_index == None:
                    self._ref_dataset_index = d_handler.d_num
                # Raise error if the reference has already been set and it's not this dataset
                elif self._ref_dataset_index != d_handler.d_num:
                    raise Exception('ALIGNED OBJECT EQUAL TO UNALIGNED OBJECT - THIS IS MOST UNLIKELY')
        self.log('\rDatasets Aligned.          ', True)

    def collect_dataset_variation_data(self):
        """Go through all of the datasets and collect lots of different characteristics of the datasets for identifying odd datasets"""

        self.log('===================================>>>', True)
        self.log('Collecting Dataset/Crystal Variation Data', True)
        self.log('===================================>>>')

        self.log('Extracting Resolutions')
        self.get_dataset_observations().add_data(data_name='resolution', data_values=[d.get_fobs_miller_array().d_max_min() for d in self.get_all_datasets()])
        self.log('Extracting Unit Cell Sizes')
        self.get_dataset_observations().add_data(data_name='cell_params', data_values=[d.get_fobs_miller_array().unit_cell().parameters() for d in self.get_all_datasets()])
        self.log('Extracting Unit Cell Volume')
        self.get_dataset_observations().add_data(data_name='cell_volume', data_values=[d.get_fobs_miller_array().unit_cell().volume() for d in self.get_all_datasets()])
        self.log('Extracting R-work, R-free')
        self.get_dataset_observations().add_data(data_name='rfree_rwork', data_values=[(d.get_input().get_r_rfree_sigma().r_free, d.get_input().get_r_rfree_sigma().r_work) for d in self.get_all_datasets()])
        self.log('Calculating Variation in Correlations between Diffraction Data')
        self.get_dataset_observations().add_data(data_name='correlation_to_unscaled_data', data_values=[d.get_scaled_fobs_miller_array().correlation(d.get_fobs_miller_array()).coefficient() if d.get_scaled_fobs_miller_array() else -1 for d in self.get_all_datasets()])

    def filter_datasets_1(self):
        """Filter out the datasets which contain different protein models (i.e. protein length, sequence, etc)"""

        self.log('===================================>>>', True)
        self.log('Filtering Datasets (Non-identical structures)). Potential Classes:', True)
        for failure_class in self._structure_mask_names:
            self.log('\t{!s}'.format(failure_class), True)
        self.log('===================================>>>', True)

        ref_counts, ref_dict = self.reference_dataset().get_structure_summary()

        # Check that the same protein structure is present in each dataset
        for d_handler in self.get_masked_datasets(mask_name='rejected - total', invert=True):

            my_counts, my_dict = d_handler.get_structure_summary()

            print '\rFiltering Dataset {!s}'.format(d_handler.d_tag),; sys.stdout.flush()
            # Check the number of chains
            if my_counts.n_chains != ref_counts.n_chains:
                self.log('\rRejecting Dataset: {!s}'.format(d_handler.d_tag))
                self.log('Different Number of Chains')
                self.log('Reference: {!s}, {!s}: {!s}'.format(ref_counts.n_chains, d_handler.d_tag, my_counts.n_chains))
                self.log('===================================>>>')
                self.get_dataset_masks().set_mask_value(mask_name='bad structure - chain counts', entry_id=d_handler.d_tag, value=True)
            # Check the ids of the chains
            elif my_counts.chain_ids != ref_counts.chain_ids:
                self.log('\rRejecting Dataset: {!s}'.format(d_handler.d_tag))
                self.log('Different Chain IDs')
                self.log('Reference: {!s}, {!s}: {!s}'.format(ref_counts.chain_ids, d_handler.d_tag, my_counts.chain_ids))
                print('===================================>>>')
                self.get_dataset_masks().set_mask_value(mask_name='bad structure - chain ids', entry_id=d_handler.d_tag, value=True)
            # Check the sequences of the chains
            elif my_dict['chain_sequences'] != ref_dict['chain_sequences']:
                self.log('\rRejecting Dataset: {!s}'.format(d_handler.d_tag))
                self.log('Different Sequences')
                for chain_id in ref_dict['chain_sequences'].keys():
                    self.log('Chain {!s}: Reference - {!s}'.format(chain_id, ref_dict['chain_sequences'][chain_id]))
                    self.log('Chain {!s}: {:>9s} - {!s}'.format(chain_id, my_dict['chain_sequences'][chain_id]))
                print('===================================>>>')
                self.get_dataset_masks().set_mask_value(mask_name='bad structure - chain sequences', entry_id=d_handler.d_tag, value=True)
            # Check the number of residues - TODO not sure it can ever get here... remove?
            elif my_counts.n_residues != ref_counts.n_residues:
                self.log('\rRejecting Dataset: {!s}'.format(d_handler.d_tag))
                self.log('Different Number of Residues')
                self.log('Reference: {!s}, {!s}: {!s}'.format(ref_counts.n_residues, d_handler.d_tag, my_counts.n_residues))
                print('===================================>>>')
                self.get_dataset_masks().set_mask_value(mask_name='bad structure - residue counts', entry_id=d_handler.d_tag, value=True)
            # Check the number of atoms
            elif my_counts.n_atoms != ref_counts.n_atoms:
                self.log('\rRejecting Dataset: {!s}'.format(d_handler.d_tag))
                self.log('Different Number of Atoms')
                self.log('Reference: {!s}, {!s}: {!s}'.format(ref_counts.n_atoms, d_handler.d_tag, my_counts.n_atoms))
                print('===================================>>>')
                self.get_dataset_masks().set_mask_value(mask_name='bad structure - atom counts', entry_id=d_handler.d_tag, value=True)
            else:
                pass

        # Combine structure_masks
        structure_reject_mask = self.get_dataset_masks().combine_masks(self._structure_mask_names)
        self.get_dataset_masks().add_mask(mask_name='rejected - structure', mask=structure_reject_mask)

        # Update the combined masks
        combined_reject_mask = self.get_dataset_masks().combine_masks(self._reject_mask_names)
        self.get_dataset_masks().add_mask(mask_name='rejected - total', mask=combined_reject_mask)

        self.log('\rDatasets Filtered.          ', True)
        self.log('Rejected Datasets (Structure): {!s}'.format(sum(self.get_dataset_masks().get_mask(mask_name='rejected - structure'))), True)
        self.log('Rejected Datasets (Total):     {!s}'.format(sum(self.get_dataset_masks().get_mask(mask_name='rejected - total'))), True)

    def filter_datasets_2(self):
        """Filter out the non-isomorphous datasets"""

        self.log('===================================>>>', True)
        self.log('Filtering Datasets (Non-isomorphous datasets). Potential Classes:', True)
        for failure_class in self._crystal_mask_names:
            self.log('\t{!s}'.format(failure_class), True)
        self.log('===================================>>>', True)

        # Check that each dataset is similar enough to be compared
        for d_handler in self.get_masked_datasets(mask_name='rejected - total', invert=True):

            d_handler.comparison_differences = self.reference_dataset().compare_structure().compare_all(d_handler.get_hierarchy())

            print '\rFiltering Dataset {!s}'.format(d_handler.d_tag),; sys.stdout.flush()
            # Check that it correlates well with itself before and after scaling
            if d_handler.get_scaled_fobs_miller_array().correlation(d_handler.get_fobs_miller_array()).coefficient() < 0.9:
                self.log('\rRejecting Dataset: {!s}'.format(d_handler.d_tag))
                self.log('Low correlation between scaled and unscaled data')
                self.log('Scaled-Unscaled Correlation: {!s}'.format(d_handler.get_scaled_fobs_miller_array().correlation(d_handler.get_fobs_miller_array()).coefficient()))
                self.log('===================================>>>')
                self.get_dataset_masks().set_mask_value(mask_name='bad crystal - data correlation', entry_id=d_handler.d_tag, value=True)
            # Check the resolution of the dataset
            elif d_handler.get_resolution() > self.get_cut_resolution():
                self.log('\rRejecting Dataset: {!s}'.format(d_handler.d_tag))
                self.log('Does not meet high-resolution cutoff')
                self.log('Cut Resolution: {!s}'.format(self.get_cut_resolution()))
                self.log('Map Resolution: {!s}'.format(d_handler.get_resolution()))
                self.log('===================================>>>')
                self.get_dataset_masks().set_mask_value(mask_name='bad crystal - resolution', entry_id=d_handler.d_tag, value=True)
            # Check the deviation from the average sites
            elif d_handler.get_calpha_sites().rms_difference(d_handler.transform_from_reference(self.reference_dataset().get_calpha_sites())) > 1:
                self.log('\rRejecting Dataset: {!s}'.format(d_handler.d_tag))
                self.log('C-alpha RMSD is too large')
                self.log('Initial (Calpha) RMSD: {!s}'.format(d_handler.get_calpha_sites().rms_difference(self.reference_dataset().get_calpha_average_sites())))
                self.log('Aligned (Calpha) RMSD: {!s}'.format(d_handler.get_calpha_sites().rms_difference(d_handler.transform_from_reference(self.reference_dataset().get_calpha_average_sites()))))
                self.log('===================================>>>')
                self.get_dataset_masks().set_mask_value(mask_name='bad crystal - isomorphous structure', entry_id=d_handler.d_tag, value=True)
            else:
                pass

        # Combine crystal masks
        crystal_reject_mask = self.get_dataset_masks().combine_masks(self._crystal_mask_names)
        self.get_dataset_masks().add_mask(mask_name='rejected - crystal', mask=crystal_reject_mask)

        # Combine all of the masks
        combined_reject_mask = self.get_dataset_masks().combine_masks(self._reject_mask_names)
        self.get_dataset_masks().add_mask(mask_name='rejected - total', mask=combined_reject_mask)

        self.log('\rDatasets Filtered.          ', True)
        self.log('Rejected Datasets (Structure): {!s}'.format(sum(self.get_dataset_masks().get_mask(mask_name='rejected - structure'))), True)
        self.log('Rejected Datasets (Crystal):   {!s}'.format(sum(self.get_dataset_masks().get_mask(mask_name='rejected - crystal'))), True)
        self.log('Rejected Datasets (Total):     {!s}'.format(sum(self.get_dataset_masks().get_mask(mask_name='rejected - total'))), True)

    def calculate_mean_structure_and_protein_masks(self, deviation_cutoff):
        """Calculate the average of all of the structures, and create masks for each protein where residues deviate from the mean by more than `deviation_cutoff`"""

        self.log('===================================>>>', True)
        self.log('Calculating Mean Structure', True)
        self.log('===================================>>>')

        # TODO Make this reject points until consensus

        # Pull all c-alpha sites for each structure
        all_sites = numpy.array([d.transform_to_reference(d.get_calpha_sites()) for d in self.get_masked_datasets(mask_name='rejected - total', invert=True)])
        # Calculate the mean x,y,z for each c-alpha
        mean_sites = numpy.mean(all_sites, axis=0)
        # Differences from the mean
        diff_sites = all_sites - mean_sites
        # Euclidean norms of the distances moved
        diff_norms = numpy.apply_along_axis(numpy.linalg.norm, axis=2, arr=diff_sites)

        # TODO MOVE THIS TO THE STRUCTURE VARIATION FUNCTION TODO

        # Create a list of masks for large-moving c-alphas
        residue_deviation_masks = []
        # Iterate by dataset, masking if the deviation of the calpha in the dataset is more than `deviation_cutoff`
        for calpha_shifts in diff_norms:
            residue_deviation_masks.append([1 if shift > deviation_cutoff else 0 for shift in calpha_shifts])

        # Save the masks
        self._average_calpha_sites = flex.vec3_double(mean_sites)
        self._residue_deviation_masks = residue_deviation_masks

        # Now calculate the variation in the structure, from the reference
        self.log('Calculating Variation in RMSD (Calphas) to Reference Structure')
        rmsds = numpy.array([-1]*self.get_number_of_datasets())
        [rmsds.put(d.d_num, d.get_calpha_sites().rms_difference(d.transform_from_reference(self.get_calpha_average_sites()))) for d in self.get_masked_datasets(mask_name='rejected - total', invert=True)]
        self.get_dataset_observations().add_data(data_name='rmsd_to_mean', data_values=rmsds.tolist())

    def collect_structure_variation_data(self):
        """Go through all of the datasets and collect lots of different structural characteristics of the datasets for identifying odd datasets"""

        self.log('===================================>>>', True)
        self.log('Calculating Structure Variation Data', True)

        backbone_atom_names = [" CA ", " C  ", " O  ", " N  "]

        # Extract the hierarchy for each of the datasets
        hierarchies = [d.get_hierarchy() for d in self.get_masked_datasets(mask_name='rejected - total', invert=True)]

        # Extract the reference residue types
        ref_residue_types = self.get_structure_observations().get_collection_info('residue types')

        self.log('===================================>>>', True)
        self.log('Calculating B-factors for residues across the datasets', True)
        self.log('===================================>>>')

        # Atomic selection cache for the reference structure - should be the same across the structures
        selection_cache = self.reference_dataset().get_hierarchy().atom_selection_cache()

        # String of 1-letter codes for the protein
        amino_sequence = ''

        for res_idx, res_id in enumerate(self.get_structure_observations().collection_ids):
            """Iterate through the residues and extract values from across the datasets"""

            res_collection = data_collection()
            res_collection.set_data_ids([d.d_tag for d in self.get_masked_datasets(mask_name='rejected - total', invert=True)])

            # Boolean mask for the selected residue
            res_bool_selection = selection_cache.selection("chain '{!s}' and resid '{!s}' and resname '{!s}'".format(*res_id))

            # Information to collect for the residues
            backbone_mean_bfactor = []
            sidechain_mean_bfactor = []
            total_mean_bfactor = []

            try:
                amino_sequence += pdb_reader.amino_acid_codes.one_letter_given_three_letter[res_id[2]]
            except:
                amino_sequence += '?'

            print '\rAnalysing Residue: {!s}-->'.format(amino_sequence),; sys.stdout.flush()

            # Iterate through the hierarchies and extract values for this residue
            for d_hierarchy in hierarchies:

                # Select the atoms for this residue
                new_root = d_hierarchy.select(res_bool_selection)

                # Calculate the sidechain and backbone b-factors for amino acids
                if ref_residue_types[res_idx] == 'common_amino_acid':

                    # Pull out the backbone atoms
                    backbone_atoms = [at for at in new_root.atoms() if at.name in backbone_atom_names]
                    if not backbone_atoms: backbone_mean_bfactor.append(None)
                    else:                  backbone_mean_bfactor.append(flex.mean(flex.double([at.b for at in backbone_atoms])))

                    # Pull out the sidechain atoms
                    sidechain_atoms = [at for at in new_root.atoms() if at.name not in backbone_atom_names]
                    if not sidechain_atoms: sidechain_mean_bfactor.append(None)
                    else:                   sidechain_mean_bfactor.append(flex.mean(flex.double([at.b for at in sidechain_atoms])))

                else:
                    # If not an amino acid, just append None
                    backbone_mean_bfactor.append(None)
                    sidechain_mean_bfactor.append(None)

                # Calculate the mean b-factor for the whole residue
                total_atoms = [at for at in new_root.atoms()]
                if not total_atoms: total_mean_bfactor.append(None)
                else:               total_mean_bfactor.append(flex.mean(flex.double([at.b for at in total_atoms])))

            res_collection.add_data(data_name='backbone_mean_bfactor', data_values=backbone_mean_bfactor)
            res_collection.add_data(data_name='sidechain_mean_bfactor', data_values=sidechain_mean_bfactor)
            res_collection.add_data(data_name='total_mean_bfactor', data_values=total_mean_bfactor)

            self.get_structure_observations().add_collection(collection_id=res_id, collection=res_collection)

        self.log('\rResidues Analysed: {!s}     '.format(amino_sequence), True)

    def load_all_map_handlers(self):
        """Load the objects for getting map values - these can't be pickled so they have to be loaded each time"""

        self.log('===================================>>>', True)
        for d_handler in self.get_masked_datasets(mask_name='rejected - total', invert=True):
            print '\rLoading Map Handlers {!s}'.format(d_handler.d_tag),; sys.stdout.flush()
            # Create maps
            d_handler.create_fft_map(miller=d_handler.get_scaled_fobs_miller_array(), map_type=self.get_map_type(), d_min=self.get_cut_resolution())
            # Scale maps
            if self.get_map_scaling() == 'none':
                pass
            elif self.get_map_scaling() == 'sigma':
                d_handler.get_map().apply_sigma_scaling()
            elif self.get_map_scaling() == 'volume':
                d_handler.get_map().apply_volume_scaling()
            d_handler.create_map_handler()
        self.log('\rMap Handlers Loaded.          ', True)

        self.update_pandda_size(tag='After Loading Map Handlers')

    def analyse_raw_maps(self):
        """Calculate the mean and std of the raw electron density maps over the masked grid points"""

        # Get masked points in reference dataset
        masked_gps = self.reference_grid().global_mask().outer_mask()
        masked_cart_ref = flex.vec3_double(masked_gps)*self.reference_grid().grid_spacing()

        # Calculate the indices of these grid points
        grid_indexer = self.reference_grid().grid_indexer()
        masked_idxs = [grid_indexer(gp) for gp in masked_gps]

        # List of statistics of map values
        map_stats = []

        self.log('===================================>>>')

        for d_handler in self.get_masked_datasets(mask_name='rejected - total', invert=True):
            # Retrieve associated files
            pdb = d_handler.get_pdb_filename()
            mtz = d_handler.get_mtz_filename()

            print '\rCalculating Maps Statistics for Dataset {!s} ({!s}/{!s}) '.format(d_handler.d_tag, d_handler.d_num+1, self.get_number_of_datasets()),; sys.stdout.flush()

            # Extract the map values at the transformed grid points
            masked_cart_new = d_handler.transform_from_reference(masked_cart_ref)
            masked_map_vals = d_handler.get_map_handler().get_cart_values(masked_cart_new)

            # Calculate the mean and std of the masked values
            ms = basic_statistics(flex.double(masked_map_vals))
            # Set the map offset and offset in the d_handler
            d_handler.set_raw_map_offset_and_scale(offset=-1.0*ms.mean, scale=1.0/ms.biased_standard_deviation)
            # Append to full list
            map_stats.append(ms)

        # Store the raw map stats for error spotting
        self.get_map_summary().set_data_length(len(map_stats))
        self.get_map_summary().add_data(data_name='raw_map_means', data_values=[m.mean for m in map_stats])
        self.get_map_summary().add_data(data_name='raw_map_stds', data_values=[m.biased_standard_deviation for m in map_stats])

        print ''
        self.log('===================================>>>', True)
        self.log('Maps Values Loaded: {!s} Datasets'.format(self.get_number_of_datasets(mask_name='rejected - total', invert=True)), True)

        self.update_pandda_size(tag='After Loading Maps')

    def calculate_map_statistics(self):
        """Take the sampled maps and calculate statistics for each grid point across the datasets"""

        # Create statistics objects for each grid point
        self.log('===================================>>>', True)
        self.log('Calculating Statistics of Grid Points', True)

        # Get masked points in reference dataset
        masked_gps = self.reference_grid().global_mask().outer_mask()
        masked_cart_ref = flex.vec3_double(masked_gps)*self.reference_grid().grid_spacing()

        # Calculate the indices of these grid points
        grid_indexer = self.reference_grid().grid_indexer()
        masked_idxs = [grid_indexer(gp) for gp in masked_gps]
        masked_num = len(masked_idxs)

        assert len(masked_idxs) == len(masked_cart_ref)

        # All dataset handlers for analysing
        all_d_handlers = self.get_masked_datasets(mask_name='rejected - total', invert=True)

        # Chunk the points into groups - Compromise between cpu time and memory usage - ~200 dataset -> chunksize of 5000
        chunk_size = 1000*int(2000/self.get_number_of_datasets(mask_name='rejected - total', invert=True))
        chunked_cart = [masked_cart_ref[i:i + chunk_size] for i in range(0, len(masked_cart_ref), chunk_size)]
        chunk_num = len(chunked_cart)

        self.log('Iterating through {!s} points in {!s} chunks'.format(masked_num, chunk_num), True)

        # Statistics objects for each of the points we're interested in
        masked_point_statistics = []

        # Calculate the statistics of the map values across the datasets
        for chunk_idx, ref_coord_chunk in enumerate(chunked_cart):
            status_bar(n=chunk_idx, n_max=chunk_num)
            p_map_vals = numpy.array([dh.get_cart_values(points=dh.transform_from_reference(flex.vec3_double(ref_coord_chunk))) for dh in all_d_handlers])

            if chunk_idx+1 < chunk_num:
                assert len(p_map_vals) == len(all_d_handlers)
                assert len(p_map_vals.T) == chunk_size
            [masked_point_statistics.append(basic_statistics(flex.double(map_vals.tolist()))) for map_vals in p_map_vals.T]

        # Calculate Mean Maps
        mean_map_vals = numpy.zeros(self.reference_grid().grid_size_1d())
        mean_map_vals.put(masked_idxs, [bs.mean for bs in masked_point_statistics])
        mean_map_vals = flex.double(mean_map_vals.tolist())
        mean_map_file = self.output_handler.get_file('reference_dataset').replace('.mtz','.mean.ccp4')
        self.write_array_to_map(output_file=mean_map_file, map_data=mean_map_vals)

        # Calculate Stds Maps
        stds_map_vals = numpy.ones(self.reference_grid().grid_size_1d())
        stds_map_vals.put(masked_idxs, [bs.biased_standard_deviation for bs in masked_point_statistics])
        stds_map_vals = flex.double(stds_map_vals.tolist())
        stds_map_file = self.output_handler.get_file('reference_dataset').replace('.mtz','.stds.ccp4')
        self.write_array_to_map(output_file=stds_map_file, map_data=stds_map_vals)

        # Calculate Skew Maps
        skew_map_vals = numpy.zeros(self.reference_grid().grid_size_1d())
        skew_map_vals.put(masked_idxs, [bs.skew for bs in masked_point_statistics])
        skew_map_vals = flex.double(skew_map_vals.tolist())
        skew_map_file = self.output_handler.get_file('reference_dataset').replace('.mtz','.skew.ccp4')
        self.write_array_to_map(output_file=skew_map_file, map_data=skew_map_vals)

        # Calculate Kurtosis Maps
        kurt_map_vals = numpy.zeros(self.reference_grid().grid_size_1d())
        kurt_map_vals.put(masked_idxs, [bs.kurtosis for bs in masked_point_statistics])
        kurt_map_vals = flex.double(kurt_map_vals.tolist())
        kurt_map_file = self.output_handler.get_file('reference_dataset').replace('.mtz','.kurt.ccp4')
        self.write_array_to_map(output_file=kurt_map_file, map_data=kurt_map_vals)

        # Calculate Bimodality Maps
        bimo_map_vals = numpy.zeros(self.reference_grid().grid_size_1d())
        bimo_map_vals.put(masked_idxs, [(bs.skew**2 + 1)/bs.kurtosis for bs in masked_point_statistics])
        bimo_map_vals = flex.double(bimo_map_vals.tolist())
        bimo_map_file = self.output_handler.get_file('reference_dataset').replace('.mtz','.bimo.ccp4')
        self.write_array_to_map(output_file=bimo_map_file, map_data=bimo_map_vals)

        # Store map vals
        self._mean_map = mean_map_vals
        self._stds_map = stds_map_vals
        self._skew_map = skew_map_vals
        self._kurt_map = kurt_map_vals
        self._bimo_map = bimo_map_vals

    def get_map(self, d_handler):
        """Extract the map values for the masked grid"""

        # Get masked points in reference dataset
        masked_gps = self.reference_grid().global_mask().outer_mask()
        masked_idxs = self.reference_grid().global_mask().outer_mask_indices()
        masked_cart_ref = flex.vec3_double(masked_gps)*self.reference_grid().grid_spacing()
        assert len(masked_gps) == len(masked_idxs), 'MASKED GPS DO NOT MATCH MASKED INDICES'

        # Transform Coordinates
        masked_cart_d = d_handler.transform_from_reference(masked_cart_ref)
        masked_vals_d = d_handler.get_cart_values(masked_cart_d)

        # Calculate the sparse vector of the masked map values
        map_sparse = sparse.vector(self.reference_grid().grid_size_1d(), dict(zip(masked_idxs, masked_vals_d)))

        return map_sparse.as_dense_vector()

    def calculate_z_map(self, map_vals, sparse=False):
        """Takes a flex.double map and calculates the z-map relative to the mean and std map"""

        assert map_vals.size() == self.reference_grid().grid_size_1d()

        # Calculate Z-values
        z_map_vals = (map_vals - self.get_mean_map())/self.get_stds_map()

        if sparse:
            # Calculate the sparse vector of the masked map values
            return sparse.vector(self.reference_grid().grid_size_1d(),
                                 dict(zip(self.reference_grid().global_mask().outer_mask_indices(),
                                          z_map_vals.select(self.reference_grid().global_mask().outer_mask_indices()))))
        else:
            return z_map_vals

    def process_z_map(self, z_map, local_mask_function=None):
        """Process the z_maps, looking for groups of high z-values"""

        # TODO DO IN PARALLEL!!!

        # Translate between 3d grid point and 1d array index
        grid_indexer = self.reference_grid().grid_indexer()

        # Set local_mask function if not given - set to rms
        if local_mask_function is None:
            # Given arguments: (map values, grid points (relative to central point))
            local_mask_function = lambda vals, gps: numpy.sqrt(numpy.mean(numpy.power(vals,2)))

        # Initialise empty sparse vector
        mod_z_map = sparse.vector(self.reference_grid().grid_size_1d())

        self.log('Post Processing Z-Map.', True)

        # Iterate through all grid points and calculate modified z-value at each point (0 at unsampled points or in buffer zone)
        for p_num, gp in enumerate(self.reference_grid().masked_grid_points()):

            status_bar(n=p_num, n_max=len(self.reference_grid().masked_grid_points()))

            # Find nearby values within local mask
            gp_masked = self.reference_grid().local_mask().apply_mask(gp)
            # Find the map values at these points
            local_map_vals = [z_map[grid_indexer(gp_m)] for gp_m in gp_masked]

            # Use the local mask function to modify the extracted map values (use local mask vectors so relative to gp)
            mod_z_val = local_mask_function(vals=local_map_vals, gps=self.reference_grid().local_mask().mask())

            # Add to new map array
            mod_z_map[grid_indexer(gp)] = mod_z_val

        # Get down-sampled map
        resamp_mod_z_map = [mod_z_map[grid_indexer(p)] for p in self.reference_grid().resampled_grid_points()]
        assert flex.product(flex.int(self.reference_grid().resampled_grid_size())) == len(resamp_mod_z_map)

        return mod_z_map.as_dense_vector(), flex.double(resamp_mod_z_map)

    def print_clustering_settings(self, z_cutoff,
                                        min_cluster_volume,
                                        clustering_cutoff,
                                        clustering_criterion = 'distance',
                                        clustering_metric    = 'euclidean',
                                        clustering_method    = 'average' ):

        # Calculate the approximate minimum volume for the cluster size
        min_cluster_size = int(min_cluster_volume/(self.reference_grid().grid_spacing()**3))

        self.log('===================================>>>', True)
        self.log('Z-Scores Clustering Scores', True)
        self.log('===================================>>>', True)
        self.log('Clustering Points with Z-Scores > {!s}'.format(z_cutoff), True)
        self.log('===================================>>>', True)
        self.log('Clustering Cutoff:           {!s}'.format(clustering_cutoff), True)
        self.log('Clustering Metric:           {!s}'.format(clustering_metric), True)
        self.log('Clustering Criterion:        {!s}'.format(clustering_criterion), True)
        self.log('Clustering Method (Linkage): {!s}'.format(clustering_method), True)
        self.log('===================================>>>', True)
        self.log('Minimum Cluster Size:        {!s}'.format(min_cluster_size), True)
        self.log('Minimum Cluster Volume:      {!s}'.format(min_cluster_volume), True)
        self.log('===================================>>>', True)

    def cluster_high_z_values(self, d_handler,
                                    z_map,
                                    z_cutoff,
                                    min_cluster_volume,
                                    clustering_cutoff,
                                    clustering_criterion = 'distance',
                                    clustering_metric    = 'euclidean',
                                    clustering_method    = 'average'):
        """Finds all the points in the z-map above `z_cutoff`, points will then be clustered into groups of cutoff `clustering_cutoff` angstroms"""

        # Check that the map is the expected size etc...
        assert len(z_map) == self.reference_grid().grid_size_1d()

        # Translate between grid coords and grid index
        grid_indexer = self.reference_grid().grid_indexer()

        # Scale the cutoff (Angstroms) into grid units
        grid_clustering_cutoff = clustering_cutoff/self.reference_grid().grid_spacing()

        # Calculate the approximate minimum volume for the cluster size
        min_cluster_size = int(min_cluster_volume/(self.reference_grid().grid_spacing()**3))

        # List of points in this dataset - TODO Revisit using global outer mask - is this the right mask?
        d_selected_points = [(gp, z_map[grid_indexer(gp)]) for gp in self.reference_grid().global_mask().outer_mask() if z_map[grid_indexer(gp)] >= z_cutoff]

        # No Cluster points found
        if not d_selected_points:
            d_handler.raw_cluster_hits = {}
            self.log('Dataset {!s}: No Clusters Found'.format(d_handler.d_tag), True)
        # Can't cluster if there are too many points
        elif len(d_selected_points) > 10000:
            self.log('===================================>>>', True)
            self.log('Dataset {!s}: Too many points to cluster: {!s} Points.'.format(d_handler.d_tag, len(d_selected_points)), True)
            d_handler.raw_cluster_hits = {}
            # This dataset is too noisy to analyse - flag!
            self.get_dataset_masks().set_mask_value(mask_name='noisy zmap', entry_id=d_handler.d_tag, value=True)

            # Still mark it as interesting!
            self.get_dataset_masks().set_mask_value(mask_name='interesting', entry_id=d_handler.d_tag, value=True)

            # Link datasets to the initial results directory
            hit_subdir = os.path.join(self.output_handler.get_dir('initial_hits'), 'Dataset-{!s}'.format(d_handler.d_tag))
            if not os.path.exists(hit_subdir):
                os.mkdir(hit_subdir)
                os.symlink(d_handler.get_pdb_filename(), os.path.join(hit_subdir, 'noisy_zmap.pdb'))
                os.symlink(d_handler.get_mtz_filename(), os.path.join(hit_subdir, 'noisy_zmap.mtz'))

            # Create maps of the high z-value points (significant points)
            highz_points = [tup[0] for tup in d_selected_points]
            highz_map_array = numpy.zeros(self.reference_grid().grid_size_1d(), dtype=int)
            highz_map_array.put(map(self.reference_grid().grid_indexer(), highz_points), [1]*len(highz_points))
            d_highz_map_file = os.path.join(hit_subdir, 'high_z_vals_mask.ccp4')
            self.write_array_to_map(d_highz_map_file, flex.double(highz_map_array.tolist()))

        # Cluster points if we have found them
        else:
            self.log('===================================>>>', True)
            self.log('Dataset {!s}: Clustering {!s} Point(s).'.format(d_handler.d_tag, len(d_selected_points)), True)
            # Points found for this cluster!
            if len(d_selected_points) == 1:
                # Only 1 point - 1 cluster! - edge case where we can't do clustering
                clust_num = 1
                # Dictionary of results
                clust_dict = {1: d_selected_points}
            else:
                # Extract only the coordinates and form an array
                point_array = numpy.array([tup[0] for tup in d_selected_points])
                # Cluster the extracted points
                t1 = time.time()
                clusts = list(scipy_cluster.hierarchy.fclusterdata(X=point_array, t=grid_clustering_cutoff, criterion=clustering_criterion, metric=clustering_metric, method=clustering_method))
                t2 = time.time()
                if t2-t1 > 30.0:
                    print('-> Clustering > Time Taken: {!s} seconds'.format(int(t2-t1)))

                # Get the number of clusters
                clust_num = max(clusts)
                # Initialise dictionary to hold the points for the different clusters (under the new indexing)
                clust_dict = dict([(i+1, []) for i in range(clust_num)])
                # Populate the clusters according to the clustering
                [clust_dict[c_idx].append(d_selected_points[p_idx]) for p_idx, c_idx in enumerate(clusts)]

            # Check that the clusters are numbered properly
            assert clust_num == max(clust_dict.keys())
            assert clust_num == len(clust_dict.keys())

            # Filter out small clusters - get numbers of clusters satisfying the minimum cluster size
            large_clusters = [c_num for c_num in range(1,clust_num+1) if len(clust_dict[c_num]) >= min_cluster_size]
            # Pull out the data for these filtered clusters and renumber them
            clust_dict = dict([(new_c_idx+1, clust_dict[old_c_num]) for new_c_idx, old_c_num in enumerate(large_clusters)])
            # Calculate the new number of clusters
            clust_num = len(large_clusters)

            if clust_num > 0:
                # Check that the clusters are numbered properly
                assert clust_num == max(clust_dict.keys())
                assert clust_num == len(clust_dict.keys())

                # Add clustered points to the dataset handler
                d_handler.raw_cluster_hits = clust_dict
                self.log('===> {!s} Cluster(s) found.'.format(clust_num), True)

                # This dataset is interesting!
                self.get_dataset_masks().set_mask_value(mask_name='interesting', entry_id=d_handler.d_tag, value=True)

                # Link datasets to the initial results directory
                hit_subdir = os.path.join(self.output_handler.get_dir('initial_hits'), 'Dataset-{!s}'.format(d_handler.d_tag))
                if not os.path.exists(hit_subdir):
                    os.mkdir(hit_subdir)
                    os.symlink(d_handler.get_pdb_filename(), os.path.join(hit_subdir, 'interesting.pdb'))
                    os.symlink(d_handler.get_mtz_filename(), os.path.join(hit_subdir, 'interesting.mtz'))

                # Create maps of the high z-value points (significant points)
                highz_points = []; [highz_points.extend([tup[0] for tup in clust_tups]) for clust_tups in clust_dict.values()]
                highz_map_array = numpy.zeros(self.reference_grid().grid_size_1d(), dtype=int)
                highz_map_array.put(map(self.reference_grid().grid_indexer(), highz_points), [1]*len(highz_points))
                d_highz_map_file = os.path.join(hit_subdir, 'high_z_vals_mask.ccp4')
                self.write_array_to_map(d_highz_map_file, flex.double(highz_map_array.tolist()))

            else:
                d_handler.raw_cluster_hits = {}
                self.log('===> No Clusters found - Minimum cluster size not reached.', True)

    def collate_all_clusters(self):
        """Collate clusters from all of the datasets"""

        self.log('===================================>>>', True)
        self.log('Collating Clusters', True)

        # List of points to be returned
        all_dataset_clusters = dict([(d.d_tag, []) for d in self.get_all_datasets()])

        for d_handler in self.get_all_datasets():

            if d_handler.raw_cluster_hits:
                all_dataset_clusters[d_handler.d_tag] = d_handler.raw_cluster_hits

        self.log('===================================>>>', True)
        # Print Cluster Summaries
        cluster_num = [(k, len(all_dataset_clusters[k])) for k in sorted(all_dataset_clusters.keys()) if all_dataset_clusters[k]]
        self.log('Total Datasets with Clusters: {!s}'.format(len(cluster_num)), True)
        cluster_total = sum([a[1] for a in cluster_num])
        self.log('Total Clusters: {!s}'.format(cluster_total), True)

        return cluster_total, cluster_num, all_dataset_clusters

    def process_z_value_clusters(self):
        """Analyse the clusters of points with high Z-values"""

        print('===================================>>>')
        print('Processing Clusters of Z-scores')

        # Summaries of the clusters
        cluster_summaries = []
        # Dataset numbers for those datasets with hits
        cluster_ids = []
        # All cluster objects
        all_clusters = []

        for d_handler in self.get_masked_datasets(mask_name='rejected - total', invert=True):

            # Check to see if there are any clustered points
            if not d_handler.raw_cluster_hits:
                continue

            print('===================================>>>')
            print('Processing Clusters in Dataset {!s}'.format(d_handler.d_tag))

            # Create cluster object from the clustered points
            cluster_obj = cluster_data(d_handler.raw_cluster_hits)
            d_handler.clustered_hits = cluster_obj
            # Add to list of all clusters
            cluster_ids.append(d_handler.d_tag)
            all_clusters.append(cluster_obj)

            # Calculate the spacings of the clusters
            dists = []
            total_c_num = len(cluster_obj.get_centroids())
            for cent_a in cluster_obj.get_centroids():
                for cent_b in cluster_obj.get_centroids():
                    if cent_a == cent_b: continue
                    dists.append((flex.double(cent_a) - flex.double(cent_b)).norm()*self.reference_grid().grid_spacing())
            if not dists: print 'Clusters, Dataset {!s}:'.format(d_handler.d_tag), '\tNum: {:4}'.format(total_c_num), '\tMin Spacing: {!s:>5} A'.format('--.--'), '\tMax Spacing: {!s:>5} A'.format('--.--')
            else:         print 'Clusters, Dataset {!s}:'.format(d_handler.d_tag), '\tNum: {:4}'.format(total_c_num), '\tMin Spacing: {:5.3} A'.format(min(dists)), '\tMax Spacing: {:5.3} A'.format(max(dists))

            # Pull out the cluster summaries
            for c_num, c_key in enumerate(cluster_obj.get_keys()):
                cluster_summaries.append([d_handler.d_tag, c_key] + cluster_obj.get_sizes([c_num]) + cluster_obj.get_means([c_num]) + cluster_obj.get_maxima([c_num]))

        # Fill out the cluster summary data
        self.get_cluster_summary().set_data_length(len(cluster_summaries))
        self.get_cluster_summary().add_data(data_name='dataset_id', data_values=[x[0] for x in cluster_summaries])
        self.get_cluster_summary().add_data(data_name='cluster_key', data_values=[x[1] for x in cluster_summaries])
        self.get_cluster_summary().add_data(data_name='cluster_size', data_values=[x[2] for x in cluster_summaries])
        self.get_cluster_summary().add_data(data_name='cluster_z_mean', data_values=[x[3] for x in cluster_summaries])
        self.get_cluster_summary().add_data(data_name='cluster_z_peak', data_values=[x[4] for x in cluster_summaries])

        # Combine all clusters into
        combined_cluster = combine_clusters(clusters=all_clusters, ids=cluster_ids)

        self.all_clustered_points = combined_cluster

        return combined_cluster

    def write_analysis_summary(self, output_file):
        """Writes an html summary of the datasets"""

        # Get template to be filled in
        template = PANDDA_HTML_ENV.get_template('output_table.html')
        # XXX 0=success, 1=none, 2=info, 3=warning, 4=failure

        dataset_observations = self.get_dataset_observations()

        # Construct the data object to populate the template
        output_data = {'PANDDA_TOP' : PANDDA_TOP}
        output_data['header'] = 'PANDDA Processing Output'
        output_data['title'] = 'PANDDA Processing Output'
        output_data['introduction'] = 'Summary of Processing of Datasets'
        output_data['table'] = {}
        output_data['table']['column_headings'] = ['Data Quality','Refinement','Isomorphous','Is Interesting']
        output_data['table']['rows'] = []
        # Add the datasets as rows
        for d in self.get_all_datasets():

            rmsd = round(dataset_observations.get_data('rmsd_to_mean')[d.d_num], 3)
            rwork, rfree = map(round, dataset_observations.get_data('rfree_rwork')[d.d_num],[3,3])

            columns = []
            overall_success = [0]
            # ------------------------------>>>
            # Test for Data Quality
            # ------------------------------>>>
            if self.get_dataset_masks().get_mask_value(mask_name='bad crystal - data quality', entry_id=d.d_tag) == True:
                columns.append({'flag':4,'message':'R-Merge: {!s}'})
            else:
                columns.append({'flag':0,'message':'OK'})
            # ------------------------------>>>
            # Test for Isomorphous - some test on crystal properties
            # ------------------------------>>>
            if self.get_dataset_masks().get_mask_value(mask_name='bad crystal - isomorphous crystal', entry_id=d.d_tag) == True:
                columns.append({'flag':4,'message':'Invalid Cell'})
            else:
                columns.append({'flag':0,'message':'OK'})
            # ------------------------------>>>
            # Test for Refinement Success - some test on r-free
            # ------------------------------>>>
            if self.get_dataset_masks().get_mask_value(mask_name='bad crystal - rfree', entry_id=d.d_tag) == True:
                columns.append({'flag':4,'message':'Free: {!s}, Work: {!s}'.format(rfree, rwork)})
            else:
                columns.append({'flag':0,'message':'Free: {!s}, Work: {!s}'.format(rfree, rwork)})
            # ------------------------------>>>
            # Test for Structure movement
            # ------------------------------>>>
            if self.get_dataset_masks().get_mask_value(mask_name='bad crystal - isomorphous structure', entry_id=d.d_tag) == True:
                columns.append({'flag':4,'message':'RMSD: {!s}'.format(rmsd)})
            else:
                columns.append({'flag':0,'message':'RMSD: {!s}'.format(rmsd)})
            # ------------------------------>>>
            # Test for if it's interesting
            # ------------------------------>>>
            if self.get_dataset_masks().get_mask_value(mask_name='interesting', entry_id=d.d_tag) == True:
                columns.append({'flag':0,'message':'Clusters Found'})
            else:
                columns.append({'flag':1,'message':''})

            # Find largest error
            overall_success = max([c['flag'] for c in columns])

            output_data['table']['rows'].append({'heading':'Dataset {!s}'.format(d.d_tag),
                                                 'success':overall_success,
                                                 'columns':columns})

        with open(os.path.join(self.outdir, self.output_handler.get_file(file_tag='summary_table')), 'w') as out_html:
            out_html.write(template.render(output_data))

    def write_array_to_map(self, output_file, map_data, grid_size=None, grid_spacing=None):
        """Takes a 1d array and writes it to a map"""
        if not grid_size:    grid_size    = self.reference_grid().grid_size()
        if not grid_spacing: grid_spacing = self.reference_grid().grid_spacing()
        self.log('> Writing Map: {!s}'.format(output_file))
        write_1d_array_as_p1_map(file_name=output_file, map_data=map_data, grid_size=grid_size, grid_spacing=grid_spacing)

    def write_grid_point_distributions(self, grid_points, output_filename):
        """Write CSV file of grid points, dataset numbers and map values"""

        self.log('===================================>>>', True)
        self.log('Writing Grid Points distributions: {!s} points'.format(len(grid_points)), True)
        self.log('Writing to {!s}'.format(output_filename), True)

        # Convert into cartesian coords
        cart_points_ref = flex.vec3_double(grid_points)*self.reference_grid().grid_spacing()

        # Iterate through wach dataset
        with open(os.path.join(self.output_handler.get_dir('manual_analyses'), output_filename), 'w') as fh:

            fh.write('grid_point, dataset_tag, dataset_num, map_val\n')
            for d_handler in self.get_masked_datasets(mask_name='rejected - total', invert=True):
                # Extract map values
                cart_points_d = d_handler.transform_from_reference(cart_points_ref)
                map_vals = d_handler.get_cart_values(cart_points_d)
                # Create list of tags for zipping
                d_tags = [d_handler.d_tag]*len(grid_points)
                d_nums = [d_handler.d_num]*len(grid_points)
                # Format and write
                out_list = [', '.join(map(str,tup)) for tup in zip(grid_points, d_tags, d_nums, map_vals)]
                out_line = '\n'.join(map(str,out_list)) + '\n'
                fh.write(out_line)

    def pickle(self, pickle_file, pickle_object, force=True):
        """Takes an object and pickles it"""
        if os.path.exists(pickle_file) and not force:
            self.log('NOT PICKLING: {!s}'.format(pickle_file))
        else:
            self.log('Pickling Object: {!s}'.format(pickle_file))
            easy_pickle.dump(pickle_file, pickle_object)

    def unpickle(self, pickle_file):
        """Takes an object and unpickles it"""
        self.log('Unpickling File: {!s}'.format(pickle_file))
        return easy_pickle.load(pickle_file)

#    def has_maps(self):
#        """Checks to see if the maps have been calculated - but does not load them if pickled"""
#        if self._maps!=[] or os.path.exists(self.pickle_handler.get_file('map_values')):
#            return True
#        return False
#    def get_maps(self, store=True):
#        """Returns maps for the used datasets - after scaling"""
#        if (self._maps==[]) and os.path.exists(self.pickle_handler.get_file('map_values')):
#            self.log('Unpickling Maps - this may take a while...', True)
#            maps = [sparse.vector(dimension=map_size, elements=map_dict) for map_size, map_dict in self.unpickle(self.pickle_handler.get_file('map_values'))]
#            if store: self._maps = maps
#            else: return maps
#        return self._maps
#    def has_z_maps(self):
#        """Checks to see if the maps have been calculated - but does not load them if pickled"""
#        if self._z_maps!=[] or os.path.exists(self.pickle_handler.get_file('z_map_values')):
#            return True
#        return False
#    def get_z_maps(self, store=True):
#        """Returns z-maps for the used datasets - loads from pickle if neccessary"""
#        if (self._z_maps==[]) and os.path.exists(self.pickle_handler.get_file('z_map_values')):
#            self.log('Unpickling Maps - this may take a while...', True)
#            maps = [sparse.vector(dimension=map_size, elements=map_dict) for map_size, map_dict in self.unpickle(self.pickle_handler.get_file('z_map_values'))]
#            if store: self._z_maps = maps
#            else: return maps
#        return self._z_maps
#    def has_modified_z_maps(self):
#        """Checks to see if the maps have been calculated - but does not load them if pickled"""
#        if self._mod_z_maps!=[] or os.path.exists(self.pickle_handler.get_file('mod_z_map_values')):
#            return True
#        return False
#    def get_modified_z_maps(self, store=True):
#        """Returns the post-processed, modified, z-maps for the used datasets"""
#        if (self._mod_z_maps==[]) and os.path.exists(self.pickle_handler.get_file('mod_z_map_values')):
#            self.log('Unpickling Maps - this may take a while...', True)
#            maps = [sparse.vector(dimension=map_size, elements=map_dict) for map_size, map_dict in self.unpickle(self.pickle_handler.get_file('mod_z_map_values'))]
#            if store: self._mod_z_maps = maps
#            else: return maps
#        return self._mod_z_maps
#
#    def _OLD_calculate_map_statistics(self):
#        """Take the sampled maps and calculate statistics for each grid point across the datasets"""
#
#        # Create statistics objects for each grid point
#        self.log('===================================>>>', True)
#        self.log('Calculating Statistics of Grid Points', True)
#        point_statistics = [basic_statistics(flex.double([map_data[i] for map_data in self.get_maps()])) for i in xrange(self.reference_grid().get_grid_size_1d())]
#        assert len(point_statistics) == self.reference_grid().get_grid_size_1d()
#        self.log('===================================>>>', True)
#
#        # Calculate Mean Maps
#        mean_map_vals = numpy.array([ps.mean for ps in point_statistics])
#        assert len(mean_map_vals) == self.reference_grid().get_grid_size_1d()
#        mean_map_file = self.output_handler.get_file('reference_dataset').replace('.mtz','.mean.ccp4')
#        self.write_array_to_map(output_file=mean_map_file, map_data=flex.double(mean_map_vals))
#
#        # Calculate Stds Maps
#        stds_map_vals = numpy.array([ps.biased_standard_deviation for ps in point_statistics])
#        assert len(stds_map_vals) == self.reference_grid().get_grid_size_1d()
#        stds_map_file = self.output_handler.get_file('reference_dataset').replace('.mtz','.stds.ccp4')
#        self.write_array_to_map(output_file=stds_map_file, map_data=flex.double(stds_map_vals))
#
#        # Calculate Skew Maps
#        skew_map_vals = numpy.array([ps.skew for ps in point_statistics])
#        assert len(skew_map_vals) == self.reference_grid().get_grid_size_1d()
#        skew_map_file = self.output_handler.get_file('reference_dataset').replace('.mtz','.skew.ccp4')
#        self.write_array_to_map(output_file=skew_map_file, map_data=flex.double(skew_map_vals))
#
#        # Calculate Kurtosis Maps
#        kurt_map_vals = numpy.array([ps.kurtosis for ps in point_statistics])
#        assert len(kurt_map_vals) == self.reference_grid().get_grid_size_1d()
#        kurt_map_file = self.output_handler.get_file('reference_dataset').replace('.mtz','.kurt.ccp4')
#        self.write_array_to_map(output_file=kurt_map_file, map_data=flex.double(kurt_map_vals))
#
#        # Calculate Bimodality Maps
#        bimo_map_vals = numpy.array([(skew_map_vals[i]**2 + 1)/kurt_map_vals[i] for i in xrange(self.reference_grid().get_grid_size_1d())])
#        assert len(bimo_map_vals) == self.reference_grid().get_grid_size_1d()
#        bimo_map_file = self.output_handler.get_file('reference_dataset').replace('.mtz','.bimo.ccp4')
#        self.write_array_to_map(output_file=bimo_map_file, map_data=flex.double(bimo_map_vals))
#
#        # Store map vals
#        self._mean_map = mean_map_vals.tolist()
#        self._stds_map = stds_map_vals.tolist()
#        self._skew_map = skew_map_vals.tolist()
#        self._kurt_map = kurt_map_vals.tolist()
#        self._bimo_map = bimo_map_vals.tolist()
#
#    def _old_extract_map_values(self):
#        """Extract map values for the filtered datasets (with scaling over the masked grid points)"""
#
#        assert self.get_maps() == [], 'Maps list is not empty!'
#
#        # Get masked points and binary grid mask (unsampled grid, outer mask only - max distance from protein)
#        masked_gps = self.reference_grid().global_mask().outer_mask()
#        binary_mask = self.reference_grid().get_global_outer_mask_binary()
#        # Translate between 3d grid point and 1d array index
#        grid_indexer = self.reference_grid().grid_indexer()
#
#        map_stats = []
#
#        for d_handler in self.get_masked_datasets(mask_name='rejected - total', invert=True):
#            # Retrieve associated files
#            pdb = d_handler.get_pdb_filename()
#            mtz = d_handler.get_mtz_filename()
#            self.log('===================================>>>')
#            self.log('Loading Maps for Dataset {!s}'.format(d_handler.d_tag))
#            self.log('===============>>>')
#            # Extract the map values at the transformed grid points
#            new_sample_points = d_handler.transform_from_reference(self.reference_grid().get_cart_points())
#            new_map_data = d_handler.get_cart_values(new_sample_points)
#
#            # Write the sampled map (raw_values)
#            new_sampled_map_file = mtz.replace('.mtz','.sampled.ccp4')
#            self.write_array_to_map(output_file=new_sampled_map_file, map_data=flex.double(new_map_data))
#
#            # Extract the map values for the masked grid to normalise maps
#            masked_map_values = [new_map_data[grid_indexer(gp)] for gp in masked_gps]
#
##            # Calculate the mean and standard deviation of the masked map
##            new_map_mean = numpy.mean(masked_map_values)
##            new_map_stdv = numpy.std(masked_map_values)
##            # Normalise the Map Values (roughly equivalent to sigma scaling, except over a subset of the map)
##            norm_map_array = normalise_array_to_z_scores(input_array=numpy.array(new_map_data), element_means=[new_map_mean]*len(new_map_data), element_stds=[new_map_stdv]*len(new_map_data), binary_mask=binary_mask)
#
#            # Sigma scale the map values relative to the masked points
#            norm_map_array = scipy_stats.zmap(scores=new_map_data, compare=masked_map_values)
#
#            # Convert to sparse vector (due to masking - lots of 0's)
#            norm_map_array_sparse = sparse.vector(len(new_map_data), dict([(i,v) for i,v in enumerate(norm_map_array) if binary_mask[i]==1]))
#            # Append the dataset, the map values to common lists
#            self._maps.append(norm_map_array_sparse)
#
#            # Calculate stats for the map
#            map_stats.append(basic_statistics(flex.double(masked_map_values)))
#
#        # Store the raw map stats for error spotting
#        self.get_map_summary().set_data_length(len(map_stats))
#        self.get_map_summary().add_data(data_name='map_means', data_values=[m.mean for m in map_stats])
#        self.get_map_summary().add_data(data_name='map_stds', data_values=[m.biased_standard_deviation for m in map_stats])
#        self.get_map_summary().add_data(data_name='map_stats', data_values=map_stats)
#
#        self.log('===================================>>>', True)
#        self.log('Maps Values Loaded: {!s} Datasets'.format(self.get_number_of_datasets(mask_name='rejected - total', invert=True)), True)
#
#        self.update_pandda_size(tag='After Loading Maps')
#
#    def old_normalise_maps_to_z_maps(self):
#        """Normalise the map values to z-values"""
#
#        assert self.get_z_maps() == [], 'Z-Maps list is not empty!'
#
#        # Get masked points (unsampled grid, outer mask only - max distance from protein)
#        masked_gps = self.reference_grid().global_mask().outer_mask()
#        binary_mask = self.reference_grid().get_global_outer_mask_binary()
#        # Translate between 3d grid point and 1d array index
#        grid_indexer = self.reference_grid().grid_indexer()
#
#        z_map_stats = []
#
#        for d_handler, map_data in zip(self.get_masked_datasets(mask_name='rejected - total', invert=True), self.get_maps()):
#            self.log('===================================>>>', True)
#            self.log('Normalising Dataset {!s}'.format(d_handler.d_tag), True)
#
#            pdb = d_handler.get_pdb_filename()
#            mtz = d_handler.get_mtz_filename()
#
#            # =====================================>>>
#            # Create Map-MeanMap Maps
#            mean_diff_array = normalise_array_to_z_scores(input_array=numpy.array(map_data.as_dense_vector()), element_means=self.get_mean_map(), element_stds=[1]*len(self.get_stds_map()), binary_mask=binary_mask)
#            # Write map
#            mean_diff_map_file = mtz.replace('.mtz','.mean_diff_values.ccp4')
#            self.write_array_to_map(output_file=mean_diff_map_file, map_data=flex.double(mean_diff_array))
#            # =====================================>>>
#
#            # Create Z-scores
#            z_array = normalise_array_to_z_scores(input_array=numpy.array(map_data.as_dense_vector()), element_means=self.get_mean_map(), element_stds=self.get_stds_map(), binary_mask=binary_mask)
#
#            # Write map
#            z_map_file = mtz.replace('.mtz','.zvalues.ccp4')
#            self.write_array_to_map(output_file=z_map_file, map_data=flex.double(z_array))
#
#            # Convert to sparse vector (due to masking - lots of 0's)
#            z_array_sparse = sparse.vector(len(z_array), dict([(i,v) for i,v in enumerate(z_array) if binary_mask[i]==1]))
#            # Store data
#            self._z_maps.append(z_array_sparse)
#
#            # Extract the z_values for the masked grid
#            masked_z_values = [z_array[grid_indexer(gp)] for gp in masked_gps]
#            # Calculate stats for the map
#            z_map_stats.append(basic_statistics(flex.double(masked_z_values)))
#
#        # Store the z map stats for error spotting
#        self.get_map_summary().add_data(data_name='z_map_mean', data_values=z_map_stats)
#
#        # No longer need maps - clear from memory after pickling
#        if (not self.keep_maps_in_memory) and (self._maps):
#            self.log('===================================>>>')
#            self.pickle(pickle_file=self.pickle_handler.get_file('map_values'), pickle_object=[(m.size, dict(m)) for m in self.get_maps()])
#            self._maps = []
#            gc.collect()
#
#    def old_post_process_z_maps(self, local_mask_function=None):
#        """Process the z_maps, looking for groups of high z-values"""
#
#        assert self.get_modified_z_maps() == [], 'Modified Z-Maps list is not empty!'
#
#        # Get the size of the down-sampled grid
#        resampled_grid_points = self.reference_grid().resampled_grid_points()
#        resampled_grid_size = self.reference_grid().resampled_grid_size()
#
#        # Get masked points (resampled grid, total mask - min+max distance from protein)
#        masked_gps = self.reference_grid().masked_grid_points()
#        # Translate between 3d grid point and 1d array index
#        grid_indexer = self.reference_grid().grid_indexer()
#
#        # Get the number of grid points in the local mask
#        local_mask_size = self.reference_grid().local_mask().size()
#
#        # Set local_mask function if not given - set to rms
#        if local_mask_function is None:
#            # Given arguments: (map values, grid points (relative to central point))
#            local_mask_function = lambda vals, gps: numpy.sqrt(numpy.mean(numpy.power(vals,2)))
#
#        # Calculate the redundancy required for the statistical tests
#        max_stat_sample_size = 344
#        redundancy = max(8, calculate_minimum_redundancy(unsampled_size=local_mask_size, max_sample_size=max_stat_sample_size))
#        self.log('===================================>>>')
#        self.log('Sampling Redundancy: {!s}'.format(redundancy))
#
#        for d_handler, z_map_data in zip(self.get_masked_datasets(mask_name='rejected - total', invert=True), self.get_z_maps()):
#            self.log('===================================>>>', True)
#            self.log('Post-Processing Dataset {!s}'.format(d_handler.d_tag), True)
#
#            pdb = d_handler.get_pdb_filename()
#            mtz = d_handler.get_mtz_filename()
#
#            mod_z_map_data = [0]*self.reference_grid().grid_size_1d()
#
#            # Iterate through all grid points and calculate modified z-value at each point (0 at unsampled points or in buffer zone)
#            for p_num, gp in enumerate(masked_gps):
#
#                status_bar(n=p_num, n_max=len(masked_gps))
#
#                # Find nearby values within local mask
#                gp_masked = self.reference_grid().local_mask().apply_mask(gp)
#                # Find the map values at these points
#                local_map_vals = [z_map_data[grid_indexer(gp_m)] for gp_m in gp_masked]
#
##                # Down-sample recorded values
##                resamp_map_vals = resample_ordered_list_of_values(vals=local_map_vals, redundancy=redundancy)
##                # Calculate significance of values
##                pval = test_significance_of_group_of_z_values(vals=resamp_map_vals)
##                # Convert to Z-score
##                mod_z_val = convert_pvalue_to_zscore(pval=pval)
#
#                # Use the local mask function to modify the extracted map values (use local mask vectors so relative to gp)
#                mod_z_val = local_mask_function(vals=local_map_vals, gps=self.reference_grid().local_mask().mask())
#
#                # Add to new map array
#                mod_z_map_data[grid_indexer(gp)] = mod_z_val
#
#            # Convert to sparse vector (due to masking - lots of 0's)
#            mod_z_map_data_sparse = sparse.vector(len(mod_z_map_data), dict([(i,v) for i,v in enumerate(mod_z_map_data) if v != 0]))
#
#            # Append to list of maps
#            self._mod_z_maps.append(mod_z_map_data_sparse)
#
#            # Write map
#            mod_z_map_file = mtz.replace('.mtz','.modzvalues.ccp4')
#            self.write_array_to_map(output_file=mod_z_map_file, map_data=flex.double(mod_z_map_data))
#
#            # Get down-sampled map
#            resamp_mod_z_map_data = [mod_z_map_data[grid_indexer(p)] for p in resampled_grid_points]
#            assert resampled_grid_size[0]*resampled_grid_size[1]*resampled_grid_size[2] == len(resamp_mod_z_map_data)
#
#            # Calculate stats for the map
#            # TODO Change this so that it's only over the masked points
##            resamp_mod_z_stats = basic_statistics(flex.double(samp_mod_z_map_data))
#
#            # Write down-sampled map
#            resamp_mod_z_map_file = mtz.replace('.mtz','.resamp.modzvalues.ccp4')
#            self.write_array_to_map(output_file=resamp_mod_z_map_file, map_data=flex.double(resamp_mod_z_map_data),
#                                    grid_size=resampled_grid_size, grid_spacing=self.reference_grid().resampled_grid_spacing())
#
#        # No longer need maps - clear from memory after pickling
#        if (not self.keep_maps_in_memory) and (self._z_maps):
#            self.log('===================================>>>')
#            self.pickle(pickle_file=self.pickle_handler.get_file('z_map_values'), pickle_object=[(m.size, dict(m)) for m in self.get_z_maps()])
#            self._z_maps = []
#            gc.collect()
#
#    def old_cluster_high_z_values(self, z_cutoff, min_cluster_size, clustering_cutoff,
#        clustering_criterion='distance', clustering_metric='euclidean', clustering_method='average'):
#        """Finds all the points in the modified z-maps above `z_cutoff`, points will then be clustered into groups of cutoff `clustering_cutoff` angstroms"""
#
#        # Translate between grid coords and grid index
#        grid_indexer = self.reference_grid().grid_indexer()
#
#        # Scale the cutoff (Angstroms) into grid units
#        grid_clustering_cutoff = clustering_cutoff/self.reference_grid().grid_spacing()
#
#        # Calculate the approximate minimum volume for the cluster size
#        min_cluster_volume = min_cluster_size * self.reference_grid().resampled_grid_spacing()**3
#
#        # List of points to be returned
#        all_dataset_clusters = dict([(d.d_num, []) for d in self.get_all_datasets()])
#
#        # Mask for the datasets that have high z-values
#        interesting_mask = [False]*self.get_number_of_datasets()
#
#        self.log('===================================>>>', True)
#        self.log('Clustering High Z-Scores', True)
#
#        self.log('===================================>>>', True)
#        self.log('Gathering Points with Z-Scores > {!s}'.format(z_cutoff), True)
#        self.log('===================================>>>', True)
#        self.log('Clustering Cutoff:           {!s}'.format(clustering_cutoff), True)
#        self.log('Clustering Metric:           {!s}'.format(clustering_metric), True)
#        self.log('Clustering Criterion:        {!s}'.format(clustering_criterion), True)
#        self.log('Clustering Method (Linkage): {!s}'.format(clustering_method), True)
#        self.log('===================================>>>', True)
#        self.log('Minimum Cluster Size:        {!s}'.format(min_cluster_size), True)
#        self.log('Minimum Cluster Volume:      {!s}'.format(min_cluster_volume), True)
#        self.log('===================================>>>', True)
#
#        # Iterate through the mod_z_maps, extract d_num, coord_num (grid_index), and mod_z_val
#        for d_handler, mz_map_data in zip(self.get_masked_datasets(mask_name='rejected - total', invert=True), self.get_modified_z_maps()):
#            # List of points in this dataset
#            d_selected_points = [(gp, mz_map_data[grid_indexer(gp)]) for gp in self.reference_grid().masked_grid_points() if mz_map_data[grid_indexer(gp)] >= z_cutoff]
#
#            # Cluster points if we have found them
#            if d_selected_points:
#                self.log('===================================>>>', True)
#                self.log('Dataset {!s}: Clustering {!s} Point(s).'.format(d_handler.d_tag, len(d_selected_points)), True)
#                # Points found for this cluster!
#                if len(d_selected_points) == 1:
#                    # Only 1 point - 1 cluster! - edge case where we can't do clustering
#                    clust_num = 1
#                    # Dictionary of results
#                    clust_dict = {1: d_selected_points}
#                else:
#                    # Extract only the coordinates and form an array
#                    point_array = numpy.array([tup[0] for tup in d_selected_points])
#                    # Cluster the extracted points
#                    clusts = list(scipy_cluster.hierarchy.fclusterdata(X=point_array, t=grid_clustering_cutoff, criterion=clustering_criterion, metric=clustering_metric, method=clustering_method))
#
#                    # Get the number of clusters
#                    clust_num = max(clusts)
#                    # Initialise dictionary to hold the points for the different clusters (under the new indexing)
#                    clust_dict = dict([(i+1, []) for i in range(clust_num)])
#                    # Populate the clusters according to the clustering
#                    [clust_dict[c_idx].append(d_selected_points[p_idx]) for p_idx, c_idx in enumerate(clusts)]
#
#                # Check that the clusters are numbered properly
#                assert clust_num == max(clust_dict.keys())
#                assert clust_num == len(clust_dict.keys())
#
#                # Filter out small clusters - get numbers of clusters satisfying the minimum cluster size
#                large_clusters = [c_num for c_num in range(1,clust_num+1) if len(clust_dict[c_num]) >= min_cluster_size]
#                # Pull out the data for these filtered clusters and renumber them
#                clust_dict = dict([(new_c_idx+1, clust_dict[old_c_num]) for new_c_idx, old_c_num in enumerate(large_clusters)])
#                # Calculate the new number of clusters
#                clust_num = len(large_clusters)
#
#                if clust_num > 0:
#                    # Check that the clusters are numbered properly
#                    assert clust_num == max(clust_dict.keys())
#                    assert clust_num == len(clust_dict.keys())
#
#                    # Add clustered points to the dataset handler
#                    d_handler.raw_cluster_hits = clust_dict
#                    # Append to the global list
#                    all_dataset_clusters[d_handler.d_num] = clust_dict
#                    self.log('===> {!s} Cluster(s) found.'.format(clust_num), True)
#                    # This dataset is interesting!
#                    interesting_mask[d_handler.d_num] = True
#                else:
#                    self.log('===> No Clusters found - Minimum cluster size not reached.', True)
#                self.log('===================================>>>', True)
#            else:
#                self.log('Dataset {!s}: No Clusters Found'.format(d_handler.d_tag), True)
#
#        # Expand the mask to include rejected datasets
#        self.get_dataset_masks().add_mask(mask_name='interesting', mask=interesting_mask)
#
#        # Print Cluster Summaries
#        cluster_num = [(k, len(all_dataset_clusters[k])) for k in sorted(all_dataset_clusters.keys()) if all_dataset_clusters[k]]
#        print('Total Datasets with Clusters: {!s}'.format(len(cluster_num)))
#        cluster_total = sum([a[1] for a in cluster_num])
#        print('Total Clusters: {!s}'.format(cluster_total))
#
#        return all_dataset_clusters

# TODO Move to Giant.Grid.Handler(?)
class grid_handler(object):
    def __init__(self, verbose=True):
        """Create and manage a grid object to be sampled across many aligned datasets"""

        self.verbose = verbose

        # Size of complete grid
        self._grid_size = None
        self._grid_spacing = None

        # Cartesian values of grid
        self._cart_max = None
        self._cart_min = None
        self._cart_size = None
        self._cart_points = None

        # Local mask for filtering/smoothing
        self._local_mask = None
        # Global mask for selecting regions
        self._global_mask = None

        # Groups of points defined by the local mask
        self._buffer_mask_points = None
        self._buffer_mask_binary = None

        # Group of points defined by sampling the grid regularly
        self._resampled_grid_points = None
        self._resampled_grid_size = None

        # Manually set group of points created by combining masks
        self._masked_grid_points = None

    def set_grid_spacing(self, spacing):
        self._grid_spacing = spacing
    def grid_spacing(self):
        return self._grid_spacing
    def grid_size(self):
        return self._grid_size
    def grid_size_1d(self):
        return self._grid_size[0]*self._grid_size[1]*self._grid_size[2]
    def grid_point_volume(self):
        return self.grid_spacing()**3

    def set_cart_extent(self, cart_min, cart_max):
        self._cart_min = cart_min
        self._cart_max = cart_max
        self._cart_size = tuple([s1 - s2 for s1,s2 in zip(cart_max, cart_min)])
    def cart_max(self):
        return self._cart_max
    def cart_min(self):
        return self._cart_min
    def cart_extent(self):
        return (self._cart_min, self._cart_max)
    def cart_size(self):
        return self._cart_size

    def cart_points(self):
        return self._cart_points
    def grid_points(self):
        return flex.nested_loop(self.grid_size())

    def grid_indexer(self):
        """Translates between 3d grid coordinates and 1d array coordinates"""
        return flex.grid(self.grid_size())

    def fake_unit_cell(self):
        """Create a unit cell as if the reference grid were a real lattice"""
        return crystal.uctbx.unit_cell('{!s} {!s} {!s} 90 90 90'.format(*self.cart_size()))

    def fake_space_group(self):
        """Create a spacegroup as if the reference grid were a real lattice"""
        return crystal.sgtbx.space_group('P1')

    def set_local_mask(self, mask):
        """Add local mask to the grid object"""

        self._local_mask = mask
        self.local_mask().show_summary()

        grid_size = self.grid_size()
        grid_jump = self.local_mask().grid_jump()
        grid_buffer = self.local_mask().buffer_size()
        grid_indexer = self.grid_indexer()

        # Resample grid points
        self._resampled_grid_points = [gp for gp in flex.nested_loop(grid_size) if not [1 for coord in gp if (coord%grid_jump != 0)]]
        self._resampled_grid_size = tuple([int(1+(g-1)/grid_jump) for g in grid_size])
        self._resampled_grid_spacing = grid_jump*self.grid_spacing()

        # Create a buffer zone at the edge of the grid
        self._buffer_mask_points = [gp for gp in flex.nested_loop(grid_size) if [1 for i_dim, coord in enumerate(gp) if (coord<grid_buffer or coord>(grid_size[i_dim]-1-grid_buffer))]]

        # Create binary mask for the buffer zone
        buffer_mask_binary = numpy.zeros(self.grid_size_1d(), int)
        [buffer_mask_binary.put(grid_indexer(gp), 1) for gp in self._buffer_mask_points]
        self._buffer_mask_binary = buffer_mask_binary.tolist()

    def local_mask(self):
        return self._local_mask

    def resampled_grid_points(self):
        """Get a down-sampled list of grid points, based on the local mask used"""
        return self._resampled_grid_points
    def resampled_grid_size(self):
        """Gets the size of the re-sampled grid"""
        return self._resampled_grid_size
    def resampled_grid_spacing(self):
        """Gets the grid spacing for the re-sampled grid"""
        return self._resampled_grid_spacing
    def resampled_grid_indexer(self):
        """Translates between 3d grid coordinates and 1d array coordinates"""
        return flex.grid(self.resampled_grid_size())

    def buffer_mask_points(self):
        """Get a list of points in the buffer zone of the map"""
        return self._buffer_mask_points
    def buffer_mask_binary(self):
        """Get a binary mask for the buffer zone"""
        return self._buffer_mask_binary

    def set_global_mask(self, mask):
        """Add a global mask to the grid object - This will create binary mask of the masked grid points"""

        self._global_mask = mask
        self.global_mask().show_summary()

    def global_mask(self):
        return self._global_mask

    def set_masked_grid_points(self, masked_points):
        self._masked_grid_points = masked_points
    def masked_grid_points(self):
        return self._masked_grid_points

    def show_summary(self):
        print '===================================>>>'
        print('Reference Grid Summary:')
        print('Grid Spacing:        {!s}'.format(round(self.grid_spacing(), 3)))
        print('Grid Point Volume:   {!s}'.format(round(self.grid_point_volume(),3)))
        print('Size of Grid (3D):   {!s}'.format(self.grid_size()))
        print('Size of Grid (1D):   {!s}'.format(self.grid_size_1d()))
        print('Min of Grid (Cart): {!s}'.format(tuple([round(x,3) for x in self.cart_min()])))
        print('Max of Grid (Cart): {!s}'.format(tuple([round(x,3) for x in self.cart_max()])))
        print('Size of Grid (Cart): {!s}'.format(tuple([round(x,3) for x in self.cart_size()])))

    def create_cartesian_grid(self, expand_to_origin):
        if expand_to_origin:
            # TODO Don't know if this will still work - raise error for now
            raise Exception('NOT CURRENTLY CHECKED')
            assert [i>0 for i in self.cart_min()], 'ALL GRID SITES MUST BE GREATER THAN 0 IF ORIGIN INCLUDED'
            assert [i>0 for i in self.cart_max()], 'ALL GRID SITES MUST BE GREATER THAN 0 IF ORIGIN INCLUDED'
            box_size, self._grid_size, self._cart_points = create_cartesian_grid(min_carts=(0,0,0),
                                                                                 max_carts=self.cart_max(),
                                                                                 grid_spacing=self.grid_spacing())
        else:
            box_size, self._grid_size, self._cart_points = create_cartesian_grid(min_carts=self.cart_min(),
                                                                                 max_carts=self.cart_max(),
                                                                                 grid_spacing=self.grid_spacing())

        print 'BOX SIZE: {!s}'.format(box_size)

        # Update max/min cart sizes as the grid will be slightly larger than the requested size
        self.set_cart_extent(cart_min=self.cart_points().min(), cart_max=self.cart_points().max())

        return self.grid_size()

    def _old_create_local_mask(self, distance_cutoff, mask='sphere'):
        """Generate a sphere of grid points within distance_cutoff of the origin"""

        print('===================================>>>')
        print('Generating Local Mask')
        if mask == 'sphere':
            self._local_mask = masking_sphere(self.get_grid_spacing(), distance_cutoff=distance_cutoff)
        else:
            raise Exception('Mask not found')

        self.local_mask().show_summary()

    def _old_create_global_mask(self, cart_sites, distance_cutoff):
        """Get the mask of all grid points within distance_cutoff of any of cart_sites"""

        print('===================================>>>')
        print('Generating Global Mask')
        self._global_mask = global_mask(cart_sites=cart_sites, grid_spacing=self.get_grid_spacing(), distance_cutoff=distance_cutoff)
        self.global_mask().show_summary()

# TODO Move to Giant.Grid.Masks
class atomic_mask(object):
    def __init__(self, cart_sites, grid_size, unit_cell, max_dist, min_dist):
        """Take a grid and calculate all grid points with a certain distance cutoff of any point in cart_sites"""

        if min_dist: assert max_dist > min_dist, 'Minimum Mask Distance must be smaller than Maximum Mask Distance'
        print('===================================>>>')

        # TODO Check that all of the points lie withing the fake unit cell
        # i.e. cart_sites + max_dist < unit_cell.parameters()

        # Store distances from masking atoms
        self._max_dist = max_dist
        self._min_dist = min_dist

        # Store grid size
        self._grid_size = grid_size
        self._grid_idxr = flex.grid(grid_size)

        # Unit cell
        self._fake_unit_cell = unit_cell

        t1 = time.time()

        # Calculate the masked indices defined by max distance from protein atoms
        self._outer_mask_indices = maptbx.grid_indices_around_sites(unit_cell=unit_cell,
                                    fft_n_real=grid_size, fft_m_real=grid_size,
                                    sites_cart=cart_sites, site_radii=flex.double(cart_sites.size(), max_dist))
        # Calculate a binary list of the selected indices
        outer_mask_binary = numpy.zeros(self._grid_idxr.size_1d(), dtype=int)
        [outer_mask_binary.put(idx, 1) for idx in self.outer_mask_indices()]
        self._outer_mask_binary = outer_mask_binary.tolist()
        # Select the masked grid points
        self._outer_mask_points = [gp for idx, gp in enumerate(flex.nested_loop(self._grid_size)) if self._outer_mask_binary[idx]==1]

        t2 = time.time()
        print('-> OUTER MASK > Time Taken: {!s} seconds'.format(int(t2-t1)))

        # Calculate the masked indices defined by min distance from protein atoms
        self._inner_mask_indices = maptbx.grid_indices_around_sites(unit_cell=unit_cell,
                                    fft_n_real=grid_size, fft_m_real=grid_size,
                                    sites_cart=cart_sites, site_radii=flex.double(cart_sites.size(), min_dist))
        # Calculate a binary list of the selected indices
        inner_mask_binary = numpy.zeros(self._grid_idxr.size_1d(), dtype=int)
        [inner_mask_binary.put(idx, 1) for idx in self.inner_mask_indices()]
        self._inner_mask_binary = inner_mask_binary.tolist()
        # Select the masked grid points
        self._inner_mask_points = [gp for idx, gp in enumerate(flex.nested_loop(self._grid_size)) if self._inner_mask_binary[idx]==1]

        t3 = time.time()
        print('-> INNER MASK > Time Taken: {!s} seconds'.format(int(t3-t2)))

        # Calculate the combination of these masks
        total_mask_binary = numpy.zeros(self._grid_idxr.size_1d(), int)
        [total_mask_binary.put(self._grid_idxr(gp), 1) for gp in self.outer_mask()]
        [total_mask_binary.put(self._grid_idxr(gp), 0) for gp in self.inner_mask()]
        self._total_mask_binary = total_mask_binary.tolist()
        # Select the masked grid indices
        self._total_mask_indices = [gp for idx, gp in enumerate(flex.nested_loop(self._grid_size)) if self._total_mask_binary[idx]==1]
        # Select the masked grid points
        self._total_mask_points = [gp for idx, gp in enumerate(flex.nested_loop(self._grid_size)) if self._total_mask_binary[idx]==1]

        t4 = time.time()
        print('-> TOTAL MASK > Time Taken: {!s} seconds'.format(int(t4-t3)))

    def total_mask(self):
        """Return the grid points allowed by the mask - combination of max_dist (allowed) and min_dist (rejected)"""
        return self._total_mask_points
    def outer_mask(self):
        """Get grid points allowed subject to max_dist"""
        return self._outer_mask_points
    def inner_mask(self):
        """Get grid points rejected subject to min_dist"""
        return self._inner_mask_points

    def total_mask_binary(self):
        return self._total_mask_binary
    def outer_mask_binary(self):
        return self._outer_mask_binary
    def inner_mask_binary(self):
        return self._inner_mask_binary

    def total_mask_indices(self):
        return self._total_mask_indices
    def outer_mask_indices(self):
        return self._outer_mask_indices
    def inner_mask_indices(self):
        return self._inner_mask_indices

    def total_size(self):
        """Returns the number of grid points in the mask"""
        return len(self.total_mask())
    def outer_size(self):
        """Returns the number of grid points inside max_dist"""
        return len(self.outer_mask())
    def inner_size(self):
        """Returns the number of grid points inside max_dist"""
        return len(self.inner_mask())

    def extent(self):
        """Returns the minimum and maximum grid points in the mask"""
        return min(self.total_mask()), max(self.total_mask())

    def show_summary(self):
        print('===================================>>>')
        print('Global Mask Summary:')
        print('Total Mask Size (1D): {!s}'.format(self.total_size()))
        print('Outer Mask Size (1D): {!s}'.format(self.outer_size()))
        print('Inner Mask Size (1D): {!s}'.format(self.inner_size()))
        print('Masked Grid Min/Max: {!s}'.format(min(self.extent())))

# TODO Move to Giant.Grid.Masks
class spherical_mask(object):
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

    def mask(self):
        return self._mask
    def size(self):
        return len(self.mask())
    def buffer_size(self):
        return self._buffer
    def grid_spacing(self):
        return self._grid_spacing
    def grid_jump(self):
        return self._grid_jump
    def radius(self):
        return self._radius
    def volume(self):
        return (4/3.0)*numpy.pi*(self.radius()**3)

    def apply_mask(self, grid_point):
        """Combine a grid point with all of the masking vectors"""
        return combine_grid_point_and_grid_vectors(start_point=grid_point, grid_vectors=self.mask())

    def show_summary(self):
        print('===================================>>>')
        print('Local Mask Summary:')
        print('Number of Mask Points:  {!s}'.format(len(self.mask())))
        print('Mask Radius (Cart):     {!s}'.format(self.radius()))
        print('Mask Volume (Cart):     {!s}'.format(round(self.volume(),3)))
        print('Largest Mask Vector:    {!s}'.format(max(self.mask())))
        print('Req. Edge Buffer Zone:  {!s}'.format(self.buffer_size()))
        print('Sampling Fraction (1D): 1/{!s}'.format(self.grid_jump()))
        print('Sampling Fraction (3D): 1/{!s}'.format(self.grid_jump()**3))



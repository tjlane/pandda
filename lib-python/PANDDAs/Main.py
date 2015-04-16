import os, sys, glob, time, re
import copy, resource, gc
import multiprocessing

import math

from scipy import spatial
from scipy import cluster as scipy_cluster

import numpy

import iotbx.pdb
import iotbx.mtz
import mmtbx.f_model as model_handler
#import iotbx.map_tools as map_tools

from libtbx import easy_mp, easy_pickle
from libtbx import phil

import cctbx.miller

from cctbx import maptbx
from cctbx import crystal
from scitbx import sparse

from scitbx.array_family import flex
from scitbx.math import superpose, basic_statistics
from scitbx.math.distributions import normal_distribution

from libtbx.math_utils import ifloor, iceil

#from iotbx.reflection_file_utils import extract_miller_array_from_file

from Bamboo.Common import compare_dictionaries
from Bamboo.Common.File import output_file_object, easy_directory
from Bamboo.Common.Data import data_collection, multiple_data_collection
from Bamboo.Common.Masks import mask_collection

#from Bamboo.Density.Edstats.Score import score_file_with_edstats
from Bamboo.Density.Edstats.Utils import score_with_edstats_to_dict

from Giant.Grid import grid_handler

from Giant.Xray.Miller.Scaling import apply_simple_scaling
from Giant.Xray.Maps.Utils import write_1d_array_as_p1_map

from Giant.Stats.Tests import test_significance_of_group_of_z_values, convert_pvalue_to_zscore
from Giant.Stats.Utils import resample_ordered_list_of_values, calculate_minimum_redundancy
from Giant.Stats.Ospina import estimate_true_underlying_sd

from Giant.Stats.Cluster import cluster_data, combine_clusters

from Giant.Structure.Align import perform_flexible_alignment

from Giant.Utils import status_bar

from PANDDAs.HTML import PANDDA_HTML_ENV
from PANDDAs.Phil import pandda_phil_def
from PANDDAs.Settings import PANDDA_TOP, PANDDA_TEXT

diff_map_hash = {'2mFo-DFc' : 'mFo-DFc',
                 'mFo'      : 'mFo-DFc',
                 'Fo'       : 'Fo-DFc'}

col_hash = { '2mFo-DFc' : {'AMP':'FWT'    , 'AMPWT':None  , 'PHA':'PHWT'    } ,
             'mFo-DFc'  : {'AMP':'DELFWT' , 'AMPWT':None  , 'PHA':'PHDELWT' } ,
             'mFo'      : {'AMP':'F'      , 'AMPWT':'FOM' , 'PHA':'FC'      }
            }

STRUCTURE_MASK_NAMES = [    'bad structure - chain counts',
                            'bad structure - chain ids',
                            'bad structure - chain sequences',
                            'bad structure - residue counts',
                            'bad structure - atom counts',
                            'bad structure - non-identical structures'  ]
CRYSTAL_MASK_NAMES   = [    'bad crystal - isomorphous crystal',
                            'bad crystal - isomorphous structure',
                            'bad crystal - space group',
                            'bad crystal - data quality',
                            'bad crystal - data correlation',
                            'bad crystal - rfree'  ]
REJECT_MASK_NAMES    = [    'rejected - total',
                            'rejected - crystal',
                            'rejected - structure',
                            'rejected - unknown'  ]
FLAG_MASK_NAMES      = [    'noisy zmap',
                            'interesting'  ]

# SET FILTERING CONSTANTS
FILTER_SCALING_CORRELATION_CUTOFF = 0.7
FILTER_GLOBAL_RMSD_CUTOFF = 1

PANDDA_VERSION = 0.1

def rel_symlink(orig, link):
    """Make a relative symlink from link to orig"""
    assert os.path.exists(orig), 'FILE DOES NOT EXIST: {!s}'.format(orig)
    assert not os.path.exists(link), 'LINK ALREADY EXISTS: {!s}'.format(link)
    orig = os.path.abspath(orig)
    link = os.path.abspath(link)
    assert not link.endswith('/'), 'LINK CANNOT END WITH /'
    os.symlink(os.path.relpath(orig, start=os.path.dirname(link)), link)

def map_handler(map_data, unit_cell):
    """Map handler for easy sampling of map"""

    basic_map = maptbx.basic_map(   maptbx.basic_map_unit_cell_flag(),
                                    map_data,
                                    map_data.focus(),
                                    unit_cell.orthogonalization_matrix(),
                                    maptbx.out_of_bounds_clamp(0).as_handle(),
                                    unit_cell
                                )

    return basic_map

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

class dataset_handler_list(object):
    """Class for grouping many dataset handlers together"""

    def __init__(self):
        self._datasets = []
        self._masks = mask_collection()

    def __getitem__(self, idx):
        if   isinstance(idx, int): return self.get(d_num = idx)
        elif isinstance(idx, str): return self.get(d_tag = idx)
        else: raise Exception('CANNOT INDEX EXCEPT BY int OR str. TYPE GIVEN: {!s}'.format(type(idx)))

    def __call__(self):
        """Return all datasets"""
        return self._datasets

    def all_masks(self):
        """Return the mask object"""
        return self._masks

    def all(self):
        """Return all datasets"""
        return self._datasets
    def all_nums(self):
        """Return the list of dataset ids"""
        return [d.d_num for d in self.all()]
    def all_tags(self):
        """Return the list of dataset tags"""
        return [d.d_tag for d in self.all()]

    def size(self, mask_name=None, invert=False):
        """Return the number of datasets in the handler (with optional mask applied)"""
        if mask_name:   return len(self.mask(mask_name=mask_name, invert=invert))
        else:           return len(self.all())

    def add(self, new_datasets):
        """Add new datasets"""

        for new_d in new_datasets:
            # Check all added datasets are the right class
            assert isinstance(new_d, dataset_handler), 'ADDED OBJECTS MUST BE OF TYPE: {!s}'.format(dataset_handler.__name__)
            # Check all added dataset id tags are strs
            assert isinstance(new_d.d_tag, str), 'DATASET TAG MUST BE str. Type given: {!s}'.format(type(new_d.d_tag))
            assert new_d.d_tag not in self.all_tags(), 'DATASET TAG ALREADY EXISTS: {!s}'.format(new_d.d_tag)
            # Check all added dataset id nums are ints
            assert isinstance(new_d.d_num, int), 'DATASET NUM MUST BE int: {!s}'.format(type(new_d.d_num))
            assert new_d.d_num not in self.all_nums(), 'DATASET NUM ALREADY EXISTS: {!s}'.format(new_d.d_num)
            # No problems, add to list
            self._datasets.append(new_d)

    def get(self, d_tag=None, d_num=None):
        """Get a dataset by tag or num"""

        assert d_num or d_tag, 'No dataset NUM or TAG given!'
        assert not (d_num and d_tag), 'BOTH dataset NUM and TAG given!'
        if d_num: matching = [d for d in self.all() if d.d_num == d_num]
        else:     matching = [d for d in self.all() if d.d_tag == d_tag]
        if len(matching) == 0: raise Exception('NO MATCHING DATASET FOUND - DNUM: {!s}, DTAG: {!s}'.format(d_num, d_tag))
        if len(matching) != 1: raise Exception('MORE THAN ONE MATCHING DATASET - DNUM: {!s}, DTAG: {!s}'.format(d_num, d_tag))
        return matching[0]

    def mask(self, mask_name, invert=False):
        """Retrieve a masked list of datasets"""
        return self._masks.mask(mask_name=mask_name, input_list=self.all(), invert=invert)

class dataset_handler(object):
    def __init__(self, dataset_number, pdb_filename, mtz_filename, dataset_tag=None):
        """Create a dataset object to allow common functions to be applied easily to a number of datasets"""

        assert os.path.exists(pdb_filename), 'PDB file does not exist!'
        assert os.path.exists(mtz_filename), 'MTZ file does not exist!'

        # Store dataset number
        self.d_num = dataset_number
        # Store the tag for the dataset
        if dataset_tag:
            self.d_tag = dataset_tag
        else:
            # If d_num < 0 - mark as a reference dataset
            if self.d_num < 0:  self.d_tag = 'REF{:05d}'.format(self.d_num)
            else:               self.d_tag = 'D{:05d}'.format(self.d_num)
        # Store a name for the dataset
        self.d_name = 'Dataset-{!s}'.format(self.d_tag)

        # Output Directories
        self.output_handler = None

        ########################################################

        # Store filenames
        self._pdb_file = os.path.abspath(pdb_filename)
        self._mtz_file = os.path.abspath(mtz_filename)

        # PDB Objects
        self._structure = self.new_structure()

        ########################################################

        # Scaled structure factors
        self.unscaled_sfs = None
        self.scaled_sfs = None
        self.native_map = None
        self.morphed_map = None

        # Initialise other variables
        self.unit_cell = None
        self.space_group = None

        # Single matrix - global alignment
        self._global_rt_transform = None
        # Multiple matrices - local alignment
        self._local_rt_transforms = None

        ########################################################

        # Map uncertainties
        self._map_uncertainty = None

        # Z-Map statistics
        self.z_map_stats = {}

        ########################################################

        # Map of the clusters in the dataset
        self.raw_cluster_hits = None
        self.clustered_hits = None

        # MANUAL VARIABLES - ONLY TO BE ACCESSED DURING TESTING
        self.number_of_sig_points = 0

        ########################################################

    def initialise_output_directory(self, outputdir):
        """Initialise a dataset output directory"""
        # Create a file and directory organiser
        self.output_handler = output_file_object(rootdir=easy_directory(outputdir))

    def get_pdb_filename(self):
        return self._pdb_file
    def get_mtz_filename(self):
        return self._mtz_file

    def get_input(self):
        return self._structure.input
    def get_hierarchy(self):
        return self._structure.hierarchy

    def new_structure(self):
        """Generate a new copy of the input-hierarchy pair, from the pdb file"""
        return iotbx.pdb.hierarchy.input(file_name=self._pdb_file)

    def reflection_data(self):
        """Return an object containing the reflection data"""
        return iotbx.mtz.object(self._mtz_file)

    def get_resolution(self):
        # TODO REMOVE - POINTLESS WRAPPER TODO
        return self.reflection_data().max_min_resolution()[1]

    def get_heavy_atom_sites(self):
        xray_structure = self.get_input().xray_structure_simple()
        return xray_structure.sites_cart().select(xray_structure.heavy_selection())
    def get_calpha_sites(self):
        xray_structure = self.get_input().xray_structure_simple()
        return xray_structure.sites_cart().select(xray_structure.backbone_selection(atom_names=['CA']))
    def get_backbone_sites(self):
        xray_structure = self.get_input().xray_structure_simple()
        return xray_structure.sites_cart().select(xray_structure.backbone_selection())
    def get_global_alignment_transform(self):
        return self._global_rt_transform
    def get_local_alignment_transforms(self):
        return self._local_rt_transforms

    def set_map_uncertainty(self, sigma):
        self._map_uncertainty = sigma
    def get_map_uncertainty(self):
        return self._map_uncertainty

    def get_pickle_copy(self):
        """Get copy of self that can be pickled - some cctbx objects cannot be pickled..."""
        return self

    def get_structure_summary(self):
        return structure_summary(hierarchy=self.get_hierarchy())

    def align_to_reference(self, ref_handler, method, sites='calpha'):
        """Calculate the rotation and translation needed to align one structure to another"""
        assert sites in ['calpha','backbone']
        assert method in ['both','local','global'], 'METHOD NOT DEFINED: {!s}'.format(method)

        if method == 'global' or method == 'both':
            if sites == 'calpha':
                my_sites = self.get_calpha_sites()
                ref_sites = ref_handler.get_calpha_sites()
            elif sites == 'backbone':
                my_sites = self.get_backbone_sites()
                ref_sites = ref_handler.get_backbone_sites()
            assert len(my_sites) == len(ref_sites)
            self._global_rt_transform = superpose.least_squares_fit(reference_sites=ref_sites, other_sites=my_sites)
        if method == 'local' or method == 'both':
            self._local_rt_transforms = perform_flexible_alignment(mov_hierarchy=self.get_hierarchy(), ref_hierarchy=ref_handler.get_hierarchy())

    def transform_from_reference(self, points, method, point_mappings=None):
        """Use alignment to map from reference frame to our frame"""
        assert method in ['local','global'], 'METHOD NOT DEFINED: {!s}'.format(method)
        if method == 'global':
            # Simple - use one rt to transform all of the points
            return self.get_global_alignment_transform().rt().inverse() * points
        elif method == 'local':
            # Use point mappings to tell us which residue rt to use for each point
            assert point_mappings, 'NO MAPPINGS GIVEN BETWEEN POINTS AND RTS'
            assert len(points) == len(point_mappings), 'POINT MAPPINGS NOT THE SAME LENGTH AS POINTS'

            # Get the set of labels to transform in groups
            lab_set = sorted(list(set(point_mappings)))
            all_idxs = []
            # Initialise output array
            rt_points = numpy.zeros(len(points), dtype=[('x',float),('y',float),('z',float)])
            for r_lab in lab_set:
                # Extract the idxs and points of points with this label in the mapping
                lab_idxs, lab_points = zip(*[(i, points[i]) for i, i_lab in enumerate(point_mappings) if r_lab==i_lab])
                # Transform all of these points at once
                lab_rt_points = self.get_local_alignment_transforms()[r_lab].rt().inverse() * flex.vec3_double(lab_points)
                # Populate the array at the appropriate place
                rt_points.put(lab_idxs, lab_rt_points)
                all_idxs.extend(lab_idxs)

            assert len(list(set(all_idxs))) == len(points)

            return flex.vec3_double(rt_points)

    def transform_to_reference(self, points, method, point_mappings=None):
        """Use alignment to map to reference frame from our frame"""
        assert method in ['local','global'], 'METHOD NOT DEFINED: {!s}'.format(method)
        if method == 'global':
            # Simple - use one rt to transform all of the points
            return self.get_global_alignment_transform().rt() * points
        elif method == 'local':
            if not point_mappings:
                # Use point mappings to tell us which residue rt to use for each point
                point_mappings = self.find_nearest_calpha(points)

            # Get the set of labels to transform in groups
            lab_set = sorted(list(set(point_mappings)))
            all_idxs = []
            # Initialise output array
            rt_points = numpy.zeros(len(points), dtype=[('x',float),('y',float),('z',float)])
            for r_lab in lab_set:
                # Extract the idxs and points of points with this label in the mapping
                lab_idxs, lab_points = zip(*[(i, points[i]) for i, i_lab in enumerate(point_mappings) if r_lab==i_lab])
                # Transform all of these points at once
                lab_rt_points = self.get_local_alignment_transforms()[r_lab].rt() * flex.vec3_double(lab_points)
                # Populate the array at the appropriate place
                rt_points.put(lab_idxs, lab_rt_points)
                all_idxs.extend(lab_idxs)

            assert len(list(set(all_idxs))) == len(points)

            # Initial brute force method
            #rt_points = [self.get_local_alignment_transforms()[r_lab].rt() * p for r_lab, p in zip(point_mappings, points)]

            return flex.vec3_double(rt_points)

    def find_nearest_calpha(self, points):
        """Returns the labels of the nearest calpha for each of the given points"""

        calpha_hierarchy = self.get_hierarchy().select(self.get_hierarchy().atom_selection_cache().selection('pepnames and name CA'))
        atom_sites, atom_labels = zip(*[(a.xyz, (a.chain_id, a.resid())) for a in calpha_hierarchy.atoms_with_labels()])

        tree = spatial.KDTree(data=atom_sites)
        nn_dists, nn_groups = tree.query(points)
        return [atom_labels[i] for i in nn_groups]

class reference_dataset_handler(dataset_handler):
    _origin_shift = (0,0,0)
    _binning = None

    def set_origin_shift(self, origin_shift):
        self._origin_shift = origin_shift
    def get_origin_shift(self):
        return self._origin_shift

class identical_structure_ensemble(object):
    """Class for collating and comparing multiple observations of the same structure"""

    def __init__(self, ref_hierarchy):
        """Initialise the comparison object"""

        self.data = multiple_data_collection()
        self.ref_hierarchy = ref_hierarchy

        # Set the collection ids to the residue labels
        residue_labels = [(rg.parent().id, rg.resid()) for rg in ref_hierarchy.residue_groups()]
        self.data.set_collection_ids(residue_labels)
        # List of just chain ids
        chain_ids = [rg.parent().id for rg in ref_hierarchy.residue_groups()]
        self.data.add_collection_info(info_name='chain ids', info_values=chain_ids)
        # List of just residue ids
        residue_ids = [rg.resid() for rg in ref_hierarchy.residue_groups()]
        self.data.add_collection_info(info_name='residue ids', info_values=residue_ids)
        # List of which residues are protein
        is_protein = [rg.parent().is_protein() for rg in ref_hierarchy.residue_groups()]
        self.data.add_collection_info(info_name='is protein', info_values=is_protein)
        # List of which residues are in multiple conformations
        has_conformers = [rg.have_conformers() for rg in ref_hierarchy.residue_groups()]
        self.data.add_collection_info(info_name='has conformers', info_values=has_conformers)
        # List of number of conformers per residue
        num_conformers = [len(rg.conformers()) for rg in ref_hierarchy.residue_groups()]
        self.data.add_collection_info(info_name='number of conformers', info_values=num_conformers)
        # List of number of unique resnames for each residue - XXX I'm not sure when this can be more than one?
        residue_names = [[s for s in rg.unique_resnames()] for rg in ref_hierarchy.residue_groups()]
        self.data.add_collection_info(info_name='residue names', info_values=residue_names)
        # List of residues types (amino, water, other, etc)
        residue_types = [iotbx.pdb.common_residue_names_get_class(rg.unique_resnames()[0]) for rg in ref_hierarchy.residue_groups()]
        self.data.add_collection_info(info_name='residue types', info_values=residue_types)
        # List of atom ids for each residue
        residue_atom_labels = [[a.id_str() for a in rg.atoms()] for rg in ref_hierarchy.residue_groups()]
        self.data.add_collection_info(info_name='atom labels', info_values=residue_atom_labels)

    def add_structures(self, new_hierarchies, hierarchy_ids=None, verbose=True):
        """Add hierarchies to the analysis object - Iterate through the residues and extract values from across the datasets"""

        # If no ids given, simply assign indexing values
        if not hierarchy_ids: hierarchy_ids = range(1, len(new_hierarchies)+1)

        # Report string to be returned
        report_string = []

        backbone_atom_names = [" CA ", " C  ", " O  ", " N  "]

        # Extract the reference residue types
        ref_residue_types = self.data.get_collection_info('residue types')
        ref_residue_names = self.data.get_collection_info('residue names')
        ref_chain_ids =     self.data.get_collection_info('chain ids')

        print('===================================>>>')
        print('Calculating B-factors for residues across the datasets')
        print('===================================>>>')

        # Atomic selection cache for the reference structure - should be the same across the structures
        selection_cache = self.ref_hierarchy.atom_selection_cache()

        # String of 1-letter codes for the protein
        amino_sequence = ''
        current_chain = ref_chain_ids[self.data.collection_ids[0]]

        for res_id in self.data.collection_ids:

            # Initialise a new object to hold the data for one residue
            res_collection = data_collection()
            res_collection.set_entry_ids(hierarchy_ids)

            # Boolean mask for the selected residue
            res_bool_selection = selection_cache.selection("chain '{!s}' and resid '{!s}'".format(*res_id))

            # Information to collect for the residues
            backbone_mean_bfactor = []
            sidechain_mean_bfactor = []
            total_mean_bfactor = []

            # Switches for tracking and reporting new chains
            if current_chain != ref_chain_ids[res_id]:
                report_string.append('\rAnalysing Chain {!s}: {!s}'.format(current_chain, amino_sequence))
                print(report_string[-1])
                report_string.append('\rChain {!s} Analysed.'.format(current_chain))
                print(report_string[-1])
                # Reset sequence and update chain
                amino_sequence = ''
                current_chain = ref_chain_ids[res_id]

            try:
                res_class = iotbx.pdb.common_residue_names_get_class(ref_residue_names[res_id][0])
                if res_class == 'common_amino_acid':
                    amino_sequence += iotbx.pdb.amino_acid_codes.one_letter_given_three_letter[ref_residue_names[res_id][0]]
                elif res_class == 'common_rna_dna':
                    amino_sequence += '<DNA/RNA>'
                elif res_class == 'ccp4_mon_lib_rna_dna':
                    amino_sequence += '<DNA/RNA>'
                elif res_class == 'common_water':
                    amino_sequence += 'o'
                elif res_class == 'common_small_molecule':
                    amino_sequence += '<MOL>'
                elif res_class == 'common_element':
                    amino_sequence += 'i'
                elif res_class == 'other':
                    amino_sequence += 'X'
                else:
                    raise Exception()
            except:
                amino_sequence += '?'

            if len(amino_sequence) >= 50:
                report_string.append('\rAnalysing Chain {!s}: {!s}...'.format(current_chain, amino_sequence))
                print(report_string[-1])
                amino_sequence = ''
            else:
                print '\rAnalysing Chain {!s}: {!s}-->'.format(current_chain, amino_sequence),; sys.stdout.flush()

            # Iterate through the hierarchies and extract values for this residue
            for i_hierarchy in new_hierarchies:

                # Select the atoms for this residue
                new_root = i_hierarchy.select(res_bool_selection)

                # Calculate the sidechain and backbone b-factors for amino acids
                if ref_residue_types[res_id] == 'common_amino_acid':

                    # Pull out the backbone atoms
                    backbone_atoms = [at for at in new_root.atoms() if at.name in backbone_atom_names]
                    if not backbone_atoms: backbone_mean_bfactor.append(numpy.nan)
                    else:                  backbone_mean_bfactor.append(flex.mean(flex.double([at.b for at in backbone_atoms])))

                    # Pull out the sidechain atoms
                    sidechain_atoms = [at for at in new_root.atoms() if at.name not in backbone_atom_names]
                    if not sidechain_atoms: sidechain_mean_bfactor.append(numpy.nan)
                    else:                   sidechain_mean_bfactor.append(flex.mean(flex.double([at.b for at in sidechain_atoms])))

                else:
                    # If not an amino acid, just append None
                    backbone_mean_bfactor.append(numpy.nan)
                    sidechain_mean_bfactor.append(numpy.nan)

                # Calculate the mean b-factor for the whole residue
                total_atoms = [at for at in new_root.atoms()]
                if not total_atoms: total_mean_bfactor.append(numpy.nan)
                else:               total_mean_bfactor.append(flex.mean(flex.double([at.b for at in total_atoms])))

            res_collection.add_data(data_name='backbone_mean_bfactor', data_values=backbone_mean_bfactor)
            res_collection.add_data(data_name='sidechain_mean_bfactor', data_values=sidechain_mean_bfactor)
            res_collection.add_data(data_name='total_mean_bfactor', data_values=total_mean_bfactor)

            self.data.add_collection(collection_id=res_id, collection=res_collection)

        report_string.append('\rAnalysing Chain {!s}: {!s}   '.format(current_chain, amino_sequence))
        print(report_string[-1])
        report_string.append('\rChain {!s} Analysed.'.format(current_chain))
        print(report_string[-1])

        return '\n'.join(report_string)

class map_list(object):
    _initialized = False
    def __init__(self, map_names):
        for m in map_names:
            self.__dict__[m] = None
        self._initialized = True

    def __setattr__(self, name, value):
        if self._initialized and (not hasattr(self, name)):
            raise AttributeError('Cannot set new attributes after initialisation: {!s}'.format(name))
        object.__setattr__(self, name, value)

class multi_dataset_analyser(object):
    """Class for the processing of datasets from a fragment soaking campaign"""

    _version = PANDDA_VERSION

    _structure_mask_names = STRUCTURE_MASK_NAMES
    _crystal_mask_names   = CRYSTAL_MASK_NAMES
    _reject_mask_names    = REJECT_MASK_NAMES
    _flag_mask_names      = FLAG_MASK_NAMES
    _custom_mask_names    = ['selected for analysis']
    _all_mask_names = _structure_mask_names + _crystal_mask_names + _reject_mask_names + _flag_mask_names + _custom_mask_names

    _map_names =    [   'mean_map',
                        'stds_map',
                        'adj_stds_map',
                        'skew_map',
                        'kurt_map',
                        'bimo_map'
                    ]

    phil = {}

    def __init__(self, args=None):
#                    outdir='./pandda', datadir='./Processing', pdb_style='*/refine.pdb', mtz_style='*/refine.mtz',
#                    ref_pdb='./reference.pdb', ref_mtz='./reference.mtz', run_mol_subst=False,
#                    verbose=True, keep_maps_in_memory=False, maxmemory=25):
        """Class for the processing of datasets from a fragment soaking campaign"""

        # Allow the program to pull from the command line if no arguments are given
        if args == None:
            args = sys.argv

        # Process the input arguments and convert to phil
        self.cmds = args
        self.master_phil, self.working_phil, self.unused_args = self.process_input_args(args)
        self.args = self.working_phil.extract().pandda
        self.params = self.args.params

        # Validate the processed parameters
        self._validate_parameters()

        # ===============================================================================>
        # INPUT FILES STUFF
        # ===============================================================================>

        self.data_dirs = os.path.abspath(self.args.input.data_dirs)

        # ===============================================================================>
        # OUTPUT FILES STUFF
        # ===============================================================================>

        self.outdir = easy_directory(os.path.abspath(self.args.output.outdir))

        self.log_file = os.path.join(self.outdir, 'pandda.log')
        self._log = ''

        # Set up the output folders and files
        self._run_directory_setup()

        # ===============================================================================>
        # SETTINGS STUFF
        # ===============================================================================>

        self._high_resolution = None
        self._low_resolution = None

        # Current and maximum sizes of the pandda in memory
        self._pandda_size = [('Zero',0)]

        # ===============================================================================>
        # DATA AND MAPS STUFF
        # ===============================================================================>

        # File names
        self._input_files = []
        # Dataset Objects
        self.datasets = dataset_handler_list()
        # Reference Objects
        self._ref_dataset_index = None
        self._ref_dataset = None
        self._ref_grid = None
        self._ref_map = None
        # Distance to displace the reference model by so that its bounding box sits on the origin
        self._ref_origin_shift = None

        # Map Statistics
        self.stat_maps = map_list(self._map_names)

        # ===============================================================================>
        # ANALYSIS OBJECTS
        # ===============================================================================>

        # Structure summary Object
        self.ensemble_summary = None
        # Dataset Summary Object
        self.datasets_summary = None



        # TODO RENAME REMOVE
        # Map and Dataset statistics
        self._map_observations = data_collection()

        # Summaries of detected points
        self._cluster_summary = data_collection()

        # Average Structure (and masks)
        self._average_calpha_sites = None
        self._residue_deviation_masks = None

        self.crystal_contact_generators = None

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

    def process_input_args(self, args):
        """Process the input arguments"""

        # Copy the args so that we can remove items from the list without affecting sys.argv etc
        args = copy.copy(args)

        assert isinstance(args, list), 'INPUT ARGUMENTS MUST BE A LIST'
        assert len(args) != 0, 'NO INPUT ARGUMENTS GIVEN'

        # Read in the master phil
        master_phil = phil.parse(pandda_phil_def)
        cmd_interpr = master_phil.command_line_argument_interpreter(home_scope="pandda")

        # Look for any eff files in the args
        eff_files = [f for f in args if (f.endswith('.eff') and os.path.isfile(f))]
        # Remove them from the original lists
        [args.remove(f) for f in eff_files]
        # Parse the 'eff' files - these should contain phils
        eff_sources = [phil.parse(open(f, 'r').read()) for f in eff_files]

        # Look for cmd line arguments
        cmd_line_args = [a for a in args if (('=' in a) and not os.path.exists(a))]
        # Remove them from the original lists
        [args.remove(a) for a in cmd_line_args]
        # Parse these arguments
        cmd_sources = [cmd_interpr.process(arg=a) for a in cmd_line_args]

        # Combine the phils
        working_phil = master_phil.fetch(
                                        sources=eff_sources+cmd_sources
                                        )

        return master_phil, working_phil, args

    def run_pandda_init(self):
        """Set up the pandda"""

        # Print logo to log
        self.log(PANDDA_TEXT.format(self._version), True)

        self.log('===================================>>>', True)
        self.log('PROCESSED ARGUMENTS', True)
        self.log('===================================>>>', True)
        self.log(self.working_phil.as_str(), True)

        # Write the used parameters to file
        with open(self.output_handler.get_file('pandda_settings'), 'w') as out_file:
            out_file.write( '\n'.join([ '# ' + 'Arguments',
                                        '# ',
                                        '# ' + '\n# '.join(self.cmds),
                                        '',
                                        '# Used Settings:',
                                        '',
                                        self.working_phil.as_str() ]))

        if self.unused_args:
            self.log('===================================>>>', True)
            self.log('UNUSED ARGUMENTS', True)
            self.log('===================================>>>', True)
            self.log('\n'.join(self.unused_args), True)

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

        # ===============================================================================>
        # LOOK FOR MATPLOTLIB TO SEE IF WE CAN GENERATE GRAPHS
        # ===============================================================================>

        try:
            import matplotlib
            # Setup so that we can write without a display connected
            matplotlib.interactive(0)
            default_backend, validate_function = matplotlib.defaultParams['backend']
            self.log('===================================>>>')
            self.log('MATPLOTLIB LOADED. Using Backend: {!s}'.format(default_backend))
            from matplotlib import pyplot
            self.log('PYPLOT loaded successfully')
        except:
            self.log('===================================>>>', True)
            self.log('>> COULD NOT IMPORT MATPLOTLIB. CANNOT GENERATE GRAPHS.', True)

        # ===============================================================================>
        # LOG THE START TIME
        # ===============================================================================>

        # Store start time and print to log
        self._init_time = time.time()
        self.log('===================================>>>', True)
        self.log('Analysis Started: {!s}'.format(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime(self._init_time))), True)

    def _validate_parameters(self):
        """Validate and preprocess the loaded parameters"""

        p = self.args

        # Input
        assert p.input.data_dirs is not None, 'pandda.input.data_dirs IS NOT DEFINED'
        assert p.input.pdb_style, 'pandda.input.pdb_style IS NOT DEFINED'
#        assert p.input.mtz_style, 'pandda.input.mtz_style IS NOT DEFINED'
        assert p.input.lig_style, 'pandda.input.lig_style IS NOT DEFINED'
        # Output
        assert p.output.outdir, 'pandda.output.outdir IS NOT DEFINED'

    def _run_directory_setup(self):
        """Initialise the pandda directory system"""

        # Create a file and directory organiser
        self.output_handler = output_file_object(rootdir=self.outdir)

        # Input parameters
        self.output_handler.add_file(file_name='pandda.eff', file_tag='pandda_settings', dir_tag='root')

        # Somewhere to flag the empty directories (failed pre-processing)
        self.output_handler.add_dir(dir_name='empty_directories', dir_tag='empty_directories', top_dir_tag='root')
        # Somewhere to flag the rejected directories (rejected from analysis)
        self.output_handler.add_dir(dir_name='rejected_datasets', dir_tag='rejected_datasets', top_dir_tag='root')
        # Somewhere to store all of the output maps
        self.output_handler.add_dir(dir_name='processed_datasets', dir_tag='processed_datasets', top_dir_tag='root')
        # Somewhere to store the interesting datasets
        self.output_handler.add_dir(dir_name='interesting_datasets', dir_tag='interesting_datasets', top_dir_tag='root')
        # Somewhere to store the noisy datasets
        self.output_handler.add_dir(dir_name='noisy_datasets', dir_tag='noisy_datasets', top_dir_tag='root')
        # Somewhere to store all of the aligned structures
        self.output_handler.add_dir(dir_name='aligned_structures', dir_tag='aligned_structures', top_dir_tag='root')
        # Somewhere to store the analyses/summaries - for me to plot graphs
        self.output_handler.add_dir(dir_name='analyses', dir_tag='analyses', top_dir_tag='root')
        self.output_handler.add_file(file_name='map_summaries.csv', file_tag='map_summaries', dir_tag='analyses')
        self.output_handler.add_file(file_name='statistical_map_summaries.csv', file_tag='stats_map_summaries', dir_tag='analyses')
        self.output_handler.add_file(file_name='dataset_summaries.csv', file_tag='dataset_summaries', dir_tag='analyses')
        self.output_handler.add_file(file_name='blob_site_summaries.csv', file_tag='blob_summaries', dir_tag='analyses')
        self.output_handler.add_file(file_name='point_distributions.csv', file_tag='point_distributions', dir_tag='analyses')
        # Somewhere to store the analysis summaries - for the user
        self.output_handler.add_dir(dir_name='results_summaries', dir_tag='output_summaries', top_dir_tag='root')
        self.output_handler.add_file(file_name='success_summary_table.html', file_tag='summary_table', dir_tag='output_summaries')
        # Somewhere to store the pickled objects
        self.output_handler.add_dir(dir_name='pickled_panddas', dir_tag='pickle', top_dir_tag='root')

        # Reference Structure and Dataset
        self.output_handler.add_dir(dir_name='reference', dir_tag='reference', top_dir_tag='root')
        self.output_handler.add_file(file_name='reference.pdb', file_tag='reference_structure', dir_tag='reference')
        self.output_handler.add_file(file_name='reference.mtz', file_tag='reference_dataset', dir_tag='reference')
        self.output_handler.add_file(file_name='reference.shifted.pdb', file_tag='reference_on_origin', dir_tag='reference')
        self.output_handler.add_file(file_name='reference.symmetry.pdb', file_tag='reference_symmetry', dir_tag='reference')

        # Statistical Maps
        self.output_handler.add_dir(dir_name='statistical_maps', dir_tag='statistical_maps', top_dir_tag='root')
        self.output_handler.add_file(file_name='mean_map.ccp4', file_tag='mean_map', dir_tag='statistical_maps')
        self.output_handler.add_file(file_name='stds_map.ccp4', file_tag='stds_map', dir_tag='statistical_maps')
        self.output_handler.add_file(file_name='stds_adj_map.ccp4', file_tag='stds_adj_map', dir_tag='statistical_maps')
        self.output_handler.add_file(file_name='skew_map.ccp4', file_tag='skew_map', dir_tag='statistical_maps')
        self.output_handler.add_file(file_name='kurt_map.ccp4', file_tag='kurt_map', dir_tag='statistical_maps')
        self.output_handler.add_file(file_name='bimo_map.ccp4', file_tag='bimo_map', dir_tag='statistical_maps')

    def _run_pickle_setup(self):
        """Initialise all of the pickle filenames"""

        # Pickle Handler
        self.pickle_handler = output_file_object(rootdir=self.output_handler.get_dir('pickle'))
        # Pickled Reference Objects
        self.pickle_handler.add_file(file_name='reference_grid.pickle', file_tag='reference_grid')
        # Pickled Datasets
        self.pickle_handler.add_file(file_name='dataset_masks.pickle', file_tag='dataset_masks')
        self.pickle_handler.add_file(file_name='dataset_objs.pickle', file_tag='dataset_objs')
        # Pickled Information
        self.pickle_handler.add_file(file_name='dataset_observations.pickle', file_tag='dataset_observations')
        self.pickle_handler.add_file(file_name='map_observations.pickle', file_tag='map_observations')
        # Pickled Stats
        self.pickle_handler.add_file(file_name='mean_map.pickle', file_tag='mean_map')
        self.pickle_handler.add_file(file_name='stds_map.pickle', file_tag='stds_map')
        self.pickle_handler.add_file(file_name='adj_stds_map.pickle', file_tag='adj_stds_map')
        self.pickle_handler.add_file(file_name='skew_map.pickle', file_tag='skew_map')
        self.pickle_handler.add_file(file_name='kurt_map.pickle', file_tag='kurt_map')
        self.pickle_handler.add_file(file_name='bimo_map.pickle', file_tag='bimo_map')

        # Pickled SELF
        self.pickle_handler.add_file(file_name='my_pandda_dict.pickle', file_tag='my_pandda_dict')

    def load_pickled_objects(self):
        """Loads any pickled objects it finds"""

        self.log('===================================>>>', True)
        self.log('Looking for Pickled Files in Input Directory: {!s}'.format(os.path.relpath(self.pickle_handler.get_dir('root'))), True)
        # Load Reference Grid
        if os.path.exists(self.pickle_handler.get_file('reference_grid')):
            self.set_reference_grid(self.unpickle(self.pickle_handler.get_file('reference_grid')))

        # Load the datasets
        if os.path.exists(self.pickle_handler.get_file('dataset_objs')):
            self.datasets.add(self.unpickle(self.pickle_handler.get_file('dataset_objs')))
            self.update_pandda_size(tag='After Unpickling Dataset Objects')

        # Load Statistical Maps
        if os.path.exists(self.pickle_handler.get_file('mean_map')):
            self.stat_maps.mean_map = self.unpickle(self.pickle_handler.get_file('mean_map'))
        if os.path.exists(self.pickle_handler.get_file('stds_map')):
            self.stat_maps.stds_map = self.unpickle(self.pickle_handler.get_file('stds_map'))
        if os.path.exists(self.pickle_handler.get_file('adj_stds_map')):
            self.stat_maps.adj_stds_map = self.unpickle(self.pickle_handler.get_file('adj_stds_map'))
        if os.path.exists(self.pickle_handler.get_file('skew_map')):
            self.stat_maps.skew_map = self.unpickle(self.pickle_handler.get_file('skew_map'))
        if os.path.exists(self.pickle_handler.get_file('kurt_map')):
            self.stat_maps.kurt_map = self.unpickle(self.pickle_handler.get_file('kurt_map'))
        if os.path.exists(self.pickle_handler.get_file('bimo_map')):
            self.stat_maps.bimo_map = self.unpickle(self.pickle_handler.get_file('bimo_map'))
            self.update_pandda_size(tag='After Unpickling Statistical Maps')

    def pickle_the_pandda(self, components=[], all=False):
        """Pickles it's major components for quick loading..."""

        if all == True:
            self.log('===================================>>>', True)
            self.log('Pickling the PanDDA', True)
        elif not components:
            self.log('===================================>>>', True)
            self.log('Pickling NOTHING', True)
            return
        else:
            self.log('===================================>>>', True)
            self.log('Selective Pickling: {!s}'.format(', '.join(components)), True)

        if all or ('grid' in components):
            self.log('===================================>>>')
            self.log('Pickling Reference Grid')
            if self._ref_grid is not None:
                self.pickle(pickle_file=self.pickle_handler.get_file('reference_grid'), pickle_object=self.reference_grid(), force=False)

        if all or ('datasets' in components):
            self.log('===================================>>>')
            self.log('Pickling Datasets')
            if self.datasets:
                self.pickle(pickle_file=self.pickle_handler.get_file('dataset_objs'), pickle_object=[d.get_pickle_copy() for d in self.datasets.all()], force=True)

        if all or ('statistical_maps' in components):
            self.log('===================================>>>')
            self.log('Pickling Statistical Maps')
            if self.stat_maps.mean_map is not None:
                self.pickle(pickle_file=self.pickle_handler.get_file('mean_map'), pickle_object=self.stat_maps.mean_map, force=True)
            if self.stat_maps.stds_map is not None:
                self.pickle(pickle_file=self.pickle_handler.get_file('stds_map'), pickle_object=self.stat_maps.stds_map, force=True)
            if self.stat_maps.adj_stds_map is not None:
                self.pickle(pickle_file=self.pickle_handler.get_file('adj_stds_map'), pickle_object=self.stat_maps.adj_stds_map, force=True)
            if self.stat_maps.skew_map is not None:
                self.pickle(pickle_file=self.pickle_handler.get_file('skew_map'), pickle_object=self.stat_maps.skew_map, force=True)
            if self.stat_maps.kurt_map is not None:
                self.pickle(pickle_file=self.pickle_handler.get_file('kurt_map'), pickle_object=self.stat_maps.kurt_map, force=True)
            if self.stat_maps.bimo_map is not None:
                self.pickle(pickle_file=self.pickle_handler.get_file('bimo_map'), pickle_object=self.stat_maps.bimo_map, force=True)

    def exit(self):
        """Exit the PANDDA, record runtime etc..."""
        self.log('===================================>>>', True)
        self.log('...FINISHED!...', True)
        self.log('===================================>>>', True)
        self._finish_time = time.time()
        self.log('Runtime: {!s}'.format(time.strftime("%H hours:%M minutes:%S seconds", time.gmtime(self._finish_time - self._init_time))))

        # Pickle myself
        self.log('===================================>>>', True)
        self.log('Pickling the PANDDA Results')
        pandda_dict =  {
                        'combined_clusters' : self.process_z_value_clusters()
                    }
        self.pickle(pickle_file=self.pickle_handler.get_file('my_pandda_dict'), pickle_object=pandda_dict, force=True)

    def log(self, message, show=False):
        """Log message to file, and mirror to stdout if verbose or force_print"""
        if not isinstance(message, str):    message = str(message)
        # Print to stdout
        if show or self.args.settings.verbose:
            print(message)
        # Remove \r from message as this spoils the log (^Ms)
        message = message.replace('\r','')
        # Store in internal string
        self._log = self._log + message + '\n'
        # Write to file
        open(self.log_file, 'a').write(message+'\n')

    def print_log(self):
        print(self._log)

    def get_pandda_size(self):
        """Returns the history of the memory consumption of the PANDDA"""
        return self._pandda_size
    def update_pandda_size(self, tag):
        pandda_size = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1000
        assert pandda_size < self.args.settings.max_memory*1024**3, 'PANDDA HAS EXCEEDED THE MAXIMUM AMOUNT OF ALLOWED MEMORY'
        self._pandda_size.append((tag, pandda_size))
        from humanize import naturalsize
        self.log(tag+': '+naturalsize(pandda_size, binary=True))

    def is_new_pandda(self):
        """Is this the first time the program has been run?"""
        return self._new_pandda

# XXX MAYBE KEEP THESE
    def set_high_resolution(self, res):
        self._high_resolution = res
    def get_high_resolution(self):
        return self._high_resolution

    def set_low_resolution(self, res):
        self._low_resolution = res
    def get_low_resolution(self):
        return self._low_resolution

    def set_cut_resolution(self, res):
        self._cut_resolution = res
    def get_cut_resolution(self):
        return self._cut_resolution
# XXX MAYBE KEEP THESE

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

    def new_files(self):
        """Get all of the files that were added on this run"""
        return self._input_files

    def get_map_observations(self):
        """Distributions of different map variables across the datasets"""
        return self._map_observations
    def get_cluster_summary(self):
        """Summaries of difference clusters across the datasets"""
        return self._cluster_summary

    def get_calpha_average_sites(self):
        return self._average_calpha_sites
    def get_calpha_deviation_masks(self):
        return self._residue_deviation_masks

    def initialise_analysis(self):
        """Add blank masks to the mask objects, based on how many datasets have been loaded"""

        self.log('===================================>>>', True)
        self.log('Initialising Dataset Masks.', True)

        # Set the dataset mask lengths and ids (dataset tags)
        self.datasets.all_masks().set_mask_length(mask_length=self.datasets.size())
        self.datasets.all_masks().set_entry_ids(entry_ids=[d.d_tag for d in self.datasets.all()])
        # Initialise standard blank masks
        for mask_name in self._all_mask_names:
            self.datasets.all_masks().add_mask(mask_name=mask_name, mask=[False]*self.datasets.size())

        self.log('===================================>>>', True)
        self.log('Initialising Dataset Observations.', True)

        # Set the dataset observations lengths
        self.datasets_summary = data_collection()
        self.datasets_summary.set_data_length(self.datasets.size())
        self.datasets_summary.set_entry_ids([d.d_tag for d in self.datasets.all()])

        # Set the map observation lengths
        self.get_map_observations().set_data_length(self.datasets.size())
        self.get_map_observations().set_entry_ids([d.d_tag for d in self.datasets.all()])

        # Initialise blank data
        self.get_map_observations().add_empty_data(data_name='masked_map_mean')
        self.get_map_observations().add_empty_data(data_name='masked_map_std')
        self.get_map_observations().add_empty_data(data_name='masked_map_uncertainty')
        self.get_map_observations().add_empty_data(data_name='z_map_mean')
        self.get_map_observations().add_empty_data(data_name='z_map_std')
        self.get_map_observations().add_empty_data(data_name='z_map_skew')
        self.get_map_observations().add_empty_data(data_name='z_map_kurtosis')

        self.log('===================================>>>', True)
        self.log('Initialising Timings.', True)

        # Set the dataset observations lengths
        self.timings.set_data_length(self.datasets.size())
        self.timings.set_entry_ids([d.d_tag for d in self.datasets.all()])

    def load_reference_dataset(self, ref_pdb, ref_mtz):
        """Set the reference dataset, to which all other datasets will be aligned and scaled"""

        self.log('===================================>>>', True)
        self.log('Loading Reference Dataset: {!s}'.format(ref_mtz), True)

        self._ref_dataset = dataset_handler(dataset_number=-1, pdb_filename=ref_pdb, mtz_filename=ref_mtz, dataset_tag='reference')

        if not os.path.exists(self.output_handler.get_file('reference_structure')):
            rel_symlink(orig=ref_pdb, link=self.output_handler.get_file('reference_structure'))
        if not os.path.exists(self.output_handler.get_file('reference_dataset')):
            rel_symlink(orig=ref_mtz, link=self.output_handler.get_file('reference_dataset'))

        # Calculate the shift required to move the reference structure into the positive quadrant
        self.set_reference_origin_shift(tuple(flex.double(3, self.params.maps.border_padding) - flex.double(self.reference_dataset().get_input().atoms().extract_xyz().min())))
        self.log('Origin Shift for reference structure: {!s}'.format(tuple([round(s,3) for s in self.get_reference_origin_shift()])))
        # Shift the reference structure by this amount so that it is aligned with the reference grid
        ref_hierarchy = self.reference_dataset().get_hierarchy()
        ref_hierarchy.atoms().set_xyz(ref_hierarchy.atoms().extract_xyz() + self.get_reference_origin_shift())

        if not os.path.exists(self.output_handler.get_file('reference_on_origin')):
            ref_hierarchy.write_pdb_file(self.output_handler.get_file('reference_on_origin'))

        self.log('===================================>>>', True)
        self.log('Initialising Ensemble Summary', True)

        # Initialise the structural ensemble object
        self.ensemble_summary = identical_structure_ensemble(ref_hierarchy=self.reference_dataset().new_structure().hierarchy)

        return self.reference_dataset()

    def create_reference_grid(self, grid_spacing, expand_to_origin, buffer=0):
        """Create a grid over the reference protein"""

        self.log('===================================>>>', True)
        self.log('Creating Reference Grid', True)

        sites_cart = self.reference_dataset().get_input().atoms().extract_xyz()

        assert (flex.vec3_double(sorted(self.reference_dataset().get_input().atoms().extract_xyz())) -
                flex.vec3_double(sorted(self.reference_dataset().get_hierarchy().atoms().extract_xyz()))).dot().norm() == 0.0, 'EH? Coordinates should be the same?'

        self._ref_grid = grid_handler(verbose=self.args.settings.verbose)
        self._ref_grid.set_grid_spacing(spacing=grid_spacing)
        self._ref_grid.set_cart_extent(cart_min=tuple([s-buffer for s in sites_cart.min()]), cart_max=tuple([s+buffer for s in sites_cart.max()]))
        self._ref_grid.create_cartesian_grid(expand_to_origin=expand_to_origin)
        self.log(self._ref_grid.summary())

        self.log('===================================>>>', True)
        self.log('Partitioning Reference Grid', True)

        # Pull out the calphas
        calpha_hierarchy = self.reference_dataset().get_hierarchy().select(self.reference_dataset().get_hierarchy().atom_selection_cache().selection('pepnames and name CA'))

        if 1 or self.params.alignment.method == 'local':
            t1 = time.time()

            # Calculate the nearest residue for each point on the grid
            self.reference_grid().create_grid_partition(atomic_hierarchy=calpha_hierarchy)
            # Partition Grid
            self.reference_grid().partition().partition_grid(cpus=self.args.settings.cpus)

            t2 = time.time()
            print('> GRID PARTITIONING > Time Taken: {!s} seconds'.format(int(t2-t1)))

    def mask_reference_grid(self, d_handler):
        """Create masks for the reference grid based on distances from atoms in the reference structure"""

        self.log('===================================>>>', True)
        self.log('Masking Reference Grid', True)
        self.log('NOT YET', True)

        # Write protein masked map
        mask_map_file = self.output_handler.get_file('reference_dataset').replace('.mtz','.totalmask.ccp4')
        self.write_array_to_map(output_file = mask_map_file,
                                map_data = flex.double(self.reference_grid().global_mask().total_mask_binary()),
                                grid_size = self.reference_grid().grid_size(),
                                grid_spacing = self.reference_grid().grid_spacing())

        # Write symmetry masked map
        mask_map_file = self.output_handler.get_file('reference_dataset').replace('.mtz','.symmask.ccp4')
        self.write_array_to_map(output_file = mask_map_file,
                                map_data = flex.double(self.reference_grid().symmetry_mask().total_mask_binary()),
                                grid_size = self.reference_grid().grid_size(),
                                grid_spacing = self.reference_grid().grid_spacing())

#    def mask_resampled_reference_grid(self):
#        """Using the local and global masks, mask the resampled grid points"""
#
#        self.log('===================================>>>', True)
#        self.log('Masking Resampled Grid (Global and Local Mask)', True)
#
#        # Get the grid points that are not masked, and points in the buffer
#        grid_indexer = self.reference_grid().grid_indexer()
#        resampled_points = self.reference_grid().resampled_grid_points()
#        buffer_mask_binary = self.reference_grid().buffer_mask_binary()
#        global_mask_binary = self.reference_grid().global_mask().total_mask_binary()
#
#        self.log('===================================>>>')
#        self.log('Resampled Grid Size (3D): {!s}'.format(self.reference_grid().resampled_grid_size()))
#        self.log('Resampled Grid Size (1D): {!s}'.format(len(resampled_points)))
#        self.log('===================================>>>')
#
#        # Remove points using the protein mask, and the mask around the edge of the grid
#        masked_resampled_points_1 = resampled_points
#        self.log('Filtering with Buffer Mask... (Edge of Cell)')
#        masked_resampled_points_2 = [p for p in masked_resampled_points_1 if buffer_mask_binary[grid_indexer(p)] == 0]
#        self.log('Filtered Points: {!s}'.format(len(masked_resampled_points_2)))
#        self.log('Filtering with Global Mask... (Protein)')
#        masked_resampled_points_3 = [p for p in masked_resampled_points_2 if global_mask_binary[grid_indexer(p)] == 1]
#        self.log('Filtered Points: {!s}'.format(len(masked_resampled_points_3)))
#        masked_resampled_points = masked_resampled_points_3
#
#        # Store points in the reference grid object
#        self.reference_grid().set_masked_grid_points(masked_resampled_points)
#
#        # Calculate the resampled mask on the original grid
#        full_binary_mask = numpy.zeros(self.reference_grid().grid_size_1d(), int)
#        [full_binary_mask.put(grid_indexer(gp), 1) for gp in masked_resampled_points]
#        # Downsample to the resampled grid
#        resampled_binary_mask = [full_binary_mask[grid_indexer(gp)] for gp in resampled_points]

    def build_input_list(self):
        """Builds a list of input files from the command line arguments passed"""

        self.log('===================================>>>', True)
        self.log('Building List of Datasets')

        dir_style = self.args.input.data_dirs.strip('./')
        pdb_style = self.args.input.pdb_style.lstrip('/')
        if self.args.input.mtz_style:
            mtz_style = self.args.input.mtz_style.lstrip('/')
        else:
            mtz_style = pdb_style.replace('.pdb','.mtz')

        # Datasets that are already added
        already_added  = [(d.get_pdb_filename(), d.get_mtz_filename()) for d in self.datasets.all()]
        new_files = []
        empty_directories = []

        for dir in sorted(glob.glob(self.data_dirs)):
            pdb_files = glob.glob(os.path.join(dir, pdb_style))
            mtz_files = glob.glob(os.path.join(dir, mtz_style))
            if not (pdb_files and mtz_files):
                print('EMPTY DIRECTORY: {!s}'.format(dir))
                empty_directories.append(dir)
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

        # Get the list of already linked empty_directories
        empty_dir_prefix = 'Dir_'
        link_old_empty_dirs = glob.glob(os.path.join(self.output_handler.get_dir('empty_directories'),'*'))
        real_old_empty_dirs = [os.path.realpath(p) for p in link_old_empty_dirs]

        # Pull out the highest current idx
        high_idx = max([0]+[int(os.path.basename(v).strip(empty_dir_prefix)) for v in link_old_empty_dirs])
        assert high_idx == len(link_old_empty_dirs), 'Numbering of empty directories is not consecutive'

        # Link the empty directories into the same directory
        for dir in empty_directories:
            # Already linked
            if dir in real_old_empty_dirs:
                continue
            # Increment counter
            high_idx += 1
            # Create new dir
            empty_dir = os.path.join(self.output_handler.get_dir('empty_directories'), empty_dir_prefix+'{:05d}'.format(high_idx))
            if not os.path.exists(empty_dir):
                os.symlink(dir, empty_dir)
            else:
                raise Exception('THIS DIRECTORY SHOULD NOT EXIST: {!s}'.format(empty_dir))

        self.log('{!s} EMPTY DIRECTORIES FOUND'.format(high_idx), True)

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
        assert self.datasets.all() == []

        # Generate arg_list for loading
        arg_list = [{'dataset_number':d_num, 'pdb_filename':pdb, 'mtz_filename':mtz, 'dataset_tag':tag} for d_num, (pdb, mtz, tag) in enumerate(self.new_files())]

        start = time.time()
        self.log('===================================>>>', True)
        print 'Loading Datasets... (using {!s} cores)'.format(self.args.settings.cpus)
        loaded_datasets = easy_mp.pool_map(fixed_func=map_func, args=arg_list, processes=self.args.settings.cpus)
        finish = time.time()
        print('> Time Taken: {!s} seconds'.format(int(finish-start)))

        lig_style = self.args.input.lig_style.strip('/')

        for d_handler in loaded_datasets:
            d_handler.initialise_output_directory(outputdir=os.path.join(self.output_handler.get_dir('processed_datasets'), d_handler.d_name))
            # These will be links to the input files
            d_handler.output_handler.add_file(file_name='{!s}-input.pdb'.format(d_handler.d_tag), file_tag='input_structure')
            d_handler.output_handler.add_file(file_name='{!s}-input.mtz'.format(d_handler.d_tag), file_tag='input_data')
            # Add links to the ligand files if they've been found
            d_handler.output_handler.add_dir(dir_name='ligand_files', dir_tag='ligand', top_dir_tag='root')
            d_handler.output_handler.add_file(file_name='{!s}-ligand.pdb'.format(d_handler.d_tag), file_tag='ligand_coordinates')
            d_handler.output_handler.add_file(file_name='{!s}-ligand.cif'.format(d_handler.d_tag), file_tag='ligand_restraints')
            d_handler.output_handler.add_file(file_name='{!s}-ligand.png'.format(d_handler.d_tag), file_tag='ligand_image')
            # Aligned structure (created maps will be aligned to this structure)
            d_handler.output_handler.add_file(file_name='{!s}-aligned.pdb'.format(d_handler.d_tag), file_tag='aligned_structure')
            # Sampled map for the aligned structure
            d_handler.output_handler.add_file(file_name='{!s}-observed.ccp4'.format(d_handler.d_tag), file_tag='sampled_map')
            d_handler.output_handler.add_file(file_name='{!s}-difference.ccp4'.format(d_handler.d_tag), file_tag='difference_map')
            # Difference from the mean map for the aligned structure
            d_handler.output_handler.add_file(file_name='{!s}-mean_diff.ccp4'.format(d_handler.d_tag), file_tag='mean_diff_map')
            # Z-map for the structure
            d_handler.output_handler.add_file(file_name='{!s}-z_map_naive.ccp4'.format(d_handler.d_tag), file_tag='z_map_naive')
            d_handler.output_handler.add_file(file_name='{!s}-z_map_naive_normalised.ccp4'.format(d_handler.d_tag), file_tag='z_map_naive_normalised')
            d_handler.output_handler.add_file(file_name='{!s}-z_map_adjusted.ccp4'.format(d_handler.d_tag), file_tag='z_map_corrected')
            d_handler.output_handler.add_file(file_name='{!s}-z_map_adjusted_normalised.ccp4'.format(d_handler.d_tag), file_tag='z_map_corrected_normalised')

            # Output images
            d_handler.output_handler.add_dir(dir_name='output_images', dir_tag='images', top_dir_tag='root')
            d_handler.output_handler.add_file(file_name='{!s}-obsv_map_dist.png'.format(d_handler.d_tag), file_tag='s_map_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name='{!s}-diff_map_dist.png'.format(d_handler.d_tag), file_tag='d_map_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name='{!s}-diff_mean_map_dist.png'.format(d_handler.d_tag), file_tag='d_mean_map_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name='{!s}-z_map_dist_naive.png'.format(d_handler.d_tag), file_tag='z_map_naive_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name='{!s}-z_map_dist_naive_normalised.png'.format(d_handler.d_tag), file_tag='z_map_naive_normalised_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name='{!s}-z_map_dist_adjusted.png'.format(d_handler.d_tag), file_tag='z_map_corrected_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name='{!s}-z_map_dist_adjusted_normalised.png'.format(d_handler.d_tag), file_tag='z_map_corrected_normalised_png', dir_tag='images')

            d_handler.output_handler.add_file(file_name='{!s}-uncertainty-qqplot.png'.format(d_handler.d_tag), file_tag='unc_qqplot_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name='{!s}-mean-v-obs-sorted-qqplot.png'.format(d_handler.d_tag), file_tag='obs_qqplot_sorted_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name='{!s}-mean-v-obs-unsorted-plot.png'.format(d_handler.d_tag), file_tag='obs_qqplot_unsorted_png', dir_tag='images')

            # Z-map mask for the structure
            d_handler.output_handler.add_file(file_name='{!s}-high_z_mask.ccp4'.format(d_handler.d_tag), file_tag='high_z_mask')

            # Edstats Scores
            d_handler.output_handler.add_file(file_name='{!s}-edstats.scores'.format(d_handler.d_tag), file_tag='edstats_scores')

            # Analysis files
            d_handler.output_handler.add_file(file_name='{!s}-z_map_peaks.csv'.format(d_handler.d_tag), file_tag='z_peaks_csv')

            # Scripts
            d_handler.output_handler.add_dir(dir_name='scripts', dir_tag='scripts', top_dir_tag='root')
            d_handler.output_handler.add_file(file_name='load_maps.pml', file_tag='pymol_script', dir_tag='scripts')
            d_handler.output_handler.add_file(file_name='ccp4mg_{!s}_{!s}.py', file_tag='ccp4mg_script', dir_tag='scripts')

            # Output blobs
            d_handler.output_handler.add_dir(dir_name='blobs', dir_tag='blobs', top_dir_tag='root')
            d_handler.output_handler.add_file(file_name='blob_{!s}_{!s}.png', file_tag='ccp4mg_png', dir_tag='blobs')

            # Link the input files to the output folder
            if not os.path.exists(d_handler.output_handler.get_file('input_structure')):
                os.symlink(d_handler.get_pdb_filename(), d_handler.output_handler.get_file('input_structure'))
            if not os.path.exists(d_handler.output_handler.get_file('input_data')):
                os.symlink(d_handler.get_mtz_filename(), d_handler.output_handler.get_file('input_data'))

            # Search for ligand files and link them to the output ligands folder
            lig_glob  = os.path.join(os.path.dirname(d_handler.get_pdb_filename()), lig_style)
            lig_files = glob.glob(lig_glob)

            for lig_file in lig_files:
                # Find all files with the same basename but allowing for different extensions. Then link to output folder.
                lig_base = os.path.splitext(lig_file)[0] + '.*'
                lig_matches = glob.glob(lig_base)
                for lig in lig_matches:
                    out_path = os.path.join(d_handler.output_handler.get_dir('ligand'), os.path.basename(lig))
                    if not os.path.exists(out_path):
                        os.symlink(lig, out_path)

        self.datasets.add(loaded_datasets)
        self.log('{!s} Datasets Loaded.          '.format(len(loaded_datasets), True))

    def select_reference_dataset(self, method='resolution', max_rfree=0.4, min_resolution=5):
        """Select dataset to act as the reference - scaling, aligning etc"""

        assert method in ['resolution','rfree'], 'METHOD FOR SELECTING THE REFERENCE DATASET NOT RECOGNISED: {!s}'.format(method)

        if self.args.input.reference.pdb and self.args.input.reference.pdb:
            self.log('===================================>>>', True)
            self.log('Reference Provided by User', True)
            return self.args.input.reference.pdb, self.args.input.reference.mtz
        else:
            self.log('===================================>>>', True)
            self.log('Selecting Reference Dataset by: {!s}'.format(method), True)
            if method == 'rfree':
                # Get RFrees of datasets (set to dummy value of 999 if resolution is too high so that it is not selected)
                r_frees = [d.get_input().get_r_rfree_sigma().r_free if (d.get_resolution() < min_resolution) else 999 for d in self.datasets.all()]
                if len(resolns) == 0: raise Exception('NO DATASETS BELOW RESOLUTION CUTOFF {!s}A - CANNOT SELECT REFERENCE DATASET'.format(min_resolution))
                self._ref_dataset_index = r_frees.index(min(r_frees))
            elif method == 'resolution':
                # Get Resolutions of datasets (set to dummy value of 999 if r-free is too high so that it is not selected)
                resolns = [d.get_resolution() if (d.get_input().get_r_rfree_sigma().r_free < max_rfree) else 999 for d in self.datasets.all()]
                if len(resolns) == 0: raise Exception('NO DATASETS BELOW RFREE CUTOFF {!s} - CANNOT SELECT REFERENCE DATASET'.format(max_rfree))
                self._ref_dataset_index = resolns.index(min(resolns))

            reference = self.datasets.all()[self._ref_dataset_index]
            assert reference.d_num == self._ref_dataset_index, 'INDEX DOES NOT MATCH DNUM'
            self.log('Reference Selected: Dataset {!s}'.format(self._ref_dataset_index+1), True)
            self.log('Resolution: {!s}, RFree: {!s}'.format(reference.get_resolution(), reference.get_input().get_r_rfree_sigma().r_free), True)

            return reference.get_pdb_filename(), reference.get_mtz_filename()

#    def scale_datasets(self):
#        """Iterate through the datasets, and scale the reflections to the reference"""
#
#        self.log('===================================>>>', True)
#        for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True):
#            print '\rScaling Dataset {!s}          '.format(d_handler.d_tag),; sys.stdout.flush()
#
#            # Scale new data to the reference dataset
#            d_handler.scale_fobs_to_reference(ref_miller=self.reference_dataset().get_fobs_miller_array())
#        self.log('\rDatasets Scaled.               ', True)

    def scale_datasets(self, ampl_label, phas_label):
        """Extract amplitudes and phases for creating map"""

        def extract_structure_factors(mtz_object, ampl_label, phas_label):

            # Get the crystal symmetry from the amplitudes' crystal
            ampl_col = mtz_object.get_column(ampl_label)
            crystal_symmetry = ampl_col.mtz_crystal().crystal_symmetry()
            # Create miller set
            mill_idx = mtz_object.extract_miller_indices()
            mill_set = cctbx.miller.set(crystal_symmetry=crystal_symmetry, indices=mill_idx)
            # Extract amplitudes and phases
            ampl_com = mtz_object.extract_complex(column_label_ampl=ampl_label, column_label_phi=phas_label)
            # Convert to miller array
            mill_sfs = mill_set.array(ampl_com.data)
            mill_sfs.set_observation_type_xray_amplitude()

            # Check it's complex
            assert mill_sfs.is_complex_array(), 'STRUCTURE FACTORS SHOULD BE COMPLEX?!'

            # Make non-anomalous
            mill_sfs = mill_sfs.as_non_anomalous_array()

            return mill_sfs

        if self.args.input.reference.ampl_label:    ref_ampl_label = self.args.input.reference.ampl_label
        else:                                       ref_ampl_label = ampl_label
        if self.args.input.reference.phas_label:    ref_phas_label = self.args.input.reference.phas_label
        else:                                       ref_phas_label = phas_label

        # TODO TRY EXCEPT ON FAILING TO EXTRACT
        ref_sfs = extract_structure_factors(self.reference_dataset().reflection_data(), ampl_label=ref_ampl_label, phas_label=ref_phas_label)
        # Extract the amplitudes for scaling
        ref_amps = ref_sfs.amplitudes()
        # Save the structure factors to the reference dataset object to create the map from later
        self.reference_dataset().unscaled_sfs = ref_sfs
        self.reference_dataset().scaled_sfs = ref_sfs

        t1 = time.time()

        self.log('===================================>>>', True)
        for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True):

            if d_handler.scaled_sfs != None:
                print '\rAlready Scaled: Dataset {!s}          '.format(d_handler.d_tag),; sys.stdout.flush()
                continue
            else:
                print '\rScaling: Dataset {!s}                 '.format(d_handler.d_tag),; sys.stdout.flush()

            # Get reflection data object
            refl_data = d_handler.reflection_data()

            # Extract miller array of structure factors
            unscaled_sfs = extract_structure_factors(refl_data, ampl_label=ampl_label, phas_label=phas_label)

            # Extract just the amplitudes for scaling
            d_amps = unscaled_sfs.amplitudes()
            d_phas = unscaled_sfs.phases()

            assert d_amps.is_real_array(), 'AMPLITUDES SHOULD BE REAL?!'
            assert d_phas.is_real_array(), 'PHASE ANGLES SHOULD BE REAL?!'

            # Scale new data to the reference dataset
            scaled_amps = apply_simple_scaling(miller=d_amps, ref_miller=ref_amps)

            # Recombine the scaled amplitudes and the phases
            scaled_sf_real = scaled_amps.data() * flex.cos(d_phas.data()*math.pi/180.0)
            scaled_sf_imag = scaled_amps.data() * flex.sin(d_phas.data()*math.pi/180.0)
            scaled_sf_com = flex.complex_double(reals=scaled_sf_real, imags=scaled_sf_imag)
            # Create a copy of the old array with the new scaled structure factors
            scaled_sfs = unscaled_sfs.array(data=scaled_sf_com)

            # Set the scaled structure factors
            d_handler.unscaled_sfs = unscaled_sfs
            d_handler.scaled_sfs = scaled_sfs
            # At the moment save the 'truncated' sfs as the whole list
            d_handler.tr_scaled_sfs = scaled_sfs

        self.log('\rDatasets Scaled.               ', True)

        t2 = time.time()
        print('> DIFFRACTION DATA SCALING > Time Taken: {!s} seconds'.format(int(t2-t1)))

    def align_datasets(self, method, sites='calpha'):
        """Align each structure the reference structure"""

        assert sites in ['calpha','backbone']
        assert method in ['local','global'], 'METHOD NOT DEFINED: {!s}'.format(method)

        # If local alignment has been chosen, also do a global alignment
        if method == 'local': method = 'both'

        t1 = time.time()

        self.log('===================================>>>', True)
        for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True):
            print '\rAligning Dataset {!s}          '.format(d_handler.d_tag),; sys.stdout.flush()

            # Align to reference structure to get mapping transform
            d_handler.align_to_reference(ref_handler=self.reference_dataset(), sites=sites, method=method)

            # Create a copy to transform
            aligned_struc = d_handler.get_hierarchy().deep_copy()

            # Always do a global alignment
            aligned_struc.atoms().set_xyz(d_handler.transform_to_reference(points=d_handler.get_hierarchy().atoms().extract_xyz(), method='global'))
            # Also write out into the aligned_structures folder
            aligned_struc.write_pdb_file(file_name=os.path.join(self.output_handler.get_dir('aligned_structures'), '{!s}-global-aligned.pdb'.format(d_handler.d_tag)))

            if method == 'both':
                # Write out the aligned structure from the local alignments
                aligned_struc.atoms().set_xyz(d_handler.transform_to_reference(points=d_handler.get_hierarchy().atoms().extract_xyz(), method='local'))
                # Also write out into the aligned_structures folder
                aligned_struc.write_pdb_file(file_name=os.path.join(self.output_handler.get_dir('aligned_structures'), '{!s}-local-aligned.pdb'.format(d_handler.d_tag)))

            # Write out next to the original files (if global, global, if local, local)
            aligned_struc.write_pdb_file(file_name=d_handler.output_handler.get_file('aligned_structure'))

            # Align to the reference to find the alignment rmsd
            mapped_ref_sites = d_handler.transform_from_reference(points=self.reference_dataset().get_hierarchy().atoms().extract_xyz(), method='global')
            alignment_rmsd = d_handler.get_hierarchy().atoms().extract_xyz().rms_difference(mapped_ref_sites)
            # Check to see if this is the reference dataset
            if alignment_rmsd == 0.0:
                # This must be the reference! Set the dataset number if it's not already set
                if self._ref_dataset_index == None:
                    self._ref_dataset_index = d_handler.d_num
                    print 'REFERENCE FOUND! {!s}                          '.format(d_handler.d_tag)
                # Raise error if the reference has already been set and it's not this dataset
                elif self._ref_dataset_index != d_handler.d_num:
                    raise Exception('ALIGNED OBJECT EQUAL TO UNALIGNED OBJECT - THIS IS MOST UNLIKELY')
        self.log('\rDatasets Aligned.               ', True)

        t2 = time.time()
        print('> STRUCTURE ALIGNMENT > Time Taken: {!s} seconds'.format(int(t2-t1)))

    def analyse_dataset_variability_1(self):
        """Go through all of the datasets and collect lots of different characteristics of the datasets for identifying odd datasets"""

        self.log('===================================>>>', True)
        self.log('Collecting Dataset/Crystal Variation Data - 1', True)
        self.log('===================================>>>')

        self.log('Extracting Resolutions')
        self.datasets_summary.add_data( data_name='high_res_limit',
                                        data_values=[d.unscaled_sfs.d_min() if d.unscaled_sfs else numpy.nan for d in self.datasets.all()]
                                        )
        self.datasets_summary.add_data( data_name='low_res_limit',
                                        data_values=[d.unscaled_sfs.d_max_min()[0] if d.unscaled_sfs else numpy.nan for d in self.datasets.all()]
                                        )
        self.log('Extracting Unit Cell Sizes')
        self.datasets_summary.add_data( data_name='cell_params',
                                        data_values=[d.unscaled_sfs.unit_cell().parameters() if d.unscaled_sfs else [numpy.nan]*6 for d in self.datasets.all()]
                                        )
        self.log('Extracting Unit Cell Volume')
        self.datasets_summary.add_data( data_name='cell_volume',
                                        data_values=[d.unscaled_sfs.unit_cell().volume() if d.unscaled_sfs else numpy.nan for d in self.datasets.all()]
                                        )
        self.log('Extracting R-work, R-free')
        self.datasets_summary.add_data( data_name='rfree',
                                        data_values=[d.get_input().get_r_rfree_sigma().r_free for d in self.datasets.all()]
                                        )
        self.datasets_summary.add_data( data_name='rwork',
                                        data_values=[d.get_input().get_r_rfree_sigma().r_work for d in self.datasets.all()]
                                        )
        self.log('Calculating Variation in Correlations between Diffraction Data')
        self.datasets_summary.add_data( data_name='correlation_to_unscaled_data',
                                        data_values=[d.scaled_sfs.amplitudes().correlation(d.unscaled_sfs.amplitudes()).coefficient() if d.scaled_sfs else numpy.nan for d in self.datasets.all()]
                                        )

    def filter_datasets_1(self):
        """Filter out the datasets which contain different protein models (i.e. protein length, sequence, etc)"""

        self.log('===================================>>>', True)
        self.log('Filtering Datasets (Non-identical structures)). Potential Classes:', True)
        for failure_class in self._structure_mask_names:
            self.log('\t{!s}'.format(failure_class), True)
        self.log('===================================>>>', True)

        ref_counts, ref_dict = self.reference_dataset().get_structure_summary()

        # Check that the same protein structure is present in each dataset - THIS MASK SHOULD INCLUDE ALL DATASETS AT FIRST
        for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True):

            my_counts, my_dict = d_handler.get_structure_summary()

            print '\rFiltering Dataset {!s}          '.format(d_handler.d_tag),; sys.stdout.flush()
            # Check the number of chains
            if my_counts.n_chains != ref_counts.n_chains:
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.d_tag))
                self.log('Different Number of Chains')
                self.log('Reference: {!s}, {!s}: {!s}'.format(ref_counts.n_chains, d_handler.d_tag, my_counts.n_chains))
                self.log('===================================>>>')
                self.datasets.all_masks().set_mask_value(mask_name='bad structure - chain counts', entry_id=d_handler.d_tag, value=True)
            # Check the ids of the chains
            elif my_counts.chain_ids != ref_counts.chain_ids:
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.d_tag))
                self.log('Different Chain IDs')
                self.log('Reference: {!s}, {!s}: {!s}'.format(ref_counts.chain_ids, d_handler.d_tag, my_counts.chain_ids))
                print('===================================>>>')
                self.datasets.all_masks().set_mask_value(mask_name='bad structure - chain ids', entry_id=d_handler.d_tag, value=True)
            # Check the sequences of the chains
            elif my_dict['chain_sequences'] != ref_dict['chain_sequences']:
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.d_tag))
                self.log('Different Sequences')
                for chain_id in ref_dict['chain_sequences'].keys():
                    self.log('Chain {!s}: Reference - {!s}'.format(chain_id, ref_dict['chain_sequences'][chain_id]))
                    self.log('Chain {!s}: {:>9s} - {!s}'.format(chain_id, my_dict['chain_sequences'][chain_id]))
                print('===================================>>>')
                self.datasets.all_masks().set_mask_value(mask_name='bad structure - chain sequences', entry_id=d_handler.d_tag, value=True)
            # Check the number of residues - TODO not sure it can ever get here... remove?
            elif my_counts.n_residues != ref_counts.n_residues:
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.d_tag))
                self.log('Different Number of Residues')
                self.log('Reference: {!s}, {!s}: {!s}'.format(ref_counts.n_residues, d_handler.d_tag, my_counts.n_residues))
                print('===================================>>>')
                self.datasets.all_masks().set_mask_value(mask_name='bad structure - residue counts', entry_id=d_handler.d_tag, value=True)
            # Check the number of atoms
            elif my_counts.n_atoms != ref_counts.n_atoms:
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.d_tag))
                self.log('Different Number of Atoms')
                self.log('Reference: {!s}, {!s}: {!s}'.format(ref_counts.n_atoms, d_handler.d_tag, my_counts.n_atoms))
                print('===================================>>>')
                self.datasets.all_masks().set_mask_value(mask_name='bad structure - atom counts', entry_id=d_handler.d_tag, value=True)
            else:
                pass

        # Combine structure_masks
        structure_reject_mask = self.datasets.all_masks().combine_masks(self._structure_mask_names)
        self.datasets.all_masks().add_mask(mask_name='rejected - structure', mask=structure_reject_mask)

        # Update the combined masks
        combined_reject_mask = self.datasets.all_masks().combine_masks(self._reject_mask_names)
        self.datasets.all_masks().add_mask(mask_name='rejected - total', mask=combined_reject_mask)

        self.log('\rDatasets Filtered.               ', True)
        self.log('Rejected Datasets (Structure): {!s}'.format(sum(self.datasets.all_masks().get_mask(mask_name='rejected - structure'))), True)
        self.log('Rejected Datasets (Total):     {!s}'.format(sum(self.datasets.all_masks().get_mask(mask_name='rejected - total'))), True)

        # Link all rejected datasets into the rejected directory
        for d_handler in self.datasets.mask(mask_name='rejected - total'):
            reject_dir = os.path.join(self.output_handler.get_dir('rejected_datasets'), d_handler.d_name)
            if not os.path.exists(reject_dir):
                rel_symlink(orig=d_handler.output_handler.get_dir('root'), link=reject_dir)

    def filter_datasets_2(self):
        """Filter out the non-isomorphous datasets"""

        self.log('===================================>>>', True)
        self.log('Filtering Datasets (Non-isomorphous datasets). Potential Classes:', True)
        for failure_class in self._crystal_mask_names:
            self.log('\t{!s}'.format(failure_class), True)
        self.log('===================================>>>', True)

        # Check that each dataset is similar enough to be compared
        for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True):

            print '\rFiltering Dataset {!s}          '.format(d_handler.d_tag),; sys.stdout.flush()
            # Check that it correlates well with itself before and after scaling
            if d_handler.scaled_sfs.amplitudes().correlation(d_handler.unscaled_sfs.amplitudes()).coefficient() < FILTER_SCALING_CORRELATION_CUTOFF:
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.d_tag))
                self.log('Low correlation between scaled and unscaled data')
                self.log('Scaled-Unscaled Correlation: {!s}'.format(d_handler.scaled_sfs.amplitudes().correlation(d_handler.unscaled_sfs.amplitudes()).coefficient()))
                self.log('===================================>>>')
                self.datasets.all_masks().set_mask_value(mask_name='bad crystal - data correlation', entry_id=d_handler.d_tag, value=True)
            # Check the deviation from the average sites
            elif d_handler.get_calpha_sites().rms_difference(d_handler.transform_from_reference(points=self.reference_dataset().get_calpha_sites(), method='global')) > FILTER_GLOBAL_RMSD_CUTOFF:
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.d_tag))
                self.log('C-alpha RMSD is too large')
                self.log('Aligned (Calpha) RMSD: {!s}'.format(d_handler.get_calpha_sites().rms_difference(d_handler.transform_from_reference(points=self.reference_dataset().get_calpha_sites(), method='global'))))
                self.log('===================================>>>')
                self.datasets.all_masks().set_mask_value(mask_name='bad crystal - isomorphous structure', entry_id=d_handler.d_tag, value=True)
            else:
                pass

        # Combine crystal masks
        crystal_reject_mask = self.datasets.all_masks().combine_masks(self._crystal_mask_names)
        self.datasets.all_masks().add_mask(mask_name='rejected - crystal', mask=crystal_reject_mask)

        # Combine all of the masks
        combined_reject_mask = self.datasets.all_masks().combine_masks(self._reject_mask_names)
        self.datasets.all_masks().add_mask(mask_name='rejected - total', mask=combined_reject_mask)

        self.log('\rDatasets Filtered.               ', True)
        self.log('Rejected Datasets (Structure): {!s}'.format(sum(self.datasets.all_masks().get_mask(mask_name='rejected - structure'))), True)
        self.log('Rejected Datasets (Crystal):   {!s}'.format(sum(self.datasets.all_masks().get_mask(mask_name='rejected - crystal'))), True)
        self.log('Rejected Datasets (Total):     {!s}'.format(sum(self.datasets.all_masks().get_mask(mask_name='rejected - total'))), True)

        # Link all rejected datasets into the rejected directory
        for d_handler in self.datasets.mask(mask_name='rejected - total'):
            reject_dir = os.path.join(self.output_handler.get_dir('rejected_datasets'), d_handler.d_name)
            if not os.path.exists(reject_dir):
                rel_symlink(orig=d_handler.output_handler.get_dir('root'), link=reject_dir)

    def filter_datasets_3(self, resolution):
        """Select datasets based on resolution"""

        self.log('===================================>>>', True)
        self.log('Selecting datasets based on resolution', True)
        self.log('===================================>>>', True)

        # Select from the datasets that haven't been rejected
        for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True):

            print '\rFiltering Dataset {!s}          '.format(d_handler.d_tag),; sys.stdout.flush()

            # Check the resolution of the dataset
            if d_handler.get_resolution() > resolution:
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.d_tag))
                self.log('Does not meet high-resolution cutoff')
                self.log('Cut Resolution: {!s}'.format(resolution))
                self.log('Map Resolution: {!s}'.format(d_handler.get_resolution()))
                self.log('===================================>>>')
                self.datasets.all_masks().set_mask_value(mask_name='selected for analysis', entry_id=d_handler.d_tag, value=False)
            else:
                self.datasets.all_masks().set_mask_value(mask_name='selected for analysis', entry_id=d_handler.d_tag, value=True)

        self.log('\rDatasets Filtered.               ', True)
        self.log('{!s} datasets selected for analysis at {!s}A'.format(sum(self.datasets.all_masks().get_mask(mask_name='selected for analysis')), resolution), True)

        # TODO LINK THESE INTO THE APPROPRIATE RESOLUTION FOLDER - TO BE CREATED

    def calculate_mean_structure_and_protein_masks(self, deviation_cutoff):
        """Calculate the average of all of the structures, and create masks for each protein where residues deviate from the mean by more than `deviation_cutoff`"""

        self.log('===================================>>>', True)
        self.log('Calculating Mean Structure', True)
        self.log('===================================>>>')

        # TODO Make this reject points until consensus

        # Pull all c-alpha sites for each structure
        all_sites = numpy.array([d.transform_to_reference(points=d.get_calpha_sites(), method='global') for d in self.datasets.mask(mask_name='rejected - total', invert=True)])
        # Calculate the mean x,y,z for each c-alpha
        mean_sites = numpy.mean(all_sites, axis=0)
        # Differences from the mean
        diff_sites = all_sites - mean_sites
        # Euclidean norms of the distances moved
        diff_norms = numpy.apply_along_axis(numpy.linalg.norm, axis=2, arr=diff_sites)

        # TODO MOVE THIS TO THE STRUCTURE VARIATION FUNCTION TODO

        # TODO CREATE A HIERARCHY FOR THE MEAN STRUCTURE (AND WITH MEAN NORMALISED B-FACTORS?)

        # Create a list of masks for large-moving c-alphas
        residue_deviation_masks = []
        # Iterate by dataset, masking if the deviation of the calpha in the dataset is more than `deviation_cutoff`
        for calpha_shifts in diff_norms:
            residue_deviation_masks.append([1 if shift > deviation_cutoff else 0 for shift in calpha_shifts])

        # Save the masks
        self._average_calpha_sites = flex.vec3_double(mean_sites)
        self._residue_deviation_masks = residue_deviation_masks

    def analyse_dataset_variability_2(self):
        """Go through all of the datasets and collect lots of different characteristics of the datasets for identifying odd datasets"""

        self.log('===================================>>>', True)
        self.log('Collecting Dataset/Crystal Variation Data - 2', True)
        self.log('===================================>>>')

        # Now calculate the variation in the structure, from the reference
        self.log('Calculating Variation in RMSD (Calphas) to Reference Structure')
        rmsds = numpy.array([numpy.nan]*self.datasets.size())
        for d in self.datasets.mask(mask_name='rejected - total', invert=True):
            rmsds.put(d.d_num, d.get_calpha_sites().rms_difference(d.transform_from_reference(points=self.get_calpha_average_sites(), method='global')))
        self.datasets_summary.add_data(data_name='rmsd_to_mean', data_values=rmsds.tolist())

    def analyse_structure_variability_1(self):
        """Go through all of the datasets and collect lots of different structural characteristics of the datasets for identifying odd datasets"""

        self.log('===================================>>>', True)
        self.log('Populating Ensemble Summary', True)

        # Extract the hierarchy for each of the datasets
        hierarchies = [d.get_hierarchy() for d in self.datasets.mask(mask_name='rejected - total', invert=True)]
        hierarchy_ids = [d.d_tag for d in self.datasets.mask(mask_name='rejected - total', invert=True)]

        self.log(self.ensemble_summary.add_structures(new_hierarchies=hierarchies, hierarchy_ids=hierarchy_ids))

    def truncate_scaled_data(self, truncation_stuff=None):
        """Truncate data at the same indices across all the datasets"""

        # Create maps for all of the datasets (including the reference dataset)
        for d_handler in self.datasets.mask(mask_name='selected for analysis'):
            continue
            d_handler.tr_scaled_sfs = None

    def load_and_morph_maps(self, overwrite=False, delete_native_map=True):
        """Create map from miller arrays. Transform map into the reference frame by sampling at the given points."""

        self.log('===================================>>>', True)
        self.log('Loading Electron Density Maps', True)
        self.log('===================================>>>')

        # Extract the points for the morphed maps (in the reference frame)
        masked_gps = self.reference_grid().global_mask().outer_mask()
        masked_idxs = self.reference_grid().global_mask().outer_mask_indices()
        masked_cart_ref = flex.vec3_double(masked_gps)*self.reference_grid().grid_spacing()
        # Mapping of grid points to rotation matrix keys (residue CA labels)
        masked_cart_mappings = self.reference_grid().partition().query_by_grid_indices(masked_idxs)

        # Create maps for all of the datasets (including the reference dataset)
        for d_handler in self.datasets.mask(mask_name='selected for analysis'):

            ################################
            #     CALCULATE NATIVE MAP     #
            ################################

            print '\rLoading Maps for Dataset {!s} ({!s}/{!s})          '.format(d_handler.d_tag, d_handler.d_num+1, self.datasets.size()),; sys.stdout.flush()

            if d_handler.morphed_map:
                continue
            elif d_handler.native_map and (not overwrite):
                pass
            else:
                # Take the scaled diffraction data for each dataset and create fft
                fft_map = d_handler.unscaled_sfs.fft_map( resolution_factor=self.params.maps.resolution_factor,
                #TODO MAKE D_MIN CUT AT MILLER INDEX INSTEAD OF ABSOLUTE
                                                        d_min=self.get_cut_resolution(),
                                                        symmetry_flags=maptbx.use_space_group_symmetry
                                                    )

                # Scale the map
                if self.params.maps.scaling == 'none':
                    pass
                elif self.params.maps.scaling == 'sigma':
                    fft_map.apply_sigma_scaling()
                elif self.params.maps.scaling == 'volume':
                    fft_map.apply_volume_scaling()

                # Save the map for future use (only part that's pickleable)
                d_handler.native_map  = fft_map.real_map()
                d_handler.unit_cell   = fft_map.unit_cell()
                d_handler.space_group = fft_map.space_group()

                # TODO PICKLE THE NATIVE MAP

            #################################
            # MORPH MAPS TO REFERENCE FRAME #
            #################################

            assert d_handler.native_map is not None, 'NO MORPHED MAP'
            if not overwrite:   assert d_handler.morphed_map is None, 'MORPHED MAP ALREADY EXISTS'

            # Create map handler and map to the reference frame
            native_map_handler = map_handler(map_data=d_handler.native_map, unit_cell=d_handler.unit_cell)

            # Transform Coordinates
            masked_cart_d = d_handler.transform_from_reference(points=masked_cart_ref, method=self.params.alignment.method, point_mappings=masked_cart_mappings)
            # Sample the map at these points
            masked_vals_d = native_map_handler.get_cart_values(masked_cart_d)

            # Normalise map
            map_mean = numpy.mean(masked_vals_d)
            map_stdv = numpy.std(masked_vals_d)
            masked_vals_d = (masked_vals_d - map_mean)/map_stdv
            # Record map means
            self.get_map_observations().set_data_value(data_name='masked_map_mean', entry_id=d_handler.d_tag, value=map_mean)
            self.get_map_observations().set_data_value(data_name='masked_map_std', entry_id=d_handler.d_tag, value=map_stdv)

            # Calculate the sparse vector of the masked map values
            morphed_map_sparse = sparse.vector(self.reference_grid().grid_size_1d(), dict(zip(masked_idxs, masked_vals_d)))
            morphed_map        = morphed_map_sparse.as_dense_vector()

            # Reshape into right shape of the grid
            grid_shape = self.reference_grid().grid_size()
            morphed_map.reshape(flex.grid(grid_shape))

#            # Write out the morphed map
#            self.write_array_to_map(output_file=d_handler.output_handler.get_file('sampled_map'),
#                                    map_data=morphed_map)

            # Save to the dataset handler
            d_handler.morphed_map = morphed_map

            # Delete the native map to save space
            if delete_native_map:
                d_handler.native_map = None
                # Collect
                gc.collect()

        self.log('\rMaps Values Loaded: {!s} Datasets'.format(self.datasets.size(mask_name='selected for analysis')), True)

        self.update_pandda_size(tag='After Loading Maps')

    def calculate_mean_map(self):
        """Calculate the mean map from all of the different observations"""

        # Create statistics objects for each grid point
        self.log('===================================>>>', True)
        self.log('Calculating Mean Map', True)

        # All dataset handlers for analysing
        all_d_handlers = self.datasets.mask(mask_name='selected for analysis')

        # Get masked points in reference dataset
        masked_idxs = self.reference_grid().global_mask().outer_mask_indices()

        # Chunk the points into groups - Compromise between cpu time and memory usage - ~200 dataset -> chunksize of 5000
        chunk_size = 1000*int(2000/len(all_d_handlers))
        chunked_indices = [masked_idxs[i:i + chunk_size] for i in range(0, len(masked_idxs), chunk_size)]
        chunk_num = len(chunked_indices)

        self.log('Iterating through {!s} points in {!s} chunks'.format(len(masked_idxs), chunk_num), True)
        t1 = time.time()

        point_means = []

        # Calculate the mean map across the datasets
        for chunk_idx, chunk in enumerate(chunked_indices):
            status_bar(n=chunk_idx-1, n_max=chunk_num)
            p_map_vals = numpy.array([dh.morphed_map.select(chunk) for dh in all_d_handlers])

            if chunk_idx+1 < chunk_num:
                assert len(p_map_vals) == len(all_d_handlers)
                assert len(p_map_vals.T) == chunk_size

            p_map_means = numpy.mean(p_map_vals, axis=0)
            point_means.extend(p_map_means.tolist())

        status_bar(n=chunk_idx, n_max=chunk_num)
        assert len(point_means) == len(masked_idxs)

        t2 = time.time()
        print('> MEAN MAP CALCULATION > Time Taken: {!s} seconds'.format(int(t2-t1)))

        # Calculate Mean Maps
        mean_map_vals = numpy.zeros(self.reference_grid().grid_size_1d())
        mean_map_vals.put(masked_idxs, point_means)
        mean_map_vals = flex.double(mean_map_vals.tolist())
        mean_map_file = self.output_handler.get_file('mean_map')
        self.write_array_to_map(output_file=mean_map_file, map_data=mean_map_vals)

        self.stat_maps.mean_map = mean_map_vals

#    def old_calculate_mean_map(self):
#        """Calculate the mean map from all of the different observations"""
#
#        # Create statistics objects for each grid point
#        self.log('===================================>>>', True)
#        self.log('Calculating Mean Map', True)
#
#        # Get masked points in reference dataset
#        masked_gps = self.reference_grid().global_mask().outer_mask()
#        masked_idxs = self.reference_grid().global_mask().outer_mask_indices()
#        masked_cart_ref = flex.vec3_double(masked_gps)*self.reference_grid().grid_spacing()
#        # Mapping of grid points to rotation matrix keys (residue CA labels)
#        masked_cart_mappings = self.reference_grid().partition().query_by_grid_indices(masked_idxs)
#
#        assert len(masked_idxs) == len(masked_cart_ref)
#
#        # All dataset handlers for analysing
#        all_d_handlers = self.datasets.mask(mask_name='rejected - total', invert=True)
#
#        # Chunk the points into groups - Compromise between cpu time and memory usage - ~200 dataset -> chunksize of 5000
#        chunk_size = 1000*int(2000/self.datasets.size(mask_name='rejected - total', invert=True))
#        chunked_cart = [masked_cart_ref[i:i + chunk_size] for i in range(0, len(masked_cart_ref), chunk_size)]
#        chunked_mapping = [masked_cart_mappings[i:i + chunk_size] for i in range(0, len(masked_cart_mappings), chunk_size)]
#        chunk_num = len(chunked_cart)
#
#        self.log('Iterating through {!s} points in {!s} chunks'.format(len(masked_idxs), chunk_num), True)
#        t1 = time.time()
#
#        point_means = []
#
#        # Calculate the mean map across the datasets
#        for chunk_idx, ref_coord_chunk in enumerate(chunked_cart):
#            status_bar(n=chunk_idx, n_max=chunk_num)
#            p_map_vals = numpy.array([dh.get_cart_values(points=dh.transform_from_reference(points=flex.vec3_double(ref_coord_chunk), method=self.get_alignment_method(), point_mappings=chunked_mapping[chunk_idx]), map_type=self.get_obs_map_type()) for dh in all_d_handlers])
#
#            if chunk_idx+1 < chunk_num:
#                assert len(p_map_vals) == len(all_d_handlers)
#                assert len(p_map_vals.T) == chunk_size
#
#            point_means.extend([numpy.mean(map_vals.tolist()) for map_vals in p_map_vals.T])
#
#        t2 = time.time()
#        print('> MEAN MAP CALCULATION > Time Taken: {!s} seconds'.format(int(t2-t1)))
#
#        # Calculate Mean Maps
#        mean_map_vals = numpy.zeros(self.reference_grid().grid_size_1d())
#        mean_map_vals.put(masked_idxs, point_means)
#        mean_map_vals = flex.double(mean_map_vals.tolist())
#        mean_map_file = self.output_handler.get_file('mean_map')
#        self.write_array_to_map(output_file=mean_map_file, map_data=mean_map_vals)
#
#        self.stat_maps.mean_map = mean_map_vals

    def calculate_map_uncertainties(self):
        """Calculate the uncertainty in each of the different maps"""

        try:
            import matplotlib
            matplotlib.interactive(0)
            from matplotlib import pyplot
            output_graphs = True
        except:
            output_graphs = False

        # Section of q-q plot to calculate slope over (-q_cut -> +q_cut)
        q_cut = 1.5

        # Mask grid points and select those near to the protein
        masked_idxs = self.reference_grid().global_mask().inner_mask_indices()

        # Extract the mean map values
        masked_mean_vals = self.stat_maps.mean_map.select(masked_idxs)
        sorted_mean_vals = sorted(masked_mean_vals)

        # Extract the theoretical quantiles that we would expect if these values were from a normal distribution
        expected_diff_vals = normal_distribution().quantiles(len(masked_idxs))
        # Select the points in the middle of the distribution
        middle_indices = (expected_diff_vals < q_cut).iselection().intersection((expected_diff_vals > -1*q_cut).iselection())
        middle_expected_diff_vals = expected_diff_vals.select(middle_indices)

        self.log('===================================>>>')

        for d_handler in self.datasets.mask(mask_name='selected for analysis'):

            if d_handler.get_map_uncertainty() != None:
                self.get_map_observations().set_data_value(data_name='masked_map_uncertainty', entry_id=d_handler.d_tag, value=d_handler.get_map_uncertainty())
                continue

            print '\rCalculating Map Uncertainty for Dataset {!s} ({!s}/{!s})          '.format(d_handler.d_tag, d_handler.d_num+1, self.datasets.size()),; sys.stdout.flush()

            # Extract the map values at the transformed grid points - FOR THE DIFFERENCE MAP
            masked_map_vals = d_handler.morphed_map.select(masked_idxs)

            # Subtract the mean map from the observed map
            diff_mean_map_vals = masked_map_vals - masked_mean_vals

            # Sort the map values
            sorted_observed_vals = sorted(masked_map_vals)
            sorted_observed_diff_vals = sorted(diff_mean_map_vals)

            # Select the middle indices
            middle_observed_diff_vals = flex.double(sorted_observed_diff_vals).select(middle_indices)
            # Calculate the slope of the centre of the graph
            fit_coeffs = numpy.polyfit(x=middle_expected_diff_vals, y=middle_observed_diff_vals, deg=1)
            map_off = fit_coeffs[1]
            map_unc = fit_coeffs[0]

            # Save the uncertainty
            d_handler.set_map_uncertainty(sigma=map_unc)
#            print '> UNCERTAINTY: {!s}'.format(fit_coeffs.round(3).tolist()[0])
            self.get_map_observations().set_data_value(data_name='masked_map_uncertainty', entry_id=d_handler.d_tag, value=d_handler.get_map_uncertainty())

            if output_graphs:

                fig = pyplot.figure()
                pyplot.title('UNSORTED MEAN v OBS SCATTER PLOT')
                pyplot.plot([-3, 10], [-3, 10], 'b--')
                pyplot.plot(masked_mean_vals, masked_map_vals, 'go')
                pyplot.xlabel('MEAN MAP VALUES')
                pyplot.ylabel('OBSV MAP VALUES')
                # Apply tight layout to prevent overlaps
                pyplot.tight_layout()
                # Save
                pyplot.savefig(d_handler.output_handler.get_file('obs_qqplot_unsorted_png'))
                pyplot.close(fig)

                fig = pyplot.figure()
                pyplot.title('SORTED MEAN v OBS Q-Q PLOT')
                pyplot.plot([-3, 10], [-3, 10], 'b--')
                pyplot.plot(sorted_mean_vals, sorted_observed_vals, 'go-')
                pyplot.xlabel('SORTED MEAN MAP VALUES')
                pyplot.ylabel('SORTED OBSV MAP VALUES')
                # Apply tight layout to prevent overlaps
                pyplot.tight_layout()
                # Save
                pyplot.savefig(d_handler.output_handler.get_file('obs_qqplot_sorted_png'))
                pyplot.close(fig)

                fig = pyplot.figure()
                pyplot.title('DIFF-MEAN-MAP Q-Q PLOT')
                pyplot.plot([map_off-5*map_unc, map_off+5*map_unc], [-5, 5], 'b--')
                pyplot.plot([-1, 1], [q_cut, q_cut], 'k-.')
                pyplot.plot([-1, 1], [-q_cut, -q_cut], 'k-.')
                pyplot.plot(sorted_observed_diff_vals, expected_diff_vals, 'go-')
                pyplot.xlabel('DIFF-MEAN-MAP QUANTILES')
                pyplot.ylabel('THEORETICAL QUANTILES')
                # Apply tight layout to prevent overlaps
                pyplot.tight_layout()
                # Save
                pyplot.savefig(d_handler.output_handler.get_file('unc_qqplot_png'))
                pyplot.close(fig)

        self.log('\rMap Uncertainties Calculated.                                         ', True)

        try:
            from ascii_graph import Pyasciigraph
            g=Pyasciigraph()
            graph_data = [(d.d_tag, round(d.get_map_uncertainty(),3)) for d in self.datasets.mask(mask_name='selected for analysis')]
            for l in g.graph(label='Sorted Map Uncertainties (Ascending Order)', data=graph_data, sort=1):
                print(l)
        except ImportError:
            print('IMPORT ERROR (ascii_graph) - CANNOT GENERATE UNCERTAINTY GRAPH')
            pass

#    def old_calculate_map_uncertainties(self):
#        """Calculate the uncertainty in each of the different maps"""
#
#        # Section of q-q plot to calculate slope over (-q_cut -> +q_cut)
#        q_cut = 1.5
#
#        # Mask grid points and select those near to the protein
#        masked_gps = self.reference_grid().global_mask().inner_mask()
#        masked_idxs = self.reference_grid().global_mask().inner_mask_indices()
#        masked_cart_ref = flex.vec3_double(masked_gps)*self.reference_grid().grid_spacing()
#        # Mapping of grid points to rotation matrix keys (residue CA labels)
#        masked_cart_mappings = self.reference_grid().partition().query_by_grid_indices(masked_idxs)
#
#        # Extract the mean map values
#        masked_mean_vals = self.stat_maps.mean_map.select(masked_idxs)
#        sorted_mean_vals = sorted(masked_mean_vals)
#
#        # Extract the theoretical quantiles that we would expect if these values were from a normal distribution
#        expected_diff_vals = normal_distribution().quantiles(len(masked_idxs))
#        # Select the points in the middle of the distribution
#        middle_indices = (expected_diff_vals < q_cut).iselection().intersection((expected_diff_vals > -1*q_cut).iselection())
#        middle_expected_diff_vals = expected_diff_vals.select(middle_indices)
#
#        try:
#            import matplotlib
#            # Setup so that we can write without a display connected
#            matplotlib.interactive(0)
#            default_backend, validate_function = matplotlib.defaultParams['backend']
#            from matplotlib import pyplot
#            output_graphs = True
#        except:
#            output_graphs = False
#
#        self.log('===================================>>>')
#
#        for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True):
#
#            print '\rCalculating Map Uncertainty for Dataset {!s} ({!s}/{!s})          '.format(d_handler.d_tag, d_handler.d_num+1, self.datasets.size()),; sys.stdout.flush()
#
#            # Extract the map values at the transformed grid points - FOR THE DIFFERENCE MAP
#            masked_cart_new = d_handler.transform_from_reference(points=masked_cart_ref, method=self.get_alignment_method(), point_mappings=masked_cart_mappings)
#            masked_map_vals = d_handler.get_cart_values(points=masked_cart_new, map_type=self.get_obs_map_type())
#
#            # Subtract the mean map from the observed map
#            diff_mean_map_vals = masked_map_vals - masked_mean_vals
#
#            # Sort the map values
#            sorted_observed_vals = sorted(masked_map_vals)
#            sorted_observed_diff_vals = sorted(diff_mean_map_vals)
#
#            # Select the middle indices
#            middle_observed_diff_vals = flex.double(sorted_observed_diff_vals).select(middle_indices)
#            # Calculate the slope of the centre of the graph
#            fit_coeffs = numpy.polyfit(x=middle_expected_diff_vals, y=middle_observed_diff_vals, deg=1)
#            map_off = fit_coeffs[1]
#            map_unc = fit_coeffs[0]
#
#            # Save the uncertainty
#            d_handler.set_map_uncertainty(sigma=map_unc)
##            print '> UNCERTAINTY: {!s}'.format(fit_coeffs.round(3).tolist()[0])
#            self.get_map_observations().set_data_value(data_name='masked_map_uncertainty', entry_id=d_handler.d_tag, value=map_unc)
#
#            if output_graphs:
#
#                fig = pyplot.figure()
#                pyplot.title('UNSORTED MEAN v OBS SCATTER PLOT')
#                pyplot.plot([-3, 10], [-3, 10], 'b--')
#                pyplot.plot(masked_mean_vals, masked_map_vals, 'go')
#                pyplot.xlabel('MEAN MAP VALUES')
#                pyplot.ylabel('OBSV MAP VALUES')
#                # Apply tight layout to prevent overlaps
#                pyplot.tight_layout()
#                # Save
#                pyplot.savefig(d_handler.output_handler.get_file('obs_qqplot_unsorted_png'))
#                pyplot.close(fig)
#
#                fig = pyplot.figure()
#                pyplot.title('SORTED MEAN v OBS Q-Q PLOT')
#                pyplot.plot([-3, 10], [-3, 10], 'b--')
#                pyplot.plot(sorted_mean_vals, sorted_observed_vals, 'go-')
#                pyplot.xlabel('SORTED MEAN MAP VALUES')
#                pyplot.ylabel('SORTED OBSV MAP VALUES')
#                # Apply tight layout to prevent overlaps
#                pyplot.tight_layout()
#                # Save
#                pyplot.savefig(d_handler.output_handler.get_file('obs_qqplot_sorted_png'))
#                pyplot.close(fig)
#
#                fig = pyplot.figure()
#                pyplot.title('DIFF-MEAN-MAP Q-Q PLOT')
#                pyplot.plot([map_off-5*map_unc, map_off+5*map_unc], [-5, 5], 'b--')
#                pyplot.plot([-1, 1], [q_cut, q_cut], 'k-.')
#                pyplot.plot([-1, 1], [-q_cut, -q_cut], 'k-.')
#                pyplot.plot(sorted_observed_diff_vals, expected_diff_vals, 'go-')
#                pyplot.xlabel('DIFF-MEAN-MAP QUANTILES')
#                pyplot.ylabel('THEORETICAL QUANTILES')
#                # Apply tight layout to prevent overlaps
#                pyplot.tight_layout()
#                # Save
#                pyplot.savefig(d_handler.output_handler.get_file('unc_qqplot_png'))
#                pyplot.close(fig)
#
#        self.log('\rMap Uncertainties Calculated.                                         ', True)
#
#        try:
#            from ascii_graph import Pyasciigraph
#            g=Pyasciigraph()
#            graph_data = [(d.d_tag, round(d.get_map_uncertainty(),3)) for d in self.datasets.mask(mask_name='rejected - total', invert=True)]
#            for l in g.graph(label='Sorted Map Uncertainties (Ascending Order)', data=graph_data, sort=1):
#                print(l)
#        except ImportError:
#            print('IMPORT ERROR (ascii_graph) - CANNOT GENERATE UNCERTAINTY GRAPH')
#            pass

    def calculate_map_statistics(self):
        """Take the sampled maps and calculate statistics for each grid point across the datasets"""

        # Create statistics objects for each grid point
        self.log('===================================>>>', True)
        self.log('Calculating Statistics of Grid Points', True)

        # Get masked points in reference dataset
        masked_idxs = self.reference_grid().global_mask().outer_mask_indices()

        # All dataset handlers for analysing
        all_d_handlers = self.datasets.mask(mask_name='selected for analysis')
        d_uncertainties = [d.get_map_uncertainty() for d in all_d_handlers]

        # Chunk the points into groups - Compromise between cpu time and memory usage - ~200 dataset -> chunksize of 5000
        chunk_size = 1000*int(2000/len(all_d_handlers))
        chunked_indices = [masked_idxs[i:i + chunk_size] for i in range(0, len(masked_idxs), chunk_size)]
        chunk_num = len(chunked_indices)

        self.log('Iterating through {!s} points in {!s} chunks'.format(len(masked_idxs), chunk_num), True)
        t1 = time.time()

        # Statistics objects for each of the points we're interested in
        masked_point_statistics = []
        masked_point_adj_sigmas = []

        # Starting guess for the underlying sigma is raw_std/guess_factor
        guess_factor = 0.001

        # Calculate the statistics of the map values across the datasets
        # Calculate the mean map across the datasets
        for chunk_idx, chunk in enumerate(chunked_indices):
            status_bar(n=chunk_idx-1, n_max=chunk_num)
            p_map_vals = numpy.array([dh.morphed_map.select(chunk) for dh in all_d_handlers])

            if chunk_idx+1 < chunk_num:
                assert len(p_map_vals) == len(all_d_handlers)
                assert len(p_map_vals.T) == chunk_size

            point_stats = [basic_statistics(flex.double(map_vals.tolist())) for map_vals in p_map_vals.T]

            # Iterate through and, using the uncertainties of the maps, calculate the underlying standard deviation of the map values
            approx_sigmas = [bs.bias_corrected_standard_deviation for bs in point_stats]
            adjusted_sigmas = [estimate_true_underlying_sd(obs_vals=map_vals.tolist(), obs_error=d_uncertainties, est_sigma=approx_sigmas[i]*guess_factor) for i, map_vals in enumerate(p_map_vals.T)]

            # Extend the complete lists
            masked_point_statistics.extend(point_stats)
            masked_point_adj_sigmas.extend(adjusted_sigmas)

            assert i+1 == len(approx_sigmas), 'LIST INDEX DOES NOT MATCH LIST LENGTH'

        status_bar(n=chunk_idx, n_max=chunk_num)
        assert len(masked_point_statistics) == len(masked_idxs)
        assert len(masked_point_adj_sigmas) == len(masked_idxs)

        t2 = time.time()
        print('> MAP POINT ANALYSIS > Time Taken: {!s} seconds'.format(int(t2-t1)))

        # Calculate Stds Maps - Set the background to be tiny but non-zero so we can still divide by it
        stds_map_vals = 1e-18*numpy.ones(self.reference_grid().grid_size_1d())
        stds_map_vals.put(masked_idxs, [bs.bias_corrected_standard_deviation for bs in masked_point_statistics])
        stds_map_vals = flex.double(stds_map_vals.tolist())
        stds_map_file = self.output_handler.get_file('stds_map')
        self.write_array_to_map(output_file=stds_map_file, map_data=stds_map_vals)

        # Calculate ADJUSTED Stds Maps - Set the background to be tiny but non-zero so we can still divide by it
        adj_stds_map_vals = 1e-18*numpy.ones(self.reference_grid().grid_size_1d())
        adj_stds_map_vals.put(masked_idxs, masked_point_adj_sigmas)
        adj_stds_map_vals = flex.double(adj_stds_map_vals.tolist())
        adj_stds_map_file = self.output_handler.get_file('stds_adj_map')
        self.write_array_to_map(output_file=adj_stds_map_file, map_data=adj_stds_map_vals)

        # Calculate Skew Maps
        skew_map_vals = numpy.zeros(self.reference_grid().grid_size_1d())
        skew_map_vals.put(masked_idxs, [bs.skew for bs in masked_point_statistics])
        skew_map_vals = flex.double(skew_map_vals.tolist())
        skew_map_file = self.output_handler.get_file('skew_map')
        self.write_array_to_map(output_file=skew_map_file, map_data=skew_map_vals)

        # Calculate Kurtosis Maps
        kurt_map_vals = numpy.zeros(self.reference_grid().grid_size_1d())
        kurt_map_vals.put(masked_idxs, [bs.kurtosis for bs in masked_point_statistics])
        kurt_map_vals = flex.double(kurt_map_vals.tolist())
        kurt_map_file = self.output_handler.get_file('kurt_map')
        self.write_array_to_map(output_file=kurt_map_file, map_data=kurt_map_vals)

        # Calculate Bimodality Maps
        bimo_map_vals = numpy.zeros(self.reference_grid().grid_size_1d())
        bimo_map_vals.put(masked_idxs, [(bs.skew**2 + 1)/bs.kurtosis for bs in masked_point_statistics])
        bimo_map_vals = flex.double(bimo_map_vals.tolist())
        bimo_map_file = self.output_handler.get_file('bimo_map')
        self.write_array_to_map(output_file=bimo_map_file, map_data=bimo_map_vals)

        # Store map vals
        self.stat_maps.stds_map     = stds_map_vals
        self.stat_maps.adj_stds_map = adj_stds_map_vals
        self.stat_maps.skew_map     = skew_map_vals
        self.stat_maps.kurt_map     = kurt_map_vals
        self.stat_maps.bimo_map     = bimo_map_vals

#    def old_calculate_map_statistics(self):
#        """Take the sampled maps and calculate statistics for each grid point across the datasets"""
#
#        # Create statistics objects for each grid point
#        self.log('===================================>>>', True)
#        self.log('Calculating Statistics of Grid Points', True)
#
#        # Get masked points in reference dataset
#        masked_gps = self.reference_grid().global_mask().outer_mask()
#        masked_idxs = self.reference_grid().global_mask().outer_mask_indices()
#        masked_cart_ref = flex.vec3_double(masked_gps)*self.reference_grid().grid_spacing()
#        # Mapping of grid points to rotation matrix keys (residue CA labels)
#        masked_cart_mappings = self.reference_grid().partition().query_by_grid_indices(masked_idxs)
#
#        assert len(masked_idxs) == len(masked_cart_ref)
#
#        # All dataset handlers for analysing
#        all_d_handlers = self.datasets.mask(mask_name='rejected - total', invert=True)
#        d_uncertainties = [d.get_map_uncertainty() for d in all_d_handlers]
#
#        # Chunk the points into groups - Compromise between cpu time and memory usage - ~200 dataset -> chunksize of 5000
#        chunk_size = 1000*int(2000/self.datasets.size(mask_name='rejected - total', invert=True))
#        chunked_cart = [masked_cart_ref[i:i + chunk_size] for i in range(0, len(masked_cart_ref), chunk_size)]
#        chunked_mapping = [masked_cart_mappings[i:i + chunk_size] for i in range(0, len(masked_cart_mappings), chunk_size)]
#        chunk_num = len(chunked_cart)
#
#        self.log('Iterating through {!s} points in {!s} chunks'.format(len(masked_idxs), chunk_num), True)
#        t1 = time.time()
#
#        # Statistics objects for each of the points we're interested in
#        masked_point_statistics = []
#        masked_point_adj_sigmas = []
#
#        # Starting guess for the underlying sigma is raw_std/guess_factor
#        guess_factor = 0.001
#
#        # Calculate the statistics of the map values across the datasets
#        for chunk_idx, ref_coord_chunk in enumerate(chunked_cart):
#            status_bar(n=chunk_idx, n_max=chunk_num)
#            p_map_vals = numpy.array([dh.get_cart_values(points=dh.transform_from_reference(points=flex.vec3_double(ref_coord_chunk), method=self.get_alignment_method(), point_mappings=chunked_mapping[chunk_idx]), map_type=self.get_obs_map_type()) for dh in all_d_handlers])
#
#            if chunk_idx+1 < chunk_num:
#                assert len(p_map_vals) == len(all_d_handlers)
#                assert len(p_map_vals.T) == chunk_size
#
#            point_stats = [basic_statistics(flex.double(map_vals.tolist())) for map_vals in p_map_vals.T]
#
#            # Iterate through and, using the uncertainties of the maps, calculate the underlying standard deviation of the map values
#            approx_sigmas = [bs.bias_corrected_standard_deviation for bs in point_stats]
#            adjusted_sigmas = [estimate_true_underlying_sd(obs_vals=map_vals.tolist(), obs_error=d_uncertainties, est_sigma=approx_sigmas[i]*guess_factor) for i, map_vals in enumerate(p_map_vals.T)]
#
#            # Extend the complete lists
#            masked_point_statistics.extend(point_stats)
#            masked_point_adj_sigmas.extend(adjusted_sigmas)
#
#            assert i+1 == len(approx_sigmas), 'LIST INDEX DOES NOT MATCH LIST LENGTH'
#
#        t2 = time.time()
#        print('> MAP POINT ANALYSIS > Time Taken: {!s} seconds'.format(int(t2-t1)))
#
#        # Calculate Stds Maps - Set the background to be tiny but non-zero so we can still divide by it
#        stds_map_vals = 1e-18*numpy.ones(self.reference_grid().grid_size_1d())
#        stds_map_vals.put(masked_idxs, [bs.bias_corrected_standard_deviation for bs in masked_point_statistics])
#        stds_map_vals = flex.double(stds_map_vals.tolist())
#        stds_map_file = self.output_handler.get_file('stds_map')
#        self.write_array_to_map(output_file=stds_map_file, map_data=stds_map_vals)
#
#        # Calculate ADJUSTED Stds Maps - Set the background to be tiny but non-zero so we can still divide by it
#        adj_stds_map_vals = 1e-18*numpy.ones(self.reference_grid().grid_size_1d())
#        adj_stds_map_vals.put(masked_idxs, masked_point_adj_sigmas)
#        adj_stds_map_vals = flex.double(adj_stds_map_vals.tolist())
#        adj_stds_map_file = self.output_handler.get_file('stds_adj_map')
#        self.write_array_to_map(output_file=adj_stds_map_file, map_data=adj_stds_map_vals)
#
#        # Calculate Skew Maps
#        skew_map_vals = numpy.zeros(self.reference_grid().grid_size_1d())
#        skew_map_vals.put(masked_idxs, [bs.skew for bs in masked_point_statistics])
#        skew_map_vals = flex.double(skew_map_vals.tolist())
#        skew_map_file = self.output_handler.get_file('skew_map')
#        self.write_array_to_map(output_file=skew_map_file, map_data=skew_map_vals)
#
#        # Calculate Kurtosis Maps
#        kurt_map_vals = numpy.zeros(self.reference_grid().grid_size_1d())
#        kurt_map_vals.put(masked_idxs, [bs.kurtosis for bs in masked_point_statistics])
#        kurt_map_vals = flex.double(kurt_map_vals.tolist())
#        kurt_map_file = self.output_handler.get_file('kurt_map')
#        self.write_array_to_map(output_file=kurt_map_file, map_data=kurt_map_vals)
#
#        # Calculate Bimodality Maps
#        bimo_map_vals = numpy.zeros(self.reference_grid().grid_size_1d())
#        bimo_map_vals.put(masked_idxs, [(bs.skew**2 + 1)/bs.kurtosis for bs in masked_point_statistics])
#        bimo_map_vals = flex.double(bimo_map_vals.tolist())
#        bimo_map_file = self.output_handler.get_file('bimo_map')
#        self.write_array_to_map(output_file=bimo_map_file, map_data=bimo_map_vals)
#
#        # Store map vals
#        self.stat_maps.stds_map     = stds_map_vals
#        self.stat_maps.adj_stds_map = adj_stds_map_vals
#        self.stat_maps.skew_map     = skew_map_vals
#        self.stat_maps.kurt_map     = kurt_map_vals
#        self.stat_maps.bimo_map     = bimo_map_vals

    def collect_map_statistics(self):
        """Collect map statistics from the datasets"""

        for d_handler in self.datasets.all():

            if d_handler.z_map_stats:
                self.get_map_observations().set_data_value(data_name='z_map_mean', entry_id=d_handler.d_tag, value=d_handler.z_map_stats['z_map_mean'])
                self.get_map_observations().set_data_value(data_name='z_map_std', entry_id=d_handler.d_tag, value=d_handler.z_map_stats['z_map_std'])
                self.get_map_observations().set_data_value(data_name='z_map_skew', entry_id=d_handler.d_tag, value=d_handler.z_map_stats['z_map_skew'])
                self.get_map_observations().set_data_value(data_name='z_map_kurtosis', entry_id=d_handler.d_tag, value=d_handler.z_map_stats['z_map_kurtosis'])

    def calculate_z_map(self, map_vals, method='naive', map_uncertainty=None, sparse=False):
        """Takes a flex.double map and calculates the z-map relative to the mean and std map"""

        # STANDARD METHOD USES THE RAW STANDARD DEVIATION OF THE MAP VALUES
        # ADJUSTED METHOD USED THE UNCERTAINTY CORRECTED STANDARD DEVIATION
        assert method in ['naive','adjusted','adjusted+uncertainty']
        if method == 'adjusted+uncertainty': assert map_uncertainty, 'UNCERTAINTY REQUIRED TO USE ADJUSTED METHOD'

        assert map_vals.size() == self.reference_grid().grid_size_1d()

        # Calculate Z-values
        if method == 'naive':
            z_map_vals = (map_vals - self.stat_maps.mean_map)/self.stat_maps.stds_map
        elif method == 'adjusted':
            z_map_vals = (map_vals - self.stat_maps.mean_map)/self.stat_maps.adj_stds_map
        elif method == 'adjusted+uncertainty':
            z_map_vals = (map_vals - self.stat_maps.mean_map)/flex.sqrt(self.stat_maps.adj_stds_map**2 + map_uncertainty**2)

        if sparse:
            # Calculate the sparse vector of the masked map values
            return sparse.vector(self.reference_grid().grid_size_1d(),
                                 dict(zip(self.reference_grid().global_mask().outer_mask_indices(),
                                          z_map_vals.select(self.reference_grid().global_mask().outer_mask_indices()))))
        else:
            return z_map_vals

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

    def print_clustering_settings(self, z_cutoff,
                                        min_cluster_volume,
                                        clustering_cutoff,
                                        clustering_criterion,
                                        clustering_metric,
                                        clustering_method   ):

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
                                    point_mask,
                                    min_cluster_volume,
                                    clustering_cutoff,
                                    clustering_criterion,
                                    clustering_metric,
                                    clustering_method   ):
        """Finds all the points in the z-map above `z_cutoff`, points will then be clustered into groups of cutoff `clustering_cutoff` angstroms"""

        self.log('===================================>>>', True)

        # Check that the map is the expected size etc...
        assert len(z_map) == self.reference_grid().grid_size_1d()

        # Translate between grid coords and grid index
        grid_indexer = self.reference_grid().grid_indexer()

        # Scale the cutoff (Angstroms) into grid units
        grid_clustering_cutoff = clustering_cutoff/self.reference_grid().grid_spacing()

        # Calculate the approximate minimum volume for the cluster size
        min_cluster_size = int(min_cluster_volume/(self.reference_grid().grid_spacing()**3))

        # List of points in this dataset - TODO Revisit using global outer mask - is this the right mask?
        d_selected_points = [(gp, z_map[grid_indexer(gp)]) for gp in point_mask if z_map[grid_indexer(gp)] >= z_cutoff]
        d_handler.number_of_sig_points = len(d_selected_points)

        # Write the output mask if points found
        if d_selected_points:
            # Create maps of the high z-value points (significant points)
            highz_points = zip(*d_selected_points)[0]
            highz_map_array = numpy.zeros(self.reference_grid().grid_size_1d(), dtype=int)
            highz_map_array.put(map(self.reference_grid().grid_indexer(), highz_points), [1]*len(highz_points))
            self.write_array_to_map(d_handler.output_handler.get_file('high_z_mask'), flex.double(highz_map_array.tolist()))

        # No Cluster points found
        if not d_selected_points:
            self.log('Dataset {!s}: No Clusters Found'.format(d_handler.d_tag), True)
        # Can't cluster if there are too many points
        elif len(d_selected_points) > 10000:
            self.log('Dataset {!s}: Too many points to cluster: {!s} Points.'.format(d_handler.d_tag, len(d_selected_points)), True)
            # This dataset is too noisy to analyse - flag!
            self.datasets.all_masks().set_mask_value(mask_name='noisy zmap', entry_id=d_handler.d_tag, value=True)

            # Link datasets to the initial results directory
            noisy_dir = os.path.join(self.output_handler.get_dir('noisy_datasets'), d_handler.d_name)
            if not os.path.exists(noisy_dir):
                rel_symlink(orig=d_handler.output_handler.get_dir('root'), link=noisy_dir)

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
                    print('> Clustering > Time Taken: {!s} seconds'.format(int(t2-t1)))

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

            if clust_num == 0:
                self.log('===> No Clusters found - Minimum cluster size not reached.', True)
            else:
                # Check that the clusters are numbered properly
                assert clust_num == max(clust_dict.keys())
                assert clust_num == len(clust_dict.keys())

                # Add clustered points to the dataset handler
                self.log('===> {!s} Cluster(s) found.'.format(clust_num), True)

                # Create a link to the interesting directories in the initial results directory
                hit_dir = os.path.join(self.output_handler.get_dir('interesting_datasets'), d_handler.d_name)
                if not os.path.exists(hit_dir):
                    rel_symlink(orig=d_handler.output_handler.get_dir('root'), link=hit_dir)

                # Return a positive result
                return clust_dict

        # Return empty - catchall if it makes it to the end of the function
        return {}

    def collate_all_clusters(self):
        """Collate clusters from all of the datasets"""

        self.log('===================================>>>', False)
        self.log('Collating Clusters', False)

        # List of points to be returned
        all_dataset_clusters = dict([(d.d_tag, []) for d in self.datasets.all()])

        for d_handler in self.datasets.all():
            if d_handler.raw_cluster_hits:
                all_dataset_clusters[d_handler.d_tag] = d_handler.raw_cluster_hits

        # Print Cluster Summaries
        cluster_num = [(k, len(all_dataset_clusters[k])) for k in sorted(all_dataset_clusters.keys()) if all_dataset_clusters[k]]
        cluster_total = sum([a[1] for a in cluster_num])

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

        for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True):

            # Check to see if there are any clustered points
            if not d_handler.raw_cluster_hits:
                continue

            # This dataset is interesting!
            self.datasets.all_masks().set_mask_value(mask_name='interesting', entry_id=d_handler.d_tag, value=True)

            cluster_obj = d_handler.clustered_hits
            # Add to list of all clusters
            cluster_ids.append(d_handler.d_tag)
            all_clusters.append(cluster_obj)

            ########################################################

            # Calculate the spacings of the clusters
#            dists = []
#            total_c_num = len(cluster_obj.get_centroids())
#            for cent_a in cluster_obj.get_centroids():
#                for cent_b in cluster_obj.get_centroids():
#                    if cent_a == cent_b: continue
#                    dists.append((flex.double(cent_a) - flex.double(cent_b)).norm()*self.reference_grid().grid_spacing())
#            if not dists: print 'Clusters, Dataset {!s}:'.format(d_handler.d_tag), '\tNum: {:4}'.format(total_c_num), '\tMin Spacing: {!s:>5} A'.format('--.--'), '\tMax Spacing: {!s:>5} A'.format('--.--')
#            else:         print 'Clusters, Dataset {!s}:'.format(d_handler.d_tag), '\tNum: {:4}'.format(total_c_num), '\tMin Spacing: {:5.3} A'.format(min(dists)), '\tMax Spacing: {:5.3} A'.format(max(dists))

            ########################################################

            # Pull out the coordinates of the peaks, in cartesian, in the frame of the dataset
            peak_sites_cart_ref = list(flex.vec3_double(cluster_obj.get_peaks())*self.reference_grid().grid_spacing())
#            peak_sites_cart = list(d_handler.transform_from_reference(points=flex.vec3_double(peak_sites_cart_ref), method=self.get_alignment_method(), point_mappings=self.reference_dataset().find_nearest_calpha(peak_sites_cart_ref)))
            peak_sites_cart = list(d_handler.transform_from_reference(points=flex.vec3_double(peak_sites_cart_ref), method=self.params.alignment.method, point_mappings=self.reference_grid().partition().query_by_grid_points(cluster_obj.get_peaks())))
            list_to_write = zip(cluster_obj.get_maxima(), cluster_obj.get_sizes(), *zip(*[p+p_ref for p,p_ref in zip(peak_sites_cart, peak_sites_cart_ref)]))

            sort_order = cluster_obj.sort(sorting_data='values', sorting_function=max, decreasing=True)
            sorted_list = [list_to_write[i] for i in sort_order]

            string_to_write = '\n'.join(['z_peak, size, x, y, z, refx, refy, refz'] + [', '.join(map(str, l)) for l in sorted_list])
            with open(d_handler.output_handler.get_file('z_peaks_csv'), 'w') as fh:
                fh.write(string_to_write)

            ########################################################

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

        # TODO XXX EITHER CREATE A NEW FUNCTION OR EXTEND THE FUNCTION ABOVE XXX TODO
        # CLUSTER THE LARGE Z-VALUES ACROSS THE DATASETS TO DETECT INTERESTING 'SITES' ON THE PROTEIN

    def write_pymol_scripts(self, d_handler):
        """Autogenerate pymol scripts"""

        # Get template to be filled in
        template = PANDDA_HTML_ENV.get_template('load_pandda_maps.pml')

        with open(d_handler.output_handler.get_file(file_tag='pymol_script'), 'w') as out_pml:
            out_pml.write(template.render({'file_dict':d_handler.output_handler.output_files}))

    def image_blob(self, script, image, d_handler, point, point_no, towards=[10,10,10]):
        """Take pictures of the maps with ccp4mg"""

        from Bamboo.Common.Command import CommandManager
        from Giant.Graphics import calculate_view_quaternion, multiply_quaternions

        # Get the template to be filled in
        template = PANDDA_HTML_ENV.get_template('ccp4mg-pic.py')

        orientation = calculate_view_quaternion(towards, point)
        rotate_1 = multiply_quaternions(orientation, (0.0, 0.5**0.5, 0.0, 0.5**0.5))
        rotate_2 = multiply_quaternions(orientation, (0.5**0.5, 0.0, 0.0, 0.5**0.5))

        for view_no, view in enumerate([orientation, rotate_1, rotate_2]):

            view_script = script.format(point_no, view_no)
            view_image  = image.format(point_no, view_no)

            print view_no, view_script, view_image

            ccp4mg_script = template.render({
                                                'view'  :{
                                                                'camera_centre' : [-1*c for c in point],
                                                                'orientation'   : list(view)
                                                            },
                                                'mol'   :{
                                                                'path'  : d_handler.output_handler.get_file('aligned_structure'),
                                                                'name'  : 'aligned_structure'
                                                            },
                                         #       'map'   :{
                                         #                       'path'    : d_handler.output_handler.get_file('sampled_map'),
                                         #                       'name'    : 'sampled_map',
                                         #                       'contour' : [1]
                                         #                   },
                                                'diff_map' :{
                                                                'path'    : d_handler.output_handler.get_file('z_map_corrected_normalised'),
                                                                'name'    : 'diff_map',
                                                            #    'neg-contour' : -3,
                                                                'pos-contour' : [3,4,5,6,7]
                                                            }
                                            })

            # Write out the ccp4mg script to the dataset's scripts folder
            with open(view_script, 'w') as fh:
                fh.write(ccp4mg_script)

            # Make the images
            c = CommandManager('ccp4mg')
            c.SetArguments(['-norestore','-picture', view_script, '-R', view_image, '-RO', """'{"size":"1000x1000"}'""", '-quit'])
            c.Run()

            if not os.path.exists(view_image):
                print 'FAILED TO MAKE IMAGES'
                print c.err

    def write_analysis_summary(self):
        """Writes an html summary of the datasets"""

        # Get template to be filled in
        template = PANDDA_HTML_ENV.get_template('output_table.html')
        # XXX 0=success, 1=none, 2=info, 3=warning, 4=failure

        d_summary = self.datasets_summary

        # Construct the data object to populate the template
        output_data = {'PANDDA_TOP' : PANDDA_TOP}
        output_data['header'] = 'PANDDA Processing Output'
        output_data['title'] = 'PANDDA Processing Output'
        output_data['introduction'] = 'Summary of Processing of Datasets'
        output_data['table'] = {}
        output_data['table']['column_headings'] = ['Data Quality', 'Identical Structure', 'Model Quality', 'RMSD to Reference', 'Overall', 'Interesting Areas']
        output_data['table']['rows'] = []
        # Add the datasets as rows
        for d in self.datasets.all():

            rmsd = round(d_summary.get_data_value(data_name='rmsd_to_mean', entry_id=d.d_tag), 3)
            rfree = round(d_summary.get_data_value(data_name='rfree', entry_id=d.d_tag), 3)
            rwork = round(d_summary.get_data_value(data_name='rwork', entry_id=d.d_tag), 3)

            columns = []
            overall_success = [0]
            # ------------------------------>>>
            # Test for Data Quality
            # ------------------------------>>>
            if self.datasets.all_masks().get_mask_value(mask_name='rejected - crystal', entry_id=d.d_tag) == True:
                columns.append({'flag':4,'message':'Rejected'.format(None)})
            else:
                columns.append({'flag':0,'message':'OK'})
            # ------------------------------>>>
            # Test for Identical Structures
            # ------------------------------>>>
            if self.datasets.all_masks().get_mask_value(mask_name='rejected - structure', entry_id=d.d_tag) == True:
                columns.append({'flag':4,'message':'Rejected'})
            else:
                columns.append({'flag':0,'message':'OK'})
            # ------------------------------>>>
            # Test for Refinement Success - some test on r-free
            # ------------------------------>>>
            if self.datasets.all_masks().get_mask_value(mask_name='bad crystal - rfree', entry_id=d.d_tag) == True:
                columns.append({'flag':4,'message':'Free: {!s}'.format(rfree)})
            else:
                columns.append({'flag':0,'message':'Free: {!s}'.format(rfree)})
            # ------------------------------>>>
            # Test for Structure movement
            # ------------------------------>>>
            if self.datasets.all_masks().get_mask_value(mask_name='bad crystal - isomorphous structure', entry_id=d.d_tag) == True:
                columns.append({'flag':4,'message':'RMSD: {!s}'.format(rmsd)})
            else:
                columns.append({'flag':0,'message':'RMSD: {!s}'.format(rmsd)})
            # ------------------------------>>>
            # Test for Structure movement
            # ------------------------------>>>
            if self.datasets.all_masks().get_mask_value(mask_name='selected for analysis', entry_id=d.d_tag) == True:
                columns.append({'flag':0,'message':'Analysed'.format(rmsd)})
            else:
                columns.append({'flag':4,'message':'Rejected'.format(rmsd)})
            # ------------------------------>>>
            # Test for if it's interesting
            # ------------------------------>>>
            if self.datasets.all_masks().get_mask_value(mask_name='interesting', entry_id=d.d_tag) == True:
                columns.append({'flag':0,'message':'Clusters Found'})
            else:
                columns.append({'flag':1,'message':''})

            # Find largest error
            overall_success = max([c['flag'] for c in columns])

            output_data['table']['rows'].append({'heading':'Dataset {!s}'.format(d.d_tag),
                                                 'success':overall_success,
                                                 'columns':columns})

        with open(self.output_handler.get_file(file_tag='summary_table'), 'w') as out_html:
            out_html.write(template.render(output_data))

    def write_map_value_distribution(self, map_vals, grid_indices, output_file, plot_normal=False):
        """Write out the value distribution for a map"""

        try:
            import matplotlib
            matplotlib.interactive(0)
            from matplotlib import pyplot
        except:
            return

        plot_vals = [map_vals[i] for i in grid_indices]

        fig = pyplot.figure()
        pyplot.title('MAP VALUE DISTRIBUTION')
        pyplot.hist(x=plot_vals, bins=30, normed=True)
        if plot_normal:
            # Plot the distribution for N(0,1)
            nd_t = normal_distribution()
            theor_x = numpy.linspace(-5,5,101)
            theor_y = [nd_t.pdf(x) for x in theor_x]
            pyplot.plot(theor_x, theor_y, c='k', ls='--', marker='o')
            # Plot the distribution for the observed distribution
            nd_o = normal_distribution(mean=numpy.mean(plot_vals), sd=numpy.std(plot_vals))
            obs_x = numpy.linspace(-5,5,101)
            obs_y = [nd_o.pdf(x) for x in obs_x]
            pyplot.plot(obs_x, obs_y, c='g', ls='-', marker='o')

        pyplot.xlabel('MAP VALUE')
        pyplot.ylabel('DENSITY')
        # Apply tight layout to prevent overlaps
        pyplot.tight_layout()
        # Save
        pyplot.savefig(output_file)
        pyplot.close(fig)

    def write_summary_graphs(self):
        """Write out graphs of dataset variables"""

        try:
            import matplotlib
            matplotlib.interactive(0)
            from matplotlib import pyplot
            assert not pyplot.isinteractive(), 'PYPLOT IS STILL IN INTERACTIVE MODE'
        except:
            return

        def filter_nans(vals):
            return [v for v in vals if v is not numpy.nan]

        # Get the output directory to write the graphs into
        out_dir = self.output_handler.get_dir('analyses')

        n_bins = 30

        ########################################################

        self.log('===================================>>>')
        self.log('Generating Summary Graphs')

        ########################################################

        self.log('=> Crystal Variation')

        # Dataset Crystal Summary
        d_summary = self.datasets_summary

        high_res = d_summary.get_data(data_name='high_res_limit')
        low_res = d_summary.get_data(data_name='low_res_limit')
        rfree = d_summary.get_data(data_name='rfree')
        rwork = d_summary.get_data(data_name='rwork')
        rmsds = d_summary.get_data(data_name='rmsd_to_mean')
        a,b,c,alpha,beta,gamma = zip(*d_summary.get_data(data_name='cell_params'))
        vols = d_summary.get_data(data_name='cell_volume')

        ########################################################

        # RESOLUTIONS
        fig = pyplot.figure()
        pyplot.title('RESOLUTION HISTOGRAMS')
        # High Resolution
        pyplot.subplot(2, 1, 1)
        pyplot.hist(x=filter_nans(high_res), bins=n_bins)
        pyplot.xlabel('HIGH RESOLUTION LIMIT (A)')
        pyplot.ylabel('COUNT')
        # Low Resolution
        pyplot.subplot(2, 1, 2)
        pyplot.hist(x=filter_nans(low_res), bins=n_bins)
        pyplot.xlabel('LOW RESOLUTION LIMIT (A)')
        pyplot.ylabel('COUNT')
        # Apply tight layout to prevent overlaps
        pyplot.tight_layout()
        # Save both
        pyplot.savefig(os.path.join(out_dir, 'd_resolutions.png'))
        pyplot.close(fig)

        # R-FACTORS
        fig = pyplot.figure()
        pyplot.title('R-FACTOR HISTOGRAMS')
        # RFree
        pyplot.subplot(2, 1, 1)
        pyplot.hist(x=filter_nans(rfree), bins=n_bins)
        pyplot.xlabel('RFREE')
        pyplot.ylabel('COUNT')
        # RWork
        pyplot.subplot(2, 1, 2)
        pyplot.hist(x=filter_nans(rwork), bins=n_bins)
        pyplot.xlabel('RWORK')
        pyplot.ylabel('COUNT')
        # Apply tight layout to prevent overlaps
        pyplot.tight_layout()
        # Save both
        pyplot.savefig(os.path.join(out_dir, 'd_rfactors.png'))
        pyplot.close(fig)

        # RMSDs
        fig = pyplot.figure()
        pyplot.title('RMSDS TO MEAN STRUCTURE HISTOGRAM')
        pyplot.hist(x=[v for v in rmsds if v > -1], bins=n_bins)
        pyplot.xlabel('RMSD (A)')
        pyplot.ylabel('COUNT')
        # Apply tight layout to prevent overlaps
        pyplot.tight_layout()
        # Save
        pyplot.savefig(os.path.join(out_dir, 'd_rmsd_to_mean.png'))
        pyplot.close(fig)

        # CELL PARAMS
        fig = pyplot.figure()
        pyplot.title('UNIT CELL PARAMS')
        # A
        pyplot.subplot(2, 3, 1)
        pyplot.hist(x=filter_nans(a), bins=n_bins)
        pyplot.xlabel('A (A)')
        pyplot.ylabel('COUNT')
        # B
        pyplot.subplot(2, 3, 2)
        pyplot.hist(x=filter_nans(b), bins=n_bins)
        pyplot.xlabel('B (A)')
        # C
        pyplot.subplot(2, 3, 3)
        pyplot.hist(x=filter_nans(c), bins=n_bins)
        pyplot.xlabel('C (A)')
        # ALPHA
        pyplot.subplot(2, 3, 4)
        pyplot.hist(x=filter_nans(alpha), bins=n_bins)
        pyplot.xlabel('ALPHA')
        pyplot.ylabel('COUNT')
        # BETA
        pyplot.subplot(2, 3, 5)
        pyplot.hist(x=filter_nans(beta), bins=n_bins)
        pyplot.xlabel('BETA')
        # GAMMA
        pyplot.subplot(2, 3, 6)
        pyplot.hist(x=filter_nans(gamma), bins=n_bins)
        pyplot.xlabel('GAMMA')
        # Apply tight layout to prevent overlaps
        pyplot.tight_layout()
        # Save both
        pyplot.savefig(os.path.join(out_dir, 'd_cell_param.png'))
        pyplot.close(fig)

        # CELL VOLUME
        fig = pyplot.figure()
        pyplot.title('UNIT CELL VOLUME HISTOGRAM')
        pyplot.hist(x=filter_nans(vols), bins=n_bins)
        pyplot.xlabel('VOLUME (A**3)')
        pyplot.ylabel('COUNT')
        # Apply tight layout to prevent overlaps
        pyplot.tight_layout()
        # Save
        pyplot.savefig(os.path.join(out_dir, 'd_cell_volume.png'))
        pyplot.close(fig)

        ########################################################

        self.log('=> Statistical Maps')

        masked_idxs = self.reference_grid().global_mask().outer_mask_indices()

        mean_map_vals = list(self.stat_maps.mean_map.select(masked_idxs))
        stds_map_vals = list(self.stat_maps.stds_map.select(masked_idxs))
        ajsd_map_vals = list(self.stat_maps.adj_stds_map.select(masked_idxs))

        ########################################################

        # STATISTICAL MAP VALUES
        fig = pyplot.figure()
        pyplot.title('STATISTICAL MAP VALUES')
        # MEAN MAP
        pyplot.hist(x=list(mean_map_vals), bins=n_bins)
        pyplot.xlabel('MEAN MAP DISTRIBUTION')
        pyplot.ylabel('COUNT')
        # Apply tight layout to prevent overlaps
        pyplot.tight_layout()
        # Save
        pyplot.savefig(os.path.join(out_dir, 'mean_map_vals.png'))
        pyplot.close(fig)

        # STATISTICAL MAP VALUES
        fig = pyplot.figure()
        pyplot.title('STATISTICAL MAP VALUES')
        # STANDARD DEVIATION MAPS
        pyplot.subplot(2, 1, 1)
        pyplot.hist(x=list(stds_map_vals), bins=n_bins)
        pyplot.xlabel('STDS MAP DISTRIBUTION')
        pyplot.ylabel('COUNT')
        # STANDARD DEVIATION MAPS
        pyplot.subplot(2, 1, 2)
        pyplot.hist(x=list(ajsd_map_vals), bins=n_bins)
        pyplot.xlabel('ADJUSTED STDS MAP DISTRIBUTION')
        pyplot.ylabel('COUNT')
        # Apply tight layout to prevent overlaps
        pyplot.tight_layout()
        # Save both
        pyplot.savefig(os.path.join(out_dir, 'stds_map_vals.png'))
        pyplot.close(fig)

        ########################################################

        self.log('=> Map Variation')

        # Dataset Map Summary
        map_summary = self.get_map_observations()

        map_uncties = map_summary.get_data(data_name='masked_map_uncertainty')
        z_map_mean  = map_summary.get_data(data_name='z_map_mean')
        z_map_std = map_summary.get_data(data_name='z_map_std')
        z_map_skew  = map_summary.get_data(data_name='z_map_skew')
        z_map_kurtosis = map_summary.get_data(data_name='z_map_kurtosis')

        ########################################################

        if filter_nans(map_uncties):

            # MAP PARAMS
            fig = pyplot.figure()
            pyplot.title('MAP STATISTICS')
            # MAP UNCERTAINTIES
            pyplot.hist(x=filter_nans(map_uncties), bins=n_bins)
            pyplot.xlabel('UNCERTAINTY OF MAP VALUES')
            pyplot.ylabel('COUNT')
            # Apply tight layout to prevent overlaps
            pyplot.tight_layout()
            # Save
            pyplot.savefig(os.path.join(out_dir, 'd_map_uncertainties.png'))
            pyplot.close(fig)

            ########################################################

            self.log('=> Scatter Plots')

            # MAP RESOLUTION V UNCERTAINTY
            fig = pyplot.figure()
            pyplot.title('HIGH RES LIMIT AGAINST MAP UNCERTAINTY')
            pyplot.scatter(x=high_res, y=map_uncties)
            pyplot.xlabel('RESOLUTION')
            pyplot.ylabel('UNCERTAINTY')
            # Apply tight layout to prevent overlaps
            pyplot.tight_layout()
            # Save
            pyplot.savefig(os.path.join(out_dir, 'resolution_v_uncertainty.png'))
            pyplot.close(fig)

            # MAP RESOLUTION V UNCERTAINTY
            fig = pyplot.figure()
            pyplot.title('HIGH RES LIMIT AGAINST RFREE')
            pyplot.scatter(x=high_res, y=rfree)
            pyplot.xlabel('RESOLUTION')
            pyplot.ylabel('RFREE')
            # Apply tight layout to prevent overlaps
            pyplot.tight_layout()
            # Save
            pyplot.savefig(os.path.join(out_dir, 'resolution_v_rfree.png'))
            pyplot.close(fig)

            # RFREE V UNCERTAINTY
            fig = pyplot.figure()
            pyplot.title('RFREE AGAINST UNCERTAINTY')
            pyplot.scatter(x=rfree, y=map_uncties)
            pyplot.xlabel('RFREE')
            pyplot.ylabel('UNCERTAINTY')
            # Apply tight layout to prevent overlaps
            pyplot.tight_layout()
            # Save
            pyplot.savefig(os.path.join(out_dir, 'rfree_v_uncertainty.png'))
            pyplot.close(fig)

        # Check to see if any map values have been added
        if filter_nans(z_map_skew):

            self.log('===================================>>>')
            self.log('Writing Z-MAP Graphs')

            # R-FACTORS
            fig = pyplot.figure()
            pyplot.title('Z-MAP DISTRIBUTION HISTOGRAMS')
            # RFree
            pyplot.subplot(2, 1, 1)
            pyplot.hist(x=filter_nans(z_map_mean), bins=n_bins)
            pyplot.xlabel('Z-MAP MEAN')
            pyplot.ylabel('COUNT')
            # RWork
            pyplot.subplot(2, 1, 2)
            pyplot.hist(x=filter_nans(z_map_std), bins=n_bins)
            pyplot.xlabel('Z_MAP_STD')
            pyplot.ylabel('COUNT')
            # Apply tight layout to prevent overlaps
            pyplot.tight_layout()
            # Save both
            pyplot.savefig(os.path.join(out_dir, 'z_map_statistics.png'))
            pyplot.close(fig)

            # Z-MAP SKEW V UNCERTAINTY
            fig = pyplot.figure()
            pyplot.title('DATASET NORMALITY')
            pyplot.scatter(x=z_map_skew, y=z_map_kurtosis)
            pyplot.xlabel('SKEW')
            pyplot.ylabel('KURTOSIS')
            # Apply tight layout to prevent overlaps
            pyplot.tight_layout()
            # Save
            pyplot.savefig(os.path.join(out_dir, 'z_map_skew_v_kurtosis.png'))
            pyplot.close(fig)

    def write_summary_csvs(self):
        """Write CSV file of dataset variables"""

        self.log('===================================>>>')
        self.log('Writing Dataset Crystal Summaries')

        d_summary = self.datasets_summary

        with open(self.output_handler.get_file('dataset_summaries'), 'w') as fh:
            fh.write('dataset_id, res_high, res_low, rfree, rwork, rmsd to mean, a, b, c, alpha, beta, gamma, cell volume\n')
            # Write out parameters for each dataset
            for d_handler in self.datasets.all():
                out_list = [  d_handler.d_tag,
                              d_summary.get_data_value(data_name='high_res_limit', entry_id=d_handler.d_tag),
                              d_summary.get_data_value(data_name='low_res_limit', entry_id=d_handler.d_tag),
                              d_summary.get_data_value(data_name='rfree', entry_id=d_handler.d_tag),
                              d_summary.get_data_value(data_name='rwork', entry_id=d_handler.d_tag),
                              d_summary.get_data_value(data_name='rmsd_to_mean', entry_id=d_handler.d_tag)  ] + \
                           list(d_summary.get_data_value(data_name='cell_params', entry_id=d_handler.d_tag)) + \
                           [  d_summary.get_data_value(data_name='cell_volume', entry_id=d_handler.d_tag)  ]
                out_line = ', '.join(map(str,out_list)) + '\n'
                fh.write(out_line)

        self.log('===================================>>>')
        self.log('Writing Dataset Map Summaries')

        map_summary = self.get_map_observations()

        with open(self.output_handler.get_file('map_summaries'), 'w') as fh:
            fh.write('dataset_id, masked_map_mean, masked_map_std, masked_map_uncertainty, z_map_mean, z_map_std, z_map_skew, z_map_kurtosis\n')
            # Write out parameters for each dataset
            for d_handler in self.datasets.mask(mask_name='selected for analysis'):
                out_list = [  d_handler.d_tag,
                              map_summary.get_data_value(data_name='masked_map_mean', entry_id=d_handler.d_tag),
                              map_summary.get_data_value(data_name='masked_map_std', entry_id=d_handler.d_tag),
                              map_summary.get_data_value(data_name='masked_map_uncertainty', entry_id=d_handler.d_tag),
                              map_summary.get_data_value(data_name='z_map_mean', entry_id=d_handler.d_tag),
                              map_summary.get_data_value(data_name='z_map_std', entry_id=d_handler.d_tag),
                              map_summary.get_data_value(data_name='z_map_skew', entry_id=d_handler.d_tag),
                              map_summary.get_data_value(data_name='z_map_kurtosis', entry_id=d_handler.d_tag)  ]

                out_line = ', '.join(map(str,out_list)) + '\n'
                fh.write(out_line)

        self.log('===================================>>>')
        self.log('Writing Statistical Map Summaries')

        gps_str = ['\"'+str(g)+'\"' for g in self.reference_grid().global_mask().outer_mask()]

        with open(self.output_handler.get_file('stats_map_summaries'), 'w') as fh:
            fh.write('gp, mean, std, adj_std, skew, kurtosis, bimodality\n')
            # Write out map values for each grid point
            for g_str, g_idx in zip(gps_str, self.reference_grid().global_mask().outer_mask_indices()):

                out_list = [  g_str,
                              self.stat_maps.mean_map[g_idx],
                              self.stat_maps.stds_map[g_idx],
                              self.stat_maps.adj_stds_map[g_idx],
                              self.stat_maps.skew_map[g_idx],
                              self.stat_maps.kurt_map[g_idx],
                              self.stat_maps.bimo_map[g_idx]  ]

                out_line = ', '.join(map(str,out_list)) + '\n'
                fh.write(out_line)

        # TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
        #self.log('===================================>>>')
        #self.log('Writing Residue Summaries - NOT YET')
        # TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO

    def write_ranked_cluster_csv(self, cluster, sorted_indices, outfile):
        """Output a combined list of all of identified blobs, for all datasets"""

        output_list = []

        # Blob Rank, Blob Index
        for c_rank, c_index in enumerate(sorted_indices):
            # Extract the label of the cluster
            d_tag, c_num = cluster.get_keys([c_index])[0]
            # Pull out the d_handler for the dataset
            d_handler = self.datasets.get(d_tag=d_tag)

            # Extract the location of the blob in the reference frame
            grid_ref = cluster.get_peaks([c_index])[0]
            cart_ref = list(flex.vec3_double([grid_ref])*self.reference_grid().grid_spacing())[0]
            # Transform it to the original frame of reference
            cart_nat = list(d_handler.transform_from_reference(
                                                                points=flex.vec3_double([cart_ref]),
                                                                method=self.params.alignment.method,
#                                                                point_mappings=self.reference_dataset().find_nearest_calpha([cart_ref])
                                                                point_mappings=self.reference_grid().partition().query_by_grid_points([grid_ref])
                                                                )
                            )[0]
            # Extract the peak value of the blob
            cluster_peak_val = cluster.get_maxima([c_index])[0]
            cluster_size_val = cluster.get_sizes([c_index])[0]

            # Combine into columns for the output
            cols_to_write = [d_tag, c_rank+1, cluster_peak_val, cluster_size_val] + list(cart_nat) + list(cart_ref) + [d_handler.get_pdb_filename(), d_handler.get_mtz_filename()]
            output_list.append(cols_to_write)

        # Form the headers, and write
        headers = 'dtag, rank, blob_peak, blob_size, x, y, z, refx, refy, refz, pdb, mtz'
        string_to_write = '\n'.join([headers] + [', '.join(map(str, l)) for l in output_list])
        with open(outfile, 'w') as fh:
            fh.write(string_to_write)

        return None

    def write_grid_point_distributions(self, grid_points, output_filename=None):
        """Write CSV file of grid points, dataset numbers and map values"""

        if not output_filename: output_filename = self.output_handler.get_file('point_distributions')

        self.log('===================================>>>', True)
        self.log('Writing Grid Points distributions: {!s} points'.format(len(grid_points)), True)
        self.log('Writing to {!s}'.format(output_filename), True)

        # Add quotations around the grid points
        grid_points_str = ['\"'+str(g)+'\"' for g in grid_points]

        # Iterate through each dataset
        with open(output_filename, 'w') as fh:

            fh.write('grid_point, dataset_tag, dataset_num, map_val\n')
            for d_handler in self.datasets.mask(mask_name='selected for analysis'):
                map_vals = [d_handler.morphed_map[g] for g in grid_points]
                # Create list of tags for zipping
                d_tags = [d_handler.d_tag]*len(grid_points)
                d_nums = [d_handler.d_num]*len(grid_points)
                # Format and write
                out_list = [', '.join(map(str,tup)) for tup in zip(grid_points_str, d_tags, d_nums, map_vals)]
                out_line = '\n'.join(map(str,out_list)) + '\n'
                fh.write(out_line)

    def write_array_to_map(self, output_file, map_data, grid_size=None, grid_spacing=None):
        """Takes a 1d array and writes it to a map"""
        if not grid_size:    grid_size    = self.reference_grid().grid_size()
        if not grid_spacing: grid_spacing = self.reference_grid().grid_spacing()
        self.log('> Writing Map: {!s}'.format(output_file))
        write_1d_array_as_p1_map(file_name=output_file, map_data=map_data, grid_size=grid_size, grid_spacing=grid_spacing)

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


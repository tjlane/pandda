from __future__ import print_function

import os, sys, glob, time, re
import copy, resource, gc
import multiprocessing
import warnings

import math

import scipy.spatial
import scipy.cluster

import numpy, pandas

import iotbx.pdb
import iotbx.mtz
import iotbx.ccp4_map
import iotbx.crystal_symmetry_from_any
#import iotbx.map_tools as map_tools

from libtbx import easy_mp, easy_pickle
from libtbx import phil

import cctbx.miller
import cctbx.maptbx

import scitbx.sparse
import scitbx.matrix

from scitbx.array_family import flex
from scitbx.math import superpose, basic_statistics
from scitbx.math.distributions import normal_distribution

from libtbx.math_utils import ifloor, iceil

#from iotbx.reflection_file_utils import extract_miller_array_from_file

from Bamboo.Common import Meta, Info
from Bamboo.Common.Logs import Log
from Bamboo.Common.File import output_file_object, easy_directory
from Bamboo.Common.Data import data_collection, multiple_data_collection
from Bamboo.Common.Masks import mask_collection

#from Bamboo.Density.Edstats.Score import score_file_with_edstats
from Bamboo.Density.Edstats.Utils import score_with_edstats_to_dict

from Giant.Grid import grid_handler
from Giant.Grid.Masks import spherical_mask, atomic_mask, grid_mask

from Giant.Xray.Miller.Scaling import apply_simple_scaling
from Giant.Xray.Maps.Utils import write_1d_array_as_p1_map
from Giant.Xray.Data import crystalSummary
from Giant.Xray.Data.utils import extract_structure_factors
from Giant.Xray.Symmetry import combine_hierarchies, generate_adjacent_symmetry_copies

from Giant.Stats.Ospina import estimate_true_underlying_sd

from Giant.Stats.Cluster import cluster_data, combine_clusters, find_connected_groups

from Giant.Structure.Align import perform_flexible_alignment

from Giant.Utils import status_bar, rel_symlink

from PANDDAs import graphs

from PANDDAs import PANDDA_VERSION
from PANDDAs.phil import pandda_phil_def
from PANDDAs.html import PANDDA_HTML_ENV
from PANDDAs.settings import PANDDA_TOP, PANDDA_TEXT

STRUCTURE_MASK_NAMES = [    'bad structure - chain counts',
                            'bad structure - chain ids',
                            'bad structure - chain sequences',
                            'bad structure - residue counts',
                            'bad structure - atom counts',
                            'bad structure - non-identical structures',
                            'bad structure - different space group'    ]
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
                            'interesting',
                            'analysed'  ]

DATASET_INFO_FIELDS  = [    'high_resolution',
                            'low_resolution',
                            'r_work',
                            'r_free',
                            'rmsd_to_mean',
                            'space_group',
                            'uc_a',
                            'uc_b',
                            'uc_c',
                            'uc_alpha',
                            'uc_beta',
                            'uc_gamma',
                            'uc_vol'                ]
DATASET_MAP_FIELDS   = [    'analysed_resolution',
                            'map_uncertainty',
                            'obs_map_mean',
                            'obs_map_rms',
                            'z_map_mean',
                            'z_map_std',
                            'z_map_skew',
                            'z_map_kurt'            ]
DATASET_EVENT_FIELDS = [    'site_idx',
                            'est_occupancy',
                            'z_peak',
                            'z_mean',
                            'cluster_size',
                            'x','y','z',
                            'refx','refy','refz'    ]
SITE_TABLE_FIELDS    = [    'centroid',
                            'num_events',
                            'nearest_residue 1',
                            'nearest_residue 2',
                            'nearest_residue 3',
                            'native_centroid',
                            'near_crystal_contacts'     ]

# SET FILTERING CONSTANTS
FILTER_SCALING_CORRELATION_CUTOFF = 0.7

class HolderList(object):
    """Class for grouping many holders together"""
    _holder_class = None

    def __init__(self, *args):
        self._holder_list = []
        self._masks = mask_collection()
        # Allow for subclassing
        self.__class_init__(args)

    def __class_init__(self, args):
        pass

    def __getitem__(self, idx):
        if   isinstance(idx, int): return self.get(num=idx)
        elif isinstance(idx, str): return self.get(tag=idx)
        else: raise Exception('CANNOT INDEX EXCEPT BY int OR str. TYPE GIVEN: {!s}'.format(type(idx)))

    def __call__(self):
        """Return all holders"""
        return self.all()

    def all(self):
        """Return all holders"""
        return self._holder_list
    def all_nums(self):
        """Return the list of holder ids"""
        return [h.num for h in self.all()]
    def all_tags(self):
        """Return the list of holder tags"""
        return [h.tag for h in self.all()]

    def all_masks(self):
        """Return the mask object"""
        return self._masks
    def mask(self, mask_name, invert=False):
        """Retrieve a masked list of datasets"""
        if not mask_name:
            return self.all()
        else:
            return self._masks.mask(mask_name=mask_name, input_list=self.all(), invert=invert)

    def size(self, mask_name=None, invert=False):
        """Return the number of holders in the list (with optional mask applied)"""
        if mask_name:   return len(self.mask(mask_name=mask_name, invert=invert))
        else:           return len(self.all())

    def add(self, new_holders):
        """Add new datasets"""

        for new_h in new_holders:
            # Check all added datasets are the right class
            if self._holder_class:
                assert isinstance(new_h, self._holder_class), 'OBJECTS MUST BE OF TYPE: {!s}\n(ADDED OF TYPE: {!s})'.format(self._holder_class, type(new_h))
            # Check all added dataset id tags are strs
            assert isinstance(new_h.tag, str), 'TAG MUST BE str. Type given: {!s}'.format(type(new_h.tag))
            assert new_h.tag not in self.all_tags(), 'HOLDER WITH TAG ALREADY EXISTS: {!s}'.format(new_h.tag)
            # Check all added dataset id nums are ints
            assert isinstance(new_h.num, int), 'NUM MUST BE int. Type given: {!s}'.format(type(new_h.num))
            assert new_h.num not in self.all_nums(), 'HOLDER WITH NUM ALREADY EXISTS: {!s}'.format(new_h.num)
            # No problems, add to list
            self._holder_list.append(new_h)

    def get(self, tag=None, num=None):
        """Get a dataset by tag or num"""

        assert [num, tag].count(None) == 1, 'Must give EITHER num OR tag'
        if num: matching = [m for m in self.all() if m.num == num]
        else:   matching = [m for m in self.all() if m.tag == tag]
        if len(matching) == 0: raise Exception('NO MATCHING HOLDER FOUND - NUM: {!s}, TAG: {!s}'.format(num, tag))
        if len(matching) != 1: raise Exception('MORE THAN ONE MATCHING HOLDER FOUND - NUM: {!s}, TAG: {!s}'.format(num, tag))
        return matching[0]

class DatasetHandler(object):
    child = None
    def __init__(self, dataset_number, pdb_filename, mtz_filename, dataset_tag=None, name_prefix=''):
        """Create a dataset object to allow common functions to be applied easily to a number of datasets"""

        assert os.path.exists(pdb_filename), 'PDB file does not exist!'
        assert os.path.exists(mtz_filename), 'MTZ file does not exist!'

        # Store dataset number
        self.num = dataset_number
        # Store the tag for the dataset
        if dataset_tag:
            self.tag = dataset_tag
        else:
            # If num < 0 - mark as a reference dataset
            if self.num < 0:  self.tag = 'REF{:05d}'.format(self.num)
            else:             self.tag = 'D{:05d}'.format(self.num)
        # Store a name for the dataset
        self.name = name_prefix + '{!s}'.format(self.tag)

        # Output Directories
        self.output_handler = None

        ########################################################

        # Store filenames
        self._pdb_file = os.path.abspath(pdb_filename)
        self._mtz_file = os.path.abspath(mtz_filename)

        # PDB Objects
        self._structure = self.new_structure()

        # Data summaries
        self.pdb_summary = None
        self.mtz_summary = crystalSummary(self._mtz_file)

        ########################################################

        # All Structure factors
        self.sfs = None
        # Truncated structure factors
        self.tr_sfs = None

        # Initialise other variables
        self.unit_cell = None
        self.space_group = None

        # Single matrix - global alignment
        self._global_rt_transform = None
        # Multiple matrices - local alignment
        self._local_rt_transforms = None

        ########################################################

        # Map of the clusters in the dataset
        self.events = []

    def initialise_output_directory(self, outputdir):
        """Initialise a dataset output directory"""
        # Create a file and directory organiser
        self.output_handler = output_file_object(rootdir=easy_directory(outputdir))

    def pdb_filename(self):
        return self._pdb_file
    def mtz_filename(self):
        return self._mtz_file

    def input(self):
        return self._structure.input
    def hierarchy(self):
        return self._structure.hierarchy

    def new_structure(self):
        """Generate a new copy of the input-hierarchy pair, from the pdb file"""
        return iotbx.pdb.hierarchy.input(file_name=self._pdb_file)

    def reflection_data(self):
        """Return an object containing the reflection data"""
        return iotbx.mtz.object(self._mtz_file)

    def get_heavy_atom_sites(self):
        xray_structure = self.input().xray_structure_simple()
        return xray_structure.sites_cart().select(xray_structure.heavy_selection())
    def get_calpha_sites(self):
        xray_structure = self.input().xray_structure_simple()
        return xray_structure.sites_cart().select(xray_structure.backbone_selection(atom_names=['CA']))
    def get_backbone_sites(self):
        xray_structure = self.input().xray_structure_simple()
        return xray_structure.sites_cart().select(xray_structure.backbone_selection())

    def set_global_alignment(self, alignment):
        self._global_rt_transform = alignment
    def global_alignment_transform(self):
        return self._global_rt_transform

    def set_local_alignments(self, alignment):
        self._local_rt_transforms = alignment
    def local_alignment_transforms(self):
        return self._local_rt_transforms

    def get_pickle_copy(self):
        """Get copy of self that can be pickled - some cctbx objects cannot be pickled..."""
        return self

    def get_structure_summary(self):
        return structure_summary(hierarchy=self.hierarchy())

    def transform_from_reference(self, points, method, point_mappings=None):
        """Use alignment to map from reference frame to our frame"""
        assert method in ['local','global'], 'METHOD NOT DEFINED: {!s}'.format(method)
        if method == 'global':
            # Simple - use one rt to transform all of the points
            return self.global_alignment_transform().inverse() * points
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
                lab_rt_points = self.local_alignment_transforms()[r_lab].inverse() * flex.vec3_double(lab_points)
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
            return self.global_alignment_transform() * points
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
                lab_rt_points = self.local_alignment_transforms()[r_lab] * flex.vec3_double(lab_points)
                # Populate the array at the appropriate place
                rt_points.put(lab_idxs, lab_rt_points)
                all_idxs.extend(lab_idxs)

            assert len(list(set(all_idxs))) == len(points)

            # Initial brute force method
            #rt_points = [self.local_alignment_transforms()[r_lab] * p for r_lab, p in zip(point_mappings, points)]

            return flex.vec3_double(rt_points)

    def find_nearest_calpha(self, points):
        """Returns the labels of the nearest calpha for each of the given points"""

        calpha_hierarchy = self.hierarchy().select(self.hierarchy().atom_selection_cache().selection('pepnames and name CA'))
        atom_sites, atom_labels = zip(*[(a.xyz, (a.chain_id, a.resid())) for a in calpha_hierarchy.atoms_with_labels()])

        tree = scipy.spatial.KDTree(data=atom_sites)
        nn_dists, nn_groups = tree.query(points)
        return [atom_labels[i] for i in nn_groups]

class ReferenceDatasetHandler(DatasetHandler):
    _origin_shift = (0,0,0)
    _binning = None

    def set_origin_shift(self, origin_shift):
        self._origin_shift = origin_shift
    def origin_shift(self):
        return self._origin_shift

    def set_map_scale(self, map_mean, map_rms):
        """Record the map mean and rms values"""
        self._map_mean = map_mean
        self._map_rms  = map_rms
    def map_scale(self):
        return (self._map_mean, self._map_rms)

class DatasetHandlerList(HolderList):
    """Class for grouping many dataset handlers together"""
    _holder_class = DatasetHandler

    def __custom_init__(self, *args):
        print(args)

class MapHolder(object):
    child = None
    parent = None
    """Class to hold map values and meta data"""
    def __init__(self, num, tag, map, unit_cell, space_group, meta, parent=None, child=None):
        assert isinstance(num, int), 'Num must be int. Type given: {!s}'.format(type(num))
        self.num = num
        assert isinstance(tag, str), 'Tag must be str. Type given: {!s}'.format(type(tag))
        self.tag = tag
        assert isinstance(map, flex.double), 'Map data must be flex.double. Type given: {!s}'.format(type(map))
        self.map = map
        assert isinstance(unit_cell, cctbx.uctbx.ext.unit_cell), 'Unit cell must be of type unit_cell. Type given: {!s}'.format(type(unit_cell))
        self.unit_cell = unit_cell
        assert isinstance(space_group, cctbx.sgtbx.ext.space_group), 'Space group must be of type space_group. Type given: {!s}'.format(type(space_group))
        self.space_group = space_group
        assert isinstance(meta, Meta), 'Meta must be dict. Type given: {!s}'.format(type(meta))
        self.meta = meta
        if parent:
            assert isinstance(parent, DatasetHandler) or isinstance(parent, MapHolder), 'parent must be of type DatasetHandler or MapHolder. Type given: {!s}'.format(type(parent))
            self.parent = parent
        if child:
            self.child = child

class MapHolderList(HolderList):
    """Class for grouping many MapHolder objects together"""

    def __custom_init__(self, *args):
        print(args)

def map_handler(map_data, unit_cell):
    """Map handler for easy sampling of map"""

    basic_map = cctbx.maptbx.basic_map( cctbx.maptbx.basic_map_unit_cell_flag(),
                                        map_data,
                                        map_data.focus(),
                                        unit_cell.orthogonalization_matrix(),
                                        cctbx.maptbx.out_of_bounds_clamp(0).as_handle(),
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

def align_dataset_to_reference(d_handler, ref_handler, method):
    """Calculate the rotation and translation needed to align one structure to another"""

    assert method in ['both','local','global'], 'METHOD NOT DEFINED: {!s}'.format(method)

    if method == 'global' or method == 'both':
        my_sites = d_handler.get_calpha_sites()
        ref_sites = ref_handler.get_calpha_sites()
        assert len(my_sites) == len(ref_sites)
        global_rt_transform = superpose.least_squares_fit(reference_sites=ref_sites, other_sites=my_sites).rt()
    else:
        global_rt_transform = None

    if method == 'local' or method == 'both':
        local_rt_transforms = perform_flexible_alignment(mov_hierarchy=d_handler.hierarchy(), ref_hierarchy=ref_handler.hierarchy())
    else:
        local_rt_transforms = None

    return global_rt_transform, local_rt_transforms

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
                print('\rAnalysing Chain {!s}: {!s}-->'.format(current_chain, amino_sequence), end=''); sys.stdout.flush()

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

class MapList(Info):
    _map_names = []
    def __init__(self, map_names=[]):
        assert self._map_names+map_names, 'No Map Names defined'
        for m in self._map_names+map_names:
            self.__dict__[m] = None
        # Intialise meta dict
        self.meta = Meta()
        # Set initialised so attributes can't be added
        self._initialized = True

class PanddaStatMapList(MapList):
    _map_names = ['mean_map','medn_map','stds_map','sadj_map','skew_map','kurt_map','bimo_map']

class PanddaMultipleStatMapList(object):
    def __init__(self):
        """Store lots of statistical maps"""
        self.map_lists = {}
    def get_resolutions(self):
        return sorted(self.map_lists.keys())
    def add(self, stat_map_list, resolution):
        assert isinstance(resolution, float), 'Resolution of map must be of type float. Type given: {!s}'.format(type(resolution))
        assert resolution not in self.map_lists.keys(), 'MAPS OF THIS RESOLUTION ALREADY ADDED'
        assert isinstance(stat_map_list, PanddaStatMapList), 'stat_map_list must be of type PanddaMultipleStatMapList. Type given: {!s}'.format(type(stat_map_list))
        self.map_lists[resolution] = stat_map_list
    def get(self, resolution):
        return self.map_lists[resolution]

class PanddaMapAnalyser(object):
    """Class to hold dataset maps, statistical maps and meta data for a set of related maps. Also holds functions for analysing the maps."""
    output_handler = None
    def __init__(self, dataset_maps, meta, statistical_maps=None, parent=None):
        # Validate the meta
        assert isinstance(meta, Meta), 'meta must be of type Meta. Type given: {!s}'.format(type(meta))
        assert hasattr(meta, 'resolution')
        assert hasattr(meta, 'grid_size')
        assert hasattr(meta, 'grid_size_1d')
        self.meta = meta
        # Validate the dataset maps
        assert isinstance(dataset_maps, MapHolderList), 'dataset_maps must be stored in a MapHolderList. Type given: {!s}'.format(type(dataset_maps))
        self.dataset_maps = dataset_maps
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
            self.log = self.parent.log
        else:
            self.parent = None
            self.log = print

    def validate_maps(self):
        """Check that all of the added maps are the same size etc..."""

        for mh in self.dataset_maps.all():
            print('Checking Map {!s}'.format(mh.name))

    def calculate_mean_map(self, masked_idxs=None, mask_name=None):
        """Calculate the mean map from all of the different observations"""

        if not masked_idxs: masked_idxs = range(0, self.meta.grid_size_1d)
        else:               assert max(masked_idxs) < self.meta.grid_size_1d, 'masked_idxs out of range of map'

        # Create statistics objects for each grid point
        self.log('===================================>>>', True)
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
            status_bar(n=i_chunk, n_max=num_chunks)
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

        status_bar(n=num_chunks, n_max=num_chunks)

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

        if output_graphs:
            try:
                import matplotlib
                matplotlib.interactive(0)
                from matplotlib import pyplot
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

        self.log('===================================>>>')
        t1 = time.time()

        for i_mh, mh in enumerate(self.dataset_maps.mask(mask_name=mask_name)):

            if mh.meta.map_uncertainty is not None:
                print('SKIPPING Dataset {!s} ({!s}/{!s})                                 '.format(mh.tag, i_mh+1, self.dataset_maps.size(mask_name=mask_name)))
                continue
            else:
                print('\rCalculating Map Uncertainty for Dataset {!s} ({!s}/{!s})          '.format(mh.tag, i_mh+1, self.dataset_maps.size(mask_name=mask_name)), end=''); sys.stdout.flush()

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

                graphs.mean_obs_scatter( f_name    = mh.parent.output_handler.get_file('obs_qqplot_unsorted_png'),
                                         mean_vals = masked_mean_vals,
                                         obs_vals  = masked_map_vals     )

                graphs.sorted_mean_obs_scatter( f_name    = mh.parent.output_handler.get_file('obs_qqplot_sorted_png'),
                                                mean_vals = sorted_mean_vals,
                                                obs_vals  = sorted_observed_vals    )

                graphs.diff_mean_qqplot( f_name   = mh.parent.output_handler.get_file('unc_qqplot_png'),
                                         map_off  = map_off,
                                         map_unc  = map_unc,
                                         q_cut    = q_cut,
                                         obs_diff = sorted_observed_diff_vals,
                                         quantile = expected_diff_vals  )

        self.log('\rMap Uncertainties Calculated.                                         ', True)

        t2 = time.time()
        self.log('> MAP UNCERTAINTY CALCULATION > Time Taken: {!s} seconds'.format(int(t2-t1)))

        return [mh.meta.map_uncertainty for mh in self.dataset_maps.mask(mask_name=mask_name)]

    def calculate_statistical_maps(self, masked_idxs=None, mask_name=None, cpus=1, ignore_warnings=True):
        """Take the sampled maps and calculate statistics for each grid point across the datasets"""

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

        # Create statistics objects for each grid point
        self.log('===================================>>>', True)
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
            status_bar(n=i_chunk, n_max=num_chunks)
            # Select map values from each map
            p_map_vals = numpy.array([mh.map.select(chunk_idxs) for mh in self.dataset_maps.mask(mask_name=mask_name)])
            if i_chunk+1 < num_chunks:
                assert len(p_map_vals) == self.dataset_maps.size(mask_name=mask_name)
                assert len(p_map_vals.T) == chunk_size
            if cpus > 1:
                # Generate arg_list for analysis
                arg_list = [{'map_vals':map_vals, 'map_uncertainties':all_uncertainties} for map_vals in p_map_vals.T]
                point_statistics = easy_mp.pool_map(fixed_func=map_statistics_map_func, args=arg_list, processes=cpus)
            else:
                point_statistics = [map_statistics_map_func({'map_vals':map_vals, 'map_uncertainties':all_uncertainties}) for map_vals in p_map_vals.T]

            masked_point_statistics.extend(point_statistics)

        status_bar(n=num_chunks, n_max=num_chunks)
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

class PanddaInfoTables(object):

    _all_dataset_fields     = DATASET_INFO_FIELDS
    _all_dataset_map_fields = DATASET_MAP_FIELDS
    _all_event_fields       = DATASET_EVENT_FIELDS
    _all_site_fields        = SITE_TABLE_FIELDS

class PanddaMultiDatasetAnalyser(object):
    """Class for the processing of datasets from a fragment soaking campaign"""

    _version = PANDDA_VERSION

    _structure_mask_names = STRUCTURE_MASK_NAMES
    _crystal_mask_names   = CRYSTAL_MASK_NAMES
    _reject_mask_names    = REJECT_MASK_NAMES
    _flag_mask_names      = FLAG_MASK_NAMES
    _custom_mask_names    = ['no_analyse', 'no_build']
    _all_mask_names       = _structure_mask_names + _crystal_mask_names + _reject_mask_names + _flag_mask_names + _custom_mask_names

    _all_dataset_fields     = DATASET_INFO_FIELDS
    _all_dataset_map_fields = DATASET_MAP_FIELDS
    _all_event_fields       = DATASET_EVENT_FIELDS
    _all_site_fields        = SITE_TABLE_FIELDS

    def __init__(self, args=None):
#                    outdir='./pandda', datadir='./Processing', pdb_style='*/refine.pdb', mtz_style='*/refine.mtz',
#                    ref_pdb='./reference.pdb', ref_mtz='./reference.mtz', run_mol_subst=False,
#                    verbose=True, keep_maps_in_memory=False, maxmemory=25):
        """Class for the processing of datasets from a fragment soaking campaign"""

        # Allow the program to pull from the command line if no arguments are given
        if args == None:
            args = sys.argv[1:]

        # ===============================================================================>
        # PROCESS INPUT ARGUMENTS
        # ===============================================================================>

        # Read in the master phil
        self.master_phil = phil.parse(pandda_phil_def)

        # Show defaults and exit
        if '--show-defaults' in args:
            self.master_phil.show()
            sys.exit()

        # Process the input arguments and convert to phil
        self.cmds = args
        self.working_phil, self.unused_args = self.process_input_args(master_phil=self.master_phil, args=args)
        # Pull out the python object of the arguments (at the `pandda` level)
        self.args = self.working_phil.extract().pandda
        # Most of the variables are contained withing the params object - create a shortcut to save typing
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

        self._log = Log(log_file=os.path.join(self.outdir, 'pandda.log'))

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
                                                            columns = self._all_dataset_fields      )
        # Record information about the created maps for each dataset
        self.tables.dataset_map_info    = pandas.DataFrame( data    = None,
                                                            index   = pandas.Index(data=[], name='dtag'),
                                                            columns = self._all_dataset_map_fields  )
        # Record the events detected in each dataset
        self.tables.event_info          = pandas.DataFrame( data    = None,
                                                            index   = pandas.MultiIndex(levels=[[],[]], labels=[[],[]], names=['dtag','event_idx']),
                                                            columns = self._all_event_fields        )
        # Record information about the clustered events (cluster of events = site)
        self.tables.site_info           = pandas.DataFrame( data    = None,
                                                            index   = pandas.Index(data=[], name='site_idx'),
                                                            columns = self._all_site_fields         )

#        # TODO --- XXX --- TODO --- XXX --- TODO --- XXX --- TODO
#        # Map and Dataset statistics - TODO DELETE TODO
#        self.map_observations = data_collection()
#        # TODO --- XXX --- TODO --- XXX --- TODO --- XXX --- TODO
#
#        # TODO --- XXX --- TODO --- XXX --- TODO --- XXX --- TODO
#        # Summaries of detected points - TODO DELETE TODO
#        self._cluster_summary = data_collection()
#        # TODO --- XXX --- TODO --- XXX --- TODO --- XXX --- TODO

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

        if os.path.exists(self.pickle_handler.get_file('dataset_meta')):
            self._new_pandda = False
        else:
            self._new_pandda = True

    def process_input_args(self, master_phil, args):
        """Process the input arguments"""

        # Copy the args so that we can remove items from the list without affecting sys.argv etc
        args = copy.copy(args)

        assert isinstance(args, list), 'INPUT ARGUMENTS MUST BE A LIST'
        assert len(args) != 0, 'NO INPUT ARGUMENTS GIVEN'

        # Build an interpreter from the master phil
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
        working_phil = master_phil.fetch(sources=eff_sources+cmd_sources)

        return working_phil, args

    def run_pandda_init(self):
        """Set up the pandda"""

        # Update Status
        self.update_status('running')

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
                                        '# ' + ' '.join(self.cmds),
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
        # PRINT SOME HELPFUL INFORMATION
        # ===============================================================================>

        self.log('===================================>>>', True)
        self.log('RUNNING FROM: {!s}'.format(sys.argv[0]), True)
        self.log('===================================>>>', True)
        self.log('READING INPUT FROM : {!s}'.format(self.data_dirs), True)
        self.log('===================================>>>', True)
        self.log('WRITING OUTPUT TO: {!s}'.format(self.outdir), True)

        # ===============================================================================>
        # REPOPULATE PANDDA FROM PREVIOUS RUNS
        # ===============================================================================>

        if os.path.exists(self.output_handler.get_file('reference_structure')) and os.path.exists(self.output_handler.get_file('reference_dataset')):
            self.log('===================================>>>', True)
            self.log('Loading Reference Dataset', True)
            self.load_reference_dataset(ref_pdb=self.output_handler.get_file('reference_structure'), ref_mtz=self.output_handler.get_file('reference_dataset'))

        # Load any objects from previous runs
        self.load_pickled_objects()

        # ===============================================================================>
        # LOOK FOR MATPLOTLIB TO SEE IF WE CAN GENERATE GRAPHS
        # ===============================================================================>

        if self.args.output.plot_graphs:
            try:
                import matplotlib
                # Setup so that we can write without a display connected
                matplotlib.interactive(0)
                default_backend, validate_function = matplotlib.defaultParams['backend']
                self.log('===================================>>>')
                self.log('MATPLOTLIB LOADED. Using Backend: {!s}'.format(default_backend))
                from matplotlib import pyplot
                assert not pyplot.isinteractive(), 'PYPLOT IS STILL IN INTERACTIVE MODE'
                self.log('PYPLOT loaded successfully')
            except:
                self.log('===================================>>>', True)
                self.log('>> COULD NOT IMPORT MATPLOTLIB. CANNOT GENERATE GRAPHS.', True)
                self.args.output.plot_graphs = False

        # ===============================================================================>
        # LOG THE START TIME
        # ===============================================================================>

        # Store start time and print to log
        self._init_time = time.time()
        self.log('===================================>>>', True)
        self.log('Analysis Started: {!s}'.format(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime(self._init_time))), True)

        # Get the size of the empty pandda
        self.log('===================================>>>', True)
        self.update_pandda_size(tag='Initialised Pandda')

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

        # ===============================================================================>
        # Global Directories that do not change from run to run
        # ===============================================================================>

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
        # Somewhere to record which dataset was processed at which resolution
        self.output_handler.add_dir(dir_name='resolutions', dir_tag='resolutions', top_dir_tag='root')
        # Somewhere to store all of the aligned structures
        self.output_handler.add_dir(dir_name='aligned_structures', dir_tag='aligned_structures', top_dir_tag='root')

        # ===============================================================================>
        # New directories that will be created for each run (so that data is not overwritten)
        # ===============================================================================>

#        # Create a label for this run of the program
#        run_label = '{!s}'.format(time.strftime("%Y-%m-%d-%H:%M", time.localtime()))
#        # Top folder for this run
#        self.output_handler.add_dir(dir_name='ANALYSIS-{!s}'.format(run_label), dir_tag='run_directory', top_dir_tag='root')

        # Input parameters
        self.output_handler.add_file(file_name='pandda.eff',     file_tag='pandda_settings', dir_tag='root')
        self.output_handler.add_file(file_name='pandda.running', file_tag='pandda_status_running', dir_tag='root')
        self.output_handler.add_file(file_name='pandda.done',    file_tag='pandda_status_done', dir_tag='root')

        # Somewhere to store the analyses/summaries - for me to plot graphs
        self.output_handler.add_dir(dir_name='analyses', dir_tag='analyses', top_dir_tag='root')
        # Somewhere to store the dataset information (general values)
        self.output_handler.add_file(file_name='dataset_info.csv', file_tag='dataset_info',  dir_tag='analyses')
        self.output_handler.add_file(file_name='dataset_map_info.csv', file_tag='dataset_map_info',  dir_tag='analyses')
        self.output_handler.add_file(file_name='dataset_combined_info.csv', file_tag='dataset_combined_info',  dir_tag='analyses')
        # Somewhere to store the dataset information (identified events + sites)
        self.output_handler.add_file(file_name='event_info.csv', file_tag='event_info',  dir_tag='analyses')
        self.output_handler.add_file(file_name='event_info_2.csv', file_tag='event_info_2',  dir_tag='analyses')
        self.output_handler.add_file(file_name='site_info.csv', file_tag='site_info',  dir_tag='analyses')
        self.output_handler.add_file(file_name='point_distributions.csv', file_tag='point_distributions', dir_tag='analyses')
        # Somewhere to store the analysis summaries - for the user
        self.output_handler.add_dir(dir_name='results_summaries', dir_tag='output_summaries', top_dir_tag='root')
        self.output_handler.add_file(file_name='success_summary_table.html', file_tag='summary_table', dir_tag='output_summaries')
        self.output_handler.add_file(file_name='site_bar_graph.png', file_tag='site_bar_graph', dir_tag='output_summaries')
        # Somewhere to store the pickled objects
        self.output_handler.add_dir(dir_name='pickled_panddas', dir_tag='pickle', top_dir_tag='root')

        # ===============================================================================>
        # Reference Structure Files (should only be needed once for writing and then only for reloading)
        # ===============================================================================>

        # Reference Structure and Dataset
        self.output_handler.add_dir(dir_name='reference', dir_tag='reference', top_dir_tag='root')
        self.output_handler.add_file(file_name='reference.pdb', file_tag='reference_structure', dir_tag='reference')
        self.output_handler.add_file(file_name='reference.mtz', file_tag='reference_dataset', dir_tag='reference')
        self.output_handler.add_file(file_name='reference.shifted.pdb', file_tag='reference_on_origin', dir_tag='reference')
        self.output_handler.add_file(file_name='reference.symmetry.pdb', file_tag='reference_symmetry', dir_tag='reference')

        # ===============================================================================>
        # Standard template files that will be populated when needed
        # ===============================================================================>

        # Statistical Maps
        self.output_handler.add_dir(dir_name='statistical_maps', dir_tag='statistical_maps', top_dir_tag='root')
        self.output_handler.add_file(file_name='{!s}A-mean_map.ccp4', file_tag='mean_map', dir_tag='statistical_maps')
        self.output_handler.add_file(file_name='{!s}A-medn_map.ccp4', file_tag='medn_map', dir_tag='statistical_maps')
        self.output_handler.add_file(file_name='{!s}A-stds_map.ccp4', file_tag='stds_map', dir_tag='statistical_maps')
        self.output_handler.add_file(file_name='{!s}A-sadj_map.ccp4', file_tag='sadj_map', dir_tag='statistical_maps')
        self.output_handler.add_file(file_name='{!s}A-skew_map.ccp4', file_tag='skew_map', dir_tag='statistical_maps')
        self.output_handler.add_file(file_name='{!s}A-kurt_map.ccp4', file_tag='kurt_map', dir_tag='statistical_maps')
        self.output_handler.add_file(file_name='{!s}A-bimo_map.ccp4', file_tag='bimo_map', dir_tag='statistical_maps')

    def _run_pickle_setup(self):
        """Initialise all of the pickle filenames"""

        # Pickle Handler
        self.pickle_handler = output_file_object(rootdir=self.output_handler.get_dir('pickle'))
        # Pickled Reference Objects
        self.pickle_handler.add_file(file_name='reference_grid.pickle', file_tag='reference_grid')
        self.pickle_handler.add_file(file_name='reference_dataset.pickle', file_tag='reference_dataset')
        # Pickled Datasets
        self.pickle_handler.add_file(file_name='dataset_masks.pickle', file_tag='dataset_masks')
        self.pickle_handler.add_file(file_name='dataset_meta.pickle', file_tag='dataset_meta')
        # Pickled Information
        self.pickle_handler.add_file(file_name='dataset_info.pickle', file_tag='dataset_info')
        self.pickle_handler.add_file(file_name='dataset_map_info.pickle', file_tag='map_info')
        # Pickled Stats
        self.pickle_handler.add_file(file_name='statistical_maps.pickle', file_tag='stat_maps')

        # Pickled SELF
        self.pickle_handler.add_file(file_name='my_pandda.pickle', file_tag='my_pandda')

    def load_pickled_objects(self):
        """Loads any pickled objects it finds"""

        self.log('===================================>>>', True)
        self.log('Looking for Pickled Files in Input Directory: {!s}'.format(os.path.relpath(self.pickle_handler.get_dir('root'))), True)

        # Record whether any pickled objects are loaded
        pickles_found = False

        # Load Reference Grid
        if os.path.exists(self.pickle_handler.get_file('reference_grid')):
            pickles_found = True
            self.log('===> Loading reference grid')
            self.set_reference_grid(self.unpickle(self.pickle_handler.get_file('reference_grid')))

        # Load Reference Dataset
        if os.path.exists(self.pickle_handler.get_file('reference_dataset')):
            pickles_found = True
            self.log('===> Loading reference dataset')
            self.set_reference_dataset(self.unpickle(self.pickle_handler.get_file('reference_dataset')))

        # Load the datasets
        if os.path.exists(self.pickle_handler.get_file('dataset_meta')):
            pickles_found = True
            # Unpickle the list of the pickled datasets from the directory structure
            self.log('===> Loading old dataset information (existing datasets)')
            self.pickled_dataset_meta = self.unpickle(self.pickle_handler.get_file('dataset_meta'))

            if self.args.method.reload_existing_datasets:
                # Extract the paths of the pickled dataset objects
                pickled_dataset_list = self.pickled_dataset_meta.dataset_pickle_list
                # Check they all exist - should be relative to the outdirectory
                for filename in pickled_dataset_list:
                    assert os.path.isfile(os.path.join(self.outdir, filename)), 'File does not exist: {!s}'.format(filename)
                # Unpickle the datasets and add them to the dataset handler list
                self.log('===> Reloading old datasets')
                self.datasets.add([self.unpickle(os.path.join(self.outdir,f)) for f in pickled_dataset_list])
                self.update_pandda_size(tag='After Unpickling Dataset Objects')
            else:
                self.log('===> Not reloading old datasets')
        else:
            # No datasets to load - this must be False
            self.args.method.reload_existing_datasets = False
            self.log('===> No old datasets found')

        # Load Statistical Maps
        if os.path.exists(self.pickle_handler.get_file('stat_maps')):
            pickles_found = True
            self.log('===> Loading old statistical maps')
            self.stat_maps = self.unpickle(self.pickle_handler.get_file('stat_maps'))

        if not pickles_found:
            self.log('===> No Pickles Found', True)

    def pickle_the_pandda(self, components=[], all=False, datasets=None):
        """Pickles it's major components for quick loading..."""

        if all == True:
            self.log('===================================>>>', True)
            self.log('Pickling the Pandda', True)
        elif not components:
            self.log('===================================>>>', True)
            self.log('Pickling NOTHING', True)
            return
        else:
            self.log('===================================>>>', True)
            self.log('Selective Pickling: {!s}'.format(', '.join(components).upper()), True)

        if all or ('grid' in components):
            self.log('===================================>>>')
            if self.reference_grid() is not None:
                self.log('Pickling Reference Grid')
                self.pickle(pickle_file   = self.pickle_handler.get_file('reference_grid'),
                            pickle_object = self.reference_grid(),
                            overwrite = False)
            else:
                self.log('No Reference Grid to Pickle')

        if all or ('datasets' in components):
            self.log('===================================>>>')

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
            self.log('===================================>>>')
            if self.stat_maps is not None:
                self.log('Pickling Statistical Maps')
                self.pickle(pickle_file   = self.pickle_handler.get_file('stat_maps'),
                            pickle_object = self.stat_maps,
                            overwrite = True)
            else:
                self.log('No Statistical Maps to Pickle')

    def update_status(self, status):
        """Set log files to indicate the status of the program"""

        assert status in ['running','done']
        running_file = self.output_handler.get_file('pandda_status_running')
        done_file =    self.output_handler.get_file('pandda_status_done')

        if status == 'running':
            if os.path.exists(done_file):        os.remove(done_file)
            if not os.path.exists(running_file): os.mknod(running_file)
        elif status == 'done':
            if os.path.exists(running_file):  os.remove(running_file)
            if not os.path.exists(done_file): os.mknod(done_file)

    def exit(self, error=False):
        """Exit the PANDDA, record runtime etc..."""
        self.log('===================================>>>', True)
        self.log('...FINISHED!...', True)
        self.log('===================================>>>', True)
        self._finish_time = time.time()
        self.log('Runtime: {!s}'.format(time.strftime("%H hours:%M minutes:%S seconds", time.gmtime(self._finish_time - self._init_time))))

        try:
            if self.args.settings.pickling.full_pickle:
                # Pickle myself
                self.log('===================================>>>', True)
                self.log('Pickling the PANDDA Results')
                self.pickle(pickle_file=self.pickle_handler.get_file('my_pandda'), pickle_object=self, overwrite=True)
        except:
            self.log('FAILED TO PICKLE MYSELF')

        try:
            # Extract meta about the datasets
            if self.pickled_dataset_meta and (not self.args.method.reload_existing_datasets):
                # Combine the old meta with the new
                self.log('Combining old dataset meta with new meta for pickle')
                number_of_datasets  = self.pickled_dataset_meta.number_of_datasets  + self.datasets.size()
                dataset_labels      = self.pickled_dataset_meta.dataset_labels      + [d.tag for d in self.datasets.all()]
                dataset_pickle_list = self.pickled_dataset_meta.dataset_pickle_list + [os.path.relpath(d.output_handler.get_file('dataset_pickle'), start=self.outdir) for d in self.datasets.all()]
            else:
                # Don't need to consider the pickle as the datasets have been reloaded
                self.log('Creating new meta for pickle')
                number_of_datasets  = self.datasets.size()
                dataset_labels      = [d.tag for d in self.datasets.all()]
                dataset_pickle_list = [os.path.relpath(d.output_handler.get_file('dataset_pickle'), start=self.outdir) for d in self.datasets.all()]
            # Create a dictionary to be stored
            dataset_meta = Meta({'number_of_datasets'    : number_of_datasets,
                                 'dataset_labels'        : dataset_labels,
                                 'dataset_pickle_list'   : dataset_pickle_list
                            })
            # Pickle the list of locations of the dataset pickles
            self.pickle(pickle_file   = self.pickle_handler.get_file('dataset_meta'),
                        pickle_object = dataset_meta,
                        overwrite = True)
        except:
            self.log('FAILED TO PICKLE META')

    def log(self, message, show=False, hide=False):
        """Log message to file, and mirror to stdout if verbose or force_print (hide overrules show)"""
        if not isinstance(message, str):    message = str(message)
        # Print to stdout
        if (not hide) and (show or self.args.settings.verbose):
            self._log.show(message)
        # Remove \r from message as this spoils the log (^Ms)
        message = message.replace('\r','')
        # Write to file
        self._log.write(message=message+'\n', mode='a')

    def print_log(self):
        self._log.show(self._log.read_all())

    def update_pandda_size(self, tag):
        pandda_size = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1000
        assert pandda_size < self.args.settings.max_memory*1024**3, 'PANDDA HAS EXCEEDED THE MAXIMUM AMOUNT OF ALLOWED MEMORY'
        self._pandda_size.append((tag, pandda_size))
        try:
            from humanize import naturalsize
            self.log(tag+': '+naturalsize(pandda_size, binary=True))
        except:
            pass

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
        self.datasets.all_masks().set_entry_ids(entry_ids=[d.tag for d in self.datasets.all()])
        # Initialise standard blank masks
        for mask_name in self._all_mask_names:
            self.datasets.all_masks().add_mask(mask_name=mask_name, mask=[False]*self.datasets.size())

        # Initialise masks for datasets that shouldn't be analysed
        if self.params.analysis.no_analyse:
            no_analyse_tags = self.params.analysis.no_analyse.split(',')
            self.log('Not analysing {!s} Datasets: {!s}'.format(len(no_analyse_tags), ', '.join(no_analyse_tags)))
            no_analyse_mask = [True if d.tag in no_analyse_tags else False for d in self.datasets.all()]
            self.datasets.all_masks().add_mask(mask_name='no_analyse', mask=no_analyse_mask)

        # Initialise mask for datasets that shouldn't be used for building
        if self.params.analysis.no_build:
            no_build_tags = self.params.analysis.no_build.split(',')
            self.log('Not building distributions from {!s} Datasets: {!s}'.format(len(no_build_tags), ', '.join(no_build_tags)))
            no_build_mask = [True if d.tag in no_build_tags else False for d in self.datasets.all()]
            self.datasets.all_masks().add_mask(mask_name='no_build', mask=no_build_mask)

        # Initialise mask for datasets that have been previously pickled
        self.datasets.all_masks().add_mask(mask_name='old datasets', mask=[False]*self.datasets.size())

        if self.pickled_dataset_meta and self.args.method.reload_existing_datasets:
            for tag in self.pickled_dataset_meta.dataset_tags:
                self.datasets.all_masks().set_mask_value(mask_name='old datasets', entry_id=tag, value=True)
            self.log('Considering {!s} datasets as "New Datasets"'.format(self.datasets.size(mask_name='old datasets', invert=True)))
            self.log('Considering {!s} datasets as "Old Datasets"'.format(self.datasets.size(mask_name='old datasets')))
        else:
            self.log('Considering all {!s} datasets as "New Datasets"'.format(self.datasets.size(mask_name='old datasets', invert=True)))
            assert self.datasets.size(mask_name='old datasets', invert=True) == self.datasets.size(), 'Total datasets should be same as total new datasets'

        self.log('===================================>>>', True)
        self.log('Initialising Dataset Information Tables.', True)

        # Add dataset tags as rows in the tables
        self.tables.dataset_info     = self.tables.dataset_info.append(pandas.DataFrame(index=[d.tag for d in self.datasets.all()]), verify_integrity=True)
        self.tables.dataset_map_info = self.tables.dataset_map_info.append(pandas.DataFrame(index=[d.tag for d in self.datasets.all()]), verify_integrity=True)

        self.log('===================================>>>', True)
        self.log('Initialising Timings.', True)

        # Set the dataset observations lengths
        self.timings.set_data_length(self.datasets.size())
        self.timings.set_entry_ids([d.tag for d in self.datasets.all()])

    def load_reference_dataset(self, ref_pdb, ref_mtz):
        """Set the reference dataset, to which all other datasets will be aligned and scaled"""

        self.log('===================================>>>', True)
        self.log('Loading Reference Dataset: {!s}'.format(ref_mtz), True)

        self.set_reference_dataset(ReferenceDatasetHandler(dataset_number=-1, pdb_filename=ref_pdb, mtz_filename=ref_mtz, dataset_tag='reference'))

        if not os.path.exists(self.output_handler.get_file('reference_structure')):
            rel_symlink(orig=ref_pdb, link=self.output_handler.get_file('reference_structure'))
        if not os.path.exists(self.output_handler.get_file('reference_dataset')):
            rel_symlink(orig=ref_mtz, link=self.output_handler.get_file('reference_dataset'))

        # Calculate the shift required to move the reference structure into the positive quadrant
        total_border_padding = self.params.maps.padding + self.params.masks.outer_mask
        self.reference_dataset().set_origin_shift(tuple(flex.double(3, total_border_padding) - flex.double(self.reference_dataset().input().atoms().extract_xyz().min())))
        self.log('Origin Shift for reference structure: {!s}'.format(tuple([round(s,3) for s in self.reference_dataset().origin_shift()])))
        # Shift the reference structure by this amount so that it is aligned with the reference grid
        ref_hierarchy = self.reference_dataset().hierarchy()
        ref_hierarchy.atoms().set_xyz(ref_hierarchy.atoms().extract_xyz() + self.reference_dataset().origin_shift())

        if not os.path.exists(self.output_handler.get_file('reference_on_origin')):
            ref_hierarchy.write_pdb_file(self.output_handler.get_file('reference_on_origin'))

        # Align reference dataset to it's shifted self
        r = scitbx.matrix.rec([1,0,0,0,1,0,0,0,1], (3,3))
        t = scitbx.matrix.rec(self.reference_dataset().origin_shift(), (3,1))
        rt = scitbx.matrix.rt((r,t))
        self.reference_dataset().set_global_alignment(alignment=rt)

        return self.reference_dataset()

    def create_reference_grid(self, grid_spacing, expand_to_origin, buffer=0):
        """Create a grid over the reference protein"""

        self.log('===================================>>>', True)
        self.log('Creating Reference Grid', True)

        sites_cart = self.reference_dataset().input().atoms().extract_xyz()

        assert (flex.vec3_double(sorted(self.reference_dataset().input().atoms().extract_xyz())) -
                flex.vec3_double(sorted(self.reference_dataset().hierarchy().atoms().extract_xyz()))).dot().norm() == 0.0, 'EH? Coordinates should be the same?'

        self._ref_grid = grid_handler(verbose=self.args.settings.verbose)
        self._ref_grid.set_grid_spacing(spacing=grid_spacing)
        self._ref_grid.set_cart_extent(cart_min=tuple([s-buffer for s in sites_cart.min()]), cart_max=tuple([s+buffer for s in sites_cart.max()]))
        self._ref_grid.create_cartesian_grid(expand_to_origin=expand_to_origin)
        self.log(self._ref_grid.summary())

        if self.params.alignment.method == 'local':
            self.log('===================================>>>', True)
            self.log('Partitioning Reference Grid', True)
            # Pull out the calphas
            calpha_hierarchy = self.reference_dataset().hierarchy().select(self.reference_dataset().hierarchy().atom_selection_cache().selection('pepnames and name CA'))
            t1 = time.time()
            # Calculate the nearest residue for each point on the grid
            self.reference_grid().create_grid_partition(atomic_hierarchy=calpha_hierarchy)
            # Partition Grid
            self.reference_grid().partition().partition_grid(cpus=self.args.settings.cpus)
            t2 = time.time()
            self.log('> GRID PARTITIONING > Time Taken: {!s} seconds'.format(int(t2-t1)))

    def mask_reference_grid(self):
        """Create masks for the reference grid based on distances from atoms in the reference structure"""

        self.log('===================================>>>', True)
        self.log('Masking Reference Grid', True)

        # ============================================================================>
        # Create neighbouring symmetry copies of the reference structures
        # ============================================================================>
        ref_sym_copies = self.generate_symmetry_copies(d_handler=self.reference_dataset(), save_operators=True)
        # Write out the symmetry sites
        if not os.path.exists(self.output_handler.get_file('reference_symmetry')):
            ref_sym_copies.write_pdb_file(self.output_handler.get_file('reference_symmetry'))
        # ============================================================================>
        # Local mask used for forming groups of points around a grid point
        # ============================================================================>
        if self.reference_grid().local_mask() is None:
            self.log('===================================>>>')
            self.log('Generating Local Mask')
            local_mask = spherical_mask(grid_spacing    = self.reference_grid().grid_spacing(),
                                        distance_cutoff = 1.2,
                                        grid_jump       = 1 )
            self.reference_grid().set_local_mask(local_mask)
        # ============================================================================>
        # Global mask used for removing points in the bulk solvent regions
        # ============================================================================>
        if self.reference_grid().global_mask() is None:
            self.log('===================================>>>')
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
            self.log('===================================>>>')
            self.log('Generating Symmetry Mask')
            # Pull out the cartesian sites of the symmetry mates
            cache = ref_sym_copies.atom_selection_cache()
            sym_sites_cart = ref_sym_copies.select(cache.selection('pepnames and not element H')).atoms().extract_xyz()
            # Generate the symmetry mask
            symmetry_mask = grid_mask(cart_sites = sym_sites_cart,
                                      grid_size  = self.reference_grid().grid_size(),
                                      unit_cell  = self.reference_grid().unit_cell(),
                                      max_dist   = self.params.masks.outer_mask,
                                      min_dist   = self.params.masks.inner_mask )
            self.reference_grid().set_symmetry_mask(symmetry_mask)

        # ============================================================================>
        # Write masked maps
        # ============================================================================>
        # Write protein masked map
        mask_map_file = self.output_handler.get_file('reference_dataset').replace('.mtz','.totalmask.ccp4')
        map_mask = flex.double(self.reference_grid().global_mask().total_mask_binary())
        map_mask.reshape(flex.grid(self.reference_grid().grid_size()))
        self.write_array_to_map(output_file=mask_map_file, map_data=map_mask)
        # Write symmetry masked map
        mask_map_file = self.output_handler.get_file('reference_dataset').replace('.mtz','.symmask.ccp4')
        map_mask = flex.double(self.reference_grid().symmetry_mask().total_mask_binary())
        map_mask.reshape(flex.grid(self.reference_grid().grid_size()))
        self.write_array_to_map(output_file=mask_map_file, map_data=map_mask)

        return

    def generate_symmetry_copies(self, d_handler, save_operators=False):
        """Generate the symmetry copies of the reference structure in the reference frame"""

        # Use symmetry operations to create the symmetry mates of the reference structure
        sym_ops, sym_hierarchies, chain_mappings = generate_adjacent_symmetry_copies(   ref_hierarchy    = d_handler.new_structure().hierarchy,
                                                                                        crystal_symmetry = d_handler.input().crystal_symmetry(),
                                                                                        buffer_thickness = self.params.masks.outer_mask+5    )
        # Record the symmetry operations that generate the crystal contacts
        if save_operators: self.crystal_contact_generators = sym_ops

        # Create a combined hierarchy of the crystal contacts
        symmetry_root = combine_hierarchies(sym_hierarchies)
        symmetry_root.atoms().set_xyz(d_handler.transform_to_reference(points=symmetry_root.atoms().extract_xyz(), method='global'))

        return symmetry_root

    def build_input_list(self):
        """Builds a list of input files from the command line arguments passed"""

        self.log('===================================>>>', True)
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

        for dir in sorted(glob.glob(self.data_dirs)):
            pdb_files = glob.glob(os.path.join(dir, pdb_style))
            mtz_files = glob.glob(os.path.join(dir, mtz_style))
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
                    pdb_regex = pdb_style.replace('*', '(.*)')
                    pdb_tag = re.findall(pdb_regex, pdb_base)
                    assert pdb_tag, 'NO PDB TAG FOUND: {!s} -> {!s}'.format(pdb_regex, pdb_base)
                    if isinstance(pdb_tag[0], tuple):
                        self.log('More than one PDB TAG found - choosing the first one of {!s}'.format(pdb_tag[0]))
                        pdb_tag = list(pdb_tag[0])[0:1]
                else: pdb_regex = pdb_tag = None

                if '*' in mtz_style:
                    mtz_base = os.path.basename(new_mtz)
                    mtz_regex = mtz_style.replace('*', '(.*)')
                    mtz_tag = re.findall(mtz_regex, mtz_base)
                    assert mtz_tag, 'NO MTZ TAG FOUND: {!s} -> {!s}'.format(mtz_regex, mtz_base)
                    if isinstance(mtz_tag[0], tuple):
                        self.log('More than one MTZ TAG found - choosing the first one of {!s}'.format(mtz_tag[0]))
                        mtz_tag = list(mtz_tag[0])[0:1]
                else: mtz_regex = mtz_tag = None

                if '*' in dir_style:
                    dir_base = os.path.dirname(pdb_files[0])
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
        if self.params.analysis.ignore_datasets:
            ignore_tags = self.params.analysis.ignore_datasets.split(',')
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
        self.log('===================================>>>', True)
        self.log('{!s} EMPTY DIRECTORIES FOUND (TOTAL)'.format(emp_num_i+emp_num_offset), True)
        self.log('{!s} EMPTY DIRECTORIES FOUND (NEW)'.format(emp_num_i), True)

        # Record total number of datasets, and total number of new datasets
        if self.pickled_dataset_meta: num_old = self.pickled_dataset_meta.number_of_datasets
        else:                         num_old = 0
        self.log('===================================>>>', True)
        self.log('{!s} DATASETS FOUND (TOTAL)'.format(len(filtered_new_files)+num_old), True)
        self.log('{!s} DATASETS FOUND (NEW)'.format(len(filtered_new_files)), True)

        return filtered_new_files

    def add_new_files(self, input_files):
        """Add (pdb, mtz) file pairs to the datasets to be processed"""

        self.log('===================================>>>', True)
        self._input_files = input_files

        self.log('{!s} Datasets Added'.format(len(input_files)), True)

    def load_new_datasets(self):
        """Read in maps for the input datasets"""

        def load_dataset_map_func(arg_dict):
            return DatasetHandler(**arg_dict)

        if not self.datasets.all() and self.is_new_pandda():
            self.log('Adding First Datasets to Pandda')
        else:
            self.log('Adding more datasets')
            self.log('{!s} already loaded'.format(self.datasets.size()))
            self.log('{!s} not loaded'.format(self.pickled_dataset_meta.number_of_datasets - self.datasets.size()))
            self.log('tags: X -> X')
            self.log('nums: Y -> Y')

        # Counting offset for dataset index
        if self.pickled_dataset_meta: n_offset = self.pickled_dataset_meta.number_of_datasets
        else:                         n_offset = 0

        # Generate arg_list for loading
        arg_list = [{'dataset_number':dnum+n_offset, 'pdb_filename':pdb, 'mtz_filename':mtz, 'dataset_tag':dtag, 'name_prefix':self.args.output.dataset_prefix} for dnum, (pdb, mtz, dtag) in enumerate(self.new_files())]

        start = time.time()
        self.log('===================================>>>', True)
        print('Adding Datasets... (using {!s} cores)'.format(self.args.settings.cpus))
        loaded_datasets = easy_mp.pool_map(fixed_func=load_dataset_map_func, args=arg_list, processes=self.args.settings.cpus)
        finish = time.time()
        self.log('> Adding Datasets > Time Taken: {!s} seconds'.format(int(finish-start)), True)

        lig_style = self.args.input.lig_style.strip('/')

        for d_handler in loaded_datasets:

            # Create a file manager object
            d_handler.initialise_output_directory(outputdir=os.path.join(self.output_handler.get_dir('processed_datasets'), d_handler.name))

            # Main input/output files
            d_handler.output_handler.add_file(file_name='{!s}-input.pdb'.format(d_handler.tag), file_tag='input_structure')
            d_handler.output_handler.add_file(file_name='{!s}-input.mtz'.format(d_handler.tag), file_tag='input_data')
            d_handler.output_handler.add_file(file_name='{!s}-info.csv'.format(d_handler.tag), file_tag='dataset_info', dir_tag='root')
            d_handler.output_handler.add_file(file_name='{!s}-aligned.pdb'.format(d_handler.tag), file_tag='aligned_structure')
            d_handler.output_handler.add_file(file_name='{!s}-sym_contacts.pdb'.format(d_handler.tag), file_tag='symmetry_copies')
            d_handler.output_handler.add_file(file_name='{!s}-observed.ccp4'.format(d_handler.tag), file_tag='sampled_map')
            d_handler.output_handler.add_file(file_name='{!s}-mean_diff.ccp4'.format(d_handler.tag), file_tag='mean_diff_map')
            d_handler.output_handler.add_file(file_name='{!s}-z_map.ccp4'.format(d_handler.tag), file_tag='z_map')
            d_handler.output_handler.add_file(file_name='{!s}-z_map_naive.ccp4'.format(d_handler.tag), file_tag='z_map_naive')
            d_handler.output_handler.add_file(file_name='{!s}-z_map_naive_normalised.ccp4'.format(d_handler.tag), file_tag='z_map_naive_normalised')
            d_handler.output_handler.add_file(file_name='{!s}-z_map_uncertainty.ccp4'.format(d_handler.tag), file_tag='z_map_uncertainty')
            d_handler.output_handler.add_file(file_name='{!s}-z_map_uncertainty_normalised.ccp4'.format(d_handler.tag), file_tag='z_map_uncertainty_normalised')
            d_handler.output_handler.add_file(file_name='{!s}-z_map_adjusted.ccp4'.format(d_handler.tag), file_tag='z_map_corrected')
            d_handler.output_handler.add_file(file_name='{!s}-z_map_adjusted_normalised.ccp4'.format(d_handler.tag), file_tag='z_map_corrected_normalised')
            d_handler.output_handler.add_file(file_name='{!s}-event_{!s}_occupancy_{!s}_map.ccp4'.format(d_handler.tag, '{!s}', '{!s}'), file_tag='occupancy_map')

            # Miscellaneous files
            d_handler.output_handler.add_file(file_name='{!s}-high_z_mask.ccp4'.format(d_handler.tag), file_tag='high_z_mask')
            d_handler.output_handler.add_file(file_name='{!s}-masked_grid.ccp4'.format(d_handler.tag), file_tag='grid_mask')

            # Links to ligand files (if they've been found)
            d_handler.output_handler.add_dir(dir_name='ligand_files', dir_tag='ligand', top_dir_tag='root')
            d_handler.output_handler.add_file(file_name='{!s}-ligand.pdb'.format(d_handler.tag), file_tag='ligand_coordinates')
            d_handler.output_handler.add_file(file_name='{!s}-ligand.cif'.format(d_handler.tag), file_tag='ligand_restraints')
            d_handler.output_handler.add_file(file_name='{!s}-ligand.png'.format(d_handler.tag), file_tag='ligand_image')

            # Native (back-rotated/transformed) maps
            d_handler.output_handler.add_file(file_name='{!s}-observed.native.ccp4'.format(d_handler.tag), file_tag='native_obs_map')
            d_handler.output_handler.add_file(file_name='{!s}-z_map.native.ccp4'.format(d_handler.tag), file_tag='native_z_map')
            d_handler.output_handler.add_file(file_name='{!s}-event_{!s}_occupancy_{!s}_map.native.ccp4'.format(d_handler.tag,'{!s}','{!s}'), file_tag='native_occupancy_map')

            # Fitted structures when modelled with pandda.inspect
            d_handler.output_handler.add_dir(dir_name='modelled_structures', dir_tag='models', top_dir_tag='root')

            # Output images
            d_handler.output_handler.add_dir(dir_name='output_images', dir_tag='images', top_dir_tag='root')
            # Smapled map
            d_handler.output_handler.add_file(file_name='{!s}-obsv_map_dist.png'.format(d_handler.tag), file_tag='s_map_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name='{!s}-diff_mean_map_dist.png'.format(d_handler.tag), file_tag='d_mean_map_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name='{!s}-z_map_dist_naive.png'.format(d_handler.tag), file_tag='z_map_naive_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name='{!s}-z_map_dist_naive_normalised.png'.format(d_handler.tag), file_tag='z_map_naive_normalised_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name='{!s}-z_map_dist_uncertainty.png'.format(d_handler.tag), file_tag='z_map_uncertainty_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name='{!s}-z_map_dist_uncertainty_normalised.png'.format(d_handler.tag), file_tag='z_map_uncertainty_normalised_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name='{!s}-z_map_dist_adjusted.png'.format(d_handler.tag), file_tag='z_map_corrected_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name='{!s}-z_map_dist_adjusted_normalised.png'.format(d_handler.tag), file_tag='z_map_corrected_normalised_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name='{!s}-z_map_dist_qq_plot.png'.format(d_handler.tag), file_tag='z_map_qq_plot_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name='{!s}-event_{!s}_occupancy_correlation.png'.format(d_handler.tag, '{!s}'), file_tag='occ_corr_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name='{!s}-uncertainty-qqplot.png'.format(d_handler.tag), file_tag='unc_qqplot_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name='{!s}-mean-v-obs-sorted-qqplot.png'.format(d_handler.tag), file_tag='obs_qqplot_sorted_png', dir_tag='images')
            d_handler.output_handler.add_file(file_name='{!s}-mean-v-obs-unsorted-plot.png'.format(d_handler.tag), file_tag='obs_qqplot_unsorted_png', dir_tag='images')

            # Edstats Scores
            d_handler.output_handler.add_file(file_name='{!s}-edstats.scores'.format(d_handler.tag), file_tag='edstats_scores')

            # Analysis files
            d_handler.output_handler.add_file(file_name='{!s}-z_map_peaks.csv'.format(d_handler.tag), file_tag='z_peaks_csv')

            # Scripts
            d_handler.output_handler.add_dir(dir_name='scripts', dir_tag='scripts', top_dir_tag='root')
            d_handler.output_handler.add_file(file_name='load_maps.pml', file_tag='pymol_script', dir_tag='scripts')
            d_handler.output_handler.add_file(file_name='ccp4mg_{!s}_{!s}.py', file_tag='ccp4mg_script', dir_tag='scripts')

            # Output blobs
            d_handler.output_handler.add_dir(dir_name='blobs', dir_tag='blobs', top_dir_tag='root')
            d_handler.output_handler.add_file(file_name='blob_{!s}_{!s}.png', file_tag='ccp4mg_png', dir_tag='blobs')

            # Pickled objects
            d_handler.output_handler.add_dir(dir_name='pickles', dir_tag='pickles', top_dir_tag='root')
            d_handler.output_handler.add_file(file_name='dataset.pickle', file_tag='dataset_pickle', dir_tag='pickles')

            ##############################################################################################################

            # Link the input files to the output folder
            if not os.path.exists(d_handler.output_handler.get_file('input_structure')):
                os.symlink(d_handler.pdb_filename(), d_handler.output_handler.get_file('input_structure'))
            if not os.path.exists(d_handler.output_handler.get_file('input_data')):
                os.symlink(d_handler.mtz_filename(), d_handler.output_handler.get_file('input_data'))

            # Search for ligand files and link them to the output ligands folder
            lig_glob  = os.path.join(os.path.dirname(d_handler.pdb_filename()), lig_style)
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
        self.log('{!s} Datasets Loaded (New).          '.format(len(loaded_datasets), True))
        self.log('{!s} Datasets Loaded (Total).        '.format(self.datasets.size(), True))

    def select_reference_dataset(self, method='resolution', max_rfree=0.4, min_resolution=5):
        """Select dataset to act as the reference - scaling, aligning etc"""

        assert method in ['resolution','rfree'], 'METHOD FOR SELECTING THE REFERENCE DATASET NOT RECOGNISED: {!s}'.format(method)

        # Filter the datasets that can be selected as the reference dataset
        filtered_datasets = self.datasets.mask(mask_name='no_build', invert=True)

        if self.args.input.reference.pdb and self.args.input.reference.pdb:
            self.log('===================================>>>', True)
            self.log('Reference Provided by User', True)
            return self.args.input.reference.pdb, self.args.input.reference.mtz
        else:
            self.log('===================================>>>', True)
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

    def load_reflection_data(self, ampl_label, phas_label):
        """Extract amplitudes and phases for creating map"""

        if self.args.input.reference.ampl_label:    ref_ampl_label = self.args.input.reference.ampl_label
        else:                                       ref_ampl_label = ampl_label
        if self.args.input.reference.phas_label:    ref_phas_label = self.args.input.reference.phas_label
        else:                                       ref_phas_label = phas_label

        # Extract reflection data for the reference dataset
        self.reference_dataset().sfs = extract_structure_factors(self.reference_dataset().reflection_data(), ampl_label=ref_ampl_label, phas_label=ref_phas_label)

        t1 = time.time()
        self.log('===================================>>>', True)
        for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True):
            # Extract miller array of structure factors
            if d_handler.sfs != None:
                print('\rAlready Loaded: Dataset {!s}          '.format(d_handler.tag), end=''); sys.stdout.flush()
            else:
                print('\rLoading: Dataset {!s}                 '.format(d_handler.tag), end=''); sys.stdout.flush()
                d_handler.sfs = extract_structure_factors(mtz_object=d_handler.reflection_data(), ampl_label=ampl_label, phas_label=phas_label)
        t2 = time.time()
        self.log('\r> Structure Factors Extracted > Time Taken: {!s} seconds'.format(int(t2-t1)), True)

#    def scale_datasets(self, ampl_label, phas_label):
#        """OLD - Scale the reflection data to the reference (now done on the maps in real space)"""
#
#        t1 = time.time()
#        self.log('===================================>>>', True)
#        for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True):
#            if d_handler.scaled_sfs != None:
#                print('\rAlready Scaled: Dataset {!s}          '.format(d_handler.tag), end=''); sys.stdout.flush()
#                continue
#            else:
#                print('\rScaling: Dataset {!s}                 '.format(d_handler.tag), end=''); sys.stdout.flush()
#
#            # Extract just the amplitudes for scaling
#            d_amps = unscaled_sfs.amplitudes()
#            d_phas = unscaled_sfs.phases()
#
#            assert d_amps.is_real_array(), 'AMPLITUDES SHOULD BE REAL?!'
#            assert d_phas.is_real_array(), 'PHASE ANGLES SHOULD BE REAL?!'
#
#            # Scale new data to the reference dataset
#            scaled_amps = apply_simple_scaling(miller=d_amps, ref_miller=ref_amps)
#            # Recombine the scaled amplitudes and the phases
#            scaled_sf_real = scaled_amps.data() * flex.cos(d_phas.data()*math.pi/180.0)
#            scaled_sf_imag = scaled_amps.data() * flex.sin(d_phas.data()*math.pi/180.0)
#            scaled_sf_com = flex.complex_double(reals=scaled_sf_real, imags=scaled_sf_imag)
#            # Create a copy of the old array with the new scaled structure factors
#            scaled_sfs = unscaled_sfs.array(data=scaled_sf_com)
#            # Set the scaled structure factors
#            d_handler.unscaled_sfs = unscaled_sfs
#            d_handler.scaled_sfs = scaled_sfs
#        self.log('\rDatasets Scaled.               ', True)
#        t2 = time.time()
#        print('> SCALING STRUCTURE FACTORS > Time Taken: {!s} seconds'.format(int(t2-t1)))

    def align_datasets(self, method):
        """Align each structure the reference structure"""

        def align_map_func(arg_dict):
            global_mx, local_mxs = align_dataset_to_reference(
                                            d_handler   = arg_dict['d_handler'],
                                            ref_handler = arg_dict['r_handler'],
                                            method      = arg_dict['method']   )
            return (arg_dict['d_handler'].tag, (global_mx, local_mxs))

        assert method in ['local','global'], 'METHOD NOT DEFINED: {!s}'.format(method)

        # If local alignment has been chosen, also do a global alignment
        if method == 'local': method = 'both'

        # Align the datasets using multiple cores if possible
        self.log('===================================>>>', True)
        print('Generating Alignments (using {!s} cores)'.format(self.args.settings.cpus))
        start = time.time()
        arg_list = [{'r_handler':self.reference_dataset(), 'd_handler':d_handler, 'method':method} for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True)]
        alignment_transforms = easy_mp.pool_map(fixed_func=align_map_func, args=arg_list, processes=self.args.settings.cpus)
        alignment_transforms = dict(alignment_transforms)
        finish = time.time()
        self.log('\r> Generating Alignments > Time Taken: {!s} seconds'.format(int(finish-start)), True)

        # Post-process the alignments (write out aligned structures etc)
        t1 = time.time()
        for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True):
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
                if self._ref_dataset_index == None:
                    self._ref_dataset_index = d_handler.num
                    print('REFERENCE FOUND! {!s}                          '.format(d_handler.tag))
                # Raise error if the reference has already been set and it's not this dataset
                elif self._ref_dataset_index != d_handler.num:
                    raise Exception('ALIGNED OBJECT EQUAL TO UNALIGNED OBJECT - THIS IS MOST UNLIKELY')

        t2 = time.time()
        self.log('\r> Aligning Structures > Time Taken: {!s} seconds'.format(int(t2-t1)), True)

    def analyse_dataset_variability_1(self):
        """Go through all of the datasets and collect lots of different characteristics of the datasets for identifying odd datasets"""

        self.log('===================================>>>', True)
        self.log('Collecting Dataset/Crystal Variation Data - 1', True)
        self.log('===================================>>>')

        self.log('Extracting Resolutions')
        for d in self.datasets.all():
            # Resolution info
            self.tables.dataset_info.set_value(d.tag, 'high_resolution', d.mtz_summary.high_res)
            self.tables.dataset_info.set_value(d.tag, 'low_resolution', d.mtz_summary.low_res)
            # Unit cell info
            self.tables.dataset_info.set_value(d.tag, ['uc_a','uc_b','uc_c','uc_alpha','uc_beta','uc_gamma'], d.mtz_summary.unit_cell.parameters())
            self.tables.dataset_info.set_value(d.tag, 'uc_vol', d.mtz_summary.unit_cell.volume())
            # Spacegroup info
            self.tables.dataset_info.set_value(d.tag, 'space_group', d.mtz_summary.space_group.info().type().lookup_symbol())
            # Quality info
            self.tables.dataset_info.set_value(d.tag, 'r_work', d.input().get_r_rfree_sigma().r_work)
            self.tables.dataset_info.set_value(d.tag, 'r_free', d.input().get_r_rfree_sigma().r_free)

    def filter_datasets_1(self):
        """Filter out the datasets which contain different protein models (i.e. protein length, sequence, etc)"""

        self.log('===================================>>>', True)
        self.log('Filtering Datasets (Non-identical structures). Potential Classes:', True)
        for failure_class in self._structure_mask_names:
            self.log('\t{!s}'.format(failure_class), True)
        self.log('===================================>>>', True)

        ref_counts, ref_dict = self.reference_dataset().get_structure_summary()
        ref_handler = self.reference_dataset()

        # Check that the same protein structure is present in each dataset - THIS MASK SHOULD INCLUDE ALL DATASETS AT FIRST
        for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True):

            my_counts, my_dict = d_handler.get_structure_summary()

            print('\rFiltering Dataset {!s}          '.format(d_handler.tag), end=''); sys.stdout.flush()
            # Check the space group of the dataset
            if d_handler.input().crystal_symmetry().space_group().info().symbol_and_number() != ref_handler.input().crystal_symmetry().space_group().info().symbol_and_number():
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.tag))
                self.log('Different Space Group')
                self.log('Reference: {!s}, {!s}: {!s}'.format(ref_handler.input().crystal_symmetry().space_group().info().symbol_and_number(),
                                                        d_handler.tag, d_handler.input().crystal_symmetry().space_group().info().symbol_and_number()))
                self.log('===================================>>>')
                self.datasets.all_masks().set_mask_value(mask_name='bad structure - different space group', entry_id=d_handler.tag, value=True)
            # Check that the hierarchies are identical
            if not d_handler.hierarchy().is_similar_hierarchy(ref_handler.hierarchy()):
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.tag))
                self.log('Non-Identical Hierarchy')
                self.log('Reference: {!s}, {!s}: {!s}'.format(ref_counts.n_chains, d_handler.tag, my_counts.n_chains))
                self.log('===================================>>>')
                self.datasets.all_masks().set_mask_value(mask_name='bad structure - non-identical structures', entry_id=d_handler.tag, value=True)
            # Check the number of chains
            elif my_counts.n_chains != ref_counts.n_chains:
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.tag))
                self.log('Different Number of Chains')
                self.log('Reference: {!s}, {!s}: {!s}'.format(ref_counts.n_chains, d_handler.tag, my_counts.n_chains))
                self.log('===================================>>>')
                self.datasets.all_masks().set_mask_value(mask_name='bad structure - chain counts', entry_id=d_handler.tag, value=True)
            # Check the ids of the chains
            elif my_counts.chain_ids != ref_counts.chain_ids:
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.tag))
                self.log('Different Chain IDs')
                self.log('Reference: {!s}, {!s}: {!s}'.format(ref_counts.chain_ids, d_handler.tag, my_counts.chain_ids))
                print('===================================>>>')
                self.datasets.all_masks().set_mask_value(mask_name='bad structure - chain ids', entry_id=d_handler.tag, value=True)
            # Check the sequences of the chains
            elif my_dict['chain_sequences'] != ref_dict['chain_sequences']:
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.tag))
                self.log('Different Sequences')
                for chain_id in ref_dict['chain_sequences'].keys():
                    self.log('Chain {!s}: Reference - {!s}'.format(chain_id, ref_dict['chain_sequences'][chain_id]))
                    self.log('Chain {!s}: {:>9s} - {!s}'.format(chain_id, my_dict['chain_sequences'][chain_id]))
                print('===================================>>>')
                self.datasets.all_masks().set_mask_value(mask_name='bad structure - chain sequences', entry_id=d_handler.tag, value=True)
            # Check the number of residues - TODO not sure it can ever get here... remove?
            elif my_counts.n_residues != ref_counts.n_residues:
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.tag))
                self.log('Different Number of Residues')
                self.log('Reference: {!s}, {!s}: {!s}'.format(ref_counts.n_residues, d_handler.tag, my_counts.n_residues))
                print('===================================>>>')
                self.datasets.all_masks().set_mask_value(mask_name='bad structure - residue counts', entry_id=d_handler.tag, value=True)
            # Check the number of atoms
            elif my_counts.n_atoms != ref_counts.n_atoms:
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.tag))
                self.log('Different Number of Atoms')
                self.log('Reference: {!s}, {!s}: {!s}'.format(ref_counts.n_atoms, d_handler.tag, my_counts.n_atoms))
                print('===================================>>>')
                self.datasets.all_masks().set_mask_value(mask_name='bad structure - atom counts', entry_id=d_handler.tag, value=True)
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
            reject_dir = os.path.join(self.output_handler.get_dir('rejected_datasets'), d_handler.name)
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

            print('\rFiltering Dataset {!s}          '.format(d_handler.tag), end=''); sys.stdout.flush()
            # Check that it correlates well with itself before and after scaling
            if d_handler.input().get_r_rfree_sigma().r_free > self.params.filtering.max_rfree:
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.tag))
                self.log('RFree is higher than cutoff: {!s}'.format(self.params.filtering.max_rfree))
                self.log('High RFree: {!s}'.format(d_handler.input().get_r_rfree_sigma().r_free))
                self.log('===================================>>>')
                self.datasets.all_masks().set_mask_value(mask_name='bad crystal - rfree', entry_id=d_handler.tag, value=True)
#            elif d_handler.scaled_sfs.amplitudes().correlation(d_handler.unscaled_sfs.amplitudes()).coefficient() < self.params.filtering.min_correlation_to_reflection_data:
#                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.tag))
#                self.log('Low correlation between scaled and unscaled data')
#                self.log('Scaled-Unscaled Correlation: {!s}'.format(d_handler.scaled_sfs.amplitudes().correlation(d_handler.unscaled_sfs.amplitudes()).coefficient()))
#                self.log('===================================>>>')
#                self.datasets.all_masks().set_mask_value(mask_name='bad crystal - data correlation', entry_id=d_handler.tag, value=True)
            # Check the deviation from the average sites
            elif d_handler.get_calpha_sites().rms_difference(d_handler.transform_from_reference(points=self.reference_dataset().get_calpha_sites(), method='global')) > self.params.filtering.max_rmsd_to_reference:
                self.log('\rRejecting Dataset: {!s}          '.format(d_handler.tag))
                self.log('C-alpha RMSD is too large')
                self.log('Aligned (Calpha) RMSD: {!s}'.format(d_handler.get_calpha_sites().rms_difference(d_handler.transform_from_reference(points=self.reference_dataset().get_calpha_sites(), method='global'))))
                self.log('===================================>>>')
                self.datasets.all_masks().set_mask_value(mask_name='bad crystal - isomorphous structure', entry_id=d_handler.tag, value=True)
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
            reject_dir = os.path.join(self.output_handler.get_dir('rejected_datasets'), d_handler.name)
            if not os.path.exists(reject_dir):
                rel_symlink(orig=d_handler.output_handler.get_dir('root'), link=reject_dir)

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

        # Create empty mask
        self.datasets.all_masks().add_mask(mask_name=analysis_mask_name, mask=[False]*self.datasets.size())
        # Select from the datasets that haven't been rejected
        for d_handler in self.datasets.mask(mask_name='rejected - total', invert=True):
            # Check the resolution of the dataset (is not too low)
            if self.tables.dataset_info.get_value(index=d_handler.tag, col='high_resolution') > high_res_large_cutoff:
                continue
            # Check the resolution of the dataset (is not too high)
            elif self.tables.dataset_info.get_value(index=d_handler.tag, col='high_resolution') < high_res_small_cutoff:
                continue
            # Check to see if this has been excluded from building
            elif self.datasets.all_masks().get_mask_value(mask_name='no_analyse', entry_id=d_handler.tag) == True:
                self.log('Rejecting Dataset {!s}: Excluded from analysis'.format(d_handler.tag))
                continue
            elif (not self.args.method.reprocess_existing_datasets) and self.datasets.all_masks().get_mask_value(mask_name='old datasets', entry_id=d_handler.tag):
                self.log('Rejecting Dataset {!s}: Already Processed (Old Dataset)'.format(d_handler.tag))
                continue
            else:
                self.datasets.all_masks().set_mask_value(mask_name=analysis_mask_name, entry_id=d_handler.tag, value=True)
        return self.datasets.all_masks().get_mask(analysis_mask_name)

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
        for d in self.datasets.mask(mask_name='rejected - total', invert=True):
            rmsd = d.get_calpha_sites().rms_difference(d.transform_from_reference(points=self.get_calpha_average_sites(), method='global'))
            self.tables.dataset_info.set_value(d.tag, 'rmsd_to_mean', rmsd)

    def truncate_scaled_data(self, dataset_handlers, truncation_stuff=None):
        """Truncate data at the same indices across all the datasets"""

        self.log('===================================>>>', True)
        self.log('Truncating Reflection Data', True)

        # Calculate which reflections are present in all datasets
        common_set = self.reference_dataset().sfs.set()
        ref_size = common_set.size()
        self.log('Number of Reflections in Reference Dataset: {!s}'.format(ref_size))
        # Create maps for all of the datasets (including the reference dataset)
        for i_dh, d_handler in enumerate(dataset_handlers):
            common_set = common_set.common_set(d_handler.sfs, assert_is_similar_symmetry=False)
#            self.log('After Dataset {!s} - Remaining Common Reflections: {!s}'.format(i_dh+1, common_set.size()))

        self.log('===================================>>>', True)
        self.log('Number of Common Reflections between Datasets: {!s} ({!s}% of reference)'.format(common_set.size(), int(100.0*common_set.size()/ref_size)))
        self.log('After Truncation - Reflections per dataset: {!s}'.format(common_set.size()))

        # Create maps for all of the datasets (including the reference dataset)
        for d_handler in [self.reference_dataset()]+dataset_handlers:
            # At the moment save the 'truncated' sfs as the whole list
            # TODO
            d_handler.tr_sfs = d_handler.sfs.common_set(common_set, assert_is_similar_symmetry=False)

        reslns = [d.tr_sfs.d_min() for d in dataset_handlers]
        min_res = min(reslns)
        max_res = max(reslns)

        self.log('After Truncation - Resolution Range: {!s}-{!s}'.format(min_res, max_res))

        # TODO WRITE HISTOGRAMS OF RESOLUTIONS?!

    def load_reference_map(self, map_resolution=0):
        """Load the reference map, and calculate some map statistics"""

        # Reference dataset handler
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
        masked_gps = list(self.reference_grid().global_mask().outer_mask())
        masked_idxs = self.reference_grid().global_mask().outer_mask_indices()
        masked_cart_ref = flex.vec3_double(masked_gps)*self.reference_grid().grid_spacing()
        # Sample the map at the masked points
        masked_vals_ref = native_map_handler.get_cart_values(masked_cart_ref)

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

        def load_maps_map_func(arg_dict):
            """Function to allow maps to be loaded in parallel"""

            d_handler            = arg_dict['d_handler']
            grid                 = arg_dict['grid']
            params               = arg_dict['params']
            map_resolution       = arg_dict['map_resolution']
            ref_map_holder       = arg_dict['ref_map_holder']
            masked_cart_ref      = arg_dict['masked_cart_ref']
            masked_cart_mappings = arg_dict['masked_cart_mappings']

            print('\rLoading Maps for Dataset {!s}'.format(d_handler.tag), end=''); sys.stdout.flush()

            # Take the scaled diffraction data for each dataset and create fft
            fft_map = d_handler.tr_sfs.fft_map( resolution_factor = params.maps.resolution_factor,
                                                d_min             = map_resolution,
                                                symmetry_flags    = cctbx.maptbx.use_space_group_symmetry  )

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
            masked_cart_d = d_handler.transform_from_reference(points=masked_cart_ref, method=params.alignment.method, point_mappings=masked_cart_mappings)
            # Sample the map at these points
            masked_vals_d = native_map_handler.get_cart_values(masked_cart_d)

            # Normalise map - Scale it to the same scale as the reference dataset on the same points
            d_map_mean = masked_vals_d.min_max_mean().mean
            d_map_rms  = masked_vals_d.standard_deviation_of_the_sample()
            masked_vals_d = (masked_vals_d - d_map_mean)*(ref_map_holder.meta.map_rms/d_map_rms) + ref_map_holder.meta.map_mean

            # Calculate the sparse vector of the masked map values
            morphed_map_sparse = scitbx.sparse.vector(grid.size_1d(), dict(zip(masked_idxs, masked_vals_d)))
            morphed_map        = morphed_map_sparse.as_dense_vector()

            # Reshape into right shape of the grid
            morphed_map.reshape(grid)

            # Create map holder
            map_holder = MapHolder( num         = d_handler.num,
                                    tag         = d_handler.tag,
                                    map         = morphed_map,
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

        self.log('===================================>>>', True)
        self.log('Loading Electron Density Maps @ {!s}A'.format(map_resolution), True)

        # Create holder for the output map objects
        map_holder_list = MapHolderList()

        # Extract the points for the morphed maps (in the reference frame)
        masked_gps = list(self.reference_grid().global_mask().outer_mask())
        masked_idxs = self.reference_grid().global_mask().outer_mask_indices()
        masked_cart_ref = flex.vec3_double(masked_gps)*self.reference_grid().grid_spacing()
        # Mapping of grid points to rotation matrix keys (residue CA labels)
        if self.params.alignment.method == 'local':
            masked_cart_mappings = self.reference_grid().partition().query_by_grid_indices(masked_idxs)
        else:
            masked_cart_mappings = None

        # Load maps using multiple cores if possible
        self.log('===================================>>>', True)
        print('Loading Maps (using {!s} cores)'.format(self.args.settings.cpus))
        start = time.time()
        arg_list = [{   'd_handler':d_handler, 'params':self.params, 'map_resolution':map_resolution, \
                        'grid':self.reference_grid().grid_indexer(), 'ref_map_holder':ref_map_holder, \
                        'masked_cart_ref':masked_cart_ref, 'masked_cart_mappings':masked_cart_mappings      } for d_handler in dataset_handlers]
        map_holders = easy_mp.pool_map(fixed_func=load_maps_map_func, args=arg_list, processes=self.args.settings.cpus, chunksize=1)
        # Append to the map holder list
        map_holder_list.add(map_holders)
        # Go through and assign parents
        for mh in map_holder_list.all():
            mh.parent = self.datasets.get(tag=mh.tag)

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

    def print_clustering_settings(self, z_cutoff,
                                        min_cluster_volume,
                                        min_cluster_z_peak,
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
        self.log('Minimum Cluster Z-Peak:      {!s}'.format(min_cluster_z_peak), True)
        self.log('===================================>>>', True)

    def cluster_high_z_values(self, d_handler,
                                    z_map,
                                    z_cutoff,
                                    point_mask,
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

        # Look for z-clusters on the mask
        d_selected_points = [(gp, z_map[grid_indexer(gp)]) for gp in point_mask if z_map[grid_indexer(gp)] >= z_cutoff]

        # No Cluster points found
        if not d_selected_points:
            self.log('Dataset {!s}: No Clusters Found'.format(d_handler.tag), True)
            return 0, {}
        # Can't cluster if there are too many points
        elif len(d_selected_points) > 10000:
            self.log('Dataset {!s}: Too many points to cluster: {!s} Points.'.format(d_handler.tag, len(d_selected_points)), True)
            num_clusters = -1
            z_clusters = {}
            # This dataset is too noisy to analyse - flag!
            self.datasets.all_masks().set_mask_value(mask_name='noisy zmap', entry_id=d_handler.tag, value=True)
            # Link datasets to the noisy results directory
            noisy_dir = os.path.join(self.output_handler.get_dir('noisy_datasets'), d_handler.name)
            if not os.path.exists(noisy_dir):
                rel_symlink(orig=d_handler.output_handler.get_dir('root'), link=noisy_dir)
        # Cluster points if we have found them
        else:
            self.log('Dataset {!s}: Clustering {!s} Point(s).'.format(d_handler.tag, len(d_selected_points)), True)
            # Points found for this cluster!
            if len(d_selected_points) == 1:
                # Only 1 point - 1 cluster! - edge case where we can't do clustering
                num_clusters = 1
                # Dictionary of results
                z_clusters = {1: d_selected_points}
            else:
                # Extract only the coordinates and form an array
                point_array = numpy.array([tup[0] for tup in d_selected_points])
                # Cluster the extracted points
                t1 = time.time()
                clusts = list(scipy.cluster.hierarchy.fclusterdata( X=point_array,
                                                                    t=grid_clustering_cutoff,
                                                                    criterion=clustering_criterion,
                                                                    metric=clustering_metric,
                                                                    method=clustering_method    )   )
                t2 = time.time()
                self.log('> Clustering > Time Taken: {!s} seconds'.format(int(t2-t1)))

                # Get the number of clusters
                num_clusters = max(clusts)
                # Initialise dictionary to hold the points for the different clusters (under the new indexing)
                z_clusters = dict([(i+1, []) for i in range(num_clusters)])
                # Populate the clusters according to the clustering
                [z_clusters[c_idx].append(d_selected_points[p_idx]) for p_idx, c_idx in enumerate(clusts)]

            # Check that the clusters are numbered properly
            assert num_clusters == max(z_clusters.keys())
            assert num_clusters == len(z_clusters.keys())

        return num_clusters, z_clusters

    def filter_z_clusters_1(self, z_clusters,
                                  min_cluster_volume,
                                  min_cluster_z_peak    ):
        """Filter the z-clusters on a variety of criteria (size, peak value)"""

        self.log('===================================>>>')
        self.log('FILTERING (STEP 1): Filtering by blob size and peak value')

        # Calculate the approximate minimum volume for the cluster size
        min_cluster_size = int(min_cluster_volume/(self.reference_grid().grid_spacing()**3))

        # Filter out small clusters - get numbers of clusters satisfying the minimum cluster size
        large_clusters = [c_num for c_num in range(1,len(z_clusters)+1) if len(z_clusters[c_num]) >= min_cluster_size]
        if len(large_clusters) == 0:
            return 0, {}

        # Filter out weak clusters - get numbers of clusters satisfying the minimum z_peak value
        strong_clusters = [c_num for c_num in large_clusters if max([c[1] for c in z_clusters[c_num]]) >= min_cluster_z_peak]
        if len(strong_clusters) == 0:
            return 0, {}

        # Renumber the filtered clusters
        filtered_z_clusters = dict([(new_c_idx+1, z_clusters[old_c_num]) for new_c_idx, old_c_num in enumerate(strong_clusters)])
        # Check that the clusters are numbered properly
        num_clusters = len(strong_clusters)
        assert num_clusters == max(filtered_z_clusters.keys())
        assert num_clusters == len(filtered_z_clusters.keys())

        self.log('Filtered {!s} Clusters to {!s} Clusters'.format(len(z_clusters), len(filtered_z_clusters)))

        return num_clusters, filtered_z_clusters

    def filter_z_clusters_2(self, z_clusters,
                                  grid_spacing,
                                  grid_origin_cart,
                                  ref_structure,
                                  min_contact_dist=6    ):
        """Find and remove clusters more than a minimum distance from the protein"""

        # min_contact_dist - blobs are rejected if they are more than this distance from the protein

        self.log('===================================>>>')
        self.log('FILTERING (STEP 2): Filtering by minimum distance from protein')

        # Extract the protein sites
        cache = ref_structure.atom_selection_cache()
        protein_bool = cache.selection('pepnames')
        ref_sites_cart = ref_structure.select(protein_bool).atoms().extract_xyz()

        # Save time - calculate the square of the contact distance
        min_contact_dist_sq = min_contact_dist**2

        # Remove any clusters that are more than min_contact_dist from the protein
        filtered_c_keys = []
        for c_key in sorted(z_clusters.keys()):
            # Extract points in cluster
            cluster_points_cart = (flex.vec3_double([c[0] for c in z_clusters[c_key]]) * grid_spacing) - grid_origin_cart
            # Calculate minimum distance to protein
            for r_site_cart in ref_sites_cart:
                diff_vecs_cart = cluster_points_cart - r_site_cart
                # Keep key if minimum distance is less than min_contact_dist
                if min(diff_vecs_cart.dot()) < min_contact_dist_sq:
                    filtered_c_keys.append(c_key)
                    break

            if self.args.settings.verbose:
                if filtered_c_keys and (filtered_c_keys[-1] == c_key):
                    print('KEEPING CLUSTER:', c_key)
                else:
                    print('REJECTING CLUSTER:', c_key, '\t', cluster_points_cart[0])

        # Check if all clusters have been rejected
        if not filtered_c_keys:
            num_clusters = 0
            filtered_z_clusters = {}
        else:
            # Renumber the filtered clusters
            filtered_z_clusters = dict([(new_c_idx+1, z_clusters[old_c_key]) for new_c_idx, old_c_key in enumerate(filtered_c_keys)])
            # Check that the clusters are numbered properly
            num_clusters = len(filtered_c_keys)
            assert num_clusters == max(filtered_z_clusters.keys())
            assert num_clusters == len(filtered_z_clusters.keys())

        self.log('Filtered {!s} Clusters to {!s} Clusters'.format(len(z_clusters), len(filtered_z_clusters)))

        return num_clusters, filtered_z_clusters

    def filter_z_clusters_3(self, z_clusters,
                                  grid_spacing,
                                  grid_origin_cart,
                                  ref_unit_cell,
                                  ref_sym_ops,
                                  ref_structure,
                                  max_contact_dist=8  ):
        """Find and remove symmetry equivalent clusters"""

        # max_contact_dist - a point contacts an atom if the atoms is within this distance of it

        if len(z_clusters) == 1:
            return 1, z_clusters
        else:
            self.log('===================================>>>')
            self.log('FILTERING (STEP 3): Filtering symmetry equivalent clusters')

        # Extract the protein sites
        cache = ref_structure.atom_selection_cache()
        protein_bool = cache.selection('pepnames')
        ref_sites_cart = ref_structure.select(protein_bool).atoms().extract_xyz()

        # Save time - calculate the square of the contact distance
        max_contact_dist_sq = max_contact_dist**2
        # Minimum distance between cart points to be equiv (due to sampling on grid)
        dist_cut_sq = 1.05 * 3 * grid_spacing * grid_spacing

        # Pull out the labels for the clusters
        cluster_keys = sorted(z_clusters.keys())

        # Cartesianise and fractionalise the points in each of the clusters
        points_frac = {}
        points_cart = {}
        for c_key in cluster_keys:
            # Calculate the cartesian coordinates of the grid points (and apply reference shift)
            points_cart[c_key] = (flex.vec3_double([c[0] for c in z_clusters[c_key]]) * grid_spacing) - grid_origin_cart
            # Fractionalise them to the unit cell of the reference structure
            points_frac[c_key] = ref_unit_cell.fractionalize(points_cart[c_key])

        # Matrix for whether they are near to each other
        equiv_sites = numpy.zeros([len(ref_sym_ops), len(cluster_keys), len(cluster_keys)], dtype=int)

        # Apply the symmetry operations to each cluster to see if it is near to other clusters
        for i_sym_op, sym_op in enumerate(ref_sym_ops):
            # Transformed clusters under this symmetry operation
            trans_points_frac = dict([(c_key, sym_op.as_rational().as_float() * points_frac[c_key]) for c_key in cluster_keys])
            # Loop through clusters
            for i_clust_1, c_key_1 in enumerate(cluster_keys):
                # Extract points for the un-transformed points
                ref_points = points_frac[c_key_1]
                # Loop through clusters again
                for i_clust_2, c_key_2 in enumerate(cluster_keys):
                    # Comparing cluster to itself - skip
                    if i_clust_1 == i_clust_2:
                        # These are cleary equivalent
                        equivalent = True
                    else:
                        # Start by assuming the clusters are not related
                        equivalent = False
                        # Extract locations for the transformed points
                        query_points = trans_points_frac[c_key_2]
                        # Loop through and see if the clusters overlap
                        for qp in query_points:
                            # Brute force - calculate the distance from all to all
                            diffs_frac = ref_points - qp
                            # Convert back to cartesian
                            diffs_cart = ref_unit_cell.orthogonalize(diffs_frac)
                            # Check if any are closer than the minimum required
                            if min(diffs_cart.dot()) < dist_cut_sq:
                                equivalent = True
                                break
                    # If the clusters overlap, they are equivalent
                    if equivalent:
                        equiv_sites[(i_sym_op, i_clust_1, i_clust_2)] = 1

        # Condense the cluster equivalence - take max over the symmetry operations and group by connected paths
        sym_cluster_groups = find_connected_groups(connection_matrix=equiv_sites.max(axis=0))
        # Number of unique clusters
        num_uniq_clusters = max(sym_cluster_groups)
        sym_cluster_zip = zip(sym_cluster_groups, cluster_keys)

        # Keys of clusters to keep
        clusters_to_keep = []
        # Iterate through the equivalent clusters and find the closest one to the protein
        for n_clust in range(1, num_uniq_clusters+1):
            # Extract the keys for the clusterss in this group
            clust_c_keys = [c[1] for c in sym_cluster_zip if c[0] == n_clust]
            # Save the number of contacts between the cluster and the protein
            c_contacts = []
            # Calculate the number of contacts each cluster has with the protein
            for c_key in clust_c_keys:
                # Initialise contact counter
                contacts = 0
                # Get the cartesian points for the cluster
                c_points_cart = points_cart[c_key]
                # Again, use the brute force all-v-all method
                for rp in ref_sites_cart:
                    diffs_cart = c_points_cart - rp
                    # Check to see if site closer to cluster than minimum
                    if min(diffs_cart.dot()) < max_contact_dist_sq:
                        contacts += 1
                # Record the number of contacts (over size of cluster)
                c_contacts.append(1.0*contacts/len(c_points_cart))
                if self.args.settings.verbose:
                    print('CLUSTER:', c_key, ', CONTACTS PER POINT:', round(c_contacts[-1],3))

            # Find the cluster with the most contacts
            max_contacts = max(c_contacts)

            if max_contacts == 0:
                raise Exception('MAX CONTACTS IS 0!')
            else:
                clusters_to_keep.append(clust_c_keys[c_contacts.index(max_contacts)])
                if self.args.settings.verbose:
                    print('KEEPING CLUSTER', clusters_to_keep[-1])

        assert len(clusters_to_keep) == num_uniq_clusters, 'NUMBER OF BLOBS AND BLOBS TO BE RETURNED NOT THE SAME'

        # Select and renumber the filtered clusters to be returned
        filtered_num_clusters = num_uniq_clusters
        filtered_z_clusters = dict([(new_c_idx+1, z_clusters[old_c_key]) for new_c_idx, old_c_key in enumerate(clusters_to_keep)])

        self.log('Filtered {!s} Clusters to {!s} Clusters'.format(len(z_clusters), len(filtered_z_clusters)))

        return filtered_num_clusters, filtered_z_clusters

    def group_clusters(self, z_clusters,
                             grid_spacing,
                             separation_cutoff=5   ):
        """Join clusters that are separated by less than max_separation"""

        if len(z_clusters) == 1:
            return 1, z_clusters
        else:
            self.log('===================================>>>')
            self.log('GROUPING: Grouping Nearby Clusters')

        # Minimum distance between grid points to be joined (squared)
        grid_spacing_cutoff_sq = (separation_cutoff/grid_spacing)**2

        # Record which clusters are to be joined
        connect_array = numpy.zeros((len(z_clusters),len(z_clusters)), dtype=int)

        cluster_keys = sorted(z_clusters.keys())

        for i_clust_1, c_key_1 in enumerate(cluster_keys):
            # Extract grid points for cluster 1
            c_gp_1 = flex.vec3_double([c[0] for c in z_clusters[c_key_1]])
            for i_clust_2, c_key_2 in enumerate(cluster_keys):
                # Skip if this is the same blob
                if i_clust_1 == i_clust_2:
                    connect_array[(i_clust_1, i_clust_2)] = 1
                    continue
                # Extract grid points for cluster 2
                c_gp_2 = flex.vec3_double([c[0] for c in z_clusters[c_key_2]])
                # Extract the minimum separation of the grid points
                min_dist_sq = min([min((c_gp_2 - gp).dot()) for gp in c_gp_1])
                # Check to see if they should be joined
                if min_dist_sq < grid_spacing_cutoff_sq:
                    connect_array[(i_clust_1, i_clust_2)] = 1

        # Find which clusters to join
        cluster_groupings = find_connected_groups(connection_matrix=connect_array)
        num_groups = max(cluster_groupings)
        grouped_zip = zip(cluster_groupings, cluster_keys)

        if self.args.settings.verbose:
            print('COMBINING', len(z_clusters), 'GROUPS INTO', num_groups, 'GROUPS')

        # Output cluster dict
        grouped_clusters = {}

        # Iterate through the groups of clusters and combine each group into one cluster
        for n_group in range(1, num_groups+1):
            # Extract the keys for the clusters in this group
            group_c_keys = [c[1] for c in grouped_zip if c[0] == n_group]
            # Extract points for these clusters
            group_points = []; [group_points.extend(z_clusters[c_key]) for c_key in group_c_keys]
            # Combine all of the points into a new cluster
            grouped_clusters[n_group] = group_points

        assert len(grouped_clusters) == num_groups
        assert sum(map(len, z_clusters.values())) == sum(map(len, grouped_clusters.values()))

        print(map(len, z_clusters.values()),'->',map(len, grouped_clusters.values()))
        self.log('Grouped {!s} Clusters together to form {!s} Clusters'.format(len(z_clusters), len(grouped_clusters)))

        return num_groups, grouped_clusters

    def collate_event_counts(self):
        """Collate eventss from all of the datasets"""

        self.log('===================================>>>', False)
        self.log('Collating Clusters', False)

        # List of points to be returned
        all_dataset_events = dict([(d.tag, d.events) for d in self.datasets.all()])

        # Print Cluster Summaries
        event_num = [(k, len(all_dataset_events[k])) for k in sorted(all_dataset_events.keys()) if all_dataset_events[k]]
        event_total = sum([a[1] for a in event_num])

        return event_total, event_num, all_dataset_events

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
                                                                'pos-contour' : [2,3,4,5]
                                                            }
                                            })

            # Write out the ccp4mg script to the dataset's scripts folder
            with open(view_script, 'w') as fh:
                fh.write(ccp4mg_script)

            # Make the images
            c = CommandManager('ccp4mg')
            c.SetArguments(['-norestore','-picture', view_script, '-R', view_image, '-RO', """'{"size":"1500x1500"}'""", '-quit'])
            c.Run()

            if not os.path.exists(view_image):
                print('FAILED TO MAKE IMAGES')
                print(c.err)

    def write_html_analysis_summary(self):
        """Writes an html summary of the datasets"""

        # Get template to be filled in
        template = PANDDA_HTML_ENV.get_template('output_table.html')
        # XXX 0=success, 1=none, 2=info, 3=warning, 4=failure

        # Extract the dataset info as a dictionary
        all_info = self.tables.dataset_info.transpose().to_dict()

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

            d_data = all_info[d.tag]
            rmsd  = round(d_data['rmsd_to_mean'], 3)
            rfree = round(d_data['r_free'], 3)
            rwork = round(d_data['r_work'], 3)

            columns = []
            overall_success = [0]
            # ------------------------------>>>
            # Test for Data Quality
            # ------------------------------>>>
            if self.datasets.all_masks().get_mask_value(mask_name='rejected - crystal', entry_id=d.tag) == True:
                columns.append({'flag':4,'message':'Rejected'.format(None)})
            else:
                columns.append({'flag':0,'message':'OK'})
            # ------------------------------>>>
            # Test for Identical Structures
            # ------------------------------>>>
            if self.datasets.all_masks().get_mask_value(mask_name='rejected - structure', entry_id=d.tag) == True:
                columns.append({'flag':4,'message':'Rejected'})
            else:
                columns.append({'flag':0,'message':'OK'})
            # ------------------------------>>>
            # Test for Refinement Success - some test on r-free
            # ------------------------------>>>
            if self.datasets.all_masks().get_mask_value(mask_name='bad crystal - rfree', entry_id=d.tag) == True:
                columns.append({'flag':4,'message':'RFree: {!s}'.format(rfree)})
            else:
                columns.append({'flag':0,'message':'RFree: {!s}'.format(rfree)})
            # ------------------------------>>>
            # Test for Structure movement
            # ------------------------------>>>
            if self.datasets.all_masks().get_mask_value(mask_name='bad crystal - isomorphous structure', entry_id=d.tag) == True:
                columns.append({'flag':4,'message':'RMSD: {!s}'.format(rmsd)})
            else:
                columns.append({'flag':0,'message':'RMSD: {!s}'.format(rmsd)})
            # ------------------------------>>>
            # Test for Structure movement
            # ------------------------------>>>
            if self.datasets.all_masks().get_mask_value(mask_name='analysed', entry_id=d.tag) == True:
                columns.append({'flag':0,'message':'Analysed'.format(rmsd)})
            else:
                columns.append({'flag':4,'message':'Rejected'.format(rmsd)})
            # ------------------------------>>>
            # Test for if it's interesting
            # ------------------------------>>>
            if self.datasets.all_masks().get_mask_value(mask_name='interesting', entry_id=d.tag) == True:
                columns.append({'flag':0,'message':'Clusters Found'})
            else:
                columns.append({'flag':1,'message':''})

            # Find largest error
            overall_success = max([c['flag'] for c in columns])

            output_data['table']['rows'].append({'heading':'Dataset {!s}'.format(d.tag),
                                                 'success':overall_success,
                                                 'columns':columns})

        with open(self.output_handler.get_file(file_tag='summary_table'), 'w') as out_html:
            out_html.write(template.render(output_data))

    def write_map_value_distribution(self, map_vals, output_file, plot_indices=None, plot_normal=False):
        """Write out the value distribution for a map"""
        if not self.args.output.plot_graphs: return
        if plot_indices: plot_vals = [map_vals[i] for i in plot_indices]
        else:            plot_vals = list(map_vals)
        graphs.map_value_distribution(f_name=output_file, plot_vals=plot_vals, plot_normal=plot_normal)

    def write_qq_plot_against_normal(self, map_vals, output_file, plot_indices=None):
        """Plot the values in map_vals against those expected from a normal distribution"""
        if not self.args.output.plot_graphs: return
        if plot_indices: plot_vals = [map_vals[i] for i in plot_indices]
        else:            plot_vals = list(map_vals)
        graphs.qq_plot_against_normal(f_name=output_file, plot_vals=plot_vals)

    def write_map_analyser_summary(self, map_analyser, analysis_mask_name):
        """Write statistical maps for a map_analyser object"""

        if self.args.output.plot_graphs:
            import matplotlib
            matplotlib.interactive(0)
            from matplotlib import pyplot

        ########################################################

        map_res = map_analyser.meta.resolution

        self.log('===================================>>>')
        self.log('=> Writing Summary of Analysis @ {!s}A'.format(map_res))

        ########################################################

        self.log('=> Writing Statistical Maps')

        self.write_array_to_map(output_file = self.output_handler.get_file('mean_map').format(map_res),
                                map_data    = map_analyser.statistical_maps.mean_map     )
        self.write_array_to_map(output_file = self.output_handler.get_file('stds_map').format(map_res),
                                map_data    = map_analyser.statistical_maps.stds_map     )
        self.write_array_to_map(output_file = self.output_handler.get_file('sadj_map').format(map_res),
                                map_data    = map_analyser.statistical_maps.sadj_map     )
        self.write_array_to_map(output_file = self.output_handler.get_file('skew_map').format(map_res),
                                map_data    = map_analyser.statistical_maps.skew_map     )
        self.write_array_to_map(output_file = self.output_handler.get_file('kurt_map').format(map_res),
                                map_data    = map_analyser.statistical_maps.kurt_map     )
        self.write_array_to_map(output_file = self.output_handler.get_file('bimo_map').format(map_res),
                                map_data    = map_analyser.statistical_maps.bimo_map     )

        ########################################################

        # Statistical Map Values
        masked_idxs = self.reference_grid().global_mask().outer_mask_indices()
        mean_map_vals = list(map_analyser.statistical_maps.mean_map.select(masked_idxs))
        try:    medn_map_vals = list(map_analyser.statistical_maps.medn_map.select(masked_idxs))
        except: medn_map_vals = mean_map_vals
        stds_map_vals = list(map_analyser.statistical_maps.stds_map.select(masked_idxs))
        sadj_map_vals = list(map_analyser.statistical_maps.sadj_map.select(masked_idxs))

        ########################################################

        # Dataset info
        d_info = self.tables.dataset_info
        # All datasets
        high_res = [d_info['high_resolution'][mh.tag] for mh in map_analyser.dataset_maps.all()]
        low_res =  [d_info['low_resolution'][mh.tag]  for mh in map_analyser.dataset_maps.all()]
        rfree =    [d_info['r_free'][mh.tag]          for mh in map_analyser.dataset_maps.all()]
        rwork =    [d_info['r_work'][mh.tag]          for mh in map_analyser.dataset_maps.all()]

        # Map info
        m_info = self.tables.dataset_map_info
        # All datasets
        map_uncties = [m_info['map_uncertainty'][mh.tag] for mh in map_analyser.dataset_maps.all()]
        # Analysed datasets only
        z_map_mean  = [m_info['z_map_mean'][mh.tag] for mh in map_analyser.dataset_maps.mask(mask_name=analysis_mask_name)]
        z_map_std   = [m_info['z_map_std'][mh.tag]  for mh in map_analyser.dataset_maps.mask(mask_name=analysis_mask_name)]
        z_map_skew  = [m_info['z_map_skew'][mh.tag] for mh in map_analyser.dataset_maps.mask(mask_name=analysis_mask_name)]
        z_map_kurt  = [m_info['z_map_kurt'][mh.tag] for mh in map_analyser.dataset_maps.mask(mask_name=analysis_mask_name)]

        ########################################################

        # Get the output directory to write the graphs into
        img_out_dir = os.path.join(self.output_handler.get_dir('analyses'), '{!s}A Maps'.format(map_analyser.meta.resolution))
        if not os.path.exists(img_out_dir): os.mkdir(img_out_dir)

        n_bins = 30

        ########################################################

        if not self.args.output.plot_graphs: return

        self.log('=> Writing Statistical Map Distributions')

        ##################################
        # STATISTICAL MAP HISTOGRAMS
        ##################################
        fig = pyplot.figure()
        pyplot.title('STATISTICAL MAP VALUES')
        # MEAN MAP
        pyplot.subplot(2, 1, 1)
        pyplot.hist(x=mean_map_vals, bins=n_bins)
        pyplot.xlabel('MEAN MAP DISTRIBUTION')
        pyplot.ylabel('COUNT')
        # MEDIAN MAP
        pyplot.subplot(2, 1, 2)
        pyplot.hist(x=medn_map_vals, bins=n_bins)
        pyplot.xlabel('MEDIAN MAP DISTRIBUTION')
        pyplot.ylabel('COUNT')
        # Apply tight layout to prevent overlaps
        pyplot.tight_layout()
        # Save both
        pyplot.savefig(os.path.join(img_out_dir, '{!s}A-mean_medn_map_vals.png'.format(map_res)))
        pyplot.close(fig)

        ##################################
        # MEAN Values v MEDIAN Values
        ##################################
        fig = pyplot.figure()
        pyplot.title('MEAN v MEDIAN MAP')
        pyplot.scatter(x=mean_map_vals, y=medn_map_vals)
        # Plot straight line between the min and max values
        min_val = min(mean_map_vals+medn_map_vals)
        max_val = max(mean_map_vals+medn_map_vals)
        pyplot.plot([min_val, max_val], [min_val, max_val], 'b--')
        # Axis labels
        pyplot.xlabel('MEAN MAP')
        pyplot.ylabel('MEDIAN MAP')
        # Apply tight layout to prevent overlaps
        pyplot.tight_layout()
        # Save
        pyplot.savefig(os.path.join(img_out_dir, '{!s}A-mean_v_median_scatter.png'.format(map_res)))
        pyplot.close(fig)

        ##################################
        # MEAN Values v MEDIAN Values
        ##################################
        graphs.med_mean_diff_hist(  f_name=os.path.join(img_out_dir, '{!s}A-mean_median_diff_hist.png'.format(map_res)),
                                    plot_vals=numpy.abs(flex.double(mean_map_vals)-flex.double(medn_map_vals)))

        ##################################
        # STATISTICAL MAP HISTOGRAMS
        ##################################
        fig = pyplot.figure()
        pyplot.title('STATISTICAL MAP VALUES')
        # STANDARD DEVIATION MAPS
        pyplot.subplot(2, 1, 1)
        pyplot.hist(x=stds_map_vals, bins=n_bins)
        pyplot.xlabel('STDS MAP DISTRIBUTION')
        pyplot.ylabel('COUNT')
        # ADJUSTED STANDARD DEVIATION MAPS
        pyplot.subplot(2, 1, 2)
        pyplot.hist(x=sadj_map_vals, bins=n_bins)
        pyplot.xlabel('ADJUSTED STDS MAP DISTRIBUTION')
        pyplot.ylabel('COUNT')
        # Apply tight layout to prevent overlaps
        pyplot.tight_layout()
        # Save both
        pyplot.savefig(os.path.join(img_out_dir, '{!s}A-stds_sadj_map_vals.png'.format(map_res)))
        pyplot.close(fig)

        ##################################
        # STD Values v ADJ STD Values
        ##################################
        fig = pyplot.figure()
        pyplot.title('RAW v ADJUSTED STDS')
        pyplot.scatter(x=stds_map_vals, y=sadj_map_vals)
        # Plot straight line between the min and max values
        min_val = min(stds_map_vals+sadj_map_vals)
        max_val = max(stds_map_vals+sadj_map_vals)
        pyplot.plot([min_val, max_val], [min_val, max_val], 'b--')
        # Axis labels
        pyplot.xlabel('RAW STDS')
        pyplot.ylabel('ADJUSTED STDS')
        # Apply tight layout to prevent overlaps
        pyplot.tight_layout()
        # Save
        pyplot.savefig(os.path.join(img_out_dir, '{!s}A-std_v_adj_std_scatter.png'.format(map_res)))
        pyplot.close(fig)

        ########################################################

        self.log('=> Writing Map Uncertainties')

        # MAP PARAMS
        fig = pyplot.figure()
        pyplot.title('MAP STATISTICS')
        # MAP UNCERTAINTIES
        pyplot.hist(x=map_uncties, bins=n_bins, range=(min(map_uncties)-0.1,max(map_uncties)+0.1))
        pyplot.xlabel('UNCERTAINTY OF MAP VALUES')
        pyplot.ylabel('COUNT')
        # Apply tight layout to prevent overlaps
        pyplot.tight_layout()
        # Save
        pyplot.savefig(os.path.join(img_out_dir, '{!s}A-d_map_uncertainties.png'.format(map_res)))
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
        pyplot.savefig(os.path.join(img_out_dir, '{!s}A-resolution_v_uncertainty.png'.format(map_res)))
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
        pyplot.savefig(os.path.join(img_out_dir, '{!s}A-resolution_v_rfree.png'.format(map_res)))
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
        pyplot.savefig(os.path.join(img_out_dir, '{!s}A-rfree_v_uncertainty.png'.format(map_res)))
        pyplot.close(fig)

        ########################################################

        self.log('=> Z-Map Distribution')

        # R-FACTORS
        fig = pyplot.figure()
        pyplot.title('Z-MAP DISTRIBUTION HISTOGRAMS')
        # RFree
        pyplot.subplot(2, 1, 1)
        pyplot.hist(x=z_map_mean, bins=n_bins, range=(min(z_map_mean)-0.1,max(z_map_mean)+0.1))
        pyplot.xlabel('Z-MAP MEAN')
        pyplot.ylabel('COUNT')
        # RWork
        pyplot.subplot(2, 1, 2)
        pyplot.hist(x=z_map_std, bins=n_bins, range=(0, max(z_map_std)+0.1))
        pyplot.xlabel('Z_MAP_STD')
        pyplot.ylabel('COUNT')
        # Apply tight layout to prevent overlaps
        pyplot.tight_layout()
        # Save both
        pyplot.savefig(os.path.join(img_out_dir, '{!s}A-z_map_statistics.png'.format(map_res)))
        pyplot.close(fig)

        # Z-MAP SKEW V UNCERTAINTY
        fig = pyplot.figure()
        pyplot.title('DATASET NORMALITY')
        pyplot.scatter(x=z_map_skew, y=z_map_kurt)
        pyplot.xlabel('SKEW')
        pyplot.ylabel('KURTOSIS')
        # Apply tight layout to prevent overlaps
        pyplot.tight_layout()
        # Save
        pyplot.savefig(os.path.join(img_out_dir, '{!s}A-z_map_skew_v_kurtosis.png'.format(map_res)))
        pyplot.close(fig)

    def write_dataset_summary_graphs(self):
        """Write out graphs of dataset variables"""

        if self.args.output.plot_graphs:
            import matplotlib
            matplotlib.interactive(0)
            from matplotlib import pyplot
        else:
            return

        def filter_nans(x):
            return [v for v in x if not numpy.isnan(v)]

        ########################################################

        self.log('===================================>>>')
        self.log('Generating Summary Graphs')

        # Extract dataset summary table
        d_info = self.tables.dataset_info
        # Get the output directory to write the graphs into
        img_out_dir = self.output_handler.get_dir('analyses')
        n_bins = 30

        ########################################################

        self.log('=> Data Quality Variation')

        # RESOLUTIONS
        fig = pyplot.figure()
        pyplot.title('RESOLUTION HISTOGRAMS')
        # High Resolution
        pyplot.subplot(2, 1, 1)
        pyplot.hist(x=d_info['high_resolution'], bins=n_bins)
        pyplot.xlabel('HIGH RESOLUTION LIMIT (A)')
        pyplot.ylabel('COUNT')
        # Low Resolution
        pyplot.subplot(2, 1, 2)
        pyplot.hist(x=d_info['low_resolution'], bins=n_bins)
        pyplot.xlabel('LOW RESOLUTION LIMIT (A)')
        pyplot.ylabel('COUNT')
        # Apply tight layout to prevent overlaps
        pyplot.tight_layout()
        # Save both
        pyplot.savefig(os.path.join(img_out_dir, 'd_resolutions.png'))
        pyplot.close(fig)

        # R-FACTORS
        fig = pyplot.figure()
        pyplot.title('R-FACTOR HISTOGRAMS')
        # RFree
        pyplot.subplot(2, 1, 1)
        pyplot.hist(x=d_info['r_free'], bins=n_bins)
        pyplot.xlabel('R-FREE')
        pyplot.ylabel('COUNT')
        # RWork
        pyplot.subplot(2, 1, 2)
        pyplot.hist(x=d_info['r_work'], bins=n_bins)
        pyplot.xlabel('R-WORK')
        pyplot.ylabel('COUNT')
        # Apply tight layout to prevent overlaps
        pyplot.tight_layout()
        # Save both
        pyplot.savefig(os.path.join(img_out_dir, 'd_rfactors.png'))
        pyplot.close(fig)

        # RMSDs
        fig = pyplot.figure()
        pyplot.title('RMSDS TO MEAN STRUCTURE HISTOGRAM')
        pyplot.hist(x=filter_nans(d_info['rmsd_to_mean']), bins=n_bins)
        pyplot.xlabel('RMSD (A)')
        pyplot.ylabel('COUNT')
        # Apply tight layout to prevent overlaps
        pyplot.tight_layout()
        # Save
        pyplot.savefig(os.path.join(img_out_dir, 'd_rmsd_to_mean.png'))
        pyplot.close(fig)

        ########################################################

        self.log('=> Crystal Variation')

        # CELL PARAMS
        fig = pyplot.figure()
        pyplot.title('UNIT CELL PARAMS')
        # A
        pyplot.subplot(2, 3, 1)
        pyplot.hist(x=d_info['uc_a'], bins=n_bins)
        pyplot.xlabel('A (A)')
        pyplot.ylabel('COUNT')
        # B
        pyplot.subplot(2, 3, 2)
        pyplot.hist(x=d_info['uc_b'], bins=n_bins)
        pyplot.xlabel('B (A)')
        # C
        pyplot.subplot(2, 3, 3)
        pyplot.hist(x=d_info['uc_c'], bins=n_bins)
        pyplot.xlabel('C (A)')
        # ALPHA
        pyplot.subplot(2, 3, 4)
        pyplot.hist(x=d_info['uc_alpha'], bins=n_bins)
        pyplot.xlabel('ALPHA')
        pyplot.ylabel('COUNT')
        # BETA
        pyplot.subplot(2, 3, 5)
        pyplot.hist(x=d_info['uc_beta'], bins=n_bins)
        pyplot.xlabel('BETA')
        # GAMMA
        pyplot.subplot(2, 3, 6)
        pyplot.hist(x=d_info['uc_gamma'], bins=n_bins)
        pyplot.xlabel('GAMMA')
        # Apply tight layout to prevent overlaps
        #pyplot.tight_layout()
        # Save both
        pyplot.savefig(os.path.join(img_out_dir, 'd_cell_param.png'))
        pyplot.close(fig)

        # CELL VOLUME
        fig = pyplot.figure()
        pyplot.title('UNIT CELL VOLUME HISTOGRAM')
        pyplot.hist(x=d_info['uc_vol'], bins=n_bins)
        pyplot.xlabel('VOLUME (A**3)')
        pyplot.ylabel('COUNT')
        # Apply tight layout to prevent overlaps
        pyplot.tight_layout()
        # Save
        pyplot.savefig(os.path.join(img_out_dir, 'd_cell_volume.png'))
        pyplot.close(fig)

    def write_output_csvs(self):
        """Write CSV file of dataset variables"""

        self.log('===================================>>>')
        self.log('Writing Dataset + Dataset Map Summary CSV')

        # Write the dataset information to csv file
        self.tables.dataset_info.to_csv(path_or_buf=self.output_handler.get_file('dataset_info'), index_label='dtag')
        self.tables.dataset_map_info.to_csv(path_or_buf=self.output_handler.get_file('dataset_map_info'), index_label='dtag')

        self.log('===================================>>>')
        self.log('Writing COMBINED Dataset Summary CSV')

        # Join the tables on the index of the main table
        comb_tab = self.tables.dataset_info.join(self.tables.dataset_map_info, how='outer')
        comb_tab.to_csv(path_or_buf=self.output_handler.get_file('dataset_combined_info'), index_label='dtag')

        self.log('===================================>>>')
        self.log('Writing Event+Site Summary CSVs')

        # Sort the event data by z-peak and write out
        sort_eve = self.tables.event_info.sort(columns=['site_idx','z_peak'], ascending=[1,0])
        sort_eve.to_csv(path_or_buf=self.output_handler.get_file('event_info_2'))
        sort_eve = sort_eve.join(comb_tab, how='right')
        sort_eve.to_csv(path_or_buf=self.output_handler.get_file('event_info'))
        # Sort the sites by number of events and write out
        sort_sit = self.tables.site_info.sort(columns=['num_events'],ascending=[0])
        sort_sit.to_csv( path_or_buf=self.output_handler.get_file('site_info'))

        # TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
        #self.log('===================================>>>')
        #self.log('Writing Residue Summaries - NOT YET')
        # TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO

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
        return

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
        self.tables.event_info.set_value(event.id, 'est_occupancy', event.info.estimated_occupancy)
        self.tables.event_info.set_value(event.id, 'z_peak', event.cluster.max)
        self.tables.event_info.set_value(event.id, 'z_mean', event.cluster.mean)
        self.tables.event_info.set_value(event.id, 'cluster_size', event.cluster.size)
        self.tables.event_info.set_value(event.id, ['refx','refy','refz'], list(flex.double(event.cluster.peak)*self.reference_grid().grid_spacing()))
        self.tables.event_info.set_value(event.id, ['x','y','z'], list(d_handler.transform_from_reference(  points=flex.vec3_double([event.cluster.peak])*self.reference_grid().grid_spacing(),
                                                                                                        method='global',
#                                                                                                        method=self.params.alignment.method,
#                                                                                                        point_mappings=self.reference_grid().partition().query_by_grid_points([grid_ref])
                                                                                                        )[0]))

    def update_event_table_site_info(self, events):
        """Update the event table for pre-existing events"""
        for e in events:
            assert e.id, 'NO ID GIVEN: {!s}'.format(e.id)
            assert e.parent, 'EVENT HAS NO PARENT: {!s}'.format(e.parent)
            self.tables.event_info.set_value(e.id, 'site_idx', e.parent.id)

    def write_grid_point_distributions(self, grid_points, map_analyser, output_filename=None):
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
            for mh in map_analyser.dataset_maps.all():
                map_vals = [mh.map[g] for g in grid_points]
                # Create list of tags for zipping
                d_tags = [mh.tag]*len(grid_points)
                d_nums = [mh.num]*len(grid_points)
                # Format and write
                out_list = [', '.join(map(str,tup)) for tup in zip(grid_points_str, d_tags, d_nums, map_vals)]
                out_line = '\n'.join(map(str,out_list)) + '\n'
                fh.write(out_line)

#    def write_array_to_map_old(self, output_file, map_data, grid_size=None, grid_spacing=None):
#        """Takes a 1d array and writes it to a map"""
#        if not grid_size:    grid_size    = self.reference_grid().grid_size()
#        if not grid_spacing: grid_spacing = self.reference_grid().grid_spacing()
#        self.log('> Writing Map (OLD METHOD): {!s}'.format(output_file))
#        write_1d_array_as_p1_map(file_name=output_file, map_data=map_data, grid_size=grid_size, grid_spacing=grid_spacing)

    def rotate_map(self, d_handler, map_data):
        """Apply an RT matrix to an array on the reference grid"""

        # For the local alignment transformation
        #point_mappings = self.find_nearest_calpha(points)
        #lab_rt_points = self.local_alignment_transforms()[r_lab].inverse() * flex.vec3_double(lab_points)

        return cctbx.maptbx.rotate_translate_map(   unit_cell          = self.reference_grid().unit_cell(),
                                                    map_data           = map_data,
                                                    rotation_matrix    = d_handler.global_alignment_transform().r.elems,
                                                    translation_vector = d_handler.global_alignment_transform().t.elems    )

    def write_array_to_map(self, output_file, map_data):
        """Take array on the reference grid and write to map"""
        iotbx.ccp4_map.write_ccp4_map(  file_name   = output_file,
                                        unit_cell   = self.reference_grid().unit_cell(),
                                        space_group = self.reference_grid().space_group(),
                                        map_data    = map_data,
                                        labels      = flex.std_string(['Map from PANDDAs'])     )

    def pickle(self, pickle_file, pickle_object, overwrite=True):
        """Takes an object and pickles it"""
        if os.path.exists(pickle_file) and not overwrite:
            self.log('NOT PICKLING: {!s}'.format(pickle_file))
        else:
            self.log('Pickling Object: {!s}'.format(pickle_file))
            easy_pickle.dump(pickle_file, pickle_object)

    def unpickle(self, pickle_file):
        """Takes an object and unpickles it"""
        self.log('Unpickling File: {!s}'.format(pickle_file))
        return easy_pickle.load(pickle_file)


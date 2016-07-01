import os

import iotbx.pdb
import iotbx.mtz
import cctbx.maptbx
import scitbx.math.superpose

from bamboo.common import Meta, Info
from bamboo.common.file import FileManager
from bamboo.common.path import easy_directory

from giant.structure import make_label
from giant.structure.align import perform_flexible_alignment, find_nearest_calphas, transform_coordinates_with_flexible_alignment
from giant.xray.data import CrystalSummary
from giant.xray.symmetry import combine_hierarchies, generate_adjacent_symmetry_copies

########################################################################################################
#
#   DATASET HANDLER CLASSES
#
########################################################################################################

class DatasetHandler(object):
    def __init__(self, dataset_number, pdb_filename, mtz_filename=None, dataset_tag=None):
        """Create a dataset object to allow common functions to be applied easily to a number of datasets"""
        if pdb_filename: assert os.path.exists(pdb_filename), 'PDB file does not exist!'
        if mtz_filename: assert os.path.exists(mtz_filename), 'MTZ file does not exist!'
        # Store dataset number
        self.num = dataset_number
        # Store the tag for the dataset
        if dataset_tag:
            self.tag = dataset_tag
        else:
            # If num < 0 - mark as a reference dataset
            if self.num < 0:  self.tag = 'REF{:05d}'.format(self.num)
            else:             self.tag = 'D{:05d}'.format(self.num)
        # Output Directories
        self.output_handler = None
        self.child = None
        self.meta = Meta()
        ########################################################
        # Store filenames
        self._pdb_file = pdb_filename
        self._mtz_file = mtz_filename
        # PDB Objects
        self._structure = self.new_structure()
        # Data summaries
        if self._pdb_file: self.pdb_summary = None
        if self._mtz_file: self.mtz_summary = CrystalSummary.from_mtz(mtz_file=self._mtz_file)
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
        self.crystal_contact_generators = None
        ########################################################
        # Map of the clusters in the dataset
        self.events = []

    #####################################################################
    #                                                                   #
    #                       UTILITY FUNCTIONS                           #
    #                                                                   #
    #####################################################################

    def initialise_output_directory(self, outputdir):
        """Initialise a dataset output directory"""
        # Create a file and directory organiser
        self.output_handler = FileManager(rootdir=easy_directory(outputdir))

    def get_pickle_copy(self):
        """Get copy of self that can be pickled - some cctbx objects cannot be pickled..."""
        return self

    def pdb_filename(self):
        return self._pdb_file
    def mtz_filename(self):
        return self._mtz_file

    #####################################################################
    #                                                                   #
    #                       HIGH LEVEL OBJECTS                          #
    #                                                                   #
    #####################################################################

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

    #####################################################################
    #                                                                   #
    #                        STRUCTURE STUFF                            #
    #                                                                   #
    #####################################################################

    def heavy_atom_sites(self):
        xray_structure = self.input().xray_structure_simple()
        return xray_structure.sites_cart().select(xray_structure.heavy_selection())
    def calpha_sites(self):
        xray_structure = self.input().xray_structure_simple()
        return xray_structure.sites_cart().select(xray_structure.backbone_selection(atom_names=['CA']))
    def backbone_sites(self):
        xray_structure = self.input().xray_structure_simple()
        return xray_structure.sites_cart().select(xray_structure.backbone_selection())

    def calphas(self):
        """Get the calphas for the structure"""
        return self.hierarchy().select(self.hierarchy().atom_selection_cache().selection('pepnames and name CA'))
    def calpha_labels(self):
        """Return the labels of the calphas of the structure"""
        return [make_label(a) for a in self.calphas().atoms_with_labels()]

    def find_nearest_calpha(self, points, hierarchy=None):
        """Returns the labels of the nearest calpha for each of the given points"""
        if hierarchy is None: hierarchy = self.hierarchy()
        return find_nearest_calphas(hierarchy, coordinates=points)

    #####################################################################
    #                                                                   #
    #                        ALIGNMENT STUFF                            #
    #                                                                   #
    #####################################################################

    def set_global_alignment(self, alignment):
        self._global_rt_transform = alignment
    def global_alignment_transform(self):
        return self._global_rt_transform

    def set_local_alignments(self, alignment):
        self._local_rt_transforms = alignment
    def local_alignment_transforms(self):
        return self._local_rt_transforms

    def transform_coordinates(self, points, method, point_mappings=None, inverse=False):
        """Transform coordinates using contained alignments"""
        assert method in ['local','global'], 'METHOD NOT DEFINED: {!s}'.format(method)
        if method == 'global':
            if inverse: return self.global_alignment_transform().inverse() * points
            else:       return self.global_alignment_transform() * points
        elif method == 'local':
            assert point_mappings is not None
            return transform_coordinates_with_flexible_alignment(   alignments  = self.local_alignment_transforms(),
                                                                    coordinates = points,
                                                                    mappings    = point_mappings,
                                                                    inverse     = inverse)

    def transform_from_reference(self, points, method, point_mappings=None):
        """Use alignment to map to reference frame from our frame"""
        return self.transform_coordinates(  points         = points,
                                            method         = method,
                                            point_mappings = point_mappings,
                                            inverse        = True   )

    def transform_to_reference(self, points, method, point_mappings=None):
        """Use alignment to map to reference frame from our frame"""
        if point_mappings is None:
            point_mappings = self.find_nearest_calpha(points)
        return self.transform_coordinates(  points         = points,
                                            method         = method,
                                            point_mappings = point_mappings,
                                            inverse        = False  )

    #####################################################################
    #                                                                   #
    #                         SYMMETRY STUFF                            #
    #                                                                   #
    #####################################################################

    def generate_symmetry_copies(self, rt_method=None, save_operators=True, buffer=10):
        """Generate the symmetry copies of the reference structure in the reference frame"""

        # Use symmetry operations to create the symmetry mates of the reference structure
        sym_ops, sym_hierarchies, chain_mappings = generate_adjacent_symmetry_copies(   ref_hierarchy    = self.new_structure().hierarchy,
                                                                                        crystal_symmetry = self.input().crystal_symmetry(),
                                                                                        buffer_thickness = buffer    )
        # Record the symmetry operations that generate the crystal contacts
        if save_operators: self.crystal_contact_generators = sym_ops

        # Create a combined hierarchy of the crystal contacts
        symmetry_root = combine_hierarchies(sym_hierarchies)
        # Transform to reference frame?
        if rt_method: symmetry_root.atoms().set_xyz(self.transform_to_reference(points=symmetry_root.atoms().extract_xyz(), method=rt_method))
        # Save coordinates
        self.symmetry_copies = symmetry_root
        return self.symmetry_copies

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

########################################################################################################
#
#   MAP HANDLER CLASSES
#
########################################################################################################

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

########################################################################################################
#
#   HANDLER FUNCTIONS
#
########################################################################################################

def align_dataset_to_reference(d_handler, ref_handler, method):
    """Calculate the rotation and translation needed to align one structure to another"""
    assert method in ['both','local','global'], 'METHOD NOT DEFINED: {!s}'.format(method)
    global_rt_transform = local_rt_transforms = None
    if method == 'global' or method == 'both':
        my_sites = d_handler.calpha_sites()
        ref_sites = ref_handler.calpha_sites()
        assert len(my_sites) == len(ref_sites)
        global_rt_transform = scitbx.math.superpose.least_squares_fit(reference_sites=ref_sites, other_sites=my_sites).rt()
    if method == 'local' or method == 'both':
        local_rt_transforms = perform_flexible_alignment(mov_hierarchy=d_handler.hierarchy(), ref_hierarchy=ref_handler.hierarchy(), cutoff_radius=10)
    return global_rt_transform, local_rt_transforms


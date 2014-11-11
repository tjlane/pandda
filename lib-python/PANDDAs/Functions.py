import os, sys, glob

import iotbx.pdb as pdb_reader
import iotbx.map_tools as map_tools

from iotbx.reflection_file_utils import extract_miller_array_from_file

from Giant.Xray.Structure.Select import get_calpha_sites
from Giant.Xray.Maps.Utils import get_fft_map_from_f_obs_and_structure
from Giant.Xray.Maps.Grid import get_bounding_box_for_structure, create_cartesian_grid


from cctbx import maptbx

class dataset_handler(object):
    def __init__(self, pdb_filename, mtz_filename):
        """Create a reference dataset object as a reference frame for the other datasets to be aligned and scaled to"""

        assert os.path.exists(pdb_filename), 'PDB file does not exist!'
        assert os.path.exists(mtz_filename), 'MTZ file does not exist!'

        # Store filenames
        self._pdb_file = os.path.abspath(pdb_filename)
        self._mtz_file = os.path.abspath(mtz_filename)
        # PDB Objects
        self._pdb_input = pdb_reader.input(source_info=None,lines=open(self._pdb_file,'r').read())
        self._pdb_struc = self._pdb_input.xray_structure_simple()

        # Store C-alpha sites
        self._calpha_sites = get_calpha_sites(input_obj=self._pdb_input, structure_obj=self._pdb_struc)

        # Extract miller array
        self._fobs_miller = extract_miller_array_from_file(self._mtz_file, 'F,SIGF', log=open(os.devnull,'a'))

        # Initialise other variables
        self._fft_map = None
        self._unit_cell = None
        self._space_group = None
        self._basic_map = None

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
    def get_fobs_miller_array(self):
        return self._fobs_miller
    def get_map(self):
        return self._fft_map
    def get_map_handler(self):
        return self._basic_map

    def get_dataset_resolution(self):
        return self.get_fobs_miller_array().d_min()
    def get_structure_min_max_sites(self):
        self._min_max_sites = get_bounding_box_for_structure(self.get_pdb_structure())
        return self._min_max_sites

    def create_fft_map(self, map_type='2mFo-DFc'):
        self._fft_map = get_fft_map_from_f_obs_and_structure(f_obs=self._fobs_miller, xray_structure=self._pdb_struc, map_type=map_type)
        self._unit_cell = self._fft_map.unit_cell()
        self._space_group = self._fft_map.space_group()
        return self._fft_map

    def create_map_handler(self):
        if self._fft_map == None: self.create_fft_map()
        self._basic_map = maptbx.basic_map(maptbx.basic_map_unit_cell_flag(), self._fft_map.real_map(), self._fft_map.real_map().focus(),
                            self._unit_cell.orthogonalization_matrix(), maptbx.out_of_bounds_clamp(0).as_handle(), self._unit_cell)
        return self._basic_map

class grid_handler(object):
    def __init__(self):
        """Create and manage a grid object to be sampled across many aligned datasets"""

        self._grid_size = None
        self._grid_spacing = None

        self._cart_min = None
        self._cart_max = None
        self._cart_size = None
        self._cart_points = None

    def get_grid_size(self):
        return self._grid_size
    def get_grid_spacing(self):
        return self._grid_spacing
    def get_min_cart(self):
        return self._cart_min
    def get_max_cart(self):
        return self._cart_max
    def get_extent_cart(self):
        return (self._cart_min, self._cart_max)
    def get_size_cart(self):
        return self._cart_size

    def get_points_cart(self):
        return self._cart_points
    def get_points_grid(self):
        return flex.nested_loop(self.get_grid_size())

    def set_grid_spacing(self, spacing):
        self._grid_spacing = spacing
    def set_extent_cart(self, cart_min, cart_max):
        self._cart_min = cart_min
        self._cart_max = cart_max

    def create_cartesian_grid(self, include_origin=True):
        if include_origin:
            box_size, self._grid_size, self._cart_points = create_cartesian_grid(min_carts=(0,0,0),
                                                                                 max_carts=self.get_max_cart(),
                                                                                 grid_spacing=self.get_grid_spacing())
        else:
            box_size, self._grid_size, self._cart_points = create_cartesian_grid(min_carts=self.get_min_cart(),
                                                                                 max_carts=self.get_max_cart(),
                                                                                 grid_spacing=self.get_grid_spacing())

        # Update max/min cart sizes as the grid will be slightly larger than the requested size
        self.set_extent_cart(cart_min=self.get_points_cart()[0], cart_max=self.get_points_cart()[-1])

        return self.get_grid_size()

def status_bar(n, n_max):
    if (n+1==n_max):
        print '\r>>', 100, '%'
    elif (n%1000==0):
        print '\r>>', int(round(100.0*(n+1)/n_max,0)), '%',
    sys.stdout.flush()




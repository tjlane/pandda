import os, copy

import numpy

import iotbx.mtz, iotbx.ccp4_map
import cctbx.maptbx, cctbx.uctbx
import cctbx_uctbx_ext
import scitbx.matrix
from scitbx.array_family import flex
from libtbx.utils import Sorry, Failure

from bamboo.common import Meta


class DensityMap(object):


    def __init__(self, map_data, unit_cell, map_indices=None, map_size=None, map_origin=(0.0,0.0,0.0), sparse=False, meta=None, parent=None, children=None):

        assert isinstance(map_data, flex.double)
        assert isinstance(unit_cell, cctbx.uctbx.unit_cell) or isinstance(unit_cell, cctbx_uctbx_ext.unit_cell)

        if sparse:
            assert map_data.nd() == 1, 'Map data must be 1-dimensional when sparse=True'
            assert [map_indices, map_size].count(None) == 0, 'Must provide map_indices and map_size when sparse=True'
            assert len(map_data) == len(map_indices), 'map_data and map_indices must be the same length when sparse=True ({} != {})'.format(len(map_data), len(map_indices))
            assert max(map_indices) < numpy.prod(map_size), 'indices are not compatible with map_size ({} > {})'.format(max(map_indices), numpy.prod(map_size))
            if not isinstance(map_indices, flex.size_t):
                map_indices = flex.size_t(map_indices)
        else:
            if map_size is None:
                assert map_data.nd() == 3, 'map_data must be 3-dimension if map_size is not given'
                map_size = map_data.all()
            assert len(map_size) == 3, 'map_size must be three dimensional'
            assert map_indices is None, 'Do not provide map_indices for non-sparse matrices'
            assert numpy.prod(map_size) == map_data.size()
            # Reshape the map data if necessary
            if map_data.nd() == 1:
                map_data = map_data.deep_copy()
                map_data.reshape(flex.grid(map_size))
            assert map_data.all() == map_size, 'map_data is not the same shape as map_size ({} != {})'.format(map_data.all(), map_size)

        self.data         = map_data
        self._map_size    = map_size
        self._map_indices = map_indices
        self._map_origin  = map_origin
        self.unit_cell = unit_cell
        self.meta      = meta if meta else Meta()
        self.parent    = parent
        self.children  = children if children else []

        assert len(self._map_size)==3, 'map_size must be tuple of length 3'
        assert sparse == self.is_sparse()

    @classmethod
    def from_fft_map(cls, fft_map):
        return cls(map_data=fft_map.real_map(), unit_cell=fft_map.unit_cell())

    def new_from_template(self, map_data, sparse=False, copy_meta=False, same_parent=False):
        """Create a new ElectronDensityMap using this map as a template"""

        # map_data is in sparse form
        if sparse:
            # Make sparse copy if template is not sparse
            if not self.is_sparse():
                self = self.copy().make_sparse()
            # Check input map data is compatible
            assert map_data.nd() == 1, 'map_data must 1-dimensional'
            assert map_data.size() == self._map_indices.size()
            # Extract parameters for sparseness
            map_size    = self._map_size
            map_indices = self._map_indices
        else:
            assert map_data.size() == numpy.prod(self._map_size)
            map_size    = self._map_size
            map_indices = None

        if copy_meta:   meta = copy.deepcopy(self.meta)
        else:           meta = None
        if same_parent: parent = self.parent
        else:           parent = None

        return ElectronDensityMap(map_data      = map_data,
                                  unit_cell     = self.unit_cell,
                                  map_indices   = map_indices,
                                  map_size      = map_size,
                                  map_origin    = self._map_origin,
                                  meta          = meta,
                                  parent        = parent,
                                  sparse        = sparse)

    def copy(self):
        return self.new_from_template(map_data    = self.data.deep_copy(),
                                      sparse      = self.is_sparse(),
                                      copy_meta   = True,
                                      same_parent = True)

    def normalised_copy(self):
        """Perform rms scaling on map data"""

        # Create output copy and make sparse (always calculate the mean and rms from the sparse values)
        result = self.copy().make_sparse()
        map_data = result.get_map_data(sparse=True)
        # Apply normalisation
        result.data = (result.data-numpy.mean(map_data)) * (1.0/numpy.std(map_data))

        # Return the modified map
        if self.is_sparse():
            return result.make_sparse()
        else:
            return result.make_dense()

    def _check_compatibility(self, other):
        assert self.is_sparse() is other.is_sparse()

    def __add__(self, other):
        if isinstance(other, ElectronDensityMap):
            self._check_compatibility(other=other)
            return self.__add__(other.data)
        else:
            return self.new_from_template(map_data=self.data+other, sparse=self.is_sparse())

    def __sub__(self, other):
        if isinstance(other, ElectronDensityMap):
            self._check_compatibility(other=other)
            return self.__sub__(other.data)
        else:
            return self.new_from_template(map_data=self.data-other, sparse=self.is_sparse())

    def __mul__(self, other):
        if isinstance(other, ElectronDensityMap):
            self._check_compatibility(other=other)
            return self.__mul__(other.data)
        else:
            return self.new_from_template(map_data=self.data*other, sparse=self.is_sparse())

    def __div__(self, other):
        if isinstance(other, ElectronDensityMap):
            self._check_compatibility(other=other)
            return self.__div__(other.data)
        else:
            return self.new_from_template(map_data=self.data*(1.0/other), sparse=self.is_sparse())

    def __rdiv__(self, other):
        return self.__div__(other)

    def is_sparse(self):
        return (self._map_indices is not None)

    def embed(self, map_data):
        """Embed map data relative to the real map origin, rather than (0,0,0)"""
        if self._map_origin == (0.0,0.0,0.0): return map_data
        return cctbx.maptbx.rotate_translate_map(unit_cell          = self.unit_cell,
                                                 map_data           = map_data,
                                                 rotation_matrix    = scitbx.matrix.rec([1,0,0,0,1,0,0,0,1], (3,3)).elems,
                                                 translation_vector = (-1.0*scitbx.matrix.rec(self._map_origin, (3,1))).elems    )

#    def as_map(self):
#        map_data = self.get_map_data(sparse=False)
#        return cctbx.maptbx.basic_map(
#                    cctbx.maptbx.basic_map_unit_cell_flag(),
#                    self.embed(map_data),
#                    map_data.focus(),
#                    self.unit_cell.orthogonalization_matrix(),
#                    cctbx.maptbx.out_of_bounds_clamp(0).as_handle(),
#                    self.unit_cell)

    def get_cart_values(self, cart_points):
        assert not self.is_sparse(), 'map must not be in sparse format for sampling'
        # Shift input points to the grid frame -- TODO implement in function so that rotations can be automated integrated
        cart_points = (cart_points - self._map_origin)
        frac_values = self.unit_cell.fractionalize(cart_points)
        # Get the map data with the correct origin
        map_data = self.get_map_data(sparse=False)
        map_vals = map(map_data.eight_point_interpolation, frac_values)
        return flex.double(map_vals)

    def to_file(self, filename, space_group):
        map_data = self.get_map_data(sparse=False)
        iotbx.ccp4_map.write_ccp4_map(
                    file_name   = filename,
                    unit_cell   = self.unit_cell,
                    space_group = space_group,
                    map_data    = self.embed(map_data),
                    labels      = flex.std_string(['Output map from giant/pandda'])     )

    def get_map_data(self, sparse):
        """Get the map data as sparse/dense without altering state of master object"""
        if sparse is not self.is_sparse():
            result = self.copy()
            if sparse:
                result.make_sparse()
            else:
                result.make_dense()
        else:
            result = self

        return result.data

    def make_sparse(self):
        """Convert the map data into sparse form"""
        if self.is_sparse(): return self

        data_flat = self.data.as_1d()
        data_mask = (data_flat != 0.0)
        sparse_idxs = data_mask.iselection()
        sparse_data = data_flat.select(data_mask)

        self.data = sparse_data
        self._map_indices = sparse_idxs

        return self

    def make_dense(self):
        """Convert the map data into dense form"""
        if not self.is_sparse(): return self

        data_bulk = numpy.zeros(numpy.product(self._map_size))
        data_bulk.put(indices=self._map_indices, values=self.data)
        data_bulk = flex.double(data_bulk)
        data_bulk.reshape(flex.grid(self._map_size))

        self.data = data_bulk
        self._map_indices = None

        return self


# For backwards compativility : TODO Delete
ElectronDensityMap = DensityMap

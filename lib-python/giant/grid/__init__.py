import os, sys, glob, time, re

from scipy import spatial
import numpy

from libtbx import easy_mp
from libtbx.math_utils import ifloor, iceil
from scitbx.array_family import flex
import iotbx.crystal_symmetry_from_any

from cctbx import crystal

from giant.grid.utils import calculate_grid_size


class Grid(object):


    def __init__(self, grid_spacing, origin, approx_max, verbose=True):
        """Create and manage a grid object"""

        self.verbose = verbose

        self._origin = tuple(origin)
        self._grid_spacing = grid_spacing
        self._grid_size = calculate_grid_size(min_carts=origin, max_carts=approx_max, grid_spacing=grid_spacing)

        self._grid_masks = {}
        self.partition = None

    def grid_spacing(self):
        return self._grid_spacing
    def grid_size(self):
        return self._grid_size
    def grid_size_1d(self):
        return numpy.product(self._grid_size)
    def grid_point_volume(self):
        return self.grid_spacing()**3

    def cart_size(self):
        return tuple([x*self._grid_spacing for x in self.grid_size()])
    def cart_origin(self):
        return self._origin

    def cart_points(self, origin=True):
        assert isinstance(origin, bool)
        return flex.vec3_double(flex.nested_loop(grid_size))*self.grid_spacing() + self.cart_origin()*origin
    def grid_points(self):
        return flex.nested_loop(self.grid_size())

    def grid2cart(self, sites_grid, origin=True):
        assert isinstance(origin, bool)
        if origin: shift = self.cart_origin()
        else:      shift = (0,0,0)
        return flex.vec3_double(sites_grid)*self.grid_spacing() + shift

    def indexer(self):
        """Translates between 3d grid coordinates and 1d array coordinates"""
        return flex.grid(self.grid_size())

    def unit_cell(self):
        """Create a unit cell as if the reference grid were a real lattice (THE UNIT CELL IS LARGER THAN cart_size() BY 1 GRID SPACING)"""
        return crystal.uctbx.unit_cell('{!s} {!s} {!s} 90 90 90'.format(*list(flex.double(self.grid_size())*self.grid_spacing())))
    def space_group(self):
        """Create a spacegroup as if the reference grid were a real lattice (Cartesian Grid == P1 SpaceGroup)"""
        return crystal.sgtbx.space_group('P1')
    def crystal_symmetry(self):
        """Create a symmetry object as if the reference grid were a real lattice (THE UNIT CELL IS LARGER THAN cart_size() BY 1 GRID SPACING)"""
        return iotbx.crystal_symmetry_from_any.from_string("{:f},{:f},{:f},90,90,90,P1".format(*list(flex.double(self.grid_size())*self.grid_spacing())))

    def create_grid_partition(self, sites_cart):
        """Partition the grid using nearest neighbour algorithm"""
        self.partition = GridPartition(grid=self, sites_cart=sites_cart)
        return self.partition

    def add_mask(self, mask_name, mask):
        """Add a an atomic mask to the reference grid"""
        assert mask_name not in self._grid_masks.keys(), 'MASK ALREADY ADDED: {!s}'.format(mask_name)
        self._grid_masks[mask_name] = mask
    def get_mask(self, mask_name):
        """Return a named atomic mask"""
        return self._grid_masks[mask_name]

    def set_global_mask(self, mask):
        """Add a global mask to the grid object - This will create binary mask of the masked grid points"""
        self._grid_masks['protein'] = mask
        print self.global_mask().summary()
    def global_mask(self):
        return self._grid_masks.get('protein', None)

    def set_symmetry_mask(self, mask):
        """Add a global mask to the grid object - This will create binary mask of the masked grid points"""
        self._grid_masks['symmetry'] = mask
        print self.symmetry_mask().summary()
    def symmetry_mask(self):
        return self._grid_masks.get('symmetry', None)

    def summary(self):
        return '\n'.join([  '----------------------------------->>>',
                            'Reference Grid Summary:',
                            'Grid Spacing:        {!s}'.format(round(self.grid_spacing(), 3)),
                            'Grid Point Volume:   {!s}'.format(round(self.grid_point_volume(),3)),
                            'Size of Grid (3D):   {!s}'.format(self.grid_size()),
                            'Size of Grid (1D):   {!s}'.format(self.grid_size_1d()),
                            'Grid Origin  (Cart):  {!s}'.format(tuple([round(x,3) for x in self.cart_origin()])),
                            'Size of Grid (Cart): {!s}'.format(tuple([round(x,3) for x in self.cart_size()]))
                        ])


class GridPartition(object):


    def __init__(self, grid, sites_cart):
        """Partition a grid based on the nearest neighbour calpha atom for each grid site"""
        assert isinstance(grid, Grid), 'grid must be of type Grid'
        self.parent=grid
        self.sites_cart=sites_cart
        self.sites_grid=(sites_cart-self.parent.cart_origin())*(1.0/self.parent.grid_spacing())

    def partition(self, mask=None, cpus=1):
        """Find the nearest neighbour for each grid point (or the subset defined by mask.outer_mask() if mask is not None)"""

        def find_sites(sites_tuple):
            ref_sites, query_sites = sites_tuple
            tree = spatial.KDTree(data=ref_sites)
            nn_dists, nn_groups = tree.query(query_sites)
            return nn_groups

        assert isinstance(cpus, int) and (cpus > 0)

        # Sites that we are partitioning
        if mask: query_sites = flex.vec3_double(mask.outer_mask())
        else:    query_sites = flex.vec3_double(self.parent.grid_points())
        # Find the nearest grid_site for each query_site (returns index of the grid site)
        if cpus == 1:
            output = [find_sites((self.sites_grid, query_sites))]
        else:
            # Chunk the points into groups
            chunk_size = iceil(1.0*len(query_sites)/cpus)
            chunked_points = [query_sites[i:i + chunk_size] for i in range(0, len(query_sites), chunk_size)]
            assert sum([len(a) for a in chunked_points]) == len(query_sites)
            assert len(chunked_points) == cpus
            # Map to cpus
            arg_list = [(self.sites_grid, chunk) for chunk in chunked_points]
            output = easy_mp.pool_map(fixed_func=find_sites, args=arg_list, processes=cpus)

        assert len(output) == cpus, '{!s} != {!s}'.format(len(output), cpus)
        # Extract the indices of the mapped points
        nn_groups = []; [nn_groups.extend(o) for o in output]
        nn_groups = numpy.array(nn_groups)
        assert len(query_sites) == len(nn_groups)
        # Reformat into full grid size
        if mask:
            self.nn_groups = -1*numpy.ones(self.parent.grid_size_1d(), dtype=int)
            self.nn_groups.put(mask.outer_mask_indices(), nn_groups)
        else:
            self.nn_groups = nn_groups

        return self

    def query_by_grid_indices(self, idxs):
        """Return the atom label for a grid site index"""
        assert self.nn_groups is not None, 'Grid not yet partitioned'
        return numpy.array([self.nn_groups[i] for i in idxs])

    def query_by_grid_points(self, gps):
        """Return the atom label for a grid point"""
        assert self.nn_groups is not None, 'Grid not yet partitioned'
        return numpy.array([self.nn_groups[self.grid_indexer(g)] for g in gps])

    def query_by_cart_points(self, sites_cart):
        """Dynamically calculate the nearest atom site to the input points"""
        tree = spatial.KDTree(data=self.sites)
        nn_dists, nn_groups = tree.query(sites_cart)
        return numpy.array(nn_groups)

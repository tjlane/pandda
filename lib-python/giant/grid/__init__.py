import os, sys, glob, time, re

from scipy import spatial
import numpy

from libtbx import easy_mp
from libtbx.math_utils import ifloor, iceil
from scitbx.array_family import flex
import iotbx.crystal_symmetry_from_any

from cctbx import crystal

from giant.structure import make_label
from giant.grid.utils import create_cartesian_grid

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

        # Atomic Masks for partitioning the grid
        self._grid_masks = {}
        # Local mask for filtering/smoothing
        self._local_mask = None

        # Partiton object for spliting the grid into sections around atoms
        self._grid_partition = None

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
        return flex.vec3_double(flex.nested_loop(grid_size)) * self.grid_spacing() + self.cart_min()
    def grid_points(self):
        return flex.nested_loop(self.grid_size())

    def grid_indexer(self):
        """Translates between 3d grid coordinates and 1d array coordinates"""
        return flex.grid(self.grid_size())

    def unit_cell(self):
        """Create a unit cell as if the reference grid were a real lattice (THE UNIT CELL IS LARGER THAN cart_size() BY 1 GRID SPACING)"""
        return crystal.uctbx.unit_cell('{!s} {!s} {!s} 90 90 90'.format(*list(flex.double(self.grid_size())*self.grid_spacing())))
    def space_group(self):
        """Create a spacegroup as if the reference grid were a real lattice (Cartesian Grid == P1 SpaceGroup)"""
        return crystal.sgtbx.space_group('P1')
    def symmetry(self):
        """Create a symmetry object as if the reference grid were a real lattice (THE UNIT CELL IS LARGER THAN cart_size() BY 1 GRID SPACING)"""
        return iotbx.crystal_symmetry_from_any.from_string("{:f},{:f},{:f},90,90,90,P1".format(*list(flex.double(self.grid_size())*self.grid_spacing())))

    def create_grid_partition(self, atomic_hierarchy):
        """Partition the grid using nearest neighbour algorithm"""
        self._grid_partition = grid_partition(  grid_size = self.grid_size(),
                                                grid_spacing = self.grid_spacing(),
                                                atomic_hierarchy = atomic_hierarchy  )
    def partition(self):
        return self._grid_partition

    def add_mask(self, mask_name, mask):
        """Add a an atomic mask to the reference grid"""
        assert mask_name not in self._grid_masks.keys(), 'MASK ALREADY ADDED: {!s}'.format(mask_name)
        self._grid_masks[mask_name] = mask
    def mask(self, mask_name):
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

    def set_local_mask(self, mask):
        """Add local mask to the grid object"""

        self._local_mask = mask
        print self.local_mask().summary()

#        grid_size = self.grid_size()
#        grid_jump = self.local_mask().grid_jump()
#        grid_buffer = self.local_mask().buffer_size()
#        grid_indexer = self.grid_indexer()
#
#        # Resample grid points
#        self._resampled_grid_points = [gp for gp in flex.nested_loop(grid_size) if not [1 for coord in gp if (coord%grid_jump != 0)]]
#        self._resampled_grid_size = tuple([int(1+(g-1)/grid_jump) for g in grid_size])
#        self._resampled_grid_spacing = grid_jump*self.grid_spacing()
#
#        # Create a buffer zone at the edge of the grid
#        buffer_mask_points = [gp for gp in flex.nested_loop(grid_size) if [1 for i_dim, coord in enumerate(gp) if (coord<grid_buffer or coord>(grid_size[i_dim]-1-grid_buffer))]]
#
#        # Create binary mask for the buffer zone
#        buffer_mask_indices = []
#        [buffer_mask_indices.append(grid_indexer(gp)) for gp in buffer_mask_points]
#
#        self._buffer_mask_indices = numpy.array(buffer_mask_indices)

    def local_mask(self):
        return self._local_mask

    def summary(self):
        return '\n'.join([  '----------------------------------->>>',
                            'Reference Grid Summary:',
                            'Grid Spacing:        {!s}'.format(round(self.grid_spacing(), 3)),
                            'Grid Point Volume:   {!s}'.format(round(self.grid_point_volume(),3)),
                            'Size of Grid (3D):   {!s}'.format(self.grid_size()),
                            'Size of Grid (1D):   {!s}'.format(self.grid_size_1d()),
                            'Min of Grid (Cart): {!s}'.format(tuple([round(x,3) for x in self.cart_min()])),
                            'Max of Grid (Cart): {!s}'.format(tuple([round(x,3) for x in self.cart_max()])),
                            'Size of Grid (Cart): {!s}'.format(tuple([round(x,3) for x in self.cart_size()]))
                        ])

    def create_cartesian_grid(self, expand_to_origin):
        if expand_to_origin:
            # TODO Don't know if this will still work - raise error for now
            raise Exception('NOT CURRENTLY CHECKED')
            assert [i>0 for i in self.cart_min()], 'ALL GRID SITES MUST BE GREATER THAN 0 IF ORIGIN INCLUDED'
            assert [i>0 for i in self.cart_max()], 'ALL GRID SITES MUST BE GREATER THAN 0 IF ORIGIN INCLUDED'
            box_size, self._grid_size, cart_points = create_cartesian_grid(min_carts=(0,0,0),
                                                                                 max_carts=self.cart_max(),
                                                                                 grid_spacing=self.grid_spacing())
        else:
            box_size, self._grid_size, cart_points = create_cartesian_grid(min_carts=self.cart_min(),
                                                                                 max_carts=self.cart_max(),
                                                                                 grid_spacing=self.grid_spacing())

        # Update max/min cart sizes as the grid will be slightly larger than the requested size
        self.set_cart_extent(cart_min=cart_points.min(), cart_max=cart_points.max())

        return self.grid_size()

class grid_partition(object):
    def __init__(self, grid_size, grid_spacing, atomic_hierarchy):
        """Partition a grid based on the nearest neighbour calpha atom for each grid site"""

        # Save inputs
        self.hierarchy = atomic_hierarchy
        self.grid_size = grid_size
        self.grid_spacing = grid_spacing
        self.grid_indexer = flex.grid(grid_size)

        # Calculate partition variables
        atoms = [at for at in atomic_hierarchy.atoms_with_labels()]
        #self.grid_sites = list(flex.nested_loop(grid_size))
        self.atom_sites_grid = numpy.array([a.xyz for a in atoms])/grid_spacing
        self.atom_sites_cart = numpy.array([a.xyz for a in atoms])
        self.atom_labels = [make_label(a) for a in atoms]

        self._atoms = atomic_hierarchy.deep_copy()

        # Distances from grid sites to nearest atom site
        self.nn_dists = None
        # Index of nearest atom to grid sites
        self.nn_groups = None
        # Label of nearest atom to grid sites
        #self.nn_atom_labels = None

    def grid_sites(self):
        return list(flex.nested_loop(self.grid_size))

    def partition_grid(self, cpus=1):
        """Find the nearest neighbour for each grid point"""

#        tree = spatial.KDTree(data=self.atom_sites_grid)
#        self.nn_dists, self.nn_groups = tree.query(self.grid_sites())
#        self.nn_atom_labels = [self.atom_labels[i] for i in self.nn_groups]

        def find_sites(sites_dict):
            ref_sites   = sites_dict['ref']
            query_sites = sites_dict['query']
            tree = spatial.KDTree(data=ref_sites)
            nn_dists, nn_groups = tree.query(query_sites)
            return {'nn_dists':nn_dists, 'nn_groups':nn_groups}

        assert isinstance(cpus, int)
        assert cpus > 0

        # Points that define neighbourhoods
        ref_sites = self.atom_sites_grid
        # Sites that we are partitioning
        query_sites = self.grid_sites()

        if cpus == 1:
            output = [find_sites({'ref':ref_sites, 'query':query_sites})]
        else:
            # Chunk the points into groups
            chunk_size = iceil(1.0*len(query_sites)/cpus)
            chunked_points = [query_sites[i:i + chunk_size] for i in range(0, len(query_sites), chunk_size)]
            assert sum([len(a) for a in chunked_points]) == len(query_sites)
            assert len(chunked_points) == cpus

            # Map to cpus
            arg_list = [{'ref':ref_sites, 'query':chunk} for chunk in chunked_points]
            output = easy_mp.pool_map(fixed_func=find_sites, args=arg_list, processes=cpus)

        assert len(output) == cpus, '{!s} != {!s}'.format(len(output), cpus)

        self.nn_dists = [];  [self.nn_dists.extend(t['nn_dists']) for t in output]
        self.nn_dists = numpy.array(self.nn_dists)
        assert len(query_sites) == len(self.nn_dists)
        self.nn_groups = []; [self.nn_groups.extend(t['nn_groups']) for t in output]
        self.nn_groups = numpy.array(self.nn_groups)
        assert len(query_sites) == len(self.nn_groups)
        #self.nn_atom_labels = [self.atom_labels[i] for i in self.nn_groups]

    def query_by_grid_indices(self, idxs):
        """Return the atom label for a grid site index"""
        if self.nn_groups is None: self.partition_grid()
        return [self.atom_labels[self.nn_groups[i]] for i in idxs]

    def query_by_grid_points(self, gps):
        """Return the atom label for a grid point"""
        if self.nn_groups is None: self.partition_grid()
        return [self.atom_labels[self.nn_groups[self.grid_indexer(g)]] for g in gps]

    def query_by_cart_points(self, sites_cart):
        """Dynamically calculate the nearest atom site to the input points"""
        tree = spatial.KDTree(data=self.atom_sites_cart)
        nn_dists, nn_groups = tree.query(sites_cart)
        return [self.atom_labels[i] for i in nn_groups]

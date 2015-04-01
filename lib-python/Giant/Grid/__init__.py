import os, sys, glob, time, re

import numpy

from scitbx.array_family import flex

from cctbx import crystal



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

        # Atomic Masks for partitioning the grid
        self._grid_masks = {}

        # Local mask for filtering/smoothing
        self._local_mask = None

#        # Groups of points defined by the local mask
#        self._buffer_mask_points = None
#        self._buffer_mask_binary = None
#
#        # Group of points defined by sampling the grid regularly
#        self._resampled_grid_points = None
#        self._resampled_grid_size = None

        # Manually set group of points created by combining masks
        self._masked_grid_points = None

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

    def create_grid_partition(self, atomic_hierarchy):
        """Partition the grid using nearest neighbour algorithm"""
        self._grid_partition = grid_partition(  grid_size = self.grid_size(),
                                                grid_spacing = self.grid_spacing(),
                                                atomic_hierarchy = atomic_hierarchy  )
    def partition(self):
        if not self._grid_partition: raise Exception('GRID NOT PARTITIONED')
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

#    def resampled_grid_points(self):
#        """Get a down-sampled list of grid points, based on the local mask used"""
#        return self._resampled_grid_points
#    def resampled_grid_size(self):
#        """Gets the size of the re-sampled grid"""
#        return self._resampled_grid_size
#    def resampled_grid_spacing(self):
#        """Gets the grid spacing for the re-sampled grid"""
#        return self._resampled_grid_spacing
#    def resampled_grid_indexer(self):
#        """Translates between 3d grid coordinates and 1d array coordinates"""
#        return flex.grid(self.resampled_grid_size())
#
#    def buffer_mask_points(self):
#        """Get a list of points in the buffer zone of the map"""
#        return self._buffer_mask_points
#    def buffer_mask_binary(self):
#        """Get a binary mask for the buffer zone"""
#        return self._buffer_mask_binary
#
#    def set_masked_grid_points(self, masked_points):
#        self._masked_grid_points = masked_points
#    def masked_grid_points(self):
#        return self._masked_grid_points

    def summary(self):
        return '\n'.join([  '===================================>>>',
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
            box_size, self._grid_size, self._cart_points = create_cartesian_grid(min_carts=(0,0,0),
                                                                                 max_carts=self.cart_max(),
                                                                                 grid_spacing=self.grid_spacing())
        else:
            box_size, self._grid_size, self._cart_points = create_cartesian_grid(min_carts=self.cart_min(),
                                                                                 max_carts=self.cart_max(),
                                                                                 grid_spacing=self.grid_spacing())

        # Update max/min cart sizes as the grid will be slightly larger than the requested size
        self.set_cart_extent(cart_min=self.cart_points().min(), cart_max=self.cart_points().max())

        return self.grid_size()


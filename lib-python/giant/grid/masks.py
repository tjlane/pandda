import os, sys, glob, re
import time

import numpy

from cctbx import maptbx
from cctbx import crystal

from scitbx.array_family import flex

from giant.grid.utils import get_grid_points_within_distance_cutoff_of_cart_sites_2
from giant.grid.utils import get_grid_points_within_distance_cutoff_of_origin, combine_grid_point_and_grid_vectors

class atomic_mask(object):
    def __init__(self, cart_sites, grid_size, unit_cell, max_dist, min_dist):
        """Take a grid and calculate all grid points with a certain distance cutoff of any point in cart_sites"""

        if min_dist: assert max_dist > min_dist, 'Minimum Mask Distance must be smaller than Maximum Mask Distance'

        # Store distances from masking atoms
        self._max_dist = max_dist
        self._min_dist = min_dist

        # Store grid size
        self._grid_size = grid_size
        self._grid_idxr = flex.grid(grid_size)

        # Unit cell
        self._fake_unit_cell = unit_cell

        # Calculate the masked indices defined by max distance from protein atoms
        self._outer_mask_indices = maptbx.grid_indices_around_sites(unit_cell  = unit_cell,
                                                                    fft_n_real = grid_size, fft_m_real=grid_size,
                                                                    sites_cart = cart_sites,
                                                                    site_radii = flex.double(cart_sites.size(), max_dist))
        outer_mask_binary = numpy.zeros(self._grid_idxr.size_1d(), dtype=bool)
        outer_mask_binary.put(self._outer_mask_indices, True)
        self._outer_mask_binary = outer_mask_binary

        # Calculate the masked indices defined by min distance from protein atoms
        self._inner_mask_indices = maptbx.grid_indices_around_sites(unit_cell  = unit_cell,
                                                                    fft_n_real = grid_size, fft_m_real = grid_size,
                                                                    sites_cart = cart_sites,
                                                                    site_radii = flex.double(cart_sites.size(), min_dist))
        inner_mask_binary = numpy.zeros(self._grid_idxr.size_1d(), dtype=bool)
        inner_mask_binary.put(self._inner_mask_indices, True)
        self._inner_mask_binary = inner_mask_binary

        # Calculate the combination of these masks
        total_mask_binary = numpy.zeros(self._grid_idxr.size_1d(), bool)
        total_mask_binary.put(self._outer_mask_indices, True)
        total_mask_binary.put(self._inner_mask_indices, False)
        self._total_mask_binary = total_mask_binary
        self._total_mask_indices = [idx for idx in xrange(self._grid_idxr.size_1d()) if self._total_mask_binary[idx]>0]

    def total_mask(self):
        """Return the grid points allowed by the mask - combination of max_dist (allowed) and min_dist (rejected)"""
        for p in flex.nested_loop(self._grid_size):
            if self._total_mask_binary[self._grid_idxr(p)]: yield p
    def outer_mask(self):
        """Get grid points allowed subject to max_dist"""
        for p in flex.nested_loop(self._grid_size):
            if self._outer_mask_binary[self._grid_idxr(p)]: yield p
    def inner_mask(self):
        """Get grid points rejected subject to min_dist"""
        for p in flex.nested_loop(self._grid_size):
            if self._inner_mask_binary[self._grid_idxr(p)]: yield p

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
        return len(self._total_mask_indices)
    def outer_size(self):
        """Returns the number of grid points inside max_dist"""
        return len(self._outer_mask_indices)
    def inner_size(self):
        """Returns the number of grid points inside min_dist"""
        return len(self._inner_mask_indices)

    def extent(self):
        """Returns the minimum and maximum grid points in the mask"""
        return min(self.total_mask()), max(self.total_mask())

    def summary(self):
        return '\n'.join([  '----------------------------------->>>',
                            'Atomic Mask Summary:',
                            'Total Mask Size (1D): {!s}'.format(self.total_size()),
                            'Outer Mask Size (1D): {!s}'.format(self.outer_size()),
                            'Inner Mask Size (1D): {!s}'.format(self.inner_size()),
                            'Masked Grid Min/Max: {!s}'.format(self.extent())
                        ])

class grid_mask(atomic_mask):
    """Creates the same mask as `atomic_mask` without allowing for lattice symmetry (rejects atoms outside of the unit cell)"""

    def __init__(self, cart_sites, grid_size, unit_cell, max_dist, min_dist):
        """Mask a grid (approximately) against a set of cartesian points"""

        if min_dist: assert max_dist >= min_dist, 'Minimum Mask Distance must be smaller than Maximum Mask Distance'

        # Store distances from masking atoms
        self._max_dist = max_dist
        self._min_dist = min_dist

        # Store grid size
        self._grid_size = grid_size
        self._grid_idxr = flex.grid(grid_size)

        # Unit cell
        self._fake_unit_cell = unit_cell

        # Filter out sites that are outside of the unit cell
        frac_sites = unit_cell.fractionalize(cart_sites)
        ax1_len_1 = (frac_sites.dot((1.0,0.0,0.0)) <  1.000001).iselection()
        ax2_len_1 = (frac_sites.dot((0.0,1.0,0.0)) <  1.000001).iselection()
        ax3_len_1 = (frac_sites.dot((0.0,0.0,1.0)) <  1.000001).iselection()
        ax1_len_2 = (frac_sites.dot((1.0,0.0,0.0)) > -0.000001).iselection()
        ax2_len_2 = (frac_sites.dot((0.0,1.0,0.0)) > -0.000001).iselection()
        ax3_len_2 = (frac_sites.dot((0.0,0.0,1.0)) > -0.000001).iselection()
        allowed_idxs_1 = ax1_len_1.intersection(ax2_len_1).intersection(ax3_len_1)
        allowed_idxs_2 = ax1_len_2.intersection(ax2_len_2).intersection(ax3_len_2)
        allowed_idxs = allowed_idxs_1.intersection(allowed_idxs_2)
        if not allowed_idxs: raise Exception('No Sites valid for masking!')
        filt_cart_sites = cart_sites.select(allowed_idxs)
#        filt_cart_sites_2 = flex.vec3_double([p for p in cart_sites if not [x for x in unit_cell.fractionalize(p) if (x>1.000001 or x<-0.000001)]])

        # Calculate the masked indices defined by max distance from protein atoms
        self._outer_mask_indices = maptbx.grid_indices_around_sites(unit_cell  = unit_cell,
                                                                    fft_n_real = grid_size, fft_m_real = grid_size,
                                                                    sites_cart = filt_cart_sites,
                                                                    site_radii = flex.double(filt_cart_sites.size(), max_dist))
        outer_mask_binary = numpy.zeros(self._grid_idxr.size_1d(), dtype=bool)
        outer_mask_binary.put(self._outer_mask_indices, True)
        self._outer_mask_binary = outer_mask_binary

        # Calculate the masked indices defined by min distance from protein atoms
        self._inner_mask_indices = maptbx.grid_indices_around_sites(unit_cell  = unit_cell,
                                                                    fft_n_real = grid_size, fft_m_real = grid_size,
                                                                    sites_cart = filt_cart_sites,
                                                                    site_radii = flex.double(filt_cart_sites.size(), min_dist))
        inner_mask_binary = numpy.zeros(self._grid_idxr.size_1d(), dtype=bool)
        inner_mask_binary.put(self._inner_mask_indices, True)
        self._inner_mask_binary = inner_mask_binary

        # Calculate the combination of these masks
        total_mask_binary = numpy.zeros(self._grid_idxr.size_1d(), bool)
        total_mask_binary.put(self._outer_mask_indices, True)
        total_mask_binary.put(self._inner_mask_indices, False)
        self._total_mask_binary  = total_mask_binary
        self._total_mask_indices = [idx for idx in xrange(self._grid_idxr.size_1d()) if self._total_mask_binary[idx]>0]

class non_symmetrical_atomic_mask(atomic_mask):
    """Creates the same mask as `atomic_mask` without allowing for unit cell shifts from the p1 box"""

    def __init__(self, cart_sites, grid_spacing, grid_size, unit_cell, max_dist, min_dist):
        """Take a grid and calculate all grid points with a certain distance cutoff of any point in cart_sites"""

        if min_dist: assert max_dist >= min_dist, 'Minimum Mask Distance must be smaller than Maximum Mask Distance'

        # Store distances from masking atoms
        self._max_dist = max_dist
        self._min_dist = min_dist

        # Store grid size
        self._grid_size = grid_size
        self._grid_idxr = flex.grid(grid_size)

        # Unit cell
        self._fake_unit_cell = unit_cell

        # First, trim off cartesian sites that won't contribute to the masking (those more than max_dist from the edge of the box)
        # Extract the maximum and minimum of the unit cell (assumed to be a p1 grid on the origin)
        uc_min = flex.double([0.0,0.0,0.0])
        uc_max = flex.double(unit_cell.parameters()[:3])
        # We'll be looking for points within max_dist of this grid
        grid_filter_minimum = uc_min - max_dist
        grid_filter_maximum = uc_max + max_dist
        # Filter out sites that are too far from the grid to be worth masking
        filt_cart_sites = flex.vec3_double([p for p in cart_sites if not (    ([i+1 for i in range(3) if p[i]<grid_filter_minimum[i]])
                                                                           or ([i+1 for i in range(3) if p[i]>grid_filter_maximum[i]])    )])

        # Mask the grid
        total_mask_points, \
        outer_mask_points, \
        inner_mask_points = get_grid_points_within_distance_cutoff_of_cart_sites_2( cart_sites=filt_cart_sites,
                                                                                    grid_spacing=grid_spacing,
                                                                                    max_dist=max_dist,
                                                                                    min_dist=min_dist   )

        # Remove points that are outside of the grid we are interested in
        total_mask_points = [gp for gp in total_mask_points if not ( ([i+1 for i in range(3) if gp[i]<=0]) or ([i+1 for i in range(3) if gp[i]>=grid_size[i]]) )]
        outer_mask_points = [gp for gp in outer_mask_points if not ( ([i+1 for i in range(3) if gp[i]<=0]) or ([i+1 for i in range(3) if gp[i]>=grid_size[i]]) )]
        inner_mask_points = [gp for gp in inner_mask_points if not ( ([i+1 for i in range(3) if gp[i]<=0]) or ([i+1 for i in range(3) if gp[i]>=grid_size[i]]) )]

        # Convert the list of points to lists of indices, and binary masks
        self._total_mask_indices = [self._grid_idxr(gp) for gp in total_mask_points]
        total_mask_binary = numpy.zeros(self._grid_idxr.size_1d(), dtype=bool)
        total_mask_binary.put(self._total_mask_indices, True)
        self._total_mask_binary = total_mask_binary

        self._outer_mask_indices = [self._grid_idxr(gp) for gp in outer_mask_points]
        outer_mask_binary = numpy.zeros(self._grid_idxr.size_1d(), dtype=bool)
        outer_mask_binary.put(self._outer_mask_indices, True)
        self._outer_mask_binary = outer_mask_binary

        self._inner_mask_indices = [self._grid_idxr(gp) for gp in inner_mask_points]
        inner_mask_binary = numpy.zeros(self._grid_idxr.size_1d(), dtype=bool)
        inner_mask_binary.put(self._inner_mask_indices, True)
        self._inner_mask_binary = inner_mask_binary

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

    def summary(self):
        return '\n'.join(['----------------------------------->>>',
                          'Local Mask Summary:',
                          'Number of Mask Points:  {!s}'.format(len(self.mask())),
                          'Mask Radius (Cart):     {!s}'.format(self.radius()),
                          'Mask Volume (Cart):     {!s}'.format(round(self.volume(),3)),
                          'Largest Mask Vector:    {!s}'.format(max(self.mask())),
                          'Req. Edge Buffer Zone:  {!s}'.format(self.buffer_size()),
                          'Sampling Fraction (1D): 1/{!s}'.format(self.grid_jump()),
                          'Sampling Fraction (3D): 1/{!s}'.format(self.grid_jump()**3)
                        ])



import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy as np


class GridMask(object):

    name = "GridMask"

    def __init__(self,
        grid_size,
        mask_indices = None,
        mask_binary = None,
        ):

        self.grid_size = grid_size

        if [(mask_indices is None), (mask_binary is None)].count(True) != 1:
            raise ValueError("must supply mask_indices OR mask_binary")

        self.mask_indices = (
            np.array(mask_indices, dtype=int).flatten() 
            if (mask_indices is not None) 
            else None
            )
        self.mask_binary = (
            np.array(mask_binary, dtype=bool).flatten() 
            if (mask_binary is not None) 
            else None
            )

        self.validate()

        # Populate derived values (after validation!)

        self.mask_size = len(self.get_mask_indices())

    def __str__(self):
        
        s_ = (
            'Object: {name}\n'
            'Grid Size: {grid_size}\n'
            'Mask Size: {mask_size}\n'
            'Mask Percentage: {mask_frac:.1%}'
            ).format(
            name = self.name,
            grid_size = self.grid_size,
            mask_size = self.mask_size,
            mask_frac = (
                float(self.mask_size) / float(np.product(self.grid_size))
                ),
            )

        return s_

    def validate(self):

        if (self.mask_binary is not None): 
            if len(self.mask_binary) != np.product(self.grid_size):
                raise ValueError("mask_binary is not the same size as grid_size")

        if (self.mask_indices is not None): 
            if max(self.mask_indices) >= np.product(self.grid_size):
                raise ValueError("Element in mask_indices is greater than grid size")
            if min(self.mask_indices) < 0: 
                raise ValueError("Element in mask_indices is negative")

    def get_mask_binary(self):

        if (self.mask_binary is not None):
            return self.mask_binary

        mask_binary = np.zeros(np.product(self.grid_size), dtype=bool)
        mask_binary[self.mask_indices] = 1

        return mask_binary

    def get_mask_indices(self):

        if (self.mask_indices is not None):
            return self.mask_indices

        mask_indices = np.where(self.mask_binary)[0]

        return mask_indices

    def get_mask_points(self):

        all_grid_points = np.transpose(
            np.indices(
                dimensions = self.grid_size,
                dtype = int,
                ), 
            axes = (1,2,3,0),
            ).reshape(-1, 3)

        return all_grid_points[self.get_mask_binary()]

    def index_on_other(self, other):
        """Create a mask for this mask relative to another mask"""

        import copy

        self_binary = self.get_mask_binary()
        other_binary = other.get_mask_binary()

        new_binary = copy.deepcopy(self_binary)
        new_binary[np.logical_not(other_binary)] = False

        new_mask = GridMask(
            grid_size = self.grid_size,
            mask_binary = new_binary,
            )

        return new_mask, self_binary[other_binary]

    def embed_data(self, data, masked_value=0.0):

        embedded = (
            masked_value * np.ones(np.product(self.grid_size))
            )

        embedded[self.get_mask_binary()] = np.array(data).flatten()
        
        embedded = embedded.reshape(self.grid_size)

        return embedded


class GetSitesMask(object):

    def __init__(self,
        grid_size,
        cart_origin,
        unit_cell,
        mask_dist,
        ):

        self.grid_size = grid_size
        self.cart_origin = cart_origin
        self.unit_cell = unit_cell

        self.mask_dist = mask_dist

    @classmethod
    def from_map_grid(cls, 
        map_grid, 
        mask_dist,
        ):

        return cls(
            grid_size = map_grid.grid_size(),
            cart_origin = map_grid.cart_origin(),
            unit_cell = map_grid.unit_cell(),
            mask_dist = mask_dist,
            )

    def __call__(self,        
        sites_cart,
        ):

        from scitbx.array_family import flex
        from cctbx import maptbx

        sites_cart = flex.vec3_double(sites_cart)

        # Calculate the masked indices defined by max distance from protein atoms
        mask_indices = maptbx.grid_indices_around_sites(
            unit_cell = self.unit_cell,
            fft_n_real = self.grid_size,
            fft_m_real = self.grid_size,
            # Masking is performed relative to the grid origin, so need to apply origin shift
            sites_cart = (sites_cart - self.cart_origin),
            site_radii = flex.double(sites_cart.size(), self.mask_dist),
            )

        return GridMask(
            grid_size = self.grid_size,
            mask_indices = np.array(mask_indices, dtype=int),
            )


class GetNonPeriodicSitesMask(object):

    def __init__(self, 
        grid_size,
        cart_origin,
        unit_cell,
        mask_dist,
        ):

        self.grid_size = grid_size
        self.cart_origin = cart_origin
        self.unit_cell = unit_cell

        self.mask_dist = mask_dist

    @classmethod
    def from_map_grid(cls, 
        map_grid, 
        mask_dist,
        ):

        return cls(
            grid_size = map_grid.grid_size(),
            cart_origin = map_grid.cart_origin(),
            unit_cell = map_grid.unit_cell(),
            mask_dist = mask_dist,
            )

    def __call__(self,
        sites_cart,
        ):

        from scitbx.array_family import flex

        sites_cart = flex.vec3_double(sites_cart)

        filtered_sites_cart = self.filter_sites(
            sites_cart = sites_cart,
            )

        mask = self.mask_sites(
            sites_cart = filtered_sites_cart,
            )

        return mask

    def filter_sites(self, sites_cart):

        from scitbx.array_family import flex

        buffer_cart = (self.mask_dist, ) * 3
        buffer_grid = self.unit_cell.fractionalize(
            flex.vec3_double([buffer_cart])
            )[0]

        # Map to grid
        sites_grid = (sites_cart - self.cart_origin)
        sites_frac = self.unit_cell.fractionalize(sites_grid)

        selection_bool = flex.bool(len(sites_cart), True)

        for filter_selection in [
            (sites_frac.dot((1.0,0.0,0.0)) > 1.0+buffer_grid[0]),
            (sites_frac.dot((0.0,1.0,0.0)) > 1.0+buffer_grid[1]),
            (sites_frac.dot((0.0,0.0,1.0)) > 1.0+buffer_grid[2]),
            (sites_frac.dot((1.0,0.0,0.0)) < 0.0-buffer_grid[0]),
            (sites_frac.dot((0.0,1.0,0.0)) < 0.0-buffer_grid[1]),
            (sites_frac.dot((0.0,0.0,1.0)) < 0.0-buffer_grid[2]),
            ]:

            selection_bool.set_selected(filter_selection, False)

        return sites_cart.select(selection_bool)

    def mask_sites(self, sites_cart):

        from scitbx.array_family import flex
        from cctbx import maptbx

        # Find out how many unit cells we need to avoid loop-round
        buffer_cart = (2.0 * self.mask_dist, ) * 3
        buffer_grid = self.unit_cell.fractionalize(
            flex.vec3_double([buffer_cart])
            )[0]
        # Number of unit cells (1 + buffer) along each axis to avoid wrap around
        n_unit_cells = tuple(
            1 + int(np.ceil(buffer_g))
            for buffer_g in buffer_grid
            )
        # Create an expanded unit cell
        from cctbx.uctbx import unit_cell
        expanded_unit_cell = unit_cell(
            tuple(
                n * p 
                for n, p in zip(
                    tuple(n_unit_cells)+(1.,1.,1.), 
                    self.unit_cell.parameters()
                    )
                )
            )
        # Number of grid points in the expanded map
        expanded_grid_size = tuple(
            n * s
            for n, s in zip(
                n_unit_cells,
                self.grid_size,
                )
            )

        # Calculate the masked indices defined by max distance from protein atoms
        mask_indices = maptbx.grid_indices_around_sites(
            unit_cell = expanded_unit_cell,
            fft_n_real = expanded_grid_size,
            fft_m_real = expanded_grid_size,
            # Masking is performed relative to the grid origin, so need to apply origin shift
            sites_cart = (sites_cart - self.cart_origin),
            site_radii = flex.double(sites_cart.size(), self.mask_dist),
            )

        # Create array for expanded array
        expanded_grid_mask = np.zeros(np.product(expanded_grid_size), dtype=bool)
        if len(mask_indices) > 0:
            expanded_grid_mask[np.array(mask_indices)] = True
        expanded_grid_mask = expanded_grid_mask.reshape(expanded_grid_size)

        # Now trim to the original grid
        grid_mask = expanded_grid_mask[
            0:self.grid_size[0],
            0:self.grid_size[1],
            0:self.grid_size[2]
            ]        

        return GridMask(
            grid_size = self.grid_size,
            mask_binary = grid_mask,
            )


def compound_grid_masks(
    positive_masks,
    negative_masks = None,
    ):

    if (positive_masks is None) or (len(positive_masks) == 0):
        raise ValueError("must provide at least one positive_mask")
    
    positive_masks = (positive_masks if (positive_masks is not None) else [])
    negative_masks = (negative_masks if (negative_masks is not None) else [])

    grid_size = positive_masks[0].grid_size

    for m in (positive_masks+negative_masks):
        assert m.grid_size == grid_size

    total_binary = np.zeros(np.product(grid_size), dtype=bool)

    for m in positive_masks: 
        total_binary[m.get_mask_binary()] = 1

    for m in negative_masks: 
        total_binary[m.get_mask_binary()] = 0

    return GridMask(
        grid_size = grid_size,
        mask_binary = total_binary,
        )

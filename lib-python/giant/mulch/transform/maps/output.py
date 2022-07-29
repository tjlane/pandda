import numpy as np

import copy, time

import cctbx.sgtbx, cctbx.uctbx, cctbx.maptbx, scitbx.matrix, iotbx.ccp4_map 

from scitbx.array_family import flex

def write_map(map_data, dataset, filename):

    iotbx.ccp4_map.write_ccp4_map(
        file_name = filename,
        unit_cell = dataset.model.crystal.unit_cell,
        space_group = dataset.model.crystal.space_group,
        map_data = flex.double(map_data),
        labels = flex.std_string(['Map from pandda'])
        )


class MapGridMapWriter(object):

    def __init__(self,
        map_grid,
        ):

        self.grid_origin = map_grid.cart_origin()
        self.grid_unit_cell = map_grid.unit_cell()

    def __call__(self,
        map_data,
        filename,
        ):

        map_data = flex.double(map_data)
        
        assert map_data.nd() == 3

        t_map_data = cctbx.maptbx.rotate_translate_map(
            unit_cell = self.grid_unit_cell,
            map_data = map_data,
            rotation_matrix = scitbx.matrix.rec([1,0,0,0,1,0,0,0,1], (3,3)).elems,
            translation_vector = (-1.0*scitbx.matrix.rec(list(self.grid_origin), (3,1))).elems,
            )
                
        iotbx.ccp4_map.write_ccp4_map(  
            file_name   = filename,
            unit_cell   = self.grid_unit_cell,
            space_group = cctbx.sgtbx.space_group('P1'),
            map_data    = t_map_data,
            labels      = flex.std_string(['Map from pandda'])     
            )


class NativeMapMaker(object):

    debug = False

    def __init__(self,
        map_grid,
        sites_mask, 
        sites_mask_radius,
        native_map_grid_spacing = None,
        ):

        if (native_map_grid_spacing is None):
            native_map_grid_spacing = map_grid.grid_spacing()

        # Used for testing
        self._write_grid_map = MapGridMapWriter(
            map_grid = map_grid,
            )

        # Map properties
        self.grid_origin = map_grid.cart_origin()
        self.grid_spacing = map_grid.grid_spacing()
        self.grid_unit_cell = map_grid.unit_cell()

        # Sites to define the map regions in the output map
        # Given in the map reference frame
        self.sites_mask = sites_mask
        self.sites_mask_radius = sites_mask_radius

        self.native_map_grid_spacing = native_map_grid_spacing

    def __call__(self, map_data, dataset, filename):
        """Process the dataset"""
        
        if self.debug:
            self._write_grid_map(
                map_data = map_data, 
                filename = filename + '-gridtest.ccp4',
                )

        native_map_data = self.embed_data(
            map_data = map_data,
            dataset = dataset,
            )

        if (filename is not None):
            write_map(
                map_data = native_map_data,
                dataset = dataset,
                filename = filename,
                )

        return native_map_data

    def embed_data(self, 
        map_data,
        dataset,
        ):

        reference_map = self.get_reference_map(
            map_data = map_data,
            )

        native_map_data = self.get_native_map_data(
            reference_map = reference_map,
            dataset = dataset,
            )

        return native_map_data

    def get_reference_map(self, map_data):

        from giant.mulch.transform.maps.utils import (
            SampleableMap,
            )

        return SampleableMap(
            map_data = map_data,
            unit_cell = self.grid_unit_cell,
            origin = self.grid_origin,
            filter_out_of_bounds_points = True, # !!! makes map non-periodic! 
            )

    def get_native_map_data(self, 
        reference_map, 
        dataset,
        ):

        native_map_data = self.embed_grid_data(
            reference_map = reference_map,
            dataset = dataset,
            )

        native_map_data = self.symmetrise_map(
            map_data = native_map_data,
            dataset = dataset,
            )

        return native_map_data

    def embed_grid_data(self, 
        reference_map,
        dataset,
        ):

        # Map sites to dataset space
        native_sites = dataset.model.alignment.ref2nat(
            self.sites_mask
            )

        # Calculate supercell required to contain sites
        sc_size = self.calculate_supercell_size(
            sites_cart = native_sites,
            dataset = dataset, 
            )

        # Create supercell
        sc_unit_cell = self.get_supercell_unit_cell(
            unit_cell = dataset.model.crystal.unit_cell,
            supercell_size = sc_size,
            )

        # Get gridding for output map
        uc_gridding = self.get_native_gridding(
            dataset = dataset,
            )

        # Expand to supercell
        sc_gridding = self.get_supercell_gridding(
            unit_cell_gridding = uc_gridding,
            supercell_unit_cell = sc_unit_cell,
            supercell_size = sc_size,
            )

        # Get vector to centre sites on supercell
        shift_vector = self.get_vector_to_centre_sites_on_unit_cell(
            sites_cart = native_sites,
            unit_cell = sc_unit_cell,
            )

        # Get points in the supercell to be sampled from reference map
        sc_mask = self.mask_unit_cell(
            sites_cart = (native_sites + shift_vector),
            gridding = sc_gridding,
            unit_cell = sc_unit_cell,
            )

        # Get the cartesian positions of the points in the supercell
        sc_sample_points = self.map_grid_points_to_cart_points(
            sites_grid = sc_mask.get_mask_points(),
            gridding = sc_gridding,
            unit_cell = sc_unit_cell,
            )

        # Shift back to the unit cell, map to grid, and sample
        sc_sample_values = reference_map.get_cart_values(
            dataset.model.alignment.nat2ref(
                sc_sample_points - shift_vector
                )
            )

        # Embed the sampled data in the supercell
        supercell_map_data = sc_mask.embed_data(
            data = sc_sample_values,
            )

        # Reverse the origin shift
        supercell_map_data = self.shift_map_data(
            map_data = supercell_map_data,
            vector = shift_vector, # maybe needs to be reversed! 
            unit_cell = sc_unit_cell,
            )

        supercell_map_data = self.roll_expand_supercell_data(
            supercell_map_data = supercell_map_data,
            supercell_size = sc_size,
            )

        uc_map_data = self.extract_subcell_of_supercell(
            supercell_map_data = supercell_map_data,
            supercell_size = sc_size,
            )

        return uc_map_data

    def calculate_supercell_size(self, sites_cart, dataset):
        """How big does the supercell need to be to enclose all of the mask points"""

        box_frac_min_max = dataset.model.crystal.unit_cell.box_frac_around_sites(
            sites_cart = sites_cart, 
            buffer = self.sites_mask_radius,
            )

        supercell_size = tuple(
            int(np.ceil(ma-mi)) for mi, ma in 
            zip(*box_frac_min_max)
            )

        return supercell_size

    def get_supercell_unit_cell(self,
        unit_cell,
        supercell_size,
        ):

        assert len(supercell_size) == 3

        uc_params = list(unit_cell.parameters())

        for i, n in enumerate(supercell_size):

            uc_params[i] = (uc_params[i] * float(n))

        return cctbx.uctbx.unit_cell(uc_params)

    def get_native_gridding(self,
        dataset,
        ):

        map_gridding = cctbx.maptbx.crystal_gridding(
            unit_cell = dataset.model.crystal.unit_cell, 
            step = self.native_map_grid_spacing,
            )

        return map_gridding

    def get_supercell_gridding(self,
        unit_cell_gridding, 
        supercell_unit_cell,
        supercell_size = (2, 2, 2),
        ):

        assert len(supercell_size) == 3

        supercell_nreal = tuple(
            int(s*n) for s,n in 
            zip(
                supercell_size, 
                unit_cell_gridding.n_real()
                )
            )

        sc_map_gridding = cctbx.maptbx.crystal_gridding(
            unit_cell = supercell_unit_cell, 
            pre_determined_n_real = supercell_nreal,
            )

        return sc_map_gridding

    def get_vector_to_centre_sites_on_unit_cell(self,
        sites_cart,
        unit_cell,
        ):
        """sites_cart + shift_vector = unit_cell_centre"""

        sites_cart = np.array(sites_cart)
        
        sites_centre = tuple(
            (
                sites_cart.max(axis=0) + 
                sites_cart.min(axis=0)
                ) / 2.0
            )

        uc_centre = unit_cell.orthogonalize((0.5, 0.5, 0.5))

        shift_vector = tuple(
            a-b for a,b in zip(uc_centre, sites_centre)
            )

        return shift_vector

    def mask_unit_cell(self,
        sites_cart,
        gridding,
        unit_cell,
        ):

        from .grid.mask import GridMask

        sites_cart = flex.vec3_double(sites_cart)

        masked_points_indices = cctbx.maptbx.grid_indices_around_sites(
            unit_cell = unit_cell,
            fft_n_real = gridding.n_real(),
            fft_m_real = gridding.n_real(),
            sites_cart = sites_cart,
            site_radii = flex.double(
                sites_cart.size(), 
                self.sites_mask_radius, # !!!
                ), 
            )

        map_mask = GridMask(
            grid_size = gridding.n_real(),
            mask_indices = np.array(masked_points_indices),
            )
        
        return map_mask

    def map_grid_points_to_cart_points(self, 
        sites_grid,
        gridding,
        unit_cell,
        ):

        g2c = cctbx.maptbx.grid2cart(
            gridding.n_real(), 
            unit_cell.orthogonalization_matrix(),
            )

        sites_cart = np.array(
            list(map(g2c, map(tuple, sites_grid)))
            )

        return sites_cart

    def shift_map_data(self,
        map_data,
        vector, 
        unit_cell,
        ):

        #vector = tuple(-v for v in vector)

        shift_map_data = cctbx.maptbx.rotate_translate_map(
            unit_cell = unit_cell,
            map_data = flex.double(map_data),
            rotation_matrix = scitbx.matrix.rec([1,0,0,0,1,0,0,0,1], (3,3)).elems,
            translation_vector = (scitbx.matrix.rec(vector, (3,1))).elems,
            )

        shift_map_data = shift_map_data.as_numpy_array()

        assert shift_map_data.ndim == 3

        return shift_map_data

    def roll_expand_supercell_data(self,
        supercell_map_data,
        supercell_size,
        ):

        n_total = supercell_map_data.shape
        n_shift = tuple(
            int(n / s) for n,s in zip(n_total, supercell_size)
            )

        total_map_data = copy.deepcopy(supercell_map_data)

        for xyz in flex.nested_loop(supercell_size):

            if xyz == (0,0,0):
                continue

            n_roll = tuple(
                n*int(i) 
                for n,i in zip(n_shift, xyz)
                )

            rolled_map_data = np.roll(
                supercell_map_data,
                shift = n_roll,
                axis = (0, 1, 2),
                )

            # Mask values that already have values in the output array
            total_mask = (total_map_data != 0.0)
            rolled_map_data[total_mask] = 0.0

            # Add to output array
            total_map_data += rolled_map_data

        return total_map_data

    def extract_subcell_of_supercell(self,
        supercell_map_data,
        supercell_size,
        ):
        
        n_all = supercell_map_data.shape
        n_one = tuple(
            int(n / s) for n,s in zip(n_all, supercell_size)
            )

        n1, n2, n3 = n_one

        one_map_data = supercell_map_data[0:n1, 0:n2, 0:n3]

        assert one_map_data.shape == n_one

        return np.array(one_map_data)

    def symmetrise_map(self,
        map_data,
        dataset,
        ):

        space_group = dataset.model.crystal.space_group
        unit_cell = dataset.model.crystal.unit_cell

        # Create copy to store output
        total_map_data = copy.deepcopy(map_data)

        map_data = flex.double(map_data)
        assert map_data.nd() == 3

        # Apply all symmetry operations to unit cell data
        for sym_op in space_group.all_ops():

            if sym_op.as_xyz() == 'x,y,z': 
                continue

            # Get the transformation matrix
            rt_mx = sym_op.as_rational().as_float()

            # Tranform the map
            rt_map_data = cctbx.maptbx.rotate_translate_map(
                    unit_cell = unit_cell,
                    map_data = map_data,
                    rotation_matrix = rt_mx.r.elems,
                    #translation_vector = native_unit_cell.orthogonalize((-1.0*rt_mx.t).elems)   )
                    translation_vector = unit_cell.orthogonalize(rt_mx.t.elems),
                    )
            rt_map_data = rt_map_data.as_numpy_array()
            assert rt_map_data.ndim == 3

            total_mask = (total_map_data != 0.0)
            rt_map_data[total_mask] = 0.0

            total_map_data += (rt_map_data)

        return total_map_data


class MaskedNativeMapMaker(object):

    def __init__(self,
        map_grid,
        map_mask,
        ):
        """Take the sparse map data and put into the reference grid"""

        self.map_mask = map_mask

        self.make_map = NativeMapMaker(
            map_grid = map_grid,
            sites_mask = map_grid.grid2cart(
                self.map_mask.get_mask_points(),
                origin_shift = True,
                ),
            # only need to include neighbouring grid points of the mask
            sites_mask_radius = 2. * map_grid.grid_spacing(),
            )
        
    def __call__(self, map_data, dataset, filename):

        map_data_dense = self.map_mask.embed_data(
            data = map_data,
            )

        return self.make_map(
            map_data = map_data_dense,
            dataset = dataset,
            filename = filename,
            )


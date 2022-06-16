import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy as np

from giant.mulch.reference import DefaultReferenceDataset


def calculate_grid_size(min_carts, max_carts, grid_spacing):
    """Calculate the number of points to be sampled for a box size and sampling distance. Returns the number of points to be sampled along each axis. Box may be larger than max_carts."""

    cart_size = tuple([
        float(max_c - min_c) 
        for min_c, max_c in zip(min_carts, max_carts)
        ])

    for c in cart_size: 
        if c <= 0.0: 
            raise ValueError('at least one dimension has size <= 0.0')
    
    grid_spacing = float(grid_spacing)

    grid_size = tuple([
        1 + int(np.ceil(c_size/grid_spacing))
        for c_size in cart_size
        ])

    return grid_size


class GridIndexer(object):

    def __init__(self, grid_size):

        self.grid_size = grid_size
        self.norm_0 = (grid_size[0] * grid_size[1] * grid_size[2])
        self.norm_1 = (grid_size[1] * grid_size[2])
        self.norm_2 = (grid_size[2])

    def __call__(self, grid_point):

        return (
            (grid_point[0] * self.norm_1) + 
            (grid_point[1] * self.norm_2) + 
            (grid_point[2])
            )


class GridInverter(GridIndexer):

    def __call__(self, grid_index):

        i1 = (grid_index) // self.norm_1
        i2 = (grid_index % self.norm_1) // self.norm_2 
        i3 = (grid_index % self.norm_1) % self.norm_2

        return (i1, i2, i3)


class MultiPointGridIndexer(GridIndexer):

    def __call__(self, grid_points):

        grid_points = np.array(grid_points)

        return (
            (grid_points[:,0] * self.norm_1) + 
            (grid_points[:,1] * self.norm_2) + 
            (grid_points[:,2])
            )


class MultiPointGridInverter(GridInverter):

    def __call__(self, grid_indices):

        grid_indices = np.array(grid_indices)

        i1 = (grid_indices) // self.norm_1
        i2 = (grid_indices % self.norm_1) // self.norm_2 
        i3 = (grid_indices % self.norm_1) % self.norm_2

        return np.transpose(
            np.array([i1,i2,i3])
            )


class MapGrid(object):

    def __init__(self, 
        grid_spacing, 
        cart_origin, 
        cart_approx_max,
        ):

        self._grid_spacing = float(grid_spacing)
        self._grid_size = calculate_grid_size(
            min_carts = cart_origin, 
            max_carts = cart_approx_max, 
            grid_spacing = float(grid_spacing),
            )
        self._cart_origin = tuple(cart_origin)

        self.masks = {}
        self.partition = None

    def __str__(self):

        s_ = (
            'Map Grid Summary:\n'
            '| Grid Origin (cartesian): {cart_origin}\n'
            '| Grid Size: {grid_size}\n'
            '| Grid Points: {grid_size_1d}\n'
            '| Grid Spacing: {grid_spacing}\n'
            '| Unit Cell: {unit_cell}\n'
            '`---->\n'
            ).format(
                cart_origin = np.round(self.cart_origin(),3),
                grid_size = self.grid_size(),
                grid_size_1d = self.grid_size_1d(),
                grid_spacing = self.grid_spacing(),
                unit_cell = str(self.unit_cell()),
            )

        for m_key, m_mask in self.masks.items():

            s_ += (
                '| Mask: {name}\n'
                '|\t{mask}\n'
                '`---->\n'
                ).format(
                name = m_key,
                mask = str(m_mask).strip().replace('\n', '\n|\t'),
                )

        if (self.partition is not None): 

            s_ += (
                '| Partition:\n'
                '|\t{partition}\n'
                '`---->'
                ).format(
                partition = str(self.partition).strip().replace('\n', '\n|\t'),
                )

        return s_.strip()

    def copy(self, lightweight=True):

        new = copy.deep_copy(self)

        if (lightweight is True):
            new.masks = None
            new.partition = None

        return new

    def grid_spacing(self): 
        return self._grid_spacing

    def grid_size(self):
        return self._grid_size

    def grid_size_1d(self):
        return np.product(self.grid_size())

    def cart_size(self):
        return tuple(s * self.grid_spacing() for s in self.grid_size())

    def cart_origin(self):
        return self._cart_origin

    def grid_points(self):

        all_grid_points = np.transpose(
            np.indices(
                dimensions = self.grid_size(),
                dtype = int,
                ), 
            axes = (1,2,3,0),
            ).reshape(-1, 3)

        return all_grid_points

    def cart_points(self, origin_shift=True):
        sites = self.grid_points().astype(np.float)
        sites *= self.grid_spacing()
        if (origin_shift is True):
            sites += self.cart_origin()
        return sites 

    def grid2cart(self, sites_grid, origin_shift=True):
        assert isinstance(origin_shift, bool)
        sites = np.array(
            sites_grid, 
            copy = True, 
            dtype = np.float,
            )
        sites *= self.grid_spacing()
        if (origin_shift is True): 
            sites += self.cart_origin()
        return sites

    def cart2grid(self, sites_cart, origin_shift=True):
        assert isinstance(origin_shift, bool)
        sites = np.array(
            sites_cart, 
            copy = True,
            dtype = np.float,
            )
        if (origin_shift is True): 
            sites -= self.cart_origin()
        sites /= self.grid_spacing()
        return sites

    def embed_data(self, data, masked_value=0.0):

        embedded = np.array(
            data, 
            copy = True,
            )
        embedded = embedded.reshape(self.grid_size())
        return embedded

    def unit_cell(self):
        """Create a unit cell as if the reference grid were a real lattice"""
        import cctbx.uctbx
        parameters = tuple(self.cart_size())+(90.,90.,90.)
        return cctbx.uctbx.unit_cell(parameters)

    # def space_group(self):
    #     """Create a spacegroup as if the reference grid were a real lattice (Cartesian Grid == P1 SpaceGroup)"""
    #     return cctbx.sgtbx.space_group('P1')

    # def crystal_symmetry(self):
    #     """Create a symmetry object as if the reference grid were a real lattice (THE UNIT CELL IS LARGER THAN cart_size() BY 1 GRID SPACING)"""
    #     return iotbx.crystal_symmetry_from_any.from_string("{:f},{:f},{:f},90,90,90,P1".format(*list(self.cart_size())))


class GetWarpedMapGrid(object):

    ReferenceDatasetClass = DefaultReferenceDataset

    def __init__(self,
        map_grid_spacing,
        outer_mask_radius,    # define the map zone (map extent)
        inner_mask_radius,    # define the internal zone (near to atoms)
        symmetry_mask_radius, # remove symmetry atoms
        mask_pdb = None, 
        align_mask_to_reference = False,
        create_grid_selection_string = "not hetero", # mask used to create the grid
        mask_grid_selection_string = "not hetero", # mask used to define regions on the grid
        partition_grid_selection_string = "pepnames and name CA and (altloc ' ' or altloc 'A')",
        processor = None,
        ):

        if (processor is None):
            from giant.processors import basic_processor
            processor = basic_processor

        self.map_grid = None
        self.map_grid_spacing = map_grid_spacing

        self.outer_mask_radius = outer_mask_radius
        self.inner_mask_radius = inner_mask_radius
        self.symmetry_mask_radius = symmetry_mask_radius

        self.mask_pdb = mask_pdb
        self.mask_dataset = None
        self.align_mask_to_reference = align_mask_to_reference
        
        self.create_grid_selection_string = create_grid_selection_string
        self.mask_grid_selection_string = mask_grid_selection_string
        self.partition_grid_selection_string = partition_grid_selection_string

        self.processor = processor

    def __call__(self, 
        reference_dataset = None,
        ):
        """Generate the grid objects for the analysis"""

        if (self.map_grid is None): 
            self.map_grid = self.make_grid(
                reference_dataset = reference_dataset,
                )
        
        return self.map_grid

    def make_grid(self, 
        reference_dataset = None,
        ):
        
        logger.heading('Creating aligned map grid')

        #logger(str(self))

        logger.subheading('Input reference dataset')
        logger(str(reference_dataset))

        # Which dataset to be used to mask the grid
        logger.subheading('Initialising masking dataset')
        mask_dataset = self.initialise_mask_dataset(
            reference_dataset = reference_dataset,
            )
        logger(str(mask_dataset))

        # Create the grid using the masking dataset (for determining size and extent of grid)
        logger.subheading("Determining map grid extent")
        map_grid = self.create_grid(
            dataset = mask_dataset, 
            selection_string = self.create_grid_selection_string,
            )
        logger(str(map_grid))

        logger.subheading("Masking map grid")
        map_grid.masks = self.mask_grid(
            map_grid = map_grid,
            dataset = reference_dataset,
            selection_string = self.mask_grid_selection_string,
            )
        logger(str(map_grid))

        # Partition the grid with the reference dataset (which grid points use which coordinate transformations)
        logger.subheading("Partitioning map grid")
        map_grid.partition = self.partition_grid(
            map_grid = map_grid,
            dataset = reference_dataset,
            selection_string = self.partition_grid_selection_string,
            )
        logger('\n'+str(map_grid.partition))

        self.mask_dataset = mask_dataset

        return map_grid

    def initialise_mask_dataset(self, 
        reference_dataset = None,
        ):

        # Initial checks
        if (reference_dataset is None) and (self.mask_pdb is None):
            raise ValueError('Must provide reference_dataset if mask_pdb is None')

        mask_dataset = None

        #####
        # Initialise mask_dataset from pdb
        if (self.mask_pdb is not None): 

            mask_dataset = self.ReferenceDatasetClass.from_file(
                model_filename = self.mask_pdb,
            ).label(
                tag = "masking_dataset",
            )

        #####
        # Select the mask dataset and align to reference if needed
        #
        if (mask_dataset is None):

            # Mask dataset is copy of reference dataset
            mask_dataset = reference_dataset.copy()

        elif (self.align_mask_to_reference is True):

            # Align input mask to reference
            if (reference_dataset is None):
                raise ValueError("Must provide reference_dataset if align_mask_to_reference is True")

            try:
                logger("Aligning mask_dataset to reference_dataset")

                alignment = mask_dataset.model.align_to(
                    other_hierarchy = reference_dataset.model.hierarchy,
                    method = "global",
                    require_hierarchies_identical = False,
                    )

            except Exception as e:
                msg = traceback.format_exc()
                msg += '\n------------------>>>'
                msg += '\n\nFailed to align masking pdb ({}) to the reference structure.'.format(
                    self.mask_pdb,
                    )
                msg += '\nIf the masking structure does not need alignment, rerun with align_mask_to_reference=False'
                raise Failure(msg)

            # Apply the alignment to the structure
            mask_dataset.model.hierarchy.atoms().set_xyz(
                alignment.nat2ref(
                    mask_dataset.model.hierarchy.atoms().extract_xyz(),
                    )
                )

        assert mask_dataset.model.hierarchy.atoms_size() > 0

        return mask_dataset

    def create_grid(self, 
        dataset,
        selection_string = None,
        ):
        """Create a grid over the given dataset"""

        mask_h = dataset.model.hierarchy

        if (selection_string is not None):

            asc = mask_h.atom_selection_cache()
            mask_h = mask_h.select(
                atom_selection = asc.selection(selection_string),
                )

            if mask_h.atoms_size() == 0:
                raise ValueError(
                    (
                        "The input mask_selection_string '{}' "
                        "does not select any atoms from the mask dataset"
                        ).format(selection_string)
                    )

        sites_cart = mask_h.atoms().extract_xyz()        

        # Calculate the extent of the grid
        grid_min = tuple(s - self.outer_mask_radius for s in sites_cart.min())
        grid_max = tuple(s + self.outer_mask_radius for s in sites_cart.max())

        # ============================================================================>
        # Create main grid object
        # ============================================================================>
        map_grid = MapGrid(
            grid_spacing = self.map_grid_spacing,
            cart_origin = grid_min,
            cart_approx_max = grid_max,
            )

        return map_grid

    def mask_grid(self, 
        map_grid, 
        dataset, 
        selection_string = None,
        ):
        """Create masks for the reference grid based on distances from atoms in the reference structure"""

        from giant.mulch.transform.maps.grid.mask import \
            GetNonPeriodicSitesMask, compound_grid_masks

        #####
        # Apply selection string to the masked hierarchy
        #

        main_h = dataset.model.hierarchy

        if (selection_string is not None):
            asc = main_h.atom_selection_cache()
            main_h = main_h.select(
                asc.selection(selection_string), 
                copy_atoms = True,
                )

        if len(main_h.atoms()) == 0: 
            raise ValueError('Zero atoms have been selected to mask the grid')

        #####
        # Apply selection string to the crystal contacts of the full input hierarchy
        #

        sym_h = dataset.model.crystal_contacts(
            distance_cutoff = (
                self.outer_mask_radius + self.symmetry_mask_radius
                ), # max radius that can affect map 
            combine_copies = True,
            )

        if (selection_string is not None):
            asc = sym_h.atom_selection_cache()
            sym_h = sym_h.select(
                asc.selection(selection_string), 
                copy_atoms = True,
                )

        # if len(sym_h.atoms()) == 0: 
        #     raise ValueError('Zero atoms have been selected to mask the grid')

        #####
        # Get coordinates for masks
        #

        ref_sites_cart = (
            main_h.atoms().extract_xyz()
            )

        sym_sites_cart = (
            sym_h.atoms().extract_xyz()
            )

        #####
        # Mask getters

        make_outer_mask = GetNonPeriodicSitesMask.from_map_grid(
            map_grid = map_grid,
            mask_dist = self.outer_mask_radius,
            )

        make_inner_mask = GetNonPeriodicSitesMask.from_map_grid(
            map_grid = map_grid,
            mask_dist = self.inner_mask_radius,
            )

        make_symmetry_mask = GetNonPeriodicSitesMask.from_map_grid(
            map_grid = map_grid,
            mask_dist = self.symmetry_mask_radius,
            )

        #####
        # Make the masks

        masks = {}

        masks['outer'] = make_outer_mask(
            sites_cart = ref_sites_cart,
            )

        masks['inner'] = make_inner_mask(
            sites_cart = ref_sites_cart,
            )

        masks['symmetry'] = make_symmetry_mask(
            sites_cart = sym_sites_cart,
            )

        masks['total'] = compound_grid_masks(
            positive_masks = [masks['outer']],
            negative_masks = [masks['inner'], masks['symmetry']],
            )

        return masks

    def partition_grid(self, 
        map_grid, 
        dataset, 
        selection_string,
        ):

        import time

        from .partition import (
            MakeVoronoiGridPartition,
            )

        make_partition = MakeVoronoiGridPartition(
            processor = self.processor,
            )

        partition = make_partition(
            map_grid = map_grid,
            dataset = dataset,
            selection_string = selection_string,
            mask_name = "outer",
        )

        return partition

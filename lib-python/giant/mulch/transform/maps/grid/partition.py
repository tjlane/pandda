import numpy as np


class GridPartition(object):

    name = "GridPartition"

    def __init__(self,
        map_grid,
        partition_sites,
        partition_mappings,
        grid_mask = None,
        ):

        self.grid_size = map_grid.grid_size()
        self.grid_mask = grid_mask # Can be None
        
        self.partition_sites = partition_sites
        self.partition_mappings = np.array(partition_mappings, dtype=int)
        self.assign_sites = AssignSites(ref_sites=self.partition_sites)

        self.validate()

        from giant.mulch.transform.maps.grid import MultiPointGridIndexer
        self.gps2idx = MultiPointGridIndexer(
            grid_size = self.grid_size,
            )

    def __str__(self):

        assigned_points = dict(
            zip(*np.unique(self.partition_mappings, return_counts=True))
            )

        point_strings = [
            'Site {idx}, {site}: {count} point(s)'.format(
                idx = i+1, 
                site = str(tuple(xyz)),
                count = assigned_points.get(i, 'no'),
                )
            for i, xyz in enumerate(self.partition_sites)
        ]

        s_ = (
            'Object: {name}\n'
            '| Grid Size: {grid_size}\n'
            '| Partition Sites: {n_partition_sites}\n'
            '| Masked Points (not assigned): {masked_points}\n'
            '| Grid site assignments:\n'
            '|\t{grid_assignments}\n'
            '`---->\n'
            ).format(
            name = self.name,
            grid_size = self.grid_size,
            n_partition_sites = len(self.partition_sites),
            masked_points = (
                (np.product(self.grid_size) - self.grid_mask.mask_size)
                if (self.grid_mask is not None)
                else "None"
                ),
            grid_assignments = (
                '\n'.join(point_strings).strip().replace('\n','\n|\t')
                if point_strings else "None"
                )
            )

        return s_.strip()

    def validate(self):

        if (self.grid_mask is not None):
            assert self.grid_size == self.grid_mask.grid_size
            assert len(self.partition_mappings) == self.grid_mask.mask_size
        else: 
            assert len(self.partition_mappings) == np.product(self.grid_size)

    def embed(self):

        if (self.grid_mask is not None):
            mappings = self.grid_mask.embed_data(
                data = self.partition_mappings,
                masked_value = -1,
                ).astype(int)
        else:
            mappings = self.partition_mappings.reshape(
                self.grid_size,
                )

        return mappings

    def query_by_grid_indices(self, idxs):
        """Return the atom label for a grid site index"""
        
        mappings = self.embed().flatten()

        return mappings[idxs]

    def query_by_grid_points(self, gps):
        """Return the atom label for a grid point"""

        return self.query_by_grid_indices(
            idxs = self.gps2idx(gps),
            )

    def query_by_cart_points(self, sites_cart):
        """Dynamically calculate the nearest atom site to the input points"""

        return self.assign_sites(sites_cart)


class AssignSites(object):

    def __init__(self, ref_sites):
        self.ref_sites = ref_sites

    def __call__(self, query_sites):
        from scipy import spatial
        tree = spatial.KDTree(data=self.ref_sites)
        nn_dists, nn_groups = tree.query(query_sites.astype(np.float))
        return np.array(nn_groups, dtype=int)


class AssignSitesChunked(object):

    chunk_size = 10000

    def __init__(self,
        assign_sites,
        processor,
        ):

        self.assign_sites = assign_sites
        self.processor = processor

    def __call__(self, query_sites):

        results = self.processor(
            self.iter_jobs(
                query_sites = query_sites,
                )
            )

        mappings = np.empty(len(query_sites), dtype=int)

        for i in np.arange(0, len(query_sites), self.chunk_size):        
            mappings[i:i+self.chunk_size] = results.pop(0)

        return mappings

    def iter_jobs(self, query_sites):

        for i in np.arange(0, len(query_sites), self.chunk_size):
            yield self.processor.make_wrapper(
                func = self.assign_sites,
                query_sites = query_sites[i:i+self.chunk_size],
                )


class MakeVoronoiGridPartition(object):

    def __init__(self,
        processor = None,
        ):

        if (processor is None):
            from giant.processors import basic_processor
            processor = basic_processor

        self.processor = processor

    def __call__(self, 
        map_grid, 
        dataset,
        selection_string,
        mask_name,
        ):

        partition_sites_cart = self.get_sites(
            hierarchy = dataset.model.hierarchy,
            selection_string = selection_string,
            )

        partition = self.get_partition(
            map_grid = map_grid,
            partition_sites_cart = partition_sites_cart,
            mask_name = mask_name,
            )

        return partition

    def get_sites(self, 
        hierarchy, 
        selection_string,
        ):

        return np.array(
            hierarchy.atoms().select(
                hierarchy.atom_selection_cache().selection(selection_string),
                ).extract_xyz()
            )

    def get_partition(self, 
        map_grid, 
        partition_sites_cart,
        mask_name = None,
        ):

        if (mask_name is not None): 
            # Calculate for masked subset of grid points
            grid_mask = map_grid.masks[mask_name]
            query_sites_grid = grid_mask.get_mask_points()
        else: 
            # All points on the grid
            grid_mask = None
            query_sites_grid = map_grid.grid_points()

        # map sites_cart to grid
        partition_sites_grid = np.array(
            map_grid.cart2grid(
                sites_cart = partition_sites_cart, 
                origin_shift = True,
                )
            )

        sites_mappings = self.get_mappings(
            ref_sites = partition_sites_grid,
            query_sites = query_sites_grid,
            )

        partition = GridPartition(
            map_grid = map_grid,
            grid_mask = grid_mask,
            partition_sites = partition_sites_cart,
            partition_mappings = sites_mappings,
            )

        return partition

    def get_mappings(self, 
        ref_sites, 
        query_sites,
        ):

        assign_sites = AssignSitesChunked(
            assign_sites = AssignSites(
                ref_sites = ref_sites,
                ),
            processor = self.processor,
            )

        partition_mappings = assign_sites(
            query_sites = query_sites,
            )

        return partition_mappings

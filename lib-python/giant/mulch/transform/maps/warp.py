import numpy as np

from giant.mulch.transform.maps.utils import SampleableMap

def closest_point_within_tolerance(query, ref_list, tol):
    dist_sq = list((ref_list - query).dot())
    dmin_sq = min(dist_sq)
    if dmin_sq > tol ** 2:
        return -1
    return dist_sq.index(dmin_sq)

def get_interpolated_mapping_between_coordinates(query_list, ref_list, tol=0.01):
    """
    Take each of query_list and find the sites in ref_list (within tolerance).
    Missing sites will be interpolated to the closest neighbouring site.
    Return list of indices mapping the site in one to the closest site in the other.
    """

    import copy
    from scitbx.array_family import flex
    
    ref_list = flex.vec3_double(ref_list)
    tmp_idxs_q_to_r = [
        closest_point_within_tolerance(
            query = q, 
            ref_list = ref_list, 
            tol = tol,
            ) 
        for q in query_list
        ]
    assert tmp_idxs_q_to_r.count(-1) != len(tmp_idxs_q_to_r), 'no matching sites found between mappings'
    out_idxs_q_to_r = copy.copy(tmp_idxs_q_to_r)
    l = len(tmp_idxs_q_to_r)
    # Populate the missing values with the nearest value
    for i in range(l):
        d_i = 0
        while out_idxs_q_to_r[i] == -1:
            d_i += 1;
            p_i = i + d_i;
            n_i = i - d_i
            if (p_i < l) and (tmp_idxs_q_to_r[p_i] != -1):
                out_idxs_q_to_r[i] = out_idxs_q_to_r[p_i]
            elif (n_i >= 0) and (tmp_idxs_q_to_r[n_i] != -1):
                out_idxs_q_to_r[i] = out_idxs_q_to_r[n_i]

    return out_idxs_q_to_r
    

class VoronoiMapWarper(object):

    debug = True

    def __init__(self, map_grid):

        self.map_grid = map_grid

    def __call__(self, dataset, fft_map):

        map_grid = self.map_grid
        map_partition = map_grid.partition
        map_mask = map_partition.grid_mask

        # Create map object for extracting map values
        native_map_true = SampleableMap.from_fft_map(fft_map)

        # Extract the map sites from the grid partition
        point_mappings_grid = map_partition.partition_mappings

        if (self.debug is True): 
            assert (point_mappings_grid == -1).sum() == 0

        # Use map mask if present.
        if (map_mask is not None):
            sites_grid = map_mask.get_mask_points()
        else: 
            sites_grid = map_grid.get_points()

        # Map the grid sites to cartesian
        sites_cart = map_grid.grid2cart(
            sites_grid,
            origin_shift = True,
            )

        # Translate the grid partition mappings to the dataset alignment mappings 
        # (there might be some missing in this dataset)
        mappings_grid2dataset = np.array(
            get_interpolated_mapping_between_coordinates(
                query_list = map_partition.partition_sites,
                ref_list = dataset.model.alignment.reference_sites,
                tol = 0.01,
                )
            )

        point_mappings_dataset = mappings_grid2dataset[point_mappings_grid]

        if (self.debug is True): 
            assert sum(point_mappings_dataset == -1) == 0

        sites_cart_d = dataset.model.alignment.ref2nat(
            coordinates = sites_cart,
            mappings = point_mappings_dataset,
            )

        morphed_map_data = native_map_true.get_cart_values(sites_cart_d)

        return morphed_map_data


class ReferenceMapGridSampler(object): 

    def __init__(self, map_grid, mask_name=None):
        
        self.map_grid = map_grid
        self.mask_name = mask_name

    def __call__(self, dataset, fft_map): 

        native_map = SampleableMap.from_fft_map(fft_map)
        
        if (self.mask_name is not None):
            sites_grid = self.map_grid.masks[self.mask_name].get_mask_points()
        else:
            sites_grid = self.map_grid.grid_points()

        sites_cart = self.map_grid.grid2cart(
            sites_grid,
            origin_shift = True,
            )

        ## Removed as unnecessary!?
        # sites_cart_map = dataset.alignment.ref2nat(
        #     coordinates = sites_cart,
        #     )

        morphed_map_data = native_map.get_cart_values(sites_cart)

        return morphed_map_data

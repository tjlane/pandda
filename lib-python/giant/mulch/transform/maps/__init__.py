import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy as np

from giant.mulch.transform.maps.grid import GetWarpedMapGrid


class WrapperToIndexTuple(object):

    def __init__(self, i, func):
        """
        Call function and place result in appropriate row of output array
        """
        self.i = i
        self.func = func

    def __call__(self):

        data = self.func()
        
        return (self.i, data)


class DatasetMapLoaderAndWarperCallback(object):

    def __init__(self, load_map, warp_map, dataset, map_resolution):
        self.load_map = load_map
        self.warp_map = warp_map
        self.dataset = dataset
        self.map_resolution = map_resolution

    def __call__(self):
        dataset_map = self.load_map(
            dataset = self.dataset,
            resolution = self.map_resolution,
            )
        warped_map = self.warp_map(
            dataset = self.dataset,
            fft_map = dataset_map,
            )
        return warped_map


class MapArrayManager(object):

    def __init__(self, dataset_keys, map_data):

        self.dataset_keys = dataset_keys
        self.map_data = map_data

    def get(self, key):
        idx = self.dataset_keys.index(key)
        return self.map_data[idx]

    # def get_all(self):
    #     return self.map_data

    def as_dict(self, delete_internal_references=False):
        output_maps = {}
        for key in self.dataset_keys:
            output_maps[key] = self.get(key)
        # Delete reference to the map data -- call gc?
        if (delete_internal_references is True):
            del self.map_data
        return output_maps


class LoadWarpAndScaleMaps(object):

    def __init__(self,
        map_grid,
        reference_dataset,
        get_reference_map,
        dataset_map_getters,
        warp_reference_map,
        warp_dataset_map,
        scale_map_array,
        processor,
        ):

        self.map_grid = map_grid
        self.reference_dataset = reference_dataset
        self.get_reference_map = get_reference_map
        self.dataset_map_getters = dataset_map_getters
        self.warp_reference_map = warp_reference_map
        self.warp_dataset_map = warp_dataset_map
        self.scale_map_array = scale_map_array
        self.processor = processor

    def __call__(self,
        datasets,
        map_resolution,
        ):

        wrappers = []
        dataset_keys = sorted(datasets.keys())

        for i, dtag in enumerate(dataset_keys):

            get_warped_map_data = DatasetMapLoaderAndWarperCallback(
                dataset = datasets[dtag],
                load_map = self.dataset_map_getters[dtag],
                warp_map = self.warp_dataset_map,
                map_resolution = map_resolution,
                )

            wrappers.append(
                WrapperToIndexTuple(
                    i = i,
                    func = get_warped_map_data,
                    )
                )

        # Populate array
        results = self.processor(funcs=wrappers)

        # Create output array and populate
        map_array = np.empty(
            shape = (
                len(datasets),
                self.map_grid.masks["outer"].mask_size,
                ),
            dtype = float,
            )

        # Extract results
        for i, data in results:
            map_array[i] = data

        # Extract reference map and scale
        reference_map = self.get_reference_map(
            dataset = self.reference_dataset,
            resolution = map_resolution,
            )

        reference_map_warped = self.warp_reference_map(
            dataset = self.reference_dataset,
            fft_map = reference_map,
            )

        # Scale
        scaled_map = self.scale_map_array(
            map_array = map_array,
            reference_map_array = reference_map_warped,
            )

        return MapArrayManager(
            map_data = scaled_map,
            dataset_keys = dataset_keys,
            )

    def initialise(self, 
        datasets,
        ):

        from giant.mulch.xray.truncate_data import (
            CommonSetMillerArrayTruncator,
            )
        
        miller_arrays = [
            self.dataset_map_getters[dkey].get_miller_array.get_data(
                dataset = dataset,
                )
            for (dkey, dataset) in datasets.items()
        ]

        data_truncator = CommonSetMillerArrayTruncator(
            miller_arrays = miller_arrays,
            )

        # Replace ALL data truncators for simplicity
        for dkey, map_getter in self.dataset_map_getters.items():

            assert hasattr(map_getter.get_miller_array, "truncate_data")
            
            map_getter.get_miller_array.truncate_data = (
                data_truncator
                )

class ScaleMapsArrayInPlace(object):

    def __init__(self, map_grid):

        self.map_grid = map_grid

    def __call__(self, map_array, reference_map_array):

        assert map_array.shape[1] == len(reference_map_array), 'incompatible arrays'

        outer_mask_binary = self.map_grid.masks['outer'].get_mask_binary()
        inner_mask_binary = self.map_grid.masks['inner'].get_mask_binary()
        scale_mask_binary = inner_mask_binary[outer_mask_binary]

        # scale_mask = np.array(
        #     self.map_grid.index_on_other(
        #         query = self.map_grid.global_mask().inner_mask_indices(),
        #         other = self.map_grid.global_mask().outer_mask_indices(),
        #         )
        #     )

        for i in range(len(map_array)):
            scaled_map_data = self.scale_map_to_reference(
                ref_values = reference_map_array,
                values = map_array[i],
                mask = scale_mask_binary,
                )
            map_array[i] = scaled_map_data

        return map_array

    def scale_map_to_reference(self, ref_values, values, mask=None):
        """Scale given map to reference over the point_mask"""

        scale_ref_values = (
            ref_values[mask]
            if (mask is not None)
            else ref_values
            )
        
        scale_values = (
            values[mask]
            if (mask is not None)
            else values
            )

        # Fit ref = a*map + b
        # a, b
        scale, offset = np.polyfit(x=scale_values, y=scale_ref_values, deg=1)

        # output to full array
        scaled_values = (values*scale) + offset

        return scaled_values


class GetWarpedMapLoader(object):

    def __init__(self,
        get_map_grid,
        processor = None,
        ):

        if (processor is None):
            from giant.processors import basic_processor
            processor = basic_processor

        self.get_map_grid = get_map_grid
        self.processor = processor

    def __call__(self,
        reference_dataset,
        reference_map_getter,
        dataset_map_getters,
        #mask_name = XXX,
        ):

        from giant.mulch.transform.maps.warp import (
            ReferenceMapGridSampler, 
            VoronoiMapWarper,
            )

        # Load and partition map grid (slow -- do only once)
        map_grid = self.get_map_grid(
            reference_dataset = reference_dataset,
            )

        # reference_dataset.set_origin(
        #     map_grid.cart_origin(),
        #     )

        load_and_warp_maps = LoadWarpAndScaleMaps(
            map_grid = map_grid,
            reference_dataset = reference_dataset,
            get_reference_map = reference_map_getter,
            dataset_map_getters = dataset_map_getters,
            warp_reference_map = ReferenceMapGridSampler(
                map_grid = map_grid,
                mask_name = "outer",
                ),
            warp_dataset_map = VoronoiMapWarper(
                map_grid = map_grid,
                ),
            scale_map_array = ScaleMapsArrayInPlace(
                map_grid = map_grid,
                ),
            processor = self.processor,
            )

        return load_and_warp_maps

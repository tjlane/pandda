import giant.logs as lg
logger = lg.getLogger(__name__)


class DataGetterPrepper(object):
    """
    Class for loading miller arrays from datasets using pre-loaded functions.

    example usage: 
        get_miller_array = DataGetterPrepper(
            get_data = get_data_func, 
            truncate_data = truncate_data_func,
            scale_data = scale_data_func, 
            )
        miller_array = get_miller_array(dataset)
    """

    def __init__(self,
        get_data,
        truncate_data = None,
        scale_data = None,
        ):

        self.get_data = get_data
        self.truncate_data = truncate_data
        self.scale_data = scale_data

    def __call__(self, dataset):

        data = self.get_data(dataset)

        if (self.truncate_data is not None):
            data = self.truncate_data(data)

        if (self.scale_data is not None):
            data = self.scale_data(data)

        return data


class MapGetterPrepper(object):
    """
    Class for loading fft maps from datasets using pre-loaded functions.

    Valid map scalings: rmsd, volume, none

    example usage: 
        get_fft_map = MapGetterPrepper(
            get_miller_array = get_miller_array_func, 
            )
        fft_map = get_fft_map(dataset)
    """

    def __init__(self, 
        get_miller_array, 
        map_scaling = "rmsd", 
        map_resolution_factor = 0.25,
        ):

        if str(map_scaling) not in ["rmsd","sigma","volume","none","None"]:
            raise ValueError(
                'invalid map scaling: {}'.format(map_scaling)
                )

        self.get_miller_array = get_miller_array
        self.map_scaling = map_scaling
        self.map_resolution_factor = map_resolution_factor

    def __call__(self, dataset, resolution):

        import cctbx.maptbx

        miller_array = self.get_miller_array(dataset)

        flex_ed_map = miller_array.fft_map(
            resolution_factor = self.map_resolution_factor,
            d_min = resolution,
            symmetry_flags = cctbx.maptbx.use_space_group_symmetry,
            )

        if self.map_scaling in ["rmsd", "sigma"]:
            flex_ed_map = flex_ed_map.apply_sigma_scaling()
        elif self.map_scaling in ["volume"]: 
            flex_ed_map = flex_ed_map.apply_volume_scaling()
        elif self.map_scaling in ["none","None",None]:
            pass
        else:
            raise ValueError(
                'invalid map scaling: {}'.format(self.map_scaling)
                )

        return flex_ed_map

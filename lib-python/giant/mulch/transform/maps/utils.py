import numpy as np
from scitbx.array_family import flex


class SampleableMap(object):

    def __init__(self,
        map_data,
        unit_cell,
        origin = (0.,0.,0.),
        filter_out_of_bounds_points = False,
        ):

        if hasattr(map_data, "shape"):
            map_data = flex.double(map_data)

        if map_data.nd() != 3: 
            raise ValueError("map_data must be 3-dimensional data")

        self.map_data = map_data
        self.unit_cell = unit_cell
        self.origin = origin

        self.filter_out_of_bounds_points = filter_out_of_bounds_points

    @classmethod
    def from_fft_map(cls, 
        fft_map,
        ):

        return cls(
            map_data = fft_map.real_map_unpadded(),
            unit_cell = fft_map.unit_cell(),
            origin = (0.,0.,0.),
            filter_out_of_bounds_points = False,
            )

    @classmethod
    def from_map_grid(cls, 
        map_data, 
        map_grid, 
        filter_out_of_bounds_points,
        ):

        return cls(
            map_data = map_data,
            unit_cell = map_grid.unit_cell(),
            origin = map_grid.cart_origin(),
            filter_out_of_bounds_points = filter_out_of_bounds_points,
            )

    def get_cart_values(self, cart_points):

        from scitbx.array_family import flex

        cart_points = (cart_points - self.origin)

        frac_points = self.unit_cell.fractionalize(
            flex.vec3_double(cart_points),
            )

        return self.get_frac_values(frac_points)

    def get_frac_values(self, frac_points):

        if (self.filter_out_of_bounds_points is True):

            frac_values = np.zeros(len(frac_points), dtype=float)

            frac_points = np.array(frac_points)

            sel_bool = np.logical_and(
                (frac_points < 1.),
                (frac_points > 0.),
                ).astype(
                int
                ).prod(
                axis = -1,
                ).astype(
                bool
                )
            
            frac_points = frac_points[sel_bool]

            if len(frac_points) == 0:
                raise ValueError('no points selected')

        map_vals = list(map(
            self.map_data.eight_point_interpolation, 
            frac_points,
            ))

        if (self.filter_out_of_bounds_points is True):

            frac_values[sel_bool] = map_vals

            return frac_values

        else: 

            return np.array(map_vals)


# class RmsdNormaliseMapDataCrude:

#     def __init__(self):
#         pass

#     def __call__(self, map_data):
#         return (map_data - numpy.mean(map_data)) * (1.0/numpy.std(map_data))


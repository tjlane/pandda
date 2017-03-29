from scitbx.math import scale_curves
from scitbx.array_family import flex
from mmtbx.scaling import absolute_scaling

from giant.stats.optimisation import LeastSquaresFitting

class RelativeBfactorScaling(LeastSquaresFitting):

    def initialise_parameters(self):
        """Initialise starting simplex"""
        v0 = 0.0    # 0th order - offset
        v1 = 0.0    # 1st order - scale
        self.starting_simplex = [    flex.double([v0,    v1-20]),
                                     flex.double([v0,    v1+20]),
                                     flex.double([v0+0.1,v1])      ]
        return [v0,v1]

    def _scale(self, values, params):
        v0,v1 = params
        out = flex.exp(v0+v1*self.x_values)*values
        return out

class IsotropicScalingFactory(object):
    """Adapted from mmtbx.scaling.relative_wilson"""

    def __init__(self, reference_miller_array, d_star_min_sq=0.0625, d_star_max_sq=2.0, n_scale=2000):
        """
        Calculate scaling between any miller array and a reference miller array.
        This does not assume that the two arrays come from the same crystal.
        Data will only be scaled over the common resolution range of the crystal.
        For best results, the two miller arrays should have high completeness
        (although this is not checked).
        """

        # Number of points for scaling
        self._npoints = n_scale

        # Convert to intensities and store
        self.ref_miller = reference_miller_array.as_intensity_array()
        self.ref_kernel = self._kernel_normalisation(miller_array=self.ref_miller)

        # Store d_star limits and then update due to reference limits
        self.d_star_min_sq, self.d_star_max_sq = d_star_min_sq, d_star_max_sq
        self.d_star_min_sq, self.d_star_max_sq = self._common_d_star_sq_range(d_star_sq=self.ref_kernel.d_star_sq_array)

        return None

    def _common_d_star_sq_range(self, d_star_sq):
        """Find the common range over current limits and those of supplied array"""

        d_star_min_sq = max(self.d_star_min_sq, d_star_sq.min_max_mean().min)
        d_star_max_sq = min(self.d_star_max_sq, d_star_sq.min_max_mean().max)
        return d_star_min_sq, d_star_max_sq

    def _kernel_normalisation(self, miller_array, n_bins=45, n_term=17):

        return absolute_scaling.kernel_normalisation(miller_array=miller_array,
                                                     auto_kernel=True,
                                                     n_bins=n_bins,
                                                     n_term=n_term)

    def calculate_scaling(self, miller_array):
        """Calculate the scaling between two arrays"""

        # Convert to intensities and extract d_star_sq
        new_miller = miller_array.as_intensity_array()
        new_kernel = self._kernel_normalisation(miller_array=new_miller)
        # Calculate new range of d_star_sq
        d_star_sq_min, d_star_sq_max = self._common_d_star_sq_range(d_star_sq=new_kernel.d_star_sq_array)

        # Create interpolator for the two arrays (new and reference)
        interpolator = scale_curves.curve_interpolator(d_star_sq_min, d_star_sq_max, self._npoints)

        # Interpolate the two curves (use full range of the two array)
        new_itpl_d_star_sq, new_itpl_mean_I, dummy, dummy = interpolator.interpolate(
                                                                x_array=new_kernel.d_star_sq_array,
                                                                y_array=new_kernel.mean_I_array)
        ref_itpl_d_star_sq, ref_itpl_mean_I, dummy, dummy = interpolator.interpolate(
                                                                x_array=self.ref_kernel.d_star_sq_array,
                                                                y_array=self.ref_kernel.mean_I_array)

        # Run optimisation on the linear scaling
        lsc = RelativeBfactorScaling(x_values = interpolator.target_x,
                                     ref_y_values = ref_itpl_mean_I,
                                     mov_y_values = new_itpl_mean_I)

        return lsc


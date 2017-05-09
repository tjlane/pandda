import numpy

from scitbx.math import scale_curves, approx_equal_relatively
from scitbx.array_family import flex
from scitbx.python_utils.robust_statistics import percentile
from mmtbx.scaling import absolute_scaling

from giant.stats.optimisation import ExponentialScaling

class IsotropicBfactorScalingFactory(object):
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

    def calculate_scaling(self, miller_array, convergence_crit_perc=0.05, convergence_reject_perc=99, max_iter=20):
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

        # Initalise convergence loop - begin by scaling over all points
        selection = flex.bool(self._npoints, True)
        # Set initial scale factor to small value
        curr_b = 1e-6
        # Percent change between iterations - convergence when delta <convergence_criterion
        n_iter = 0
        while n_iter < max_iter:
            print 'ITER: '+str(n_iter)
            # Run optimisation on the linear scaling
            lsc = ExponentialScaling(x_values = interpolator.target_x,
                                     ref_values = ref_itpl_mean_I,
                                     scl_values = new_itpl_mean_I,
                                     weights = selection.as_double())
            # Calculate scaling B-factor
            lsc.scaling_b_factor = -0.5 * list(lsc.optimised_values)[0]
            # Break if fitted to 0
            if approx_equal_relatively(0.0, lsc.scaling_b_factor, 1e-6):
                break
            # Calculate percentage change
            print 'Curr/New: '+str(curr_b)+'\t'+str(lsc.scaling_b_factor)
            delta = abs((curr_b-lsc.scaling_b_factor)/curr_b)
            print 'Delta: '+str(delta)
            if delta < convergence_crit_perc: break
            # Update selection
            print 'Curr Selection Size: '+str(sum(selection))
            abs_diffs = flex.abs(flex.log(lsc.ref_values)-flex.log(lsc.out_values))
            sel_diffs = abs_diffs.select(selection)
            rej_val = numpy.percentile(sel_diffs, convergence_reject_perc)
            print 'Percentile: '+str(convergence_reject_perc)+'\t'+str(rej_val)
            selection.set_selected(abs_diffs>rej_val, False)
            print 'New Selection Size: '+str(sum(selection))
            # Update loop params
            curr_b = lsc.scaling_b_factor
            n_iter += 1

        lsc.unscaled_ln_rmsd = (flex.log(lsc.ref_values)-flex.log(lsc.scl_values)).norm()/(lsc.ref_values.size()**0.5)
        lsc.scaled_ln_rmsd   = (flex.log(lsc.ref_values)-flex.log(lsc.out_values)).norm()/(lsc.ref_values.size()**0.5)

        lsc.unscaled_ln_dev = flex.sum(flex.abs(flex.log(lsc.ref_values)-flex.log(lsc.scl_values)))
        lsc.scaled_ln_dev   = flex.sum(flex.abs(flex.log(lsc.ref_values)-flex.log(lsc.out_values)))

        return lsc


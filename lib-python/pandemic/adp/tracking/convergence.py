import giant.logs as lg
logger = lg.getLogger(__name__)

import copy, collections
import numpy
from libtbx import adopt_init_args, group_args

from pandemic.adp import constants


class PandemicConvergenceChecker(object):

    def __init__(self,
        parent,
        max_rmsd_b = None,
        max_delta_b = None,
        delta_b_window_frac = 0.05,
        delta_b_window_min = 5,
        eps_b = 0.01,
        ):
        """
        Applies multiple checks to find if the model has converged.
        max_rmsd_b: rmsd of model and target must be less than this for convergence
        max_delta_b: change in all B-factors over the last delta_n_cycles must be less than this,
            where delta_n_cycles is less given by max(delta_cycle_min, delta_cycle_frac*n_cycle_total)
        eps_b: mean B-factor of model must be greater than this
        """

        convergence_info = collections.OrderedDict()

        # Cycle number where the model became non-zero
        effective_n_start = None

        adopt_init_args(self, locals())

    def store_values(self,
        n_cycle,
        **kw_args
        ):

        for key, value in kw_args.items():
            key_dict = self.convergence_info.setdefault(key, collections.OrderedDict())
            key_dict[n_cycle] = value

    def update(self):

        # Extract data from parent
        n_cycle = int(self.parent.n_cycle)
        uij_current = self.parent.uij_history.get(n_cycle=n_cycle)

        assert uij_current is not None

        # Calculate iso-B (LEVEL, DATASET, ATOM)
        b_current = constants.EIGHTPISQ * uij_current[..., 0:3].mean(axis=-1)

        # Check if model is still zero and record cycle number if not
        non_zero_b = bool((b_current > self.eps_b).any())
        if (self.effective_n_start is None) and (non_zero_b is True):
            self.effective_n_start = n_cycle

        # Calculate delta B
        largest_delta_b = 0.0
        checking_from = None
        # First cycle where model is non-zero: compare to zero
        if (self.effective_n_start == n_cycle):
            largest_delta_b = b_current.max()
        # Compare over previous N cycles
        elif (self.effective_n_start is not None):
            # Calculate number of cycles to check the convergence over
            n_non_zero_cycles = max(1, n_cycle - self.effective_n_start)
            # Size of the window (minimum size or as fraction of cycles)
            n_check_start_delta = max(
                int(self.delta_b_window_min),
                int(numpy.ceil(self.delta_b_window_frac * n_non_zero_cycles)),
            )
            # Find the start cycle - must be at least the first cycle
            n_check_start = max(1, n_cycle - n_check_start_delta)
            checking_from = n_check_start

            # Calculate changes between selected previous cycles and current cycle
            for nn_cyc in range(n_check_start, n_cycle):

                # Extract the eigenvalues of the change over the last cycle
                uij_eigenvalues = self.parent.uij_history.get_delta_eigenvalues(
                    n_cycle_2 = n_cycle,
                    n_cycle_1 = nn_cyc,
                    )
                # Extract the largest eigenvalue (change)
                max_change = numpy.abs(uij_eigenvalues).mean(axis=-1).max()

                # maximum change since last cycle
                largest_delta_b = max(
                    largest_delta_b,
                    constants.EIGHTPISQ * max_change,
                    )

        self.store_values(
            n_cycle = n_cycle,
            non_zero = non_zero_b,
            mean_b = b_current.sum(axis=0).mean(),
            checking_from = checking_from,
            delta_b = largest_delta_b,
            rmsd_b = self.parent.table['rmsd'].iloc[-1],
        )

    def show(self):

        ci = self.convergence_info
        n_cyc = self.parent.n_cycle

        non_zero      = ci['non_zero'][n_cyc]
        mean_b        = ci['mean_b'][n_cyc]
        checking_from = ci['checking_from'][n_cyc]
        delta_b       = ci['delta_b'][n_cyc]
        rmsd_b        = ci['rmsd_b'][n_cyc]

        s = ''

        s += 'Level B-factor Changes:\n'
        b_eigenvalues = constants.EIGHTPISQ * self.parent.uij_history.get_delta_eigenvalues()
        for i_level, b_eigs in enumerate(b_eigenvalues):
            s += '> Level {level}: \n'.format(level=i_level+1)
            s += '\tMinimum: {minimum:+f}\n'.format(minimum=b_eigs.min())
            s += '\tAverage: {average:+f}\n'.format(average=b_eigs.mean())
            s += '\tMaximum: {maximum:+f}\n'.format(maximum=b_eigs.max())

        s += '\n'
        s += 'Convergence Checker Summary:\n'
        s += '> Model is approximately zero: {}\n'.format(
            'no' if non_zero else 'yes',
        )
        s += '> Average B-factor of model: {}\n'.format(
            mean_b,
        )
        s += '> Maximum change over recent cycles (since {}): {} (B-factor)\n'.format(
            checking_from,
            delta_b,
        )
        s += '    (cutoff for convergence: {})\n'.format(
            self.max_delta_b,
        )
        s += '> RMSD from target and fitted Uij : {} (B-factor)\n'.format(
            rmsd_b,
        )
        s += '    (cutoff for convergence: {})\n'.format(
            self.max_rmsd_b,
        )

        logger(s)

    def is_converged(self):

        # Initalise to false -- require only one success to set to true
        converged = False

        ci = self.convergence_info
        n_cyc = self.parent.n_cycle

        non_zero      = ci['non_zero'][n_cyc]
        mean_b        = ci['mean_b'][n_cyc]
        checking_from = ci['checking_from'][n_cyc]
        delta_b       = ci['delta_b'][n_cyc]
        rmsd_b        = ci['rmsd_b'][n_cyc]

        # Check if model is still zero -- never converged if this is the case
        if (bool(non_zero) is False):
            logger('Model is zero -- not converged')
            converged = False
            return converged

        if (self.max_rmsd_b is not None) and (rmsd_b < self.max_rmsd_b):
            logger('RMSD is below threshold -- converged')
            converged = True

        # Check if the change in B is less than tolerance
        if (self.max_delta_b is not None) and (delta_b < self.max_delta_b):
            logger('Delta B is below threshold -- converged')
            converged = True

        return converged

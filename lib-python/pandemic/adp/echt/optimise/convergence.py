import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy as np

from pandemic.adp.tracking.uijs import UijHistory


class EchtModelAmplitudeConvergenceChecker(object):

    def __init__(self,
        max_amplitude_change = 0.01,
        ):

        self.max_amplitude_change = max_amplitude_change
        self.previous_amplitudes = None

    def __call__(self, model_object):

        current_amplitudes = self.get_amplitudes(model_object)

        if self.previous_amplitudes is None: 

            self.previous_amplitudes = current_amplitudes

            return False

        diff_amplitudes = (
            current_amplitudes - self.previous_amplitudes
            )

        max_change = np.abs(diff_amplitudes).max()

        converged = (
            max_change < self.max_amplitude_change
            )

        logger.subheading(
            'Max amplitude change {max_amplitude_change:.4f}, model is {converged}converged'.format(
                max_amplitude_change = max_change,
                converged = ('' if converged else 'not '),
                )
            )

        self.previous_amplitudes = current_amplitudes

        return converged

    def get_amplitudes(self, model_object):

        array = np.concatenate([
            [
                [
                    t.amplitudes.get() 
                    for t in g.tls_parameters
                    ]
                for g in model_object.tls_objects[i]
                ]
            for i in range(model_object.n_tls_levels)
        ])

        array = array.flatten()

        return array

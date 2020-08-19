import numpy as np


def rms(values, axis=None, weights=None):
    """Calculate the rms of an array with optional weights"""

    if (weights is not None): 
        assert weights.shape == values.shape

    values_sq = np.power(values, 2)
    values_msq = np.average(values_sq, axis=axis, weights=weights)
    values_rms = np.sqrt(values_msq)

    return values_rms


import numpy as np

def finite(x):
    """Returns finite values of the given array"""
    return np.array(x)[np.isfinite(x)]

def round_no_fail(a, decimals=0):
    try:
        return np.round(a, decimals)
    except:
        return None

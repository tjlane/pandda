
import numpy

def finite(x):
    """Returns finite values of the given array"""
    return numpy.array(x)[numpy.isfinite(x)]

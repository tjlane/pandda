from bamboo.plot.radar import Radar

import numpy

def finite(x):
    """Returns finite values of the given array"""
    return x[numpy.isfinite(x)]

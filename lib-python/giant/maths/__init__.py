
import math
import numpy

DEG2RAD = math.pi / 180.
DEG2RADSQ = DEG2RAD ** 2.

RAD2DEG = 1. / DEG2RAD
RAD2DEGSQ = 1. / DEG2RADSQ

EIGHTPISQ = 8. * math.pi * math.pi

def finite(x):
    """Returns finite values of the given array"""
    return numpy.array(x)[numpy.isfinite(x)]

def round_no_fail(a, decimals=0):
    try:
        return numpy.round(a, decimals)
    except:
        return None

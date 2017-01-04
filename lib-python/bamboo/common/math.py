import numpy

def round_no_fail(a, decimals=0):
    try:    return numpy.round(a, decimals)
    except: return None

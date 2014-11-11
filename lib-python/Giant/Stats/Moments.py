import numpy

def skew(data):
    """Calculate the third moment of the data"""

    data_len = 1.0*len(data)
    data_mean = numpy.mean(data)
    data_cubed = [(x-data_mean)**3 for x in data]

    return sum(data_cubed)/(data_len)

def kurtosis(data):
    """Calculate the fourth moment of the data"""

    data_len = 1.0*len(data)
    data_mean = numpy.mean(data)
    data_quad = [(x-data_mean)**4 for x in data]

    return sum(data_quad)/(data_len)


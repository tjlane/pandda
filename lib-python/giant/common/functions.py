from libtbx import adopt_init_args

import numpy


class _ScalingFunction(object):

    def __call__(self, value):
        return self.map(value)

    def map(self, value):
        raise Exception('Not implemented!')

    def set(self, **kw_args):
        """Set function parameters"""
        for k in kw_args:
            assert k in self.__dict__
            self.__setattr__(k, kw_args[k])


class Linear(_ScalingFunction):
    name = 'Linear'

    def __init__(self, y_scale=1.0, y_intercept=0.0):
        adopt_init_args(self, locals())

    def map(self, x):
        return self.y_scale * x + self.y_intercept


class Sigmoid(_ScalingFunction):
    name = 'Sigmoid'

    def __init__(self, y_scale=1.0, x_width=1.0, x_offset=0.0, y_offset=0.0):
        adopt_init_args(self, locals())

    def rescale(self, x):
        return ( x - self.x_offset ) / self.x_width

    def map(self, x):
        rescaled = self.rescale(x)
        # Exponential function cannot take values larger than X (and does not matter as always gives 1. over this value)
        rescaled = numpy.clip(rescaled, -700.0, 700.0)
        return self.y_offset + self.y_scale / ( numpy.exp( - rescaled ) + 1.0 )


class Logarithm(_ScalingFunction):
    name = 'Logartihmic'

    def __init__(self, y_scale=1.0, x_width=1.0, x_offset=0.0):
        adopt_init_args(self, locals())

    def rescale(self, x):
        return ( x - self.x_offset ) / self.x_width

    def map(self, x):
        return self.y_scale * numpy.log( self.rescale(x) )


class PiecewiseLinear(_ScalingFunction):
    name = 'Piecewise Linear'

    def __init__(self, y_scale=0.0, y_discontinuity=0.0, x_discontinuity_location=0.0):
        adopt_init_args(self, locals())

    def rescale(self, x):
        return ( x - self.x_discontinuity_location )

    def map(self, x):
        rescaled = self.rescale(x)
        return (self.y_scale*rescaled + self.y_discontinuity) * (rescaled > 0.0)


def get_function(form, **kw_args):
    """Get and initialise one of several mathematical functions"""

    if form == 'linear':
        keys = ['y_scale','y_intercept']
        func = Linear
    elif form == 'piecewise_linear':
        keys = ['y_scale','y_discontinuity','x_discontinuity_location']
        func = PiecewiseLinear
    elif form == 'sigmoid':
        keys = ['y_scale', 'x_width', 'x_offset', 'y_offset']
        func = Sigmoid
    elif form == 'logarithm':
        keys = ['y_scale', 'x_width', 'x_offset']
        func = Logarithm
    else:
        raise Exception('Function not found: {}'.format(form))

    # Extract relevant keys and initialise function class
    filt_kw_args = {k:kw_args[k] for k in keys if k in kw_args}
    return func(**filt_kw_args)



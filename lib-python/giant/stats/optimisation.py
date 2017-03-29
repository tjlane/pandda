from scitbx.array_family import flex
from scitbx import simplex


class LeastSquaresFitting(object):

    def __init__(self, x_values, ref_y_values, mov_y_values, weights=None):
        """Calculate the optimum scaling of mov -> ref for a target function"""
        self.x_values = flex.double(x_values)
        self.ref_y_values = flex.double(ref_y_values)
        self.mov_y_values = flex.double(mov_y_values)
        # Store or create weightings
        if weights is not None:
            self.weight_array = flex.double(weights)
        else:
            self.weight_array = flex.double(self.ref_y_values.size(), 1.0/self.ref_y_values.size())
        # Validate inputs
        assert self.x_values.size() == self.mov_y_values.size()
        assert self.x_values.size() == self.ref_y_values.size()
        assert self.x_values.size() == self.weight_array.size()
        # Initialise simplex
        self.starting_values = self.initialise_parameters()
        # Calculate scaling
        self.run(initial_simplex=self.starting_simplex)

    def run(self, initial_simplex):
        """Calculate scaling"""
        # Optimise the simplex
        self.optimised = simplex.simplex_opt(dimension = 2,
                                             matrix    = initial_simplex,
                                             evaluator = self)
        # Extract solution
        self.optimised_values = self.optimised.get_solution()
        return self.optimised_values

    def target(self, vector):
        """Target function for the simplex optimisation"""
        scaled = self.transform(values=self.mov_y_values, params=vector)
        diff = (scaled-self.ref_y_values)
        diff_sq = diff*diff
        result = flex.sum(self.weight_array*diff_sq)
        return result

    def transform(self, values, params=None):
        """Function defining how the fitting parameters are used to transform the input vectors"""
        if params is None: params=self.optimised_values
        values = flex.double(values)
        assert values.size() == self.x_values.size()
        return self._scale(values=values, params=params)

    def new_x_values(self, x_values):
        self.x_values = flex.double(x_values)
        self.ref_y_values = None
        self.mov_y_values = None
        self.weight_array = None

class LinearFitting(LeastSquaresFitting):

    def initialise_parameters(self):
        """Initialise starting simplex"""
        v0 = 0.0    # 0th order - offset
        v1 = 1.0    # 1st order - scale
        self.starting_simplex = [    flex.double([v0,    v1]),
                                     flex.double([v0,    v1+0.1]),
                                     flex.double([v0+0.1,v1])      ]
        return [v0,v1]

    def _scale(self, values, params):
        v0,v1 = params
        out = v0 + values*v1
        return out


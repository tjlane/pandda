from libtbx import adopt_init_args
from scitbx.array_family import flex
from pandemic.adp.echt.sanitise.decompose_u import get_positive_u_component


class SanitiseUijValue(object):

    def __init__(self,
        tolerance,
        eps,
        ):

        # Create a near-zero uij value
        eps_value = (eps, eps, eps, 0.0, 0.0, 0.0)
        null_value = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        # Create validation function
        from pandemic.adp.echt.validate.uij import ValidateUijValues
        validate = ValidateUijValues(uij_tol = tolerance)

        adopt_init_args(self, locals())

        self.validate_input()

    def validate_input(self):
        assert self.tolerance >= 0.0, 'uij tolerance cannot be negative: {}'.format(self.tolerance)
        assert self.eps >= 0.0, 'uij eps cannot be negative: {}'.format(self.eps)

    def __call__(self, u):

        eigenvalues = self.validate.eigenvalues(u)

        # If the Uij is "null" then return the eps value
        if (eigenvalues < self.eps).all_eq(True):
            return self.eps_value

        # If the eigenvalues are too negative, then return only the valid component
        if flex.min(eigenvalues) < self.tolerance:
            u = get_positive_u_component(u)

        return u

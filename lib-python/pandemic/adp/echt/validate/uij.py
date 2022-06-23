from libtbx import adopt_init_args
from scitbx.array_family import flex
from giant.structure.uij import sym_mat3_eigenvalues


class ValidateUijValues(object):

    def __init__(self,
            uij_tol,
            ):
        adopt_init_args(self, locals())
        self.validate_input()

    def validate_input(self):
        assert self.uij_tol >= 0.0, 'uij tolerance cannot be negative: {}'.format(self.uij_tol)

    def __call__(self,
            u,
            ):
        eigenvalues = self.eigenvalues(u)
        # invalid if has large negative eigenvalues, or negative average eigenvalue?
        is_valid = (
            # no large negative eigenvalues
            (eigenvalues < -self.uij_tol).all_eq(False) and \
            # average eigenvalue (B-factor) is positive
            (flex.mean(eigenvalues) >= 0.0)
            )
        return is_valid

    @staticmethod
    def eigenvalues(u):
        return sym_mat3_eigenvalues(tuple(u))

from libtbx import adopt_init_args


class ValidateTLSMode(object):

    def __init__(self,
            matrix_tol = -1,
            ):
        adopt_init_args(self, locals())

    def __call__(self,
            tls_mode,
            ):
        return tls_mode.matrices.is_valid(self.matrix_tol)


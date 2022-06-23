import operator
from libtbx import adopt_init_args
from scitbx.array_family import flex
from mmtbx.tls.utils import uij_eigenvalues


class TargetTerm_UijLeastSquares(object):

    def __init__(self):
        adopt_init_args(self, locals())

    def __call__(self,
            uij_deltas,
            ):
        uij_deltas_1d = uij_deltas.as_1d().as_double()
        lsq = flex.sum(uij_deltas_1d*uij_deltas_1d*self.uij_weights)
        return self.scale * lsq

    def set_uij_weights(self, uij_weights):
        self.uij_weights = uij_weights.as_1d().matrix_outer_product(flex.double(6, 1)).as_1d()
        # Divide through by total number of atoms
        self.scale = 1.0
        #self.scale = 1.0 / float(reduce(operator.mul, uij_weights.all()))

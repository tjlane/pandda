import collections
import numpy

from libtbx import adopt_init_args


class UijHistory(object):

    def __init__(self):
        uijs = collections.OrderedDict()
        adopt_init_args(self, locals())

    def n_max(self):
        return max(self.uijs.keys())

    def add(self, n_cycle, uijs):
        self.uijs[n_cycle] = uijs

    def get(self, n_cycle, default=None):
        return self.uijs.get(n_cycle, default)

    def get_delta(self, n_cycle_2=None, n_cycle_1=None, default=0.0):
        if n_cycle_2 is None:
            n_cycle_2 = self.n_max()
        if n_cycle_1 is None:
            n_cycle_1 = n_cycle_2 - 1
        return self.get(n_cycle_2, default=default) - self.get(n_cycle_1, default=default)

    def get_delta_eigenvalues(self, n_cycle_2=None, n_cycle_1=None, default=0.0):
        deltas = self.get_delta(n_cycle_2=n_cycle_2, n_cycle_1=n_cycle_1, default=default)
        eigenvalues = self._uij_eigenvalues(deltas)
        return eigenvalues

    def _uij_eigenvalues(self, u_values):
        from mmtbx.tls.utils import uij_eigenvalues
        from scitbx.array_family import flex
        sh = u_values.shape
        # Reshape for calculating eigenvalues
        sh1 = (numpy.product(sh[:-2]), sh[-2], 6)
        # Reshape back to input
        sh2 = sh[:-1] + (3, )
        eig_values = numpy.array(list(map(uij_eigenvalues, map(flex.sym_mat3_double, u_values.reshape(sh1)))))
        return eig_values.reshape(sh2)

import time
import numpy as np
import functools

from libtbx import adopt_init_args
from scitbx.array_family import flex
from giant.structure.uij import sym_mat3_eigenvalues
from mmtbx.tls.utils import uij_eigenvalues
from pandemic.adp.echt.sanitise.decompose_u import get_negative_u_component


class SanitiseTLSGroup(object):

    max_cycles = 10

    def __init__(self,
        matrix_eps,
        amplitude_eps,
        ):
        adopt_init_args(self, locals())
        self.validate_input()

    def validate_input(self):
        assert self.matrix_eps >= 0.0, 'matrix eps cannot be negative: {}'.format(self.matrix_eps)
        assert self.amplitude_eps >= 0.0, 'amplitude eps cannot be negative: {}'.format(self.amplitude_eps)

    def __call__(self, tls_group):
        """Fix a TLS group object"""

        for i, mode in enumerate(tls_group.tls_parameters):

            # Reset the matrices of null modes
            if not (mode.amplitudes.get() > self.amplitude_eps).all_eq(True):
                mode.matrices.reset()
                mode.matrices.set(values=(1.,1.,1.,0.,0.,0.), component_string='T')

            # Reset small amplitudes (leave larger ones alone - change this?)
            if not (mode.amplitudes.get() >= self.amplitude_eps).all_eq(True):
                sel = (mode.amplitudes.get() < self.amplitude_eps).iselection()
                mode.amplitudes.set(values=[self.amplitude_eps]*sel.size(), selection=sel)

        # Make sure the produced Uijs have no negative eigenvalues
        self.remove_negative_eigenvalues(tls_group)

    def remove_negative_eigenvalues(self, tls_group):
        """Make sure that tls matrices do not generate negative eigenvalues"""

        for _ in range(self.max_cycles):

            unchanged = True

            for i_m, (amps, uijs) in enumerate(tls_group.uijs_unmultiplied()):
                # Sanity checks
                assert len(uijs) == tls_group.n_datasets
                assert uijs[0].size() == tls_group.n_atoms
                # Reshape to 1d
                uijs = np.array(uijs)
                uijs = uijs.reshape((np.product(uijs.shape[:-1]), 6))
                # Get eigenvalues for each dataset
                neg_comps = np.array(list(map(get_negative_u_component, uijs)))
                neg_magnitudes = np.abs(neg_comps[...,0:3].mean(axis=-1))
                max_mag = neg_magnitudes.max()

                # No negative magnitudes? Great!
                if max_mag <= 0.0:
                    continue

                # Something going to be changed
                unchanged = False

                # The largest negative component of a uij
                i_neg = np.where(neg_magnitudes == max_mag)[0][0]
                neg_uij = neg_comps[i_neg]

                # SUBTRACT this component from the T-matrix
                mode = tls_group.tls_parameters[i_m]
                new_t = [t-dt for t, dt in zip(mode.matrices.T, neg_uij)]
                mode.matrices.set(values=flex.double(new_t), component_string='T')

            if unchanged is True:
                break
        return

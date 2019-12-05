from libtbx import adopt_init_args
from scitbx.array_family import flex
from bamboo.common.logs import Log

import numpy


class RepairEchtModel:


    def __init__(self, 
        minimimum_uij_eigenvalue=0.0,
        verbose = False, 
        log = None,
        ):
        if log is None: log = Log()
        adopt_init_args(self, locals())

    def __call__(self, 
        model_object,
        ):

        for i_l, l_tls_objects in enumerate(model_object.tls_objects):
            for tls_group in l_tls_objects:
                self.repair_tls_group(tls_group=tls_group)

        model_object.adp_values = self.repair_adp_values(model_object.adp_values)

    def repair_tls_group(self, 
        tls_group,
        ):

        from mmtbx.tls.utils import uij_eigenvalues
        for i_m, (amps, uijs) in enumerate(tls_group.uijs_unmultiplied()):
            assert len(uijs) == tls_group.n_datasets
            assert uijs[0].size() == tls_group.n_atoms
            eigs = [uij_eigenvalues(u).as_double() for u in uijs]
            min_eig = min([flex.min(e) for e in eigs])
            if min_eig < self.minimimum_uij_eigenvalue:
                delta = abs(min_eig - self.minimimum_uij_eigenvalue)
                mode = tls_group.tls_parameters[i_m]
                current_t = list(mode.matrices.T)
                for i in range(3): current_t[i] += delta
                mode.matrices.set(values=tuple(current_t), component_string="T") 

    def repair_adp_values(self, 
        adp_values,
        ):
        from mmtbx.tls.utils import uij_eigenvalues
        eigs = numpy.array(uij_eigenvalues(adp_values))
        assert eigs.shape == adp_values.all() + (3,)
        min_eigs = eigs.min(axis=1)
        assert min_eigs.shape == adp_values.all()
        deltas = (min_eigs - self.minimimum_uij_eigenvalue)
        deltas[deltas > 0.0] = 0.0
        deltas = numpy.abs(deltas)
        deltas_triplets = deltas.reshape((deltas.size,1)).repeat(3, axis=1)
        zeroes_triplets = numpy.zeros((deltas.size, 3))
        deltas_sym = numpy.concatenate([deltas_triplets, zeroes_triplets], axis=1)
        deltas_flex = flex.sym_mat3_double(deltas_sym)
        return adp_values + deltas_flex
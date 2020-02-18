from libtbx import adopt_init_args
from scitbx.array_family import flex
from bamboo.common.logs import Log

import numpy

class SanitiseEchtModel:


    def __init__(self,
        tls_matrices_eps = 1e-6,
        tls_amplitudes_eps = 1e-6,
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
                self.sanitise_tls_group(tls_group=tls_group)

    def sanitise_tls_group(self,
        tls_group,
        ):

        for i_m, mode in enumerate(tls_group.tls_parameters):

            # Reset "null" groups
            if mode.is_null(
                    matrices_tolerance = self.tls_matrices_eps,
                    amplitudes_tolerance = self.tls_amplitudes_eps,
                    ):

                self.log('Resetting group: "{}" (Mode {})'.format(tls_group.label, i_m+1))

                mode.matrices.reset()
                mode.matrices.set(values=(1.,1.,1.,0.,0.,0.), component_string='T')
                mode.amplitudes.zero_values()

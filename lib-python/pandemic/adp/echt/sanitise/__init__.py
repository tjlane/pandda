import giant.logs as lg
logger = lg.getLogger(__name__)

from libtbx import adopt_init_args
from scitbx.array_family import flex

import numpy


class SanitiseEchtModel(object):

    def __init__(self,
        tls_parameters_dict,
        uij_parameters_dict,
        ):

        # Initialise sub-classes
        sanitise_tls_group = sanitise_uij_value = None
        if tls_parameters_dict is not None:
            from pandemic.adp.echt.sanitise.tls import SanitiseTLSGroup
            sanitise_tls_group = SanitiseTLSGroup(**tls_parameters_dict)
        if uij_parameters_dict is not None:
            from pandemic.adp.echt.sanitise.uij import SanitiseUijValue
            sanitise_uij_value = SanitiseUijValue(**uij_parameters_dict)

        adopt_init_args(self, locals())

    def __call__(self, model_object):
        for i_l, l_tls_objects in enumerate(model_object.tls_objects):
            self.sanitise_tls_groups(l_tls_objects)
        self.sanitise_uij_values(model_object.adp_values)

    def sanitise_tls_groups(self, tls_groups):
        """Sanitise a set of tls groups in place"""
        if self.sanitise_tls_group is None:
            return
        for tls_group in tls_groups:
            self.sanitise_tls_group(tls_group)

    def sanitise_uij_values(self, uij_values):
        """Sanitise a set of uij values in place"""
        if self.sanitise_uij_value is None:
            return
        for i, u in enumerate(uij_values):
            uij_values[i] = self.sanitise_uij_value(u)


import math
import numpy

import iotbx.pdb

from bamboo.common.logs import LogStream

from scitbx import matrix, linalg
from scitbx.array_family import flex
from mmtbx.tls import tlso, uaniso_from_tls_one_group, analysis

deg_to_rad_scale = math.pi/180

############

def extract_tls_from_pdb(pdb_file):
    ih = iotbx.pdb.hierarchy.input(pdb_file)
    tls_params = ih.input.extract_tls_params(ih.hierarchy)
    return tls_params

def uij_from_tls_vector_and_origin(xyz, tls_vector, origin):
    assert len(tls_vector) == 21
    tls_obj = tlso(t=tls_vector[0:6], l=tls_vector[6:12], s=tls_vector[12:21], origin=origin)
    return uij_from_tls_object(xyz=xyz, tls_obj=tls_obj)

def uij_from_tls_object(xyz, tls_obj):
    uij_tls = uaniso_from_tls_one_group(tlso       = tls_obj,
                                        sites_cart = flex.vec3_double(xyz),
                                        zeroize_trace = False)
    return numpy.array(uij_tls)

############

def get_t_l_s_from_vector(vals):
    assert len(vals) == 21
    return vals[0:6], vals[6:12], vals[12:21]

def tls_str_to_n_params(s, include_szz=True):
    assert not set(s).difference('TLS')
    return 6*('T' in s) + 6*('L' in s) + (8+include_szz)*('S' in s)

############

def validate_tls_object(tls_object, eps=1.e-6):
    return validate_tls_params(tls_vector=tls_object.t+tls_object.l+tls_object.s, eps=eps)

def validate_tls_params(tls_vector, eps=1.e-6, silent=False):
    """Validate a set of TLS matrices using the approach from phenix.tls_analysis. Returns True or False."""
    t, l, s = get_t_l_s_from_vector(tls_vector)
    try:
        checker = analysis.run(T = matrix.sym(sym_mat3=flex.double(t)),
                               L = matrix.sym(sym_mat3=flex.double(l))*(deg_to_rad_scale**2),
                               S = matrix.sqr(elems=flex.double(s))*deg_to_rad_scale,
                               log = LogStream(silent=silent),
                               eps = eps)
    except Exception as e:
        return e

    return True


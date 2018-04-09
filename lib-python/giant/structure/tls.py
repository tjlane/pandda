
import math
import numpy

import iotbx.pdb

from bamboo.common.logs import LogStream

from scitbx import matrix, linalg
from scitbx.array_family import flex
from mmtbx.tls import tlso, uaniso_from_tls_one_group
from mmtbx.tls.decompose import decompose_tls_matrices

deg_to_rad_scale = math.pi/180

############

def extract_tls_from_pdb(pdb_file):
    ih = iotbx.pdb.hierarchy.input(pdb_file)
    tls_params = ih.input.extract_tls_params(ih.hierarchy)
    return tls_params

def uij_from_tls_vector_and_origin(xyz, tls_vector, origin):
    assert len(tls_vector) == 21
    t,l,s = get_t_l_s_from_vector(tls_vector)
    tls_obj = tlso(t=t, l=l, s=s, origin=origin)
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


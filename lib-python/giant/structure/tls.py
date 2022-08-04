import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy

from giant.exceptions import Sorry, Failure

def phenix_find_tls_groups(filename=None, hierarchy=None):

    assert [filename, hierarchy].count(None) == 1, 'Must supply either filename or hierarchy'

    if hierarchy is not None:
        import tempfile
        tmp_file = tempfile.mktemp(suffix='.pdb')
        hierarchy.write_pdb_file(tmp_file)
        filename = tmp_file
    else:
        tmp_file = None

    from giant.dispatcher import Dispatcher
    cmd = Dispatcher('phenix.find_tls_groups')
    cmd.append_arg(filename)
    logger(cmd.as_string())
    result = cmd.run()

    if result['exitcode'] != 0:
        logger(str(result.stdout))
        logger(str(result.stderr))
        raise Failure('Failed to determine TLS groups: {}'.format(cmd.program))

    import re
    regex = re.compile(r"refinement\.refine\.adp \{([\s\S]*?)\}")
    tls_command = regex.findall(str(cmd.result.stdout))[0]
    tls_selections = [s.strip() for s in tls_command.split('tls =') if s.strip()]

    if (tmp_file is not None):
        import os
        os.remove(tmp_file) # called by tmp_file just in case...

    return tls_selections

############

def extract_tls_from_pdb(pdb_file):
    import iotbx.pdb
    ih = iotbx.pdb.hierarchy.input(pdb_file)
    tls_params = ih.input.extract_tls_params(ih.hierarchy)
    return tls_params

def uij_from_tls_vector_and_origin(xyz, tls_vector, origin):
    from mmtbx.tls import tlso
    assert len(tls_vector) == 21
    t,l,s = get_t_l_s_from_vector(tls_vector)
    tls_obj = tlso(t=t, l=l, s=s, origin=origin)
    return uij_from_tls_object(xyz=xyz, tls_obj=tls_obj)

def uij_from_tls_object(xyz, tls_obj):
    from mmtbx.tls import uaniso_from_tls_one_group
    from scitbx.array_family import flex
    uij_tls = uaniso_from_tls_one_group(
        tlso = tls_obj,
        sites_cart = flex.vec3_double(xyz),
        zeroize_trace = False,
    )
    return numpy.array(uij_tls)

def get_t_l_s_from_vector(vals):
    assert len(vals) == 21
    return vals[0:6], vals[6:12], vals[12:21]

def tls_str_to_n_params(s, include_szz=True):
    assert not set(s).difference('TLS')
    return 6*('T' in s) + 6*('L' in s) + (8+include_szz)*('S' in s)

def make_tls_isotropic(t, l, s, l_and_s_in_degrees=True):

    from scitbx import matrix, linalg
    from mmtbx.tls import decompose, analysis

    d = decompose.decompose_tls_matrices(
        T = t,
        L = l,
        S = s,
        l_and_s_in_degrees = l_and_s_in_degrees,
    )

    # Make vibrational amplitudes isotropic
    t = numpy.sqrt(numpy.mean(numpy.square(d.v_amplitudes)))

    tls = analysis.tls_from_motions(
        dx = d.l_amplitudes[0],
        dy = d.l_amplitudes[1],
        dz = d.l_amplitudes[2],
        l_x = matrix.rec(d.l_axis_directions[0], (3,1)),
        l_y = matrix.rec(d.l_axis_directions[1], (3,1)),
        l_z = matrix.rec(d.l_axis_directions[2], (3,1)),
        sx = d.s_amplitudes[0],
        sy = d.s_amplitudes[1],
        sz = d.s_amplitudes[2],
        tx = t, # set all the same
        ty = t,
        tz = t,
        v_x = matrix.rec((1,0,0), (3,1)), # don't need to provide actual axes as no longer defined
        v_y = matrix.rec((0,1,0), (3,1)),
        v_z = matrix.rec((0,0,1), (3,1)),
        w_M_lx = matrix.rec(d.l_axis_intersections[0], (3,1)),
        w_M_ly = matrix.rec(d.l_axis_intersections[1], (3,1)),
        w_M_lz = matrix.rec(d.l_axis_intersections[2], (3,1)),
    )

    T = tls.T_M
    L = tls.L_M
    S = tls.S_M

    # Should return in degrees also
    if l_and_s_in_degrees is True:
        from giant.common import RAD2DEG, RAD2DEGSQ
        T = T
        L = L * RAD2DEGSQ
        S = S * RAD2DEG

    return (T.as_sym_mat3(), L.as_sym_mat3(), S.as_mat3())


import copy, math
import numpy

from giant.structure.tls import tls_str_to_n_params
from scitbx import simplex, linalg, matrix


class _Simplex(object):
    _steps = {}

    def __init__(self, **kw_args):
        """Initialise Simplex"""
        self.set(**kw_args)

    def current_values(self):
        return copy.copy(self._steps)

    def get(self, key):
        assert self._steps.has_key(key), 'Key {} not one of {}'.format(key,self._steps.keys())
        return self._steps[key]

    def set(self, **kw_args):
        for k, v in kw_args.items():
            assert self._steps.has_key(k), 'Key {} not one of {}'.format(k,self._steps.keys())
            self._steps[k] = float(v)

    def get_simplex(self, start, variables, **kw_args):
        delta = self.get_deltas(variables, **kw_args)
        assert len(start)+1==len(delta) # e.g. 3d space (start) - 4 simplex points (delta)
        start_simplex = numpy.repeat([start], len(start)+1, axis=0)
        assert start_simplex.shape == delta.shape
        start_simplex += delta
        return start_simplex


class TLSSimplex(_Simplex):
    _steps = {
            'vibration_delta_min'   : 1e-3,  # minimum absolute change in matrix values
            'libration_delta_min'   : 1e-3,  # minimum absolute change in matrix values
            'amplitude_delta_min'   : 1e-6,  # minimum absolute change in amplitude values
            'vibration_delta_frac'  : 1e-2,  # fractional change in matrix values
            'libration_delta_frac'  : 1e-2,  # fractional change in matrix values
            'amplitude_delta_frac'  : 1e-2,  # fractional change in amplitude values
            }

    _DEG2RAD = math.pi/180.0
    _RAD2DEG = 1.0 / _DEG2RAD
    _RAD2DEGSQ = _RAD2DEG*_RAD2DEG

    _I_sqr = (1.,0.,0.,
              0.,1.,0.,
              0.,0.,1.)

    def get_deltas(self, variables, tls_matrices=None, tls_amplitudes=None): #, **kw_args

        p = variables
        # matrices, amplitudes, components, n_mdl, n_dst, tls_matrices, i_dst, i_mdl

        p_mdl = 0
        p_amp = 0
        if p.optimise_matrices is True:
            p_mdl += tls_str_to_n_params(p.components, include_szz=False)*p.n_mdl
        if p.optimise_amplitudes is True:
            p_amp += p.n_mdl*p.n_dst

        # Populate the delta array
        # Total number of parameters to vary
        n_all = p_mdl+p_amp
        # Begin indexing at 0 -- naturally
        i_this = 0
        deltas = numpy.zeros((n_all, n_all))
        # Deal with models
        if p.optimise_matrices is True:
            assert tls_matrices is not None, 'must provide tls_matrices when params.optimise_matrices is True'
            assert p.components is not None, 'must provide components when optimise_matrices is True'
            assert set('TLS').intersection(p.components), 'components must contain T, L or S'
            for i_mdl in p.i_mdl:
                # Extract decomposition
                m = tls_matrices[i_mdl]
                d = m.decompose()
                # If decomposition is invalid (fixing mode), obviously can't use the decomposition
                if not d.is_valid():
                    d = None
                # Generate base for each component
                for c in 'TLS':
                    # Skip if not selected
                    if c not in p.components:
                        continue
                    c_deltas = None
                    if (c=='T'):
                        n_c = 6
                        c_deltas = self.get_sym_mat3_deltas(
                                mat=m.T,
                                delta_frac=self.get('vibration_delta_frac'),
                                delta_min=self.get('vibration_delta_min'),
                                rot=(d.R_ML if d else None),
                                )
                    elif (c=='L'):
                        n_c = 6
                        c_deltas = self.get_sym_mat3_deltas(
                                mat=m.L,
                                delta_frac=self.get('libration_delta_frac'),
                                delta_min=self.get('libration_delta_min'),
                                rot=(d.R_ML if d else None),
                                )
                    elif (c=='S'):
                        n_c = 8
                        c_deltas = self.get_mat3_deltas(
                                mat=m.S,
                                delta_frac=self.get('libration_delta_frac'),
                                delta_min=self.get('libration_delta_min'),
                                rot=(d.R_ML if d else None),
                                )
                        # Must truncate this to 8 values (don't need Szz)
                        c_deltas = c_deltas[0:n_c,0:n_c]
                    else:
                        assert 0, 'Non-valid components given?'
                    # Check to see if any deltas extracted
                    assert c_deltas is not None
                    assert c_deltas.shape == (n_c, n_c), (c, n_c, c_deltas, c_deltas.shape, d)
                    # Apply deltas to simplex
                    for i, dels in enumerate(c_deltas):
                        deltas[i_this+i, i_this:i_this+n_c] = dels
                    i_this += n_c
        assert i_this == p_mdl
        # Deal with amplitudes
        if p.optimise_amplitudes is True:
            assert tls_amplitudes is not None, 'must provide amplitudes when params.optimise_amplitudes is True'
            assert len(tls_amplitudes) == p_amp
            amp_frac = self.get('amplitude_delta_frac')
            amp_min = self.get('amplitude_delta_min')
            for i_mdl in p.i_mdl:
                for i_dst in p.i_dst:
                    delt = abs(amp_frac * tls_amplitudes[i_mdl][i_dst])
                    deltas[i_this,i_this] = max(amp_min, delt)
                    i_this +=1
        assert i_this == n_all, 'something has gone wrong during simplex-delta generation...'
        # Prepend a list of zeros so that start point is also included
        deltas = numpy.append(numpy.zeros((1, n_all)), deltas, axis=0)
        assert deltas.shape == (n_all+1, n_all)
        return deltas

    def get_sym_mat3_deltas(self, mat, delta_frac, delta_min, rot=None):
        """Get deltas caused by changes to a symmetric matrix in a rotated frame"""

        if rot is None:
            rot = self._I_sqr

        # Create rotation matrix
        R = matrix.sqr(rot)
        R_t = R.transpose()

        # Transform input matrix
        mat_R = (R * matrix.sym(sym_mat3=mat) * R_t).as_sym_mat3()

        deltas = []
        n = 6

        mat_deltas = [
                (1., 1., 1., 0., 0., 0.),
                (1.,-1., 0., 0., 0., 0.),
                (0., 1.,-1., 0., 0., 0.),
                (0., 0., 0., 1., 0., 0.),
                (0., 0., 0., 0., 1., 0.),
                (0., 0., 0., 0., 0., 1.),
                ]

        for i in range(n):
            # Form the delta matrix in the rotated frame
            d_R = list(mat_deltas[i])
            for i_el in range(n):
                d_R[i_el] *= delta_min
            # Transform back
            d = (R_t * matrix.sym(sym_mat3=tuple(d_R)) * R).as_sym_mat3()
            # Append to output
            deltas.append(d)
        deltas = numpy.array(deltas)
        # remark: This check is necessary because providing a numpy float can cause
        # remark: the matrix multiplication to become a numpy array of matrices
        assert (deltas.shape == (n,n))
        return deltas

    def get_mat3_deltas(self, mat, delta_frac, delta_min, rot=None):
        """Get deltas caused by changes to a non-symmetric matrix in a rotated frame"""

        if rot is None:
            rot = self._I_sqr

        # Create rotation matrix
        R = matrix.sqr(rot)
        R_t = R.transpose()

        # Transform input matrix
        mat_R = (R * matrix.sqr(elems=mat) * R_t).as_mat3()

        deltas = []
        n = 9
        for i in range(n):
            # Generate the delta
            v = mat_R[i] * delta_frac
            if abs(v) < delta_min:
                v = delta_min
            # Form the delta matrix in the rotated frame
            d_R = [0.0,]*n
            d_R[i] = v
            # Transform back
            d = (R_t * matrix.sqr(tuple(d_R)) * R).as_mat3()
            # Append to output
            deltas.append(d)
        deltas = numpy.array(deltas)
        # remark: This check is necessary because providing a numpy float can cause
        # remark: the matrix multiplication to become a numpy array of matrices
        assert (deltas.shape == (n,n))
        return deltas

    def summary(self):
        s = 'Simplex optimisation deltas'
        for a in ['amplitude','vibration','libration']:
            s += '\n{} Fractional Change:   {}'.format(a.capitalize(), self.get(a+'_delta_frac'))
            s += '\n{} Minimum Change:   {}'.format(a.capitalize(), self.get(a+'_delta_min'))
        return s

class UijSimplex(_Simplex):
    _steps = {'uij_delta': 1e-6}

    def get_deltas(self, variables):
        n_uij = 6
        deltas = numpy.zeros((n_uij+1,n_uij))
        step = self.get('uij_delta')
        for i in range(n_uij):
            deltas[i+1,i] = step
        return deltas

    def summary(self):
        s = 'Simplex optimisation deltas'
        s += '\nUij Step Size:  {}'.format(self.get('uij'))
        return s


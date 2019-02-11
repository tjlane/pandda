import copy, math
import numpy

from giant.structure.tls import tls_str_to_n_params
from scitbx import simplex, linalg, matrix


class _Simplex(object):
    _steps = {}

    def __init__(self, **kw_args):
        """Initialise Simplex"""
        self.set_step_sizes(**kw_args)

    def current_values(self):
        return copy.copy(self._steps)

    def step_size(self, key):
        assert self._steps.has_key(key), 'Key {} not one of {}'.format(key,self._steps.keys())
        return self._steps[key]

    def set_step_sizes(self, **kw_args):
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
    _steps = {'vibration'   : 0.1,  # rms vibration magnitude change (angstrom)
              'libration'   : 1.0,  # rms libration magnitude change (degrees)
              'angle'       : 1.0,  # vector direction change (degrees)
              'amplitude'   : 0.1,  # multiplier (dimensionless)
             }

    _DEG2RAD = math.pi/180.0
    _RAD2DEG = 1.0 / _DEG2RAD
    _RAD2DEGSQ = _RAD2DEG*_RAD2DEG

    _I_sqr = (1.,0.,0.,
              0.,1.,0.,
              0.,0.,1.)
    _I_diag = (1.,1.,1.)

    _eigenvalue_change_basis = ((+1., 0., 0.,0.,0.,0.),
                                ( 0.,+1., 0.,0.,0.,0.),
                                ( 0., 0.,+1.,0.,0.,0.))
    _eigenvector_change_basis = (( 0., 1., 0.,
                                  -1., 0., 0.,
                                   0., 0., 0.),
                                 ( 0., 0., 1.,
                                   0., 0., 0.,
                                  -1., 0., 0.),
                                 ( 0., 0., 0.,
                                   0., 0., 1.,
                                   0.,-1., 0.))
    # 8 basis matrices with zero trace
    _s_change_basis = ((1.,0.,0.,
                        0.,0.,0.,
                        0.,0.,0.),
                       (0.,0.,0.,
                        0.,1.,0.,
                        0.,0.,0.),
                       (0.,1.,0.,
                        0.,0.,0.,
                        0.,0.,0.),
                       (0.,0.,1.,
                        0.,0.,0.,
                        0.,0.,0.),
                       (0.,0.,0.,
                        0.,0.,1.,
                        0.,0.,0.),
                       (0.,0.,0.,
                        1.,0.,0.,
                        0.,0.,0.),
                       (0.,0.,0.,
                        0.,0.,0.,
                        1.,0.,0.),
                       (0.,0.,0.,
                        0.,0.,0.,
                        0.,1.,0.))

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
        n_this = 0
        deltas = numpy.zeros((n_all, n_all))
        # Deal with models
        if p.optimise_matrices is True:
            assert tls_matrices is not None, 'must provide tls_matrices when params.optimise_matrices is True'
            assert p.components is not None, 'must provide components when optimise_matrices is True'
            assert set('TLS').intersection(p.components), 'components must contain T, L or S'
            for i_mdl in xrange(p.n_mdl):
                # Extract decomposition
                d = tls_matrices[i_mdl].decompose()
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
                        c_deltas = self.get_t_deltas(decomposition=d)
                    elif (c=='L'):
                        n_c = 6
                        c_deltas = self.get_l_deltas(decomposition=d)
                    elif (c=='S'):
                        n_c = 8
                        c_deltas = self.get_s_deltas(decomposition=d)
                    else:
                        assert 0, 'Non-valid components given?'
                    # Check to see if any deltas extracted
                    assert c_deltas is not None
                    assert c_deltas.shape == (n_c, n_c), (c, n_c, c_deltas, c_deltas.shape, d)
                    # Apply deltas to simplex
                    for i, dels in enumerate(c_deltas):
                        deltas[n_this+i, n_this:n_this+n_c] = dels
                    n_this += n_c
        # Deal with amplitudes
        if p.optimise_amplitudes is True:
            #TODO assert tls_ampltidues is not None, 'must provide amplitudes when params.optimise_amplitudes is True'
            # TODO make a percentage of the current amplitude!
            amp_step = self.step_size('amplitude')
            for i in range(n_this, n_this+p_amp):
                deltas[i,i] = amp_step
            n_this += p_amp
        assert n_this == n_all, 'something has gone wrong during simplex-delta generation...'
        # Prepend a list of zeros so that start point is also included
        deltas = numpy.append(numpy.zeros((1, n_all)), deltas, axis=0)
        assert deltas.shape == (n_all+1, n_all)
        return deltas

    def get_eigenvalue_deltas(self, R, scale):
        """Get the matrix contributions for a change in the eigenvalues of a basis"""
        assert scale > 0.0
        # Calculate the contributions to a matrix R_t * D * R, (where D is
        # diagonal) from a change to the diagonal elements of D
        R_t = R.transpose()
        result = []
        for dm_prime_vals in self._eigenvalue_change_basis:
            dm_prime = matrix.sym(sym_mat3=dm_prime_vals) * float(scale)
            dm = R_t*dm_prime*R
            result.append(dm)
        return result

    def get_eigenvector_deltas(self, R, eigenvalues, scale):
        """Get the matrix contributions for a change in the eigenvectors of a basis"""
        assert scale > 0.0
        # If all eigenvalues are the same size, slightly modify them to
        # become bigger/smaller, otherwise all of the output is zero
        if len(set(eigenvalues)) == 1:
            assert eigenvalues[0] != 0.0
            eigenvalues = list(eigenvalues)
            eigenvalues[0] *= 1.05
            eigenvalues[2] /= 1.05
        # Convert to diagonal matrix (in the eigenbasis defined by R)
        m_prime = matrix.diag(map(float, eigenvalues))
        R_t = R.transpose()
        result = []
        # Calculate the contributions to the matrix R_t * m_prime * R from
        # a change in basis by a small rotation around an eigen-axis
        for dr_prime_vals in self._eigenvector_change_basis:
            dr_prime = matrix.sqr(elems=dr_prime_vals) * float(scale)
            comp1 = dr_prime.transpose() * m_prime
            comp2 = comp1.transpose()
            comp3 = comp1 * dr_prime
            dm_prime = comp1+comp2+comp3
            dm = R_t*dm_prime*R
            result.append(dm)
        return result

    def get_t_deltas(self, decomposition):
        """Transform modifications to the vibrational components to the M-basis"""
        # If there is no rotation, act as if there is no decomposition
        if (decomposition is not None) and sum(decomposition.v_amplitudes) < 1.e-16:
            decomposition = None
        # Extract coordinate frame
        if decomposition is not None:
            R = decomposition.R_MV
            eigs = tuple(map(float, numpy.power(decomposition.v_amplitudes, 2)))
        # Else use simple basis
        else:
            R = self._I_sqr
            #eigs = tuple(self._I_diag)
            eigs = (self.step_size('vibration')**2.0,)*3
        assert len(R) == 9
        assert len(eigs) == 3
        R = matrix.sqr(R)
        # Calculate step size for eigenvalues (sq. angstroms) and vectors (degrees)
        eig_scale = self.step_size('vibration') ** 2.0
        vec_scale = self.step_size('angle')*self._DEG2RAD
        # First three modifications to the T-matrix are changes to vibrational eigenvalues
        d_eig_val = self.get_eigenvalue_deltas(R=R, scale=eig_scale)
        # Second three changes are rotations of the principal axes
        d_eig_vec = self.get_eigenvector_deltas(R=R, eigenvalues=eigs, scale=vec_scale)
        # Convert to sym_mat and return
        #deltas = [m.as_sym_mat3() for m in d_eig_val+d_eig_vec]
        deltas = []
        for m in d_eig_val+d_eig_vec:
            deltas.append(m.as_sym_mat3())
        deltas = numpy.array(deltas)
        # remark: This check is necessary because providing a numpy float can cause
        # remark: the matrix multiplication to become a numpy array of matrices
        assert (deltas.shape == (6,6))
        return deltas

    def get_l_deltas(self, decomposition):
        """Transform modifications to the librational components to the M-basis"""
        # If there is no rotation, act as if there is no decomposition
        if (decomposition is not None) and sum(decomposition.l_amplitudes) < 1.e-16:
            decomposition = None
        # Extract coordinate frame
        if decomposition is not None:
            R = decomposition.R_ML
            eigs = tuple(map(float, self._RAD2DEGSQ*numpy.power(decomposition.l_amplitudes, 2)))
        # Else use simple basis
        else:
            R = self._I_sqr
            #eigs = tuple(self._I_diag)
            eigs = (self.step_size('libration')**2.0,)*3
        assert len(R) == 9
        assert len(eigs) == 3
        R = matrix.sqr(map(float, R))
        # Calculate step size for eigenvalues (sq. degrees) and vectors (degrees)
        eig_scale = self.step_size('libration') ** 2.0
        vec_scale = self.step_size('angle')*self._DEG2RAD
        # First three modifications to the T-matrix are changes to vibrational eigenvalues
        d_eig_val = self.get_eigenvalue_deltas(R=R, scale=eig_scale)
        # Second three changes are rotations of the principal axes
        d_eig_vec = self.get_eigenvector_deltas(R=R, eigenvalues=eigs, scale=vec_scale)
        # Convert to sym_mat and return
        #deltas = [m.as_sym_mat3() for m in d_eig_val+d_eig_vec]
        deltas = []
        for m in d_eig_val+d_eig_vec:
            deltas.append(m.as_sym_mat3())
        deltas = numpy.array(deltas)
        # remark: This check is necessary because providing a numpy float can cause
        # remark: the matrix multiplication to become a numpy array of matrices
        assert (deltas.shape == (6,6))
        return deltas

    def get_s_deltas(self, decomposition):
        """Transform modifications to the librational components to the M-basis"""
        if (decomposition is not None) and sum(decomposition.v_amplitudes+decomposition.l_amplitudes) < 1.e-16:
            decomposition = None
        # Extract coordinate frame
        if decomposition is not None:
            R = decomposition.R_ML
        # Else use simple basis
        else:
            R = self._I_sqr
        R = matrix.sqr(R)
        R_t = R.transpose()
        # Calculate step size (degrees * angstroms)
        scale = self.step_size('libration') * self.step_size('vibration')
        deltas = []
        for dm_prime_vals in self._s_change_basis:
            # Apply scaling
            dm_prime = scale * matrix.sqr(dm_prime_vals)
            # Calculate contribution in proper coordinate frame
            dm = R_t*dm_prime*R
            # Only need the first 8 components of each matrix
            deltas.append(dm.as_mat3()[:8])
            assert len(deltas[-1]) == 8
        assert len(deltas) == 8
        return numpy.array(deltas)

    def summary(self):
        s = 'Simplex optimisation deltas'
        s += '\nVibration (RMS) Step Size: {}'.format(self.step_size('vibration'))
        s += '\nLibration (RMS) Step Size: {}'.format(self.step_size('libration'))
        s += '\nVector Change Step Size:   {}'.format(self.step_size('angle'))
        s += '\nTLS Amplitude Step Size:   {}'.format(self.step_size('amplitude'))
        return s

class UijSimplex(_Simplex):
    _steps = {'uij': 0.001}

    def get_deltas(self, variables):
        n_uij = 6
        deltas = numpy.zeros((n_uij+1,n_uij))
        step = self.step_size('uij')
        for i in range(n_uij):
            deltas[i+1,i] = step
        return deltas

    def summary(self):
        s = 'Simplex optimisation deltas'
        s += '\nUij Step Size:  {}'.format(self.step_size('uij'))
        return s

class AmplitudeSimplex(_Simplex):
    _steps = {
            'absolute' : 0.01,
            'relative' : 0.01,
            }

    def _group_hash(self, variables, i_level, i_group, i_mode, i_dataset):
        v_ = variables
        assert i_level < len(v_.n_groups_per_level)
        assert i_group < v_.n_groups_per_level[i_level]
        assert i_mode < v_.n_tls_per_group
        assert i_dataset < v_.n_dst
        n_other = sum([v_.n_var_per_group * n for n in v_.n_groups_per_level[:i_level]])
        n_this = i_group * v_.n_var_per_group + \
                i_mode * v_.n_dst + \
                i_dataset
        return n_other + n_this

    def _residual_hash(self, variables, i_atm):
        v_ = variables
        n_other = sum([v_.n_var_per_group * n for n in v_.n_groups_per_level])
        return n_other + i_atm

    def get_deltas(self, variables):

        v = variables

        n_tot = 0
        for l,g in v.level_group_pairs:
            if l != 'X':
                n_tot += v.n_tls * v.n_dst
            else:
                n_tot += 1

        deltas = numpy.zeros((n_tot+1,n_tot))
        step = self.step_size('relative')

        i = 0
        for l,g in v.level_group_pairs:
            if l=='X':
                mult = -1.
                deltas[i+1,i] = mult*step
                i += 1
            else:
                mult = (-1.)**((l-1)%2)
                for _ in range(v.n_tls * v.n_dst):
                    deltas[i+1,i] = mult*step
                    i += 1
        assert i == n_tot

        #for i in range(n_tot):
        #    deltas[i+1,i] = step

        return deltas

        #p_lvl = sum(v.n_groups_per_level) * v.n_tls_per_group * v.n_dst
        #p_res = v.n_atm
        #p_tot = p_lvl + p_res

        #deltas = numpy.zeros((p_tot+1,p_tot))

        #i_row = 1

        #l_hash = dict([(l.index, i) for i, l in enumerate(levels)])

        ## Iter datasets
        #for i_d in xrange(v.n_dst):
        #    # Iter modes
        #    for i_m in xrange(v.n_tls_per_group):

        #        # Iter levels
        #        for i_l1, l1 in enumerate(levels):

        #            print '============================='

        #            # For the first level, create independent "things"
        #            if (i_l1 == 0):
        #                print 'First level'
        #                # Create simplex point for each group
        #                for i_g1 in xrange(l1.n_groups()):
        #                    h = self._group_hash(
        #                            variables = v,
        #                            i_level   = i_l1,
        #                            i_group   = i_g1,
        #                            i_mode    = i_m,
        #                            i_dataset = i_d,
        #                            )
        #                    print i_row, i_l1, i_g1, i_m, i_d, h
        #                    deltas[i_row, h] = - self._steps.get('absolute')
        #                    i_row += 1

        #            print '-------------'
        #            print 'i_row, i_l1, i_g1, i_m, i_d, h'

        #            # Get map for the lower level
        #            l1_dict = group_hash.get(l1.index)

        #            # Iterate through groups in this level
        #            for i_g1, (g1_n, g1_sel, g1) in enumerate(l1):
        #                # Get the link information for this group in this level
        #                g1_dict = l1_dict.get(g1_n)
        #                if g1_dict is None: continue

        #                # Create positive contribution for group in lower level
        #                g1_h = self._group_hash(
        #                        variables = v,
        #                        i_level   = i_l1,
        #                        i_group   = i_g1,
        #                        i_mode    = i_m,
        #                        i_dataset = i_d,
        #                        )
        #                print '============================='
        #                print 'Lower level'
        #                print i_row, i_l1, i_g1, i_m, i_d, g1_h

        #                print '-------------'

        #                # Need to create len(g2_ns) simplex points
        #                idxs = []
        #                steps = []

        #                # Iterate through the levels connected to this group in this level
        #                # Reformat this into 1-d [(l,g),(l,g),...]
        #                for l2_n, g2_ns in g1_dict.items():
        #                    # Extract the hash of the group
        #                    if l2_n == 'X':
        #                        s = -self._steps.get('relative')
        #                        g2_hash = dict([(a_n, i_a) for i_a, (a_n, a_s, a) in enumerate(residual)])
        #                        for g2_n in g2_ns:
        #                            i_g2 = g2_hash[g2_n]
        #                            g2_h = self._residual_hash(
        #                                    variables = v,
        #                                    i_atm = i_g2,
        #                                    )
        #                            idxs.append(g2_h)
        #                            steps.append(s)
        #                    else:
        #                        s = -self._steps.get('absolute')
        #                        i_l2 = l_hash[l2_n]
        #                        g2_hash = dict([(g_n, i_g) for i_g, (g_n, g_s, g) in enumerate(levels[i_l2])])
        #                        for g2_n in g2_ns:
        #                            i_g2 = g2_hash[g2_n]
        #                            g2_h = self._group_hash(
        #                                    variables = v,
        #                                    i_level   = i_l2,
        #                                    i_group   = i_g2,
        #                                    i_mode    = i_m,
        #                                    i_dataset = i_d,
        #                                    )
        #                            idxs.append(g2_h)
        #                            steps.append(s)

        #                    try:
        #                        print 'Upper level ({})'.format(l2_n)
        #                        print i_row, i_l2, i_g2, i_m, i_d, ','.join(map(str,idxs))
        #                        deltas[i_row, idxs] = -s # NEGATIVE
        #                    except:
        #                        print 'Upper level (residual)'
        #                        print i_row, '-', i_g2, '-', '-', ','.join(map(str,idxs))
        #                        deltas[i_row, idxs] = -s # NEGATIVE

        #                # First point - level below goes up, level above goes down
        #                deltas[i_row, g1_h] = self._steps.get('absolute')   # POSITIVE
        #                deltas[i_row, idxs] = steps                         # NEGATIVE
        #                # Finished with this point -- increment
        #                i_row += 1

        #                # Create another N-1 points from internal fluctuations of the group
        #                for i_cut in range(1, len(idxs)):
        #                    deltas[i_row, idxs[:i_cut]] = steps[:i_cut]
        #                    deltas[i_row, idxs[i_cut:]] = [-val for val in steps[i_cut:]]
        #                    print deltas[i_row, numpy.where(deltas[i_row] != 0.0)[0]]
        #                    i_row += 1

        #assert i_row == len(deltas)

        #return deltas


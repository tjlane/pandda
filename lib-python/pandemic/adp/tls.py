import copy, operator
import numpy

from giant.structure.tls import uij_from_tls_vector_and_origin

############################################################################

def get_t_l_s_from_vector(vals):
    assert len(vals) == 21
    return vals[0:6], vals[6:12], vals[12:21]

def str_to_n_params(s, include_szz=True):
    assert not set(s).difference('TLS')
    return 6*('T' in s) + 6*('L' in s) + (8+include_szz)*('S' in s)

############################################################################

class TLSModel(object):

    _decimals = 3
    _n_prm = 21

    def __init__(self, values=None):
        self.values = numpy.array(values) if (values is not None) else numpy.zeros(self._n_prm)
        assert len(self.values) == self._n_prm

    def __add__(self, other):
        return self.add(other)

    def _str_to_idx(self, components, include_szz=True):
        """Map T-L-S to parameter indexes"""
        assert isinstance(components, str), 'not a string: {}'.format(components)
        assert not set(components).difference('TLS'), 'components can only contain letters "TLS"'
        idx = (range(00,06)                    if 'T' in components else []) + \
              (range(06,12)                    if 'L' in components else []) + \
              (range(12,20)+[20]*(include_szz) if 'S' in components else [])
              # One of the S parameters is not free so may be deselected with include_szz=False
        return idx

    def add(self, other):
        assert len(other.values) == len(self.values)
        self = self.copy()
        self.values += other.values
        return self

    def any(self, components='TLS', tol=1e-6):
        return (numpy.abs(self.get(components=components)) > tol).any()

    def copy(self):
        return copy.deepcopy(self)

    def get(self, components, include_szz=True):
        idx = self._str_to_idx(components=components, include_szz=include_szz)
        return self.values[idx]

    def multiply(self, amplitudes):
        assert len(amplitudes) == self._n_prm
        self = self.copy()
        self.values = amplitudes*self.values
        return self

    def n_params(self, free=True, non_zero=False):
        if non_zero is False:
            return self._n_prm - 1*(free)
        else:
            return self.any('T')*6 + self.any('L')*6 + self.any('S')*(9-int(free))
            #t,l,s = get_t_l_s_from_vector(vals=self.values)
            #return (t!=0.0).any()*(t.size) + (l!=0.0).any()*(l.size) + (s!=0.0).any()*(s.size-free)

    def reset(self, components):
        self.set(vals=0.0, components=components, include_szz=True)

    def round(self, decimals):
        """Round TLS values to a set number of decimal places"""
        self.values.round(decimals, out=self.values)

    def set(self, vals, components, include_szz=True):
        idx = self._str_to_idx(components=components, include_szz=include_szz)
        self.values[idx] = vals
        self.round(decimals=self._decimals)
        if (include_szz is False) and ('S' in components):
            self.set_szz_value_from_sxx_and_syy(target_trace=0.0)

    def set_szz_value_from_sxx_and_syy(self, target_trace=0.0):
        """Update Szz so that Sxx+Syy+Szz=target_trace"""
        # S stored as S11, S12, S13, S21, S22, S23, S31, S32, S33
        # So this is: Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz
        # So indexed:  12,  13,  14,  15,  16,  17,  18,  19,  20
        self.values[20] = target_trace - (self.values[12] + self.values[16])

    def summary(self):
        """Print summary of the TLS components"""
        r = '> TLS parameters'
        t,l,s = get_t_l_s_from_vector(vals=self.values)
        r += '\n\tT: '+', '.join(['{:8.3f}'.format(v) for v in t])
        r += '\n\tL: '+', '.join(['{:8.3f}'.format(v) for v in l])
        r += '\n\tS: '+', '.join(['{:8.3f}'.format(v) for v in s])
        return r

    def uij(self, xyz, origin):
        return uij_from_tls_vector_and_origin(xyz=xyz, tls_vector=self.values, origin=origin)

class _TLSAmplitudeSet(object):

    _model = None
    _n_amp = None
    _terms = None
    _terms_ordered = None

    _decimals = 3

    def __init__(self, n_dst=1):
        assert len(self._terms) == self._n_amp, 'length of amplitude names not equal to number of amplitudes'
        self.values = numpy.ones((n_dst, self._n_amp))
        self._n_prm = numpy.product(self.values.shape)
        self._terms_ordered = list(self._terms)

    def __getitem__(self, key):
        if isinstance(key, list) or isinstance(key, tuple):
            assert isinstance(key[0], int)
        else:
            assert isinstance(key, int)
            key = [key]
        return self.values[key]

    def _str_to_idx(self, components):
        """Map T-L-S string to amplitude parameter indexes"""
        if isinstance(components, str):
            components = [components]
        idxs = [self._terms.index(c) for c in components]
        assert len(idxs) > 0, 'invalid components {} (should be one of {})'.format(components, self._terms)
        return idxs

    def copy(self):
        return copy.deepcopy(self)

    def component_sets(self):
        """Get a list of all of the "modes" of this amplitude set (the set of matrices that each amplitude affects, e.g. [TLS] or [TS,LS])"""
        return list(self._terms_ordered)

    def expand(self, datasets=None):
        """Convert n-element vector into 21 element vector for TLS multiplication"""
        # Get the multipliers for each of the TLS matrices
        t_mult, l_mult, s_mult = self.t_l_s_multipliers(datasets)
        n_dst = len(datasets) if datasets is not None else len(self.values)
        # Expand the multipliers to 21-length parameters for multiplying onto the TLS model
        t_amps = numpy.repeat(t_mult, 6, axis=0).reshape((n_dst, 6))
        l_amps = numpy.repeat(l_mult, 6, axis=0).reshape((n_dst, 6))
        s_amps = numpy.repeat(s_mult, 9, axis=0).reshape((n_dst, 9))
        # Combine into one matrix
        exp_amps = numpy.concatenate([t_amps, l_amps, s_amps], axis=1)
        assert exp_amps.shape == (n_dst,21)
        return exp_amps

    def get(self, components, datasets=None):
        idx = self._str_to_idx(components)
        if datasets is not None:    return self.values[datasets, idx]
        else:                       return self.values[:, idx]

    def n_amplitudes(self):
        return self._n_amp

    def n_params(self, non_zero=False):
        if non_zero is False:
            return self._n_prm
        else:
            n_dst = len(self.values)
            return sum([(self.values[:,i]!=1.0).any()*(n_dst) for i in range(self._n_amp)])

    def normalise(self):
        raise Exception('Not implemented for base class!')

    def reset(self, components, datasets=None):
        self.set(vals=1.0, components=components, datasets=datasets)

    def round(self, decimals):
        """Round TLS values to a set number of decimal places"""
        self.values.round(decimals, out=self.values)

    def set(self, vals, components, datasets=None):
        idx = self._str_to_idx(components)
        if datasets is not None:    self.values[datasets, idx] = vals
        else:                       self.values[:, idx] = vals
        self.round(decimals=self._decimals)

    def set_component_sets_order(self, sequence):
        """Change the order of the component sets as returned by component_sets"""
        assert not set(self._terms).symmetric_difference(sequence), 'terms must contain each of {} (currently {})'.format(self._terms, sequence)
        self._terms_ordered = sequence

    def summary(self):
        """Print summary of the TLS amplitudes"""
        r = '> TLS amplitudes ({})'.format(self._model)
        for i, vals in enumerate(self.values):
            r += '\n\tDataset {:4}: '.format(i+1)+' '.join(['{:8.3f} ({})'.format(v,g) for v,g in zip(vals, self._terms)])
        return r

    def t_l_s_multipliers(self, datasets):
        raise Exception('Not implemented for base class!')

class SimpleTLSAmplitudeSet(_TLSAmplitudeSet):

    _model = 'one amplitude per TLS group'
    _terms = ['TLS']
    _n_amp = 1

    def normalise(self):
        """Normalise the amplitudes and return multipliers to be applied to the TLS model object"""
        # Calculate the average of all the amplitudes
        amp_mean = numpy.mean(self.get(components='TLS'))
        # Abort normalisation if average amplitude is approx zero
        if amp_mean < 1e-6:
            return (None, None, None)
        # Normalise amplitudes and return multipliers
        self.set(vals=self.get(components='TLS')*(1.0/amp_mean), components='TLS')
        # Multipliers are all the same under this model
        return amp_mean, amp_mean, amp_mean

    def t_l_s_multipliers(self, datasets):
        """Return T, L, and S matrix multipliers for each dataset"""
        amps = self.values
        if datasets is not None:
            amps = amps[datasets]
        n_dst = len(amps)
        assert amps.shape == (n_dst,self._n_amp)
        # Extract the only columns as the multipliers
        mult = amps[:,0]
        return (mult, mult, mult)

class TS_LS_TLSAmplitudeSet(_TLSAmplitudeSet):
    """Test-class for more complicated TLS-amplitude models"""

    _model = 'individual T- and L-amplitudes per TLS group'
    _terms = ['TS', 'LS']
    _n_amp = 2

    def normalise(self):
        raise Exception('Not yet implemented for this class!')

    def t_l_s_multipliers(self, datasets):
        """Return T, L, and S matrix multipliers for each dataset"""
        amps = self.values
        if datasets is not None:
            amps = amps[datasets]
        n_dst = len(amps)
        assert amps.shape == (n_dst,self._n_amp)
        t_mult = amps[:,0] * amps[:,0]
        l_mult = amps[:,1] * amps[:,1]
        s_mult = amps[:,0] * amps[:,1]
        return (t_mult, l_mult, s_mult)

class MultiDatasetTLSModel(object):

    _mdl_class = TLSModel
    _amp_class = SimpleTLSAmplitudeSet

    def __init__(self, n_dst=1, index=0, log=None):
#        if log is None: log = Log(verbose=False)
        self.log = log
        self.index = index
        self.n_dst = n_dst
        self.model = self._mdl_class()
        self.amplitudes = self._amp_class(n_dst=n_dst)

        # Block class function - should not be used after initialisation
        self.set_amplitude_model = None

    @classmethod
    def set_amplitude_model(cls, amplitude_class):
        cls._amp_class = amplitude_class
        return cls

    def component_sets(self):
        return self.amplitudes.component_sets()

    def expand(self, datasets=None):
        """Generate amplitude-multiplied TLS models for selected datasets"""
        exp_amps = self.amplitudes.expand(datasets=datasets)
        return [self.model.multiply(amplitudes=a) for a in exp_amps]

    def n_params(self, non_zero=False):
        return self.model.n_params(non_zero=non_zero) + self.amplitudes.n_params(non_zero=non_zero)

    def normalise(self):
        """Normalise the TLS amplitudes to an average of 1"""
        # Normalise amplitudes and extract multipliers for model
        multipliers = self.amplitudes.normalise()
        # Iterate through each matrix
        for cpt, mult in zip('TLS', multipliers):
            # Skip if multiplier is None
            if mult is None: continue
            self.log('> Normalising {} matrix of model {} -- applying multiplier of {:5.3f}'.format(cpt, self.index+1, mult))
            # Apply normalisation to TLS model
            self.model.set(vals=self.model.get(components=cpt)*(1.0*mult), components=cpt)

    def summary(self):
        """Print summary of the TLS model"""
        r = '> TLS Model {}'.format(self.index)
        r += '\n\t'+self.model.summary().replace('\n','\n\t')
        r += '\n\t'+self.amplitudes.summary().replace('\n','\n\t')
        return r

    def uij(self, xyzs, origins, datasets=None):
        """Generate Uijs from TLS models"""
        n_dst = len(datasets) if datasets is not None else self.n_dst
        assert len(xyzs) == len(origins) == n_dst
        exp_models = self.expand(datasets=datasets)
        assert len(exp_models) == n_dst
        uij = numpy.array([m.uij(xyz=xyzs[i], origin=origins[i]) for i,m in enumerate(exp_models)])
        return uij

    def zero_negative_amplitudes(self, error_tol=0.01):
        """Reset small negative values of amplitudes (raise errors for values less than -error_tol)"""

        # Make sure error_cut negative
        error_cut = -1.0*abs(error_tol)

        # Iterate through the amplitude sets (e.g. ['TLS'] or ['TS', 'LS'])
        for cpts in self.component_sets():
            self.log('Check for negative {} amplitudes of model {}'.format(cpts, self.index+1))

            amp_vals = self.amplitudes.get(components=cpts)
            if (amp_vals < error_cut).any():
                raise Failure('Negative amplitudes < {} obtained for model {} component {}.\n'.format(error_cut, self.index+1, cpts) + \
                              '\n'.join(['> Dataset {}, value: {}'.format(d+1, a) for d,a in enumerate(amp_vals) if a<error_cut]))
            reset_sel = (amp_vals < 0.0)
            if reset_sel.any():
                reset_dst = numpy.where(reset_sel)[0]
                self.log('Resetting negative amplitudes for {} datasets in model {}, component {} (values {})'.format(len(reset_dst), self.index+1, cpts, str(amp_vals[reset_dst])))
                self.amplitudes.reset(components=cpts, datasets=reset_dst)

    def zero_null_modes(self, mdl_tol=1e-9, amp_tol=1e-6):
        """Identify if either model or amplitudes have zero-values, and reset accordingly"""

        # Iterate through the amplitude sets (e.g. ['TLS'] or ['TS', 'LS'])
        for cpts in self.component_sets():
            self.log('Check for zero modes of {} of model {}'.format(cpts, self.index+1))

            # Check if all of the TLS model is approx zeroes
            mdl_vals = self.model.get(components=cpts)
            mdl_avge = numpy.mean(numpy.abs(mdl_vals))
            if (mdl_avge < mdl_tol):
                self.log('Zero-value model: resetting model {}, component(s) {}, values {}'.format(self.index+1, cpts, str(mdl_vals)))
                self.model.reset(components=cpts)
                self.amplitudes.reset(components=cpts)
                continue
            # Calculate the average amplitude
            amp_vals = self.amplitudes.get(components=cpts)
            amp_avge = numpy.mean(numpy.abs(amp_vals))
            if (amp_avge < amp_tol):
                self.log('Zero-value average amplitude: resetting model {}, component(s) {}, values {}'.format(self.index+1, cpts, str(amp_vals)))
                self.model.reset(components=cpts)
                self.amplitudes.reset(components=cpts)
                continue

class MultiDatasetTLSModelList(object):

    def __init__(self, n_mdl=1, n_dst=1, log=None):
#        if log is None: log = Log(verbose=False)
        self.log = log
        self.n_dst = n_dst
        self.n_mdl = n_mdl
        self._list = [MultiDatasetTLSModel(n_dst=self.n_dst, index=i, log=self.log) for i in xrange(self.n_mdl)]

    def __iter__(self):
        return iter(self._list)

    def get(self, index):
        if isinstance(index, int):
            if index > self.n_mdl: raise Failure('index {} too high (max {})'.format(index, self.n_mdl))
            return self._list[index]
        elif isinstance(index, list):
            return (self.get(index=i) for i in index)
        else:
            raise Failure('Must provide integer or list of integers')

    def expand(self, datasets=None):
        """Generate amplitude-multiplied TLS models for selected datasets and sum over all models"""
        n_dst = len(datasets) if (datasets is not None) else self.n_dst
        exp_models = [m.expand(datasets=datasets) for m in self]
        assert len(exp_models) == self.n_mdl
        assert len(exp_models[0]) == n_dst
        dst_models = [reduce(operator.add, m) for m in zip(*exp_models)]
        assert len(dst_models) == n_dst
        return dst_models

    def n_params(self, non_zero=False):
        return sum([m.n_params(non_zero=non_zero) for m in self])

    def normalise_amplitudes(self):
        """Normalise amplitudes model by model to average of 1"""
        for m in self:
            m.normalise()

    def reset_amplitudes(self, models, datasets=None, components=None):
        """Reset amplitudes for selected models to 1"""
        for m in self.get(models):
            cpts = components if (components is not None) else m.component_sets()
            m.amplitudes.reset(components=cpts, datasets=datasets)

    def uij(self, xyzs, origins, datasets=None):
        """Convert a set of parameter vectors to a set of uijs"""
        n_dst = len(datasets) if datasets is not None else self.n_dst
        n_atm = xyzs.shape[1]
        assert xyzs.shape == (n_dst, n_atm, 3)
        assert origins.shape == (n_dst, 3)
        uijs = numpy.sum([m.uij(xyzs=xyzs, origins=origins, datasets=datasets) for m in self], axis=0)
        assert uijs.shape == (n_dst, n_atm, 6)
        return uijs

    def summary(self):
        """Print summary of the TLS collection"""
        r = '> Multi-TLS Model summary'
        for m in self:
            r += '\n\t'+m.summary().replace('\n','\n\t')
        return r

    def zero_amplitudes(self, models, datasets=None, components=None):
        """Reset amplitudes for selected models to 1"""
        for m in self.get(models):
            cpts = components if (components is not None) else m.component_sets()
            m.amplitudes.set(vals=1.0, components=cpts, datasets=datasets)

    def zero_matrices(self, models, components=None):
        """Reset selected TLS matrices to 0"""
        for m in self.get(models):
            cpts = components if (components is not None) else m.component_sets()
            m.model.reset(components=cpts)

    def zero_negative_amplitudes(self, error_tol=0.01):
        """Reset negative ampltidues (raise error if amplitude < -1*error_tol)"""
        for i, m in enumerate(self):
            m.zero_negative_amplitudes(error_tol=error_tol)

    def zero_null_modes(self, mdl_tol=1e-9, amp_tol=1e-6):
        """Reset models and amplitudes that refine to zero"""
        for i, m in enumerate(self):
            m.zero_null_modes(mdl_tol=mdl_tol, amp_tol=amp_tol)


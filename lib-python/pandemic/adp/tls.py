import copy, operator
import numpy

from libtbx.utils import Sorry, Failure

from giant.structure.uij import uij_eigenvalues
from giant.structure.tls import uij_from_tls_vector_and_origin, get_t_l_s_from_vector
from mmtbx.tls.decompose import decompose_tls_matrices

############################################################################

class TLSModel(object):

    _dtype = float
    _decimals = 3
    _n_prm = 21

    _tolerance = 1.e-6

    def __init__(self, values=None, log=None):
#        if log is None: log = Log(verbose=False)
        self.log = log
        self.values = numpy.array(values, dtype=self._dtype) if (values is not None) else numpy.zeros(self._n_prm, dtype=self._dtype)
        assert len(self.values) == self._n_prm

        # Store precision, etc as instance variables so that not changed if class updated
        self._decimals = self._decimals
        self.set_precision = None
        self.set_tolerance = None

    def __add__(self, other):
        return self.add(other)

    def _str_to_idx(self, components, include_szz=True):
        """Map T-L-S to parameter indexes"""
        if components == 'all':
            components = 'TLS'
        assert isinstance(components, str), 'not a string: {}'.format(components)
        assert not set(components).difference('TLS'), 'components can only contain letters "TLS"'
        idx = (range(00,06)                    if 'T' in components else []) + \
              (range(06,12)                    if 'L' in components else []) + \
              (range(12,20)+[20]*(include_szz) if 'S' in components else [])
              # One of the S parameters is not free so may be deselected with include_szz=False
        return idx

    @classmethod
    def set_precision(cls, decimals):
        cls._decimals = decimals

    @classmethod
    def set_tolerance(cls, tolerance):
        cls._tolerance = tolerance

    def get_precision(self):
        return self._decimals

    def get_tolerance(self):
        return self._tolerance

    def add(self, other):
        assert len(other.values) == len(self.values)
        self = self.copy()
        self.values += other.values
        self.round()
        return self

    def any(self, components='TLS', tol=1.e-6):
        return (numpy.abs(self.get(components=components)) > tol).any()

    def copy(self):
        return copy.deepcopy(self)

    def get(self, components, include_szz=True):
        idx = self._str_to_idx(components=components, include_szz=include_szz)
        return self.values[idx]

    def decompose(self, tol=None):
        """Decompose TLS matrices into physical motions"""
        if tol is None: tol = self._tolerance
        T,L,S = get_t_l_s_from_vector(vals=self.values)
        try:
            result = decompose_tls_matrices(T=T, L=L, S=S,
                                            l_and_s_in_degrees=True,
                                            tol=tol)
        except Exception as e:
            self.log(e)
            return None
        return result

    def is_valid(self, tol=None):
        """Check whether the TLS matrices are physically valid"""
        if tol is None: tol = self._tolerance
        result = self.decompose(tol=tol)
        if result is None: return False
        return result.is_valid()

    def multiply(self, amplitudes):
        assert len(amplitudes) == self._n_prm
        self = self.copy()
        self.values = amplitudes*self.values
        self.round()
        return self

    def normalise(self, xyz, origin, target=1.0, tolerance=1e-16):
        """Given a set of xyz, scale matrices to give Uijs with approx dimensions of <target>"""
        # Extract Uijs for the supplied coordinates
        uijs = numpy.array(self.uij(xyz=xyz, origin=origin))
        # Extract the maximum axis length of each uij
        eigs = numpy.apply_along_axis(uij_eigenvalues, axis=1, arr=uijs)
        assert eigs.shape == xyz.shape
        maxs = numpy.max(eigs, axis=1)
        assert len(maxs) == len(xyz)
        # Calculate average of the maxima
        mean_max = numpy.mean(maxs)
        # Abort normalisation if max eigenvalue is approximately zero
        if mean_max < tolerance:
            return None
        # Calculate the scaling that modifies the mean to correspond to target
        mult = mean_max / target
        # Report
        self.log('> Normalise Uij values from TLS matrices to approx {}A -- dividing matrix values by {}'.format(target, mult))
        # Apply scaling and set szz value
        self.scale(multiplier=(1.0/mult))
        self.set_szz_value_from_sxx_and_syy()
        # Multipliers are all the same under this model
        return mult

    def n_params(self, free=True, non_zero=False):
        if non_zero is False:
            return self._n_prm - 1*(free)
        else:
            return self.any('T')*6 + self.any('L')*6 + self.any('S')*(9-int(free))
            #t,l,s = get_t_l_s_from_vector(vals=self.values)
            #return (t!=0.0).any()*(t.size) + (l!=0.0).any()*(l.size) + (s!=0.0).any()*(s.size-free)

    def reset(self, components):
        self.set(vals=0.0, components=components, include_szz=True)

    def round(self, decimals=None):
        """Round TLS values to a set number of decimal places"""
        if decimals is None: decimals = self._decimals
        self.values.round(decimals, out=self.values)

    def scale(self, multiplier):
        """Multiply the TLS matric values by a constant"""
        # Extract model values
        model_vals = self.get(components='TLS', include_szz=True)
        # Find out which values are zero before normalisation to ensure they are zero after normalisation
        approx_zero = numpy.abs(model_vals) < 10.**(-self.get_precision())
        # Apply normalisation to TLS matrices
        model_vals = model_vals*multiplier
        # Make all of the zero values are zero
        model_vals[approx_zero] = 0.0
        # Put the values back
        self.set(vals=model_vals, components='TLS')

    def set(self, vals, components, include_szz=True):
        idx = self._str_to_idx(components=components, include_szz=include_szz)
        self.values[idx] = vals
        if (include_szz is False) and ('S' in components):
            self.set_szz_value_from_sxx_and_syy(target_trace=0.0)
        self.round()

    def set_szz_value_from_sxx_and_syy(self, target_trace=0.0):
        """Update Szz so that Sxx+Syy+Szz=target_trace"""
        # S stored as S11, S12, S13, S21, S22, S23, S31, S32, S33
        # So this is: Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz
        # So indexed:  12,  13,  14,  15,  16,  17,  18,  19,  20
        self.values[20] = target_trace - (self.values[12] + self.values[16])

    def summary(self):
        """Print summary of the TLS components"""
        # Create format function for output
        format_func = numpy.vectorize('{{:{}.{}f}}'.format(self._decimals+5, self._decimals).format)
        # Generate summary string
        r = '> TLS parameters'
        t,l,s = get_t_l_s_from_vector(vals=self.values)
        r += '\n\tT: '+', '.join(format_func(t))
        r += '\n\tL: '+', '.join(format_func(l))
        r += '\n\tS: '+', '.join(format_func(s))
        return r

    def uij(self, xyz, origin):
        return uij_from_tls_vector_and_origin(xyz=xyz, tls_vector=self.values, origin=origin)

class TLS_AmplitudeSet(object):
    """Single-amplitude TLS-amplitude model"""

    _dtype = float
    _decimals = 2

    _model = 'single amplitude per TLS group'
    _terms = ['TLS']
    _n_amp = 1

    _description = """
    TLS amplitude model:
        One amplitude (a) per TLS model.
        All TLS matrices are coupled together.

        T -> a * T
        L -> a * L
        S -> a * S
    """

    def __init__(self, n_dst=1, log=None):
#        if log is None: log = Log(verbose=False)
        self.log = log
        assert len(self._terms) == self._n_amp, 'length of amplitude names not equal to number of amplitudes'
        self.values = numpy.ones((n_dst, self._n_amp), dtype=self._dtype)
        self._n_prm = numpy.product(self.values.shape)
        self._terms_ordered = list(self._terms)

        # Store precision, etc as instance variables so that not changed if class updated
        self._decimals = self._decimals
        self.set_precision = None

    @classmethod
    def description(cls):
        return cls._description

    @classmethod
    def set_precision(cls, decimals):
        cls._decimals = decimals

    def get_precision(self):
        return self._decimals

    def __getitem__(self, key):
        if isinstance(key, list) or isinstance(key, tuple):
            assert isinstance(key[0], int)
        else:
            assert isinstance(key, int)
            key = [key]
        return self.values[key]

    def _str_to_idx(self, components):
        """Map T-L-S string to amplitude parameter indexes"""
        # If all given -- reset all ampltidues
        if components == 'all':
            return range(self._n_amp)
        # If string, put in list
        if isinstance(components, str):
            components = [components]
        # Extract the indices for each term
        idxs = [self._terms.index(c) for c in components]
        assert len(idxs) > 0, 'invalid components {} (should be one of {})'.format(components, self._terms)
        return idxs

    def all_components(self):
        """Return all possible matrices that can be varied with this model (e.g. T or TL or TLS)"""
        return ''.join(sorted(set(''.join(self.component_sets())), key=lambda x: 'TLS'.index(x)))

    def component_sets(self):
        """Get a list of all of the "modes" of this amplitude set (the set of matrices that each amplitude affects, e.g. [TLS] or [TS,LS])"""
        return list(self._terms_ordered)

    def copy(self):
        return copy.deepcopy(self)

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

    def normalise(self, target=1.0):
        """Normalise the amplitudes and return multipliers to be applied to the TLS model object"""
        # Calculate the average of all the amplitudes
        amp_mean = numpy.round(numpy.mean(self.get(components='TLS')), self._decimals)
        # Abort normalisation if average amplitude is approx zero
        if amp_mean < 1e-6:
            return (None, None, None)
        # Modify the ampltiude mean to correspond to target
        mult = amp_mean / target
        # Normalise amplitudes and return multipliers
        self.log('> Normalising all amplitudes to average of {} -- dividing amplitudes by {}'.format(target, mult))
        self.scale(multiplier=(1.0/mult))
        # Multipliers are all the same under this model
        return (mult, mult, mult)

    def reset(self, components, datasets=None):
        self.set(vals=1.0, components=components, datasets=datasets)

    def round(self, decimals=None):
        """Round TLS values to a set number of decimal places"""
        if decimals is None: decimals = self._decimals
        self.values.round(decimals, out=self.values)

    def scale(self, multiplier):
        """Scale the amplitudes by a fixed constant"""
        self.log('> Scaling all amplitudes by constant -- multiplying amplitude values by {}'.format(multiplier, multiplier))
        for cpt in self.component_sets():
            old_vals = self.get(components=cpt)
            self.set(vals=old_vals*multiplier, components=cpt)

    def set(self, vals, components, datasets=None):
        idx = self._str_to_idx(components)
        if datasets is not None:    self.values[datasets, idx] = vals
        else:                       self.values[:, idx] = vals
        self.round()

    def set_component_sets_order(self, sequence):
        """Change the order of the component sets as returned by component_sets"""
        assert not set(self._terms).symmetric_difference(sequence), 'terms must contain each of {} (currently {})'.format(self._terms, sequence)
        self._terms_ordered = sequence

    def summary(self):
        """Print summary of the TLS amplitudes"""
        # Create format function for output
        format_func = numpy.vectorize('{{:{}.{}f}}'.format(self._decimals+5, self._decimals).format)
        # Generate output summary string
        r = '> TLS amplitudes ({})'.format(self._model)
        for i, vals in enumerate(format_func(self.values)):
            r += '\n\tDataset {:4}: '.format(i+1)+' '.join(['{} ({})'.format(v,g) for v,g in zip(vals, self._terms)])
        return r

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

class T_AmplitudeSet(TLS_AmplitudeSet):
    """Single-amplitude TLS-amplitude model"""

    _model = 'single translational amplitude per T-matrix. L-matrix and S-matrix set to zero.'
    _terms = ['T']
    _n_amp = 1

    _description = """
    T amplitude model:
        One amplitude (a) per TLS model.
        The amplitude of translational disorder (a) is varied.
        The L- and S-matrices are set to zero.

        T -> a * T
        L -> 0
        S -> 0
    """

    def normalise(self, target=1.0):
        """Normalise the amplitudes and return multipliers to be applied to the TLS model object"""
        # Calculate the average of all the amplitudes
        amp_mean = numpy.round(numpy.mean(self.get(components='T')), self._decimals)
        # Abort normalisation if average amplitude is approx zero
        if amp_mean < 1e-6:
            return (None, None, None)
        # Modify the ampltiude mean to correspond to target
        mult = amp_mean / target
        # Normalise amplitudes and return multipliers
        self.log('> Normalising T amplitudes to average of {} -- dividing amplitudes by {}'.format(target, mult))
        self.set(vals=self.get(components='T')*(1.0/mult), components='T')
        # Multipliers are all the same under this model
        return (mult, None, None)

    def t_l_s_multipliers(self, datasets):
        """Return T, L, and S matrix multipliers for each dataset"""
        amps = self.values
        if datasets is not None:
            amps = amps[datasets]
        n_dst = len(amps)
        assert amps.shape == (n_dst,self._n_amp)
        # Extract the only columns as the multipliers
        t_mult = amps[:,0]
        l_mult = s_mult = numpy.zeros(n_dst, dtype=self._dtype)
        return (t_mult, l_mult, s_mult)

class T_L_AmplitudeSet(TLS_AmplitudeSet):
    """Two-amplitude TLS-amplitude model"""

    _model = 'individual translational and librational amplitudes per TLS group, with no screw matrix.'
    _terms = ['T', 'L']
    _n_amp = 2

    _description = """
    T/L amplitude model:
        Two amplitudes (a,b) per TLS model.
        The amplitude of translational (a) and librational (b) motion are varied independently.
        The resulting T and L matrices are proportional to these amplitudes (a and b, respectively).
        The S matrix is set to zero.

        T -> a * T
        L -> b * L
        S -> 0
    """

    def normalise(self, target=1.0):
        """Normalise the amplitudes and return multipliers to be applied to the TLS model object"""
        # Calculate the average of all the amplitudes
        t_mean = numpy.round(numpy.mean(self.get(components='T')), self._decimals)
        l_mean = numpy.round(numpy.mean(self.get(components='L')), self._decimals)
        # Start with unitary multipliers
        t_mult = l_mult = None
        # Normalise T-S amplitude
        if t_mean > 1e-6:
            # Calculate normalisation for matrix values
            t_mult = t_mean / target
            # Normalise amplitudes and return multipliers
            self.log('> Normalising T amplitudes to average of {} -- dividing amplitudes by {}'.format(target, t_mult))
            self.set(vals=self.get(components='T')*(1.0/t_mult), components='T')
        # Normalise L-S amplitude
        if l_mean > 1e-6:
            # Calculate normalisation for matrix values
            l_mult = l_mean / target
            # Normalise amplitudes and return multipliers
            self.log('> Normalising L amplitudes to average of {} -- dividing amplitudes by {}'.format(target, l_mult))
            self.set(vals=self.get(components='L')*(1.0/l_mult), components='L')
        # Return the multipliers for the T-L-S matrices
        return (t_mult, l_mult, None)

    def t_l_s_multipliers(self, datasets):
        """Return T, L, and S matrix multipliers for each dataset"""
        amps = self.values
        if datasets is not None:
            amps = amps[datasets]
        n_dst = len(amps)
        assert amps.shape == (n_dst,self._n_amp)
        t_mult = amps[:,0]
        l_mult = amps[:,1]
        s_mult = numpy.zeros(n_dst, dtype=self._dtype)
        return (t_mult, l_mult, s_mult)

class MultiDatasetTLSModel(object):

    _mdl_classes = {'TLS': TLSModel}
    _amp_classes = {'TLS':      TLS_AmplitudeSet,
                    'T':        T_AmplitudeSet,
                    'T/L':      T_L_AmplitudeSet,
                   }

    def __init__(self, n_dst=1, index=0, model='TLS', amplitudes='TLS', log=None):
#        if log is None: log = Log(verbose=False)
        self.log = log
        self.index = index
        self.n_dst = n_dst
        self.model      = self._mdl_classes.get(model)(log=self.log)
        self.amplitudes = self._amp_classes.get(amplitudes)(n_dst=n_dst, log=self.log)

        # Block class function and delete attributes - should not be used after initialisation
        self.get_amplitude_model_class = None
        self._mdl_classes = None
        self._amp_classes = None

    @classmethod
    def get_amplitude_model_class(cls, amplitude_model):
        if amplitude_model not in cls._amp_classes:
            raise Sorry('Amplitude model "{}" not found. Possible choices: {}'.format(amplitude_model, sorted(cls._amp_classes.keys())))
        return cls._amp_classes.get(amplitude_model)

    def all_components(self):
        return self.amplitudes.all_components()

    def component_sets(self):
        return self.amplitudes.component_sets()

    def copy(self):
        return copy.deepcopy(self)

    def expand(self, datasets=None):
        """Generate amplitude-multiplied TLS models for selected datasets"""
        exp_amps = self.amplitudes.expand(datasets=datasets)
        return [self.model.multiply(amplitudes=a) for a in exp_amps]

    def n_params(self, non_zero=False):
        return self.model.n_params(non_zero=non_zero) + self.amplitudes.n_params(non_zero=non_zero)

    def normalise_by_amplitudes(self, target=1.0):
        """Normalise the TLS amplitudes to an average of 1"""
        # Normalise amplitudes and extract multipliers for model
        multipliers = self.amplitudes.normalise(target=target)
        # Iterate through each matrix
        for cpt, mult in zip('TLS', multipliers):
            # Skip if multiplier is None
            if mult is None: continue
            # Report
            fmt_mult = '{{:5.{}f}}'.format(self.amplitudes.get_precision()).format(mult)
            self.log('> Normalising {} matrix of model {} -- applying multiplier of {}'.format(cpt, self.index+1, fmt_mult))
            # Extract model values
            model_vals = self.model.get(components=cpt, include_szz=True)
            # Find out which values are zero before normalisation to ensure they are zero after normalisation
            approx_zero = numpy.abs(model_vals) < 10.**(-self.model.get_precision())
            # Apply normalisation to TLS model
            model_vals = (model_vals*mult)
            # Make all of the zero values zero
            model_vals[approx_zero] = 0.0
            # Put the valuse back
            self.model.set(vals=model_vals, components=cpt)
            self.model.set_szz_value_from_sxx_and_syy()
        self.log.bar(True, True)
        self.log(self.summary())
        self.log.bar(True, False)

    def normalise_by_matrices(self, xyzs, origins, target=1.0):
        """Normalise the TLS model uijs to have an average maximum eigenvalue of <target>"""
        # Validate input
        assert xyzs.ndim == 3
        assert origins.ndim == 2
        assert xyzs.shape[0] == origins.shape[0]
        assert xyzs.shape[2] == origins.shape[1] == 3
        # Convert xyzs and origins to one list of xyzs to pass to model
        n_dst, n_atm = xyzs.shape[:2]
        origins = origins.reshape((n_dst, 1, 3)).repeat(n_atm, axis=1)
        assert xyzs.shape == origins.shape
        xyzs = (xyzs-origins).reshape((n_dst*n_atm, 3))
        # Get the scaling for the model
        mult = self.model.normalise(xyz=xyzs, origin=(0,0,0), target=target)
        # Skip if multiplier is None
        if mult is None:
            self.log('Not able to normalise')
            return
        # Apply scaling to the amplitudes
        self.amplitudes.scale(multiplier=mult)
        if not self.model.is_valid():
            self.log('WARNING: TLS matrices are now invalid after normalisation.\nWARNING: These need to be fixed before optimisation can be continued')
        self.log.bar(True, True)
        self.log(self.summary())
        self.log.bar(True, False)

    def reset(self, components='all'):
        """Reset model back to zero-value matrices and unitary amplitudes"""
        # Model requires string (containing T-L-S letters only)
        if isinstance(components, list):
            self.model.reset(components=''.join(components))
        else:
            self.model.reset(components=components)
        # Amplitudes can reset either string or list
        self.amplitudes.reset(components=components)

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
            # Extract amplitudes
            amp_vals = self.amplitudes.get(components=cpts)
            # Error for any negative amplitudes below cutoff
            if (amp_vals < error_cut).any():
                raise Failure('Negative amplitudes < {} obtained for model {} component {}.\n'.format(error_cut, self.index+1, cpts) + \
                              '\n'.join(['> Dataset {}, value: {}'.format(d+1, a) for d,a in enumerate(amp_vals) if a<error_cut]))
            # Reset any other negative amplitudes
            reset_sel = (amp_vals < 0.0)
            if reset_sel.any():
                reset_dst = numpy.where(reset_sel)[0]
                self.log('Resetting negative amplitudes for {} datasets in model {}, component {} (values {})'.format(len(reset_dst), self.index+1, cpts, str(amp_vals[reset_dst])))
                self.amplitudes.set(vals=0.0, components=cpts, datasets=reset_dst)

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
                self.reset(components=cpts)
                continue
            # Calculate the average amplitude
            amp_vals = self.amplitudes.get(components=cpts)
            amp_avge = numpy.mean(numpy.abs(amp_vals))
            if (amp_avge < amp_tol):
                self.log('Zero-value average amplitude: resetting model {}, component(s) {}, values {}'.format(self.index+1, cpts, str(amp_vals)))
                self.reset(components=cpts)
                continue

class MultiDatasetTLSModelList(object):

    def __init__(self, n_mdl=1, n_dst=1, amplitudes='TLS', log=None):
#        if log is None: log = Log(verbose=False)
        self.log = log
        self.n_dst = n_dst
        self.n_mdl = n_mdl
        self._list = [MultiDatasetTLSModel(n_dst=self.n_dst, index=i, amplitudes=amplitudes, log=self.log) for i in xrange(self.n_mdl)]

    def __iter__(self):
        return iter(self._list)

    def any(self):
        """Return True if any model has non-zero values, else False"""
        for m in self:
            if m.model.any(): return True
        return False

    def copy(self):
        return copy.deepcopy(self)

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

    def normalise_by_amplitudes(self, target=1.0):
        """
        Normalise model by model so that the averages of
        the model amplitudes is <target>
        """
        self.log.subheading('Normalising TLS amplitudes')
        for m in self:
            m.normalise_by_amplitude(target=target)

    def normalise_by_matrices(self, xyzs, origins, target=1.0):
        """
        Normalise model by model so that the average of
        the maximum eigenvalue of each uij of the
        un-multiplied TLS matrices is <target>.
        """
        self.log.subheading('Normalising TLS matrices')
        for m in self:
            m.normalise_by_matrices(xyzs=xyzs, origins=origins, target=target)

    def reset(self, models=None, components=None):
        """Reset all matrices to 0 and amplitudes to 1"""
        models = models if (models is not None) else range(self.n_mdl)
        for m in self.get(models):
            cpts = components if (components is not None) else m.component_sets()
            m.reset(components=cpts)

    def reset_amplitudes(self, models=None, datasets=None, components=None):
        """Reset amplitudes for selected models to 1"""
        models = models if (models is not None) else range(self.n_mdl)
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

    def zero_amplitudes(self, models=None, datasets=None, components=None):
        """Reset amplitudes for selected models to 1"""
        for m in self.get(models):
            cpts = components if (components is not None) else m.component_sets()
            m.amplitudes.set(vals=0.0, components=cpts, datasets=datasets)

    def zero_matrices(self, models=None, components=None):
        """Reset selected TLS matrices to 0"""
        models = models if (models is not None) else range(self.n_mdl)
        for m in self.get(models):
            cpts = components if (components is not None) else 'TLS'
            m.model.reset(components=cpts)

    def zero_negative_amplitudes(self, error_tol=0.01):
        """Reset negative ampltidues (raise error if amplitude < -1*error_tol)"""
        for i, m in enumerate(self):
            m.zero_negative_amplitudes(error_tol=error_tol)

    def zero_null_modes(self, mdl_tol=1e-9, amp_tol=1e-6):
        """Reset models and amplitudes that refine to zero"""
        for i, m in enumerate(self):
            m.zero_null_modes(mdl_tol=mdl_tol, amp_tol=amp_tol)


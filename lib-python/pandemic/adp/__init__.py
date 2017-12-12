import os, sys, copy, traceback
import math, itertools, operator, collections

import scipy.cluster
import numpy, pandas

import libtbx.phil, libtbx.easy_mp, libtbx.easy_pickle
import iotbx.pdb
import mmtbx.tls.tools

import multiprocessing, multiprocessing.pool

from libtbx.utils import Sorry, Failure
from scitbx.array_family import flex
from scitbx import simplex, linalg

from mmtbx.secondary_structure import dssp

from bamboo.common import Meta, Info, ListStream
from bamboo.common.logs import Log, LogStream
from bamboo.common.path import easy_directory, rel_symlink
from bamboo.common.command import CommandManager

from giant.manager import Program
from giant.dataset import CrystallographicModel
from giant.structure.select import protein
from giant.structure.tls import uij_from_tls_vector_and_origin
from giant.structure.formatting import ShortLabeller, PhenixSelection, PymolSelection
from giant.structure.select import backbone, sidechains
from giant.structure.pymol import auto_residue_images, auto_chain_images, selection_images
from giant.xray.crystal import CrystalSummary
from giant.xray.refine import refine_phenix
from giant.xray.tls import phenix_find_tls_groups

from giant.jiffies import multi_table_ones

from pandemic.adp import html as pandemic_html

try:
    import matplotlib
    matplotlib.interactive(False)
    from matplotlib import pyplot, patches
    pyplot.switch_backend('agg')
    pyplot.style.use('ggplot')
    pyplot.interactive(0)
    pyplot.rc('font', family='monospace')
except Exception as e:
    print e

numpy.set_printoptions(threshold=numpy.nan)

from IPython import embed

EIGHT_PI_SQ = 8*math.pi*math.pi

############################################################################

DaemonicPool = multiprocessing.pool.Pool

class NonDaemonicProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

class NonDaemonicPool(multiprocessing.pool.Pool):
    Process = NonDaemonicProcess

############################################################################

PROGRAM = 'pandemic.adp'

DESCRIPTION = """
    Fit a consensus B-factor model to a series of datasets
"""

############################################################################

blank_arg_prepend = {'.pdb':'pdb=', '.cif':'cif=', '.pickle':'input.pickle='}

master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .help = "input pdb files - with isotropic/anisotropic b-factors"
        .multiple = True
        .type = str
    labelling = filename *foldername
        .type = choice
        .multiple = False
    pickle = None
        .type = path
        .multiple = False
}
output {
    out_dir = pandemic-adp
        .help = "output directory"
        .type = str
    pickle = True
        .type = bool
        .multiple = False
    diagnostics = False
        .help = "Write diagnostic graphs -- adds a significant amount of runtime"
        .type = bool
    pymol_images = False
        .help = "Write residue-by-residue images of the output B-factors"
        .type = bool
    group_graphs = False
        .help = "Write distribution graphs for each TLS group"
        .type = bool
}
levels {
    auto_levels = *chain *auto_group secondary_structure *residue *backbone *sidechain atom
        .type = choice(multi=True)
    custom_level
        .multiple = True
    {
        depth = None
            .help = "where to insert this level into the hierarchy? (after auto_levels have been determined)"
            .type = int
            .multiple = False
        selection = None
            .help = "list of selections that define groups for this level"
            .type = str
            .multiple = True
        mask_atoms = False
            .help = "remove atoms with high residual Uijs from this level?"
            .type = bool
            .multiple = False
    }
}
fitting {
    number_of_macro_cycles = 2
        .help = 'how many fitting cycles to run (over all levels)'
        .type = int
    number_of_micro_cycles = 3
        .help = 'how many fitting cycles to run (for each level)'
        .type = int
    number_of_modes_per_group = 1
        .help = 'how many TLS models to fit per group of atoms?'
        .type = int
    optimisation {
        max_datasets = None
            .help = 'takes up to this number of datasets for TLS parameter optimisation'
            .type = int
        max_resolution = None
            .help = 'resolution limit for dataset to be used for TLS optimisation'
            .type = float
    }
    penalties
        .help = 'penalties during optimisation. penalty function is 0 if p is 0, else (a + b*p), i.e. (p>0)*(a+b*p).'
    {
        invalid_tls_values
            .help = 'penalties for invalid TLS models. penalty p is number of T L or S matrices that is unphysical.'
            .expert_level = 3
        {
            a = 1e2
                .type = float
                .help = 'fixed penalty for positive penalty values'
            b = 1e2
                .type = float
                .help = 'multiplicative penalty for positive penalty values'
        }
        invalid_amplitudes
            .help = 'penalties for negative TLS amplitudes. penalty p is sum of amplitudes less than zero.'
            .expert_level = 3
        {
            a = 1e2
                .type = float
                .help = 'fixed penalty for positive penalty values'
            b = 1e2
                .type = float
                .help = 'multiplicative penalty for positive penalty values'
        }
        invalid_uij_values
            .help = 'penalties for invalid Uij values. penalty p is total number of negative eigenvalues for each Uij.'
            .expert_level = 3
        {
            a = 1e2
                .type = float
                .help = 'fixed penalty for positive penalty values'
            b = 1e3
                .type = float
                .help = 'multiplicative penalty for positive penalty values'
        }
        over_target_values
            .help = 'penalties when the fitted Uij is greater than the target Uij. penalty p is number of atoms with a Uij(fitted) > Uij(target). calculated as number of negative eigenvalues of tensor Uij(fitted)-Uij(target).'
        {
            a = 0.0
                .type = float
                .help = 'fixed penalty for positive penalty values'
            b = 1e3
                .type = float
                .help = 'multiplicative penalty for positive penalty values'
        }
    }
}
refine {
    refine_output_structures = False
        .help = "Refine the structures after fitting (coordinates and occupancies)"
        .type = bool
}
table_ones_options {
    include scope giant.jiffies.multi_table_ones.options_phil
}
settings {
    cpus = 1
        .type = int
        .multiple = False
    verbose = False
        .type = bool
    dry_run = False
        .type = bool
}
""", process_includes=True)

############################################################################

def get_t_l_s_from_vector(vals):
    return vals[0:6], vals[6:12], vals[12:21]

def str_to_n_params(s):
    return 6*('T' in s) + 6*('L' in s) + 9*('S' in s)

def rms(vals, axis=None):
    return numpy.sqrt(numpy.mean(numpy.power(vals,2), axis=axis))

def uij_to_b(uij):
    return EIGHT_PI_SQ*numpy.mean(uij[:,0:3],axis=1)

def default_secondary_structure_selections(hierarchy):
    return dssp.dssp(hierarchy).get_annotation().as_atom_selections()

############################################################################

def wrapper_plot_histograms(args):
    MultiDatasetUijPlots.multi_histogram(**args)

def wrapper_run(arg):
    return arg.run()

def wrapper_optimise(args):
    fitter = args
    try:
        #print '> Start optimisation: Fitter {} - Datasets {}'.format(fitter.label, ','.join(map(str,fitter._opt_dsts)))
        fitter._optimise(running_summary=False)
        #print '> Finished optimisation: Fitter {} - Datasets {}'.format(fitter.label, ','.join(map(str,fitter._opt_dsts)))
        return fitter
    except:
        return traceback.format_exc()

def wrapper_fit(args):
    fitter, kw_args = args
    try:
        fitter.optimise(**kw_args)
        return fitter
    except:
        return traceback.format_exc()

############################################################################

class TLSModel(object):

    _n_prm = 21

    def __init__(self, values=None):
        self.values = numpy.array(values) if (values is not None) else numpy.zeros(self._n_prm)
        assert len(self.values) == self._n_prm

    def __add__(self, other):
        return self.add(other)

    def _str_to_idx(self, components):
        """Map T-L-S to parameter indexes"""
        components = components.upper()
        assert not set(components).difference('TLS'), 'components can only contain letters "TLS"'
        idx = (range(00,06) if 'T' in components else []) + \
              (range(06,12) if 'L' in components else []) + \
              (range(12,21) if 'S' in components else [])
        return idx

    def copy(self):
        return copy.deepcopy(self)

    def get(self, components):
        idx = self._str_to_idx(components=components)
        return self.values[idx]

    def set(self, vals, components):
        idx = self._str_to_idx(components=components)
        self.values[idx] = vals

    def multiply(self, amplitudes):
        assert len(amplitudes) == self._n_prm
        self = self.copy()
        self.values = amplitudes*self.values
        return self

    def add(self, other):
        assert len(other.values) == len(self.values)
        self = self.copy()
        self.values += other.values
        return self

    def reset(self, components):
        self.set(vals=0.0, components=components)

    def uij(self, xyz, origin):
        return uij_from_tls_vector_and_origin(xyz=xyz, tls_vector=self.values, origin=origin)

    def n_params(self, non_zero=False):
        if non_zero is False:
            return self._n_prm
        else:
            t,l,s = get_t_l_s_from_vector(vals=self.values)
            return (t!=0.0).any()*(t.size) + (l!=0.0).any()*(l.size) + (s!=0.0).any()*(s.size)

    def summary(self):
        """Print summary of the TLS components"""
        r = '> TLS parameters'
        t,l,s = get_t_l_s_from_vector(vals=self.values)
        r += '\n\tT: '+', '.join(['{:8.3f}'.format(v) for v in t])
        r += '\n\tL: '+', '.join(['{:8.3f}'.format(v) for v in l])
        r += '\n\tS: '+', '.join(['{:8.3f}'.format(v) for v in s])
        return r

class TLSAmplitudeSet(object):

    def __init__(self, n_dst=1):
        self.values = numpy.ones((n_dst, 3))
        self._n_prm = numpy.product(self.values.shape)

    def __getitem__(self, key):
        if isinstance(key, list) or isinstance(key, tuple):
            assert isinstance(key[0], int)
        else:
            assert isinstance(key, int)
            key = [key]
        return self.values[key]

    def _str_to_idx(self, components):
        """Map T-L-S to parameter indexes"""
        components = components.upper()
        assert not set(components).difference('TLS'), 'components can only contain letters "TLS"'
        idx = ([0] if 'T' in components else []) + \
              ([1] if 'L' in components else []) + \
              ([2] if 'S' in components else [])
        return idx

    def copy(self):
        return copy.deepcopy(self)

    def get(self, components, datasets=None):
        idx = self._str_to_idx(components)
        if datasets is not None:    return self.values[datasets, idx]
        else:                       return self.values[:, idx]

    def set(self, vals, components, datasets=None):
        idx = self._str_to_idx(components)
        if datasets is not None:    self.values[datasets, idx] = vals
        else:                       self.values[:, idx] = vals

    def reset(self, components, datasets=None):
        self.set(vals=1.0, components=components, datasets=datasets)

    def expand(self, datasets=None):
        """Convert 3-element vector into 21 element vector for TLS multiplication"""
        amps = self.values
        if datasets is not None:
            amps = amps[datasets]
        n_dst = len(amps)
        assert amps.shape == (n_dst,3)
        t_amps = numpy.repeat(amps[:,0], 6, axis=0).reshape((n_dst, 6))
        l_amps = numpy.repeat(amps[:,1], 6, axis=0).reshape((n_dst, 6))
        s_amps = numpy.repeat(amps[:,2], 9, axis=0).reshape((n_dst, 9))
        exp_amps = numpy.concatenate([t_amps, l_amps, s_amps], axis=1)
        assert exp_amps.shape == (n_dst,21)
        return exp_amps

    def n_params(self, non_zero=False):
        if non_zero is False:
            return self._n_prm
        else:
            n_dst = len(self.values)
            return sum([(self.values[:,i]!=1.0).any()*(n_dst) for i in range(3)])

    def summary(self):
        """Print summary of the TLS amplitudes"""
        r = '> TLS amplitudes'
        for i, vals in enumerate(self.values):
            r += '\n\tDataset {:4}: '.format(i+1)+'{:8.3f} (T) {:8.3f} (L) {:8.3f} (S)'.format(*vals)
        return r

class MultiDatasetTLSModel(object):

    def __init__(self, n_dst=1, index=0, log=None):
#        if log is None: log = Log(verbose=False)
        self.log = log
        self.index = index
        self.n_dst = n_dst
        self.model = TLSModel()
        self.amplitudes = TLSAmplitudeSet(n_dst=n_dst)

    def n_params(self, non_zero=False):
        return self.model.n_params(non_zero=non_zero) + self.amplitudes.n_params(non_zero=non_zero)

    def normalise(self):
        """Normalise the TLS ampltidues to an average of 1"""
        for cpt in 'TLS':
            amp_mean = numpy.mean(self.amplitudes.get(components=cpt))
            if abs(amp_mean) < 1e-6: continue
            # Apply normalisation to amplitudes and TLS model
            self.model.set(vals=self.model.get(components=cpt)*(1.0*amp_mean), components=cpt)
            self.amplitudes.set(vals=self.amplitudes.get(components=cpt)*(1.0/amp_mean), components=cpt)

    def expand(self, datasets=None):
        """Generate amplitude-multiplied TLS models for selected datasets"""
        exp_amps = self.amplitudes.expand(datasets=datasets)
        return [self.model.multiply(amplitudes=a) for a in exp_amps]

    def uij(self, xyzs, origins, datasets=None):
        """Generate Uijs from TLS models"""
        n_dst = len(datasets) if datasets is not None else self.n_dst
        assert len(xyzs) == len(origins) == n_dst
        exp_models = self.expand(datasets=datasets)
        assert len(exp_models) == n_dst
        uij = numpy.array([m.uij(xyz=xyzs[i], origin=origins[i]) for i,m in enumerate(exp_models)])
        return uij

    def summary(self):
        """Print summary of the TLS model"""
        r = '> TLS Mode {}'.format(self.index)
        r += '\n\t'+self.model.summary().replace('\n','\n\t')
        r += '\n\t'+self.amplitudes.summary().replace('\n','\n\t')
        return r

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
        """Generate amplitude-multiplied TLS models for selected datasets and sum over all modes"""
        n_dst = len(datasets) if (datasets is not None) else self.n_dst
        exp_models = [m.expand(datasets=datasets) for m in self]
        assert len(exp_models) == self.n_mdl
        assert len(exp_models[0]) == n_dst
        dst_models = [reduce(operator.add, modes) for modes in zip(*exp_models)]
        assert len(dst_models) == n_dst
        return dst_models

    def n_params(self, non_zero=False):
        return sum([m.n_params(non_zero=non_zero) for m in self])

    def normalise_amplitudes(self):
        for m in self:
            m.normalise()

    def reset_models(self, components, modes):
        for m in self.get(modes):
            m.model.reset(components=components)

    def reset_amplitudes(self, components, modes, datasets=None):
        for m in self.get(modes):
            m.amplitudes.reset(components=components, datasets=datasets)

    def reset_zero_value_modes(self, mdl_tol=1e-3, amp_tol=1e-6):
        """Reset models and amplitudes that refine to zero"""
        for i, m in enumerate(self):
            for c in 'TLS':
                # Calculate the average model value
                mdl_vals = m.model.get(components=c)
                mdl_avge = numpy.mean(numpy.abs(mdl_vals))
                if (mdl_avge < mdl_tol):
                    self.log('Zero-value model: resetting model {}, component {}, values {}'.format(i+1, c, str(mdl_vals)))
                    m.model.reset(components=c)
                    m.amplitudes.reset(components=c)
                    continue
                # Calculate the average amplitude
                amp_vals = m.amplitudes.get(components=c)
                amp_avge = numpy.mean(numpy.abs(amp_vals))
                if (amp_avge < amp_tol):
                    self.log('Zero-value average amplitude: resetting model {}, component {}, values {}'.format(i+1, c, str(amp_vals)))
                    m.model.reset(components=c)
                    m.amplitudes.reset(components=c)
                    continue

    def reset_negative_amplitudes(self, error_tol=0.01):
        """Reset negative ampltidues (raise error if amplitude < -1*error_tol)"""
        error_cut = -1.0*abs(error_tol)
        for i, m in enumerate(self):
            for c in 'TLS':
                amp_vals = m.amplitudes.get(components=c)
                if (amp_vals < error_cut).any():
                    raise Failure('Negative amplitudes < {} obtained for model {} component {}.\n'.format(error_cut, i+1, c) + \
                                  '\n'.join(['> Dataset {}, value: {}'.format(d+1, a) for d,a in enumerate(amp_vals) if a<error_cut]))
                reset_sel = (amp_vals < -0.0)
                if reset_sel.any():
                    reset_dst = numpy.where(reset_sel)[0]
                    self.log('Resetting negative amplitudes for {} datasets in model {}, component {} (values {})'.format(len(reset_dst), i+1, c, str(amp_vals[reset_dst])))
                    m.amplitudes.reset(components=c, datasets=reset_dst)

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
        r = '> Multi-TLS Mode summary'
        for m in self:
            r += '\n\t'+m.summary().replace('\n','\n\t')
        return r

class UijPenalties(object):
    _tol = 1e-6
    _defaults = {'mdl_fixed' : 1e2,
                 'mdl_slope' : 1e2,
                 'amp_fixed' : 1e1,
                 'amp_slope' : 1e1,
                 'uij_fixed' : 1e2,
                 'uij_slope' : 1e3,
                 'ovr_fixed' : 1e0,
                 'ovr_slope' : 1e3}

    def __init__(self, **args):
        self.values = Info(self._defaults)
        self.set_weights(**args)

    def set_test_xyz(self, xyz, com):
        self._tst_xyz = xyz
        self._tst_com = com

    def set_weights(self, **args):
        """Set penalties for parameters to be invalid"""
        for k in args:
            assert k in self.values.__dict__.keys()
            self.values.__setattr__(k, args[k])
        return self.values.summary()

    def tls_params(self, values):
        """Return penalty for having unphysical TLS models"""
        assert len(values) == 21
        t,l,s = get_t_l_s_from_vector(vals=values)
        # Calculate eigenvalues of T+L matrices, and penalise negative eigenvalues (proportional to size of negative eigenvalues)
        t_penalty = self._sym_mat3_penalty(self._sym_mat3_eigenvalues(t))
        l_penalty = self._sym_mat3_penalty(self._sym_mat3_eigenvalues(l))
        # Only calculate S-penalties if S is non-zero
        if numpy.sum(numpy.abs(s)) > 0.0:
            # Calculate test Uijs for S-matrices
            s_uij_values = uij_from_tls_vector_and_origin(xyz=self._tst_xyz, tls_vector=[0.0]*12+list(s), origin=self._tst_com)
            # Only take the penalty of the largest invalid Uij (should be enough!)
            s_penalty = numpy.max([self._sym_mat3_penalty(self._sym_mat3_eigenvalues(uij)) for uij in s_uij_values])
        else:
            s_penalty = 0.0
        # Calculate and return total penalty
        penalty = numpy.sum([t_penalty, l_penalty, s_penalty])
        return self._standard_penalty_function(penalty       = penalty,
                                               fixed_penalty = self.values.mdl_fixed,
                                               slope_penalty = self.values.mdl_slope)

    def amplitudes(self, values):
        """Return penalty for having amplitudes with negative values"""
        penalty = abs(numpy.sum(values[values<-1.0*self._tol]))
        return self._standard_penalty_function(penalty       = penalty,
                                               fixed_penalty = self.values.amp_fixed,
                                               slope_penalty = self.values.amp_slope)

    def uij_valid(self, values):
        assert len(values) == 6
        penalty = self._standard_negative_eigenvalue_penalty(sym_mat3_vals=values)
        return self._standard_penalty_function(penalty       = penalty,
                                               fixed_penalty = self.values.uij_fixed,
                                               slope_penalty = self.values.uij_slope)

    def uij_size(self, fitted, target):
        """Add penalty for having fitted B-factors greater than observed"""
        if (self.values.ovr_fixed == 0.0) and (self.values.ovr_slope == 0.0):
            return 0.0
        values = [a-b for a,b in zip(target,fitted)]
        penalty = self._standard_negative_eigenvalue_penalty(sym_mat3_vals=values)
        return self._standard_penalty_function(penalty       = penalty,
                                               fixed_penalty = self.values.ovr_fixed,
                                               slope_penalty = self.values.ovr_slope)

    def _standard_penalty_function(self, penalty, fixed_penalty, slope_penalty):
        """Simple step function for non-zero penalty + slope to provide minimisation gradient"""
        assert penalty > -1.0*self._tol, 'penalties cannot be negative!'
        return (penalty>self._tol)*(fixed_penalty + slope_penalty*penalty)

    def _standard_negative_eigenvalue_penalty(self, sym_mat3_vals):
        eigenvalues = self._sym_mat3_eigenvalues(vals=sym_mat3_vals)
        return self._sym_mat3_penalty(eigenvalues=eigenvalues)

    def _sym_mat3_eigenvalues(self, vals):
        assert len(vals) == 6
        return linalg.eigensystem_real_symmetric(vals).values()

    def _sym_mat3_penalty(self, eigenvalues):
        """Calculate penalty for eigenvalues"""
        return -1.0*flex.sum(eigenvalues.select(eigenvalues<(-1.0*self._tol)))

class _Simplex(object):
    def get_simplex(self, start, **kw_args):
        delta = self.get_deltas(**kw_args)
        assert len(start)==len(delta)
        start_simplex = numpy.repeat([start], len(start)+1, axis=0)
        for i in xrange(len(start)):
            start_simplex[i+1][i] += delta[i]
        return start_simplex

class TLSSimplex(_Simplex):
    _del_mdl = 0.05
    _del_amp = 0.05

    def __init__(self, mdl_delta=None, amp_delta=None):
        """Initialise TLS Simplex"""
        self.set_deltas(mdl_delta=mdl_delta, amp_delta=amp_delta)

    def get_deltas(self, model, amplitudes, components, n_cpt, n_mdl, n_dst, **kw_args):
        p_mdl = 0
        p_amp = 0
        if model is True:
            p_mdl += str_to_n_params(components)*n_mdl
        if amplitudes is True:
            p_mdl += n_mdl*n_cpt*n_dst
        return numpy.array((self._del_mdl,)*p_mdl + (self._del_amp,)*p_amp)

    def set_deltas(self, mdl_delta=None, amp_delta=None):
        if mdl_delta is not None: self._del_mdl = mdl_delta
        if amp_delta is not None: self._del_amp = amp_delta

    def summary(self):
        s = 'Simplex optimisation deltas'
        s += '\nTLS Model Step Size:     {}'.format(self._del_mdl)
        s += '\nTLS Amplitude Step Size: {}'.format(self._del_amp)
        return s

class UijSimplex(_Simplex):
    _del_uij = 0.001

    def __init__(self, uij_delta=None):
        """Initialise AtomicUij Simplex"""
        self.set_deltas(uij_delta=uij_delta)

    def get_deltas(self):
        return numpy.array((self._del_uij,)*6)

    def set_deltas(self, uij_delta=None):
        if uij_delta is not None: self._del_uij = uij_delta

    def summary(self):
        s = 'Simplex optimisation deltas'
        s += '\nUij Residual Step Size:  {}'.format(self._del_uij)
        return s

############################################################################

class MultiDatasetUijParameterisation(Program):

    master_phil = master_phil

    def __init__(self,
                 models,
                 params,
                 levels,
                 level_labels=None,
                 log=None):
        """Object for fitting a series of TLS models to a set of structures"""

#        if log is None: log = Log(verbose=False)
        self.log = log

        # List of non-fatal errors from the program (to be reported at the end)
        self.errors = []

        self.params = params

        self._n_cpu = params.settings.cpus
        self._n_opt = params.fitting.optimisation.max_datasets

        self._allow_isotropic = True

        self._opt_datasets_res_limit = params.fitting.optimisation.max_resolution
        self._opt_datasets_selection = []

        self.models = models
        self.levels = levels
        self.level_labels = level_labels if level_labels else range(1,len(levels)+1)
        self.fitter = None

        # Misc files
        self.cifs = None

        # Create plot object
        self.plot = MultiDatasetUijPlots

        # Validate and add output paths, etc.
        self._init_file_system()
        self._init_input_models()
        self._init_level_groups()
        self._init_tables()
        self._init_fitter()

        #self.write_running_parameters_to_log(params=params)

    def _init_file_system(self):

        self.out_dir = easy_directory(self.params.output.out_dir)
        fm = self.initialise_file_manager(rootdir=self.out_dir)
        for d in ['logs','model','graphs','structures','table_ones']:
            fm.add_dir(dir_name=d, dir_tag=d)

        if self.params.output.diagnostics is True:
            for d in ['diagnostics']:
                fm.add_dir(dir_name=d, dir_tag=d)

    def _init_input_models(self):
        """Prepare the input models"""

        self.log.subheading('Processing input models')

        # Use the first hierarchy as the reference
        self.log('Using {} as reference structure'.format(self.models[0].tag))
        self.master_h = self.models[0].hierarchy.deep_copy()

        errors = []

        # Extract the TableOne columns as a set to compare with the columns in the MTZs
        table_one_cols = set(self.params.table_ones_options.column_labels.split(',')+[self.params.table_ones_options.r_free_label])

        for i_m, m in enumerate(self.models):
            # Check that all of the structures are the same
            if not self.master_h.is_similar_hierarchy(m.hierarchy):
                errors.append(Failure("Structures are not all the same. Model {}. File: {}".format(i_m, m.filename)))
                continue

            # Directory for this model
            m_dir = easy_directory(os.path.join(self.file_manager.get_dir('structures'), m.tag))

            # Copy input files to folder
            m.i_pdb = os.path.join(m_dir, m.tag+'.input.pdb')
            m.i_mtz = os.path.join(m_dir, m.tag+'.input.mtz')
            m.o_pdb = None
            m.o_mtz = None
            m.r_pdb = None
            m.r_mtz = None

            # Write the input pdb to the output folder (without the header as this can affect R-factors calculations)
            m.hierarchy.write_pdb_file(m.i_pdb, crystal_symmetry=m.crystal_symmetry)
            # Link the input mtz to the output folder
            rel_symlink(m.filename.replace('.pdb','.mtz'), m.i_mtz)

            assert os.path.exists(m.i_pdb), 'PDB does not exist: {}'.format(m.i_pdb)
            assert os.path.exists(m.i_mtz), 'MTZ does not exist: {}'.format(m.i_mtz)

            cs = CrystalSummary.from_mtz(m.i_mtz)
            if (self._opt_datasets_res_limit is None) or (cs.high_res < self._opt_datasets_res_limit):
                self._opt_datasets_selection.append(i_m)

            if table_one_cols.difference(cs.column_labels):
                errors.append(Failure("MTZ {} does not contain the correct columns.".format(m.filename) + \
                                      "\n\tLooking for: {}".format(','.join(table_one_cols)) + \
                                      "\n\tMTZ contains: {}".format(','.join(cs.column_labels)) + \
                                      "\n\tCan't find: {}".format(','.join(table_one_cols.difference(cs.column_labels))) + \
                                      "\nChange required columns with {}".format("table_ones_options.column_labels or table_ones_options.r_free_label")))

        # Check for errors - Structures not all the same?
        if errors:
            self.log.heading('Errors found while processing input files')
            for e in errors:
                self.log(str(e))
                self.log('')
            self.log.bar(False, True)
            raise Failure("Some structures are invalid or missing information.")
        # Check for errors - No high resolution structures?
        if len(self._opt_datasets_selection) == 0:
            raise Sorry('No datasets above resolution cutoff: {}'.format(self._opt_datasets_res_limit))

        # Limit the number of datasets for optimisation
        if self._n_opt is not None:
            self.log('Limiting list of datasets for TLS optimisation to {} datasets'.format(self._n_opt))
            self._opt_datasets_selection = self._opt_datasets_selection[:self._n_opt]
        # Print optimisation datasets
        self.log('Using {} datasets for TLS and residual parameterisation'.format(len(self._opt_datasets_selection)))
        for i_m in self._opt_datasets_selection:
            self.log('\t'+self.models[i_m].tag)

    def _init_level_groups(self):
        """Create the selection array that builds up the hierarchy of the fitting"""

        self.log.subheading('Processing input levels')

        # Extract the atoms for each tls group
        atom_cache = self.master_h.atom_selection_cache()
        # Array of which atoms are in which group at which level
        level_array = -1 * numpy.ones((len(self.levels), self.master_h.atoms().size()), dtype=int)
        # List of any selections that result in no atoms
        failures = []
        # Each column is one level
        for i_level, selections in enumerate(self.levels):
            self.log('\n> Level {} ({})\n'.format(i_level+1, self.level_labels[i_level]))
            for i_group, sel_group in enumerate(selections):
                sel = numpy.array(atom_cache.selection(sel_group))
                self.log('\tgroup {:<5d} - {:50}: {:>5d} atoms'.format(i_group+1, sel_group, sum(sel)))
                if sum(sel) == 0:
                    failures.append('Group "{}": {} atoms'.format(sel_group, sum(sel)))
                level_array[i_level, sel] = i_group+1
        if failures: raise Failure('One or more group selections do not select any atoms: \n\t{}'.format('\n\t'.join(failures)))

        # Find all atoms that are selected in at least one level of the mask
        self.atom_mask = (level_array!=-1).any(axis=0)
        # Filter and store array
        self.level_array = level_array[:,self.atom_mask]

        self.log.subheading('Level atom counts')
        self.log('> {} levels for {} atoms'.format(*self.level_array.shape))
        for i_level, label in enumerate(self.level_labels):
            self.log('\tLevel {} ({}): {} atoms'.format(i_level+1, label, sum(self.level_array[i_level]!=-1)))

    def _init_tables(self):
        """Create tables"""

        self.tables = Meta(['statistics', 'tracking'])
        self.tables.tracking = pandas.DataFrame(columns=['cycle', 'level', #'group',
                                                         'rmsd',
                                                         'u_iso (level)',
                                                         'b_iso (level)',
                                                         'u_iso (total)',
                                                         'b_iso (total)'])

    def _init_fitter(self):
        """Create fitting object"""

        # Extract all atoms from all datasets
        atoms = [m.hierarchy.atoms() for m in self.models]
        assert len(atoms) > 0, 'No models have been used?!'
        assert len(atoms[0]) > 0, 'No atoms have been extracted from models'

        # Extract uij and xyz from all datasets (and only select atoms we're interested in)
        observed_uij = numpy.array([a.extract_uij() for a in atoms])[:,self.atom_mask]
        observed_xyz = numpy.array([a.extract_xyz() for a in atoms])[:,self.atom_mask]

        # Check all uijs are present
        if (observed_uij==-1).all() and (self._allow_isotropic is True):
            # Set model type to isotropic
            self.disorder_model = 'isotropic'
            # Extract isotropic component from atoms (again only for the atoms we're interested in)
            self.log('All atoms are missing anisotropic displacement parameters -- using the isotropic parameters instead')
            observed_uij = numpy.zeros_like(observed_uij)
            observed_uij[:,:,0] = numpy.array([a.extract_b()/EIGHT_PI_SQ for a in atoms])[:,self.atom_mask]
            observed_uij[:,:,2] = observed_uij[:,:,1] = observed_uij[:,:,0]
        elif (observed_uij==-1).any():
            raise Failure('Some atoms for fitting (but not all) do not have anisotropically-refined B-factors -- either all atoms for fitting must have anistropic atoms or none')
        else:
            # Set model type to anisotropic
            self.disorder_model = 'anisotropic'

        # Create the fitting object
        self.fitter = MultiDatasetHierarchicalUijFitter(observed_uij = observed_uij,
                                                        observed_xyz = observed_xyz,
                                                        level_array  = self.level_array,
                                                        level_labels = self.level_labels,
                                                        params = self.params.fitting,
                                                        log = self.log)
        # Set whether input is anisotropic or isotropic
        self.fitter.set_input_info(disorder_model=self.disorder_model)
        # Select the datasets for be used for TLS optimisation
        self.fitter.set_optimisation_datasets(self._opt_datasets_selection)

        # Table to be populated during optimisation
        self.fitter.set_tracking(table=self.tables.tracking, out_dir=self.file_manager.get_dir('graphs'))

        # Write summary of the fitted model (groups & levels)
        model_dir = self.file_manager.get_dir('model')
        #self.log.heading('Writing summary of the hierarchical model')
        #self.hierarchy_summary(out_dir=model_dir)

        # Calculate Parameter-observed data ratio
        self.log.subheading('Data-Parameter Ratios')
        n_o = numpy.product(self.fitter.observed_uij.shape)
        n_p = self.fitter.n_params(non_zero=False)
        self.log('Number of model parameters: {}'.format(n_p))
        self.log('Number of input uij values: {}'.format(n_o))
        self.log('Parameter ratio Gain (should be <1): {}'.format(float(n_p)/float(n_o)))

    def show_errors(self):
        """Print non-fatal errors accumulated during running"""
        if len(self.errors) == 0:
            self.log('Program completed with 0 errors')
            return
        self.log.subheading('{} non-fatal errors occurred during the program'.format(len(self.errors)))
        for i, e in enumerate(self.errors):
            self.log.bar()
            self.log('Error {} of {}'.format(i+1, len(self.errors)))
            self.log.bar()
            self.log(e)
        self.log.bar()
        self.log.subheading('{} non-fatal errors occurred during the program'.format(len(self.errors)))

    def blank_master_hierarchy(self):
        h = self.master_h.deep_copy()
        h.atoms().set_uij(flex.sym_mat3_double(h.atoms().size(), [0.0]*6))
        h.atoms().set_b(flex.double(h.atoms().size(), 0.0))
        return h

    def custom_master_hierarchy(self, uij=None, iso=None, mask=None):
        """Copy of the master hierarchy with custom ADPs"""
        m_h = self.blank_master_hierarchy()
        m_a = m_h.atoms()
        if mask is not None:
            m_a = m_a.select(mask)
        # Apply new uijs
        if uij is not None:
            m_a.set_uij(flex.sym_mat3_double(uij))
        if iso is not None:
            m_a.set_b(flex.double(list(iso)))
        return m_h

    def partition_boundaries_for_level(self, i_level):
        """Find the boundaries between the level partitions"""
        l = self.level_array[i_level]
        mask = self.atom_mask.tolist()
        return self.partition_boundaries_custom(atom_labels=l, mask=mask)

    def partition_boundaries_custom(self, atom_labels, mask=None):
        """Find the boundaries for labelled regions"""
        atom_labels = numpy.array(atom_labels)
        mask=flex.bool(mask)
        b = self.blank_master_hierarchy()
        for i in numpy.unique(atom_labels):
            # Set the b-factors to the membership of the group
            h = self.blank_master_hierarchy()
            h.select(mask).atoms().set_b(flex.double((atom_labels==i).tolist()))
            # Where are the boundaries of the group?
            sel_int = numpy.array(h.atoms().extract_b(), dtype=int)
            boundaries = numpy.array(b.atoms()[:-1].extract_b(), dtype=bool) + numpy.array((sel_int[:-1]*(1-sel_int[1:]))+(sel_int[1:]*(1-sel_int[:-1])), dtype=bool)
            b.atoms()[:-1].set_b(flex.double(boundaries.tolist()))
        return b

    def fit_hierarchical_uij_model(self):
        """Optimise the TLS+residual model against the input data"""

        self.log.heading('Fitting hierarchical B-factor model', spacer=True)

        n_macro_cycles = self.params.fitting.number_of_macro_cycles
        n_micro_cycles = self.params.fitting.number_of_micro_cycles
        self.log('Macro-cycles: {}'.format(n_macro_cycles))
        self.log('Micro-cycles: {}'.format(n_micro_cycles))
        self.log('')

        # Run the fitting
        self.fitter.fit(n_cpus = self._n_cpu,
                        n_macro_cycles = n_macro_cycles,
                        n_micro_cycles = n_micro_cycles)

        self.log.heading('Parameterisation complete')

    def process_results(self):
        """Extract output and generate summaries"""

        self.log.heading('Processing results', spacer=True)

        #------------------------------------------------------------------------------#
        #---#                            Extract uijs                              #---#
        #------------------------------------------------------------------------------#

        self.log('Extracting output Uijs from parameterised model...')

        # Extract the different components of the fitted uij values
        uij_lvl = self.fitter.extract_tls(sum_levels=False)
        uij_tls = uij_lvl.sum(axis=0)
        uij_res = self.fitter.residual.extract()
        uij_all = uij_tls + uij_res
        # Extract the input uij values
        uij_inp = self.fitter.observed_uij

        #------------------------------------------------------------------------------#
        #---#                     Write output structures                          #---#
        #------------------------------------------------------------------------------#

        self.log.subheading('Writing various sets of output structures')

        # Output structure (will contain folders for each dataset)
        structure_dir = self.file_manager.get_dir('structures')

        # Create TLS headers for each dataset
        self.log('Creating TLS headers for each structure')
        tls_headers = self.make_tls_headers(datasets=None, levels=None)
        # Fully parameterised structures
        self.log.bar()
        self.log('Writing fitted structures (full uij models)')
        self.log.bar()
        pdbs = self.output_structures(uij = uij_all,
                                      iso = map(uij_to_b, uij_all),
                                      headers = tls_headers,
                                      out_dir = structure_dir,
                                      model_suffix = '.all.pdb')
        # Add these for models for table ones
        for i, mdl in enumerate(self.models):
            mdl.o_pdb = pdbs[i]
            mdl.o_mtz = pdbs[i].replace('.pdb','.mtz')
            if not os.path.exists(mdl.o_mtz):
                rel_symlink(mdl.i_mtz, mdl.o_mtz)
        # TLS-parameterised structures
        self.log.bar()
        self.log('Writing fitted structures (TLS components only)')
        self.log.bar()
        self.output_structures(uij = uij_tls,
                               iso = map(uij_to_b, uij_tls),
                               headers = tls_headers,
                               out_dir = structure_dir,
                               model_suffix = '.tls.pdb')
        # Level by level TLS-parameterised structures
        for i_level in xrange(len(self.levels)):
            self.log.bar()
            self.log('Writing fitted structures (TLS components only - level {})'.format(i_level+1))
            self.log.bar()
            self.output_structures(uij = uij_lvl[i_level],
                                   iso = map(uij_to_b, uij_lvl[i_level]),
                                   headers = self.make_tls_headers(datasets=None, levels=[i_level]),
                                   out_dir = structure_dir,
                                   model_suffix = '.tls-level-{:04}.pdb'.format(i_level))
        # TLS-subtracted structures ("residual structures")
        self.log.bar()
        self.log('Writing fitted structures (TLS subtracted from input Uijs)')
        self.log.bar()
        self.output_structures(uij = uij_inp-uij_tls,
                               iso = map(uij_to_b, uij_inp-uij_tls),
                               out_dir = structure_dir,
                               model_suffix = '.tls_subtracted.pdb')

        #------------------------------------------------------------------------------#
        #---#                         Model summaries                              #---#
        #------------------------------------------------------------------------------#

        model_dir = self.file_manager.get_dir('model')

        # Write average uijs for each level
        self.log.heading('Summaries of the level-by-level TLS fittings')
        self.tls_level_summary(out_dir=model_dir)

        # Write residuals
        self.log.heading('Summaries of the residual atomic components')
        self.residuals_summary(out_dir=model_dir)

        # Distributions of the uijs for groups
        self.log.heading('Calculating distributions of uijs over the model')
        self.uij_distribution_summary(uij_lvl=uij_lvl,
                                      uij_res=uij_res,
                                      #uij_fit=uij_all,
                                      uij_inp=uij_inp,
                                      out_dir=model_dir)

        #------------------------------------------------------------------------------#
        #---#        Compare the fittend and input uijs for all structures         #---#
        #------------------------------------------------------------------------------#

        if self.params.output.diagnostics is True:

            self.file_manager.add_dir(dir_name='diagnostics', dir_tag='diagnostics')
            fit_dir = self.file_manager.get_dir('diagnostics')

            # Dataset-by-dataset and atom-by-atom fit rmsds
            self.log.heading('Calculating rmsd metrics between input and fitted Uijs')
            self.fit_rms_distributions(uij_fit=uij_all,
                                       uij_inp=uij_inp,
                                       out_dir=fit_dir)

            # Correlations between the residual atomic component and the
            self.log.heading('Calculating correlations between atomic residuals and input-fitted Uijs')
            self.fit_residual_correlations(uij_diff=(uij_inp-uij_all),
                                           out_dir=fit_dir)

        #------------------------------------------------------------------------------#
        #---#              Refine the output models if requested                   #---#
        #------------------------------------------------------------------------------#

        if self.params.refine.refine_output_structures:
            self.log.heading('Refining output structures')
            self.refine_fitted_dataset_models()

        #------------------------------------------------------------------------------#
        #---#             Analyse the quality of the output models                 #---#
        #------------------------------------------------------------------------------#

        # Calculate new R-frees, etc.
        self.log.heading('Generating Table Ones for all structures')
        self.generate_fitted_table_ones()

        # Write output CSV of... everything
        self.log.heading('Writing output csvs')
        self.write_combined_csv(uij_fit=uij_all,
                                uij_inp=uij_inp,
                                out_dir=self.file_manager.get_dir('root'))

        #------------------------------------------------------------------------------#
        #                                                                              #
        #------------------------------------------------------------------------------#

    def output_structures(self,  uij, iso=None, headers=None, out_dir='./', model_suffix='.pdb'):
        """Write sets of uij to models."""

        # Make sure ouput directory exists
        easy_directory(out_dir)

        # Validate AnoUijs
        uij = numpy.array(uij)
        assert uij.shape == (len(self.models), sum(self.atom_mask), 6)
        # Validate IsoBs
        if iso is not None:
            iso = numpy.array(iso)
            assert iso.shape == (len(self.models), sum(self.atom_mask))
        # Validate header lines
        if headers is not None:
            if isinstance(headers[0], str):
                headers = [headers]
            assert set(map(len,headers)) == {len(self.models)}
            open_append = True
        else:
            open_append = False

        # Mask to allow us to apply uijs back to structure
        sel = flex.bool(self.atom_mask.tolist())
        # List of output pdbs
        pdbs = []

        # Apply to each model and write out
        for i, mdl in enumerate(self.models):
            # Create model paths
            mdl_d = easy_directory(os.path.join(out_dir, mdl.tag))
            mdl_f = os.path.join(mdl_d, mdl.tag+model_suffix)
            self.log('{} > {}'.format(mdl.tag, mdl_f))
            try:
                # Create new copy of hierarchy and extract atoms
                h = mdl.hierarchy.deep_copy()
                a = h.atoms().select(sel)
                # Apply new uijs
                a.set_uij(flex.sym_mat3_double(uij[i]))
                # Apply iso if given or set to 0
                if iso is not None:
                    a.set_b(flex.double(iso[i]))
                else:
                    a.set_b(flex.double(a.size(), 0))
                # Write headers
                if headers is not None:
                    with open(mdl_f, 'w') as fh:
                        [fh.write(l[i].strip()+'\n') for l in headers]
                h.write_pdb_file(mdl_f, open_append=open_append, crystal_symmetry=mdl.crystal_symmetry)
                pdbs.append(mdl_f)
            except:
                self.log.bar()
                self.log('ERROR: Failed to write structure for dataset {}: ({})'.format(mdl.tag, mdl_f))
                self.log(traceback.format_exc())
                self.log.bar()
                continue

        return pdbs

    def hierarchy_summary(self, out_dir='./'):
        """Write out the composition of the hierarchical model"""

        out_dir = easy_directory(out_dir)

        # Extract the global mask and convert to flex
        global_sel = flex.bool(self.atom_mask.tolist())

        # Write out the groups for each level
        for i_level, level in enumerate(self.fitter.levels):
            self.log('Writing partition groups for level {}'.format(i_level+1))
            # Add all structures to the same root
            level_root = None
            for i_group, sel, fitter in level:
                sel = flex.bool(sel.tolist())
                g = self.blank_master_hierarchy()
                g.atoms().select(global_sel).select(sel).set_b(flex.double(sum(sel), 1))
                if level_root is None: level_root = g
                else: [level_root.append_model(m.detached_copy()) for m in g.models()]
            filename = os.path.join(out_dir, 'level-{:04d}_atoms.pdb'.format(i_level+1))
            self.log('\t> {} ({} models)'.format(filename, len(level_root.models())))
            level_root.write_pdb_file(filename)

        # Generate hierarchy for each level with groups as b-factors
        self.log('Writing partition figures for each chain')
        level_hierarchies = []
        for g_vals in self.level_array:
            h = self.custom_master_hierarchy()
            h = self.custom_master_hierarchy(uij=None, iso=g_vals, mask=global_sel)
            level_hierarchies.append(h)
        # Write hierarchy plot for each chain
        b_h = self.blank_master_hierarchy()
        b_c = b_h.atom_selection_cache()
        for c in b_h.chains():
            chain_sel = b_c.selection('chain {}'.format(c.id))
            hierarchies = [h.select(chain_sel, copy_atoms=True) for h in level_hierarchies]
            # Skip if no partitions in this chain
            if not (numpy.array([h.atoms().extract_b() for h in hierarchies]).sum() > 1e-3):
                continue
            filename = os.path.join(out_dir, 'level-partitioning-chain-{}.png'.format(c.id))
            self.log('\t> {}'.format(filename))
            self.plot.level_plots(filename=filename, hierarchies=hierarchies, title='chain {}'.format(c.id))

    def make_tls_headers(self, levels=None, datasets=None):
        """Create header lines for each pdb file containing the TLS matrices"""

        # If no levels given take all levels
        if levels is None: levels = range(len(self.levels))
        # If no datasets given take all datasets
        n_dst = len(datasets) if (datasets is not None) else len(self.models)

        # Validate levels and ensure they're sorted!
        assert isinstance(levels, list)
        assert set(map(type, levels)) == {int}
        levels = sorted(levels)
        n_levels = len(levels)

        # Extract the centres of mass for each dataset
        coms = self.fitter.observed_com if (datasets is None) else self.fitter.observed_com[datasets]

        # Extract the selected level objects
        sel_level_objs = [self.fitter.levels[i] for i in levels]
        # Extract the selection strings for all levels
        sel_level_selections = [self.levels[i] for i in levels]
        # Extract the group assignments for the selected levels
        sel_level_array = self.level_array[levels]
        # Extract the unique group assignments for the selected levels
        unique_combinations = sorted(set(map(tuple,sel_level_array.T.tolist())))

        # List of tlso objects for each dataset
        dataset_objs = [[] for i in xrange(n_dst)]

        # Iterate through the different UNIQUE grouping combinations
        for sel_groups in unique_combinations:

            # Groups are indexed from 1 but level list is indexed from 0
            sel_indexs = [g-1 if g>0 else None for g in sel_groups]
            # All Nones -> no contribution for this part for this level
            if sel_indexs.count(None) == len(sel_indexs):
                continue

            # Combine the selection strings for this combination
            selection_string = '('+') and ('.join([sel_level_selections[i_l][i_g] for i_l, i_g in enumerate(sel_indexs) if (i_g is not None)])+')'

            # Extract the fitters for each level
            fitters = [sel_level_objs[i_l].fitters.get(g, None) for i_l, g in enumerate(sel_groups)]
            # Extract the expanded TLS models for each dataset for each level
            tls_models = zip(*[f.parameters().expand(datasets=datasets) for f in fitters if (f is not None)])
            assert len(tls_models) == n_dst
            assert len(tls_models[0]) == len(fitters)-fitters.count(None)

            # Iterate through the datasets and create TLSOs for each dataset
            for i_dst, tls_list in enumerate(tls_models):
                # Add all the levels for each dataset
                total_model = reduce(operator.add, tls_list)
                # Extract TLS remark lines from the
                tls_obj = mmtbx.tls.tools.tlso(t=total_model.get('T'), l=total_model.get('L'), s=total_model.get('S'), origin=coms[i_dst])
                # Append to the tlso list
                dataset_objs[i_dst].append((tls_obj, selection_string))

        # Create header for each dataset
        dataset_headers = []
        for tls_tuples in dataset_objs:
            # Turn list of tlso-string tuples into separate lists
            tlsos, sel_strs = zip(*tls_tuples)
            # Store output in list stream
            l = ListStream()
            # Get the remark records
            mmtbx.tls.tools.remark_3_tls(tlsos=tlsos, selection_strings=sel_strs, out=l)
            # And finally append
            dataset_headers.append(str(l))
        assert len(dataset_headers) == n_dst

        return dataset_headers

    def tls_level_summary(self, out_dir='./'):
        """Write the various TLS uijs to the master hierarchy structure"""

        out_dir = easy_directory(out_dir)
        csv_dir = easy_directory(os.path.join(out_dir, 'csvs'))
        png_dir = easy_directory(os.path.join(out_dir, 'graphs'))
        pml_dir = easy_directory(os.path.join(out_dir, 'pymol'))

        self.log.subheading('Writing TLS models and amplitudes for each level')
        # Iterate through the levels
        for level in self.fitter.levels:
            self.log('Level {}'.format(level.index))
            # Table for TLS model components
            mdl_filename = os.path.join(csv_dir, 'tls_models_level_{:04d}.csv'.format(level.index))
            mdl_table = pandas.DataFrame(columns=["group", "mode",
                                                  "T11","T22","T33","T12","T13","T23",
                                                  "L11","L22","L33","L12","L13","L23",
                                                  "S11","S12","S13","S21","S22","S23","S31","S32","S33"])
            # Create amplitude table
            amp_filename = os.path.join(csv_dir, 'tls_amplitudes_level_{:04d}.csv'.format(level.index))
            amp_table = pandas.DataFrame(columns=["group", "mode", "cpt"]+[mdl.tag for mdl in self.models])
            # Create list of plot arguments and generate all at end
            #plot_args = []
            # Iterate through the groups in this level
            for i_group, sel, fitter in level:
                tls_model, tls_amps = fitter.result()
                assert tls_model.shape == (self.params.fitting.number_of_modes_per_group, 21)
                assert tls_amps.shape  == (self.params.fitting.number_of_modes_per_group, len(self.models), 3)

                # Add to model and amplitudes tables
                for i_tls in xrange(tls_model.shape[0]):
                    # Add model values to last row of table
                    mdl_table.loc[len(mdl_table.index)] = numpy.concatenate([[i_group, i_tls], tls_model[i_tls]])
                    # Add amplitudes to last row of table
                    for i_cpt, cpt in enumerate('TLS'):
                        amp_table.loc[len(amp_table.index)] = numpy.concatenate([[i_group, i_tls, cpt], tls_amps[i_tls,:,i_cpt]])

                # Write histograms of amplitudes -- only for non-zero models
                if (tls_model.sum() > 0.0) and self.params.output.group_graphs:
                    x_vals = []; [[x_vals.append(tls_amps[i_m,:,i_c]) for i_c in xrange(tls_amps.shape[2])] for i_m in xrange(tls_amps.shape[0])]
                    filename = os.path.join(png_dir, 'tls-model-amplitudes-level-{}-group-{}.png'.format(level.index, i_group))
                    self.log('\t> {}'.format(filename))
                    self.plot.multi_histogram(filename  = filename,
                                              x_vals    = x_vals,
                                              titles    = numpy.concatenate(['T (mode {a})-L (mode {a})-S (mode {a})'.format(a=i_t+1).split('-') for i_t in xrange(tls_amps.shape[0])]),
                                              x_labs    = ['']*tls_amps.shape[0]*tls_amps.shape[2],
                                              rotate_x_labels = True,
                                              shape     = (tls_amps.shape[0], tls_amps.shape[2]),
                                              n_bins    = 30)
            # Write tables
            mdl_table.to_csv(mdl_filename)
            amp_table.to_csv(amp_filename)
            # Generate plots
            #self.log('Generating {} histogram plots of TLS mode amplitudes'.format(len(plot_args)))
            #libtbx.easy_mp.pool_map(processes=self._n_cpu, func=wrapper_plot_histograms, args=plot_args)

        self.log.subheading('Writing T-L-S Uij components for each level')
        # Extract the atom mask to apply b-factors
        sel = flex.bool(self.atom_mask.tolist())
        # Iterate through the levels and plot T-L-S contributions for each mode for each level
        for i_level, level in enumerate(self.fitter.levels):
            self.log('Level {}'.format(level.index))
            # Cumulative uijs (for this level)
            str_lvl = {}
            uij_lvl = {'T':0,'L':0,'S':0}
            # Boundaries for this level
            boundaries = self.partition_boundaries_for_level(i_level=i_level)
            # Iterate through the different tls models
            for i_tls in xrange(self.params.fitting.number_of_modes_per_group):
                # Prefix for structures & graphs of uijs for this level and mode
                prefix = os.path.join(out_dir, 'level_{}-mode_{}'.format(i_level+1, i_tls+1))
                # Cumulative uijs (for this mode)
                str_mode = {}
                uij_mode = {'T':0,'L':0,'S':0}
                # Extract the T-L-S for this mode
                for cpt in 'TLS':
                    # Boolean selections
                    t = 'T' in cpt
                    l = 'L' in cpt
                    s = 'S' in cpt
                    self.log('\tExtracting {}{}{} components of mode {} of level {}'.format('T' if t else '_','L' if l else '_','S' if s else '_', i_tls+1, i_level+1))
                    # Create copy of the level for resetting T-L-S components
                    l_copy = copy.deepcopy(level)
                    # Reset the required TLS of this mode
                    l_copy.reset_models(modes=[i_tls], t=(not t), l=(not l), s=(not s))
                    # Reset the TLS for all other modes
                    if self.params.fitting.number_of_modes_per_group > 1:
                        other_modes = range(self.params.fitting.number_of_modes_per_group)
                        other_modes.remove(i_tls)
                        l_copy.reset_models(modes=other_modes, t=True, l=True, s=True)
                    # Extract uijs
                    uij = l_copy.extract(average=True)

                    # Create structure for this level and mode
                    m_h = self.custom_master_hierarchy(uij=uij, iso=uij_to_b(uij), mask=sel)
                    # Only write if more than one mode
                    if self.params.fitting.number_of_modes_per_group > 1:
                        m_f = prefix+'-{}.pdb'.format(cpt)
                        self.log('\t> {}'.format(m_f))
                        m_h.write_pdb_file(m_f)

                    # Store this structure for plotting
                    str_mode[cpt] = m_h
                    uij_mode[cpt] = uij
                    # Add to cumulative values
                    uij_lvl[cpt] += uij_mode[cpt]

                # Write output for this MODE
                if self.params.fitting.number_of_modes_per_group > 1:
                    # Write structure with combined TLS values
                    uij_mode_tls = uij_mode['T'] + uij_mode['L'] + uij_mode['S']
                    m_h = self.custom_master_hierarchy(uij=uij_mode_tls, iso=uij_to_b(uij_mode_tls), mask=sel)
                    m_f = prefix+'-TLS.pdb'
                    self.log('\t> {}'.format(m_f))
                    m_h.write_pdb_file(m_f)
                    # Write stacked bar plot of TLS for this mode
                    self.log('\t> {}xxx.png'.format(prefix+'-TLS'))
                    self.plot.stacked_bar(prefix=prefix+'-TLS',
                                          hierarchies=[str_mode[cpt] for cpt in 'TLS'],
                                          legends=list('TLS'),
                                          title='Level {}, Mode {} - TLS Contributions'.format(i_level+1, i_tls+1),
                                          v_line_hierarchy=boundaries)

            self.log('\tWriting structures')
            # Prefix for structures & graphs of uijs for this level and mode
            prefix = os.path.join(out_dir, 'level_{}-all-modes'.format(i_level+1))
            # Dictionary to hold the TLS components for this mode+level
            str_lvl = {}
            # Write out structure & graph of uijs (this level)
            for cpt in 'TLS':
                uij = uij_lvl[cpt]
                # Create structure for this level and mode
                m_h = self.custom_master_hierarchy(uij=uij, iso=uij_to_b(uij), mask=sel)
                m_f = prefix+'-{}.pdb'.format(cpt)
                self.log('\t> {}'.format(m_f))
                m_h.write_pdb_file(m_f)
                # Store this structure for plotting
                str_lvl[cpt] = m_h
            # Write structure with combined TLS values
            uij_lvl_tls = uij_lvl['T'] + uij_lvl['L'] + uij_lvl['S']
            m_h = self.custom_master_hierarchy(uij=uij_lvl_tls, iso=uij_to_b(uij_lvl_tls), mask=sel)
            m_f = prefix+'-TLS.pdb'
            self.log('\t> {}'.format(m_f))
            m_h.write_pdb_file(m_f)
            # Write pymol images of each chain and residue
            if self.params.output.pymol_images:
                self.log('\tGenerating pymol images')
                pml_prefix = os.path.join(pml_dir, os.path.basename(prefix))
                self.log('\t> {}xxx.png'.format(pml_prefix))
                # Images for each chain
                auto_chain_images(structure_filename = m_f,
                                  output_prefix = pml_prefix,
                                  style = 'lines+ellipsoids',
                                  width=1000, height=750)
                # Images for each group
                blank_h = self.blank_master_hierarchy()
                cache_h = blank_h.atom_selection_cache()
                selection_images(structure_filename = m_f,
                                 output_prefix = pml_prefix+'-group_',
                                 selections = [PymolSelection.join_or([PymolSelection.format(a) for a in blank_h.atoms().select(cache_h.selection(s))]) for s in self.levels[i_level]],
                                 style = 'lines+ellipsoids',
                                 width = 250 if level.n_groups()>50 else 1000,
                                 height= 250 if level.n_groups()>50 else 750)
            # Write stacked bar plot of TLS for this level
            self.log('\tWriting stacked bar plot of TLS contributions for this level')
            self.log('\t> {}xxx.png'.format(prefix+'-TLS'))
            self.plot.stacked_bar(prefix=prefix+'-TLS',
                                  hierarchies=[str_lvl[cpt] for cpt in 'TLS'],
                                  legends=list('TLS'),
                                  title='Level {} - TLS Contributions'.format(i_level+1),
                                  v_line_hierarchy=boundaries if (sum(boundaries.atoms().extract_b()) < 0.5*len(list(boundaries.select(sel).residue_groups()))) else None)
            #embed()

        return

    def residuals_summary(self, out_dir='./'):
        """Write the residual uijs to the master hierarchy structure"""

        out_dir = easy_directory(out_dir)
        pml_dir = easy_directory(os.path.join(out_dir, 'pymol'))

        uij = self.fitter.residual.extract()

        # Write out structure & graph of uijs (this level)
        self.log('Writing structure and plot of residual ADPs')
        prefix = os.path.join(out_dir, 'all-residual')
        m_h = self.custom_master_hierarchy(uij=uij, iso=uij_to_b(uij), mask=flex.bool(self.atom_mask.tolist()))
        m_f = prefix+'.pdb'
        self.log('\t> {}'.format(m_f))
        m_h.write_pdb_file(m_f)
        self.log('\t> {}xxx.png'.format(prefix))
        self.plot.stacked_bar(prefix=prefix,
                              hierarchies=[m_h],
                              legends=['Residual'],
                              title='Residual Level - fixed contributions')
        # Write pymol images of each chains and residue
        if self.params.output.pymol_images:
            self.log('\tGenerating pymol images')
            pml_prefix = os.path.join(pml_dir, 'residual')
            self.log('\t> {}xxx.png'.format(pml_prefix))
            auto_chain_images(structure_filename = m_f,
                              output_prefix = pml_prefix,
                              style = 'lines+ellipsoids',
                              width=1000, height=750)
            auto_residue_images(structure_filename = m_f,
                                output_prefix = pml_prefix,
                                style = 'lines+ellipsoids',
                                width=250, height=250)

        # TODO Output CSV of all residual components

        return

    def uij_distribution_summary(self, uij_lvl, uij_res, uij_inp, out_dir='./'):
        """Write the distributions of each level TLS contributions over the structure"""

        out_dir = easy_directory(out_dir)

        self.log('Writing distribution of TLS Uijs over the model')

        uij_hierarchies = []
        for i_level in range(len(uij_lvl)):
            uij_lvl_av = numpy.mean(uij_lvl[i_level], axis=0)
            h_lvl =self.custom_master_hierarchy(iso=uij_to_b(uij_lvl_av), mask=flex.bool(self.atom_mask.tolist()))
            uij_hierarchies.append(h_lvl)
        prefix = os.path.join(out_dir, 'all-stacked')
        self.log('\t> {}xxx.png'.format(prefix))
        self.plot.stacked_bar(prefix=prefix,
                              hierarchies=uij_hierarchies,
                              legends=['Level {}'.format(i+1) for i in range(len(uij_lvl))],
                              title='TLS Contributions',
                              reverse_legend_order=True)

    def fit_rms_distributions(self, uij_fit, uij_inp, out_dir='./', max_x_width=25):
        """Analyse the dataset-by-dataset and residue-by-residue and atom-by-atom fit qualities"""

        out_dir = easy_directory(out_dir)

        #------------------------------------------------------------------------------#
        # Create series of output directories
        #------------------------------------------------------------------------------#
        dst_dir = easy_directory(os.path.join(out_dir, 'dataset_by_dataset'))
        grp_dir = easy_directory(os.path.join(out_dir, 'group_by_group'))
        res_dir = easy_directory(os.path.join(out_dir, 'residue_by_residue'))
        atm_dir = easy_directory(os.path.join(out_dir, 'atom_by_atom'))

        #------------------------------------------------------------------------------#
        # Extract series of labels for datasets, atoms ...
        #------------------------------------------------------------------------------#

        dst_labels = numpy.array([m.tag for m in self.models])

        # Labels for each group in each level
        grp_labels = [[g.label for i,s,g in level] for level in self.fitter.levels]

        # Extract residue and atom labels for fitted atoms from the master
        m_h = self.blank_master_hierarchy().select(flex.bool(self.atom_mask.tolist()), copy_atoms=True)
        atm_labels = numpy.array([ShortLabeller.format(a) for a in m_h.atoms()])
        res_labels = numpy.array([ShortLabeller.format(a.parent().parent()) for a in m_h.atoms()])

        # Extract selections for each residue
        _, unq_res_idxs = numpy.unique(res_labels, return_index=True)
        unq_res_labs = res_labels[sorted(unq_res_idxs)]
        all_res_sels = [res_labels==l for l in unq_res_labs]

        #------------------------------------------------------------------------------#
        # Measure differences between input and fitted uijs
        #------------------------------------------------------------------------------#

        uij_diff = uij_inp-uij_fit

        # Calculate the rms of the difference (n_dsts x n_atoms)
        rmsds_all = rms(uij_diff, axis=2)
        rmsds_min = numpy.min(rmsds_all)
        rmsds_max = numpy.max(rmsds_all)

        #------------------------------------------------------------------------------#
        # Dataset-by-dataset "error" boxplots
        #------------------------------------------------------------------------------#

        self.log.subheading('Calculating dataset-by-dataset rmsds to input Uijs')
        for i_level, level in enumerate(self.fitter.levels):
            self.log('Level {}'.format(level.index))
            for i_group, sel, fitter in level:
                # Extract the uij rmsds for this group
                rmsds_sel = rmsds_all[:,sel]
                # Create boxplots of the rmsds by dataset
                for i_d in xrange(0, len(self.models), max_x_width):
                    filename = os.path.join(dst_dir, 'rmsds_level-{}_group-{:04d}_datasets-{:04d}-{:04d}.png'.format(level.index, i_group, i_d+1, i_d+max_x_width))
                    self.log('\t> {}'.format(filename))
                    self.plot.violinplot(filename=filename,
                                         y_vals=rmsds_sel.T[:,i_d:i_d+max_x_width],
                                         x_labels=dst_labels[i_d:i_d+max_x_width].tolist(),
                                         title='rmsd of fitted and refined B-factors (by dataset)',
                                         x_lab='dataset', y_lab='rmsd', rotate_x_labels=True,
                                         x_lim=(0,max_x_width+1), y_lim=(rmsds_min, rmsds_max))

        #------------------------------------------------------------------------------#
        # Atom-by-atom and residue-by-residue "error" boxplots
        #------------------------------------------------------------------------------#

        self.log.subheading('Calculating atom-by-atom and residue-by-residue rmsds to input Uijs')
        # Use this for the atom-by-atom graphs with atom_label=res_labels
        residue_boundaries = self.partition_boundaries_custom(atom_labels=res_labels, mask=self.atom_mask.tolist())
        residue_boundaries = residue_boundaries.select(flex.bool(self.atom_mask.tolist()))
        # Create boxplots of the rmsds by atom (for the first level only!)
        for i_group, sel, fitter in self.fitter.levels[0]:
            self.log.bar()
            self.log('Group {}'.format(i_group))
            self.log.bar()
            # Extract the uij rmsds for this group
            rmsds_sel = rmsds_all[:,sel]
            # Extract the residue selections for this group
            res_labels_sel = res_labels[sel]
            _, u_res_idxs = numpy.unique(res_labels_sel, return_index=True)
            g_res_labs = res_labels_sel[sorted(u_res_idxs)]
            g_res_sels = [res_labels_sel==l for l in g_res_labs]
            # Create boxplots of the rmsds by atom
            self.log('Atom-by-atom plots')
            for i_a in xrange(0, sum(sel), max_x_width):
                filename = os.path.join(atm_dir, 'rmsds_level-0_group_{:02d}_atoms-{:06d}-{:06d}.png'.format(i_group, i_a+1, i_a+max_x_width))
                self.log('\t> {}'.format(filename))
                self.plot.violinplot(filename=filename,
                                     y_vals=rmsds_sel[:,i_a:i_a+max_x_width],
                                     x_labels=atm_labels[sel].tolist()[i_a:i_a+max_x_width],
                                     title='rmsd of fitted and refined B-factors (by atom)',
                                     x_lab='atom', y_lab='rmsd', rotate_x_labels=True,
                                     x_lim=(0,max_x_width+1), y_lim=(rmsds_min, rmsds_max),
                                     vlines=(numpy.where(numpy.array(list(residue_boundaries.select(flex.bool(sel.tolist())).atoms().extract_b()[i_a:i_a+max_x_width]), dtype=bool))[0] + 1.5))
            self.log('Residue-by-residue plots')
            # Create boxplots of the rmsds by residue
            for i_r in xrange(0, len(g_res_sels), max_x_width):
                filename = os.path.join(res_dir, 'rmsds_level-0_group-{:02d}_residues-{:06d}-{:06d}.png'.format(i_group, i_r+1, i_r+max_x_width))
                self.log('\t> {}'.format(filename))
                self.plot.violinplot(filename=filename,
                                     y_vals=[rmsds_sel[:,a_sel].flatten() for a_sel in g_res_sels[i_r:i_r+max_x_width]],
                                     x_labels=g_res_labs[i_r:i_r+max_x_width],
                                     title='rmsd of fitted and refined B-factors (by residue)',
                                     x_lab='residue', y_lab='rmsd', rotate_x_labels=True,
                                     x_lim=(0,max_x_width+1), y_lim=(rmsds_min, rmsds_max))

        #------------------------------------------------------------------------------#
        # Atom-by-atom and residue-by-residue "error" boxplots
        #------------------------------------------------------------------------------#

        self.log('Writing average rmsd for each atom (average quality of fit over all datasets)')
        # calculate dataset-averaged rmsds
        rmsd_avg = numpy.mean(rmsds_all, axis=0)
        prefix = os.path.join(out_dir, 'rmsds_average')
        m_h = self.custom_master_hierarchy(uij=None, iso=rmsd_avg.tolist(), mask=flex.bool(self.atom_mask.tolist()))
        m_f = prefix+'.pdb'
        self.log('\t> {}'.format(m_f))
        m_h.write_pdb_file(m_f)
        self.log('\t> {}...png'.format(prefix))
        # Output as b-factor plot
        self.plot.residue_by_residue(prefix=prefix,
                                     hierarchy=m_h)

        self.log('Writing IQR rmsd for each atom (variability of quality of fit over all datasets)')
        # calculate dataset-IQR'd rmsds
        rmsd_iqr = numpy.subtract(*numpy.percentile(rmsds_all, [75, 25],axis=0))
        # Output as structure
        prefix = os.path.join(out_dir, 'rmsds_average')
        m_h = self.custom_master_hierarchy(uij=None, iso=rmsd_iqr.tolist(), mask=flex.bool(self.atom_mask.tolist()))
        m_f = prefix+'.pdb'
        self.log('\t> {}'.format(m_f))
        m_h.write_pdb_file(m_f)
        self.log('\t> {}...png'.format(prefix))
        # Output as b-factor plot
        self.plot.residue_by_residue(prefix=prefix,
                                     hierarchy=m_h)

    def fit_residual_correlations(self, uij_diff, out_dir='./', max_x_width=25):
        """Calculate correlations between the atomic residuals and the fitted-input differences"""

        cor_dir = easy_directory(os.path.join(out_dir, 'error_correlations'))

        #------------------------------------------------------------------------------#
        # Calculate correlations between the residual uij and the fitting error
        #------------------------------------------------------------------------------#
        self.log.subheading('Calculating correlations between residual uij and fitting errors')

        # Extract residual Uijs
        res = self.fitter.residual.extract()
        assert res.shape == uij_diff.shape[1:]
        # Calculate correlations -- iterate through atoms
        corr = numpy.zeros(uij_diff.shape[:2])
        for i in xrange(res.shape[0]):
            # Calculate correlations
            corr_vals = numpy.corrcoef(uij_diff[:,i,:], res[i])[-1,:-1]
            # If nans, set to tiny noise
            if numpy.isnan(corr_vals).all():
                corr_vals = numpy.random.randn(len(corr_vals)) * 0.0000001
            # Add to array
            corr[:,i] = corr_vals

        # Extract atom labels and model labels
        m_h = self.blank_master_hierarchy().select(flex.bool(self.atom_mask.tolist()), copy_atoms=True)
        atm_labels = numpy.array([ShortLabeller.format(a) for a in m_h.atoms()])
        res_labels = numpy.array([ShortLabeller.format(a.parent().parent()) for a in m_h.atoms()])
        mdl_labels = numpy.array([mdl.tag for mdl in self.models])

        # Extract residue boundaries
        residue_boundaries = self.partition_boundaries_custom(atom_labels=res_labels, mask=self.atom_mask.tolist())
        residue_boundaries = residue_boundaries.select(flex.bool(self.atom_mask.tolist()))
        # Make violin plot of the correlations
        for i_group, sel, fitter in self.fitter.levels[0]:
            self.log.bar()
            self.log('Group {}'.format(i_group))
            self.log.bar()
            # Extract the uij rmsds for this group
            corr_sel = corr[:,sel]
            # Extract the residue selections for this group
            res_labels_sel = res_labels[sel]
            _, u_res_idxs = numpy.unique(res_labels_sel, return_index=True)
            g_res_labs = res_labels_sel[sorted(u_res_idxs)]
            g_res_sels = [res_labels_sel==l for l in g_res_labs]
            # Create boxplots of the rmsds by atom
            self.log('Atom-by-atom plots')
            for i_a in xrange(0, sum(sel), max_x_width):
                filename = os.path.join(cor_dir, 'residual-correction_level-0_group_{:02d}_atoms-{:06d}-{:06d}.png'.format(i_group, i_a+1, i_a+max_x_width))
                self.log('\t> {}'.format(filename))
                self.plot.violinplot(filename=filename,
                                     y_vals=corr_sel[:,i_a:i_a+max_x_width],
                                     x_labels=atm_labels[sel].tolist()[i_a:i_a+max_x_width],
                                     title='correlations of residual and input-fitted Uijs (by atom)',
                                     x_lab='atom', y_lab='correlation', rotate_x_labels=True,
                                     x_lim=(0,max_x_width+1), y_lim=(-1, 1),
                                     vlines=(numpy.where(numpy.array(list(residue_boundaries.select(flex.bool(sel.tolist())).atoms().extract_b()[i_a:i_a+max_x_width]), dtype=bool))[0] + 1.5))
        # Create as table and write
        corr_table = pandas.DataFrame(data=corr, index=mdl_labels, columns=atm_labels)
        filename = os.path.join(cor_dir, 'correlation_uij_residual_to_fitting_residual.csv')
        self.log('\nWriting: {}'.format(filename))
        corr_table.to_csv(filename)

    def refine_fitted_dataset_models(self, suffix='-refined'):
        """Refine coordinates of the fitted structures"""

        refine_phenix.auto = False

        proc_args = []
        for mdl in self.models:
            if not os.path.exists(mdl.o_mtz): rel_symlink(mdl.i_mtz, mdl.o_mtz)
            obj = refine_phenix(pdb_file=mdl.o_pdb, mtz_file=mdl.o_mtz, cif_files=self.cifs,
                                out_prefix=os.path.splitext(mdl.o_pdb)[0]+suffix,
                                strategy='individual_sites+occupancies', n_cycles=5)
            obj.tag = mdl.tag
            proc_args.append(obj)

        # Refine all of the models
        self.log('Running refinements')
        refined = libtbx.easy_mp.pool_map(processes=self._n_cpu, func=wrapper_run, args=proc_args, chunksize=1)

        for mdl, ref in zip(self.models, proc_args):
            assert mdl.tag == ref.tag
            mdl.r_pdb = ref.out_pdb_file
            mdl.r_mtz = ref.out_mtz_file
            assert os.path.exists(mdl.r_pdb), '{} does not exist'.format(mdl.r_pdb)
            assert os.path.exists(mdl.r_mtz), '{} does not exist'.format(mdl.r_mtz)

    def generate_fitted_table_ones(self):
        """Write table-ones for each structure before and after fitting"""

        table_one_csv_orig = self.file_manager.add_file(file_name='table_one_input.csv',   file_tag='table_one_input',   dir_tag='table_ones')
        table_one_csv_fitd = self.file_manager.add_file(file_name='table_one_fitted.csv',  file_tag='table_one_fitted',  dir_tag='table_ones')
        table_one_csv_refd = self.file_manager.add_file(file_name='table_one_refined.csv', file_tag='table_one_refined', dir_tag='table_ones')

        for mdl in self.models:
            if not os.path.exists(mdl.o_mtz): rel_symlink(mdl.i_mtz, mdl.o_mtz)

            assert  os.path.exists(mdl.i_pdb) and \
                    os.path.exists(mdl.i_mtz) and \
                    os.path.exists(mdl.o_pdb) and \
                    os.path.exists(mdl.o_mtz)

        output_eff_orig = table_one_csv_orig.replace('.csv', '.eff')
        output_eff_fitd = table_one_csv_fitd.replace('.csv', '.eff')
        output_eff_refd = table_one_csv_refd.replace('.csv', '.eff')

        # Populate table one phil
        phil = multi_table_ones.master_phil.extract()
        phil.input.dir        = []
        phil.options          = self.params.table_ones_options
        phil.settings.cpus    = self.params.settings.cpus
        phil.settings.verbose = False

        # Run 1
        phil.input.pdb = [mdl.i_pdb for mdl in self.models]
        phil.input.labelling  = self.params.input.labelling
        phil.output.parameter_file = output_eff_orig
        phil.output.output_basename = os.path.splitext(output_eff_orig)[0]
        multi_table_ones.run(params=phil)

        # Run 2
        phil.input.pdb = [mdl.o_pdb for mdl in self.models]
        phil.input.labelling = 'foldername'
        phil.output.parameter_file = output_eff_fitd
        phil.output.output_basename = os.path.splitext(output_eff_fitd)[0]
        multi_table_ones.run(params=phil)

        # Run 3
        if self.models[0].r_pdb is not None:
            phil.input.pdb = [mdl.r_pdb for mdl in self.models]
            phil.input.labelling = 'foldername'
            phil.output.parameter_file = output_eff_refd
            phil.output.output_basename = os.path.splitext(output_eff_refd)[0]
            multi_table_ones.run(params=phil)

        self.log('Running phenix.table_one to calculate R-factors')
        for f in [output_eff_orig,output_eff_fitd,output_eff_refd]:
            if not os.path.exists(f): continue
            self.log.bar()
            cmd = CommandManager('phenix.table_one')
            cmd.add_command_line_arguments([f])
            cmd.print_settings()
            cmd.run()
            self.log('')
            cmd.write_output(f.replace('.eff','.log'))
        self.log.bar()

        # Clear all of the symlinks
        #for mdl in self.models:
        #    if os.path.islink(mdl.o_mtz):
        #        os.remove(mdl.o_mtz)

        assert os.path.exists(table_one_csv_orig)
        assert os.path.exists(table_one_csv_fitd)
        if self.models[0].r_pdb is not None:
            assert os.path.exists(table_one_csv_refd)

    def write_combined_csv(self, uij_fit, uij_inp, out_dir='./'):
        """Add data to CSV and write"""

        out_dir = easy_directory(out_dir)

        # Create output table
        main_table = pandas.DataFrame(index=[m.tag for m in self.models],
                                      columns=['High Resolution Limit',
                                               'Low Resolution Limit',
                                               'R-free Change (Fitted-Input)',
                                               'R-work Change (Fitted-Input)',
                                               'R-gap Change (Fitted-Input)',
                                               'R-free Change (Refined-Input)',
                                               'R-work Change (Refined-Input)',
                                               'R-gap Change (Refined-Input)',
                                               'Average B-factor Change (Fitted-Input)',
                                               'Average B-factor (fitted atoms) Change (Fitted-Input)',
                                     ])

        # columns to be extracted from each file
        input_cols = ['High Resolution Limit', 'Low Resolution Limit', 'Unique reflections','Completeness (%)','Wilson B-factor']
        other_cols = ['R-work','R-free','Average B-factor']
        other_cols_sort = ['R-free','R-work','R-gap','Average B-factor']

        # Extract data from the table one CSVs
        self.log.subheading('Looking for table one data')
        for suff, csv in [(' (Input)',      self.file_manager.get_file('table_one_input')),
                          (' (Fitted)',     self.file_manager.get_file('table_one_fitted')),
                          (' (Refined)',    self.file_manager.get_file('table_one_refined'))]:
            if not os.path.exists(csv):
                if 'refined' in suff.lower():
                    continue
                raise Exception('Cannot read: {}'.format(csv))
            self.log('Reading: {}'.format(csv))
            table_one = pandas.read_csv(csv, index_col=0, dtype=str).transpose()
            # Separate high- and low-resolution limits
            table_one['Low Resolution Limit'], table_one['High Resolution Limit'] = zip(*table_one['Resolution range'].apply(lambda x: x.split('(')[0].split('-')))
            # Select columns for export
            if 'input' in suff.lower():
                table_one = table_one[input_cols+other_cols]    # Pull many columns out of the first file
            else:
                table_one = table_one[other_cols]               # Only a subset of columns from the other files
            # Remove "high-resolution shell" statistics
            for col in table_one.columns:
                self.log('> Formatting col: {}'.format(col))
                table_one[col] = table_one[col].apply(lambda x: x.split('(')[0])
            # Redetermine data types
            table_one = table_one.apply(lambda x: pandas.to_numeric(x, errors='coerce'))
            # If first file, transfer constant statistics
            if 'input' in suff.lower():
                # Transfer selected columns directly to main table
                for col in input_cols:
                    main_table[col] = table_one[col]
                # Reselect the remaining columns
                table_one = table_one[other_cols]
            # Calculate additional columns
            table_one['R-gap'] = table_one['R-free'] - table_one['R-work']
            # Resort columns
            table_one = table_one[other_cols_sort]
            # Add suffix
            table_one.columns = table_one.columns + suff
            # Transfer data to other
            main_table = main_table.join(table_one, how="outer")
            self.log('')
        # Create columns for the deltas between variables
        self.log('Creating additonal columns:')
        for delta_col, col in [('R-work Change', 'R-work'),
                               ('R-free Change', 'R-free'),
                               ('R-gap Change', 'R-gap'),
                               ('Average B-factor Change', 'Average B-factor')]:
            for suff1, suff2 in [(' (Input)', ' (Fitted)'),
                                 (' (Input)', ' (Refined)')]:
                if (col+suff1 not in main_table.columns) or (col+suff2 not in main_table.columns):
                    continue
                suff = ' ({}-{})'.format(suff2.strip(' ()'), suff1.strip(' ()'))
                self.log('> Creating col "{}" = "{}" - "{}"'.format(delta_col+suff, col+suff2, col+suff1))
                main_table[delta_col+suff] = main_table[col+suff2] - main_table[col+suff1]

        self.log.bar(True, False)
        self.log('Calculating mean and median fitting rmsds by dataset')
        self.log.bar()
        # Calculate rmsd between input and fitted uijs
        uij_rmsd = rms(uij_inp-uij_fit, axis=2)
        # Extract mean/median dataset-by-dataset RMSDs
        dset_medn_rmsds = numpy.median(uij_rmsd, axis=1)
        dset_mean_rmsds = numpy.mean(uij_rmsd, axis=1)
        main_table['Mean Fitting RMSD']   = dset_mean_rmsds
        main_table['Median Fitting RMSD'] = dset_medn_rmsds
        for i in xrange(0, min(10, len(self.models))):
            self.log('Model {:10}: {:6.3f} (mean), {:6.3f} (median)'.format(self.models[i].tag, dset_mean_rmsds[i], dset_medn_rmsds[i]))

        self.log.bar(True, False)
        self.log('Calculating isotropic ADPs for input and fitted ADPs by dataset')
        self.log.bar()
        # Calculate isotropic ADPs for input and fitted uijs
        uij_inp_iso = numpy.array(map(uij_to_b, uij_inp))
        uij_fit_iso = numpy.array(map(uij_to_b, uij_fit))
        # Calculate mean/median ADPs for each atom
        dset_mean_inp_iso = numpy.mean(uij_inp_iso, axis=1)
        dset_mean_fit_iso = numpy.mean(uij_fit_iso, axis=1)
        main_table['Average B-factor (fitted atoms) (Input)']               = dset_mean_inp_iso
        main_table['Average B-factor (fitted atoms) (Fitted)']              = dset_mean_fit_iso
        main_table['Average B-factor (fitted atoms) Change (Fitted-Input)'] = dset_mean_fit_iso - dset_mean_inp_iso
        for i in xrange(0, min(10, len(self.models))):
            self.log('Model {:10}: {:6.3f} (input) -> {:6.3f} (fitted)'.format(self.models[i].tag, dset_mean_inp_iso[i], dset_mean_fit_iso[i]))

        # Remove unpopulated columns
        main_table = main_table.dropna(axis='columns', how='all')

        # Write output csv
        filename = os.path.join(out_dir, 'dataset_scores.csv')
        self.log('')
        self.log('Writing output csv: {}'.format(filename))
        main_table.to_csv(filename)

        # Store table for HTML output
        self.tables.statistics = main_table

        # Make graphs for the table
        self.table_one_summaries(table=main_table, out_dir=os.path.join(out_dir,'graphs'))
        self.table_one_graphs(table=main_table,    out_dir=os.path.join(out_dir,'graphs'))

    def table_one_summaries(self, table, out_dir='./'):
        """Look at pre- and post-fitting statistics"""

        out_dir = easy_directory(out_dir)

        self.log.subheading('Model improvement summary')

        # Extract Column averages (means and medians)
        table_means = table.mean().round(3)
        out_str = '{:>35s} | {:>15} | {:>15} | {:>15}'

        self.log.bar()
        self.log('Dataset Averages:')
        self.log.bar()
        # Columns without old/new prefix
        for lab in table_means.index:
            if lab.lower().endswith('-input)'): continue
            if lab.lower().endswith('(input)'): continue
            if lab.lower().endswith('(fitted)'): continue
            if lab.lower().endswith('(refined)'): continue
            self.log(out_str.format(lab, table_means[lab], '', ''))
        self.log.bar()
        # Columns with old/new prefix
        for comp_lab in ['Fitted','Refined']:
            comp_suff = ' ({})'.format(comp_lab)
            # Check whether the columns are here - check the R-free as always should be present
            if 'R-free ({})'.format(comp_lab) not in table.columns:
                continue
            self.log('Comparing Input and {} Uijs'.format(comp_lab))
            self.log(out_str.format('', 'Input', comp_lab, 'Difference'))
            for lab in table_means.index:
                if not lab.endswith(comp_suff):
                    continue
                lab_base = lab[:lab.index(comp_suff)]
                self.log(out_str.format(lab_base,
                                        table_means[lab_base+' ({})'.format('Input')],
                                        table_means[lab_base+' ({})'.format(comp_lab)],
                                        table_means[lab_base+' Change ({}-{})'.format(comp_lab, 'Input')]))
            self.log.bar()

        return

    def table_one_graphs(self, table, out_dir='./'):
        """Look at pre- and post-fitting graphs"""

        out_dir = easy_directory(out_dir)

        self.log.subheading('Model improvement graphs')

        labs = ['Input']   * ('R-free (Input)'   in table.columns) + \
               ['Fitted']  * ('R-free (Fitted)'  in table.columns) + \
               ['Refined'] * ('R-free (Refined)' in table.columns)

        # ------------------------------------------------------------>
        # Raw scores
        # ------------------------------------------------------------>
        # Histogram of the R-FREE for input+fitted+refined B-factors
        filename = os.path.join(out_dir, 'r-free-by-resolution.png')
        self.log('> {}'.format(filename))
        self.plot.binned_boxplot(filename = filename,
                                 x = table['High Resolution Limit'],
                                 y_vals = [100*table['R-free ({})'.format(l)] for l in labs],
                                 legends = labs,
                                 title = 'R-free for {} and {} fitted B-factors'.format(', '.join(labs[:-1]), labs[-1]),
                                 x_lab = 'Resolution (A)',
                                 y_lab = 'R-free (%)',
                                 rotate_x_labels = True,
                                 min_bin_width = 0.1)
        # Histogram of the R-WORK for input+fitted+refined B-factors
        filename = os.path.join(out_dir, 'r-work-by-resolution.png')
        self.log('> {}'.format(filename))
        self.plot.binned_boxplot(filename = filename,
                                 x = table['High Resolution Limit'],
                                 y_vals = [100*table['R-work ({})'.format(l)] for l in labs],
                                 legends = labs,
                                 title = 'R-work for {} and {} B-factors'.format(', '.join(labs[:-1]), labs[-1]),
                                 x_lab = 'Resolution (A)',
                                 y_lab = 'R-work (%)',
                                 rotate_x_labels = True,
                                 min_bin_width = 0.1)
        # Histogram of the R-GAP for input+fitted+refined B-factors
        filename = os.path.join(out_dir, 'r-gap-by-resolution.png')
        self.log('> {}'.format(filename))
        self.plot.binned_boxplot(filename = filename,
                                 x = table['High Resolution Limit'],
                                 y_vals = [100*table['R-gap ({})'.format(l)] for l in labs],
                                 legends = labs,
                                 title = 'R-gap for {} and {} B-factors'.format(', '.join(labs[:-1]), labs[-1]),
                                 x_lab = 'Resolution (A)',
                                 y_lab = 'R-gap (%)',
                                 rotate_x_labels = True,
                                 min_bin_width = 0.1)
        # ------------------------------------------------------------>
        # Difference scores
        # ------------------------------------------------------------>
        # Histogram of the R-GAP for input and fitted B-factors
        filename = os.path.join(out_dir, 'r-values-change-by-resolution.png')
        self.log('> {}'.format(filename))
        self.plot.binned_boxplot(filename = filename,
                                 x = table['High Resolution Limit'],
                                 y_vals = [100*table['R-free Change (Fitted-Input)'],
                                           100*table['R-work Change (Fitted-Input)'],
                                           100*table['R-gap Change (Fitted-Input)']],
                                 legends = ['R-free change','R-work change','R-gap change'],
                                 title = 'R-value changes between input and fitted B-factors',
                                 x_lab = 'Resolution (A)',
                                 y_lab = 'R-value change (%)',
                                 rotate_x_labels = True,
                                 min_bin_width = 0.1,
                                 hlines = [0.0])
        # ------------------------------------------------------------>
        # Correlations between variables
        # ------------------------------------------------------------>
        # Scatter of resolution and fit RMSDs
        filename = os.path.join(out_dir, 'resolution-vs-rmsds-and-bfactors.png')
        self.log('> {}'.format(filename))
        self.plot.multi_scatter(filename = filename,
                                x = table['High Resolution Limit'],
                                y_vals = [table['Mean Fitting RMSD'],
                                          table['Median Fitting RMSD'],
                                          table['Average B-factor (fitted atoms) (Input)'],
                                          table['Average B-factor (fitted atoms) (Fitted)']],
                                x_lab = 'Resolution (A)',
                                y_lab = ['Mean RMSD',
                                         'Median RMSD',
                                         'Mean Input B-factor',
                                         'Mean Fitted B-factor'],
                                title = 'Resolution v. fitting metrics',
                                shape = (2,2))
        # Scatter of R-works and fit RMSDs
        filename = os.path.join(out_dir, 'r-values-vs-rmsds-and-bfactors.png')
        self.log('> {}'.format(filename))
        self.plot.multi_scatter(filename = filename,
                                x_vals = [table['Mean Fitting RMSD'],
                                          table['Median Fitting RMSD'],
                                          table['Average B-factor (fitted atoms) (Input)'],
                                          table['Average B-factor (fitted atoms) (Fitted)']],
                                y_vals = [table['R-free Change (Fitted-Input)'],
                                          table['R-free Change (Fitted-Input)'],
                                          table['R-free Change (Fitted-Input)'],
                                          table['R-free Change (Fitted-Input)'],
                                          table['R-gap Change (Fitted-Input)'],
                                          table['R-gap Change (Fitted-Input)'],
                                          table['R-gap Change (Fitted-Input)'],
                                          table['R-gap Change (Fitted-Input)']],
                                x_lab = ['Mean RMSD',
                                         'Median RMSD',
                                         'Mean Input B-factor',
                                         'Mean Fitted B-factor'],
                                y_lab = ['R-change'],
                                legends = ['R-free Change',
                                           'R-free Change',
                                           'R-free Change',
                                           'R-free Change',
                                           'R-gap Change',
                                           'R-gap Change',
                                           'R-gap Change',
                                           'R-gap Change'],
                                shape = (2,2))
        return

class MultiDatasetUijPlots(object):

    @staticmethod
    def failure_graph(filename, exception=None, title=None, **args):
        fig = pyplot.figure()
        if title is None: title=filename
        if exception is not None: title += '\n---\n'+str(exception) + '\n---'
        pyplot.title('Failed to make {}'.format( title ))
        fig.savefig(filename)
        pyplot.close(fig)

    @staticmethod
    def bin_x_values(data, n_bins=10):
        """Generate bins for data and return binned values and indices of binning"""

        data = numpy.array(data)
        bins = numpy.linspace(data.min(), data.max(), n_bins+1)
        # To avoid extra bin for the max value
        bins[-1] += 0.001*abs(bins[-1])
        # Assign data to bins
        indices = numpy.digitize(data, bins)
        # Generate bin labels
        bin_labels = ['{:.2f} - {:.2f}'.format(bins[i],bins[i+1]) for i in xrange(n_bins)]
        return bins, indices, bin_labels

    @classmethod
    def binned_boxplot(cls, filename,
                       x, y=None, y_vals=None, legends=None,
                       title='', x_lab='x', y_lab='y',
                       rotate_x_labels=True,
                       max_bins=10, min_bin_width=None,
                       hlines=[], vlines=[]):
        """Generate a binned boxplot from data (or array of data)"""

        assert [(y is None),(y_vals is None)].count(False) == 1, 'must provide y OR y_vals'
        if y is not None: y_vals = [y]
        y_vals = [numpy.array(y) for y in y_vals]
        n_y = len(y_vals)

        if legends is not None:
            assert len(legends) == n_y

        assert isinstance(max_bins, int) and max_bins>0

        if min_bin_width is not None:
            assert min_bin_width > 0
            min_x = numpy.min(x)
            max_x = numpy.max(x)
            n_bins = int(min(max_bins, 1+((max_x-min_x)//min_bin_width)))
        else:
            n_bins = max_bins
        assert n_bins >= 1

        # Stupid list of things to not cut off of the graph
        extra_artists = []

        # Generate binning for the x axis
        bins, indices, bin_labels = cls.bin_x_values(data=x, n_bins=n_bins)
        # Sort the y_vals into the bins
        binned_y = [[y[indices==i] for i in xrange(1, n_bins+1)] for y in y_vals]
        # Width of the boxplot bars
        bar_width = 0.5/n_y
        # Colors of each of y_vals
        colors = pyplot.cm.rainbow(numpy.linspace(0,1,n_y))

        # Create figures
        fig, axis = pyplot.subplots(nrows=1, ncols=1)
        axis.set_title(title)
        # Draw horizontal/vertical lines (do first so they're at the bottom)
        for v in hlines: axis.axhline(y=v, linewidth=2, zorder=1)
        for v in vlines: axis.axvline(x=v, linewidth=2, zorder=1)
        # Store plot objects
        plot_dicts = []
        for i_y, y in enumerate(binned_y):
            # Offset for this bar set relative to (n-1)
            x_offset = 0.5 + (1+2*i_y)*bar_width
            positions = numpy.arange(n_bins) + x_offset
            plt = axis.boxplot(y, positions=positions,
                               widths=bar_width,
                               showmeans=False,
                               patch_artist=True)
            plot_dicts.append(plt)
            # Color the boxplots
            c = colors[i_y]
            for patch in plt['boxes']:
                patch.set_facecolor(c)
            for line in plt['medians']:
                line.set_color('k')

        # Make labels
        bin_centres = numpy.arange(n_bins)+1
        axis.set_xticks(bin_centres)
        axis.set_xticklabels(bin_labels)
        # Axis labels
        axis.set_xlabel(x_lab)
        axis.set_ylabel(y_lab)
        # Axis limits
        axis.set_xlim((0, n_bins+1))
        # Plot v_lines
        if n_y > 1:
            for i in xrange(n_bins+1):
                axis.axvline(i+0.5, ls="dotted")
        # X-axis rotations
        if rotate_x_labels:
            pyplot.setp(axis.get_xticklabels(), rotation=45)
        # Make legend
        if legends is not None:
            handles = [patches.Patch(color=colors[i_y], label=legends[i_y]) for i_y in xrange(n_y)]
            lgd = axis.legend(handles=handles,
                              bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
                              #bbox_to_anchor=(0.0, 1.02, 1.0, 0.102),
                              #loc=3, ncol=n_y, mode="expand", borderaxespad=0.0)
            extra_artists.append(lgd)
        fig.tight_layout()
        fig.savefig(filename,
                    bbox_extra_artists=extra_artists,
                    bbox_inches='tight',
                    dpi=300)
        pyplot.close(fig)

        return

    @classmethod
    def multi_scatter(cls, filename,
                      x=None, x_vals=None, y=None, y_vals=None,
                      legends=None,
                      title='', x_lab='x', y_lab='y',
                      x_lim=None, y_lim=None,
                      shape=None,
                      hlines=[], vlines=[]):
        """Generate a scatter plot from data (or list of data)"""

        assert [(x is None),(x_vals is None)].count(False) == 1, 'must provide x OR x_vals'
        assert [(y is None),(y_vals is None)].count(False) == 1, 'must provide y OR y_vals'
        if x is not None: x_vals = [x]
        if y is not None: y_vals = [y]
        x_vals = [numpy.array(x) for x in x_vals]
        y_vals = [numpy.array(y) for y in y_vals]
        n_x = len(x_vals)
        n_y = len(y_vals)
        # Share X- or Y-axes?
        share_x = (n_x==1)
        share_y = (n_y==1)
        # Allow for differently sized lists if one is of unitary length
        #assert (n_x==1) or (n_y==1) or (n_x==n_y), 'number of x_vals not equal to number of y_vals'
        # Create cycles for simplicity
        x_vals = itertools.cycle(x_vals)
        y_vals = itertools.cycle(y_vals)
        # Create cycles of the x/y-labels
        x_lab = itertools.cycle([x_lab]) if isinstance(x_lab, str) else itertools.cycle(x_lab)
        y_lab = itertools.cycle([y_lab]) if isinstance(y_lab, str) else itertools.cycle(y_lab)

        # Create one plot if no shape specified
        if shape is None: shape = (1,1)

        n_line = max(n_x, n_y)
        n_plot = numpy.product(shape)

        if legends is not None:
            assert len(legends) == n_line

        # Colors of each of y_vals
        if n_plot == n_line:
            colors = ['b']*n_line
        elif (n_line % n_plot) == 0:
            colors = sorted(numpy.concatenate([pyplot.cm.rainbow(numpy.linspace(0,1,n_line//n_plot))]*n_plot).tolist())
        else:
            colors = pyplot.cm.rainbow(numpy.linspace(0,1,n_line))
        assert len(colors)==n_line

        # Create figures
        fig, axes = pyplot.subplots(nrows=shape[0], ncols=shape[1],
                                    sharex=share_x, sharey=share_y)
        # Create list if only one plot
        if shape == (1,1): axes = [axes]
        else: axes = numpy.array(axes).flatten()
        # Set the figure title
        fig.suptitle(title)
        # Draw horizontal/vertical lines (do first so they're at the bottom)
        for axis in axes:
            for v in hlines: axis.axhline(y=v, linewidth=2, zorder=1)
            for v in vlines: axis.axvline(x=v, linewidth=2, zorder=1)
        # Store plot objects
        plot_dicts = []
        axes_cycle = itertools.cycle(axes)
        for i in xrange(n_line):
            axis = axes_cycle.next()
            plt = axis.scatter(x=x_vals.next(),
                               y=y_vals.next(),
                               c=colors[i])
            if legends: plt.set_label(legends[i])
            plot_dicts.append(plt)
        # Format axes
        for axis in axes:
            # Axis labels
            axis.set_xlabel(x_lab.next())
            axis.set_ylabel(y_lab.next())
            # Make axis legends
            if legends:
                axis.legend()
        # Axis limits
        if x_lim: axis.set_xlim(x_lim)
        if y_lim: axis.set_ylim(y_lim)
        # X-axis rotations
#        if rotate_x_labels:
#            pyplot.setp(axis.get_xticklabels(), rotation=45)
        fig.tight_layout()
        fig.savefig(filename)#, dpi=300)
        pyplot.close(fig)

        return

    @staticmethod
    def multi_histogram(filename, x_vals, titles, x_labs, rotate_x_labels=True, shape=None, n_bins=30):
        """Generate standard histogram"""

        if shape is not None:
            nrow, ncol = shape
        else:
            nrow, ncol = (1,len(x_vals))

        xlim = [numpy.min(x_vals), numpy.max(x_vals)]
        if (xlim[1] - xlim[0]) < 0.1:
            xlim[0] -= 0.5
            xlim[1] += 0.5

        fig, axes = pyplot.subplots(nrows=nrow, ncols=ncol, sharey=False)
        for i, axis in enumerate(axes.flatten()):
            axis.set_title(titles[i])
            axis.hist(x=x_vals[i], bins=n_bins, range=xlim)
            axis.set_xlabel(x_labs[0])
            axis.set_ylabel('Count')
            axis.set_xlim(xlim)
            if rotate_x_labels:
                pyplot.setp(axis.get_xticklabels(), rotation=90)
        fig.tight_layout()
        fig.savefig(filename)#, dpi=300)
        pyplot.close(fig)

        return

    @staticmethod
    def boxplot(filename, y_vals, x_labels,
                title, x_lab='x', y_lab='y',
                x_lim=None, y_lim=None,
                rotate_x_labels=True,
                hlines=[], vlines=[]):
        """Generate standard boxplot"""

        fig, axis = pyplot.subplots(nrows=1, ncols=1)
        axis.set_title(title)
        axis.boxplot(y_vals, labels=x_labels, showmeans=True)
        for v in hlines: axis.axhline(v)
        for v in vlines: axis.axvline(v)
        axis.set_xlabel(x_lab)
        axis.set_ylabel(y_lab)
        axis.set_xlim(x_lim)
        axis.set_ylim(y_lim)
        if rotate_x_labels:
            pyplot.setp(axis.get_xticklabels(), rotation=90)
        fig.tight_layout()
        fig.savefig(filename)#, dpi=300)
        pyplot.close(fig)

        return

    @staticmethod
    def violinplot(filename,
                   y_vals,
                   x_labels,
                   title,
                   x_lab='x',
                   y_lab='y',
                   x_lim=None,
                   y_lim=None,
                   rotate_x_labels=True,
                   hlines=[],
                   vlines=[]):
        """Generate standard violin plot"""

        fig, axis = pyplot.subplots(nrows=1, ncols=1)
        axis.set_title(title)
        axis.violinplot(y_vals, showmeans=True)
        axis.set_xticks(range(1,len(x_labels)+1))
        axis.set_xticklabels(x_labels)
        for v in hlines: axis.axhline(v)
        for v in vlines: axis.axvline(v)
        axis.set_xlabel(x_lab)
        axis.set_ylabel(y_lab)
        axis.set_xlim(x_lim)
        axis.set_ylim(y_lim)
        if rotate_x_labels:
            pyplot.setp(axis.get_xticklabels(), rotation=90)
        fig.tight_layout()
        fig.savefig(filename)#, dpi=300)
        pyplot.close(fig)

        return

    @staticmethod
    def residue_by_residue(prefix, hierarchy, title=None, v_line_hierarchy=None):
        """Write out residue-by-residue b-factor graphs"""

        h = protein(hierarchy)
        for chain_id in [c.id for c in h.chains()]:
            sel = h.atom_selection_cache().selection('chain {}'.format(chain_id))
            sel_h = h.select(sel)
            y_vals = numpy.array([numpy.mean(rg.atoms().extract_b()) for rg in sel_h.residue_groups()])
            if not y_vals.any(): continue # Skip chains with no Bs
            # Create x-values for each residue starting from 1
            x_vals = numpy.array(range(len(list(sel_h.residue_groups()))))+1
            x_labels = ['']+[ShortLabeller.format(rg) for rg in sel_h.residue_groups()]
            filename = prefix + '-chain_{}.png'.format(chain_id)
            fig, axis = pyplot.subplots(nrows=1, ncols=1, sharey=True)
            if title is not None: axis.set_title(label=str(title))
            axis.set_xlabel('Residue')
            axis.set_ylabel('Isotropic B ($\AA^2$)')
            #axis.plot(x_vals, y_vals, '-ko', markersize=2)
            axis.bar(left=(x_vals-0.5), height=y_vals, width=1.0)
            axis.set_xticklabels([x_labels[int(i)] if (i<len(x_labels)) and (float(int(i))==i) else '' for i in axis.get_xticks()])
            pyplot.setp(axis.get_xticklabels(), rotation=90)
            # Plot boundaries
            if v_line_hierarchy:
                v_lines = numpy.where(numpy.array([max(rg.atoms().extract_b()) for rg in v_line_hierarchy.select(sel).residue_groups()], dtype=bool))[0] + 1.5
                for val in v_lines: axis.axvline(x=val, ls='dotted')
            # Format and save
            fig.tight_layout()
            fig.savefig(filename)#, dpi=300)
            pyplot.close(fig)

    @staticmethod
    def level_plots(filename, hierarchies, title, rotate_x_labels=True):
        """Plot a schematic representation of the hierarchical partitioning"""

        fig, axis = pyplot.subplots()

        for i_h, h in enumerate(hierarchies):
            # Extract B-factors and atom labels
            b_vals = numpy.array(h.atoms().extract_b())
            # Iterate through b values and draw boxes for each
            for b in numpy.unique(b_vals):
                # Skip zero or negative values
                if b <= 0: continue
                # Indices of values in this group
                idx = numpy.where((b_vals==b))[0]
                # Cluster the indices to create bars
                if len(idx)==1:
                    cluster = [1]
                else:
                    cluster = scipy.cluster.hierarchy.fclusterdata(X=idx.reshape((len(idx),1)), t=1.1, criterion='distance', metric='euclidean', method='single')
                # Iterate through clusters and draw
                plot_vals = []
                for g in numpy.unique(cluster):
                    min_i = numpy.min(idx[g==cluster])
                    max_i = numpy.max(idx[g==cluster])
                    plot_vals.append((min_i-0.25, max_i-min_i+0.5))
                # Plot
                axis.broken_barh(plot_vals, (i_h+0.5, 1.0))

        # Set axes, etc.
        if title is not None: axis.set_title(label=str(title))
        axis.set_xlabel('Atom')
        axis.set_ylabel('')
        a_labels = numpy.array([ShortLabeller.format(a.parent().parent()) for a in hierarchies[0].atoms()])
#        axis.set_xticklabels([a_labels[int(i)] if ((i>0) and (i<len(a_labels)) and (float(int(i))==i)) else '' for i in axis.get_xticks()])
        axis.set_xticks(range(0, len(a_labels), max(1, len(a_labels)//40)))
        axis.set_xticklabels([a_labels[i] for i in axis.get_xticks()])
        axis.set_yticks(range(1, len(hierarchies)+1))
        axis.set_yticklabels(['Level {}'.format(i) for i in xrange(1, len(hierarchies)+1)])
        axis.set_xlim((0, len(a_labels)+1))
        pyplot.setp(axis.get_xticklabels(), rotation=90)
        # Invert y-axis
        axis.invert_yaxis()
        # Format and save
        fig.tight_layout()
        fig.savefig(filename)#, dpi=300)
        pyplot.close(fig)

    @staticmethod
    def stacked_bar(prefix,
                    hierarchies,
                    legends,
                    title,
                    y_lab='Isotropic B ($\AA^2$)',
                    y_lim=None,
                    v_line_hierarchy=None,
                    rotate_x_labels=True,
                    reverse_legend_order=False):
        """Plot stacked bar plots for a series of hierarchies (plotted values are the B-factors of the hierarchies)"""

        legends = list(legends)
        assert len(hierarchies) == len(legends)

        m_h = hierarchies[0]

        # Check all hierarchies are the same
        for h in hierarchies:
            assert m_h.is_similar_hierarchy(h)

        # Create a plot for each chain
        for chain_id in [c.id for c in m_h.chains()]:

            # Create a selection for each chain and select it in each hierarchy
            sel = m_h.atom_selection_cache().selection('chain {}'.format(chain_id))
            sel_hs = [h.select(sel) for h in hierarchies]
            sel_mh = sel_hs[0]

            # Filename!
            filename = prefix + '-chain_{}.png'.format(chain_id)

            # Create x-values for each residue starting from 1
            x_vals = numpy.array(range(len(list(sel_mh.residue_groups()))))+1
            x_labels = ['']+[ShortLabeller.format(rg) for rg in sel_mh.residue_groups()]
            # Cumulative y-values (for bottoms of bars)
            cuml_y = None

            # Create colors + hatches
            colors = pyplot.cm.rainbow(numpy.linspace(0,1,len(sel_hs)))
            hatchs = itertools.cycle(['//', 'x', '\\'])

            # Create the output figure
            fig, axis = pyplot.subplots(nrows=1, ncols=1)
            if title is not None: axis.set_title(label=str(title))
            axis.set_xlabel('Residue')
            axis.set_ylabel(y_lab)

            # Iterative though hierarchies
            handles = []
            for i_h, h in enumerate(sel_hs):
                # Extract b-factors from this hierarchy
                y_vals = numpy.array([numpy.mean(rg.atoms().extract_b()) for rg in h.residue_groups()])
                # Skip chains with no Bs
                if not y_vals.any():
                    continue
                # Initialise cumulative object
                if cuml_y is None:
                    cuml_y = numpy.zeros_like(y_vals)
                # Plot the bar
                hdl = axis.bar(left=(x_vals-0.5),
                               height=y_vals,
                               bottom=cuml_y,
                               width=1.0,
                               color=colors[i_h],
                               label=legends[i_h],
                               hatch=hatchs.next(),
                               edgecolor='black')
                handles.append(hdl)

                # Append to cumulative y
                cuml_y += y_vals

            # No lines if cuml_y is None
            if cuml_y is None:
                continue

            # Legends (reverse the plot legends!)
            if reverse_legend_order is True: handles.reverse()
            lgd = axis.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            # Axis labels
            axis.set_xticklabels([x_labels[int(i)] if (i<len(x_labels)) and (float(int(i))==i) else '' for i in axis.get_xticks()])
            # Rotate axis labels
            if rotate_x_labels: pyplot.setp(axis.get_xticklabels(), rotation=90)

            # Plot boundaries
            if v_line_hierarchy:
                v_lines = numpy.where(numpy.array([max(rg.atoms().extract_b()) for rg in v_line_hierarchy.select(sel).residue_groups()], dtype=bool))[0] + 1.5
                for val in v_lines: axis.axvline(x=val, ls='dotted')

            # Format and save
            fig.tight_layout()
            fig.savefig(filename,
                        bbox_extra_artists=[lgd],
                        bbox_inches='tight',
                        dpi=300)
            pyplot.close(fig)

    @staticmethod
    def tracking_plots(table, filename):

        fig, axes = pyplot.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
        # Create list if only one plot
        axes = numpy.array(axes).flatten()
        # Set the figure title
        fig.suptitle('Tracking data', y=0.99)

        # Create xpositions and labels
        x_vals = numpy.arange(0, len(table.index))
        x_max = max(x_vals)
        x_labs = [v if v in ['start','reset','residual'] else 'level {}'.format(v) for v in table['level'].values.astype(str).tolist() ]

        # FIRST AXIS
        ax = axes[0]
        handles = []
        # Create RMSD plot
        hdl = ax.plot(x_vals, table['rmsd'].tolist(), 'bo-', label='model')
        handles.extend(hdl)
        # Axis stuff
        ax.xaxis.set_ticks_position('bottom')
        ax.set_ylabel('rmsd to input ($\AA^2$)')
        ax.set_ylim(bottom=0.0)
        # Create legend for axis
        lgd0 = ax.legend(handles=handles, bbox_to_anchor=(0.98, 0.95), loc=1, borderaxespad=0.)

        # SECOND AXIS
        ax = axes[1]
        handles = []
        # Create B-iso lines for all of the levels in the table
        levels = [l for l in table['level'][table['cycle']==table['cycle'][0]] if l not in ['start','reset']]
        for l in levels:
            l_sel = (table['level'] == l)
            l_x_vals = x_vals[l_sel]
            l_y_vals = table['b_iso (level)'][l_sel]
            hdl = ax.bar(left=l_x_vals-0.2, height=l_y_vals, width=0.4, label='level')
            handles.append(hdl)
        # Only want to print once
        handles = handles[:1]
        # Create an overall B-iso line
        hdl = ax.plot(x_vals, table['b_iso (total)'], 'ko-', label='total')
        handles.extend(hdl)
        # Axis stuff
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xlabel('Optimisation stage')
        ax.set_ylabel('B-iso Equivalent ($\AA^2$)')
        ax.set_xticks(x_vals)
        ax.set_xticklabels(x_labs, rotation=90, ha='center')
        ax.set_xlim(left=-0.5, right=max(x_vals)+0.5)
        ax.set_ylim(bottom=0.0)
        # Create legend for axis
        lgd1 = ax.legend(handles=handles, bbox_to_anchor=(0.02, 0.95), loc=2, borderaxespad=0.)

        # BOTH AXES -- Add vertical lines between macro-cycles
        curr_v = -0.5
        for i, v in enumerate(x_vals[(table['level']=='residual')] + 0.5):
            # Dashed lines to separate cycles
            if v < x_max:
                for ax in axes:
                    ax.axvline(x=v, linewidth=1, zorder=1, color='k', linestyle='--')
            # Text to label each cycle
            axes[0].text(x=(v+curr_v)/2.0,
                         y=0.95*axes[0].get_ylim()[0] + 0.05*axes[0].get_ylim()[1],
                         s='cycle {}'.format(i+1),
                         horizontalalignment='center',
                         verticalalignment='bottom',
                        )
            curr_v = v

        fig.tight_layout()
        fig.savefig(filename,
                    bbox_extra_artists=[lgd0,lgd1],
                    bbox_inches='tight')
        pyplot.close(fig)

class MultiDatasetHierarchicalUijFitter(object):

    def __init__(self,
                 observed_uij,
                 observed_xyz,
                 level_array,
                 level_labels=None,
                 params=None,
                 log=None):

#        if log is None: log = Log(verbose=False)
        self.log = log
        self.params = params

        assert observed_uij.shape[1]  == level_array.shape[1]
        assert observed_uij.shape[:2] == observed_xyz.shape[:2]
        assert observed_uij.shape[2]  == 6
        assert observed_xyz.shape[2]  == 3

        # Store observed values (needed for later)
        self.observed_uij = observed_uij
        self.observed_xyz = observed_xyz
        self.observed_com = numpy.mean(observed_xyz, axis=1)

        # Masks to exclude datasets or atoms from optimisation
        self.dataset_mask = numpy.ones(observed_uij.shape[0], dtype=bool)
        self.atomic_mask  = numpy.ones(observed_uij.shape[1], dtype=bool)

        # Level labels, and grouping for each level
        self.level_labels = level_labels if level_labels else ['Level {}'.format(i) for i in xrange(1, len(level_array)+1)]
        self.level_array = level_array

        # Series of objects to fit the Uij TLS groups
        self.levels = []
        for idx, (lab, group_idxs) in enumerate(zip(self.level_labels, self.level_array)):
            # Create fitter for each level
            self.levels.append(MultiDatasetUijTLSGroupLevel(observed_uij = observed_uij,
                                                            observed_xyz = observed_xyz,
                                                            observed_com = self.observed_com,
                                                            group_idxs = group_idxs,
                                                            n_tls = self.params.number_of_modes_per_group,
                                                            index = idx+1,
                                                            label = lab,
                                                            penalties = self.params.penalties,
                                                            log = self.log))

        # One object to fit all the Uij residuals
        self.residual = MultiDatasetUijResidualLevel(observed_uij = observed_uij,
                                                     index = len(self.levels)+1,
                                                     label = 'residual',
                                                     penalties = self.params.penalties,
                                                     log = self.log)

        assert len(self.level_labels) == len(self.level_array)
        assert len(self.level_labels) == len(self.levels)

        self.apply_masks()

    def __iter__(self):
        for i_level, level in enumerate(self.levels):
            yield (i_level+1, self.level_labels[i_level], level)

    def _target_uij(self, fitted_uij_by_level, i_level):
        arr = fitted_uij_by_level.copy()
        arr[i_level] = 0.0
        arr_sum = numpy.sum(arr, axis=0)
        assert arr_sum.shape == self.observed_uij.shape
        return self.observed_uij - arr_sum

    def n_levels(self):
        return len(self.levels)

    def n_params(self, non_zero=False):
        return sum([l.n_params(non_zero=non_zero) for l in self.levels+[self.residual]])

    def n_params_per_atom(self, non_zero=False):
        return self.n_params(non_zero=non_zero) / (self.observed_uij.shape[0] * self.observed_uij.shape[1])

    def n_input_params_per_atom(self):
        return self._input_params_per_atom

    def n_input_values(self):
        return self.observed_uij.shape[0] * self.observed_uij.shape[1] * self.n_input_params_per_atom()

    def parameter_ratio_gain(self, non_zero=False):
        """Total number of input values divided by number of model parameters (used or total)"""
        return self.n_params(non_zero=non_zero) / self.n_input_values()

    def set_input_info(self, disorder_model):
        self.input_disorder_model = disorder_model
        if disorder_model == 'anisotropic':
            self._input_params_per_atom = 6
        elif disorder_model == 'isotropic':
            self._input_params_per_atom = 3
        else:
            raise Exception('Invalid disorder model?!')

    def set_optimisation_datasets(self, dataset_indices):
        self.dataset_mask[:] = False
        self.dataset_mask[dataset_indices] = True

    def set_tracking(self, table, out_dir):
        self.tracking_data = table
        self._tracking_dir = out_dir

    def update_tracking(self, uij_lvl, i_cycle, i_level=None):
        """Update the tracking table"""

        self.log.subheading('Updating tracking...')

        # Extract uijs for all of the levels for all datasets
        uij_tot = uij_lvl.sum(axis=0)
        # Calculate the rms between fitted and input
        rmsd = rms(self.observed_uij-uij_tot, axis=None)

        # Average over all datasets
        uij_tot = uij_tot.mean(axis=0)
        uij_lvl = uij_lvl.mean(axis=1)

        assert uij_lvl.shape == (len(self.levels)+1,) + self.observed_uij.shape[1:]
        assert uij_tot.shape == self.observed_uij.shape[1:]

        # Extract the Uij for the selected level
        if isinstance(i_level, int):
            uij_sel = uij_lvl[i_level]
            level_lab = str(i_level + 1)
        elif i_level == 'residual':
            uij_sel = uij_lvl[-1]
            level_lab = str(i_level)
        else:
            uij_sel = None
            level_lab = str(i_level)

        # Calculate U-iso & B-iso for selected level
        if uij_sel is not None:
            assert uij_sel.shape == self.observed_uij.shape[1:]
            b_iso_sel = numpy.mean(uij_to_b(uij_sel))
            u_iso_sel = b_iso_sel / EIGHT_PI_SQ
        else:
            b_iso_sel = None
            u_iso_sel = None

        # Calculate U-iso & B-iso for complete model
        b_iso_tot = numpy.mean(uij_to_b(uij_tot))
        u_iso_tot = b_iso_tot / EIGHT_PI_SQ

        # Create human-readable cycle number
        cycle_lab = str(i_cycle+1)

        # Add to tracking table
        self.tracking_data.loc[len(self.tracking_data.index)] = {'cycle' : cycle_lab,
                                                                 'level' : level_lab,
                                                                 'rmsd'  : rmsd,
                                                                 'u_iso (level)' : u_iso_sel,
                                                                 'b_iso (level)' : b_iso_sel,
                                                                 'u_iso (total)' : u_iso_tot,
                                                                 'b_iso (total)' : b_iso_tot}

        self.log(self.tracking_data.loc[len(self.tracking_data.index)-1].to_string())
        # Dump to csv
        self.tracking_data.to_csv(os.path.join(self._tracking_dir, 'tracking_data.csv'))

        # Make plots
        MultiDatasetUijPlots.tracking_plots(table=self.tracking_data, filename=os.path.join(self._tracking_dir, 'tracking_data.png'))

    def apply_masks(self):
        self.log.subheading('Setting atom and datasets masks')
        self.log('> Applying dataset masks')
        for level in self.levels+[self.residual]:
            self.log('...level {}'.format(level.index))
            level.set_dataset_mask(self.dataset_mask)
        self.log('> Applying atomic masks')
        for level in self.levels:
            self.log('...level {}'.format(level.index))
            level.set_atomic_mask(self.atomic_mask)

    # TODO make this z-score based on a half-normal distribution?
    def update_dataset_mask(self, percentile, observed_fitted_differences):
        """Identify the datasets with the worst fit and remove from optimisation"""
        assert observed_fitted_differences.shape == self.observed_uij.shape
        # Calculate the average absolute deviation
        fit_err = numpy.abs(observed_fitted_differences).mean(axis=(1,2))
        mask_thresh = numpy.percentile(fit_err, percentile)
        mask = (fit_err < mask_thresh)
        if (sum(mask) > 0) and (sum(mask*self.dataset_mask)>0):
          self.log('> {} datasets removed with mask'.format(sum(self.dataset_mask)-sum(self.dataset_mask*mask)))
          self.dataset_mask *= mask
        else:
          self.log('> No datasets removed with mask')
        return mask

    # TODO make this z-score based on a half-normal distribution?
    def update_atomic_mask(self, percentile, atomic_uij):
        assert atomic_uij.shape == self.observed_uij.shape[1:]
        # Calculate the isotropic equivalent for each atom
        iso_uij = numpy.mean(numpy.abs(atomic_uij[:,:3]),axis=1)
        mask_thresh = numpy.percentile(iso_uij, percentile)
        mask = (iso_uij < mask_thresh)
        if (sum(mask) > 0) and (sum(mask*self.atomic_mask)>0):
          self.log('> {} atoms removed with mask'.format(sum(self.atomic_mask)-sum(self.atomic_mask*mask)))
          self.atomic_mask *= mask
        else:
          self.log('> No atoms removed with mask')
        return mask

    def fit(self, n_cpus=1, n_macro_cycles=1, n_micro_cycles=1):
        """Run macro-cycles of parameter optimisation"""

        self.log('Fitting TLS models for {} levels (+ residual)'.format(len(self.levels)))

        # Cumulative fitted Uij from the different levels (+1 level for the residuals!)
        # Indexing: [level, dataset, atom, uij]
        fitted_uij_by_level = numpy.zeros((len(self.levels)+1,)+self.observed_uij.shape)

        # Run macro cycles of optimisation
        for i_macro in xrange(n_macro_cycles):
            self.log.heading('Macrocycle {} of {}'.format(i_macro+1, n_macro_cycles), spacer=True)

            # Update masks at beginning of each cycles
            if i_macro == 0:
                # Initialise tracking at the beginning of the run
                self.update_tracking(uij_lvl=fitted_uij_by_level, i_cycle=i_macro, i_level='start')
            else:
                self.log.subheading('Updating parameters for next iteration')
                self.log('Removing atoms with high residual uij from TLS optimisation')
                self.update_atomic_mask(percentile=90, atomic_uij=self.residual.extract())
                # XXX Currently disabled -- need to think of better method XXX
                #self.log('Removing datasets with high fit rmsds from TLS optimisation')
                #self.update_dataset_mask(percentile=95, observed_fitted_differences=self.observed_uij-fitted_uij_by_level.sum(axis=0))
                self.log('Resetting residual Uijs to zero for next cycle')
                fitted_uij_by_level[-1] = 0.0
                # Re-initialise tracking at the beginning of each subsequent cycle
                self.update_tracking(uij_lvl=fitted_uij_by_level, i_cycle=i_macro, i_level='reset')

            # Ensure the masks are up-to-date
            self.apply_masks()

            # Iterate through the TLS levels of the fitting
            for i_level, fitter in enumerate(self.levels):
                self.log.subheading('Macrocycle {} of {}: '.format(i_macro+1, n_macro_cycles)+'Fitting TLS Groups (level {} - {})'.format(fitter.index, fitter.label))
                # Update the target uij by subtracting contributions from other levels
                self.log('Updating target Uijs for optimisation')
                fitter.set_target_uij(target_uij=self._target_uij(fitted_uij_by_level=fitted_uij_by_level, i_level=i_level))
                # Optimise
                fitted_uij_by_level[i_level] = fitter.run(n_cpus=n_cpus, n_cycles=n_micro_cycles)
                # Update tracking
                self.update_tracking(uij_lvl=fitted_uij_by_level, i_cycle=i_macro, i_level=i_level)

            # Fit the residuals
            self.log.subheading('Macrocycle {} of {}: '.format(i_macro+1, n_macro_cycles)+'Fitting residual atomic Uijs')
            # Update the target uij by subtracting contributions from other levels
            self.residual.set_target_uij(target_uij=self._target_uij(fitted_uij_by_level=fitted_uij_by_level, i_level=-1))
            # Update fitters and optimise -- always run two cycles of this
            fitted_uij_by_level[-1] = self.residual.run(n_cpus=n_cpus, n_cycles=3)
            # Update tracking
            self.update_tracking(uij_lvl=fitted_uij_by_level, i_cycle=i_macro, i_level='residual')

        return self.extract()

    def extract(self, average=False):
        """Extract the fitted Uij for all structures or the average structure (averaged TLS)"""
        return self.extract_tls(sum_levels=True, average=average) + self.residual.extract()

    def extract_tls(self, sum_levels=False, average=False):
        """Extract the fitted TLS Uij for all structure or the average structure (averaged TLS)"""
        uij = numpy.array([f.extract(average=average) for f in self.levels])
        if sum_levels is True: uij=uij.sum(axis=0)
        return uij

    def summary(self, show=False):
        s = ''
        s += '\n> Input Data:'
        s += '\nNumber of datasets: {}'.format(len(self.dataset_mask))
        s += '\nNumber of atoms:    {}'.format(len(self.atomic_mask))
        s += '\nNumber of input values: {}'.format(self.n_input_values())
        s += '\nInput Uij Type: {}'.format(self.input_disorder_model)
        s += '\n'
        s += '\n> Hierarchical Model:'
        s += '\nNumber of levels (incl. residual): {}'.format(self.n_levels()+1)
        s += '\nNumber of TLS modes per group: {}'.format(self.params.number_of_modes_per_group)
        s += '\nNumber of model parameters: {}'.format(self.n_params())
        s += '\n'
        s += '\n> Optimisation:'
        s += '\nDatasets used for TLS optimisation: {} of {}'.format(sum(self.dataset_mask), len(self.dataset_mask))
        s += '\nAtoms used for TLS optimisation:    {} of {}'.format(sum(self.atomic_mask), len(self.atomic_mask))
        s += '\n'
        s += '\n> Parameterisation:'
        s += '\nNumber of parameters per atom:'
        s += '\n... in input data: {:2.0f}'.format(self.n_input_params_per_atom())
        s += '\n... in fitted model (possible): {:5.2f}'.format(self.n_params_per_atom(non_zero=False))
        s += '\n... in fitted model (non-zero): {:5.2f}'.format(self.n_params_per_atom(non_zero=True))
        s += '\nNon-zero model parameters / number of input data: {:5.2%}'.format(self.parameter_ratio_gain(non_zero=True))
        s += '\n'
        if show: self.log(s)
        return s

#        s += '\n\tUsing {}/{} atoms'.format(len(self.get_atomic_mask()), self._n_atm)
#        s += '\n\tUsing {}/{} datasets'.format(len(self.get_dataset_mask()), self._n_dst)

class _MultiDatasetUijLevel(object):

    _uij_shape = None

    def n_params(self, non_zero=False):
        return sum([f.n_params(non_zero=non_zero) for i,g,f in self])

    def run(self, n_cycles=1, n_cpus=1):
        # Check to see if multiple cpus are available per job
        n_cpus_per_job = max(1, n_cpus//self._n_obj)
        # Create job objects
        jobs = [(fitter, {'n_cycles':n_cycles, 'n_cpus':n_cpus_per_job}) for (i, sel, fitter) in self]
        # Drop CPUs if not required
        n_cpus = min(len(jobs), n_cpus)
        # Only do multiprocessing if actually needed
        if n_cpus > 1:
            # Run jobs in parallel
            self.log('Running {} job(s) in {} process(es) [and {} cpu(s) per process]'.format(len(jobs), n_cpus, n_cpus_per_job))
            workers = NonDaemonicPool(n_cpus)
            finished_jobs = workers.map(func=wrapper_fit, iterable=jobs)
            workers.close()
            #finished_jobs = libtbx.easy_mp.pool_map(processes=n_cpus, func=wrapper_fit, args=jobs, chunksize=1)
        else:
            self.log('Running {} job(s) [with {} cpu(s)]'.format(len(jobs), n_cpus_per_job))
            finished_jobs = [wrapper_fit(j) for j in jobs]
        self.log('')
        self.log('Optimisation complete')
        self.log.subheading('Optimised values')
        # Record list of errors and raise all at end
        errors = []
        for i_iter, (i, sel, fitter) in enumerate(self):
            ret_job = finished_jobs[i_iter]
            if isinstance(ret_job, str) or (ret_job is None):
                errors.append((fitter,ret_job))
                continue
            # Write out the summary to its own log
            ret_job.log.write_to_log(clear_data=True)
            # Print heading for every group for other levels
            if self.label != 'residual':
                self.log.subheading('Results - level {} ({}) - group {}'.format(self.index, self.label, ret_job.label))
            # Print summaries (in main log)
            self.log(ret_job.summary(show=False))
            # Store the returned job
            self.fitters[i] = ret_job
        # Report errors
        if errors:
            for fitter, e in errors:
                self.log.bar()
                self.log('Level "{}", Group "{}": error returned'.format(self.label, fitter.label))
                self.log.bar()
                self.log(e)
            self.log.bar()
            raise Failure('Errors raised during optimisation. Messages printed above')
        # Return the fitted B-factors
        return self.extract()

    def set_target_uij(self, target_uij):
        assert target_uij.shape == self._uij_shape
        for i, sel, fitter in self:
            fitter.set_target_uij(target_uij=target_uij[:,sel])

    def set_atomic_mask(self, mask):
        for i, sel, fitter in self:
            # Apply the mask if contains more than one atoms
            if sum(mask[sel]) > 1:
                fitter.set_atomic_mask(mask[sel])
            # Else do nothing
            else:
                self.log('Level "{}", Group "{}": Attempting to apply mask of one atom or less -- not applying mask'.format(self.label, fitter.label))

    def set_dataset_mask(self, mask):
        for i, sel, fitter in self:
            fitter.set_dataset_mask(mask)

class MultiDatasetUijTLSGroupLevel(_MultiDatasetUijLevel):

    def __init__(self,
                 observed_uij,
                 observed_xyz,
                 observed_com,
                 group_idxs,
                 n_tls=1,
                 index=0,
                 label=None,
                 penalties=None,
                 log=None):

#        if log is None: log = Log(verbose=False)
        self.log = log

        assert observed_uij.shape[:2] == observed_xyz.shape[:2]
        assert observed_uij.shape[2]  == 6
        assert observed_xyz.shape[2]  == 3
        assert observed_uij.shape[1]  == len(group_idxs)
        assert observed_com.shape[0] == observed_uij.shape[0]
        assert observed_com.shape[1] == 3

        self.index = index
        self.label = label if label else 'Level {}'.format(index)

        self.group_idxs = group_idxs

        self._n_groups = sum(numpy.unique(self.group_idxs)!=-1)
        self._n_obj = self._n_groups

        self._uij_shape = observed_uij.shape

        self.fitters = {}
        for i, sel, f in self:
            assert f is None
            # Create a separate log for each fitter
            ls = LogStream(log_file = os.path.splitext(self.log.log_file)[0] + '-level{:04d}-group{:06d}.log'.format(self.index, i),
                           verbose  = (self.log.verbose) or (self._n_groups==1),
                           silent   = not(self._n_groups==1))
            # Create fitter object
            self.fitters[i] = MultiDatasetUijTLSOptimiser(target_uij = observed_uij[:,sel],
                                                          atomic_xyz = observed_xyz[:,sel],
                                                          atomic_com = observed_com,
                                                          n_tls = n_tls,
                                                          label = '{:4d} of {:4d}'.format(i, self._n_groups),
                                                          penalties = penalties,
                                                          log = ls)

    def __iter__(self):
        for i in numpy.unique(self.group_idxs):
            if i==-1: continue
            yield (i, (self.group_idxs==i), self.fitters.get(i, None))

    def n_groups(self):
        return self._n_groups

    def extract(self, average=False):
        fitted_uij = numpy.zeros(self._uij_shape)
        for i, sel, fitter in self:
            fitted_uij[:,sel] = fitter.extract()
        if average:
            fitted_uij = fitted_uij.mean(axis=0)
        return fitted_uij

    def reset_models(self, modes, t=False, l=False, s=False):
        components = 'T'*t + 'L'*l + 'S'*s
        assert components is not ''
        # Iterate through each group/partition
        for i, sel, fitter in self:
            fitter._parameters.reset_models(components=components, modes=modes)

    def reset_amplitudes(self, modes, t=False, l=False, s=False):
        components = 'T'*t + 'L'*l + 'S'*s
        assert components is not ''
        # Iterate through each group/partition
        for i, sel, fitter in self:
            fitter._parameters.reset_amplitudes(components=components, modes=modes)

    def summary(self, show=True):
        num_model_params = self.n_params(non_zero=False)
        num_used_params = self.n_params(non_zero=True)
        n_dst, n_atm, _ = self._uij_shape
        s = '> Level {} - {}'.format(self.index, self.label)
        s += '\n\tNumber of Groups: {}'.format(self.n_groups())
        s += '\n\tNumber of model parameters (total): {}'.format(num_model_params)
        s += '\n\tNumber of model parameters (used): {}'.format(num_used_params)
        s += '\n\tNumber of model parameters per atom (total): {:5.3f}'.format(num_model_params/(n_dst*n_atm))
        s += '\n\tNumber of model parameters per atom (used):  {:5.3f}'.format(num_used_params/(n_dst*n_atm))
        s += '\n'
        if show: self.log(s)
        return s

class MultiDatasetUijResidualLevel(_MultiDatasetUijLevel):

    def __init__(self,
                 observed_uij,
                 index=-1,
                 label=None,
                 penalties=None,
                 log=None):

#        if log is None: log = Log(verbose=True)
        self.log = log

        self.index = index
        self.label = label if label else 'Level {} (residual)'.format(index)

        self._uij_shape = observed_uij.shape

        self._n_atm = self._uij_shape[1]
        self._n_obj = self._n_atm

        # Create one log for all fitters
        ls = LogStream(log_file = os.path.splitext(self.log.log_file)[0]+'-levelXXXX.log',
                       verbose  = False,
                       silent   = False)
        # Create a fitter for each atom
        self.fitters = {}
        for i, sel, f in self:
            assert f is None
            # Create fitter object
            self.fitters[i] = MultiDatasetUijAtomOptimiser(target_uij = observed_uij[:,sel],
                                                           label = 'atom {:5d} of {:5d}'.format(i,self._n_atm),
                                                           penalties = penalties,
                                                           log = ls)

    def __iter__(self):
        for i in xrange(self._n_atm):
            yield (i+1, i, self.fitters.get(i+1, None))

    def extract(self, expanded=False):
        fitted_uij = numpy.zeros(self._uij_shape[1:])
        for i, sel, fitter in self:
            fitted_uij[sel] = fitter.extract()
        if expanded:
            fitted_uij = fitted_uij.reshape((1,)+fitted_uij.shape).repeat(self._uij_shape[0], axis=0)
        return fitted_uij

    def summary(self, show=True):
        num_model_params = self.n_params(non_zero=False)
        num_used_params = self.n_params(non_zero=True)
        n_dst, n_atm, _ = self._uij_shape
        s = '> Level {} - {}'.format(self.index, self.label)
        s += '\n\tNumber of Atoms: {}'.format(n_atm)
        s += '\n\tNumber of model parameters (total): {}'.format(num_model_params)
        s += '\n\tNumber of model parameters (used): {}'.format(num_used_params)
        s += '\n\tNumber of model parameters per atom (total): {:5.3f}'.format(num_model_params/(n_dst*n_atm))
        s += '\n\tNumber of model parameters per atom (used):  {:5.3f}'.format(num_used_params/(n_dst*n_atm))
        s += '\n'
        if show: self.log(s)
        return s

class _UijOptimiser(object):

    def __init__(self,
                 target_uij,
                 atomic_xyz=None,
                 label='',
                 penalties=None,
                 log=None):

        assert penalties is not None

#        if log is None: log = Log(verbose=False)
        self.log = log

        # Store verboseness and silence as overridden in some functions
        self.verbose = self.log.verbose
        self.silent  = self.log.silent

        self._n_prm = 0
        self._n_dst = 0
        self._n_atm = 0

        self.target_uij = target_uij
        self.atomic_xyz = atomic_xyz

        self.label = label

        self.penalties = penalties
        self.penalty = UijPenalties(mdl_fixed=penalties.invalid_tls_values.a,
                                    mdl_slope=penalties.invalid_tls_values.b,
                                    amp_fixed=penalties.invalid_amplitudes.a,
                                    amp_slope=penalties.invalid_amplitudes.b,
                                    uij_fixed=penalties.invalid_uij_values.a,
                                    uij_slope=penalties.invalid_uij_values.b,
                                    ovr_fixed=penalties.over_target_values.a,
                                    ovr_slope=penalties.over_target_values.b)

        self.optimisation_rmsd = numpy.inf
        self.optimisation_penalty = numpy.inf

        # The object to be used for evaluating the target function
        self.evaluator = self

    #===========================================+>
    # Private Functions
    #===========================================+>

    def _blank_atom_selection(self):
        return numpy.zeros(self._n_atm, dtype=bool)

    def _blank_dataset_selection(self):
        return numpy.zeros(self._n_dst, dtype=bool)

    def _blank_parameter_selection(self):
        return numpy.zeros(self._n_prm, dtype=bool)

    def _optimise(self, running_summary=False):
        """Run the optimisation"""

        # Prep variables for target function
        self._n_call = 0
        self._running_summary = running_summary
        # Apply the masks to the target uij
        self._update_target()
        # Initialise the RMSD measure
        self.optimisation_rmsd = 1e6
        self.optimisation_penalty = 0.0
        # Get the variables selected for optimisatiom, and their values
        sel_dict = self._select()
        sel_vals = self._values(selection=sel_dict)
        # Create simplex for these parameters
        opt_simplex = self.simplex.get_simplex(start=sel_vals, **sel_dict)
        # Optimise these parameters
        optimised = simplex.simplex_opt(dimension = len(opt_simplex[0]),
                                        matrix    = map(flex.double, opt_simplex),
                                        evaluator = self.evaluator,
                                        tolerance = 1e-04)
        # Extract and update current values
        self._inject(values=optimised.get_solution(), selection=sel_dict)
        if self._running_summary:
            self.log.bar()

    #===========================================+>
    # Public Functions
    #===========================================+>

    def copy(self):
        return copy.deepcopy(self)

    def get_atomic_mask(self):
        return self._mask_atom

    def get_dataset_mask(self):
        return self._mask_dset

    def set_target_uij(self, target_uij):
        assert target_uij.shape == self.target_uij.shape
        self.target_uij = target_uij

    def set_atomic_mask(self, mask):
        if not isinstance(mask, list):
            mask = list(numpy.where(mask)[0])
        assert len(mask) > 0, 'No atoms in this mask!'
        self._mask_atom = mask

    def set_dataset_mask(self, mask):
        if not isinstance(mask, list):
            mask = list(numpy.where(mask)[0])
        self._mask_dset = mask

    #=======================+>

    def parameters(self):
        return self._parameters

    def target(self, parameters):
        """Target function for the simplex optimisation"""

        # Increment counter
        self._n_call += 1
        # Combine the optimising parameters in the complete parameter set
        self._inject(values=parameters, selection=None)
        # Calculate physical penalties - reject this set if model is not physical
        ppen = self._parameter_penalties()
        # Print info line if necessary
        if self._running_summary:
            if self._n_call%20==1:
                header = '[{}] -> ({:^10}, {:^10})'.format(', '.join(['{:>10}'.format('parameter') for r in parameters]), 'fit/rmsd', 'penalty')
                line = '-'*len(header)
                self.log(line, show=False)
                self.log(header, show=False)
                self.log(line, show=False)
        # Return now if physical penalty if non-zero to save time
        if ppen > 0.0:
            if self._running_summary:
                self.log('[{}] -> ({:>10}, {:10.0f})'.format(', '.join(['{:+10.5f}'.format(r) for r in parameters]), 'UNPHYSICAL', ppen), show=False)
            return ppen
        # Get the fitted uijs (including masked atoms)
        self._update_fitted()
        # Calculate RMSD
        rmsd = numpy.sqrt(numpy.mean(numpy.power(self._target_uij-self._fitted_uij, 2)))
        # Calculate fitting penalties (add to rmsd)
        fpen = self._fitting_penalties(uij_fit=self._fitted_uij, uij_obs=self._target_uij)
        # Total target function value
        tot_val = rmsd + fpen
        # Update minima
        if tot_val < self.optimisation_rmsd+self.optimisation_penalty:
            self.optimisation_rmsd    = rmsd
            self.optimisation_penalty = fpen
        if self._running_summary:
            self.log('[{}] -> ({:10f}, {:10.3f})'.format(', '.join(['{:+10.5f}'.format(r) for r in parameters]), rmsd, fpen), show=False)
        return tot_val

class MultiDatasetUijAtomOptimiser(_UijOptimiser):

    def __init__(self,
                 target_uij,
                 label='',
                 penalties=None,
                 log=None):
        super(MultiDatasetUijAtomOptimiser, self).__init__(target_uij=target_uij,
                                                           label=label,
                                                           penalties=penalties,
                                                           log=log)

        # Should be n_dataset observations of 6 parameters
        assert len(self.target_uij.shape) == 2
        assert self.target_uij.shape[1] == 6

        # Number of datasets
        self._n_dst = self.target_uij.shape[0]
        # Number of parameters
        self._n_prm = 6

        # Output list of Uijs
        self._parameters = numpy.zeros(self._n_prm)

        # Initialise simplex generator
        self.simplex = UijSimplex()
        # Initialse penalty weights
        self.penalty.set_weights(ovr_fixed=0.0, ovr_slope=0.0)

        # Initialise the mask to all datasets
        self.set_dataset_mask(range(self._n_dst))

    #===========================================+>
    # Private Functions - common to parent class
    #===========================================+>

    def _fitting_penalties(self, uij_fit, uij_obs):
        return 0.0

    def _parameter_penalties(self):
        return self.penalty.uij_valid(values=self.result())

    #=======================+>

    def _update_fitted(self):
        self._fitted_uij = self.extract()

    def _update_target(self):
        self._target_uij = self.target_uij[self.get_dataset_mask()]

    #=======================+>

    def _select(self, **kw_args):
        """Define which variables are to be optimised and return selection dictionary"""
        return {}

    def _values(self, **kw_args):
        return self._parameters

    def _inject(self, values, **kw_args):
        """Insert a set of parameters into the complete parameter set"""
        assert len(values) == self._n_prm
        self._parameters[:] = values

    #===========================================+>
    # Public Functions - common to parent class
    #===========================================+>

    def n_params(self, non_zero=False):
        if non_zero is False:
            return numpy.product(self._parameters.shape)
        else:
            return (self.parameters()!=0.0).any()*(self.n_params(non_zero=False))

    def optimise(self, n_cycles=1, n_cpus=1):
        """Optimise the residual for a set of atoms"""
        # If only 1 dataset provided, set current values to target values and return
        if self.target_uij.shape[0] == 1:
            self._inject(values=self.target_uij[0])
            return
        # Otherwise optimise
        for i_cycle in xrange(n_cycles):
            self._optimise()

    def result(self):
        """Return the fitted parameters (same as extract for this class)"""
        return tuple(self._parameters)

    def extract(self):
        """Return the fitted uijs - for all atoms"""
        return tuple(self._parameters)

    def summary(self, show=True):
        """Print the number of parameters/input data"""
        uij = self.result()
        s = 'Uij ({}): '.format(self.label)+', '.join(['{:8.3f}'.format(v) for v in uij])
        if show: self.log(s)
        return s

class MultiDatasetUijTLSOptimiser(_UijOptimiser):

    def __init__(self,
                 target_uij,
                 atomic_xyz,
                 atomic_com,
                 n_tls=1,
                 tls_params=None,
                 label='',
                 penalties=None,
                 log=None):
        super(MultiDatasetUijTLSOptimiser, self).__init__(target_uij=target_uij,
                                                          atomic_xyz=atomic_xyz,
                                                          label=label,
                                                          penalties=penalties,
                                                          log=log)

        # Store the centre of mass of the atoms (for the rotation/screw components)
        self.atomic_com = atomic_com

        # Should be n_dataset observations of n_atm with 6 parameters
        assert len(self.target_uij.shape) == 3
        assert self.target_uij.shape[2] == 6
        assert self.atomic_xyz is not None
        assert self.atomic_com is not None

        self._n_dst = self.target_uij.shape[0]
        self._n_atm = self.target_uij.shape[1]

        assert self.target_uij.shape == (self._n_dst, self._n_atm, 6)
        assert self.atomic_xyz.shape == (self._n_dst, self._n_atm, 3)
        assert self.atomic_com.shape == (self._n_dst, 3)

        # Allow for supplied TLS, or initialise new
        if tls_params is not None:
            inp_tls = tls_params
            num_tls = len(tls_params)
            vec_tls = [p.t+p.l+p.s for p in tls_params]
        elif n_tls is not None:
            inp_tls = None
            num_tls = n_tls
            vec_tls = [numpy.mean(self.target_uij, axis=(0,1)).tolist()+[0.0]*15]+[[0.0]*21]*(num_tls-1)
        else:
            raise Sorry('No TLS models provided')
        assert len(vec_tls) == num_tls
        assert set(map(len, vec_tls)) == {21}

        # Number of models and parameters
        self._n_tls = num_tls
        self._n_prm_mdl = 21 * self._n_tls
        self._n_prm_amp = 03 * self._n_tls * self._n_dst

        # Model-Amplitude parameter sets
        self._parameters = MultiDatasetTLSModelList(n_mdl=self._n_tls,
                                                    n_dst=self._n_dst,
                                                    log = self.log)

        # Extract number of parameters
        self._n_prm = self._parameters.n_params()

        # Optimisation ranges
        self._select(optimise_model = True,
                     optimise_amplitudes = True,
                     components = 'TLS',
                     models   = range(0, self._n_tls),
                     datasets = range(0, self._n_dst),
                     atoms    = range(0, self._n_atm))

        # Initialise simplex generator
        self.simplex = TLSSimplex(mdl_delta=0.1, amp_delta=0.1)

        # Initialse penalty set of test points (for identifying physically-valid TLS matrices)
        box_size = (numpy.min(self.atomic_xyz, axis=(0,1)),
                    numpy.max(self.atomic_xyz, axis=(0,1)))
        box_edge = numpy.array([(box_size[i][0],box_size[j][1],box_size[k][2]) for i,j,k in flex.nested_loop((2,2,2))])
        self.penalty.set_test_xyz(xyz=box_edge, com=self.atomic_com.mean(axis=0))

        # Initialise the masks
        self.set_atomic_mask(range(self._n_atm))
        self.set_dataset_mask(range(self._n_dst))

    #===========================================+>
    # Private Functions - common to parent class
    #===========================================+>

    def _fitting_penalties(self, uij_fit, uij_obs):
        ovr_penalties = []; [ovr_penalties.extend([self.penalty.uij_size(*vv) for vv in zip(*v)]) for v in zip(uij_fit,uij_obs)]
        # Normalise penalties by the number of datasets!
        return numpy.sum(ovr_penalties)*1.0/self._opt_dict['n_dst']

    def _parameter_penalties(self):
        tls_mdl, tls_amp = self.result()
        # Calculate model penalties (if required -- not required if not optimising model...)
        if self._opt_dict['model']:         tls_penalties = [self.penalty.tls_params(values=v) for v in tls_mdl]
        else:                               tls_penalties = []
        if self._opt_dict['amplitudes']:    amp_penalties = [self.penalty.amplitudes(values=v) for v in tls_amp]
        else:                               amp_penalties = []
        return numpy.sum(tls_penalties+amp_penalties)

    #=======================+>

    def _update_fitted(self):
        self._fitted_uij = self._extract(mask_datasets=self.get_dataset_mask(), mask_atoms=self.get_atomic_mask())

    def _update_target(self):
        self._target_uij = self.target_uij[self.get_dataset_mask()][:,self.get_atomic_mask()]

    #=======================+>

    def _select(self, optimise_model=None, optimise_amplitudes=None, components=None, models=None, datasets=None, atoms=None):
        """Set optimisation parameters:
            optimise_model: True or False
                should TLS model components be optimised?
            optimise_amplitudes: True or False
                should TLS amplitudes be optimised?
            components: string
                which model/amplitudes elements should be optimised.
                can only contain letters T or L or S
            models: list of integers
                which TLS models should be optmised?
            datasets: list of integers
                for/using which datasets should parameters be optimised?
            atoms: list of integers
                for/using which atoms should parameters be optimised?
        """

        if optimise_model is not None:
            assert isinstance(optimise_model, bool)
            self._optimise_model = optimise_model
        if optimise_amplitudes is not None:
            assert isinstance(optimise_amplitudes, bool)
            self._optimise_amplitudes = optimise_amplitudes
        if components is not None:
            assert not set(components).difference('TLS')
            self._opt_cpts = components
        if models is not None:
            assert min(models) >= 0
            assert max(models) < self._n_tls
            self._opt_mdls = models
        if datasets is not None:
            assert min(datasets) >= 0
            assert max(datasets) < self._n_dst
            self.set_dataset_mask(datasets)
        if atoms is not None:
            assert min(atoms) >= 0
            assert max(atoms) < self._n_atm
            self.set_atomic_mask(atoms)

        # Create a dictionary that FULLY defines the optimisation cycle
        self._opt_dict = {'model'      : self._optimise_model,
                          'amplitudes' : self._optimise_amplitudes,
                          'components' : self._opt_cpts,
                          'n_cpt' : len(self._opt_cpts),
                          'n_mdl' : len(self._opt_mdls),
                          'n_dst' : len(self.get_dataset_mask()),
                          'i_mdl' : self._opt_mdls,
                          'i_dst' : self.get_dataset_mask()}
#                          'i_atm' : self._opt_atms

        return self._opt_dict

    def _values(self, selection=None):
        """Extract a set of values for a set of variables"""

        # Use default selection if not provided
        if selection is None:
            selection = self._opt_dict
        # Create object from dict
        s = Info(selection)
        # Iterate through objects and extract requested parameter values
        values = []
        if s.model is True:
            for i in s.i_mdl:
                mdl = self._parameters.get(index=i)
                values.append(mdl.model.get(components=s.components))
        if s.amplitudes is True:
            for i in s.i_mdl:
                mdl = self._parameters.get(index=i)
                values.append(mdl.amplitudes.get(components=s.components, datasets=s.i_dst).flatten())
        return numpy.concatenate(values)

    def _inject(self, values, selection=None):
        """Change a set of values for a set of variables"""

        # Use default selection if not provided
        if selection is None:
            selection = self._opt_dict
        # Create object from dict
        s = Info(selection)
        # Iterate through objects and extract requested parameter values
        values = collections.deque(values)
        if s.model is True:
            n_tls_params = str_to_n_params(s.components)
            for i in s.i_mdl:
                mdl = self._parameters.get(index=i)
                mdl.model.set(vals = [values.popleft() for _ in xrange(n_tls_params)],
                              components = s.components)
        if s.amplitudes is True:
            for i in s.i_mdl:
                mdl = self._parameters.get(index=i)
                mdl.amplitudes.set(vals = numpy.array([values.popleft() for _ in xrange(s.n_dst*s.n_cpt)]).reshape((s.n_dst, s.n_cpt)),
                                   components = s.components,
                                   datasets = s.i_dst)
        assert len(values) == 0, 'not all values have been popped'

    #===========================================+>
    # Private Functions - custom for this class
    #===========================================+>

    def _extract(self, mask_datasets=None, mask_atoms=None):
        """Calculate the TLS components for a selection of atoms and datasets"""
        # Get the atomic coordinates
        xyzs = self.atomic_xyz
        coms = self.atomic_com
        # Apply masks - datasets
        if mask_datasets is not None:
            assert isinstance(mask_datasets, list)
            n_dst = len(mask_datasets)
            xyzs = xyzs[mask_datasets]
            coms = coms[mask_datasets]
        else:
            n_dst = self._n_dst
        # Apply masks - atoms
        if mask_atoms is not None:
            assert isinstance(mask_atoms, list)
            n_atm = len(mask_atoms)
            xyzs = xyzs[:,mask_atoms]
        else:
            n_atm = self._n_atm
        # Validate
        assert xyzs.shape == (n_dst, n_atm, 3)
        # Calculate the TLS components
        uij = self._parameters.uij(xyzs=xyzs, origins=coms, datasets=mask_datasets)
        assert uij.shape == (n_dst, n_atm, 6)
        return uij

    #===========================================+>
    # Public Functions
    #===========================================+>

    def n_params(self, non_zero=False):
        return self._parameters.n_params(non_zero=non_zero)

    def optimise(self, n_cycles=1, n_cpus=1):
        """Optimise a (series of) TLS model(s) against the target data"""

        # Extract the masks (so that can be reapplied if changed)
        opt_dset_mask = self.get_dataset_mask()
        opt_atom_mask = self.get_atomic_mask()

        # Cumulative objects
        mdl_cuml = []

        # Workers for multiprocessing
        if n_cpus > 1:
            workers = DaemonicPool(n_cpus)

        # Optimise!
        for i_cycle in xrange(n_cycles):
            self.log.subheading('Group {} - Optimisation cycle {} of {}'.format(self.label, i_cycle+1, n_cycles))

            for i_tls in xrange(self._n_tls):
                self.log.subheading('Optimising TLS model {} of {}'.format(i_tls+1, self._n_tls))
                self.log('Optimising using {} atoms'.format(len(opt_atom_mask)))
                self.log('Optimising using {} datasets'.format(len(opt_dset_mask)))

                # Append to running list
                if i_tls not in mdl_cuml: mdl_cuml.append(i_tls)

                # ---------------------------------->
                # Optimise each T, L, S component separately
                for cpt in 'TLS':
                    self.log.subheading('Optimising {} parameters for TLS mode {}'.format(cpt, i_tls+1))
                    # Set penalty weights
                    self.log.bar()
                    self.log('Updating penalty weights --->')
                    self.log.bar()
                    self.log(self.penalty.set_weights(ovr_fixed=self.penalties.over_target_values.a))
                    self.log.bar(False, True)
                    # Select variables for optimisation -- model only
                    self._select(optimise_model      = True,
                                 optimise_amplitudes = False,
                                 components = cpt,
                                 models     = [i_tls],
                                 datasets   = opt_dset_mask,
                                 atoms      = opt_atom_mask)
                    # Run optimisation
                    self._optimise(running_summary=True)
                    # Log model summary
                    self.log(self._parameters.get(index=i_tls).model.summary())
                    # Log current RMSD and penalty
                    self.log.bar()
                    self.log('> Optimisation RMSD: {:5.3f}'.format(self.optimisation_rmsd))
                    self.log('> Optimisation Pen.: {:5.3f}'.format(self.optimisation_penalty))
                    self.log.bar()

                # ---------------------------------->
                # Optimise TLS amplitude parameters (all amplitudes!)
                self.log.subheading('Optimising TLS amplitudes for all datasets')
                self.log.bar()
                self.log('Updating penalty weights --->')
                self.log.bar()
                self.log(self.penalty.set_weights(ovr_fixed=10.0))
                self.log.bar()
                # If it's the first iteration, assert all amplitudes are zero (in case amplitudes are optimised for zero-value TLS models)
                if i_cycle == 0:
                    self._parameters.reset_amplitudes(components='TLS', modes=[i_tls])
                # Reset amplitudes to zero to prevent overparameterisation
                self.log('Resetting all amplitudes to zero')
                self._parameters.get(index=i_tls).amplitudes.set(vals=0.0, components='TLS')
                # Optimise amplitudes
                self.log('Setting up T-L-S amplitude optimisation')
                # Select variables for optimisation -- amplitudes only
                self._select(optimise_model      = False,
                             optimise_amplitudes = True,
                             components = 'TLS',
                             models     = mdl_cuml,
                             datasets   = None,
                             # Optimise with all atoms to prevent over-sizing
                             atoms      = range(self._n_atm))
                # Optimise all amplitudes dataset-by-dataset
                jobs = []
                for i_dst in xrange(self._n_dst):
                    # Select this dataset
                    self._select(datasets=[i_dst])
                    # Append to job list (create copy) or run optimisation
                    if n_cpus > 1:
                        jobs.append(self.copy())
                    else:
                        self._optimise(running_summary=False)
                        self.log('> dataset {} of {} (rmsd {:.3f}; penalty {:.1f})'.format(i_dst+1, self._n_dst, self.optimisation_rmsd, self.optimisation_penalty))
                # Run parallel jobs and inject results
                if n_cpus > 1:
                    self.log.subheading('Running {} jobs using {} cpus'.format(len(jobs), min(n_cpus,len(jobs))))
                    finished_jobs = workers.map(func=wrapper_optimise, iterable=jobs)
                    for i_dst in xrange(self._n_dst):
                        job = finished_jobs[i_dst]
                        if isinstance(job, str):
                            self.log(job)
                            raise Failure('error returned')
                        self._inject(values=job._values(), selection=job._select())
                        self.log('> dataset {} of {} (rmsd {:.3f}; penalty {:.1f})'.format(i_dst+1, self._n_dst, job.optimisation_rmsd,  job.optimisation_penalty))
                # Log model summary
                self.log.bar(blank_before=True, blank_after=True)
                self.log(self._parameters.get(index=i_tls).amplitudes.summary())

            # Reapply atom and dataset masks
            self.set_dataset_mask(opt_dset_mask)
            self.set_atomic_mask(opt_atom_mask)

            # End of cycle house-keeping
            self.log('')
            self.log('Looking for negative TLS amplitudes')
            self._parameters.reset_negative_amplitudes()
            self.log('')
            self.log('Normalising TLS amplitudes')
            self._parameters.normalise_amplitudes()
            self.log('')
            self.log('Looking for zero-value modes')
            self._parameters.reset_zero_value_modes()

            self.log.subheading('End-of-cycle summary')
            self.summary(show=True)

        # Let's be good
        if n_cpus > 1:
            workers.close()

    def extract(self):
        """Extract the fitted uij"""
        return self._extract(mask_datasets=None, mask_atoms=None)

    def result(self):
        """Extract the fitted parameters"""
        tls_mdls = numpy.array([p.model.values      for p in self._parameters])
        tls_amps = numpy.array([p.amplitudes.values for p in self._parameters])
        return (tls_mdls,tls_amps)

    def summary(self, show=True):
        """Print the number of parameters/input data"""
        s = self.log._bar()+'\nTLS Group Fit Summary: {}\n'.format(self.label)+self.log._bar()
        if self._optimise_model and (self.optimisation_rmsd is not numpy.inf):
            s += '\n> Optimisation Summary'
            s += '\nOptimisation RMSD:    {:5.3f}'.format(self.optimisation_rmsd)
            s += '\nOptimisation Penalty: {:5.3f}'.format(self.optimisation_penalty)
        for i_tls in xrange(self._n_tls):
            s += '\n> TLS model {}'.format(i_tls+1)
            mode = self._parameters.get(index=i_tls)
            s += '\n\t' + mode.model.summary().replace('\n','\n\t')
            s += '\n\t' + mode.amplitudes.summary().replace('\n','\n\t')
        if show: self.log(s)
        return s

############################################################################

def build_levels(model, params, log=None):
    """Build the levels for the hierarchical fitting"""

    log.subheading('Building hierarchy selections')
    levels = []; labels=[];

    # FIXME Only run on the protein for the moment FIXME
    filter_h = protein(model.hierarchy)

    if 'chain' in params.levels.auto_levels:
        log('Level {}: Creating level with groups for each chain'.format(len(levels)+1))
        levels.append([PhenixSelection.format(c) for c in filter_h.chains()])
        labels.append('chain')
    if 'auto_group' in params.levels.auto_levels:
        log('Level {}: Creating level with groups determined by phenix.find_tls_groups'.format(len(levels)+1))
        levels.append([s.strip('"') for s in phenix_find_tls_groups(model.filename)])
        labels.append('groups')
    if 'secondary_structure' in params.levels.auto_levels:
        log('Level {}: Creating level with groups based on secondary structure'.format(len(levels)+1))
        levels.append([s.strip('"') for s in default_secondary_structure_selections(model.hierarchy)])
        labels.append('secondary structure')
    if 'residue' in params.levels.auto_levels:
        log('Level {}: Creating level with groups for each residue'.format(len(levels)+1))
        levels.append([PhenixSelection.format(r) for r in filter_h.residue_groups()])
        labels.append('residue')
    if 'backbone' in params.levels.auto_levels:
        log('Level {}: Creating level with groups for each residue backbone'.format(len(levels)+1))
        levels.append([PhenixSelection.format(r)+' and (name C or name CA or name N or name O)'     for r in backbone(filter_h).atom_groups() if (r.resname not in ['ALA','GLY','PRO'])])
        labels.append('backbone')
    if 'sidechain' in params.levels.auto_levels:
        log('Level {}: Creating level with groups for each residue sidechain'.format(len(levels)+1))
        levels.append([PhenixSelection.format(r)+' and not (name C or name CA or name N or name O)' for r in sidechains(filter_h).atom_groups() if (r.resname not in ['ALA','GLY','PRO'])])
        labels.append('sidechain')
    if 'atom' in params.levels.auto_levels:
        log('Level {}: Creating level with groups for each atom'.format(len(levels)+1))
        levels.append([PhenixSelection.format(a) for a in filter_h.atoms()])
        labels.append('atom')

    # TODO IF statement for if custom levels are defined TODO
    # TODO allow ability to insert
    # TODO Take any custom groups and insert them here TODO (levels.insert)

    # Report
    log('\n> {} levels created\n'.format(len(levels)))
    for i_l, level in enumerate(levels):
        log.bar()
        log('Level {}'.format(i_l+1))
        log.bar()
        for l in level: log('\t'+l)

    return levels, labels

def run(params):

    easy_directory(params.output.out_dir)

    assert params.table_ones_options.column_labels
    assert params.table_ones_options.r_free_label

    log_dir = easy_directory(os.path.join(params.output.out_dir, 'logs'))
    log = Log(os.path.join(log_dir, '_fitting.log'), verbose=params.settings.verbose)

    # Report parameters
    log.heading('Input command')
    log(' \\\n\t'.join(['pandemic.adp'] + sys.argv[1:]))
    log.heading('Processed parameters')
    log(master_phil.format(params).as_str())

    # Load input pickle object
    if (params.input.pickle is None):
        log.heading('Running setup')

        # Load input structures
        log.subheading('Building model list'.format(len(params.input.pdb)))
        if params.input.labelling == 'foldername':
            label_func = lambda f: os.path.basename(os.path.dirname(f))
        else:
            label_func = lambda f: os.path.basename(os.path.splitext(f)[0])
        models = [CrystallographicModel.from_file(f).label(tag=label_func(f)) for f in params.input.pdb]#[:10]
        models = sorted(models, key=lambda m: m.tag)
        log('{} models loaded'.format(len(models)))

        # Construct the groups of atoms for each level
        levels, labels = build_levels(model=models[0], params=params, log=log)

        # Run the program main
        log.heading('Beginning parameterisation')
        p = MultiDatasetUijParameterisation(models = models,
                                            levels = levels,
                                            level_labels = labels,
                                            params = params,
                                            log = log)

        if params.settings.dry_run:
            log.heading('Exiting after initialisation: dry_run=True')
            sys.exit()

        # Fit TLS models
        p.fit_hierarchical_uij_model()

        # Process output
        p.process_results()

        # Pickle output object
        if p.params.output.pickle is True:
            log.subheading('Pickling output')
            p.pickle(pickle_file    = os.path.join(p.params.output.out_dir, 'pandemic.pickle'),
                     pickle_object  = p,
                     overwrite      = True)

    else:
        # Load previously-pickled parameterisation
        assert os.path.exists(params.input.pickle), 'Input pickle file does not exist: {}'.format(params.input.pickle)
        log.heading('Loading previous ADP parameterisation')
        p = libtbx.easy_pickle.load(params.input.pickle)
        assert isinstance(p, MultiDatasetUijParameterisation), 'Input pickle contains type: {}'.format(type(p))

        # Update unpickled object
        p.params = params
        #p.out_dir = os.path.abspath(params.output.out_dir)

    # Print all errors
    p.show_errors()

    # Write HTML
    p.log.heading('Writing HTML output')
    pandemic_html.write_adp_summary(parameterisation=p, out_dir=params.output.out_dir)


############################################################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
        program             = PROGRAM,
        description         = DESCRIPTION)



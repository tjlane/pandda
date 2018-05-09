import os, sys, glob, copy, traceback
import math, itertools, operator, collections

import scipy.cluster
import numpy, pandas

import libtbx.phil, libtbx.easy_mp, libtbx.easy_pickle
import iotbx.pdb
import mmtbx.tls.tools

import multiprocessing, multiprocessing.pool

from libtbx.utils import Sorry, Failure
from scitbx.array_family import flex
from scitbx import simplex, linalg, matrix

from bamboo.common import Meta, Info, ListStream
from bamboo.common.logs import Log, LogStream
from bamboo.common.path import easy_directory, rel_symlink
from bamboo.common.command import CommandManager

from giant.manager import Program
from giant.dataset import CrystallographicModel
from giant.structure.tls import tls_str_to_n_params, get_t_l_s_from_vector
from giant.structure.uij import uij_positive_are_semi_definite, \
        uij_eigenvalues, uij_to_b, calculate_uij_anisotropy_ratio, \
        scale_uij_to_target_by_selection
from giant.structure.formatting import ShortLabeller, PhenixSelection, PymolSelection
from giant.structure.select import protein, backbone, sidechains, default_secondary_structure_selections_filled
from giant.structure.pymol import auto_residue_images, auto_chain_images, selection_images
from giant.xray.crystal import CrystalSummary
from giant.xray.refine import refine_phenix
from giant.xray.tls import phenix_find_tls_groups

from giant.jiffies import multi_table_ones

from pandemic.adp import html as pandemic_html
from pandemic.adp.tls import TLSModel, TLS_AmplitudeSet, MultiDatasetTLSModelList, MultiDatasetTLSModel

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

numpy.set_printoptions(linewidth=numpy.inf, threshold=numpy.nan)

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
    Fit a consensus anisotropic B-factor model to a (series of) molecular structure(s)
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
        .type = choice(multi=False)
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
    html = True
        .help = "Write all output data in an html file"
        .type = bool
    images {
        all = False
            .help = "Generate all graphical output - may significantly increase runtime"
            .type = bool
        pymol = none *chain all
            .help = "Write residue-by-residue images of the output B-factors"
            .type = choice(multi=False)
        distributions = False
            .help = "Write distribution graphs for each TLS group"
            .type = bool
    }
    diagnostics = False
        .help = "Write diagnostic graphs -- adds a significant amount of runtime"
        .type = bool
    clean_up_files = *compress_logs *delete_mtzs
        .help = "Delete unnecessary output files (MTZs)"
        .type = choice(multi=True)
}
levels {
    overall_selection = None
        .type = str
    auto_levels = *chain auto_group *secondary_structure *residue *backbone *sidechain atom
        .type = choice(multi=True)
    custom_level
        .multiple = True
    {
        depth = None
            .help = "Where to insert this level into the hierarchy? (after auto_levels have been determined). Inserting as '1' will make this the first level and '2' will add as the second level, etc. Any auto-generated levels will be shifted down appropriately."
            .type = int
            .multiple = False
        selection = None
            .help = "list of selections that define groups for this level"
            .type = str
            .multiple = True
    }
}
fitting {
    tls_models_per_tls_group = 1
        .help = 'how many TLS models to fit per group of atoms?'
        .type = int
    tls_amplitude_model = *TLS
        .help = 'which amplitudes should be used for multiplying each TLS group in each dataset? \n\tTLS = one amplitude for each TLS group.'
        .type = choice(multi=False)
    number_of_macro_cycles = 2
        .help = 'how many fitting cycles to run (over all levels) -- should be more than 0'
        .type = int
    number_of_micro_cycles = 2
        .help = 'how many fitting cycles to run (for each level) -- should be more than 1'
        .type = int
    optimisation {
        dataset_weights = none inverse_resolution inverse_resolution_squared *inverse_resolution_cubed
            .help = 'control how datasets are weighted during optimisation?'
            .type = choice(multi=False)
        max_datasets = None
            .help = 'takes up to this number of datasets for TLS parameter optimisation'
            .type = int
        max_resolution = None
            .help = 'resolution limit for dataset to be used for TLS optimisation'
            .type = float
        step_size
            .help = 'set the various step sizes taken during simplex optimisation'
        {
            vibration = 0.1
                .help = 'RMS Vibration magnitude step size (angstroms)'
                .type = float
            libration = 0.1
                .help = 'RMS Libration magnitude step size (degrees)'
                .type = float
            angle = 1.0
                .help = 'Vector orientation change step size (degrees)'
                .type = float
            amplitude = 0.1
                .help = 'Amplitude magnitude change step size (dimensionless)'
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
                b = 0
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
    precision
        .help = 'set various levels of precision for calculations'
    {
        tls_tolerance = 1e-12
            .help = "tolerance for validating TLS matrices (cutoff for defining zero when calculating negative eigenvalues, etc)"
            .type = float
        uij_tolerance = 1e-1
            .help = "tolerance for validating Uij values (maximum allowed negative eigenvalues of Uijs)"
            .type = float
        tls_model_decimals = 6
            .help = "how many decimals of precision to be used for TLS matrices (PDB maximum is three, but more needed here due to small Uij-values)"
            .type = int
        tls_amplitude_decimals = 9
            .help = "how many decimals of precision to be used for TLS amplitudes"
            .type = int
    }
}
analysis {
    refine_output_structures = True
        .help = "Refine the structures after fitting (coordinates and occupancies)"
        .type = bool
    calculate_r_factors = True
        .help = "Recalculate r-factors for the fitted B-factors"
        .type = bool
}
refinement {
    cif = None
        .help = "Cif files required for refinement"
        .type = str
        .multiple = True
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
#                             Utility Functions                            #
############################################################################

def rms(vals, axis=None):
    return numpy.sqrt(numpy.mean(numpy.power(vals,2), axis=axis))

############################################################################
#                        Multi-processing functions                        #
############################################################################

def wrapper_plot_histograms(args):
    MultiDatasetUijPlots.multi_histogram(**args)

def wrapper_run(arg):
    return arg.run()

def wrapper_optimise(args):
    fitter, tolerance = args
    try:
        fitter._optimise(tolerance=tolerance,
                         running_summary=False)
        return fitter
    except Exception as e:
        return traceback.format_exc()

def wrapper_fit(args):
    fitter, kw_args = args
    s = kw_args.pop('silent', True)
    try:
        msg = 'Error'
        fitter.optimise(**kw_args)
        msg = 'Finished'
        return fitter
    except Exception as e:
        tr = traceback.format_exc()
        fitter.log(tr)
        return tr
    finally:
        print '{} ({}). Log written to {}'.format(msg, fitter.label, fitter.log.log_file)
        fitter.log.write_to_log(clear_data=True)

############################################################################

class UijPenalties(object):
    _base_tol = 1e-6
    _uij_tolerance = _base_tol
    _tls_tolerance = _base_tol
    _amp_tolerance = _base_tol
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

    @classmethod
    def set_tls_tolerance(cls, tolerance):
        cls._tls_tolerance = tolerance

    def tls_tolerance(self):
        return self._tls_tolerance

    @classmethod
    def set_uij_tolerance(cls, tolerance):
        cls._uij_tolerance = tolerance

    @classmethod
    def set_amp_tolerance(cls, tolerance):
        cls._amp_tolerance = tolerance

    def set_test_xyz(self, xyz, com):
        self._tst_xyz = flex.vec3_double(xyz)
        self._tst_com = com

    def set_weights(self, **args):
        """Set penalties for parameters to be invalid"""
        for k in args:
            assert k in self.values.__dict__.keys()
            self.values.__setattr__(k, args[k])
        return self.values.summary()

    def tls_params(self, model):
        """Return penalty for having unphysical TLS models"""
        # Convert to float for penalty
        penalty = 0.0 if model.is_valid(tol=self.tls_tolerance()) else 1.0
        # Calculate and return total penalty
        return self._standard_penalty_function(penalty       = penalty,
                                               fixed_penalty = self.values.mdl_fixed,
                                               slope_penalty = self.values.mdl_slope)

    def amplitudes(self, values):
        """Return penalty for having amplitudes with negative values"""
        penalty = abs(numpy.sum(values[values<-1.0*self._amp_tolerance]))
        return self._standard_penalty_function(penalty       = penalty,
                                               fixed_penalty = self.values.amp_fixed,
                                               slope_penalty = self.values.amp_slope)

    def uij_valid(self, values):
        assert len(values) == 6
        penalty = self._standard_negative_eigenvalue_penalty(sym_mat3_vals=values, tolerance=self._uij_tolerance)
        return self._standard_penalty_function(penalty       = penalty,
                                               fixed_penalty = self.values.uij_fixed,
                                               slope_penalty = self.values.uij_slope)

    def uij_size(self, fitted, target):
        """Add penalty for having fitted B-factors greater than observed"""
        if (self.values.ovr_fixed == 0.0) and (self.values.ovr_slope == 0.0):
            return 0.0
        values = [a-b for a,b in zip(target,fitted)]
        penalty = self._standard_negative_eigenvalue_penalty(sym_mat3_vals=values, tolerance=self._base_tol)
        return self._standard_penalty_function(penalty       = penalty,
                                               fixed_penalty = self.values.ovr_fixed,
                                               slope_penalty = self.values.ovr_slope)

    def _standard_penalty_function(self, penalty, fixed_penalty, slope_penalty):
        """Simple step function for non-zero penalty + slope to provide minimisation gradient"""
        assert penalty > -1.0*self._base_tol, 'penalties cannot be negative!'
        return (penalty>self._base_tol)*(fixed_penalty + slope_penalty*penalty)

    def _standard_negative_eigenvalue_penalty(self, sym_mat3_vals, tolerance):
        eigenvalues = self._sym_mat3_eigenvalues(vals=sym_mat3_vals)
        return self._sym_mat3_penalty(eigenvalues=eigenvalues, tolerance=tolerance)

    def _sym_mat3_eigenvalues(self, vals):
        assert len(vals) == 6
        return linalg.eigensystem_real_symmetric(vals).values()

    def _sym_mat3_penalty(self, eigenvalues, tolerance=None):
        """Calculate penalty for eigenvalues"""
        return -1.0*flex.sum(eigenvalues.select(eigenvalues<(-1.0*tolerance)))

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

    def get_simplex(self, start, **kw_args):
        delta = self.get_deltas(**kw_args)
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

    _eigenvalue_change_basis = ((+1.,+1.,+1.,0.,0.,0.),
                                (+1.,-1.,+0.,0.,0.,0.),
                                (+0.,+1.,-1.,0.,0.,0.))
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
    _s_change_basis = ((1., 0., 0.,
                        0.,-1., 0.,
                        0., 0., 0.),
                       (0., 0., 0.,
                        0., 1., 0.,
                        0., 0.,-1.),
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

    def get_deltas(self, model, amplitudes, components, n_cpt, n_mdl, n_dst, tls_models, i_dst, i_mdl): #, **kw_args
        p_mdl = 0
        p_amp = 0
        if model is True:
            p_mdl += tls_str_to_n_params(components, include_szz=False)*n_mdl
        if amplitudes is True:
            p_amp += n_mdl*n_cpt*n_dst

        # Populate the delta array
        # Total number of parameters
        n_all = p_mdl+p_amp
        # Begin indexing at 0 -- naturally
        n_this = 0
        deltas = numpy.zeros((n_all, n_all))
        # Deal with models
        if model is True:
            assert components is not None, 'must provide components when model is True'
            assert set('TLS').intersection(components), 'components must contain T, L or S'
            for i_mdl in xrange(n_mdl):
                # Extract decomposition
                d = tls_models[i_mdl].decompose()
                # If decomposition is invalid (fixing mode), obviously can't use the decomposition
                if not d.is_valid():
                    d = None
                # Generate base for each component
                for c in 'TLS':
                    # Skip if not selected
                    if c not in components:
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
        if amplitudes is True:
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

    def get_deltas(self):
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
        self.warnings = []

        self.params = params

        self._n_cpu = params.settings.cpus
        self._n_opt = params.fitting.optimisation.max_datasets

        self._allow_isotropic = True

        self._opt_datasets_res_limit = params.fitting.optimisation.max_resolution
        self._opt_datasets_selection = []
        self.dataset_weights = []

        self.models = models
        self.levels = levels
        self.level_labels = level_labels if level_labels else range(1,len(levels)+1)
        self.fitter = None

        # Misc files
        if self.params.refinement.cif:
            self.cifs = self.params.refinement.cif
        else:
            self.cifs = None

        # Create plot object
        self.plot = MultiDatasetUijPlots

        # Validate and add output paths, etc.
        self._init_checks()
        self._init_file_system()
        self._init_input_models()
        self._init_level_groups()
        self._init_tables()
        self._init_fitter()

        #self.write_running_parameters_to_log(params=params)

    def _init_checks(self):
        dependencies = ['phenix.refine', 'phenix.table_one']
        if self.params.output.images.pymol != 'none':
            dependencies.append('pymol')
        self.check_programs_are_available(dependencies)

    def _init_file_system(self):
        """Prepare the file manager object"""

        self.out_dir = easy_directory(self.params.output.out_dir)
        fm = self.initialise_file_manager(rootdir=self.params.output.out_dir)
        # Top level folders
        for d in ['logs','model','structures','results','analysis']:
            fm.add_dir(dir_name=d, dir_tag=d)
        # Analysis sub-folders
        for d in ['table_ones']:
            fm.add_dir(dir_name=d, dir_tag=d, top_dir_tag='results')

    def _init_input_models(self):
        """Prepare the input models"""

        self.log.subheading('Processing input models')

        # Use the first hierarchy as the reference
        self.log('Using {} as reference structure'.format(self.models[0].tag))
        self.master_h = self.models[0].hierarchy.deep_copy()

        errors = []

        # Extract the TableOne columns as a set to compare with the columns in the MTZs
        table_one_cols = set(self.params.table_ones_options.column_labels.split(',')+[self.params.table_ones_options.r_free_label])

        # Create dataset weight function
        if self.params.fitting.optimisation.dataset_weights == 'none':
            dataset_weight = lambda r: 1.0
        elif self.params.fitting.optimisation.dataset_weights == 'inverse_resolution':
            dataset_weight = lambda r: r**(-1.0)
        elif self.params.fitting.optimisation.dataset_weights == 'inverse_resolution_squared':
            dataset_weight = lambda r: r**(-2.0)
        elif self.params.fitting.optimisation.dataset_weights == 'inverse_resolution_cubed':
            dataset_weight = lambda r: r**(-3.0)
        else:
            raise Sorry('Invalid parameter provided for fitting.optimisation.dataset_weights: {}'.format(self.params.fitting.optimisation.dataset_weights))

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

            # Calculate dataset weight
            high_res = m.input.get_r_rfree_sigma().high
            if (self.params.fitting.optimisation.dataset_weights != 'none') and (high_res is None):
                errors.append(Failure("PDB does not contain resolution information in the REMARK 3 records.\n" + \
                                      "This is normally present in structures that have come from refinement.\n" + \
                                      "Are these structures from refinement?"))
            else:
                self.dataset_weights.append(dataset_weight(high_res))

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
            self.log('Dataset resolutions:')
            for m in self.models:
                self.log('> {}: {}'.format(m.tag, m.input.get_r_rfree_sigma().high))
            raise Sorry('No datasets above resolution cutoff: {}'.format(self._opt_datasets_res_limit))

        # Limit the number of datasets for optimisation
        if self._n_opt is not None:
            self.log('Limiting list of datasets for TLS optimisation to {} datasets'.format(self._n_opt))
            self._opt_datasets_selection = self._opt_datasets_selection[:self._n_opt]
        # Print optimisation datasets
        self.log('Using {} datasets for TLS and residual parameterisation'.format(len(self._opt_datasets_selection)))
        for i_m in self._opt_datasets_selection:
            self.log('\t'+self.models[i_m].tag)
        # Print optimisation weights
        self.log('Optimisation weights for TLS and residual parameterisation')
        for i_m, m in enumerate(self.models):
            self.log('\t{} (model {}) -> {:+10.6f} (resolution {})'.format(m.tag, i_m, self.dataset_weights[i_m], m.input.get_r_rfree_sigma().high))

    def _init_level_groups(self):
        """Create the selection array that builds up the hierarchy of the fitting"""

        self.log.subheading('Processing input levels')

        # Extract the atoms for each tls group
        atom_cache = self.master_h.atom_selection_cache()
        # Array of which atoms are in which group at which level
        level_array = -1 * numpy.ones((len(self.levels), self.master_h.atoms().size()), dtype=int)
        # List of any selections that result in no atoms
        warnings = []
        # Each column is one level
        for i_level, selections in enumerate(self.levels):
            no_atoms = 0
            self.log('\n> Level {} ({})\n'.format(i_level+1, self.level_labels[i_level]))
            for i_group, sel_group in enumerate(selections):
                sel = numpy.array(atom_cache.selection(sel_group))
                self.log('\tgroup {:<5d} - {:50}: {:>5d} atoms'.format(i_group+1, sel_group, sum(sel)))
                if sum(sel) == 0:
                    warnings.append('Level {}, Group "{}": {} atoms'.format(i_level+1, sel_group, sum(sel)))
                    no_atoms += 1
                    continue
                level_array[i_level, sel] = i_group+1
            # Sanity check that the expected number of groups have been produced
            n_groups = len([v for v in numpy.unique(level_array[i_level]) if v>0])
            assert n_groups+no_atoms == len(selections)
        if warnings:
            msg = 'WARNING: One or more group selections do not select any atoms: \n\t{}'.format('\n\t'.join(warnings))
            self.log(msg)
            self.warnings.append(msg)

        # Apply the overall filter if provided
        if self.params.levels.overall_selection:
            ivt_selection = numpy.logical_not(atom_cache.selection(self.params.levels.overall_selection))
            level_array[:, ivt_selection] = -1
            # Renumber the selections from 1
            for i_level, level in enumerate(self.levels):
                old_grp_nums = [v for v in numpy.unique(level_array[i_level]) if v>0]
                for i_new_grp, old_grp in enumerate(old_grp_nums):
                    group_sel = (level_array[i_level] == old_grp)
                    level_array[i_level, group_sel] = i_new_grp+1

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
                                                         'b_min (level)',
                                                         'b_max (level)',
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
                                                        dataset_weights = numpy.array(self.dataset_weights),
                                                        params = self.params.fitting,
                                                        log = self.log)
        # Set whether input is anisotropic or isotropic
        self.fitter.set_input_info(disorder_model=self.disorder_model)
        # Select the datasets for be used for TLS optimisation
        self.fitter.set_optimisation_datasets(self._opt_datasets_selection)

        # Table to be populated during optimisation
        self.fitter.set_tracking(table    = self.tables.tracking,
                                 csv_path = self.file_manager.add_file(file_name='tracking_data.csv', file_tag='tracking_csv', dir_tag='results'),
                                 png_path = self.file_manager.add_file(file_name='tracking_data.png', file_tag='tracking_png', dir_tag='results'))

        # Write summary of the fitted model (groups & levels)
        self.log.heading('Writing summary of the hierarchical model')
        self.file_manager.add_dir(dir_name='partitions', dir_tag='partitions', top_dir_tag='model')
        self.hierarchy_summary(out_dir_tag='partitions')

        # Calculate Parameter-observed data ratio
        self.log.subheading('Data-Parameter Ratios')
        n_o = numpy.product(self.fitter.observed_uij.shape)
        n_p = self.fitter.n_params(non_zero=False)
        self.log('Number of model parameters: {}'.format(n_p))
        self.log('Number of input uij values: {}'.format(n_o))
        self.log('Parameter ratio Gain (should be <1): {}'.format(float(n_p)/float(n_o)))

    def show_warnings(self):
        """Print non-fatal errors accumulated during running"""
        if len(self.warnings) == 0:
            self.log('Program completed with 0 errors')
            return
        self.log.subheading('{} non-fatal errors occurred during the program (see below)'.format(len(self.warnings)))
        for i, e in enumerate(self.warnings):
            self.log.bar()
            self.log('Error {} of {}'.format(i+1, len(self.warnings)))
            self.log.bar()
            self.log(e)
        self.log.bar()
        self.log.subheading('{} non-fatal errors occurred during the program (see above)'.format(len(self.warnings)))

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

    def extract_results(self):
        """Return the fitted uij for each level (+residual)"""
        uij_lvl = self.fitter.extract_tls(sum_levels=False)
        uij_res = self.fitter.residual.extract()
        return uij_lvl, uij_res

    def write_output_structures(self):
        """Extract output and generate summaries"""

        #------------------------------------------------------------------------------#
        #---#                            Extract uijs                              #---#
        #------------------------------------------------------------------------------#

        self.log('Extracting output Uijs from parameterised model...')

        # Extract the different components of the fitted uij values
        uij_lvl, uij_res = self.extract_results()
        # Calculate cumulative combinations
        uij_tls = uij_lvl.sum(axis=0)
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
        # TLS-parameterised structures (total TLS)
        self.log.bar()
        self.log('Writing fitted structures (total TLS components - all levels)')
        self.log.bar()
        self.output_structures(uij = uij_tls,
                               iso = map(uij_to_b, uij_tls),
                               headers = tls_headers,
                               out_dir = structure_dir,
                               model_suffix = '.tls-all-levels.pdb')
        # Level by level TLS-parameterised structures (single level contribution)
        for i_level in xrange(len(self.levels)):
            self.log.bar()
            self.log('Writing fitted structures (individual TLS components - level {})'.format(i_level+1))
            self.log.bar()
            self.output_structures(uij = uij_lvl[i_level],
                                   iso = map(uij_to_b, uij_lvl[i_level]),
                                   headers = self.make_tls_headers(datasets=None, levels=[i_level]),
                                   out_dir = structure_dir,
                                   model_suffix = '.tls-level-{:04}.pdb'.format(i_level+1))
        # Level by level TLS-parameterised structures (cumulative level contributions)
        for i_level in xrange(1, len(self.levels)):
            self.log.bar()
            self.log('Writing fitted structures (cumulative TLS components - level {})'.format(i_level+1))
            self.log.bar()
            cuml_uij = uij_lvl[0:i_level+1].sum(axis=0)
            self.output_structures(uij = cuml_uij,
                                   iso = map(uij_to_b, cuml_uij),
                                   headers = self.make_tls_headers(datasets=None, levels=range(i_level+1)),
                                   out_dir = structure_dir,
                                   model_suffix = '.tls-level-0001-to-{:04}.pdb'.format(i_level+1))
        # TLS-subtracted structures ("residual structures")
        self.log.bar()
        self.log('Writing fitted structures (total TLS subtracted from input Uijs)')
        self.log.bar()
        self.output_structures(uij = uij_inp-uij_tls,
                               iso = map(uij_to_b, uij_inp-uij_tls),
                               out_dir = structure_dir,
                               model_suffix = '.tls_subtracted.pdb')

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
            except Exception as e:
                self.log.bar()
                self.log('ERROR: Failed to write structure for dataset {}: ({})'.format(mdl.tag, mdl_f))
                self.log(traceback.format_exc())
                self.log.bar()
                continue

        return pdbs

    def hierarchy_summary(self, out_dir_tag):
        """Write out the composition of the hierarchical model"""

        # Extract the global mask and convert to flex
        atom_sel = flex.bool(self.atom_mask.tolist())

        # Generate hierarchy for each level with groups as b-factors
        self.log('Writing partition figures for each chain')
        level_hierarchies = []
        f_template = self.file_manager.add_file(file_name = 'level-{:04d}_atoms.pdb',
                                                file_tag  = 'level-partition-template',
                                                dir_tag   = out_dir_tag)
        for i_level, g_vals in enumerate(self.level_array):
            h = self.custom_master_hierarchy(uij=None, iso=g_vals, mask=atom_sel)
            # Write the structure
            filename = f_template.format(i_level+1)
            self.log('\t> {}'.format(filename))
            h.write_pdb_file(filename)
            # Append to list for plotting
            level_hierarchies.append(h)
        # Write hierarchy plot for each chain
        b_h = self.blank_master_hierarchy()
        b_c = b_h.atom_selection_cache()
        f_template = self.file_manager.add_file(file_name = 'level-partitions-chain-{}.png',
                                                file_tag  = 'level-partitions-by-chain-template',
                                                dir_tag   = out_dir_tag)
        # Generate png of partitioning
        for c in b_h.chains():
            chain_sel = b_c.selection('chain {}'.format(c.id))
            hierarchies = [h.select(chain_sel, copy_atoms=True) for h in level_hierarchies]
            # Skip if no partitions in this chain
            if not (numpy.array([h.atoms().extract_b() for h in hierarchies]).sum() > 1e-3):
                continue
            filename = f_template.format(c.id)
            self.log('\t> {}'.format(filename))
            self.plot.level_plots(filename=filename, hierarchies=hierarchies, title='chain {}'.format(c.id))
        # Generate images of each chain of each level coloured by group
        lvl_pml_prt_template = self.file_manager.add_file(file_name='partition-level_{}-chain_{}.png',
                                                          file_tag='pml-level-partition-template',
                                                          dir_tag=out_dir_tag)
        lvl_pml_prt_prefix = lvl_pml_prt_template.replace('-chain_{}.png','')
        if self.params.output.images.pymol != 'none':
            for i_level in xrange(len(self.level_array)):
                # Images for each chain (of the partitions) - coloured by group
                self.log('\t> {}'.format(lvl_pml_prt_template.format(i_level+1,'*')))
                s_f = self.file_manager.get_file('level-partition-template').format(i_level+1)
                # Choose the style based on whether interested in atoms or larger regions
                styles = 'cartoon' if self.level_labels[i_level] in ['chain','groups','secondary structure','residue'] else 'lines+spheres'
                auto_chain_images(structure_filename = s_f,
                                  output_prefix = lvl_pml_prt_prefix.format(i_level+1),
                                  style = styles,
                                  colour_by = 'bfactor',
                                  settings = [('cartoon_oval_length', 0.5),
                                              ('cartoon_discrete_colors', 'on'),
                                              ('sphere_scale', 0.25)],
                                  width=1000, height=750)
                assert glob.glob(lvl_pml_prt_template.format(i_level+1,'*')), 'no files have been generated!'

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
        if (datasets is None):
            tls_origins = self.fitter.tls_origins
        else:
            tls_origins = self.fitter.tls_origins[:,datasets,:]    # (n_partitions, n_datasets, 3)
        tls_origins_hash = self.fitter.tls_origins_hash

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

            # Get the origin for this dataset
            group_origins = [f for f in fitters if (f is not None)][0].atomic_com

            # Iterate through the datasets and create TLSOs for each dataset
            for i_dst, tls_list in enumerate(tls_models):
                # Add all the levels for each dataset
                total_model = reduce(operator.add, tls_list)
                # Extract TLS remark lines from the
                tls_obj = mmtbx.tls.tools.tlso(t = total_model.get('T'),
                                               l = total_model.get('L'),
                                               s = total_model.get('S'),
                                               origin = group_origins[i_dst])
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

    def write_tls_level_summary(self, out_dir_tag):
        """Write the various TLS uijs to the master hierarchy structure"""

        fm = self.file_manager
        out_dir = fm.get_dir(out_dir_tag)
        csv_dir = fm.add_dir(dir_name='csvs',   dir_tag='tls-csvs',   top_dir_tag=out_dir_tag)
        pdb_dir = fm.add_dir(dir_name='pdbs',   dir_tag='tls-pdbs',   top_dir_tag=out_dir_tag)
        png_dir = fm.add_dir(dir_name='graphs', dir_tag='tls-graphs', top_dir_tag=out_dir_tag)
        pml_dir = fm.add_dir(dir_name='pymol',  dir_tag='tls-pymol',  top_dir_tag=out_dir_tag)

        # ------------------------
        self.log.subheading('Writing TLS models and amplitudes for each level')
        # ------------------------
        # Output file names -- standard templates to be used for many things...
        # ------------------------
        mdl_template = fm.add_file(file_name='tls_models_level_{:04d}.csv',                file_tag='csv-tls-mdl-template', dir_tag='tls-csvs')
        amp_template = fm.add_file(file_name='tls_amplitudes_level_{:04d}.csv',            file_tag='csv-tls-amp-template', dir_tag='tls-csvs')
        dst_template = fm.add_file(file_name='tls-model-amplitudes-level-{}-group-{}.png', file_tag='png-tls-amp-dist-template', dir_tag='tls-graphs')
        # ------------------------
        # Iterate through the levels
        for level in self.fitter.levels:
            self.log('Level {}'.format(level.index))
            # Table for TLS model components
            mdl_filename = mdl_template.format(level.index)
            mdl_table = pandas.DataFrame(columns=["group", "model",
                                                  "T11","T22","T33","T12","T13","T23",
                                                  "L11","L22","L33","L12","L13","L23",
                                                  "S11","S12","S13","S21","S22","S23","S31","S32","S33"])
            # Create amplitude table
            amp_filename = amp_template.format(level.index)
            amp_table = pandas.DataFrame(columns=["group", "model", "cpt"]+[mdl.tag for mdl in self.models])
            # Iterate through the groups in this level
            for i_group, sel, fitter in level:
                tls_model, tls_amps = fitter.result()
                assert tls_model.shape == (self.params.fitting.tls_models_per_tls_group, 21)
                assert tls_amps.shape  == (self.params.fitting.tls_models_per_tls_group, len(self.models), self.params.fitting.tls_amplitude_model.count('/')+1)

                # Extract parameters
                p = fitter.parameters()

                # Add to model and amplitudes tables
                for i_tls in xrange(tls_model.shape[0]):
                    # Add model values to last row of table
                    mdl_table.loc[len(mdl_table.index)] = numpy.concatenate([[i_group, i_tls], tls_model[i_tls]])
                    # Add amplitudes to last row of table
                    for i_cpt, cpt in enumerate(p.get(index=i_tls).component_sets()):
                        amp_table.loc[len(amp_table.index)] = numpy.concatenate([[i_group, i_tls, cpt], tls_amps[i_tls,:,i_cpt]])

                # Write histograms of amplitudes -- only for non-zero models
                if (tls_model.sum() > 0.0) and self.params.output.images.distributions:
                    filename = dst_template.format(level.index, i_group)
                    self.log('\t> {}'.format(filename))
                    titles = numpy.concatenate([['Model {}: {}'.format(i_t+1, c) for c in p.get(index=i_t).component_sets()] for i_t in xrange(tls_amps.shape[0])])
                    x_vals = numpy.concatenate([[tls_amps[i_t,:,i_c]             for i_c in xrange(tls_amps.shape[2])]       for i_t in xrange(tls_amps.shape[0])])
                    self.plot.multi_histogram(filename  = filename,
                                              x_vals    = x_vals,
                                              titles    = titles,
                                              x_labs    = ['']*tls_amps.shape[0]*tls_amps.shape[2],
                                              rotate_x_labels = True,
                                              shape     = (tls_amps.shape[0], tls_amps.shape[2]),
                                              n_bins    = 30, x_lim=[0, None])
            # Write tables
            self.log('\t> {}'.format(mdl_filename))
            mdl_table.to_csv(mdl_filename)
            self.log('\t> {}'.format(amp_filename))
            amp_table.to_csv(amp_filename)

        # ------------------------
        self.log.subheading('Writing contribution of each TLS model for each level')
        # ------------------------
        # Output file names
        # ------------------------
        sgl_pdb_template = fm.add_file(file_name='level_{}-mode_{}.pdb',
                                       file_tag='pdb-tls-one-mode-template',
                                       dir_tag='tls-pdbs')
        lvl_pdb_template = fm.add_file(file_name='level_{}-all.pdb',
                                       file_tag='pdb-tls-one-level-template',
                                       dir_tag='tls-pdbs')
        # Stacked images
        lvl_plt_template = fm.add_file(file_name='uij-profile-level_{}-chain_{}.png',
                                       file_tag='png-tls-profile-template',
                                       dir_tag=out_dir_tag)
        lvl_plt_prefix = lvl_plt_template.replace('-chain_{}.png','')
        # Uij anisotropy
        lvl_aniso_plt_template = fm.add_file(file_name='uij-anisotropy-level_{}-chain_{}.png',
                                             file_tag='png-tls-anisotropy-template',
                                             dir_tag=out_dir_tag)
        lvl_aniso_plt_prefix = lvl_aniso_plt_template.replace('-chain_{}.png','')
        # Pymol templates -- create even if not used so that html can populate with dummy files
        lvl_pml_chn_template = fm.add_file(file_name='level_{}-chain_{}.png',  file_tag='pml-level-chain-template',  dir_tag='tls-pymol')
        lvl_pml_grp_template = fm.add_file(file_name='level_{}-group_{}.png',  file_tag='pml-level-group-template',  dir_tag='tls-pymol')
        lvl_pml_scl_template = fm.add_file(file_name='level_{}-scaled_{}.png', file_tag='pml-level-scaled-template', dir_tag='tls-pymol')
        lvl_pml_chn_prefix = lvl_pml_chn_template.replace('-chain_{}.png','')
        lvl_pml_grp_prefix = lvl_pml_grp_template.replace('{}.png','')
        lvl_pml_scl_prefix = lvl_pml_scl_template.replace('{}.png','')
        # ------------------------
        # Extract the atom mask to apply b-factors
        atom_sel = flex.bool(self.atom_mask.tolist())
        # Iterate through the levels and plot TLS contributions for each mode for each model for each level
        for i_level, level in enumerate(self.fitter.levels):
            self.log('Level {}'.format(level.index))
            # Boundaries for this level
            boundaries = self.partition_boundaries_for_level(i_level=i_level).select(atom_sel)
            # Cumulative uijs
            uij_all = []
            # Extract Uijs for individual modes/models
            for i_tls in xrange(self.params.fitting.tls_models_per_tls_group):
                self.log('- Extracting TLS contributions of model {} of level {}'.format(i_tls+1, i_level+1))
                # Create copy of the level for resetting T-L-S components
                l_copy = copy.deepcopy(level)
                # Reset the TLS for all other models
                if self.params.fitting.tls_models_per_tls_group > 1:
                    other_models = range(self.params.fitting.tls_models_per_tls_group)
                    other_models.remove(i_tls)
                    l_copy.zero_amplitudes(models=other_models)
                # Extract uijs
                uij_this = l_copy.extract(average=True)
                # Append to total list
                uij_all.append((i_tls, uij_this))
                # Output structure for MODEL - XXX set to True for developing XXX
                if True or (self.params.fitting.tls_models_per_tls_group > 1):
                    uij = numpy.sum([t[-1] for t in uij_all if (t[0]==i_tls)], axis=0)
                    m_h = self.custom_master_hierarchy(uij=uij, iso=uij_to_b(uij), mask=atom_sel)
                    m_f = sgl_pdb_template.format(i_level+1, i_tls+1)
                    self.log('\t> {}'.format(m_f))
                    m_h.write_pdb_file(m_f)

            # ------------------------
            # Write stacked bar plot of TLS for this level
            # ------------------------
            self.log('- Writing stacked bar plot of TLS contributions for this level')
            mdls, uijs = zip(*uij_all)
            self.log('\t> {}'.format(lvl_plt_template.format(i_level+1,'*')))
            self.plot.stacked_bar(prefix        = lvl_plt_prefix.format(i_level+1),
                                  hierarchies   = [self.custom_master_hierarchy(uij=u, iso=uij_to_b(u), mask=atom_sel).select(atom_sel) for u in uijs],
                                  legends       = ['TLS (Model {})'.format(i+1) for i in mdls],
                                  title         = 'Level {} - individual TLS model contributions'.format(i_level+1),
                                  v_line_hierarchy = boundaries if (sum(boundaries.atoms().extract_b()) < 0.5*len(list(boundaries.residue_groups()))) else None)
            assert glob.glob(lvl_plt_template.format(i_level+1,'*')), 'no files have been generated!'

            # ------------------------
            # Write out structure & graph of combined uijs
            # ------------------------
            # Output structure for LEVEL
            self.log('- Combining contributions for all TLS modes of level {}'.format(i_level+1))
            uij = numpy.sum([t[-1] for t in uij_all], axis=0)
            m_h = self.custom_master_hierarchy(uij=uij, iso=uij_to_b(uij), mask=atom_sel)
            m_f = lvl_pdb_template.format(i_level+1)
            self.log('\t> {}'.format(m_f))
            m_h.write_pdb_file(m_f)
            # Write pymol images of each chain
            if self.params.output.images.pymol != 'none':
                # Generate pymol images
                self.log('- Generating pymol images')
                # Images for each chain
                self.log('\t> {}'.format(lvl_pml_chn_template.format(i_level+1,'*')))
                auto_chain_images(structure_filename = m_f,
                                  output_prefix = lvl_pml_chn_prefix.format(i_level+1),
                                  style = 'lines+ellipsoids',
                                  width=1000, height=750)
                assert glob.glob(lvl_pml_chn_template.format(i_level+1,'*')), 'no files have been generated!'
            # Write pymol images of each group
            if self.params.output.images.pymol == 'all':
                # Output structure for LEVEL (scaled!) - for html only
                m_h_scl = scale_uij_to_target_by_selection(hierarchy=m_h,
                                                           selections=self.levels[i_level],
                                                           target=0.5)
                m_f_scl = m_f.replace('.pdb', '-scaled.pdb')
                self.log('\t> {}'.format(m_f_scl))
                m_h_scl.write_pdb_file(m_f_scl)
                # Images for each group
                blank_h = self.blank_master_hierarchy()
                cache_h = blank_h.atom_selection_cache()
                self.log('\t> {}'.format(lvl_pml_grp_template.format(i_level+1,'*')))
                selections = [PymolSelection.join_or([PymolSelection.format(a) for a in blank_h.atoms().select(cache_h.selection(s))]) for s in self.levels[i_level]]
                selection_images(structure_filename = m_f,
                                 output_prefix = lvl_pml_grp_prefix.format(i_level+1),
                                 selections = selections,
                                 style = 'lines+ellipsoids',
                                 width = 250 if level.n_groups()>50 else 1000,
                                 height= 250 if level.n_groups()>50 else 750)
                assert glob.glob(lvl_pml_grp_template.format(i_level+1,'*')), 'no files have been generated!'
                # Images for each group
                blank_h = self.blank_master_hierarchy()
                cache_h = blank_h.atom_selection_cache()
                self.log('\t> {}'.format(lvl_pml_scl_template.format(i_level+1,'*')))
                selections = [PymolSelection.join_or([PymolSelection.format(a) for a in blank_h.atoms().select(cache_h.selection(s))]) for s in self.levels[i_level]]
                selection_images(structure_filename = m_f_scl,
                                 output_prefix = lvl_pml_scl_prefix.format(i_level+1),
                                 selections = selections,
                                 style = 'lines+ellipsoids',
                                 width = 250 if level.n_groups()>50 else 1000,
                                 height= 250 if level.n_groups()>50 else 750)
                assert glob.glob(lvl_pml_scl_template.format(i_level+1,'*')), 'no files have been generated!'
            # ------------------------
            # Write out structure & graph of anisotropy
            # ------------------------
            self.log('- Calculating anistropy of Uijs')
            ani = 1.0 - calculate_uij_anisotropy_ratio(uij=uij)
            m_h = self.custom_master_hierarchy(uij=None, iso=ani, mask=atom_sel)
            # Make plot of the uij anisotropy over the chain
            self.log('\t> {}'.format(lvl_aniso_plt_template.format(i_level+1,'*')))
            self.plot.stacked_bar(prefix=lvl_aniso_plt_prefix.format(i_level+1),
                                  hierarchies=[m_h.select(atom_sel)],
                                  legends=['Anisotropy  '],
                                  title='Anisotropy of Level {}\n[fully isotropic = 0; fully anisotropic = 1]'.format(i_level+1),
                                  y_lim=(0.0,1.0),
                                  y_lab='Anisotropy of Uij ($1 - \\frac{E_{min}}{E_{max}}$)')
            assert glob.glob(lvl_aniso_plt_template.format(i_level+1,'*')), 'no files have been generated!'


        return

    def write_residuals_summary(self, out_dir_tag):
        """Write the residual uijs to the master hierarchy structure"""

        fm = self.file_manager
        out_dir = fm.get_dir(out_dir_tag)
        pdb_dir = fm.add_dir(dir_name='pdbs',   dir_tag='res-pdbs',   top_dir_tag=out_dir_tag)
        #png_dir = fm.add_dir(dir_name='graphs', dir_tag='res-graphs', top_dir_tag=out_dir_tag)
        pml_dir = fm.add_dir(dir_name='pymol',  dir_tag='res-pymol',  top_dir_tag=out_dir_tag)

        self.log('Writing structure and plots of residual ADPs')
        # ------------------------
        # Output file names
        # ------------------------
        res_pdb = fm.add_file(file_name='residual.pdb',
                              file_tag='pdb-residual',
                              dir_tag='tls-pdbs')
        # Uij profiles
        res_plt_template = fm.add_file(file_name='uij-profile-residual-chain_{}.png',
                                       file_tag='png-residual-profile-template',
                                       dir_tag=out_dir_tag)
        res_plt_prefix = res_plt_template.replace('-chain_{}.png','')
        # Uij anisotropy
        res_aniso_plt_template = fm.add_file(file_name='uij-anisotropy-residual-chain_{}.png',
                                             file_tag='png-residual-anisotropy-template',
                                             dir_tag=out_dir_tag)
        res_aniso_plt_prefix = res_aniso_plt_template.replace('-chain_{}.png','')
        # Pymol templates -- create even if not used so that html can populate with dummy files
        res_pml_chn_template = fm.add_file(file_name='residual-chain_{}.png',
                                           file_tag='pml-residual-chain-template',
                                           dir_tag='res-pymol')
        res_pml_grp_template = fm.add_file(file_name='residual-residue_{}.png',
                                           file_tag='pml-residual-group-template',
                                           dir_tag='res-pymol')
        res_pml_chn_prefix = res_pml_chn_template.replace('-chain_{}.png','')
        res_pml_grp_prefix = res_pml_grp_template.replace('-residue_{}.png','')
        # ------------------------
        # Extract residual uijs
        # ------------------------
        uij = self.fitter.residual.extract()
        # Extract the global mask and convert to flex
        atom_sel = flex.bool(self.atom_mask.tolist())
        # ------------------------
        # Write out structure & graph of uijs
        # ------------------------
        m_h = self.custom_master_hierarchy(uij=uij, iso=uij_to_b(uij), mask=atom_sel)
        m_f = res_pdb
        self.log('\t> {}'.format(m_f))
        m_h.write_pdb_file(m_f)
        # Make plot of the uij profile over each chain
        self.log('\t> {}'.format(res_plt_template.format('*')))
        self.plot.stacked_bar(prefix=res_plt_prefix,
                              hierarchies=[m_h.select(atom_sel)],
                              legends=['Residual    '],
                              title='Uij Profile of Residual Level',
                              colour_space=(1,1))
        assert glob.glob(res_plt_template.format('*')), 'no files have been generated!'
        # Write pymol images for each chain
        if self.params.output.images.pymol != 'none':
            self.log('- Generating pymol images')
            self.log('\t> {}'.format(res_pml_chn_template.format('*')))
            auto_chain_images(structure_filename = m_f,
                              output_prefix = res_pml_chn_prefix,
                              style = 'lines+ellipsoids',
                              width=1000, height=750)
            assert glob.glob(res_pml_chn_template.format('*')), 'no files have been generated!'
        # Write pymol images for each group
        if self.params.output.images.pymol == 'all':
            self.log('\t> {}'.format(res_pml_grp_template.format('*')))
            auto_residue_images(structure_filename = m_f,
                                output_prefix = res_pml_grp_prefix,
                                style = 'lines+ellipsoids',
                                width=250, height=250)
            assert glob.glob(res_pml_grp_template.format('*')), 'no files have been generated!'
        # ------------------------
        # Write out structure & graph of anisotropy
        # ------------------------
        self.log('- Calculating anistropy of Uijs')
        ani = 1.0 - calculate_uij_anisotropy_ratio(uij=uij)
        m_h = self.custom_master_hierarchy(uij=None, iso=ani, mask=atom_sel)
        # Make plot of the uij anisotropy over the chain
        self.log('\t> {}'.format(res_aniso_plt_template.format('*')))
        self.plot.stacked_bar(prefix=res_aniso_plt_prefix,
                              hierarchies=[m_h.select(atom_sel)],
                              legends=['Anisotropy  '],
                              title='Anisotropy of Residual Level\n[fully isotropic = 0; fully anisotropic = 1]',
                              y_lim=(0.0,1.0),
                              y_lab='Anisotropy of Uij ($1 - \\frac{E_{min}}{E_{max}}$)')
        assert glob.glob(res_aniso_plt_template.format('*')), 'no files have been generated!'

        # TODO Output CSV of all residual components

        return

    def write_combined_summary_graphs(self, out_dir_tag):
        """Write the distributions of each level TLS contributions over the structure"""

        fm = self.file_manager
        out_dir = fm.get_dir(out_dir_tag)
        stacked_template = fm.add_file(file_name='uij-profile-stacked-chain_{}.png', file_tag='png-combined-profile-template', dir_tag=out_dir_tag)
        stacked_prefix = stacked_template.replace('-chain_{}.png','')

        self.log('Writing profiles of TLS and residual Uijs over the model')

        # Extract the different components of the fitted uij values
        uij_lvl, uij_res = self.extract_results()
        # Extract the global mask and convert to flex
        atom_sel = flex.bool(self.atom_mask.tolist())
        # List of hierarchies for each level
        uij_hierarchies = []
        # One hierarchy for each level
        for i_level in range(len(uij_lvl)):
            uij_lvl_av = numpy.mean(uij_lvl[i_level], axis=0)
            h_lvl = self.custom_master_hierarchy(iso=uij_to_b(uij_lvl_av), mask=atom_sel)
            uij_hierarchies.append(h_lvl)
        # One hierarchy for the residual
        h_res = self.custom_master_hierarchy(iso=uij_to_b(uij_res), mask=atom_sel)
        uij_hierarchies.append(h_res)
        # Create stacked plot
        self.log('\t> {}'.format(stacked_template.format('*')))
        self.plot.stacked_bar(prefix=stacked_prefix,
                              hierarchies=[h.select(atom_sel) for h in uij_hierarchies],
                              legends=['Level {}     '.format(i+1) for i in range(len(uij_lvl))]+['Residual    '],
                              title='TLS and residual contributions',
                              reverse_legend_order=True)
        assert glob.glob(stacked_template.format('*')), 'no files have been generated!'

    def run_diagnostics(self):
        """Compare the fitted and input uijs for all structures"""

        fit_dir = self.file_manager.add_dir(dir_name='diagnostics', dir_tag='diagnostics')

        # Extract the different components of the fitted uij values
        uij_lvl, uij_res = self.extract_results()
        # Calculate cumulative combinations
        uij_tls = uij_lvl.sum(axis=0)
        uij_all = uij_tls + uij_res
        # Extract the input uij values
        uij_inp = self.fitter.observed_uij

        # Dataset-by-dataset and atom-by-atom fit rmsds
        self.log.heading('Calculating rmsd metrics between input and fitted Uijs')
        self.analyse_fit_rms_distributions(uij_fit=uij_all,
                                           uij_inp=uij_inp,
                                           out_dir=fit_dir)

        # Correlations between the residual atomic component and the
        self.log.heading('Calculating correlations between atomic residuals and input-fitted Uijs')
        self.analyse_fit_residual_correlations(uij_diff=(uij_inp-uij_all),
                                               out_dir=fit_dir)

    def analyse_rms_fit_distributions(self, uij_fit, uij_inp, out_dir='./', max_x_width=25):
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
        self.plot.multi_hierarchy_bar(prefix=prefix,
                                      hierarchies=[m_h],
                                      y_labs=['RMSD Mean (fitted uij - input uij) ($\AA$)'])

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
        self.plot.multi_hierarchy_bar(prefix=prefix,
                                      hierarchies=[m_h],
                                      y_labs=['RMSD IQR (fitted uij - input uij) ($\AA$)'])

    def analyse_residual_correlations(self, uij_diff, out_dir='./', max_x_width=25):
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

    def process_output_structures(self):
        """Refine the output models and calculate new R-free/R-work values"""

        #------------------------------------------------------------------------------#
        #---#              Refine the output models if requested                   #---#
        #------------------------------------------------------------------------------#

        if self.params.analysis.refine_output_structures:
            self.log.heading('Refining output structures')
            self.refine_fitted_dataset_models()

        #------------------------------------------------------------------------------#
        #---#             Analyse the quality of the output models                 #---#
        #------------------------------------------------------------------------------#

        if self.params.analysis.calculate_r_factors:
            self.log.heading('Generating Table Ones for all structures')
            self.generate_fitted_table_ones()

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
            if os.path.exists(obj.out_pdb_file):
                raise Failure('Refined PDB already exists! (model {})'.format(mdl.tag))
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

        if os.path.exists(table_one_csv_orig) or \
           os.path.exists(table_one_csv_fitd) or \
           os.path.exists(table_one_csv_refd):
            raise Failure('Output table ones already exist.')

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

        try:
            assert os.path.exists(table_one_csv_orig)
            assert os.path.exists(table_one_csv_fitd)
            if self.models[0].r_pdb is not None:
                assert os.path.exists(table_one_csv_refd)
        except Exception as e:
            raise Sorry('phenix.table_one did not generate output for one or more sets of structures.\n' + \
                        'Please look at the log files in this folder for information:\n' + \
                        os.path.dirname(os.path.realpath(table_one_csv_orig)))

    def write_combined_csv(self, out_dir_tag):
        """Add data to CSV and write"""

        fm = self.file_manager
        out_dir = fm.get_dir(out_dir_tag)
        out_csv = fm.add_file(file_name='output_data.csv', file_tag='output-csv', dir_tag=out_dir_tag)

        # Extract uij values
        uij_lvl, uij_res = self.extract_results()
        uij_fit = uij_lvl.sum(axis=0) + uij_res
        uij_inp = self.fitter.observed_uij

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
        for suff, csv in [(' (Input)',   self.file_manager.get_file('table_one_input')),
                          (' (Fitted)',  self.file_manager.get_file('table_one_fitted')),
                          (' (Refined)', self.file_manager.get_file('table_one_refined'))]:
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
        self.log('')
        self.log('Writing output csv: {}'.format(out_csv))
        main_table.to_csv(out_csv)

        # Store table for HTML output
        self.tables.statistics = main_table

        # Make graphs for the table
        self.show_table_one_summary(table=main_table)
        self.write_table_one_graphs(table=main_table,
                                    results_dir_tag='results',
                                    analysis_dir_tag='analysis')

    def show_table_one_summary(self, table):
        """Look at pre- and post-fitting statistics"""

        self.log.subheading('Model improvement summary')

        # Extract Column averages (means and medians)
        table_means = table.mean().round(3)
        shrt_str = '{:>35s} | {:>15} |'
        long_str = '{:>35s} | {:>15} | {:>15} | {:>15}'

        self.log.bar()
        self.log('Dataset Averages:')
        self.log.bar()
        # Columns without old/new prefix
        for lab in table_means.index:
            if lab.lower().endswith('-input)'): continue
            if lab.lower().endswith('(input)'): continue
            if lab.lower().endswith('(fitted)'): continue
            if lab.lower().endswith('(refined)'): continue
            self.log(shrt_str.format(lab, table_means[lab]))
        self.log.bar()
        # Columns with old/new prefix
        for comp_lab in ['Fitted','Refined']:
            comp_suff = ' ({})'.format(comp_lab)
            # Check whether the columns are here - check the R-free as always should be present
            if 'R-free ({})'.format(comp_lab) not in table.columns:
                continue
            self.log('Comparing Input and {} Uijs'.format(comp_lab))
            self.log(long_str.format('', 'Input', comp_lab, 'Difference'))
            for lab in table_means.index:
                if not lab.endswith(comp_suff):
                    continue
                lab_base = lab[:lab.index(comp_suff)]
                self.log(long_str.format(lab_base,
                                         '{:.3f}'.format(table_means[lab_base+' ({})'.format('Input')]),
                                         '{:.3f}'.format(table_means[lab_base+' ({})'.format(comp_lab)]),
                                         '{:.3f}'.format(table_means[lab_base+' Change ({}-{})'.format(comp_lab, 'Input')])))
            self.log.bar()

        return

    def write_table_one_graphs(self, table, results_dir_tag, analysis_dir_tag):
        """Look at pre- and post-fitting graphs"""

        fm = self.file_manager
        result_dir = fm.get_dir(results_dir_tag)
        analys_dir = fm.get_dir(analysis_dir_tag)

        self.log.subheading('Model improvement graphs')

        labs = ['Input']   * ('R-free (Input)'   in table.columns) + \
               ['Fitted']  * ('R-free (Fitted)'  in table.columns) + \
               ['Refined'] * ('R-free (Refined)' in table.columns)

        # ------------------------------------------------------------>
        # Raw scores
        # ------------------------------------------------------------>
        # Histogram of the R-FREE for input+fitted+refined B-factors
        filename = fm.add_file(file_name='r-free-by-resolution.png', file_tag='rfree_v_resolution', dir_tag=results_dir_tag)
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
        filename = fm.add_file(file_name='r-work-by-resolution.png', file_tag='rwork_v_resolution', dir_tag=results_dir_tag)
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
        filename = fm.add_file(file_name='r-gap-by-resolution.png', file_tag='rgap_v_resolution', dir_tag=results_dir_tag)
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
        filename = fm.add_file(file_name='r-values-change-by-resolution.png', file_tag='all_rvalues_v_resolution', dir_tag=results_dir_tag)
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
        filename = fm.add_file(file_name='resolution-vs-rmsds-and-bfactors.png', file_tag='resolution_v_rmsds_scatter', dir_tag=analysis_dir_tag)
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
        filename = fm.add_file(file_name='r-values-vs-rmsds-and-bfactors.png', file_tag='rvalues_v_rmsds_scatter', dir_tag=analysis_dir_tag)
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

    def tidy_output_folder(self, actions):
        """Remove MTZ files, etc"""
        if not actions: return
        self.log.heading('Tidying/Cleaning output folder')
        if 'delete_mtzs' in actions:
            self.log.subheading('Deleting unnecessary MTZ files')
            for i_mod, m in enumerate(self.models):
                for f in [m.o_mtz, m.r_mtz]:
                    if (f is not None) and os.path.exists(f):
                        self.log('Removing file: {}'.format(f))
                        os.remove(f)
        if 'compress_logs' in actions:
            self.log.subheading('Compressing log files')
            from bamboo.common.file import compress_file
            main_log = os.path.abspath(self.log.log_file)
            for log in glob.glob(os.path.join(self.file_manager.get_dir('logs'),'*.log')):
                # Do not compress main log file
                if os.path.abspath(log) == main_log:
                    continue
                self.log('Compressing: {}'.format(log))
                zip_file = compress_file(log, delete_original=True)

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
                patch.set_edgecolor('k')
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
    def multi_histogram(filename, x_vals, titles, x_labs, rotate_x_labels=True, shape=None, n_bins=30, x_lim=None):
        """Generate standard histogram"""

        if shape is not None:
            nrow, ncol = shape
        else:
            nrow, ncol = (len(x_vals), 1)
        assert nrow*ncol == len(x_vals)

        fig, axes = pyplot.subplots(nrows=nrow, ncols=ncol, sharey=False)

        # Fix for 1 subplot
        if nrow==ncol==1: axes = numpy.array([axes])

        for i, axis in enumerate(axes.flatten()):

            # Select x-axis limits
            if x_lim is not None:
                assert len(x_lim) == 2
                i_x_lim = [min(numpy.min(x_vals[i]), x_lim[0]) if (x_lim[0] is not None) else numpy.min(x_vals[i]),
                           max(numpy.max(x_vals[i]), x_lim[1]) if (x_lim[1] is not None) else numpy.max(x_vals[i])]
            else:
                i_x_lim = [numpy.min(x_vals[i]), numpy.max(x_vals[i])]
            # Make sure does not have zero width
            if (i_x_lim[1] - i_x_lim[0]) < 1e-6:
                i_x_lim[0] -= 0.01*abs(i_x_lim[0])
                i_x_lim[1] += 0.01*abs(i_x_lim[1])

            axis.set_title(titles[i])
            axis.hist(x=x_vals[i], bins=n_bins, range=i_x_lim)
            axis.set_xlabel(x_labs[0])
            axis.set_ylabel('Count')
            axis.set_xlim(i_x_lim)
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
        axis.boxplot(y_vals,
                     labels=x_labels,
                     edgecolor='black',
                     showmeans=True)
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
                    reverse_legend_order=False,
                    colour_space=(0,1)):
        """Plot stacked bar plots for a series of hierarchies (plotted values are the average B-factors of each residue of the hierarchies)"""

        # TODO MERGE stacked_bar and multi_bar with mode='stacked' mode='grid'

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
            colors = pyplot.cm.rainbow(numpy.linspace(colour_space[0],colour_space[1],len(sel_hs)))
            hatchs = itertools.cycle(['//', 'x', '\\'])

            # Create the output figure
            fig, axis = pyplot.subplots(nrows=1, ncols=1)
            if title is not None: axis.set_title(label=str(title))
            axis.set_xlabel('Residue')
            axis.set_ylabel(y_lab)
            if y_lim:
                axis.set_ylim(y_lim)

            # Iterative though hierarchies
            handles = []
            for i_h, h in enumerate(sel_hs):
                # Extract b-factors from this hierarchy
                y_vals = numpy.array([numpy.mean(rg.atoms().extract_b()) for rg in h.residue_groups()])
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
            # Axis ticks & labels
            x_ticks = numpy.arange(1, len(x_labels)+1, max(1.0, numpy.floor(len(x_labels)/20)))
            #x_ticks = numpy.unique(numpy.floor(numpy.linspace(1,len(x_labels)+1,20)))
            axis.set_xticks(x_ticks)
            axis.set_xticklabels([x_labels[int(i)] if (i<len(x_labels)) and (float(int(i))==i) else '' for i in axis.get_xticks()])
            # Rotate axis labels
            if rotate_x_labels: pyplot.setp(axis.get_xticklabels(), rotation=90)

            # Plot boundaries
            if v_line_hierarchy is not None:
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
    def multi_hierarchy_bar(prefix,
                            hierarchies,
                            titles,
                            shape=None,
                            y_labs='Isotropic B ($\AA^2$)',
                            v_line_hierarchy=None,
                            rotate_x_labels=True):
        """Plot multiple bar plots for a series of hierarchies (plotted values are the average B-factors of each residue of the hierarchies)"""

        # Check titles same lengths
        titles = list(titles)
        assert len(hierarchies) == len(titles)
        # Check y-labs same lengths
        if isinstance(y_labs, str):
            y_labs = [y_labs]*len(hierarchies)
        assert len(hierarchies) == len(y_labs)
        # Get and check the shape of the subplots
        if shape is not None:
            nrow, ncol = shape
        else:
            nrow, ncol = (len(hierarchies), 1)
        assert len(hierarchies) == nrow*ncol

        # Extract one of the hierarchies as the "master"
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

            # Extract y-vals for each structure
            y_vals = []
            for i,h in enumerate(sel_hs):
                y_vals.append([numpy.mean(rg.atoms().extract_b()) for rg in h.residue_groups()])
            y_vals = numpy.array(y_vals)
            # Skip chains with no Bs
            if not y_vals.any():
                continue

            # Filename!
            filename = prefix + '-chain_{}.png'.format(chain_id)

            # Create x-values for each residue starting from 1
            x_vals = numpy.array(range(len(list(sel_mh.residue_groups()))))+1
            x_labels = ['']+[ShortLabeller.format(rg) for rg in sel_mh.residue_groups()]

            # Create the output figure
            fig, axes = pyplot.subplots(nrows=nrow, ncols=ncol)
            # Fix for 1 subplot
            if nrow==ncol==1: axes = numpy.array([axes])

            # Iterative though axes and plot
            handles = []
            for i, axis in enumerate(axes):
                # Set labels & title
                axis.set_ylabel(y_labs[i])
                axis.set_title(label=titles[i])
                # Plot the bar
                hdl = axis.bar(left=(x_vals-0.5),
                               height=y_vals[i],
                               width=1.0,
                               edgecolor='black')
                # Plot boundaries
                if v_line_hierarchy is not None:
                    v_lines = numpy.where(numpy.array([max(rg.atoms().extract_b()) for rg in v_line_hierarchy.select(sel).residue_groups()], dtype=bool))[0] + 1.5
                    for val in v_lines: axis.axvline(x=val, ls='dotted')

            # Make sure only the last axes have x-labels
            for a in axes[-ncol:]:
                a.set_xlabel('Residue')
                a.set_xticklabels([x_labels[int(i)] if (i<len(x_labels)) and (float(int(i))==i) else '' for i in axis.get_xticks()])

            # Format and save
            fig.tight_layout()
            fig.savefig(filename)#, dpi=300)
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
        # Create an overall B-iso TOTAL line
        hdl = ax.bar(left   = x_vals,
                     height = table['b_iso (total)'],
                     width  = 0.4,
                     align  = 'center',
                     label  = 'total')
        handles.append(hdl)
        # Create B-iso lines for each LEVEL
        hdl = ax.plot(x_vals,
                      table['b_iso (level)'].values,
                      'ko:',
                      label = 'level')
        handles.extend(hdl)
        ax.errorbar(x     = x_vals,
                    y     = table['b_iso (level)'].values,
                    fmt   = 'none',
                    yerr  = ((table['b_iso (level)'] - table['b_min (level)']).values,
                             (table['b_max (level)'] - table['b_iso (level)']).values),
                    ecolor     = 'k',
                    elinewidth = 1,
                    capthick   = 1,
                    capsize    = int(10./len(x_vals)))
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
                 dataset_weights=None,
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

        # Masks to exclude datasets or atoms from optimisation
        self.dataset_mask = numpy.ones(observed_uij.shape[0], dtype=bool)
        self.atomic_mask  = numpy.ones(observed_uij.shape[1], dtype=bool)

        # Weights for each dataset
        self.dataset_weights = dataset_weights

        # Level labels, and grouping for each level
        self.level_labels = level_labels if level_labels else ['Level {}'.format(i) for i in xrange(1, len(level_array)+1)]
        self.level_array = level_array

        assert len(self.level_labels) == len(self.level_array)

        # Calculate the tls origins based on the optimial segmentation
        self.tls_origins_hash, self.tls_origins = self.get_non_segmenting_tls_origins_and_partitions()
        assert self.tls_origins_hash.shape == (observed_xyz.shape[1],)                # n_atoms
        assert self.tls_origins.shape[0] == len(numpy.unique(self.tls_origins_hash))  # n_partitions
        assert self.tls_origins.shape[1] == observed_xyz.shape[0]                     # n_datasets
        assert self.tls_origins.shape[2] == 3                                         # ... n_dim

        # Series of objects to fit the Uij TLS groups
        self.levels = []
        for idx, (lab, group_idxs) in enumerate(zip(self.level_labels, self.level_array)):
            # Create fitter for each level
            self.levels.append(MultiDatasetUijTLSGroupLevel(observed_uij = observed_uij,
                                                            observed_xyz = observed_xyz,
                                                            tls_origins      = self.tls_origins,
                                                            tls_origins_hash = self.tls_origins_hash,
                                                            group_idxs = group_idxs,
                                                            n_tls = self.params.tls_models_per_tls_group,
                                                            index = idx+1,
                                                            label = lab,
                                                            weights = dataset_weights,
                                                            amplitudes = self.params.tls_amplitude_model,
                                                            params = self.params.optimisation,
                                                            log = self.log))

        # One object to fit all the Uij residuals
        self.residual = MultiDatasetUijResidualLevel(observed_uij = observed_uij,
                                                     index = len(self.levels)+1,
                                                     label = 'residual',
                                                     weights = dataset_weights,
                                                     params = self.params.optimisation,
                                                     log = self.log)

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

    def get_non_segmenting_tls_origins_and_partitions(self):
        """Calculate the largest groups that are possible without any partitioning on any level existing between two groups"""

        self.log.subheading('Identifying partitioning of atoms that does not split any TLS groups')

        array = self.level_array
        labels = self.level_labels

        atom_xyzs = self.observed_xyz
        atom_hash = None
        atom_coms = []

        # Find the levels with the largest coverage (number of atoms in groups in this level)
        self.log('Calculating the number of atoms that are in groups for each level:')
        n_atoms = array.shape[1]
        n_atoms_per_level = [(l>=0).sum() for l in array]
        self.log('Atoms covered per level: \n\t{}'.format('\n\t'.join(['Level {} -> {} atoms'.format(i+1, n) for i,n in enumerate(n_atoms_per_level)])))
        levels_with_max_coverage = [i_l for i_l, n in enumerate(n_atoms_per_level) if n==n_atoms]
        self.log('Only considering partitionings that cover all atoms ({} atoms).'.format(n_atoms))
        self.log('Considering partitionings from levels: \n\t{}'.format('\n\t'.join(['{} ({})'.format(l+1, labels[l]) for l in levels_with_max_coverage])))
        # Find how many groups per level, and find level indices with fewest groups
        self.log('Looking for which of these levels has the smallest number of groups')
        n_groups_per_level = [len(numpy.unique(l[l>=0])) for l in array[levels_with_max_coverage]]
        self.log('Groups per level: \n\t{}'.format('\n\t'.join(['Level {} -> {}'.format(levels_with_max_coverage[i]+1, n) for i,n in enumerate(n_groups_per_level)])))
        levels_with_min_ngroups = [levels_with_max_coverage[i] for i, n in enumerate(n_groups_per_level) if n==min(n_groups_per_level)]
        self.log('Only considering partitionings with the smallest number of groups ({} groups).'.format(min(n_groups_per_level)))
        self.log('Considering partitionings from levels: \n\t{}'.format('\n\t'.join(['{} ({})'.format(i_l+1, labels[i_l]) for i_l in levels_with_min_ngroups])))
        # Set to False just because (even though it is immediately overwritten)
        found_valid_level = False
        # Test each level to see if it splits tls groups in other levels
        for i_l in levels_with_min_ngroups:
            self.log.bar()
            self.log('Checking to see if the partitioning on level {} ({}) cuts any TLS groups on other levels'.format(i_l+1, labels[i_l]))
            # Assume level is valid until proven otherwise
            found_valid_level = True
            # Extract groups for this level
            level_groups = array[i_l]
            # Check to see if groups on other levels occupy the non_sel positions (the "null space" of this level)
            level_sel = (level_groups >= 0)
            have_no_group_this_level = array[:,numpy.bitwise_not(level_sel)]
            have_group_other_levels = (have_no_group_this_level >= 0)
            # group members outside selection
            if have_group_other_levels.sum():
                found_valid_level = False
                self.log('> There are atoms in groups on other levels that do not have a group on this level.')
                self.log('> This partitioning is not appropriate.')
                continue
            else:
                self.log('> All atoms in other levels are contained by groups on this level.')
            # Get the group idxs (labels/names) for this level
            g_idxs = numpy.unique(level_groups)
            # Iterate through groups and get selection for each
            self.log('> Checking that every TLS group on other levels is contained within only one group on this level')
            for g in g_idxs:
                if g < 0: continue
                self.log('  ...checking this condition for atoms in group {} of this level'.format(g))
                # Get the selection for this group
                g_sel = (level_groups==g)
                # Get the chunk of values for other levels -- for this selection, and NOT this selection
                g_this = array[:,g_sel]
                g_other = array[:,numpy.bitwise_not(g_sel)]
                # Check condition for all other levels
                for i_l_test in range(len(array)):
                    # Skip checking against itself
                    if i_l==i_l_test: continue
                    self.log('     ...checking this condition for the partitioning on level {}.'.format(i_l_test+1))
                    # Check if groups in this selection in this level are present in the level outside this selection
                    if set(g_this[i_l_test]).intersection(set(g_other[i_l_test])):
                        self.log('      ...there is a group on level {} whose atoms are split between this group and another'.format(i_l_test+1))
                        self.log('> This partitioning is not appropriate.')
                        found_valid_level = False
                        break
                # If invalid it's not necessary to check the rest
                if found_valid_level is False:
                    break
            # All checks passed - return this
            if found_valid_level is True:
                self.log.bar()
                self.log('This level -- level {} ({}) -- has a partitioning that is compatible with all other levels.'.format(i_l+1,labels[i_l]))
                self.log('Using this partitioning.')
                atom_hash = copy.copy(array[i_l])
                break
        # No valid levels -- create over-arching partition
        if found_valid_level is False:
            self.log('No suitable level found -- using one partition containing all atoms.')
            atom_hash = numpy.ones(array.shape[1], dtype=int)
        # Renumber the array to start from 0
        for i,v in enumerate(sorted(numpy.unique(atom_hash))):
            if v < 0: continue
            atom_hash[atom_hash==v] = i
        # Iterate through and calculate origins for each dataset for each group
        self.log.bar(True, False)
        self.log('Calculating Centres-of-mass for each partition')
        self.log.bar()
        for i in sorted(numpy.unique(atom_hash)):
            if i<0: continue
            self.log('\nCentres-of-mass for atoms in Partition {}'.format(i+1))
            # Get the selection for this group
            i_sel = (atom_hash == i)
            # Get the atoms in this group
            xyzs = atom_xyzs[:,i_sel,:]
            # Calculate centre of mass for each dataset
            com = numpy.mean(xyzs, axis=1)
            self.log('')
            for j, c in enumerate(com):
                self.log('> Dataset {:>3d}: {}'.format(j+1, c))
            atom_coms.append(com)
        # Convert to numpy array and transpose
        atom_coms = numpy.array(atom_coms)

        return atom_hash, atom_coms

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

    def set_tracking(self, table, csv_path, png_path):
        self.tracking_data = table
        self.tracking_csv = csv_path
        self.tracking_png = png_path

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
            # average values
            b_iso_sel = numpy.mean(uij_to_b(uij_sel))
            u_iso_sel = b_iso_sel / EIGHT_PI_SQ
            # min/max values
            b_min_sel = numpy.min(uij_to_b(uij_sel))
            b_max_sel = numpy.max(uij_to_b(uij_sel))
        else:
            # average values
            b_iso_sel = 0.0
            u_iso_sel = 0.0
            # min/max values
            b_min_sel = numpy.nan
            b_max_sel = numpy.nan

        # Calculate U-iso & B-iso for complete model
        b_iso_tot = numpy.mean(uij_to_b(uij_tot))
        u_iso_tot = b_iso_tot / EIGHT_PI_SQ

        # Create human-readable cycle number
        cycle_lab = str(i_cycle+1)

        # Add to tracking table
        self.tracking_data.loc[len(self.tracking_data.index)] = {'cycle' : cycle_lab,
                                                                 'level' : level_lab,
                                                                 'rmsd'  : numpy.round(rmsd,3),
                                                                 'u_iso (level)' : numpy.round(u_iso_sel,3),
                                                                 'b_iso (level)' : numpy.round(b_iso_sel,3),
                                                                 'b_min (level)' : numpy.round(b_min_sel,3),
                                                                 'b_max (level)' : numpy.round(b_max_sel,3),
                                                                 'u_iso (total)' : numpy.round(u_iso_tot,3),
                                                                 'b_iso (total)' : numpy.round(b_iso_tot,3)}

        self.log(self.tracking_data.loc[len(self.tracking_data.index)-1].to_string())
        # Dump to csv
        self.tracking_data.to_csv(self.tracking_csv)
        # Make plots
        MultiDatasetUijPlots.tracking_plots(table    = self.tracking_data,
                                            filename = self.tracking_png)

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
                fitted_uij = fitter.run(n_cpus=n_cpus, n_cycles=n_micro_cycles)
                # Validate and store
                for i_dst, dataset_uij in enumerate(fitted_uij):
                    uij_valid = uij_positive_are_semi_definite(uij=dataset_uij, tol=self.params.precision.uij_tolerance)
                    if uij_valid.sum() > 0.0:
                        err_msg = '\n\t'.join(['Atom {:>5d} (group {}): {}'.format(i+1, fitter.get_atom_group(index=i).label, dataset_uij[i]) for i in numpy.where(uij_valid)[0]])
                        raise Failure('Uijs for dataset {} are not positive-semi-definite ({} atoms)\n\t{}'.format(i_dst, uij_valid.sum(), err_msg))
                # Store in output array
                fitted_uij_by_level[i_level] = fitted_uij
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
        s += '\nNumber of TLS models TLS per group: {}'.format(self.params.tls_models_per_tls_group)
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
    _chunksize = None

    def n_params(self, non_zero=False):
        return sum([f.n_params(non_zero=non_zero) for i,g,f in self])

    def run(self, n_cycles=1, n_cpus=1):
        # Check to see if multiple cpus are available per job
        n_cpus_per_job = max(1, n_cpus//self._n_obj)
        # Create job objects
        jobs = [(fitter, {'n_cycles':n_cycles, 'n_cpus':n_cpus_per_job, 'silent':(n_cpus>1)}) for (i, sel, fitter) in self]
        # Drop CPUs if not required
        n_cpus = min(len(jobs), n_cpus)
        # Only do multiprocessing if actually needed
        if n_cpus > 1:
            # Run jobs in parallel
            self.log('Running {} job(s) in {} process(es) [and {} cpu(s) per process]'.format(len(jobs), n_cpus, n_cpus_per_job))
            self.log('')
            workers = NonDaemonicPool(n_cpus)
            finished_jobs = workers.map(func=wrapper_fit, iterable=jobs, chunksize=self._chunksize)
            workers.close()
        else:
            self.log('Running {} job(s) [with {} cpu(s)]'.format(len(jobs), n_cpus_per_job))
            self.log('')
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
            # Print heading for every group for other levels
            if self.label != 'residual':
                self.log.subheading('Results - level {} ({}) - group {}'.format(self.index, self.label, ret_job.label))
            # Print summaries (in main log)
            self.log(ret_job.summary(show=False))
            # Store the returned job
            self.fitters[i] = ret_job
        # Report errors
        if errors:
            self.log.heading('Fatal errors returned during optimisation')
            for fitter, e in errors:
                self.log.bar()
                self.log('Level "{}", Group "{}": error returned'.format(self.label, fitter.label))
                self.log.bar()
                self.log(e)
                self.log.bar()
                self.log('Log file written to {}'.format(fitter.log.log_file))
            self.log.bar()
            raise Failure('Errors raised during optimisation. Messages printed above')
        # Return the fitted B-factors
        return self.extract()

    def set_target_uij(self, target_uij):
        assert target_uij.shape == self._uij_shape
        for i, sel, fitter in self:
            fitter.set_target_uij(target_uij=target_uij[:,sel])

    def set_atomic_mask(self, mask):
        n_min = 4
        for i, sel, fitter in self:
            # Apply the mask if contains more than n atoms
            if sum(mask[sel]) > n_min:
                fitter.set_atomic_mask(mask[sel])
            # Else do nothing
            else:
                self.log('Level "{}", Group "{}": Attempting to apply mask of {} atoms or less -- not applying mask'.format(self.label, fitter.label, n_min))

    def set_dataset_mask(self, mask):
        for i, sel, fitter in self:
            fitter.set_dataset_mask(mask)

class MultiDatasetUijTLSGroupLevel(_MultiDatasetUijLevel):

    _chunksize = 1

    def __init__(self,
                 observed_uij,
                 observed_xyz,
                 tls_origins,
                 tls_origins_hash,
                 group_idxs,
                 n_tls=1,
                 index=0,
                 label=None,
                 weights=None,
                 amplitudes=None,
                 params=None,
                 log=None):

#        if log is None: log = Log(verbose=False)
        self.log = log

        assert observed_uij.shape[:2] == observed_xyz.shape[:2]
        assert observed_xyz.shape[2]  == 3
        assert observed_uij.shape[1]  == len(group_idxs)
        assert observed_uij.shape[2]  == 6
        assert tls_origins.shape[1] == observed_uij.shape[0]
        assert tls_origins.shape[2] == 3
        assert tls_origins_hash.shape[0] == observed_xyz.shape[1]

        if weights is not None:
            assert weights.shape == (observed_uij.shape[0],)

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
                           #silent   = False)
                           silent   = not (self._n_groups==1))
            # Decide which tls_origins to use fo this
            i_origin = numpy.unique(tls_origins_hash[sel])
            assert len(i_origin) == 1
            observed_com = tls_origins[i_origin[0]]
            assert observed_com.shape == (observed_xyz.shape[0], 3)
            # Create fitter object
            self.fitters[i] = MultiDatasetUijTLSOptimiser(target_uij = observed_uij[:,sel],
                                                          atomic_xyz = observed_xyz[:,sel],
                                                          atomic_com = observed_com,
                                                          n_tls = n_tls,
                                                          label = '{:4d} of {:4d}'.format(i, self._n_groups),
                                                          weights = weights,
                                                          amplitudes = amplitudes,
                                                          params = params,
                                                          log = ls)

    def __iter__(self):
        for i in numpy.unique(self.group_idxs):
            if i==-1: continue
            yield (i, (self.group_idxs==i), self.fitters.get(i, None))

    def get_atom_group(self, index):
        """Get the group of an atom (by index)"""
        group_n = self.group_idxs[index]
        return self.fitters.get(group_n)

    def n_groups(self):
        return self._n_groups

    def extract(self, average=False):
        fitted_uij = numpy.zeros(self._uij_shape)
        for i, sel, fitter in self:
            fitted_uij[:,sel] = fitter.extract()
        if average:
            fitted_uij = fitted_uij.mean(axis=0)
        return fitted_uij

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

    def zero_amplitudes(self, models, components=None):
        # Iterate through each group/partition
        for i, sel, fitter in self:
            fitter.parameters().zero_amplitudes(models=models, components=components)

    def zero_matrices(self, models, components=None):
        # Iterate through each group/partition
        for i, sel, fitter in self:
            fitter.parameters().zero_matrices(models=models, components=components)

class MultiDatasetUijResidualLevel(_MultiDatasetUijLevel):

    def __init__(self,
                 observed_uij,
                 index=-1,
                 label=None,
                 weights=None,
                 params=None,
                 log=None):

#        if log is None: log = Log(verbose=True)
        self.log = log

        if weights is not None:
            assert weights.shape == (observed_uij.shape[0],)

        self.index = index
        self.label = label if label else 'Level {} (residual)'.format(index)

        self._uij_shape = observed_uij.shape

        self._n_atm = self._uij_shape[1]
        self._n_obj = self._n_atm

        # Create one log for all fitters
        ls = LogStream(log_file = os.path.splitext(self.log.log_file)[0]+'-level-residual.log',
                       verbose  = False,
                       silent   = False)
        # Create a fitter for each atom
        self.fitters = {}
        for i, sel, f in self:
            assert f is None
            # Create fitter object
            self.fitters[i] = MultiDatasetUijAtomOptimiser(target_uij = observed_uij[:,sel],
                                                           label = 'atom {:5d} of {:5d}'.format(i,self._n_atm),
                                                           weights = weights,
                                                           params = params,
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
                 label='',
                 weights=None,
                 params=None,
                 log=None):

        assert params.penalties is not None

#        if log is None: log = Log(verbose=False)
        self.log = log

        # Store verboseness and silence as overridden in some functions
        self.verbose = self.log.verbose
        self.silent  = self.log.silent

        self._n_prm = 0
        self._n_dst = 0
        self._n_atm = 0

        self.target_uij = target_uij

        self.label = label

        if weights is None:
            self.weights = numpy.ones(target_uij.shape[0])
        else:
            self.weights = numpy.array(weights)
        assert self.weights.shape == (target_uij.shape[0],)

        self.penalties = params.penalties
        self.penalty = UijPenalties(mdl_fixed=params.penalties.invalid_tls_values.a,
                                    mdl_slope=params.penalties.invalid_tls_values.b,
                                    amp_fixed=params.penalties.invalid_amplitudes.a,
                                    amp_slope=params.penalties.invalid_amplitudes.b,
                                    uij_fixed=params.penalties.invalid_uij_values.a,
                                    uij_slope=params.penalties.invalid_uij_values.b,
                                    ovr_fixed=params.penalties.over_target_values.a,
                                    ovr_slope=params.penalties.over_target_values.b)

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

    def _optimise(self, tolerance=1e-3, running_summary=False):
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
        # Get additional input required for the simplex generation
        simplex_kw = self._simplex_input()
        simplex_kw.update(sel_dict)
        # Create simplex for these parameters
        opt_simplex = self.simplex.get_simplex(start=sel_vals, **simplex_kw)
        #
        if self._running_summary:
            fmt = '{:>+10.5f}'
            self.log.subheading('Starting simplex')
            self.log('Starting Point:')
            self.log(','.join([fmt.format(v) for v in sel_vals]))
            self.log('Points on Simplex (relative to starting point):')
            for point in opt_simplex:
                self.log(','.join([fmt.format(v) for v in (point-sel_vals)]))
            self.log.subheading('Running optimisation')
        # Optimise these parameters
        assert self._mask_wgts.shape == (self._target_uij.shape[0],)
        optimised = simplex.simplex_opt(dimension = len(sel_vals),
                                        matrix    = map(flex.double, opt_simplex),
                                        evaluator = self.evaluator,
                                        tolerance = tolerance)
        # Extract and update current values
        self._inject(values=optimised.get_solution(), selection=sel_dict)
        if self._running_summary:
            self.log.bar()

    def _simplex_input(self):
        """Default function -- no additional variables"""
        return {}

    #===========================================+>
    # Public Functions
    #===========================================+>

    def copy(self):
        return copy.deepcopy(self)

    def get_atomic_mask(self):
        return self._mask_atom

    def get_dataset_mask(self):
        return self._mask_dset

    def get_dataset_weights(self):
        return self._mask_wgts

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
        self._mask_wgts = self.weights[mask]

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
                header = '[{}] -> ({:^12}, {:^10})'.format(', '.join(['{:>7}'.format('param') for r in parameters]), 'fit/rmsd', 'penalty')
                line = '-'*len(header)
                self.log(line, show=False)
                self.log(header, show=False)
                self.log(line, show=False)
        # Return now if physical penalty if non-zero to save time
        if ppen > 0.0:
            if self._running_summary:
                self.log('[{}] -> ({:>12}, {:10.3f})'.format(', '.join(['{:+7.4f}'.format(r) for r in parameters]), 'UNPHYSICAL', ppen), show=False)
            return ppen
        # Get the fitted uijs (including masked atoms)
        self._update_fitted()
        # Calculate RMSD
        sd = numpy.power(self._target_uij-self._fitted_uij, 2)
        msd = numpy.average(sd, axis=0, weights=self._mask_wgts).mean()
        rmsd = numpy.sqrt(msd)
        # Calculate fitting penalties (add to rmsd)
        fpen = self._fitting_penalties(uij_fit=self._fitted_uij, uij_obs=self._target_uij)
        # Total target function value
        tot_val = rmsd + fpen
        # Update minima
        if tot_val < self.optimisation_rmsd+self.optimisation_penalty:
            self.optimisation_rmsd    = rmsd
            self.optimisation_penalty = fpen
        if self._running_summary:
            self.log('[{}] -> ({:12.8f}, {:10.6f})'.format(', '.join(['{:+7.4f}'.format(r) for r in parameters]), rmsd, fpen), show=False)
        return tot_val

class MultiDatasetUijAtomOptimiser(_UijOptimiser):

    def __init__(self,
                 target_uij,
                 label='',
                 weights=None,
                 params=None,
                 log=None):
        super(MultiDatasetUijAtomOptimiser, self).__init__(target_uij=target_uij,
                                                           label=label,
                                                           weights=weights,
                                                           params=params,
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
                 weights=None,
                 amplitudes=None,
                 params=None,
                 log=None):
        super(MultiDatasetUijTLSOptimiser, self).__init__(target_uij=target_uij,
                                                          label=label,
                                                          weights=weights,
                                                          params=params,
                                                          log=log)

        # Store the centre of mass of the atoms (for the rotation/screw components)
        self.atomic_xyz = atomic_xyz
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
        self._parameters = MultiDatasetTLSModelList(n_mdl = self._n_tls,
                                                    n_dst = self._n_dst,
                                                    amplitudes = amplitudes,
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
        self.simplex = TLSSimplex(vibration = params.step_size.vibration,
                                  libration = params.step_size.libration,
                                  angle     = params.step_size.angle,
                                  amplitude = params.step_size.amplitude)

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

        tls_penalties = []
        amp_penalties = []

        for mode in self._parameters:

            # Calculate model penalties (if required -- not required if not optimising model...)
            if self._opt_dict['model']:
                # Penalise physically-invalid TLS models
                tls_penalties.append(self.penalty.tls_params(model=mode.model))

            # Calculate dataset penalties
            if self._opt_dict['amplitudes']:
                # Penalties for combinations of amplitudes that lead to non-physical TLS models
                datasets = self._opt_dict['i_dst']
                # Penalties for negative amplitudes, etc
                amp_penalties.append(self.penalty.amplitudes(values=mode.amplitudes.values))
                # XXX remark: next line commented out because floating point/rounding errors
                # XXX remark: prevent amplitude optimisation for small TLS values
                #amp_penalties.extend([self.penalty.tls_params(model=model) for model in mode.expand(datasets=datasets)])

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
            assert not set(''.join(components) if isinstance(components, list) else components).difference('TLS')
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
                values.append(mdl.model.get(components=s.components, include_szz=False))
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
            n_tls_params = tls_str_to_n_params(s.components, include_szz=False)
            for i in s.i_mdl:
                mdl = self._parameters.get(index=i)
                mdl.model.set(vals = [values.popleft() for _ in xrange(n_tls_params)],
                              components = s.components,
                              include_szz = False)
        if s.amplitudes is True:
            for i in s.i_mdl:
                mdl = self._parameters.get(index=i)
                mdl.amplitudes.set(vals = numpy.array([values.popleft() for _ in xrange(s.n_dst*s.n_cpt)]).reshape((s.n_dst, s.n_cpt)),
                                   components = s.components,
                                   datasets = s.i_dst)
        assert len(values) == 0, 'not all values have been popped'

    #===========================================+>

    def _simplex_input(self):
        """Return the current values of the TLS matrices"""
        tls_models = [m.model for m in self.parameters()]
        return {'tls_models':tls_models}

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

        # Cumulative i_tls indices
        mdl_cuml = []

        # Optimise!
        for i_cycle in xrange(n_cycles):
            self.log.subheading('Group {} - Optimisation cycle {} of {}'.format(self.label, i_cycle+1, n_cycles))

            # Optimise one TLS "mode" at time
            for i_tls in xrange(self._n_tls):

                # How should T-L-S be optimised - if all values are zero, optimise all separately then as those for each amplitude/component set
                all_cpts = self.parameters().get(index=i_tls).all_components()
                # Are all matrix values zero?
                null_model = not self.parameters().get(index=i_tls).model.any(components=all_cpts)
                # If null model then optimise each matrix individually (for speed) and then optimise all together
                cpt_groups = list(all_cpts)*null_model + [all_cpts]*((not null_model) or (len(all_cpts) > 1))

                self.log.subheading('Optimising TLS model {} of {}'.format(i_tls+1, self._n_tls))
                self.log('Optimising using {} atoms'.format(len(opt_atom_mask)))
                self.log('Optimising using {} datasets'.format(len(opt_dset_mask)))
                self.log.bar(True, True)
                self.log('Optimisation order for TLS matrix values: \n\t{}'.format('\n\t'.join(['{}) {}'.format(i+1, c) for i, c in enumerate(cpt_groups)])))
                self.log('Optimising all TLS amplitude values: \n\t{}'.format('\n\t'.join(['{}) {}'.format(i+1, c) for i, c in enumerate(self.parameters().get(index=i_tls).component_sets())])))

                # If starting from zero values, set translation matrix to sensible starting values
                if (null_model is True):
                    self.log.bar(True,True)
                    self.log('Identifing a starting value for the T-matrix')
                    # Get the unfitted Uij (target uij minus uij already fitted from other models)
                    uij_vals = self.target_uij[opt_dset_mask][:,opt_atom_mask] - self._extract(mask_datasets=opt_dset_mask, mask_atoms=opt_atom_mask)
                    assert uij_vals.shape == (len(opt_dset_mask), len(opt_atom_mask),6)
                    uij_eigs = numpy.apply_along_axis(uij_eigenvalues, axis=2, arr=uij_vals)
                    assert uij_eigs.shape == (len(opt_dset_mask),len(opt_atom_mask),3)
                    uij_eigs_min = numpy.min(uij_eigs, axis=2)
                    assert uij_eigs_min.shape == (len(opt_dset_mask),len(opt_atom_mask))
                    #uij_vols = numpy.array([matrix.sym(sym_mat3=u).determinant() for u in uij_vals])
                    min_eig = numpy.min(uij_eigs_min)
                    if min_eig > 0.0:
                        dst, atm = zip(*numpy.where(uij_eigs_min==min_eig))[0]
                        uij_start = tuple(uij_vals[dst,atm].tolist())
                        self.log('Minimum uij eigenvalue: {}'.format(round(min_eig,6)))
                        self.log('Atom with minimal uij eigenvalue: Dataset {}, Atom {}'.format(dst, atm))
                        self.log('Atomic uij: {}'.format(tuple([round(v,6) for v in uij_start])))
                        self.log('Initialising the T-matrix with an isotropic uij of this size')
                        uij_start = (round(min_eig,6),)*3 + (0.0,)*3
                        self.log('Setting initial T-matrix values to {}'.format(uij_start))
                        self.parameters().get(index=i_tls).model.set(vals=uij_start, components='T')
                        self.log(self.parameters().get(index=i_tls).model.summary())
                        # Normalise matrices to give Uijs of approximately xA
                        self.parameters().normalise_by_matrices(xyzs=self.atomic_xyz, origins=self.atomic_com, target=1.0)
                    else:
                        self.log('There is an atom with negative eigenvalues (value {})'.format(min_eig))
                        self.log('Starting with a T-matrix with zero values')
                        self.log('Looking for an appropriate scale for the amplitudes instead')
                        for perc in range(10,110,10):
                            scale = numpy.percentile(uij_eigs_min, q=perc)
                            if scale > 0.0:
                                self.log('Using the {}th percentile ({})'.format(perc,scale))
                                break
                        if scale <= 0.0:
                            self.log('No minimal eigenvalues are positive.')
                            scale = 1e-3
                            self.log('Scaling amplitudes to minimal size ({})'.format(scale))
                        self.log('Scaling the amplitudes...')
                        self.parameters().get(index=i_tls).amplitudes.scale(multiplier=scale)
                        self.log(self.parameters().get(index=i_tls).summary())

                # Append to running list
                if i_tls not in mdl_cuml: mdl_cuml.append(i_tls)

                # ---------------------------------->
                # Optimise each set of selected components separately
                # ---------------------------------->
                for cpts in cpt_groups:
                    # Skip S optimisation if T and L are zeros
                    if (cpts == 'S') and not (self.parameters().get(index=i_tls).model.any('T') and \
                                              self.parameters().get(index=i_tls).model.any('L')):
                        self.log('T and L matrices are zero -- not optimising S-matrix')
                        continue
                    self.log.subheading('Optimising {} parameters for TLS mode {}'.format(cpts, i_tls+1))
                    # Set penalty weights
                    self.log.bar()
                    self.log('Penalty weights:')
                    self.log.bar()
                    self.log(self.penalty.set_weights(ovr_fixed=self.penalties.over_target_values.a))
                    self.log.bar(False, True)
                    # Run optimisation
                    self.optimise_models(
                        models      = [i_tls],
                        components  = cpts,
                        datasets    = opt_dset_mask,
                        atoms       = opt_atom_mask)
                    # Log model summary
                    self.log(self.parameters().get(index=i_tls).model.summary())
                    # Log current RMSD and penalty
                    self.optimisation_summary()

                # Normalise matrices to give Uijs of approximately xA
                self.parameters().normalise_by_matrices(xyzs=self.atomic_xyz, origins=self.atomic_com, target=1.0)
                # Check that normalisation has not made any of the models invalid and correct if possible
                self.optimise_invalid_models()

                # ---------------------------------->
                # Optimise TLS amplitude parameters (all amplitudes!)
                # ---------------------------------->
                # Only optimise ampltidues if model values are non-zero
                if not self.parameters().get(index=i_tls).model.any(all_cpts):
                    self.log('{} matrices are all zero -- not optimising amplitudes'.format(all_cpts))
                    self.parameters().get(index=i_tls).reset()
                else:
                    # Report
                    self.log.subheading('Optimising TLS amplitudes for all datasets')
                    # Reset amplitudes to zero to prevent overparameterisation
                    #self.log('Resetting all amplitudes to zero')
                    #self.parameters().zero_amplitudes(models=[i_tls])
                    # Run optimisation
                    self.optimise_amplitudes(
                        models      = mdl_cuml,
                        components  = self.parameters().get(index=i_tls).component_sets(),
                        datasets    = None,
                        atoms       = range(self._n_atm),   # Optimise with all atoms to prevent over-sizing
                        n_cpus = n_cpus)
                    # Log amplitudes summary
                    self.log.bar(blank_before=True, blank_after=True)
                    self.log(self.parameters().get(index=i_tls).amplitudes.summary())

            # Reapply atom and dataset masks
            self.set_dataset_mask(opt_dset_mask)
            self.set_atomic_mask(opt_atom_mask)

            # End of cycle house-keeping
            self.log.bar()
            # Break optimisation if all matrices are zero -- not worth running following optimisation cycles
            if not self.parameters().any():
                self.log('All matrices have refined to zero -- optimisation finished.')
                self.log.bar()
                # Reset parameters just to be sure
                self.parameters().reset()
                break
            # Check to see if any amplitudes are negative
            self.log('Looking for negative TLS amplitudes')
            self.parameters().zero_negative_amplitudes()
            self.log.bar()
            # Identify modes that are zero and reset these
            self.log('Looking for zero-value modes')
            self.parameters().zero_null_modes()
            self.log.bar()

            # Show summary at end of cycle...
            self.log.subheading('End-of-cycle summary')
            self.summary(show=True)
            self.log('')

    def optimise_models(self, models, components, datasets=None, atoms=None, silent=False, tolerance=1e-3):
        """Optimise the matrices for combinations of models/componenets/datasets"""

        # Take all datasets/atoms if not provided
        if datasets is None: datasets = range(self._n_dst)
        if atoms is None:    atoms = range(self._n_atm)

        # Select variables for optimisation -- model only
        self._select(optimise_model      = True,
                     optimise_amplitudes = False,
                     models     = models,
                     components = components,
                     datasets   = datasets,
                     atoms      = atoms)

        # Run optimisation
        self._optimise(tolerance=tolerance,
                       running_summary=(not silent))

    def optimise_amplitudes(self, models, components, datasets=None, atoms=None, n_cpus=1, silent=False, tolerance=1e-3):
        """Optimise the amplitudes for combinations of models/components/datasets"""

        # Take all datasets/atoms if not provided
        if datasets is None: datasets = range(self._n_dst)
        if atoms is None:    atoms = range(self._n_atm)

        self.log('Setting up TLS amplitude optimisation')
        # Select variables for optimisation -- amplitudes only
        self._select(optimise_model      = False,
                     optimise_amplitudes = True,
                     models     = models,
                     components = components,
                     datasets   = None,         # This will be set individually
                     atoms      = atoms)
        # ---------------------------------->
        self.log('Running amplitude optimisations: --->')
        # Optimise all amplitudes dataset-by-dataset
        jobs = []
        for i_dst in datasets:
            # Select this dataset
            self._select(datasets=[i_dst])
            # Append to job list (create copy) or run optimisation
            if n_cpus > 1:
                jobs.append((self.copy(), tolerance))
            else:
                self._optimise(tolerance=tolerance,
                               running_summary=(not silent))
                self.log('> dataset {} of {} (rmsd {:.3f}; penalty {:.1f})'.format(i_dst+1, self._n_dst, self.optimisation_rmsd, self.optimisation_penalty))
        # Run parallel jobs and inject results
        if n_cpus > 1:
            # Workers for multiprocessing
            workers = DaemonicPool(n_cpus)
            # Report and map to workers
            self.log.subheading('Running {} jobs using {} cpus'.format(len(jobs), min(n_cpus,len(jobs))))
            finished_jobs = workers.map(func=wrapper_optimise, iterable=jobs)
            for i_dst in datasets:
                job = finished_jobs[i_dst]
                if isinstance(job, str):
                    self.log(job)
                    raise Failure('error returned')
                self._inject(values=job._values(), selection=job._select())
                self.log('> dataset {} of {} (rmsd {:.3f}; penalty {:.1f})'.format(i_dst+1, self._n_dst, job.optimisation_rmsd,  job.optimisation_penalty))
            # Let's be good
            workers.close()

    def optimise_invalid_models(self, n_cpus=1):
        """Check if models are invalid and attempt to correct them"""

        self.log.subheading('Checking for invalid TLS matrices')

        # Extract the current model delta
        orig_step_values = self.simplex.current_values()

        for i_tls, model in enumerate(self.parameters()):

            self.log.bar()
            self.log('Checking physicality of TLS model {}'.format(i_tls+1))
            self.log.bar()

            # Check that the core model is valid
            if model.model.is_valid():
                self.log('Core model is valid')
                continue

            # Model invalid -- try to fix
            self.log('Core model is invalid\n...attempting to re-optimise:')
            # Create copy of the model first
            model_cpy = model.copy()
            self.log('Before re-optimisation: \n{}'.format(model_cpy.model.summary()))

            # Try resetting S if T or L is zero
            if (numpy.abs(model.model.get('T')) < model.model.get_precision()).all() or \
               (numpy.abs(model.model.get('L')) < model.model.get_precision()).all():
                self.log('...T or L matrix are zero -- resetting S-matrix')
                model.model.set(vals=0.0, components='S')
            if model.model.is_valid():
                self.log('...core model is now valid')
                continue

            # Determine the range that the model values exist over (precision -> maximum value)
            log_min_model_delta     = int(numpy.ceil(-model.model.get_precision() / 2.0))
            log_max_model_delta     = int(numpy.floor(numpy.log10(numpy.max(numpy.abs(model.model.values)))/2.0))
            self.log('Min model step size:  {:e}'.format(10.**log_min_model_delta))
            self.log('Max model step size:  {:e}'.format(10.**log_max_model_delta))
            # Iteratively increase the model step until model becomes valid
            for step in 10.**(numpy.arange(log_min_model_delta, log_max_model_delta+1)):
                self.log('...re-optimising with model steps of {}'.format(step))
                step = float(step)
                # Set simplex step
                self.simplex.set_step_sizes(vibration=step,
                                            libration=step,
                                            angle=step)
                # Run optimisations with the small model step
                self.optimise_models(models     = [i_tls],
                                     components = 'TLS',
                                     datasets   = self.get_dataset_mask(),
                                     atoms      = self.get_atomic_mask(),
                                     silent     = False)
                # Check again if the core model is valid
                model = self.parameters().get(index=i_tls)
                if model.model.is_valid():
                    break
            if model.model.is_valid():
                self.log('...core model is now valid (after model step: {})'.format(step))
                self.log('Before re-optimisation: \n{}'.format(model_cpy.model.summary()))
                self.log('After re-optimisation: \n{}'.format(model.model.summary()))
                continue

            # Reset to zero and re-optimise
            self.log('...core model is still invalid: \n{}'.format(model.model.summary()))
            self.log('...resetting and redetermining TLS values')
            # Reset TLS matrices to zero
            model.model.reset(components='TLS')
            # Set model step to original value
            self.simplex.set_step_sizes(**orig_step_values)
            # Iterate through and optimise as originally
            for cpts in ['T','L','S','TLS']:
                self.optimise_models(models     = [i_tls],
                                     components = cpts,
                                     datasets   = self.get_dataset_mask(),
                                     atoms      = self.get_atomic_mask(),
                                     silent     = False)
            model = self.parameters().get(index=i_tls)
            if model.model.is_valid():
                self.log('...core model now valid (after resetting and reoptimising)')
                self.log('Before re-optimisation: \n{}'.format(model_cpy.model.summary()))
                self.log('After re-optimisation: \n{}'.format(model.model.summary()))
                rms_values = rms(model.model.values-model_cpy.model.values)
                self.log('...rms between initial model and "fixed" model: {}'.format(rms_values))
                continue

            # XXX XXX XXX
            raise Failure('Failed to fix core model. \n{}'.format(model_cpy.summary()))
            # XXX XXX XXX

        # Reset the old model delta
        self.simplex.set_step_sizes(**orig_step_values)

    def extract(self):
        """Extract the fitted uij"""
        return self._extract(mask_datasets=None, mask_atoms=None)

    def result(self):
        """Extract the fitted parameters"""
        tls_mdls = numpy.array([p.model.values      for p in self._parameters])
        tls_amps = numpy.array([p.amplitudes.values for p in self._parameters])
        return (tls_mdls,tls_amps)

    def optimisation_summary(self, show=True):
        """Print the optimisation values and weights"""
        s = self.log._bar()+'\nOptimisation Summary: {}\n'.format(self.label)+self.log._bar()
        s += '\nOptimisation RMSD:    {:5.3f}'.format(self.optimisation_rmsd)
        s += '\nOptimisation Penalty: {:5.3f}'.format(self.optimisation_penalty)
        s += '\n'+self.log._bar()
        s += '\nModels used for optimisation (and weights):'
        for i,w in zip(self.get_dataset_mask(), self.get_dataset_weights()):
            s += '\n\t{:>3d} -> {:.3f}'.format(i,w)
        s += '\n'+self.log._bar()
        if show: self.log(s)
        return s

    def summary(self, show=True):
        """Print the number of parameters/input data"""
        s = self.log._bar()+'\nTLS Group Fit Summary: {}\n'.format(self.label)+self.log._bar()
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

    # Filter the hierarchy by the overall selection
    if params.levels.overall_selection:
        cache = filter_h.atom_selection_cache()
        choice = cache.selection(params.levels.overall_selection)
        filter_h = filter_h.select(choice)

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
        levels.append([s.strip('"') for s in default_secondary_structure_selections_filled(model.hierarchy)])
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

    # Load input pickle object
    if (params.input.pickle is None):

        if os.path.exists(params.output.out_dir):
            raise Sorry('Output directory already exists: {}\nPlease delete the directory or provide a new output directory'.format(params.output.out_dir))

        assert params.table_ones_options.column_labels
        assert params.table_ones_options.r_free_label

        easy_directory(params.output.out_dir)
        log_dir = easy_directory(os.path.join(params.output.out_dir, 'logs'))
        log = Log(os.path.join(log_dir, 'fitting.log'), verbose=params.settings.verbose)

        # Report parameters
        log.heading('Input command')
        log(' \\\n\t'.join(['pandemic.adp'] + sys.argv[1:]))
        log.heading('Input parameters')
        log(master_phil.format(params).as_str())

        log.heading('Processing input parameters')
        # Process/report on input parameters
        log.bar(True, False)
        log('Selected TLS amplitude model: {}'.format(params.fitting.tls_amplitude_model))
        log(MultiDatasetTLSModel.get_amplitude_model_class(params.fitting.tls_amplitude_model).description())
        log.bar(False, True)
        # Set precisions, etc
        log.bar()
        log('Setting TLS Model and Amplitude Precision')
        # Model precision
        log('TLS Model precision     -> {:10} / {} decimals'.format(10.**-params.fitting.precision.tls_model_decimals, params.fitting.precision.tls_model_decimals))
        TLSModel.set_precision(params.fitting.precision.tls_model_decimals)
        # Amplitude precision
        log('TLS Amplitude precision -> {:10} / {} decimals'.format(10.**-params.fitting.precision.tls_amplitude_decimals, params.fitting.precision.tls_amplitude_decimals))
        TLS_AmplitudeSet.set_precision(params.fitting.precision.tls_amplitude_decimals)
        log.bar()
        log('Setting TLS Model tolerance')
        # Model validation tolerance
        log('TLS Model tolerance     -> {:10}'.format(params.fitting.precision.tls_tolerance))
        TLSModel.set_tolerance(params.fitting.precision.tls_tolerance)
        # Model penalty tolerance
        UijPenalties.set_tls_tolerance(params.fitting.precision.tls_tolerance)
        log.bar()
        # Output images
        if params.output.images.all:
            log('params.output.images.all = {}'.format(params.output.images.all))
            # ----------------
            new = 'all'
            log('Setting params.output.images.pymol to {}'.format(new))
            params.output.images.pymol = new
            # ----------------
            new = True
            log('Setting params.output.images.distributions to {}'.format(new))
            params.output.images.distributions = new
            # ----------------
            log.bar()

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
        # Exit flags
        if params.settings.dry_run:
            log.heading('Exiting after initialisation: dry_run=True')
            sys.exit()

        #------------------------------------------------------------------------------#
        #---#                         Do the fitting                               #---#
        #------------------------------------------------------------------------------#

        # Fit TLS models
        p.fit_hierarchical_uij_model()

        # Write fitted structures
        p.write_output_structures()

        #------------------------------------------------------------------------------#
        #---#                         Model summaries                              #---#
        #------------------------------------------------------------------------------#

        # Write average uijs for each level
        p.log.heading('Summaries of the level-by-level TLS fittings')
        p.write_tls_level_summary(out_dir_tag='model')

        # Write residuals
        p.log.heading('Summaries of the residual atomic components')
        p.write_residuals_summary(out_dir_tag='model')

        # Distributions of the uijs for groups
        p.log.heading('Calculating distributions of uijs over the model')
        p.write_combined_summary_graphs(out_dir_tag='results')

        if p.params.output.diagnostics is True:
            p.run_diagnostics()

        #------------------------------------------------------------------------------#
        #---#                Refine models & calculate R-values                    #---#
        #------------------------------------------------------------------------------#

        # Refine output structures and generate table ones
        p.process_output_structures()

        #------------------------------------------------------------------------------#
        #---#                      Output CSVs and pickles                         #---#
        #------------------------------------------------------------------------------#

        # Write output CSV of... everything
        p.log.heading('Writing output csvs')
        p.write_combined_csv(out_dir_tag='root')

        # Pickle output object
        if p.params.output.pickle is True:
            log.subheading('Pickling output')
            p.pickle(pickle_file    = os.path.join(p.params.output.out_dir, 'pandemic.pickle'),
                     pickle_object  = p,
                     overwrite      = True)

        # Reduce file space
        p.tidy_output_folder(p.params.output.clean_up_files)

    else:
        log = Log()
        # Load previously-pickled parameterisation
        assert os.path.exists(params.input.pickle), 'Input pickle file does not exist: {}'.format(params.input.pickle)
        log.heading('Loading previous ADP parameterisation')
        p = libtbx.easy_pickle.load(params.input.pickle)
        assert isinstance(p, MultiDatasetUijParameterisation), 'Input pickle contains type: {}'.format(type(p))

        # Update unpickled object
        p.params = params

    # Print all errors
    p.show_warnings()

    # Write HTML
    p.log.heading('Writing HTML output')
    pandemic_html.write_adp_summary(parameterisation=p)


############################################################################

if __name__=='__main__':

    from giant.jiffies import run_default
    from pandemic import module_info
    run_default._module_info = module_info
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
        program             = PROGRAM,
        description         = DESCRIPTION)


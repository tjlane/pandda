import copy, traceback
import collections, operator
import numpy

from scitbx.array_family import flex
from scitbx import simplex, lbfgs
from libtbx.math_utils import iceil

from libtbx import adopt_init_args
from libtbx.utils import Sorry, Failure
from bamboo.common import Info, dict_from_class, logs
from bamboo.maths.functions import rms, get_function

from giant.structure.tls import tls_str_to_n_params, get_t_l_s_from_vector
from giant.structure.uij import sym_mat3_eigenvalues
from mmtbx.tls.utils import uij_eigenvalues, TLSMatricesAndAmplitudesList
from mmtbx.tls.optimise_amplitudes import MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator

from pandemic.adp.simplex import TLSSimplex, UijSimplex

from IPython import embed

############################################################################

import multiprocessing
DaemonicPool = multiprocessing.pool.Pool

def _wrapper_optimise(args):
    fitter, options = args
    stdout = options.pop('stdout', False)
    try:
        if hasattr(fitter, 'log'):
            fitter.log.toggle(stdout)
        fitter._optimise()
        return fitter
    except Exception as e:
        return traceback.format_exc()

############################################################################


class BaseOptimiser(object):


    def __init__(self,
                 target_uij,
                 convergence_tolerance=1e-3,
                 label='',
                 weights=None,
                 params=None,
                 verbose=False,
                 log=None):
        adopt_init_args(self, locals())
        #if self.log is None: log = Log()

        self._n_prm = None
        self._n_dst = None
        self._n_atm = None

        self._op = None

        self._mask_dset = None
        self._mask_atom = None

        # Uij weights
        if self.weights is None:
            self.weights = numpy.ones(target_uij.shape[:-1])
        assert self.weights.shape == target_uij.shape[:-1]

        self.optimisation_lst_sqr = numpy.inf
        self.optimisation_penalty = numpy.inf

        self.use_fitting_penalties = True


class BaseGroupOptimiser(BaseOptimiser):


    def _optimise(self):
        """Run the optimisation for the current set of variables"""

        # Prep variables for target function
        self._n_call = 0
        # Initialise the target metrics
        self.optimisation_lst_sqr = 1e6      # arbitrarily high number
        self.optimisation_penalty = 0.0
        # Apply the masks to the target uij
        self._update_target()
        # Get the variables selected for optimisatiom, and their values
        sel_vals = self._extract_values(variables=self.get_optimisation_variables())
        # Get additional input required for the simplex generation
        simplex_kw = self._simplex_input()
        # Create simplex for these variables
        opt_simplex = self.simplex.get_simplex(start=sel_vals, variables=self.get_optimisation_variables(), **simplex_kw)
        #
        if self.verbose:
            fmt = '{:>+10.5f}'
            self.log.subheading('Starting simplex')
            self.log('Starting Point:')
            self.log(','.join([fmt.format(v) for v in sel_vals]))
            self.log('Points on Simplex (relative to starting point):')
            for point in opt_simplex:
                self.log(','.join([fmt.format(v) for v in (point-sel_vals)[:30]]))
            self.log.subheading('Running optimisation')
        # Optimise these parameters
        optimised = simplex.simplex_opt(dimension = len(sel_vals),
                                        matrix    = map(flex.double, opt_simplex),
                                        evaluator = self,
                                        tolerance = self.convergence_tolerance)
        # Extract and update current values
        self._inject_values(values=optimised.get_solution(), variables=self.get_optimisation_variables())
        if self.verbose:
            self.log.bar()
            self.log('\n...optimisation complete ({} iterations)\n'.format(self._n_call))

    def _simplex_input(self):
        """Default function -- no additional variables"""
        return {}

    #===========================================+>
    # Public Functions
    #===========================================+>

    def copy(self):
        return copy.deepcopy(self)

    def n_atoms(self):
        return self._n_atm

    def n_datasets(self):
        return self._n_dst

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
        assert min(mask) >= 0
        assert max(mask) < self._n_atm
        self._mask_atom = mask

    def set_dataset_mask(self, mask):
        if not isinstance(mask, list):
            mask = list(numpy.where(mask)[0])
        assert len(mask) > 0, 'No atoms in this mask!'
        assert min(mask) >= 0
        assert max(mask) < self._n_dst
        # Store mask
        self._mask_dset = mask

    #=======================+>

    def parameters(self):
        return self._parameters

    def _msqr_delta_u(self, delta_u):
        delta_u_1d = delta_u.as_1d().as_double()
        assert delta_u_1d.size() == self._op.wgts_6.size()
        # Normalise by the NUMBER OF DATASETS * NUMBER OF ATOMS * 3 (number of eigenvalues of U)
        msqr_delta_u = flex.sum(delta_u_1d*delta_u_1d*self._op.wgts_6)/(self._op.n_dst*self._op.n_atm)
        return msqr_delta_u

    def _mean_delta_u(self, delta_u):
        delta_u_eigs_1d = uij_eigenvalues(delta_u).as_1d().as_double()
        assert delta_u_eigs_1d.size() == self._op.wgts_3.size()
        # Normalise by the NUMBER OF DATASETS * NUMBER OF ATOMS * 3 (number of eigenvalues of U)
        mean_delta_u = abs(flex.sum(delta_u_eigs_1d*self._op.wgts_3)/(self._op.n_dst*self._op.n_atm))
        return mean_delta_u

    def _mabs_delta_u(self, delta_u):
        delta_u_1d = delta_u.as_1d().as_double()
        assert delta_u_1d.size() == self._op.wgts_6.size()
        # Normalise by the NUMBER OF DATASETS * NUMBER OF ATOMS * 3 (number of eigenvalues of U)
        mabs_delta_u = flex.sum(flex.abs(delta_u_1d)*self._op.wgts_6)/(self._op.n_dst*self._op.n_atm)
        return mabs_delta_u

    def target(self, values):
        """Target function for the simplex optimisation"""

        # Increment counter
        self._n_call += 1

        # Print info line if necessary
        if self.verbose:
            if self._n_call%20==1:
                header = '[{}] -> ({:^12}, {:^10})'.format(', '.join(['{:>10}'.format('param') for r in values]), 'target', 'penalty')
                line = '-'*len(header)
                self.log(line)
                self.log(header)
                self.log(line)

        # Combine the optimising values in the complete parameter set
        self._inject_values(values=values)

        # Calculate physical penalties - reject this set if model is not physical
        p_valid = self.is_valid()
        # Return now if physical penalty if non-zero to save time
        if p_valid is False:
            if self.verbose:
                self.log('[{}] -> ({:>12}, {:>10})'.format(', '.join(['{:+10.7f}'.format(r) for r in values]), 'UNPHYSICAL', ''))
            # Ensure that these parameters always have a higher values than the current best solution
            return 2.0*(self.optimisation_lst_sqr+self.optimisation_penalty)

        # Get the fitted uijs (including masked atoms)
        self._update_fitted()

        # Calculate weighted sum of squares
        # !!!IMPORTANT!!!
        # MUST be (fitted-target) so that eigenvalues of delta_u are positive when fitted greater than target!
        # !!!IMPORTANT!!!
        delta_u = self._fitted_uij - self._target_uij
        msqr_delta_u = self._msqr_delta_u(delta_u=delta_u)
        target = msqr_delta_u

        # Calculate fitting penalties (add to target)
        if self.use_fitting_penalties:
            penalty = self.additional_target_terms(delta_u=delta_u, msqr_delta_u=msqr_delta_u)
        else:
            penalty = 0.0

        # Total target function value
        tot_val = target + penalty
        # Update minima
        if tot_val < self.optimisation_lst_sqr+self.optimisation_penalty:
            self.optimisation_lst_sqr = target
            self.optimisation_penalty = penalty
        if self.verbose:
            self.log('[{}] -> ({:12.9f}, {:10.7f})'.format(', '.join(['{:+10.7f}'.format(r) for r in values]), target, penalty))
        return tot_val


class MultiDatasetTLSFitter(BaseGroupOptimiser):


    _amp_tolerance = 1e-12

    def __init__(self,
                 target_uij,
                 atomic_xyz,
                 atomic_com,
                 n_tls=1,
                 label='',
                 weights=None,
                 params=None,
                 verbose=False,
                 log=None):
        super(MultiDatasetTLSFitter, self).__init__(
                target_uij=target_uij,
                convergence_tolerance=params.simplex_convergence,
                label=label,
                weights=weights,
                params=params,
                verbose=verbose,
                log=log)

        # Store the centre of mass of the atoms (for the rotation/screw components)
        self.atomic_xyz = atomic_xyz
        self.atomic_com = atomic_com

        # Should be n_dataset observations of n_atm with 6 parameters
        assert len(self.target_uij.shape) == 3
        assert self.target_uij.shape[-1] == 6
        assert self.atomic_xyz is not None
        assert self.atomic_com is not None

        self._n_dst = self.target_uij.shape[0]
        self._n_atm = self.target_uij.shape[1]

        assert self.target_uij.shape == (self._n_dst, self._n_atm, 6)
        assert self.atomic_xyz.shape == (self._n_dst, self._n_atm, 3)
        assert self.atomic_com.shape == (self._n_dst, 3)

        # Number of models and parameters
        self._n_tls = n_tls

        # Model-Amplitude parameter sets
        self._parameters = TLSMatricesAndAmplitudesList(length = self._n_tls,
                                                        n_amplitudes = self._n_dst)
        # XXX Set precisions here? XXX

        # Extract number of parameters
        self._n_prm = self._parameters.n_params()

        # Initialise the masks
        self.set_atomic_mask(range(self._n_atm))
        self.set_dataset_mask(range(self._n_dst))

        # Optimisation ranges
        self.set_variables(optimise_matrices = True,
                           optimise_amplitudes = True,
                           components = 'TLS',
                           modes = range(self._n_tls))

        # Initialise simplex generator
        self.simplex = TLSSimplex(
                vibration_delta = params.simplex.vibration_delta,
                libration_delta = params.simplex.libration_delta,
                amplitude_delta_min = params.simplex.amplitude_delta_min,
                amplitude_delta_frac = params.simplex.amplitude_delta_frac,
                )

        # Fitting penalty functions
        self._model_input_penalty_function = get_function(**dict_from_class(params.penalties.over_target_values))

    #===========================================+>
    # Private Functions - common to parent class
    #===========================================+>

    def additional_target_terms(self, delta_u, msqr_delta_u):
        """Calculate penalties to enforce priors about fitting (e.g. model should not be greater than input)"""

        # Calculate the average delta_u
        mabs_delta_u = self._mean_delta_u(delta_u)
        # Get penalty for each dataset
        penalties = self._model_input_penalty_function(uij_eigenvalues(delta_u).as_1d().as_double())
        assert penalties.size == self._op.wgts_3.size()

        # Average penalties, with weights for each dataset
        total_penalty =  mabs_delta_u * numpy.sum(penalties) #/ 3.0 # Divide by the number of dimensions

        # Normalise by the NUMBER OF DATASETS
        total_penalty /= (self._op.n_dst * self._op.n_atm)

        return total_penalty

    def is_valid(self):
        """Calculate penalties to enforce that the parameters are physical (e.g. positive amplitudes)"""
        penalties = []
        for mode in self._parameters:
            if self._op.optimise_matrices:
                if not mode.matrices.is_valid():
                    return False
            if self._op.optimise_amplitudes:
                values = mode.amplitudes.get()
                if abs(flex.sum(values.select(values<0.0))):
                    return False
        return True

    #=======================+>

    def _update_fitted(self):
        op = self._op
        uij = self._parameters.uijs(sites_carts=op.xyzs, origins=op.coms, selection=op.i_dst)
        assert uij.all() == (op.n_dst, op.n_atm,)
        self._fitted_uij = uij

    def _update_target(self):
        op = self._op
        uij = self.target_uij[op.i_dst][:,op.i_atm]
        uij = flex.sym_mat3_double(uij.reshape((op.n_dst*op.n_atm,6)))
        uij.reshape(flex.grid((op.n_dst,op.n_atm)))
        assert uij.all() == (op.n_dst, op.n_atm,)
        self._target_uij = uij

    #=======================+>

    def get_optimisation_variables(self):
        return self._op

    def set_variables(self, optimise_matrices=False, optimise_amplitudes=False, components=None, modes=None):
        """Set optimisation variables:
            optimise_matrices: True or False
                should TLS matrix components be optimised?
            optimise_amplitudes: True or False
                should TLS amplitudes be optimised?
            components: string
                which matrices/amplitudes elements should be optimised.
                can only contain letters T or L or S
            modes: list of integers
                which TLS modes should be optmised?
        """

        assert isinstance(optimise_matrices, bool)
        assert isinstance(optimise_amplitudes, bool)

        if components is not None:
            assert not set(''.join(components) if isinstance(components, list) else components).difference('TLS')
        else:
            components = 'TLS'

        if modes is not None:
            assert min(modes) >= 0
            assert max(modes) < self._n_tls
        else:
            modes = range(self._n_tls)

        # Logic checks
        if optimise_matrices is True:
            assert components is not None
        if optimise_amplitudes is True:
            assert modes is not None

        # Extract masks
        dst_mask = self.get_dataset_mask()
        atm_mask = self.get_atomic_mask()

        # Extract masked coordinates
        xyzs = self.atomic_xyz[dst_mask][:,atm_mask]
        coms = self.atomic_com[dst_mask]
        wgts = self.weights[dst_mask][:,atm_mask]
        # Renormalise weights
        wgts = wgts / wgts.mean()
        # Extract sizes
        n_dst = len(dst_mask)
        n_atm = len(atm_mask)
        # Check shapes
        assert xyzs.shape == (n_dst, n_atm, 3)
        assert coms.shape == (n_dst, 3)
        assert wgts.shape == (n_dst, n_atm)

        # Convert to flex
        xyzs = flex.vec3_double(xyzs.reshape(n_dst*n_atm,3))
        xyzs.reshape(flex.grid(n_dst,n_atm))
        coms = flex.vec3_double(coms)

        # Do with wgts as n_dst*n_atm (but as 1d array for speed!)
        wgts_1 = flex.double(wgts.reshape((n_dst*n_atm,)).tolist())

        # Do with wgts as n_dst*n_atm*3 (but as 1d array for speed!)
        wgts_3 = numpy.repeat(wgts.reshape(wgts.shape+(1,)), repeats=3, axis=2)
        wgts_3 = flex.double(wgts_3.reshape((n_dst*n_atm*3,)).tolist())

        # Do with wgts as n_dst*n_atm*6 (but as 1d array for speed!)
        wgts_6 = numpy.repeat(wgts.reshape(wgts.shape+(1,)), repeats=6, axis=2)
        wgts_6 = flex.double(wgts_6.reshape((n_dst*n_atm*6,)).tolist())

        # Make sure all weights are normalised!
        #wgts_1 = wgts_1 / flex.mean(wgts_1)
        #wgts_3 = wgts_3 / flex.mean(wgts_3)
        #wgts_6 = wgts_6 / flex.mean(wgts_6)
        assert abs(flex.mean(wgts_1)-1.0) < 1e-3, 'weights should be normalised!'
        assert abs(flex.mean(wgts_3)-1.0) < 1e-3, 'weights should be normalised!'
        assert abs(flex.mean(wgts_6)-1.0) < 1e-3, 'weights should be normalised!'

        # Create a dictionary that FULLY defines the optimisation cycle
        self._op = Info({
            'optimise_matrices'   : optimise_matrices,
            'optimise_amplitudes' : optimise_amplitudes,
            'components'          : components,
            # Lists of indices
            'i_mdl' : modes,
            'i_dst' : flex.size_t(dst_mask),
            'i_atm' : flex.size_t(atm_mask),
            # Lengths of above lists
            'n_mdl' : len(modes),
            'n_dst' : n_dst,
            'n_atm' : n_atm,
            # Masked coordinates
            'xyzs' : xyzs,
            'coms' : coms,
            # Masked Uij weights
            'wgts'   : wgts,
            'wgts_1' : wgts_1,
            'wgts_3' : wgts_3,
            'wgts_6' : wgts_6,
            })

        return self._op

    def _extract_values(self, variables=None):
        """Extract a set of values for a set of variables"""

        # Use default selection if not provided
        if variables is None:
            variables = self._op
        # Create object from dict
        s = variables
        # Iterate through objects and extract requested parameter values
        values = []
        if s.optimise_matrices:
            for i in s.i_mdl:
                mdl = self._parameters.get(index=i)
                values.append(mdl.matrices.get(component_string=s.components, include_szz=False))
        if s.optimise_amplitudes:
            for i in s.i_mdl:
                mdl = self._parameters.get(index=i)
                values.append(mdl.amplitudes.get(selection=s.i_dst))
        return numpy.concatenate(values)

    def _inject_values(self, values, variables=None):
        """Change a set of values for a set of variables"""

        # Use default selection if not provided
        if variables is None:
            variables = self._op
        # Create object from dict
        s = variables
        # Iterate through objects and extract requested parameter values
        values = collections.deque(values)
        if s.optimise_matrices:
            n_tls_params = tls_str_to_n_params(s.components, include_szz=False)
            for i in s.i_mdl:
                mdl = self._parameters.get(index=i)
                mdl.matrices.set(values = [values.popleft() for _ in xrange(n_tls_params)],
                                 component_string = s.components,
                                 include_szz = False)
        if s.optimise_amplitudes:
            for i in s.i_mdl:
                mdl = self._parameters.get(index=i)
                mdl.amplitudes.set(values = [values.popleft() for _ in xrange(s.n_dst)],
                                   selection = s.i_dst)
                mdl.amplitudes.zero_negative_values()
        assert len(values) == 0, 'not all values have been popped'

    #===========================================+>

    def _simplex_input(self):
        """Return the current values of the TLS matrices"""
        tls_matrices = [m.matrices for m in self.parameters()]
        tls_amplitudes = [m.amplitudes.values for m in self.parameters()]
        return {'tls_matrices':tls_matrices, 'tls_amplitudes':tls_amplitudes}

    #===========================================+>
    # Public Functions
    #===========================================+>

    def n_params(self, non_zero=False):
        return self._parameters.n_params(free=True, non_zero=non_zero)

    def initialise_matrices_and_amplitudes(self, i_tls, dataset_mask, atomic_mask, minimum_amplitude=1e-6, n_cpus=1):
        """Set the TLS matrices & amplitudes to sensible starting values"""

        n_dst, n_atm = len(dataset_mask), len(atomic_mask)

        self.log.bar(True,True)
        self.log('Resetting mode {}'.format(i_tls+1))
        self.parameters().get(index=i_tls).reset()

        self.log('Looking for an appropriate scale for the amplitudes')
        # Get the unfitted Uij (target uij minus uij already fitted from other modes)
        uij_vals = (self.target_uij-self.extract())[dataset_mask][:,atomic_mask]
        assert uij_vals.shape == (n_dst, n_atm, 6)
        uij_vals_flex = flex.sym_mat3_double(uij_vals.reshape((n_dst*n_atm,6)))
        uij_eigs_flex = uij_eigenvalues(uij_vals_flex)
        uij_eigs = numpy.array(uij_eigs_flex).reshape((n_dst,n_atm,3))
        assert uij_eigs.shape == (n_dst, n_atm, 3)
        # Calculate the minimum eigenvalues for each dataset
        uij_eigs_min = uij_eigs.min(axis=(1,2))
        assert uij_eigs_min.shape == (len(dataset_mask),)
        min_eig = numpy.min(uij_eigs_min)
        # Set minimum eigenvalues as the amplitudes
        amp_scales = numpy.zeros(self._n_dst)
        amp_scales[dataset_mask] = uij_eigs_min
        amp_scales[amp_scales<0.0] = 0.0
        if amp_scales.any():
            self.log('Scaling amplitudes by minimum eigenvalues of input Uijs')
            self.log('Amplitude range: {} - {}'.format(amp_scales.min(), amp_scales.max()))
        else:
            self.log('No minimal eigenvalues of input Uijs are positive.')
            self.log('Scaling amplitudes to a fixed minimal size ({})'.format(minimum_amplitude))
            amp_scales[dataset_mask] = minimum_amplitude
        self.log('Setting amplitudes')
        mode = self.parameters().get(index=i_tls)
        mode.amplitudes.set(values=amp_scales.tolist())

        self.log('Initialising the T-matrix with an isotropic uij')
        uij_start = (1.0,)*3 + (0.0,)*3
        self.log('Setting initial T-matrix values to {}'.format(uij_start))
        self.parameters().get(index=i_tls).matrices.set(values=uij_start, component_string='T')
        self.log('Starting TLS-matrix values:')
        self.log(self.parameters().get(index=i_tls).matrices.summary())

        self.log('Optimising Amplitudes')
        self.optimise_amplitudes(
            modes  = [i_tls],
            n_cpus = n_cpus)
        self.log('Optimised TLS-amplitude values:')
        self.log(self.parameters().get(index=i_tls).amplitudes.summary())

        self.log('Setting zero value amplitudes to minimum value')
        amps = self.parameters().get(index=i_tls).amplitudes.get()
        amps.set_selected(amps<minimum_amplitude, minimum_amplitude)
        self.parameters().get(index=i_tls).amplitudes.set(amps)

        self.log.bar()
        self.log('Starting amplitude values:')
        self.log(self.parameters().get(index=i_tls).amplitudes.summary())

        #self.log.bar()
        #self.log('Resetting Matrices to zero')
        #self.log(self.parameters().get(index=i_tls).matrices.reset())
        #self.log('Setting L Matrices:')
        #self.parameters().get(index=i_tls).matrices.set(values=uij_start, component_string='L')

        self.log.bar()
        self.log('Starting values:')
        self.log(self.parameters().get(index=i_tls).summary())

        return self.parameters().get(index=i_tls)

    def optimise(self, n_cycles=1, n_cpus=1):
        """Optimise a (series of) TLS model(s) against the target data"""

        # Extract the masks (so that can be reapplied if changed)
        opt_dset_mask = self.get_dataset_mask()
        opt_atom_mask = self.get_atomic_mask()

        # Cumulative i_tls indices
        mdl_cuml = []

        # Did this optimisation start from zero-valued model?
        global_null = False #(not self.parameters().any())

        # Optimise!
        for i_cycle in xrange(n_cycles):
            self.log.subheading('Group {} - Optimisation cycle {} of {}'.format(self.label, i_cycle+1, n_cycles))

            # Optimise one TLS "mode" at time
            for i_tls in xrange(self._n_tls):

                # How should T-L-S be optimised - if all values are zero, optimise all separately then as those for each amplitude/component set
                all_cpts = "TLS"
                # Are all matrix values zero?
                null_mode = not self.parameters().get(index=i_tls).matrices.any(component_string=all_cpts)
                # If null mode then optimise each matrix individually (for speed) and then optimise all together
                #cpt_groups = list(all_cpts)*null_mode + [all_cpts]*((not null_mode) or (len(all_cpts) > 1))
                cpt_groups = ["TLS"]

                self.log.subheading('Optimising TLS mode {} of {}'.format(i_tls+1, self._n_tls))
                self.log('Optimising using {} atoms'.format(len(opt_atom_mask)))
                self.log('Optimising using {} datasets'.format(len(opt_dset_mask)))
                self.log.bar(True, True)
                self.log('Optimisation order for TLS matrix values: \n\t{}'.format('\n\t'.join(['({}) {}'.format(i+1, c) for i, c in enumerate(cpt_groups)])))

                # If starting from zero values, set translation matrix to sensible starting values
                if (null_mode is True):
                    self.initialise_matrices_and_amplitudes(
                            i_tls=i_tls,
                            dataset_mask=range(self._n_dst),
                            atomic_mask=opt_atom_mask,
                            minimum_amplitude=self.params.simplex.amplitude_delta_min,
                            n_cpus=n_cpus,
                            )
                    self.log(self.parameters().get(index=i_tls).summary())

                # Append to running list
                if i_tls not in mdl_cuml: mdl_cuml.append(i_tls)

                # ---------------------------------->
                # Optimise each set of selected components separately
                # ---------------------------------->
                for cpts in cpt_groups:
                    # Skip S optimisation if T and L are zeros
                    if (cpts == 'S') and not (self.parameters().get(index=i_tls).matrices.any('T') and \
                                              self.parameters().get(index=i_tls).matrices.any('L')):
                        self.log('T and L matrices are zero -- not optimising S-matrix')
                        continue
                    self.log.subheading('Optimising {} parameters for TLS mode {}'.format(cpts, i_tls+1))
                    self.log.bar(False, True)
                    # Run optimisation
                    self.optimise_matrices(
                        modes       = [i_tls],
                        components  = cpts)
                    # Log model summary
                    self.log(self.parameters().get(index=i_tls).matrices.summary()+'\n')
                    # Log current target and penalty
                    self.optimisation_summary()
                    # Apply multiplier to amplitudes of less than unity if null mode to allow for "wiggling"
                    if null_mode and (cpts in list('TLS')):
                        multiplier = 0.95
                        self.log.subheading('Scaling amplitudes after individual component optimisation')
                        self.log('Applying scaling to amplitudes of {} to reduce restrictions on subsequent simplex optimisations'.format(multiplier))
                        self.parameters().get(index=i_tls).amplitudes.multiply(scalar=multiplier)

                # If TLS Matrices optimise to negligible values, move on to next
                if not self.parameters().get(index=i_tls).matrices.any(all_cpts, tolerance=self.params.normalised_tls_eps):
                    self.log('{} matrices are all zero -- not optimising amplitudes'.format(all_cpts))
                    self.parameters().get(index=i_tls).reset()
                    # Set to isotropic so can still be used in inter-level optimisation
                    if i_tls == 0:
                        self.parameters().get(index=i_tls).matrices.set((1,1,1,0,0,0), "T")
                        self.parameters().get(index=i_tls).amplitudes.zero_values()
                    continue

                # Normalise matrices to give Uijs of approximately xA
                xyzs = flex.vec3_double(self.atomic_xyz.reshape(self._n_dst*self._n_atm,3))
                xyzs.reshape(flex.grid(self._n_dst,self._n_atm))
                coms = flex.vec3_double(self.atomic_com)
                self.parameters().get(index=i_tls).normalise_by_matrices(sites_carts=xyzs, origins=coms, target=self.params.normalised_tls_scale)
                self.optimise_invalid_modes() # Check that normalisation has not made any of the matrices invalid and correct if possible

                # If TLS Matrices optimise to negligible values, move on to next
                if not self.parameters().get(index=i_tls).matrices.any(all_cpts, tolerance=self.params.normalised_tls_eps):
                    self.log('{} matrices are all zero -- not optimising amplitudes'.format(all_cpts))
                    self.parameters().get(index=i_tls).reset()
                    # Set to isotropic so can still be used in inter-level optimisation
                    if i_tls == 0:
                        self.parameters().get(index=i_tls).matrices.set((1,1,1,0,0,0), "T")
                        self.parameters().get(index=i_tls).amplitudes.zero_values()
                    continue

                # ---------------------------------->
                # Optimise TLS amplitude parameters (all amplitudes!)
                # ---------------------------------->
                # Report
                self.log.subheading('Optimising TLS amplitudes for all datasets')
                # Reset amplitudes to zero to prevent overparameterisation
                #self.log('Resetting all amplitudes to zero')
                #self.parameters().zero_amplitudes(models=[i_tls])
                # Run optimisation
                self.optimise_amplitudes(
                    modes  = mdl_cuml,
                    n_cpus = n_cpus)
                # Log amplitudes summary
                self.log.bar(blank_before=True, blank_after=True)
                self.log(self.parameters().get(index=i_tls).amplitudes.summary())

            # Reapply atom and dataset masks
            self.set_dataset_mask(opt_dset_mask)
            self.set_atomic_mask(opt_atom_mask)

            # End of cycle house-keeping
            self.log.bar(True, False)
            # Break optimisation if all matrices are zero -- not worth running following optimisation cycles
            if self.parameters().is_null():
                self.log('All matrices have refined to zero -- optimisation finished.')
                # Reset parameters just to be sure
                #self.parameters().reset()
                self.log.bar()
                break
            # Check to see if any amplitudes are negative
            self.log('Looking for negative TLS amplitudes')
            self.parameters().zero_negative_amplitudes()
            self.log.bar()
            ## Identify modes that are zero and reset these
            #self.log('Resetting zero-value modes')
            #self.parameters().reset_null_modes()
            #self.log.bar()

            # Show summary at end of cycle...
            self.log.subheading('End-of-cycle summary')
            self.summary(show=True)
            self.log('')

        self.log.subheading('End of optimisation summary')
        self.summary(show=True)

        return

    def optimise_matrices(self, modes, components):
        """Optimise the matrices for combinations of modes/componenets/datasets"""

        # Select variables for optimisation -- matrices only
        self.set_variables(optimise_matrices   = True,
                           optimise_amplitudes = False,
                           modes      = modes,
                           components = components)
        # Run optimisation
        self._optimise()

    def optimise_amplitudes(self, modes, n_cpus=1):
        """Optimise the amplitudes for combinations of modes"""

        # Extract masks
        original_d_mask = self.get_dataset_mask()
        original_a_mask = self.get_atomic_mask()

        # Select all atoms
        self.set_atomic_mask(range(self._n_atm))

        self.log('Running amplitude optimisations: --->')
        # Optimise all amplitudes dataset-by-dataset
        jobs = []
        for i_dst in xrange(self._n_dst):
            # Select this dataset
            self.set_dataset_mask([i_dst])
            self.set_variables(optimise_matrices   = False,
                               optimise_amplitudes = True,
                               modes      = modes)
            # Append to job list (create copy) or run optimisation
            if n_cpus > 1:
                jobs.append((self.copy(), {'stdout':False}))
            else:
                self._optimise()
                self.log('> dataset {} of {} (target {:.3f}; penalty {:.1f})'.format(i_dst+1, self._n_dst, self.optimisation_lst_sqr, self.optimisation_penalty))
        # Run parallel jobs and inject results
        if n_cpus > 1:
            # Workers for multiprocessing
            workers = DaemonicPool(n_cpus)
            # Report and map to workers
            self.log.subheading('Running {} jobs using {} cpus'.format(len(jobs), min(n_cpus,len(jobs))))
            finished_jobs = workers.map(func=_wrapper_optimise, iterable=jobs)
            # Inject results from jobs
            self.log('Optimisation Results:')
            for i_dst in xrange(self._n_dst):
                job = finished_jobs[i_dst]
                if isinstance(job, str):
                    self.log(job)
                    raise Failure('error returned')
                self._inject_values(values=job._extract_values(), variables=job.get_optimisation_variables())
                self.log('> dataset {} of {} (target {:.3f}; penalty {:.1f})'.format(i_dst+1, self._n_dst, job.optimisation_lst_sqr,  job.optimisation_penalty))
            # Let's be good
            workers.close()

        self.set_dataset_mask(original_d_mask)
        self.set_atomic_mask(original_a_mask)

    def optimise_invalid_modes(self, n_cpus=1):
        """Check if models are invalid and attempt to correct them"""

        self.log.subheading('Checking for invalid TLS matrices')

        for i_tls, mode in enumerate(self.parameters()):

            self.log.bar()
            self.log('Checking physicality of TLS mode {}'.format(i_tls+1))
            self.log.bar()

            # Check that the core matrices are valid
            if mode.matrices.is_valid():
                self.log('Core TLS matrices are valid')
                continue

            # TODO Check if it has almost zero-values and just reset

            # Mode invalid -- try to fix
            self.log('Core matrices are invalid\n...attempting to re-optimise:')
            # Create copy of the matrices first
            mode_cpy = mode.copy()
            self.log('Before re-optimisation: \n{}'.format(mode_cpy.matrices.summary()))

            # Try resetting S if T or L is zero
            if not (mode.matrices.any("T") and mode.matrices.any("L")):
                self.log('...T or L matrix are zero -- resetting S-matrix')
                mode.matrices.set(values=[0.0]*9, component_string='S')
                if mode.matrices.is_valid():
                    self.log('...core TLS matrices are now valid')
                    continue

            # Try resetting and optimising the TLS Matrices
            # Reset to zero and re-optimise
            self.log('...core matrices are still invalid: \n{}'.format(mode.matrices.summary()))
            self.log('...resetting and redetermining TLS values')
            # Reset TLS matrices to zero
            mode.matrices.reset()
            # Iterate through and optimise as originally
            self.optimise_matrices(modes      = [i_tls],
                                   components = 'TLS')
            mode = self.parameters().get(index=i_tls)
            if mode.matrices.is_valid():
                self.log('...core matrices now valid (after resetting and reoptimising)')
                self.log('Before re-optimisation: \n{}'.format(mode_cpy.matrices.summary()))
                self.log('After re-optimisation: \n{}'.format(mode.matrices.summary()))
                rms_values = rms(mode.matrices.get()-mode_cpy.matrices.get())
                self.log('...rms between initial matrices and "fixed" matrices: {}'.format(rms_values))
                continue

            raise Failure('Failed to fix core model. \n{}'.format(mode_cpy.summary()))

    def extract(self, sum_modes=True):
        """Extract the fitted uij"""
        # Get the atomic coordinates
        xyzs = self.atomic_xyz
        coms = self.atomic_com
        # Apply masks - datasets
        n_dst = self._n_dst
        n_atm = self._n_atm
        n_tls = self._n_tls
        # Validate
        assert xyzs.shape == (n_dst, n_atm, 3)
        # Convert arrays to flex
        xyzs = flex.vec3_double(xyzs.reshape(n_dst*n_atm,3))
        xyzs.reshape(flex.grid(n_dst,n_atm))
        coms = flex.vec3_double(coms)
        # Calculate the TLS components
        if sum_modes is True:
            uij = self.parameters().uijs(sites_carts=xyzs, origins=coms)
            assert uij.all() == (n_dst, n_atm,)
            uij = uij.as_1d().as_double().as_numpy_array().reshape((n_dst,n_atm,6))
        else:
            uij = [m.uijs(sites_carts=xyzs, origins=coms) for m in self.parameters()]
            assert uij[0].all() == (n_dst, n_atm,)
            uij = numpy.array([u.as_1d().as_double().as_numpy_array().reshape((n_dst,n_atm,6)) for u in uij])
            assert uij.shape == (n_tls, n_dst, n_atm, 6)
        return uij

    def extract_core(self):
        """Extract the fitted uij (but not multiplied by the dataset amplitudes)"""
        # Get the atomic coordinates
        xyzs = self.atomic_xyz
        coms = self.atomic_com
        # Apply masks - datasets
        n_dst = self._n_dst
        n_atm = self._n_atm
        n_tls = self._n_tls
        # Validate
        assert xyzs.shape == (n_dst, n_atm, 3)
        # Calculate the TLS components
        uij = [[m.matrices.uijs(sites_cart=flex.vec3_double(xyzs[i]), origin=tuple(coms[i])).as_1d().as_double().as_numpy_array() for i in xrange(n_dst)] for m in self.parameters()]
        uij = numpy.array(uij).reshape((n_tls, n_dst, n_atm, 6))
        assert uij.shape == (n_tls, n_dst, n_atm, 6)
        return uij

    def result(self):
        """Extract the fitted parameters"""
        tls_mdls = numpy.array([p.matrices.get()   for p in self.parameters()])
        tls_amps = numpy.array([p.amplitudes.get() for p in self.parameters()])
        return (tls_mdls,tls_amps)

    def optimisation_summary(self, show=True):
        """Print the optimisation values and weights"""
        s = self.log._bar()+'\nOptimisation Summary: {}\n'.format(self.label)+self.log._bar()
        s += '\nOptimisation Least Squares: {:5.3f}'.format(self.optimisation_lst_sqr)
        s += '\nOptimisation Penalty: {:5.3f}'.format(self.optimisation_penalty)
        s += '\nOptimisation Target: {:5.3f}'.format(self.optimisation_lst_sqr+self.optimisation_penalty)
        s += '\n'+self.log._bar()
        s += '\nModels used for optimisation (and total weights):'
        wts = self._op.wgts.mean(axis=1)
        assert len(wts) == len(self.get_dataset_mask())
        for i_rel, i_abs in enumerate(self.get_dataset_mask()):
            s += '\n\t{:>3d} -> {:.3f}'.format(i_abs,wts[i_rel])
        s += '\n'+self.log._bar()
        if show: self.log(s)
        return s

    def summary(self, show=True):
        """Print the number of parameters/input data"""
        s = self.log._bar()+'\nTLS Group Fit Summary: {}\n'.format(self.label)+self.log._bar()
        for i_tls in xrange(self._n_tls):
            s += '\n> TLS model {}'.format(i_tls+1)
            mode = self.parameters().get(index=i_tls)
            s += '\n\t' + mode.matrices.summary().replace('\n','\n\t')
            s += '\n\t' + mode.amplitudes.summary().replace('\n','\n\t')
        if show: self.log(s)
        return s


class MultiDatasetUijFitter(BaseGroupOptimiser):


    _uij_tolerance = 1e-6

    def __init__(self,
                 target_uij,
                 label='',
                 weights=None,
                 params=None,
                 verbose=False,
                 log=None):
        super(MultiDatasetUijFitter, self).__init__(
                target_uij=target_uij,
                convergence_tolerance=params.simplex_convergence,
                label=label,
                weights=weights,
                params=params,
                verbose=verbose,
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
        self.simplex = UijSimplex(uij_delta=params.simplex.uij_delta)

        # Initialise the mask to all datasets
        self.set_dataset_mask(range(self._n_dst))

    #===========================================+>
    # Private Functions - common to parent class
    #===========================================+>

    def additional_target_terms(self, delta_u, msqr_delta_u):
        return 0.0

    def is_valid(self):
        return self.uij_is_valid(self.result())

    def uij_is_valid(self, u):
        eigenvalues = sym_mat3_eigenvalues(u)
        penalty = abs(flex.sum(eigenvalues.select(eigenvalues<(-1.0*self._uij_tolerance))))
        if penalty > 0.0:
            return False
        return True

    #=======================+>

    def _update_fitted(self):
        self._fitted_uij = flex.sym_mat3_double([self.extract()]*self._op.n_dst)

    def _update_target(self):
        self._target_uij = flex.sym_mat3_double(self.target_uij[self._op.i_dst])

    #=======================+>

    def get_optimisation_variables(self):
        return self._op

    def set_variables(self, mode='all'):
        """Define which variables are to be optimised and return selection dictionary"""
        # Check mode
        assert mode in ['all', 'multiplier']
        # Extract mask and mask length
        dst_mask = self.get_dataset_mask()
        n_dst = len(dst_mask)
        # Extract weights
        wgts = self.weights[dst_mask]
        # Renormalise the extrated weights
        wgts = wgts / wgts.mean()
        # Do with wgts as n_dst*6
        wgts_6 = numpy.repeat(wgts.reshape(wgts.shape+(1,)), repeats=6, axis=1)
        wgts_6 = flex.double(wgts_6.reshape((n_dst*6,)).tolist())
        # Make sure all weights are normalised!
        assert abs(flex.mean(wgts_6)-1.0) < 1e-3, 'weights should be normalised!'

        self._op = Info({
            # Mode of optimisation
            'mode'  : mode,
            'start' : self.parameters().copy(),
            # Lists of indices
            'i_dst' : flex.size_t(dst_mask),
            # Lengths of above lists
            'n_dst' : n_dst,
            'n_atm' : 1,
            # Masked coordinates
            'wgts'   : wgts,
            'wgts_6' : wgts_6,
            })
        return self._op

    def _extract_values(self, variables=None):
        # Use default selection if not provided
        if variables is None:
            variables = self._op
        # Either get values, or divide by start values
        if variables.mode == 'all':
            return self._parameters
        elif variables.mode == 'multiplier':
            return (self._parameters / variables.start)[0]
        else:
            raise Exception('Invalid mode: {}'.format(variables.mode))

    def _inject_values(self, values, variables=None):
        """Insert a set of values into the complete parameter set"""
        # Use default selection if not provided
        if variables is None:
            variables = self._op
        # Either set values, or multiply from start values
        if variables.mode == 'all':
            assert len(values) == self._n_prm
            self._parameters[:] = values
        elif variables.mode == 'multiplier':
            assert len(values) == 1
            self._parameters[:] = (values[0] * variables.start)
        else:
            raise Exception('Invalid mode: {}'.format(variables.mode))

    #===========================================+>
    # Public Functions - common to parent class
    #===========================================+>

    def n_params(self, non_zero=False):
        if non_zero is False:
            return numpy.product(self.parameters().shape)
        else:
            return (self.parameters()!=0.0).any()*(self.n_params(non_zero=False))

    def optimise(self, n_cycles=None, n_cpus=None): # Options not used in optimisation
        """Optimise the residual for a set of atoms"""
        # Initialise/update optimisation info
        self.set_variables(mode='all')
        # If only 1 dataset provided, set current values to target values and return
        if (self.target_uij.shape[0] == 1) and self.uij_is_valid(self.target_uij[0]):
            self._inject_values(values=self.target_uij[0])
            return
        self._optimise()

    def result(self):
        """Return the fitted parameters (same as extract for this class)"""
        return tuple(self.parameters())

    def extract(self):
        """Return the fitted uijs - for all atoms"""
        return tuple(self.parameters())

    def summary(self, show=True):
        """Print the number of parameters/input data"""
        uij = self.result()
        s = 'Uij ({}): '.format(self.label)+', '.join(['{:8.3f}'.format(v) for v in uij])
        if show: self.log(s)
        return s

    #===========================================+>
    # Penalty functions - unique to this class
    #===========================================+>

    @classmethod
    def set_uij_tolerance(cls, tolerance):
        cls._uij_tolerance = tolerance


class OptimiseAmplitudes:


    def __init__(self,
            target_uijs,
            target_weights,
            base_amplitudes,
            base_uijs,
            base_atom_indices,
            dataset_hash,
            residual_uijs,
            residual_amplitudes=None,
            residual_mask=None,
            convergence_tolerance=1e-5):

        # Compatibility of parameters is checked in the optimiser so not checking here

        # Mask for residual optimisation
        if (residual_mask is not None):
            residual_mask = flex.bool(map(bool, residual_mask))

        self.target_uijs        = target_uijs
        self.target_weights     = target_weights
        self.base_amplitudes    = base_amplitudes
        self.base_uijs          = base_uijs
        self.base_atom_indices  = base_atom_indices
        self.dataset_hash       = dataset_hash
        self.residual_uijs      = residual_uijs
        self.residual_amplitudes = residual_amplitudes
        self.residual_mask      = residual_mask
        self.convergence_tolerance = convergence_tolerance

        self.initial = None
        self.result = None

    def run(self):
        self._optimise()
        return self

    def _optimise(self):

        self.log = logs.Log()

        # Setup
        self.f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
                self.target_uijs,
                self.base_amplitudes,
                self.base_uijs,
                self.base_atom_indices,
                self.dataset_hash,
                self.residual_uijs,
                self.target_weights)
        # Which datasets should be used for residual optimisation
        if self.residual_mask is not None:
            self.f_g_calculator.set_residual_mask(self.residual_mask)

        # Apply residual amplitudes if given
        if self.residual_amplitudes is not None:
            current_amps = self.f_g_calculator.x
            isel = flex.size_t_range(current_amps.size()-self.residual_uijs.size(), current_amps.size())
            assert len(isel) == len(self.residual_uijs)
            assert len(isel) == len(self.residual_amplitudes)
            assert max(isel) == current_amps.size() - 1
            current_vals = current_amps.select(isel)
            assert current_vals.all_eq(1.0)
            self.f_g_calculator.x.set_selected(isel, self.residual_amplitudes)

        self.initial = copy.deepcopy(self.f_g_calculator.x)

        # Optimise
        t_params = lbfgs.termination_parameters(
                traditional_convergence_test_eps=self.convergence_tolerance,
                )
        minimizer = lbfgs.run(
                target_evaluator=self.f_g_calculator,
                termination_params=t_params,
                )

        # Extract
        self.n_iter = minimizer.iter()
        self.result = numpy.array(self.f_g_calculator.x)
        self.result[self.result < 0.0] = 0.0 # Report this!?!

        # Tidy up
        del self.f_g_calculator # Delete after usage as not pickle-enabled

        return self


class InterLevelAmplitudeOptimiser(BaseOptimiser):


    def __init__(self,
                 target_uij,
                 levels,
                 residual,
                 group_tree,
                 weights=None,
                 params=None,
                 verbose=False,
                 log=None):
        super(InterLevelAmplitudeOptimiser, self).__init__(
                target_uij=target_uij,
                convergence_tolerance=params.optimisation.gradient_convergence,
                label=None,
                weights=weights,
                params=params,
                verbose=verbose,
                log=log)

        self.n_lvl = len(levels)
        self.n_tls = self.params.n_tls_modes_per_tls_group
        self.n_dst = target_uij.shape[0]
        self.n_atm = target_uij.shape[1]

        self.group_tree = group_tree

        self._ready = False

        # Dictionaries for selecting levels/index of levels in list from level.index (which counts from 1)
        # Names of the levels (confusingly also called index)
        self.level_labels = [l.index for l in levels]
        # Map between l.index and i_level
        self.level_index_hash  = {l.index:i for i,l in enumerate(levels)}
        self.level_index_hash['X'] = len(levels)
        # Map between l.index and number of groups in level
        self.level_n_groups_hash = {l.index:l.n_groups() for l in levels}
        # Mapping between (level, group) and atom selection
        self.level_group_selection_hash = {}
        for l in levels:
            for g_n, g_s, g_f in l:
               self.level_group_selection_hash.setdefault(l.index,{})[g_n] = g_s

        # Which datasets to be used for optimising residuals
        self.set_residual_mask(numpy.ones(self.n_dst, dtype=bool))

        # Extract the current uijs for levels and residual
        self.refresh(levels=levels, residual=residual)

    def set_residual_mask(self, mask):
        assert mask.size == self.n_dst
        self.residual_mask = mask

    def refresh(self, levels, residual):
        """Update the input uijs and the group multipliers"""
        # Extract level uijs for all levels and atoms
        self.level_uijs = numpy.array([l.extract_core() for l in levels])
        assert self.level_uijs.shape == (self.n_lvl, self.n_tls, self.n_dst, self.n_atm, 6)
        self.residual_uijs = residual.extract(expanded=False)
        # Amplitudes for each of the groups (what will be optimised in this object)
        self.group_multipliers = self._extract_multipliers(levels, residual)
        self._ready = True

    def _extract_multipliers(self, levels, residual):
        mults = [numpy.zeros(shape=(l.n_groups(), self.n_tls, self.n_dst)) for l in levels] + \
                [numpy.ones(residual.n_atoms())] # These are relative to the input
        for i_l, l in enumerate(levels):
            for i_g, (n_g, g_s, g_f) in enumerate(l):
                for i_tls in range(self.n_tls):
                    m = g_f.parameters().get(i_tls)
                    mults[i_l][i_g, i_tls, :] = m.amplitudes.values
        return mults

    def apply_multipliers(self, levels, residual):
        """Apply the optimised multipliers to the input levels"""
        # Iterate through levels, fitters and apply multiplier to each set of tls amplitudes
        for i_lvl, level in enumerate(levels):
            for i_g, (g_n, g_s, g_f) in enumerate(level):
                mults = self.group_multipliers[i_lvl][i_g]
                for i_tls in range(self.n_tls):
                    mode = g_f.parameters().get(index=i_tls)
                    mode.amplitudes.set(mults[i_tls])
        # And for the residual level
        for a_n, a_s, a_f in residual:
            mults = self.group_multipliers[-1][a_s]
            vals = a_f.parameters()
            vals[:] = mults * vals

        self._ready = False

    def _generate_level_group_sets(self, tree, max_recursions=None):

        def get_linked_level_group_pairs_recursive(tree, level, group, max_recursions=None):
            t_out = []
            linked_groups = tree.get(level, {}).get(group, None)
            if linked_groups is None:
                return []
            # Get the linked groups on this level
            for l in sorted(linked_groups.keys()):
                if l == 'X': continue
                for g in linked_groups.get(l):
                    t_out.append((l,g))
            # Get the linked groups from the next level
            r_out = []
            if (max_recursions is None) or (max_recursions > 0):
                if max_recursions is not None:
                    max_recursions -= 1
                for l,g in t_out:
                    r_out.extend(get_linked_level_group_pairs_recursive(tree, l, g, max_recursions))
            return t_out + r_out

        job_list = []
        n_links = len(tree.keys()) # Number of links between levels
        # Need to generate sets for the first n_levels - recursions
        for i_lvl, l1 in enumerate(sorted(tree.keys())):
            # Check if already optimised all
            if (i_lvl>0) and ((max_recursions is None) or (i_lvl+max_recursions > n_links)):
                break
            l_jobs = []
            l_tree = tree[l1]
            for g1 in sorted(l_tree.keys()):
                lg_pairs = [(l1,g1)]
                if (max_recursions is None) or (max_recursions > 0):
                    if max_recursions is not None:
                        max_recursions -= 1
                    linked_pairs = get_linked_level_group_pairs_recursive(tree, l1, g1, max_recursions)
                    lg_pairs.extend(linked_pairs)
                l_jobs.append(lg_pairs)
            job_list.append(l_jobs)

        return job_list

    def _generate_optimisation_jobs_from_sets(self, level_group_set, optimise_residual=True):
        """Generate control list for optimisation -- if datasets are independent, can be optimised separately, etc.."""

        dataset_indices = range(self.n_dst)

        output = []
        for s in level_group_set:
            # Generate output groups for each dataset, or combined, depending on whether contains residual
            if optimise_residual:
                output.append((dataset_indices, s))                 # One optimisation job
            else:
                output.extend([([i], s) for i in dataset_indices])  # One optimsation job for each dataset
        return output

    #=======================+>

    def get_optimisation_variables(self):
        return self._op

    def _extract_target_for_optimise(self):
        """Subtract fitted uijs for groups not in exclude_level_group_pairs"""

        op = self._op

        l_mults = numpy.zeros_like(self.level_uijs)

        for i_l in xrange(self.n_lvl):
            l_n = self.level_labels[i_l]
            for i_g in xrange(self.level_n_groups_hash[l_n]):
                g_n = i_g + 1
                # Skip the groups being optimised
                if op.optimisation_groups_mask[i_l][i_g]:
                    continue
                g_s = self.level_group_selection_hash[l_n][g_n]
                # Skip groups outside the atom selection
                if not (g_s * op.b_atm).any():
                    continue
                # Extract multipliers
                for i_t in xrange(self.n_tls):
                    for i_d in op.i_dst: # Only need to iterate over selected datasets
                        l_mults[i_l, i_t, i_d, g_s] = self.group_multipliers[i_l][i_g, i_t, i_d]

        # Multiply the fitted uijs and subtract them from target
        o_fitted = (l_mults * self.level_uijs).sum(axis=(0,1))
        if not op.optimise_residual:
            o_fitted = o_fitted + self.residual_uijs

        # Extract the ones required
        t_uij = (self.target_uij - o_fitted)[op.i_dst][:, op.b_atm]
        assert t_uij.shape == (op.n_dst, op.n_atm, 6)

        # Reshape for output
        t_uij_flex = flex.sym_mat3_double(t_uij.reshape((op.n_dst*op.n_atm,6)))
        t_uij_flex.reshape(flex.grid((op.n_dst,op.n_atm)))

        # Extract optimisation weights
        wgts = self.weights[op.i_dst][:,op.b_atm]
        # Renormalise weights
        wgts = wgts / wgts.mean()
        # Do with wgts as n_dst*n_atm (but as 1d array for speed!)
        wgts = flex.double(wgts.reshape((op.n_dst*op.n_atm,)).tolist())
        wgts.reshape(flex.grid((op.n_dst, op.n_atm)))

        return t_uij_flex, wgts

    def _extract_fitted_for_optimise(self):
        op = self._op

        muls_base = flex.double(op.n_base, 1.0)
        uijs_base = []
        sels_base = []

        # Calculate new atomic indexes of the truncated set of atoms from the full set of indices
        i_atm_rel = -1*numpy.ones(self.n_atm, dtype=int)
        i_atm_rel[op.b_atm] = numpy.arange(op.n_atm)

        i = 0
        for (l, g) in op.optimisation_groups_list:
            # Map labels to indices
            i_l = self.level_index_hash[l]
            i_g = g - 1
            # Extract atom selection for group
            g_s = self.level_group_selection_hash[l][g]

            # Combined selection (relative to the full set of atoms)
            # Not sure of any case where the whole group would not be selected, but ...
            g_s_mult = g_s*op.b_atm
            # Map to indices on the masked array (relative to the truncated set)
            new_i_sel = flex.size_t(i_atm_rel[g_s_mult].tolist())
            n_sel = len(new_i_sel)
            # Extract uijs for tls modes and datasets
            for i_t in op.i_tls:
                for i_d in op.i_dst:
                    u = self.level_uijs[i_l, i_t, i_d, g_s_mult]
                    uijs_base.append(flex.sym_mat3_double(u.reshape((n_sel,6))))
                    sels_base.append(new_i_sel)
                    muls_base[i] = self.group_multipliers[i_l][i_g, i_t, i_d]
                    i += 1
        assert i == op.n_base

        # Extract residual to be optimised
        if op.optimise_residual:
            uijs_res = flex.sym_mat3_double(self.residual_uijs[op.b_atm])
            muls_res = flex.double(self.group_multipliers[-1][op.b_atm].tolist())
        else:
            uijs_res = flex.sym_mat3_double(op.n_atm, (0.,0.,0.,0.,0.,0.))
            muls_res = None

        # Mark which elements belong to which datasets
        dataset_hash = flex.size_t(range(op.n_dst)*op.n_grp*op.n_tls)

        return (muls_base, uijs_base, sels_base, uijs_res, muls_res, dataset_hash)

    def set_variables(self, level_group_pairs, dataset_indices, optimise_residual):
        """Prepare for optimisation"""

        # Generate combined atom mask for all groups -- which atoms are included in optimisation
        atom_selection_bool = numpy.zeros(self.n_atm, dtype=bool)
        for (l,g) in level_group_pairs:
            sel = self.level_group_selection_hash[l][g]
            atom_selection_bool[sel] = True
        atom_selection_isel = flex.size_t(numpy.where(atom_selection_bool)[0].tolist())

        # Numbers for checks and iterators
        n_grp = len(level_group_pairs) # Number of non-residual groups
        n_tls = self.n_tls
        n_dst = len(dataset_indices)
        n_atm = len(atom_selection_isel)
        # Number of parameters in optimisation
        n_base = (n_grp * n_tls * n_dst) # Total number of amplitudes
        n_resl = n_atm * int(optimise_residual)

        # Create a boolean mask for the levels/groups/datasets being optimised
        level_group_mask = [numpy.zeros(level_row.shape[0], dtype=bool) for level_row in self.group_multipliers]
        # Create mapping between optimisation values and the level/group/affected uijs
        for i_grp, (l,g) in enumerate(level_group_pairs):
            i_l = self.level_index_hash[l]
            i_g = (g - 1) # XXX Possibly dangerous shortcut XXX
            # Record this group as being optimised
            level_group_mask[i_l][i_g] = True

        self._op = Info({
            'i_dst' : dataset_indices,      # Indices of datasets in optimisation (relative to full list)
            'i_tls' : range(n_tls),         # Indices of tls modes in optimisation (relative to full list)
            'i_atm' : atom_selection_isel,  # Atom mask relative to the full atom array
            'b_atm' : atom_selection_bool,  # Atom mask relative to the full atom array
            'n_grp' : n_grp,                # Number of groups being optimised
            'n_tls' : n_tls,                # Number of tls modes in optimisation (currently all)
            'n_dst' : n_dst,                # Number of datasets in optmisation
            'n_atm' : n_atm,                # Total number of atoms in uij array
            'n_base' : n_base,              # Number of amplitudes from base elements
            'n_resl' : n_resl,              # Number of amplitudes from residual
            'n_amps' : n_base + n_resl,     # Total number of amplitudes
            # Arrays for the optimisation variables
            'optimisation_groups_list' : level_group_pairs,
            'optimisation_groups_mask' : level_group_mask,
            # Flag
            'optimise_residual': optimise_residual,
            })

    def _extract_values(self, variables=None):
        """Extract a set of values for a set of variables"""
        # Use default selection if not provided
        if variables is None:
            variables = self._op
        # Iterate through objects and inject requested parameter values
        values = []
        for (l, g) in variables.optimisation_groups_list:
            i_l = self.level_index_hash[l]
            i_g = g - 1
            if (l != 'X'):
                values.append(self.group_multipliers[i_l][i_g, variables.i_tls, variables.i_dst].flatten())
            else:
                values.append([self.group_multipliers[-1][i_g]])
        return numpy.concatenate(values)

    def _inject_values(self, values, variables=None):
        """Change a set of values for a set of variables"""
        # Use default selection if not provided
        if variables is None:
            variables = self._op
        # Iterate through objects and extract requested parameter values
        values = collections.deque(values)
        for (l, g) in variables.optimisation_groups_list:
            i_l = self.level_index_hash[l]
            i_g = g - 1
            for i_tls in variables.i_tls:
                for i_dst in variables.i_dst:
                    self.group_multipliers[i_l][i_g, i_tls, i_dst] = values.popleft()

        # Remaining values should be the amplitudes of the residuals (regardless of whether optimised)
        assert len(values) == variables.n_atm

        if variables.optimise_residual:
            self.group_multipliers[-1][variables.b_atm] = numpy.array(values)
        else:
            assert (numpy.array(values) == 1.0).all()

        return

    #===========================================+>

    def optimise(self, n_cpus=1, max_recursions=None, optimise_residual=True):
        """Optimise the amplitudes of a series of TLS/Uij Levels"""
        assert self._ready is True, 'Not ready for optimise -- use .refresh() first'

        # Optimisation sets of which amplitudes should be co-optimised
        level_group_sets = self._generate_level_group_sets(self.group_tree, max_recursions=max_recursions)

        self.log.heading('Running Amplitude Optimisations')

        for lg_set in level_group_sets:
            jobs_input = self._generate_optimisation_jobs_from_sets(lg_set, optimise_residual=optimise_residual)
            # Perform jobs sequentially
            chunk_size = min(len(jobs_input), 100*n_cpus)
            for i in xrange(0, len(jobs_input), chunk_size):
                self.log('Running jobs {} -> {} (of {})'.format(1+i, min(1+i+chunk_size, len(jobs_input)), len(jobs_input)))
                j_chunk = jobs_input[i:i+chunk_size]
                self.optimise_level_amplitudes(jobs_input=j_chunk, n_cpus=n_cpus, optimise_residual=optimise_residual)

    def optimise_level_amplitudes(self, jobs_input, n_cpus=1, optimise_residual=True):
        # Optimise all amplitudes dataset-by-dataset
        jobs = []

        # Generate jobs (and process if cpus==1)
        n_cpus = min(n_cpus, len(jobs_input))
        self.log.subheading('Preparing optimisations...')
        for i, (dst_indices, lvl_groups) in enumerate(jobs_input):
            self.log.bar()
            self.log('Optimisation {}'.format(i+1))
            self.log('Datasets: {}'.format(str(dst_indices)))
            i = 0
            srt_groups = sorted(lvl_groups)
            while i < len(srt_groups):
                lvl, grp = srt_groups[i]
                all_grps = []
                cur_lvl = lvl
                cur_grp = grp
                while (cur_lvl == lvl):
                    all_grps.append(cur_grp)
                    i += 1
                    if i == len(srt_groups):
                        break
                    cur_lvl, cur_grp = srt_groups[i]
                self.log('Groups from Level {}:'.format(lvl))
                for j in range(0,len(all_grps), 10):
                    self.log('\t' + ', '.join(['{:>5}'.format(v) for v in all_grps[j:j+10]]))
            # Select this dataset
            self.set_variables(
                    dataset_indices=dst_indices,
                    level_group_pairs=lvl_groups,
                    optimise_residual=optimise_residual,
                    )
            # Append to job list (create copy) or run optimisation
            if n_cpus > 1:
                jobs.append((self.get_optimiser(), {'stdout':False}))
            else:
                self._optimise()
                #print self._extract_values()

        # Run parallel jobs and inject results
        if n_cpus > 1:
            # Workers for multiprocessing
            workers = DaemonicPool(n_cpus)
            # Report and map to workers
            self.log.subheading('Running {} jobs using {} cpus'.format(len(jobs), min(n_cpus,len(jobs))))
            finished_jobs = workers.map(func=_wrapper_optimise, iterable=jobs)
            # Inject results from jobs
            #self.log('Optimisation Results:')
            for job in finished_jobs:
                if isinstance(job, str):
                    self.log(job)
                    raise Failure('error returned')
                self._inject_values(values=job.result, variables=job.variables)
            # Let's be good
            workers.close()

    def _optimise(self):
        """Run the optimisation for the current set of variables"""
        # Create optimisation object
        optimiser = self.get_optimiser()._optimise()
        # Extract and update current values
        self._inject_values(values=optimiser.result, variables=self.get_optimisation_variables())

    def get_optimiser(self):
        # Extract fitted uijs for each mode
        muls_base, uijs_base, sels_base, uijs_res, muls_res, dataset_hash = self._extract_fitted_for_optimise()
        # Extract target uijs for the relevant atoms
        t_uij, wgts = self._extract_target_for_optimise()
        # Optimise these parameters
        optimiser = OptimiseAmplitudes(
                target_uijs         = t_uij,
                target_weights      = wgts,
                base_amplitudes     = muls_base,
                base_uijs           = uijs_base,
                base_atom_indices   = sels_base,
                dataset_hash        = dataset_hash,
                residual_uijs       = uijs_res,
                residual_amplitudes = muls_res,
                residual_mask       = self.residual_mask[self._op.i_dst],
                convergence_tolerance = self.convergence_tolerance,
                )
        return optimiser

    def optimisation_summary(self, show=True):
        """Print the optimisation values and weights"""
        return ""
    def summary(self, show=True):
        """Print the number of parameters/input data"""
        return ""


#!/usr/bin/env ccp4-python

import os, sys, re, glob, shutil, copy, tempfile, gc
import math, re, time

import pandas, numpy

import libtbx.phil, libtbx.easy_mp
import iotbx.pdb
import mmtbx.tls.tools

from libtbx.utils import Sorry, Failure
from scitbx.array_family import flex
from scitbx import simplex

from bamboo.common.logs import Log

from giant.dataset import CrystallographicModel
from giant.structure.select import protein
from giant.structure.tls import uij_from_tls_vector_and_origin, extract_tls_from_pdb

import matplotlib
matplotlib.interactive(False)
from matplotlib import pyplot
#pyplot.switch_backend('agg')
pyplot.interactive(0)

numpy.set_printoptions(threshold=numpy.nan)

from IPython import embed
log = Log(verbose=True)

############################################################################

PROGRAM = 'giant.datasets.b_factor_fitting'

DESCRIPTION = """
    Analyse the variation/conservartion of B-factors (under different models) of a set of related structures.
"""

############################################################################

blank_arg_prepend = {'.pdb':'pdb=', '.cif':'cif='}

master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .help = "input pdb files - with anisotropic b-factors"
        .multiple = True
        .type = str
    labelling = filename *foldername
        .type = choice
        .multiple = False
    tls_group = None
        .type = str
        .multiple = True
}
output {
    out_dir = multi-dataset-b-factor-fitting
        .help = "output directory"
        .type = str
    log_file = multi-dataset-b-factor-fitting.log
        .type = str
}
settings {
    cpus = 48
        .type = int
        .multiple = False
}
""")

############################################################################

def proc_wrapper(arg):
    arg._optimise()
    return arg

class MultiDatasetTLSFitter(object):

    def __init__(self, observed_xyz, observed_uij, tls_params=None, num_tls_models=None, cpus=48):

        self._test = False
        self._iter = 0
        self._cpus = cpus

        self.optimisation_rmsd = numpy.inf

        # ---------------------------->
        # Input data
        # ---------------------------->
        self.observed_uij = numpy.array(observed_uij)
        self.observed_xyz = numpy.array(observed_xyz)
        # TODO Make this variable over the datasets ? TODO
        self.observed_com = numpy.mean(self.observed_xyz, axis=(0,1))
        # ---------------------------->
        # Process input TLS options
        # ---------------------------->
        if tls_params is not None:
            inp_tls = tls_params
            num_tls = len(tls_params)
            vec_tls = [p.t+p.l+p.s for p in tls_params]
        elif num_tls_models is not None:
            inp_tls = None
            num_tls = num_tls_models
            vec_tls = [numpy.mean(self.observed_uij, axis=(0,1)).tolist()+[0.0]*15]+[[0.0]*21]*(num_tls-1)
        else:
            raise Sorry('No TLS models provided')

        assert len(vec_tls) == num_tls
        assert set(map(len, vec_tls)) == {21}

        # ---------------------------->
        # Extract variables and validate
        # ---------------------------->
        # Extract size of objects
        self._n_dst = self.observed_xyz.shape[0]
        self._n_atm = self.observed_xyz.shape[1]
        self._n_tls = num_tls
        # Number of observed data
        self._n_obs_total = numpy.product(self.observed_uij.shape)
        # Calculate number of fitting parameters: 21 parameters per tls model + 3 amplitudes per TLS model per dataset + 6 residual uij per atom
        self._n_prm_tls_mdl = 21 * self._n_tls
        self._n_prm_tls_amp = 3  * self._n_dst * self._n_tls
        self._n_prm_uij_res = 6  * self._n_atm
        self._n_prm_total = self._n_prm_tls_mdl + self._n_prm_tls_amp + self._n_prm_uij_res

        # Check we're not in dangerous territory...
        assert self._n_obs_total > self._n_prm_total

        # ---------------------------->
        # Initial parameter values
        # ---------------------------->
        self.initial_tls_parameters = numpy.concatenate(vec_tls)
        self.initial_tls_amplitudes = numpy.array([1.0]*self._n_prm_tls_amp)
        self.initial_uij_residuals  = numpy.array([0.0]*self._n_prm_uij_res)
        # ---------------------------->
        # Output variables (initialise to None)
        # ---------------------------->
        self.fitted_tls_parameters = None
        self.fitted_tls_amplitudes = None
        self.fitted_uij_residuals  = None

        # ---------------------------->
        # Parameter group selections
        # ---------------------------->
        # All TLS parameters
        self._sel_tls_mdl = numpy.array([1]*self._n_prm_tls_mdl + [0]*self._n_prm_tls_amp + [0]*self._n_prm_uij_res, dtype=bool)
        # All TLS amplitudes
        self._sel_tls_amp = numpy.array([0]*self._n_prm_tls_mdl + [1]*self._n_prm_tls_amp + [0]*self._n_prm_uij_res, dtype=bool)
        # All residual uijs
        self._sel_uij_res = numpy.array([0]*self._n_prm_tls_mdl + [0]*self._n_prm_tls_amp + [1]*self._n_prm_uij_res, dtype=bool)
        # ---------------------------->
        # Make selections for each component of the TLS model (parameters and amplitudes)
        # ---------------------------->
        self._sel_t = self._blank_parameter_selection()
        self._sel_l = self._blank_parameter_selection()
        self._sel_s = self._blank_parameter_selection()
        # Any "T" component of a TLS parameter set (regularly spaced)
        self._sel_t[self._sel_tls_mdl] = ([1]*6 + [0]*6 + [0]*9)*self._n_tls
        self._sel_l[self._sel_tls_mdl] = ([0]*6 + [1]*6 + [0]*9)*self._n_tls
        self._sel_s[self._sel_tls_mdl] = ([0]*6 + [0]*6 + [1]*9)*self._n_tls
        # Any "T" amplitude of a dataset parameter set (regularly spaced)
        self._sel_t[self._sel_tls_amp] = [1,0,0]*self._n_tls*self._n_dst
        self._sel_l[self._sel_tls_amp] = [0,1,0]*self._n_tls*self._n_dst
        self._sel_s[self._sel_tls_amp] = [0,0,1]*self._n_tls*self._n_dst
        # ---------------------------->
        # (list) Make amplitude selections for each dataset
        # ---------------------------->
        self._sel_dst = [self._blank_parameter_selection() for i in range(self._n_dst)]
        for i in range(self._n_dst):
            self._sel_dst[i][self._sel_tls_amp] = [0]*3*self._n_tls*i + [1]*3*self._n_tls + [0]*3*self._n_tls*(self._n_dst-i-1)
        # ---------------------------->
        # (list) Make uij selections for each atom
        # ---------------------------->
        self._sel_atm = [self._blank_parameter_selection() for i in range(self._n_atm)]
        for i in range(self._n_atm):
            self._sel_atm[i][self._sel_uij_res] = [0]*6*i + [1]*6 + [0]*6*(self._n_atm-i-1)
        # ---------------------------->
        # (list) Make selections for each TLS model (parameters and amplitudes)
        # ---------------------------->
        self._sel_tls = [self._blank_parameter_selection() for i in range(self._n_tls)]
        for i in range(self._n_tls):
            self._sel_tls[i][self._sel_tls_mdl] = ([0]*21*i + [1]*21 + [0]*21*(self._n_tls-i-1) )
            self._sel_tls[i][self._sel_tls_amp] = ([0]*3*i  + [1]*3  + [0]*3*(self._n_tls-i-1)  )*self._n_dst
        # ---------------------------->
        # Check each selection has the right number of parameters
        # ---------------------------->
        assert set(map(numpy.sum, self._sel_dst)) == {3*self._n_tls}
        assert set(map(numpy.sum, self._sel_atm)) == {6}
        assert set(map(numpy.sum, self._sel_tls)) == {21+3*self._n_dst}

        # ---------------------------->
        # Atom and dataset Masks
        # ---------------------------->
        self._mask_dsets = numpy.ones(self._n_dst, dtype=bool)
        self._mask_atoms = numpy.ones(self._n_atm, dtype=bool)

        # ---------------------------->
        # Minimisation variables
        # ---------------------------->
        # Amount of expected variation for each data type (axes of simplex)
        self._del_tls_mdl = 0.1
        self._del_tls_amp = 0.1
        self._del_uij_res = 0.1
        self._del_simplex = numpy.array([self._del_tls_mdl]*self._n_prm_tls_mdl +
                                        [self._del_tls_amp]*self._n_prm_tls_amp +
                                        [self._del_uij_res]*self._n_prm_uij_res )
        # ---------------------------->
        # Loop variables during optimisation
        # ---------------------------->
        # Running parameter vector (Changes during target function during optimisation cycles)
        self._var_current = numpy.concatenate([self.initial_tls_parameters,
                                               self.initial_tls_amplitudes,
                                               self.initial_uij_residuals])
        # Updated at the end of each optimisation cycle
        self._var_optimised = None
        # Select which parameters to optimise
        self._var_current_sel = None
        # ---------------------------->
        # Report
        # ---------------------------->
        self.input_summary()
        self.parameter_summary()

    ################################################################################################
    ###
    ###                             Miscellaneous
    ###
    ################################################################################################

    def copy(self):
        return copy.deepcopy(self)

    def _prep_for_mp(self):
        """Clear parts of the object that are not needed for optimisation after selection - save memory when pickling"""
        self.initial_tls_parameters = None
        self.initial_tls_amplitudes = None
        self.initial_uij_residuals  = None
        return self

    ################################################################################################
    ###
    ###                             Summaries
    ###
    ################################################################################################

    def input_summary(self):
        """Print the number of parameters/input data"""
        log.subheading('Input summary')
        log('input uij parameters: {}'.format(self.observed_uij.shape))
        log('input xyz parameters: {}'.format(self.observed_xyz.shape))
        log('Centre of mass: {}'.format(self.observed_com))
        log.bar()
        log('Number of parameters for TLS fitting: {}'.format(self._n_prm_total))
        log('Number of observations: {}'.format(self._n_obs_total))
        log('Data/parameter ratio is {:.3f}'.format(self._n_obs_total*1.0/self._n_prm_total))
        log.bar()
        log('Number of datasets:   {}'.format(self._n_dst))
        log('Number of atoms:      {}'.format(self._n_atm))
        log('Number of TLS models: {}'.format(self._n_tls))
        log.bar()

    def optimisation_summary(self, full=True):
        """Print the fitted parameters"""
        tls_prms, tls_amps, res_uij = self._extract(vector=self._var_optimised)

        log.subheading('Optimisation Summary')

        tls_refined = numpy.sum(self._sel_tls_mdl*self._var_current_sel)
        amp_refined = numpy.sum(self._sel_tls_amp*self._var_current_sel)
        uij_refined = numpy.sum(self._sel_uij_res*self._var_current_sel)

        log('Optimisation used data from {} datasets'.format(len(self._cur_datasets)))
        log('Optimisation used data for {} atoms'.format(len(self._cur_atoms)))
        log.bar()
        log('Optimised {} parameters'.format(numpy.sum(self._var_current_sel)))
        log('... {} TLS parameters,'.format(tls_refined))
        log('... {} TLS amplitudes,'.format(amp_refined))
        log('... {} uij parameters.'.format(uij_refined))
        log.bar()
        log('Optimisation RMSD: {}'.format(self.optimisation_rmsd))
        log.bar()
        log('')
        if tls_refined or full:
            self.parameter_summary(tls_params=tls_prms)
        if amp_refined or full:
            self.parameter_summary(tls_amplitudes=tls_amps)
        if uij_refined or full:
            self.parameter_summary(uij_residual=res_uij)

    def parameter_summary(self, tls_params=None, tls_amplitudes=None, uij_residual=None):
        if tls_params is tls_amplitudes is uij_residual is None:
            tls_params, tls_amplitudes, uij_residual = self._extract(vector=self._var_current)
        log.bar()
        if tls_params is not None:
            log('TLS Parameters:')
            for i, vals in enumerate(tls_params):
                log.bar()
                log('Model {:2}'.format(i+1))
                log('T: '+', '.join(['{:8.3f}'.format(v) for v in vals[00:06]]))
                log('L: '+', '.join(['{:8.3f}'.format(v) for v in vals[06:12]]))
                log('S: '+', '.join(['{:8.3f}'.format(v) for v in vals[12:21]]))
            log.bar()
        if tls_amplitudes is not None:
            log('TLS Amplitudes:')
            for i, vals in enumerate(tls_amplitudes):
                log.bar()
                log('Dataset {:4}'.format(i+1))
                for i, amps in enumerate(vals):
                    log('Model {:2}:'.format(i+1)+' {:8.3f} (T) {:8.3f} (L) {:8.3f} (S)'.format(*amps))
            log.bar()
        if uij_residual is not None:
            log('Residual Uijs:')
            for i, vals in enumerate(uij_residual):
                log('Atom {:5}: '.format(i+1) + ', '.join(['{:8.3f}'.format(v) for v in vals]))
            log.bar()

    ################################################################################################
    ###
    ###                             Update internal variables or masks
    ###
    ################################################################################################

    def _blank_parameter_selection(self):
        return numpy.zeros(self._n_prm_total, dtype=bool)
    def _blank_dataset_selection(self):
        return numpy.zeros(self._n_dst, dtype=bool)
    def _blank_atom_selection(self):
        return numpy.zeros(self._n_atm, dtype=bool)

    def reset_tls_amplitudes(self):
        self._var_current[self._sel_tls_amp] = 1.0
    def reset_uij_residual(self):
        self._var_current[self._sel_uij_res] = 0.0
    def reset_current_selection(self):
        self._var_current_sel = self._blank_parameter_selection()

    def _update_atom_masks(self):
        _, _, uij_res = self._extract()
        uij_max = numpy.max(numpy.abs(uij_res[:,:3]),axis=1)
        thresh = numpy.percentile(uij_max, 90)
        self._mask_atoms *= (uij_max < thresh)

    def _update_dataset_masks(self):
        pass

    def _extract(self, vector=None):
        """Convert 1d vector into objects"""
        if vector is None: vector=self._var_current
        tls_mdls = vector[self._sel_tls_mdl].reshape((self._n_tls, 21)             )
        tls_amps = vector[self._sel_tls_amp].reshape((self._n_dst, self._n_tls, 3) )
        uij_resl = vector[self._sel_uij_res].reshape((self._n_atm, 6)              )
        return (tls_mdls,tls_amps,uij_resl)

    ################################################################################################
    ###
    ###                             Update internal variables after optimising/adopting
    ###
    ################################################################################################

    def _update_output(self):
        # Copy running parameters to output holder
        self._var_optimised = self._var_current.copy()

    ################################################################################################
    ###
    ###                             Change the parameters of the fitter
    ###
    ################################################################################################

    def _adopt(self, sub_vector, selection=None):
        """Insert a set of parameters into the complete parameter set"""
        if selection is None: selection = self._var_current_sel
        assert len(sub_vector) == numpy.sum(selection)
        self._var_current[selection] = sub_vector

    def _adopt_from_others(self, others):
        # Updating from others, so set own selection to none
        self.reset_current_selection()
        for n in others:
            # Add this selection to our selection
            self._var_current_sel += n._var_current_sel
            # Copy the values across
            self._adopt(sub_vector=n._var_current[n._var_current_sel], selection=n._var_current_sel)
        self._update_output()

    ################################################################################################
    ###
    ###                             Main methods for running/optimisation
    ###
    ################################################################################################

    def run(self, macro_cycles=3):
        """Run macro-cycles of parameter optimisation"""
        for i_cyc in range(macro_cycles):

            #########################################
            log.heading('Macrocycle {}'.format(i_cyc+1), spacer=True)
            #########################################
            # Remove outlier atoms from minimisation
            #########################################
            if i_cyc > 0:
                log.subheading('Removing atoms with high residual uij from TLS optimisation')
                self._update_atom_masks()
                log('Using {} atoms for TLS optimisation'.format(numpy.sum(self._mask_atoms)))

                #self._update_dataset_masks()

            #########################################
            log('Resetting amplitudes and uij residuals')
            self.reset_tls_amplitudes()
            self.reset_uij_residual()

            #########################################
            for i_tls in range(self._n_tls):
                for c_name, c_sel in [('T',self._sel_t), ('L',self._sel_l), ('S',self._sel_s)]:
                    log.heading('Optimising {} parameters for TLS model {}'.format(c_name, i_tls+1))
                    self._select(parameter_selection = self._sel_tls_mdl*self._sel_tls[i_tls]*c_sel,
                                 datasets = None,
                                 atoms    = None)
                    self._optimise()
                    self.optimisation_summary(False)
                    #########################################
                    log.heading('Optimising {} amplitudes for TLS model {}'.format(c_name, i_tls+1))
                    proc_args = []
                    for i_dst in range(self._n_dst):
                        n = self.copy()
                        n._select(parameter_selection = self._sel_dst[i_dst]*self._sel_tls[i_tls]*c_sel,
                                  datasets = [i_dst],
                                  atoms    = None)
                        proc_args.append(n._prep_for_mp())
                    self._adopt_from_others(libtbx.easy_mp.pool_map(processes=self._cpus, func=proc_wrapper, args=proc_args))
                    self.optimisation_summary(False)
                #########################################
                log.heading('Optimising all parameters for TLS model {}'.format(i_tls+1))
                self._select(parameter_selection = self._sel_tls_mdl*self._sel_tls[i_tls],
                             datasets = None,
                             atoms    = None)
                self._optimise()
                self.optimisation_summary(False)
            #########################################
            log.heading('Optimising residual Uijs')
            proc_args = []
            for i in range(self._n_atm):
                n = self.copy()
                n._select(parameter_selection = self._sel_atm[i],
                          datasets = None,
                          atoms    = [i])
                proc_args.append(n._prep_for_mp())
            self._adopt_from_others(libtbx.easy_mp.pool_map(processes=self._cpus, func=proc_wrapper, args=proc_args))
            self.optimisation_summary(False)
            #########################################
            log.subheading('End of macrocycle {}'.format(i_cyc+1))
            self.optimisation_summary()

            #self.optimisation_summary()

        return self

    def _select(self, parameter_selection, datasets=None, atoms=None):
        """Select variables for optimisation"""

        # Save the parameter selection
        self._var_current_sel = parameter_selection
        # Extract selections for datasets and atoms
        self._cur_datasets, self._cur_atoms = datasets, atoms
        # Override the dataset and atom selection
        if datasets is None:
            self._cur_datasets = list(numpy.where(self._mask_dsets)[0])
        if atoms is None:
            self._cur_atoms = list(numpy.where(self._mask_atoms)[0])

    def _optimise(self):
        """Run the optimisation"""

        # Initialise the RMSD measure
        self.optimisation_rmsd = 1e6
        # Create simplex for these parameters
        cur_simplex = self._get_simplex(self._var_current_sel)
        # Optimise these parameters
        optimised = simplex.simplex_opt(dimension = len(cur_simplex[0]),
                                        matrix    = map(flex.double, cur_simplex),
                                        evaluator = self)
        # Extract and update current values
        self._adopt(optimised.get_solution())
        self._update_output()

    def target(self, sub_vector):
        """Target function for the simplex optimisation"""
        # Combine the optimising parameters in the complete parameter set
        self._adopt(sub_vector)
        # Get the fitted and the observed uijs
        uij_fit = self.extract_fitted_uij(datasets=self._cur_datasets, atoms=self._cur_atoms)
        uij_obs = self.extract_observed_uij(datasets=self._cur_datasets, atoms=self._cur_atoms)
        # Calculate RMSD
        rmsd = numpy.sqrt(numpy.mean(numpy.power(uij_obs-uij_fit, 2)))
        self.optimisation_rmsd = min(self.optimisation_rmsd, rmsd)

        if self._test: embed()

        return rmsd

    def _get_simplex(self, selection):
        starting_values = self._var_current[selection]
        starting_deltas = self._del_simplex[selection]
        starting_simplex = numpy.repeat([starting_values], len(starting_values)+1, axis=0)
        for i in range(len(starting_values)):
            starting_simplex[i+1][i] += starting_deltas[i]
        return starting_simplex

    ################################################################################################
    ###
    ###                             Public Helper functions
    ###
    ################################################################################################

    def tls_uij(self, xyz, tls_vectors, origin):
        """Convert a set of parameter vectors to a set of uijs"""
        return numpy.sum([uij_from_tls_vector_and_origin(xyz=xyz, tls_vector=v, origin=origin) for v in tls_vectors], axis=0)

    def expand_tls_amplitudes(self, tls_amps):
        """Convert 3-element vector into 21 element vector for TLS multiplication"""
        cur_n_dst = len(tls_amps)
        assert tls_amps.shape == (cur_n_dst,self._n_tls,3)
        t_amps = numpy.repeat(tls_amps[:,:,0], 6, axis=1).reshape((cur_n_dst, self._n_tls, 6))
        l_amps = numpy.repeat(tls_amps[:,:,1], 6, axis=1).reshape((cur_n_dst, self._n_tls, 6))
        s_amps = numpy.repeat(tls_amps[:,:,2], 9, axis=1).reshape((cur_n_dst, self._n_tls, 9))
        exp_tls_amps = numpy.concatenate([t_amps, l_amps, s_amps], axis=2)
        assert exp_tls_amps.shape == (cur_n_dst,self._n_tls,21)
        return exp_tls_amps

    ################################################################################################
    ###
    ###                             Public Output functions
    ###
    ################################################################################################

    def get_dataset_and_atom_selection_bool(self, datasets=[], atoms=[]):
        datasets_bool = self._blank_dataset_selection()
        if (datasets is None) or (len(datasets) == 0):  datasets_bool += True   # select all
        else:                                           datasets_bool.put(datasets, True)
        atoms_bool = self._blank_atom_selection()
        if (atoms is None) or (len(atoms) == 0):        atoms_bool += True      # select all
        else:                                           atoms_bool.put(atoms, True)
        return datasets_bool, atoms_bool

    def selection_filter(self, datasets=None, atoms=None):
        """If either datasets or atoms is None, returns slice over all datasets/atoms"""
        if datasets is None: datasets = range(self._n_dst)
        if atoms    is None: atoms    = range(self._n_atm)
        return datasets, atoms

    def extract_observed_xyz(self, datasets=None, atoms=None):
        datasets, atoms = self.selection_filter(datasets=datasets, atoms=atoms)
        return self.observed_xyz[datasets][:,atoms]

    def extract_observed_uij(self, datasets=None, atoms=None):
        datasets, atoms = self.selection_filter(datasets=datasets, atoms=atoms)
        return self.observed_uij[datasets][:,atoms]

    def extract_fitted_uij(self, datasets=None, atoms=None):
        """Extract total fitted uijs for a subset of datasets or atoms"""
        datasets, atoms = self.selection_filter(datasets=datasets, atoms=atoms)
        # Extract the optimised values
        tls_p, tls_a, uij_r = self._extract()
        # Extract relevant coordinates
        xyz = self.extract_observed_xyz(datasets=datasets, atoms=atoms)
        assert xyz.shape == (len(datasets), len(atoms), 3)
        # Extract only those that we're interested in
        tls_a = tls_a[datasets]
        assert tls_a.shape == (len(datasets), self._n_tls, 3)
        uij_r = uij_r[atoms]
        assert uij_r.shape == (len(atoms), 6)
        # Multiply tls amplitudes and models
        tls_f = self.expand_tls_amplitudes(tls_amps=tls_a) * tls_p
        assert tls_f.shape == (len(datasets), self._n_tls, 21)

        assert len(xyz) == len(tls_f)
        # Calculate the tls component of uij
        uij_fit = numpy.array([self.tls_uij(xyz=xyz[i], tls_vectors=tls_f[i], origin=self.observed_com) for i in range(len(datasets))])
        assert uij_fit.shape == (len(datasets), len(atoms), 6)
        # Add the residual uij
        uij_fit += uij_r

        return uij_fit

class MultiDatasetTLSFit(object):

    def __init__(self, u_residual, tls_selections, tls_params):

        self.u_residual = u_residual
        self.tls_sels = tls_selections
        self.tls_objs = tls_objects

    def extract(self, hierarchy):
        pass

############################################################################

def run(params):

    print 'Building model list'
    models = [CrystallographicModel.from_file(f).label(tag=os.path.basename(os.path.dirname(f))) for f in params.input.pdb]

    # Check that all models are the same
    # TODO

#    # Extract the TLS model of the first dataset
#    log('No TLS models provided. Trying to extract from first crystallographic model')
#    tls_params = extract_tls_from_pdb(models[0].filename).tls_params
#    # Print the found models
#    for i_tls, tls in enumerate(tls_params):
#        log.subheading('TLS MODEL: {}'.format(i_tls+1))
#        log('T: {}'.format(tuple(tls.t)))
#        log('L: {}'.format(tuple(tls.l)))
#        log('S: {}'.format(tuple(tls.s)))
#    assert len(tls_params) != 0, 'No TLS models found in structure.'

    h=models[0].hierarchy
    for m in models:
        assert h.is_similar_hierarchy(m.hierarchy)

    for tls_sel in params.input.tls_group:
        # Extract the atoms for this tls selection
        atm_sel = h.atom_selection_cache().selection(tls_sel)

        log.heading('Calculating TLS fit for selection: {}'.format(tls_sel))

        # Extract filtered hierarchies
        atoms = [m.hierarchy.atoms().select(atm_sel) for m in models]

        # Extract uij and xyz
        obs_uij = numpy.array([a.extract_uij() for a in atoms])
        obs_xyz = numpy.array([a.extract_xyz() for a in atoms])

        log.subheading('Fitting TLS models to data')
        fitter = MultiDatasetTLSFitter(observed_uij=obs_uij,
                                       observed_xyz=obs_xyz,
                                       #tls_params=tls_params,
                                       num_tls_models=1,
                                       cpus=params.settings.cpus)

        # Calculate scaling
        fitter.run(2)

        # Extract the fitted output for each dataset
        for i_mdl, mdl in enumerate(models):
            log('Extracting B-factors for model: {}'.format(mdl.filename))
            mdl_uij = fitter.extract_fitted_uij(datasets=[i_mdl])
            h=mdl.hierarchy.deep_copy()
            h.atoms().select(atm_sel).set_uij(flex.sym_mat3_double(mdl_uij[0]))
            h.write_pdb_file(mdl.filename.replace('.pdb', '.mda.pdb'))

    #embed()



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



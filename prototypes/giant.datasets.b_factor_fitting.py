#!/usr/bin/env ccp4-python

import os, sys, re, glob, shutil, copy, tempfile, gc
import math, re

import mdp
import pandas, numpy

import libtbx.phil
import libtbx.easy_mp
import iotbx.pdb
import mmtbx.tls.tools

from libtbx.utils import Sorry, Failure
from scitbx.array_family import flex
from scitbx import simplex

from bamboo.common.path import easy_directory
from bamboo.common.command import CommandManager
from bamboo.common.logs import Log
from bamboo.stats.cluster import generate_group_idxs

from giant.dataset import CrystallographicModel
from giant.structure.b_factors import occupancy_weighted_average_b_factor
from giant.structure.select import protein, backbone, sidechains

import matplotlib
matplotlib.interactive(False)
from matplotlib import pyplot
#pyplot.switch_backend('agg')
pyplot.interactive(0)

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d

numpy.set_printoptions(threshold=numpy.nan)

EIGHT_PI_SQ = 8*math.pi*math.pi

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

from giant.stats.optimisation import _LeastSquaresFitter

############################################################################

def fit_adps_to_multiple_datasets(models):

    # Check that all models are the same
    # TODO

    ref = models[0].hierarchy
    cac = ref.atom_selection_cache()

    hierarchies = [protein(m.hierarchy) for m in models]

    # Extract uij
    obs_uij = numpy.array([h.atoms().extract_uij() for h in hierarchies])
    obs_xyz = numpy.array([h.atoms().extract_xyz() for h in hierarchies])

    fitter = MultiDatasetTLSFitter(observed_uij=obs_uij, observed_xyz=obs_xyz)
    fitter.run()

    from IPython import embed; embed()

class MultiDatasetTLSFitter(object):

    def __init__(self, observed_uij, observed_xyz):

        self.observed_uij = numpy.array(observed_uij)
        self.observed_xyz = numpy.array(observed_xyz)
        # TODO Make this variable over the datasets TODO
        self.observed_com = numpy.mean(self.observed_xyz[0])

        print 'UIJ', self.observed_uij[0].shape
        print 'XYZ', self.observed_xyz[0].shape
        print 'COM', self.observed_com

        self.n_obs = len(self.observed_xyz)
        self.n_ats = self.observed_xyz[0].shape[0]
        self.n_tls = 1

        self.observed_selection = flex.bool(self.n_ats, True)

        # Total params = 6 residual uij per atom + 21 parameters per tls + 1 per dataset per TLS
        self.total_params = 6*self.n_ats + 21*self.n_tls + self.n_obs*self.n_tls

        self.fitted_uij_tot = None
        self.fitted_uij_tls = None
        self.fitted_uij_res = None

        # Initialise simplex
        self.starting_values = self.initialise_parameters()
        # Calculate scaling
        self.run(initial_simplex=self.starting_simplex)

    def run(self, initial_simplex):
        """Calculate scaling"""
        # Optimise the simplex
        self.optimised = simplex.simplex_opt(dimension = len(initial_simplex[0]),
                                             matrix    = initial_simplex,
                                             evaluator = self)
        # Extract solution
        self.optimised_values = self.optimised.get_solution()

        # Scaled input values
        self.out_values = self.transform(values=self.scl_values)

        return self.optimised_values

    def extract_parameters(self, vector):
        """Convert 1d vector into objects"""

        n, m, p = (self.n_obs, self.n_ats, self.n_tls)

        n_offset     = 0
        uij_residual = [vector[n_offset + 6*i  : n_offset + 6*(i+1)  ] for i in range(m)]
        n_offset     = n_offset + 6*m
        tls_ampltdes = [vector[n_offset + p*i  : n_offset + p*(i+1)  ] for i in range(n)]
        n_offset     = n_offset + n*p
        tls_params   = [vector[n_offset + 21*i : n_offset + 21*(i+1) ] for i in range(p)]
        n_offset     = n_offset + 21*p
        assert n_offset == self.total_params

        return map(numpy.array, (uij_residual, tls_params, tls_ampltdes))

    def parameters_to_uijs(self, observed_xyz, tls_params, tls_amplitudes, uij_residual):
        """Convert a vector of parameters to a set of uijs"""

        assert len(observed_xyz) == len(uij_residual)

        tls_sels = [flex.bool(len(observed_xyz), True)]*len(tls_params)

        tls_p = numpy.array(tls_params)
        tls_a = numpy.array(tls_amplitudes)

            this_tls_amps = tls_ampltdes[i_obs]
            this_tls_prms = tls_params * tls_ampltdes[i_obs]
            this_tls_objs = [mmtbx.tls.tools.tlso(t=p[0:6], l=p[6:12], s=p[12:21], origin=self.observed_com) for p in this_tls_prms]
            all_uij_tls.append(mmtbx.tls.tools.u_cart_from_tls(sites_cart = flex.vec3_double(xyz),
                                                               selections = [self.observed_selections]*self.n_tls,
                                                               tlsos      = this_tls_objs))
        all_uij_tls = numpy.array(all_uij_tls)



    def target(self, vector):
        """Target function for the simplex optimisation"""
        uij_obs = self.observed_uij
        uij_residual, tls_params, tls_ampltdes = self.extract_parameters(vector=vector)

        for i_obs in range(self.n_obs):



        scaled = self.transform(values=self.scl_values, params=vector)
        diff = (scaled-self.ref_values)
        diff_sq = diff*diff
        result = flex.sum(self.weight_array*diff_sq)
        return result

    def transform(self, values, params=None):
        """Function defining how the fitting parameters are used to transform the input vectors"""
        if params is None:
            params=self.optimised_values
        return self.parameter_to_values(values=values, params=params)

    def initialise_parameters(self):
        """Initialise starting simplex"""

        starting_values = []
        for i in range(self.total_params)
        v0 = 0.0    # 0th order - offset
        v1 = 1.0    # 1st order - scale
        self.starting_simplex = [    flex.double([v0,    v1]),
                                     flex.double([v0,    v1+0.1]),
                                     flex.double([v0+0.1,v1])      ]
        return [v0,v1]

    def _scale(self, values, params):
        v0,v1 = params
        out = v0 + values*v1
        return out

class MultiDatasetTLSFit(object):

    def __init__(self, u_residual, tls_selections, tls_params):

        self.u_residual = u_residual
        self.tls_sels = tls_selections
        self.tls_objs = tls_objects

    def extract(self, hierarchy):
        pass

############################################################################

def run(params):

    all_models = [CrystallographicModel.from_file(f).label(tag=os.path.basename(os.path.dirname(f))) for f in params.input.pdb]

    result = fit_adps_to_multiple_datasets(models)


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



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

    fitter = MultiDatasetTLSFitter(observed_uij=obs_uij, observed_xyz=obs_xyz
    fitter.run()


class MultiDatasetTLSFitter():

    def __init__(self, observed_uij, observed_xyz):

        self.observed_uij = observed_uij
        self.observed_xyz = observed_xyz

        print self.observed_uij.shape
        print self.observed_xyz.shape

        self.observed_u = observed_u
        self.n_dataset

        self.fitted_u_residual = flex.double

        if tls_selections is None:
            tls_selections = ['pepnames']

        self.models = models

        self.n_datasets = len(self.models)
        self.m_atoms    = len(self.models[0].hierarchy.atoms())

        self.i_tls_params     = None
        self.i_tls_selections = tls_selections

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



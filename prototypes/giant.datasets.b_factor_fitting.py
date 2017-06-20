#!/usr/bin/env ccp4-python

import os, sys, re, glob, shutil, copy, tempfile, gc
import math, re, time

import scipy.stats
import numpy, pandas

import libtbx.phil, libtbx.easy_mp
import iotbx.pdb
import mmtbx.tls.tools

from libtbx.utils import Sorry, Failure
from scitbx.array_family import flex
from scitbx import simplex, matrix

from bamboo.common.logs import Log
from bamboo.common.path import easy_directory

from giant.dataset import CrystallographicModel
from giant.structure.select import protein
from giant.structure.tls import uij_from_tls_vector_and_origin, extract_tls_from_pdb
from giant.structure.formatting import ShortLabeller

import matplotlib
matplotlib.interactive(False)
from matplotlib import pyplot
pyplot.switch_backend('agg')
pyplot.interactive(0)

numpy.set_printoptions(threshold=numpy.nan)

from IPython import embed

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
fitting {
    number_of_macro_cycles = 1
        .help = 'how many fitting cycles to run'
        .type = int
    refine_tls_parameters = True
        .help = "refine TLS parameters or only apply to input structures?"
        .type = bool
    refine_uij_residual = True
        .help = "refine UIJ residual or only apply to input structures?"
        .type = bool
    tls {
        number_of_models_per_group = 1
            .help = 'how many superposed TLS models should be fit to the same group of atoms?'
            .type = int
    }
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

def rms(vals, axis=None):
    return numpy.sqrt(numpy.mean(numpy.power(vals,2), axis=axis))

class MultiDatasetTLSParameterisation(object):

    def __init__(self, models, groups=None, n_cpu=1, n_tls=1, log=None):

        if log is None: log = Log(verbose=True)
        self.log = log

        self._n_cpu=n_cpu
        self._n_tls=n_tls

        self.models = models
        self.groups = groups
        self.fits = {}

        # Use the first hierarchy as the reference
        self.master_h = models[0].hierarchy.deep_copy()
        # Extract the atoms for each tls group
        self.atom_selections = dict([(g, self.master_h.atom_selection_cache().selection(g)) for g in self.groups])

        self.atom_selection_all = flex.bool(self.master_h.atoms().size(), False)
        for sel in self.atom_selections.values():
            self.atom_selection_all.set_selected(sel, True)

        self.validate_input()

    def validate_input(self):

        for m in self.models:
            assert self.master_h.is_similar_hierarchy(m.hierarchy)

    def blank_master_hierarchy(self):
        h = self.master_h.deep_copy()
        h.atoms().set_uij(flex.sym_mat3_double(h.atoms().size(), [0.0]*6))
        h.atoms().set_b(flex.double(h.atoms().size(), 0.0))
        return h

    def fit(self, n_cycles=1):

        for group in self.groups:

            self.log.heading('Parameterising Uijs for selection: {}'.format(group))

            # Get the selection for this group
            atom_sel = self.atom_selections[group]
            # Get all atoms for this group
            atoms = [m.hierarchy.atoms().select(atom_sel) for m in self.models]
            # Extract uij and xyz
            obs_uij = numpy.array([a.extract_uij() for a in atoms])
            obs_xyz = numpy.array([a.extract_xyz() for a in atoms])

            self.log.subheading('Fitting TLS models to data')
            fitter = MultiDatasetTLSFitter(observed_uij = obs_uij,
                                           observed_xyz = obs_xyz,
                                           n_tls = self._n_tls,
                                           n_cpu = self._n_cpu,
                                           log   = self.log)
            # Calculate scaling
            fitter.fit(n_cycles)

            self.fits[group] = fitter

        self.log.heading('Parameterisation complete')

    def write_fitted_uij_structures(self, out_dir='./'):
        """Write residual B-factors to master hierarchy."""

        easy_directory(out_dir)

        # ----------------------------------------------------
        # Apply the residual B-factors to the master h
        # ----------------------------------------------------
        self.log.subheading('Writing fitted residual atomic B-factors')
        h = self.blank_master_hierarchy()
        for group in self.groups:
            fit = self.fits[group]
            sel = self.atom_selections[group]
            uij = fit.fitted_uij_residual()
            h.atoms().select(sel).set_b(flex.double(EIGHT_PI_SQ*numpy.mean(uij[:,0:3],axis=1)))
            h.atoms().select(sel).set_uij(flex.sym_mat3_double(uij))
        h.write_pdb_file(os.path.join(out_dir, 'fit_uij_residual.pdb'))
        # ----------------------------------------------------
        # Apply the TLS B-factors to the master h
        # ----------------------------------------------------
        self.log.subheading('Writing fitted TLS atomic B-factors')
        h = self.blank_master_hierarchy()
        t = self.blank_master_hierarchy()
        l = self.blank_master_hierarchy()
        s = self.blank_master_hierarchy()
        tl = self.blank_master_hierarchy()
        for group in self.groups:
            fit = self.fits[group]
            sel = self.atom_selections[group]
            xyz = h.atoms().select(sel).extract_xyz()
            # ALL TLS contributions
            uij = fit.fitted_uij_tls(xyz=xyz)
            h.atoms().select(sel).set_b(flex.double(EIGHT_PI_SQ*numpy.mean(uij[:,0:3],axis=1)))
            h.atoms().select(sel).set_uij(flex.sym_mat3_double(uij))
            # Only T contributions
            uij = fit.fitted_uij_tls(xyz=xyz, t=True, l=False, s=False)
            t.atoms().select(sel).set_b(flex.double(EIGHT_PI_SQ*numpy.mean(uij[:,0:3],axis=1)))
            t.atoms().select(sel).set_uij(flex.sym_mat3_double(uij))
            # Only L contributions
            uij = fit.fitted_uij_tls(xyz=xyz, t=False, l=True, s=False)
            l.atoms().select(sel).set_b(flex.double(EIGHT_PI_SQ*numpy.mean(uij[:,0:3],axis=1)))
            l.atoms().select(sel).set_uij(flex.sym_mat3_double(uij))
            # Only S contributions
            uij = fit.fitted_uij_tls(xyz=xyz, t=False, l=False, s=True)
            s.atoms().select(sel).set_b(flex.double(EIGHT_PI_SQ*numpy.mean(uij[:,0:3],axis=1)))
            s.atoms().select(sel).set_uij(flex.sym_mat3_double(uij))
            # Combine T and L contributions
            tl.atoms().select(sel).set_b(t.atoms().select(sel).extract_b() + l.atoms().select(sel).extract_b())
            tl.atoms().select(sel).set_uij(t.atoms().select(sel).extract_uij() + l.atoms().select(sel).extract_uij())
        # Write out structures
        h.write_pdb_file(os.path.join(out_dir, 'fit_uij_tls.pdb'))
        t.write_pdb_file(os.path.join(out_dir, 'fit_uij_tls_t_only.pdb'))
        l.write_pdb_file(os.path.join(out_dir, 'fit_uij_tls_l_only.pdb'))
        s.write_pdb_file(os.path.join(out_dir, 'fit_uij_tls_s_only.pdb'))

        for chain_id in [c.id for c in h.chains()]:
            sel = h.atom_selection_cache().selection('chain {}'.format(chain_id))
            sel_h = h.select(sel)
            sel_t = t.select(sel)
            sel_l = l.select(sel)
            sel_s = s.select(sel)
            sel_tl = tl.select(sel)
            # Write graphs of individual and cumulative Tls contributions
            filename = os.path.join(out_dir, 'tls_contributions-chain-{}.png'.format(chain_id))
            x_vals   = numpy.array(range(len(list(sel_h.residue_groups()))))+1
            x_labels = ['0']+[ShortLabeller.format(rg) for rg in sel_h.residue_groups()]
            l_styles = ['r-', 'g-', 'b-', 'k-']
            l_labels = ['T', 'L', 'S', 'TLS']
            # Make plot
            fig, axes = pyplot.subplots(nrows=2, ncols=1, sharey=True)
            # Individual lines
            tls_vals = [numpy.array([numpy.mean(rg.atoms().extract_b()) for rg in sel_t.residue_groups()]),
                        numpy.array([numpy.mean(rg.atoms().extract_b()) for rg in sel_l.residue_groups()]),
                        numpy.array([numpy.mean(rg.atoms().extract_b()) for rg in sel_s.residue_groups()]),
                        numpy.array([numpy.mean(rg.atoms().extract_b()) for rg in sel_h.residue_groups()])]
            axes[0].set_title('Individual contributions')
            axes[0].set_xlabel('Residue')
            axes[0].set_ylabel('Isotropic-ised B')
            for i in range(4):
                axes[0].plot(x_vals, tls_vals[i], l_styles[i], label=l_labels[i])
            axes[0].legend()
            # Cumulative contributions
            tls_vals = [numpy.array([numpy.mean(rg.atoms().extract_b()) for rg in sel_t.residue_groups() ]),
                        numpy.array([numpy.mean(rg.atoms().extract_b()) for rg in sel_tl.residue_groups()]),
                        numpy.array([numpy.mean(rg.atoms().extract_b()) for rg in sel_h.residue_groups() ])]
            axes[1].set_title('Cumulative contributions')
            axes[1].set_xlabel('Residue')
            axes[1].set_ylabel('Isotropic-ised B')
            for i in range(3):
                axes[1].plot(x_vals, tls_vals[i], l_styles[i], label=''.join(l_labels[:i+1]))
            xticks = axes[1].get_xticks()
            axes[1].set_xticklabels([x_labels[i] if i<len(x_labels) else '' for i in xticks.astype(int)])
            labels = axes[1].get_xticklabels()
            pyplot.setp(labels, rotation=90)
            axes[1].legend()
            # Format and save
            pyplot.tight_layout()
            pyplot.savefig(filename)
            pyplot.close(fig)

    def write_parameterised_structures(self):
        """Write fitted B-factors to output structures."""

        # Extract the fitted output for each dataset
        self.log.subheading('Exporting parameterised B-factors')
        for i_mdl, mdl in enumerate(self.models):
            self.log('Writing structure for model: {}'.format(mdl.filename))
            h = mdl.hierarchy.deep_copy()
            for group in self.groups:
                sel = self.atom_selections[group]
                uij = self.fits[group].extract_fitted_uij(datasets=[i_mdl])[0]
                h.atoms().select(sel).set_uij(flex.sym_mat3_double(uij))
            h.write_pdb_file(mdl.filename.replace('.pdb', '.mda.pdb'))

    def write_summary_statistics_csv(self, out_dir='./'):
        """Add data to CSV and write"""

        easy_directory(out_dir)

        table = pandas.DataFrame(index=[m.tag for m in self.models])
        all_dataset_rmsds = numpy.array([numpy.concatenate([self.fits[g].uij_fit_obs_all_rmsds()[i] for g in self.groups]) for i in range(len(self.models))])
        medn_rmsds = numpy.median(all_dataset_rmsds, axis=1)
        table['median_rmsds'] = medn_rmsds
        mean_rmsds = numpy.mean(all_dataset_rmsds, axis=1)
        table['mean_rmsds'] = mean_rmsds
        table.to_csv(os.path.join(out_dir, 'dataset_scores.csv'))

        for group in self.groups:
            table = pandas.DataFrame(data    = self.fits[g].uij_fit_obs_all_rmsds(),
                                     index   = [m.tag for m in self.models],
                                     columns = [ShortLabeller.format(a) for a in self.master_h.atoms().select(self.atom_selections[group])])
            table.to_csv(os.path.join(out_dir, 'all_rmsd_scores.csv'))

    def write_summary_statistics_graphs_and_structures(self, out_dir='./'):
        """Write atom-by-atom and dataset-by-dataset graphs"""

        easy_directory(out_dir)

        # ----------------------------------------------------
        # Write the rmsds from the refined uijs
        # ----------------------------------------------------
        self.log.subheading('Calculating atom-by-atom rmsds to refined B-factors')
        h = self.blank_master_hierarchy()
        for i_g, group in enumerate(self.groups):
            fit = self.fits[group]
            sel_h = h.select(self.atom_selections[group])
            # Apply rmsd values to structure
            sel_h.atoms().set_b(flex.double(fit.uij_fit_obs_atom_averaged_rmsds()))
            sel_h.atoms().set_uij(flex.sym_mat3_double(fit.uij_fit_obs_atom_averaged_differences()))
            # Extract the atoms in this group as list for indexing
            sel_a = list(sel_h.atoms())
            rmsds = fit.uij_fit_obs_all_rmsds()
            # ------------------------
            # Write boxplots for all atoms in each dataset
            # ------------------------
            for i_d in range(0, len(self.models), 50):
                self.boxplot(filename=os.path.join(out_dir, 'dataset-by-dataset-rmsds-group-{:02d}-datasets-{:04d}-{:04d}.png'.format(i_g+1,i_d+1,i_d+50)),
                             y_vals=rmsds.T[:,i_d:i_d+50], x_labels=numpy.arange(i_d,min(i_d+50,len(self.models)))+1,
                             title='rmsd of fitted and refined B-factors (by dataset)',
                             x_lab='dataset', y_lab='rmsd', rotate_x_labels=True,
                             x_lim=(0,51), y_lim=(numpy.min(rmsds), numpy.max(rmsds)))
            # ------------------------
            # Write boxplots for all of the atoms in this group
            # ------------------------
            # Breaks between residues
            i_brk = numpy.array([float(i)+1.5 for i, (a1,a2) in enumerate(zip(sel_a,sel_a[1:])) if a1.parent().parent().resid()!=a2.parent().parent().resid()])
            for i_a in range(0, len(sel_a), 50):
                self.boxplot(filename=os.path.join(out_dir, 'atom-by-atom-rmsds-group-{:02d}-atoms-{:06d}-{:06d}.png'.format(i_g+1,i_a+1,i_a+50)),
                             y_vals=rmsds[:,i_a:i_a+50], x_labels=[ShortLabeller.format(a) for a in sel_a[i_a:i_a+50]],
                             title='rmsd of fitted and refined B-factors (by atom)',
                             x_lab='atom', y_lab='rmsd', rotate_x_labels=True,
                             x_lim=(0,51), y_lim=(numpy.min(rmsds), numpy.max(rmsds)), vlines=i_brk-i_a)
            # ------------------------
            # Write boxplots for atoms in each residue group separately
            # ------------------------
            for rg in sel_h.residue_groups():
                self.boxplot(filename=os.path.join(out_dir, 'atom-by-atom-rmsds-group-{:02d}-residue-{}.png'.format(i_g+1,ShortLabeller.format(rg))),
                             y_vals=rmsds[:, [sel_a.index(a) for a in rg.atoms()]],
                             x_labels=[ShortLabeller.format(a) for a in rg.atoms()],
                             title='average rmsd of fitted and refined B-factors: {}'.format(ShortLabeller.format(rg)),
                             x_lab='atom', y_lab='rmsd', rotate_x_labels=True,
                             y_lim=(numpy.min(rmsds), numpy.max(rmsds)))
        # ----------------------------------------------------
        # Write out structure
        # ----------------------------------------------------
        h.write_pdb_file(os.path.join(out_dir, 'fit_rmsd_uij.pdb'))
        sel_h = h.select(self.atom_selection_all)
        # ----------------------------------------------------
        # Write averaged rmsds for all atoms by residue
        # ----------------------------------------------------
        self.log.subheading('Calculating residue-by-residue rmsds to refined B-factors')
        all_bs = numpy.concatenate([list(rg.atoms().extract_b()) for rg in sel_h.residue_groups()])
        min_b, max_b = numpy.min(all_bs), numpy.max(all_bs)
        for sel_c in sel_h.chains():
            all_rgs = list(sel_c.residue_groups())
            for i_rg in range(0, len(all_rgs), 50):
                this_rgs = all_rgs[i_rg:i_rg+50]
                self.boxplot(filename=os.path.join(out_dir, 'residue-by-residue-rmsds-chain-{}-residues-{:04d}-{:04d}.png'.format(sel_c.id,i_rg+1,i_rg+50)),
                             y_vals=[numpy.array(rg.atoms().extract_b()) for rg in this_rgs],
                             x_labels=[ShortLabeller.format(rg) for rg in this_rgs],
                             title='averaged rmsds of fitted and refined B-factors (for each residue)',
                             x_lab='residue', y_lab='rmsd', rotate_x_labels=True,
                             x_lim=(0,51), y_lim=(numpy.min(all_bs), numpy.max(all_bs)))
        # ----------------------------------------------------
        # Write average rmsds for each dataset
        # ----------------------------------------------------
        self.log.subheading('Calculating dataset-by-dataset rmsds to refined B-factors')
        all_rmsds = numpy.array([self.fits[g].uij_fit_obs_all_rmsds() for g in self.groups])
        for i_m in range(0, len(self.models), 50):
            m_idxs = numpy.arange(i_m,min(i_m+50,len(self.models)))
            self.boxplot(filename=os.path.join(out_dir, 'dataset-by-dataset-rmsds-datasets-{:04d}-{:04d}.png'.format(i_m+1, i_m+50)),
                         y_vals=[all_rmsds[:,i].flatten() for i in m_idxs], x_labels=m_idxs+1,
                         title='rmsds for each dataset of fitted and refined B-factors',
                         x_lab='dataset', y_lab='rmsd', rotate_x_labels=True,
                         x_lim=(0,51), y_lim=(numpy.min(all_rmsds), numpy.max(all_rmsds)))
        # ----------------------------------------------------
        # Write distributions of amplitudes for tls models
        # ----------------------------------------------------
        for i_g, group in enumerate(self.groups):
            fit = self.fits[group]
            amps = fit.fitted_tls_amplitudes()
            x_vals = []; [[x_vals.append(amps[:,i_m,i_a]) for i_a in range(3)] for i_m in range(amps.shape[1])]
            self.histograms(filename=os.path.join(out_dir, 'tls-model-amplitudes-group-{}.png'.format(i_g+1)), x_vals=x_vals,
                            titles=numpy.concatenate(['T (group {a})-L (group {a})-S (group {a})'.format(a=i_m+1).split('-') for i_m in range(amps.shape[1])]),
                            x_labs=['']*numpy.product(amps.shape[1:]), rotate_x_labels=True, shape=amps.shape[1:], n_bins=30)

    def boxplot(self, filename, y_vals, x_labels, title, x_lab='x', y_lab='y', x_lim=None, y_lim=None, rotate_x_labels=True, vlines=None):

        self.log('Writing: {}'.format(filename))

        fig = pyplot.figure()
        pyplot.rc('font', family='monospace')
        pyplot.title(title)
        pyplot.boxplot(y_vals, labels=x_labels, showmeans=True)
        if (vlines is not None) and (y_lim is not None):
            for v in vlines:
                pyplot.vlines(v, y_lim[0], y_lim[1])
        pyplot.xlabel(x_lab)
        pyplot.ylabel(y_lab)
        pyplot.xlim(x_lim)
        pyplot.ylim(y_lim)
        if rotate_x_labels:
            locs, labels = pyplot.xticks()
            pyplot.setp(labels, rotation=90)
        pyplot.tight_layout()
        pyplot.savefig(filename)
        pyplot.close(fig)

        return

    def histograms(self, filename, x_vals, titles, x_labs, rotate_x_labels=True, shape=None, n_bins=30):

        self.log('Writing: {}'.format(filename))

        if shape is not None:
            nrow, ncol = shape
        else:
            nrow, ncol = (1,len(x_vals))

        fig, axes = pyplot.subplots(nrows=nrow, ncols=ncol, sharey=True)
        for i, axis in enumerate(axes.flatten()):
            axis.set_title(titles[i])
            axis.hist(x=x_vals[i], bins=n_bins)
            axis.set_xlabel(x_labs[0])
            axis.set_ylabel('Count')
            if rotate_x_labels:
                labels = axis.get_xticklabels()
                pyplot.setp(labels, rotation=90)
        pyplot.tight_layout()
        pyplot.savefig(filename)
        pyplot.close(fig)

        return

class MultiDatasetTLSFitter(object):

    _tls_weight = 1.0
    _amp_weight = 1.0
    _uij_weight = 1.0
    _ovr_weight = 1.0

    def __init__(self, observed_xyz, observed_uij, tls_params=None, n_tls=None, n_cpu=1, log=None):

        if log is None: log = Log(verbose=True)
        self.log = log

        self._test = False
        self._iter = 0

        self._n_cpu = n_cpu

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
        elif n_tls is not None:
            inp_tls = None
            num_tls = n_tls
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
        self.initial_tls_amplitudes = numpy.ones(self._n_prm_tls_amp)
        self.initial_uij_residuals  = numpy.zeros(self._n_prm_uij_res)
        # ---------------------------->
        # Output variables (initialise to None)
        # ---------------------------->
        #self.fitted_tls_parameters = None
        #self.fitted_tls_amplitudes = None
        #self.fitted_uij_residuals  = None

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
#        self._mask_dsets[:20] = True
        self._mask_atoms = numpy.ones(self._n_atm, dtype=bool)

        # ---------------------------->
        # Minimisation variables
        # ---------------------------->
        # Amount of expected variation for each data type (axes of simplex)
        self._del_tls_mdl = 0.25
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

    def _adopt(self, sub_vector, selection=None):
        """Insert a set of parameters into the complete parameter set"""
        if selection is None: selection = self._var_current_sel
        assert len(sub_vector) == numpy.sum(selection)
        self._var_current[selection] = sub_vector

    def _adopt_from_others(self, others):
        # Updating from others, so set own selection to none
        self._reset_current_selection()
        for n in others:
            # Add this selection to our selection
            self._var_current_sel += n._var_current_sel
            # Copy the values across
            self._adopt(sub_vector=n._var_current[n._var_current_sel], selection=n._var_current_sel)

    def _prep_for_mp(self):
        """Clear parts of the object that are not needed for optimisation after selection - save memory when pickling"""
        self.initial_tls_parameters = None
        self.initial_tls_amplitudes = None
        self.initial_uij_residuals  = None
        return self

    ################################################################################################
    ###
    ###                             Internal functions
    ###
    ################################################################################################

    def _blank_parameter_selection(self):
        return numpy.zeros(self._n_prm_total, dtype=bool)
    def _blank_dataset_selection(self):
        return numpy.zeros(self._n_dst, dtype=bool)
    def _blank_atom_selection(self):
        return numpy.zeros(self._n_atm, dtype=bool)

    def _reset_tls_amplitudes(self):
        self._var_current[self._sel_tls_amp] = 1.0
    def _reset_uij_residual(self):
        self._var_current[self._sel_uij_res] = 0.0
    def _reset_current_selection(self):
        self._var_current_sel = self._blank_parameter_selection()

    def _extract_parameters(self, vector=None):
        """Convert 1d vector into objects"""
        if vector is None: vector=self._var_current
        tls_mdls = vector[self._sel_tls_mdl].reshape((self._n_tls, 21)             )
        tls_amps = vector[self._sel_tls_amp].reshape((self._n_dst, self._n_tls, 3) )
        uij_resl = vector[self._sel_uij_res].reshape((self._n_atm, 6)              )
        return (tls_mdls,tls_amps,uij_resl)

    def _unpack_tls_parameters(self, vals):
        return vals[0:6], vals[6:12], vals[12:21]

    def _expand_tls_amplitudes(self, tls_amps):
        """Convert 3-element vector into 21 element vector for TLS multiplication"""
        cur_n_dst = len(tls_amps)
        assert tls_amps.shape == (cur_n_dst,self._n_tls,3)
        t_amps = numpy.repeat(tls_amps[:,:,0], 6, axis=1).reshape((cur_n_dst, self._n_tls, 6))
        l_amps = numpy.repeat(tls_amps[:,:,1], 6, axis=1).reshape((cur_n_dst, self._n_tls, 6))
        s_amps = numpy.repeat(tls_amps[:,:,2], 9, axis=1).reshape((cur_n_dst, self._n_tls, 9))
        exp_tls_amps = numpy.concatenate([t_amps, l_amps, s_amps], axis=2)
        assert exp_tls_amps.shape == (cur_n_dst,self._n_tls,21)
        return exp_tls_amps

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

    def _get_simplex(self, selection):
        starting_values = self._var_current[selection]
        starting_deltas = self._del_simplex[selection]
        starting_simplex = numpy.repeat([starting_values], len(starting_values)+1, axis=0)
        for i in range(len(starting_values)):
            starting_simplex[i+1][i] += starting_deltas[i]
        return starting_simplex

    def _get_dataset_and_atom_selection_bool(self, datasets=[], atoms=[]):
        datasets_bool = self._blank_dataset_selection()
        if (datasets is None) or (len(datasets) == 0):  datasets_bool += True   # select all
        else:                                           datasets_bool.put(datasets, True)
        atoms_bool = self._blank_atom_selection()
        if (atoms is None) or (len(atoms) == 0):        atoms_bool += True      # select all
        else:                                           atoms_bool.put(atoms, True)
        return datasets_bool, atoms_bool

    def _selection_filter(self, datasets=None, atoms=None):
        """If either datasets or atoms is None, returns slice over all datasets/atoms"""
        if datasets is None: datasets = range(self._n_dst)
        if atoms    is None: atoms    = range(self._n_atm)
        return datasets, atoms

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
        # Calculate penalties
        pen = self._all_penalties(vector=self._var_current, uij_fit=uij_fit, uij_obs=uij_obs)
        return rmsd+pen

    ################################################################################################
    ###
    ###                             Constraints / Restraints
    ###
    ################################################################################################

    def set_penalty_weights(self, tls_weight=None, amp_weight=None, uij_weight=None, ovr_weight=None):
        """Set penalties for parameters to be invalid"""
        if tls_weight is not None: self._tls_weight = tls_weight
        if amp_weight is not None: self._amp_weight = amp_weight
        if uij_weight is not None: self._uij_weight = uij_weight
        if ovr_weight is not None: self._ovr_weight = ovr_weight

    def _all_penalties(self, vector, uij_fit, uij_obs):
        tls_mdl, tls_amp, uij_res = self._extract_parameters(vector=vector)
        tls_penalties = [self._tls_penalty(values=v) for v in tls_mdl]
        amp_penalties = [self._amp_penalty(values=v) for v in tls_amp]
        uij_penalties = [self._uij_penalty(values=v) for v in uij_res]
        fit_penalties = []; [fit_penalties.extend([self._uij_penalty(values=vv) for vv in v]) for v in uij_fit]
        ovr_penalties = []; [ovr_penalties.extend([self._ovr_penalty(*vv) for vv in zip(*v)]) for v in zip(uij_fit,uij_obs)]
        return numpy.sum(tls_penalties+amp_penalties+uij_penalties+fit_penalties+ovr_penalties)

    def _tls_penalty(self, values):
        assert len(values) == 21
        t,l,s = self._unpack_tls_parameters(vals=values)
        t_eig_values = numpy.linalg.eigvals(matrix.sym(sym_mat3=t).as_numpy_array())
        t_penalty = numpy.sum(t_eig_values<0.0)
        l_eig_values = numpy.linalg.eigvals(matrix.sym(sym_mat3=l).as_numpy_array())
        l_penalty = numpy.sum(l_eig_values<0.0)
        return self._tls_weight*numpy.sum([t_penalty, l_penalty])

    def _amp_penalty(self, values):
        return self._amp_weight*numpy.sum(values<0.0)

    def _uij_penalty(self, values):
        assert len(values) == 6
        eig_values = numpy.linalg.eigvals(matrix.sym(sym_mat3=values).as_numpy_array())
        return self._uij_weight*numpy.sum(eig_values<0.0)

    def _ovr_penalty(self, fit, obs):
        """Add penalty for having fitted B-factors greater than observed"""
        eig_values_fit = numpy.linalg.eigvals(matrix.sym(sym_mat3=fit).as_numpy_array())
        eig_values_obs = numpy.linalg.eigvals(matrix.sym(sym_mat3=obs).as_numpy_array())
        return self._ovr_weight*(max(eig_values_fit)>max(eig_values_obs))

    ################################################################################################
    ###
    ###                             Update functions
    ###
    ################################################################################################

    def _update_atom_masks(self):
        _, _, uij_res = self._extract_parameters()
        uij_max = numpy.max(numpy.abs(uij_res[:,:3]),axis=1)
        thresh = numpy.percentile(uij_max, 90)
        self._mask_atoms *= (uij_max < thresh)

    def _update_dataset_masks(self):
        thresh = 3.0
        d_rmsd = self.uij_fit_obs_dataset_averaged_rmsds()
        zscore = scipy.stats.zscore(d_rmsd)
        self._mask_dsets *= (zscore < thresh)

    ################################################################################################
    ###
    ###                             Main methods for running/optimisation
    ###
    ################################################################################################

    def fit(self, n_cycles=3):
        """Run macro-cycles of parameter optimisation"""
        for i_cyc in range(n_cycles):

            #########################################
            self.log.heading('Macrocycle {}'.format(i_cyc+1), spacer=True)
            #########################################
            # Remove outlier atoms from minimisation
            #########################################
            if i_cyc > 0:
                self.log.subheading('Removing atoms with high residual uij from TLS optimisation')
                self._update_atom_masks()
                self.log('Using {} atoms for TLS optimisation'.format(numpy.sum(self._mask_atoms)))

                self.log.subheading('Removing datasets with high fit rmsds from TLS optimisation')
                self._update_dataset_masks()
                self.log('Using {} datasets for TLS optimisation'.format(numpy.sum(self._mask_dsets)))


            #########################################
            #self.log('Resetting tls amplitudes')
            #self._reset_tls_amplitudes()
            self.log('Resetting uij residuals')
            self._reset_uij_residual()

            #########################################
            # Moving selection for amplitude optimisation
            prev_opt = self._blank_parameter_selection()
            for i_tls in range(self._n_tls):
                for c_name, c_sel in [('T',self._sel_t), ('L',self._sel_l), ('S',self._sel_s)]:
                    #########################################
                    # Add these parameters to the parameters for optimisation
                    prev_opt += self._sel_tls[i_tls]*c_sel
                    #########################################
                    # Update penalty weight for overfitting
                    if (c_name=='T') and (i_tls==0): self.set_penalty_weights(ovr_weight=1e6)
                    else:                            self.set_penalty_weights(ovr_weight=0.0)
                    #########################################
                    self.log.heading('Optimising {} parameters for TLS model {}'.format(c_name, i_tls+1))
                    self._select(parameter_selection = self._sel_tls[i_tls]*self._sel_tls_mdl*c_sel,
                                 datasets = None,
                                 atoms    = None)
                    self._optimise()
                    self.optimisation_summary(False)
                    #########################################
                    self.set_penalty_weights(ovr_weight=0.1)
                    #########################################
                    self.log.heading('Optimising all amplitudes for all TLS models for all datasets'.format(i_tls+1))
                    proc_args = []
                    for i_dst in range(self._n_dst):
                        n = self.copy()
                        n._select(parameter_selection = self._sel_dst[i_dst]*prev_opt,
                                  datasets = [i_dst],
                                  atoms    = None)
                        proc_args.append(n._prep_for_mp())
                    self._adopt_from_others(libtbx.easy_mp.pool_map(processes=self._n_cpu, func=proc_wrapper, args=proc_args))
                    self.optimisation_summary(False)
                #########################################
                self.set_penalty_weights(ovr_weight=10.0)
                #########################################
                self.log.heading('Optimising all parameters and amplitudes for TLS model {}'.format(i_tls+1))
                self._select(parameter_selection = self._sel_tls[i_tls],
                             datasets = None,
                             atoms    = None)
                self._optimise()
                self.optimisation_summary(False)
            #########################################
            self.set_penalty_weights(ovr_weight=0.0)
            #########################################
            self.log.heading('Optimising residual Uijs')
            proc_args = []
            for i in range(self._n_atm):
                n = self.copy()
                n._select(parameter_selection = self._sel_atm[i],
                          datasets = None,
                          atoms    = [i])
                proc_args.append(n._prep_for_mp())
            self._adopt_from_others(libtbx.easy_mp.pool_map(processes=self._n_cpu, func=proc_wrapper, args=proc_args))
            self.optimisation_summary(False)
            #########################################
            self.log.subheading('End of macrocycle {}'.format(i_cyc+1))
            self.optimisation_summary()

        return self

    def extract_observed_xyz(self, datasets=None, atoms=None):
        datasets, atoms = self._selection_filter(datasets=datasets, atoms=atoms)
        return self.observed_xyz[datasets][:,atoms]

    def extract_observed_uij(self, datasets=None, atoms=None):
        datasets, atoms = self._selection_filter(datasets=datasets, atoms=atoms)
        return self.observed_uij[datasets][:,atoms]

    def tls_uij(self, xyz, tls_vectors, origin):
        """Convert a set of parameter vectors to a set of uijs"""
        return numpy.sum([uij_from_tls_vector_and_origin(xyz=xyz, tls_vector=v, origin=origin) for v in tls_vectors], axis=0)

    def extract_fitted_uij(self, datasets=None, atoms=None):
        """Extract total fitted uijs for a subset of datasets or atoms"""
        datasets, atoms = self._selection_filter(datasets=datasets, atoms=atoms)
        # Extract the optimised values
        tls_p, tls_a, uij_r = self._extract_parameters()
        # Extract relevant coordinates
        xyz = self.extract_observed_xyz(datasets=datasets, atoms=atoms)
        assert xyz.shape == (len(datasets), len(atoms), 3)
        # Extract only those that we're interested in
        tls_a = tls_a[datasets]
        assert tls_a.shape == (len(datasets), self._n_tls, 3)
        uij_r = uij_r[atoms]
        assert uij_r.shape == (len(atoms), 6)
        # Multiply tls amplitudes and models
        tls_f = self._expand_tls_amplitudes(tls_amps=tls_a) * tls_p
        assert tls_f.shape == (len(datasets), self._n_tls, 21)

        assert len(xyz) == len(tls_f)
        # Calculate the tls component of uij
        uij_fit = numpy.array([self.tls_uij(xyz=xyz[i], tls_vectors=tls_f[i], origin=self.observed_com) for i in range(len(datasets))])
        assert uij_fit.shape == (len(datasets), len(atoms), 6)
        # Add the residual uij
        uij_fit += uij_r

        return uij_fit

    def fitted_tls_models(self):
        v,_,_ = self._extract_parameters()
        return v
    def fitted_tls_amplitudes(self, datasets=None):
        _,v,_ = self._extract_parameters()
        if datasets is not None: return v[datasets]
        return v
    def fitted_uij_residual(self, atoms=None):
        _,_,v = self._extract_parameters()
        if atoms is not None: return v[atoms]
        return v
    def fitted_uij_tls(self, xyz, models=None, t=True, l=True, s=True):
        tls = self.fitted_tls_models()
        if models is not None: tls=tls[models]
        if t is not True: tls[:,00:06] = 0.0
        if l is not True: tls[:,06:12] = 0.0
        if s is not True: tls[:,12:21] = 0.0
        return self.tls_uij(xyz=xyz, tls_vectors=tls, origin=self.observed_com)

    ################################################################################################
    ###
    ###                             Summaries
    ###
    ################################################################################################

    def input_summary(self):
        """Print the number of parameters/input data"""
        self.log.subheading('Input summary')
        self.log('input uij parameters: {}'.format(self.observed_uij.shape))
        self.log('input xyz parameters: {}'.format(self.observed_xyz.shape))
        self.log('Centre of mass: {}'.format(self.observed_com))
        self.log.bar()
        self.log('Number of parameters for TLS fitting: {}'.format(self._n_prm_total))
        self.log('Number of observations: {}'.format(self._n_obs_total))
        self.log('Data/parameter ratio is {:.3f}'.format(self._n_obs_total*1.0/self._n_prm_total))
        self.log.bar()
        self.log('Number of datasets:   {}'.format(self._n_dst))
        self.log('Number of atoms:      {}'.format(self._n_atm))
        self.log('Number of TLS models: {}'.format(self._n_tls))
        self.log.bar()

    def optimisation_summary(self, full=True):
        """Print the fitted parameters"""
        tls_prms, tls_amps, res_uij = self._extract_parameters()

        self.log.subheading('Optimisation Summary')

        tls_refined = numpy.sum(self._sel_tls_mdl*self._var_current_sel)
        amp_refined = numpy.sum(self._sel_tls_amp*self._var_current_sel)
        uij_refined = numpy.sum(self._sel_uij_res*self._var_current_sel)

        self.log('Optimisation used data from {} datasets'.format(len(self._cur_datasets)))
        self.log('Optimisation used data for {} atoms'.format(len(self._cur_atoms)))
        self.log.bar()
        self.log('Optimised {} parameters'.format(numpy.sum(self._var_current_sel)))
        self.log('... {} TLS parameters,'.format(tls_refined))
        self.log('... {} TLS amplitudes,'.format(amp_refined))
        self.log('... {} uij parameters.'.format(uij_refined))
        self.log.bar()
        self.log('Optimisation RMSD: {}'.format(self.optimisation_rmsd))
        self.log.bar()
        self.log('')
        if tls_refined or full:
            self.parameter_summary(tls_params=tls_prms)
        if amp_refined or full:
            self.parameter_summary(tls_amplitudes=tls_amps)
        if uij_refined or full:
            self.parameter_summary(uij_residual=res_uij)

    def parameter_summary(self, tls_params=None, tls_amplitudes=None, uij_residual=None):
        if tls_params is tls_amplitudes is uij_residual is None:
            tls_params, tls_amplitudes, uij_residual = self._extract_parameters(vector=self._var_current)
        self.log.bar()
        if tls_params is not None:
            self.log('TLS Parameters:')
            for i, vals in enumerate(tls_params):
                self.log.bar()
                self.log('Model {:2}'.format(i+1))
                t,l,s = self._unpack_tls_parameters(vals=vals)
                self.log('T: '+', '.join(['{:8.3f}'.format(v) for v in t]))
                self.log('L: '+', '.join(['{:8.3f}'.format(v) for v in l]))
                self.log('S: '+', '.join(['{:8.3f}'.format(v) for v in s]))
            self.log.bar()
        if tls_amplitudes is not None:
            self.log('TLS Amplitudes:')
            for i, vals in enumerate(tls_amplitudes):
                self.log.bar()
                self.log('Dataset {:4}'.format(i+1))
                for i, amps in enumerate(vals):
                    self.log('Model {:2}:'.format(i+1)+' {:8.3f} (T) {:8.3f} (L) {:8.3f} (S)'.format(*amps))
            self.log.bar()
        if uij_residual is not None:
            self.log('Residual Uijs:')
            for i, vals in enumerate(uij_residual):
                self.log('Atom {:5}: '.format(i+1) + ', '.join(['{:8.3f}'.format(v) for v in vals]))
            self.log.bar()

    def uij_fit_obs_differences(self):
        uij_obs = self.extract_observed_uij()
        uij_fit = self.extract_fitted_uij()
        return uij_obs-uij_fit

    def uij_fit_obs_all_rmsds(self):
        return rms(self.uij_fit_obs_differences(), axis=2)

    def uij_fit_obs_dataset_averaged_rmsds(self):
        return numpy.mean(self.uij_fit_obs_all_rmsds(), axis=1)

    def uij_fit_obs_dataset_averaged_differences(self):
        return numpy.mean(self.uij_fit_obs_differences(), axis=1)

    def uij_fit_obs_atom_averaged_rmsds(self):
        return numpy.mean(self.uij_fit_obs_all_rmsds(), axis=0)

    def uij_fit_obs_atom_averaged_differences(self):
        return numpy.mean(self.uij_fit_obs_differences(), axis=0)


############################################################################

def run(params):

    log = Log(params.output.log_file, verbose=True)

    log('Building model list')
    if params.input.labelling == 'foldername':
        label_func = lambda f: os.path.basename(os.path.dirname(f))
    else:
        label_func = lambda f: os.path.basename(os.path.splitext(f)[0])
    models = [CrystallographicModel.from_file(f).label(tag=label_func(f)) for f in params.input.pdb]#[:10]

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


    p = MultiDatasetTLSParameterisation(models = models,
                                        groups = params.input.tls_group,
                                        n_tls  = params.fitting.tls.number_of_models_per_group,
                                        n_cpu  = params.settings.cpus,
                                        log = log)
    p.fit(params.fitting.number_of_macro_cycles)
    p.write_summary_statistics_csv(out_dir=params.output.out_dir)
    p.write_summary_statistics_graphs_and_structures(out_dir=params.output.out_dir)
    p.write_fitted_uij_structures(out_dir=params.output.out_dir)
    p.write_parameterised_structures()

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



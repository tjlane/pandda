#!/usr/bin/env ccp4-python

import os, sys, copy
import math

import scipy.stats
import numpy, pandas, json

import libtbx.phil, libtbx.easy_mp
import iotbx.pdb
import mmtbx.tls.tools

from libtbx.utils import Sorry, Failure
from scitbx.array_family import flex
from scitbx import simplex, matrix, linalg

from bamboo.common.logs import Log
from bamboo.common.path import easy_directory, rel_symlink
from bamboo.common.command import CommandManager

from giant.dataset import CrystallographicModel
from giant.structure.select import protein
from giant.structure.tls import uij_from_tls_vector_and_origin, extract_tls_from_pdb
from giant.structure.formatting import ShortLabeller
from giant.xray.crystal import CrystalSummary
from giant.xray.refine import refine_phenix

from giant.jiffies import multi_table_ones

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
    Fit a consensus B-factor model to a series of datasets
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
}
fitting {
    optimise = *tls_models *uij_residuals
        .help = "Which parameters should be optimised across the datasets?"
        .type = choice(multi=True)
    tls {
        number_of_models_per_group = 1
            .help = 'how many superposed TLS models should be fit to the same group of atoms?'
            .type = int
        number_of_datasets_for_optimisation = None
            .help = 'how many datasets should be used for optimising the TLS parameters?'
            .type = int
        resolution_limit_for_fitting = 1.8
            .help = 'resolution limit for dataset to be used for TLS optimisation'
            .type = float
    }
    number_of_macro_cycles = 2
        .help = 'how many fitting cycles to run'
        .type = int
    number_of_micro_cycles = 2
        .help = 'how many fitting cycles to run (for each T-L-S component)'
        .type = int
}
refine {
    refine_output_structures = True
        .help = "Refine the structures after fitting (coordinates and occupancies)"
        .type = bool
}
table_ones_options {
    include scope giant.jiffies.multi_table_ones.options_phil
}
settings {
    cpus = 48
        .type = int
        .multiple = False
}
""", process_includes=True)

############################################################################

def wrapper_optimise(arg):
    arg._optimise(verbose=False)
    return arg

def wrapper_run(arg):
    return arg.run()

def rms(vals, axis=None):
    return numpy.sqrt(numpy.mean(numpy.power(vals,2), axis=axis))

class MultiDatasetBFactorParameterisation(object):

    def __init__(self, models, groups, params, n_cpu=1, log=None):

        if log is None: log = Log(verbose=True)
        self.log = log

        self.params = params

        self.out_dir = params.output.out_dir

        self._n_cpu = n_cpu
        self._n_tls = params.fitting.tls.number_of_models_per_group
        self._n_opt = params.fitting.tls.number_of_datasets_for_optimisation

        self._allow_isotropic = True

        self._optimise_tls_models = ('tls_models' in params.fitting.optimise)
        self._optimise_uij_resdls = ('uij_residuals' in params.fitting.optimise)

        self._opt_datasets_res_limit = params.fitting.tls.resolution_limit_for_fitting
        self._opt_datasets_selection = []

        self.models = models
        self.groups = groups
        self.fits = {}

        # Misc files
        self.cifs = None

        # Validate and add output paths, etc.
        self._init_input_models()
        self._init_atom_selections()

        self.table = None
        self.table_one_csv_input   = None
        self.table_one_csv_fitted  = None
        self.table_one_csv_refined = None

    def _init_input_models(self):

        # Use the first hierarchy as the reference
        self.master_h = self.models[0].hierarchy.deep_copy()

        s_dir = easy_directory(os.path.join(self.out_dir, 'parameterised_structures'))

        for i_m, m in enumerate(self.models):
            assert self.master_h.is_similar_hierarchy(m.hierarchy)

            m.i_pdb = m.filename
            m.i_mtz = m.i_pdb.replace('.pdb', '.mtz')
            m.o_pdb = None
            m.o_mtz = None
            m.r_pdb = None
            m.r_mtz = None

            assert os.path.exists(m.i_pdb), 'PDB does not exist: {}'.format(m.i_pdb)
            assert os.path.exists(m.i_mtz), 'MTZ does not exist: {}'.format(m.i_mtz)

            cs = CrystalSummary.from_mtz(m.i_mtz)
            if cs.high_res < self._opt_datasets_res_limit:
                self._opt_datasets_selection.append(i_m)

        if len(self._opt_datasets_selection) == 0:
            raise Exception('No datasets above resolution cutoff: {}'.format(self._opt_datasets_res_limit))
        if self._n_opt is not None:
            self._opt_datasets_selection = self._opt_datasets_selection[:self._n_opt]

    def _init_atom_selections(self):
        # Extract the atoms for each tls group
        self.atom_selections = dict([(g, self.master_h.atom_selection_cache().selection(g)) for g in self.groups])
        # Find which atoms in the structure are part of any group
        self.atom_selection_all = flex.bool(self.master_h.atoms().size(), False)
        for sel in self.atom_selections.values():
            self.atom_selection_all.set_selected(sel, True)

    def blank_master_hierarchy(self):
        h = self.master_h.deep_copy()
        h.atoms().set_uij(flex.sym_mat3_double(h.atoms().size(), [0.0]*6))
        h.atoms().set_b(flex.double(h.atoms().size(), 0.0))
        return h

    def iterate_groups(self):
        return ((group, self.fits.get(group,None), self.atom_selections[group]) for group in self.groups)
    def enumerate_groups(self):
        return enumerate(self.iterate_groups())

    def fit(self):
        """Optimise the TLS+residual model against the input data"""

        n_macro_cycles = self.params.fitting.number_of_macro_cycles
        n_micro_cycles = self.params.fitting.number_of_micro_cycles

        self.log.heading('Fitting B-factors', spacer=True)
        self.log.subheading('TLS Groups')
        for i, group in enumerate(self.groups):
            self.log('Group {}: {}'.format(i+1, group))
        self.log.subheading('Settings')
        self.log('Macro-cycles: {}'.format(n_macro_cycles))
        self.log('Micro-cycles: {}'.format(n_micro_cycles))

        for group, _, sel in self.iterate_groups():

            self.log.heading('Parameterising Uijs for selection: {}'.format(group))

            # Get all atoms for this group
            atoms = [m.hierarchy.atoms().select(sel) for m in self.models]
            # Extract uij and xyz
            obs_uij = numpy.array([a.extract_uij() for a in atoms])
            obs_xyz = numpy.array([a.extract_xyz() for a in atoms])

            # Check all uijs are present
            if (obs_uij==-1).all() and (self._allow_isotropic is True):
                self.log('All atoms are missing anisotropic displacement parameters -- using the isotropic parameters instead')
                obs_uij = numpy.zeros(obs_uij.shape)
                obs_uij[:,:,0] = numpy.array([a.extract_b()/EIGHT_PI_SQ for a in atoms])
                obs_uij[:,:,2] = obs_uij[:,:,1] = obs_uij[:,:,0]
            elif (obs_uij==-1).any():
                raise Failure('Some atoms in {} do not have anisotropically-refined B-factors'.format(group))

            fitter = MultiDatasetTLSFitter(observed_uij = obs_uij,
                                           observed_xyz = obs_xyz,
                                           n_tls = self._n_tls,
                                           n_cpu = self._n_cpu,
                                           optimisation_datasets = self._opt_datasets_selection,
                                           log   = self.log)
            # Calculate scaling
            fitter.fit(n_macro_cycles, n_micro_cycles)

            self.fits[group] = fitter

        self.log.heading('Parameterisation complete')

    def write_output(self):
        """Write output and summaries"""

        self.log.heading('Writing parameterisation summaries')
        self.write_parameterisation_summary(out_dir=self.out_dir)
        self.write_parameterisation_models(out_dir=self.out_dir)
        self.write_parameterisation_analysis(out_dir=os.path.join(self.out_dir,'parameterisation_quality'))

        self.log.heading('Writing TLS-subtracted models')
        self.write_tls_subtracted_models(out_dir=os.path.join(self.out_dir,'tls_subtracted_models'))

        self.log.heading('Outputting fitted structures for each dataset')
        self.write_fitted_dataset_models(out_dir=os.path.join(self.out_dir, 'fitted_structures'))
        if self.params.refine.refine_output_structures:
            self.refine_fitted_dataset_models()

        self.log.heading('Generating Table Ones for all structures')
        self.generate_fitted_table_ones(out_dir=os.path.join(self.out_dir, 'table_ones'))

        self.log.heading('Writing output csvs')
        self.write_combined_csv(out_dir=self.out_dir)

    def write_parameterisation_summary(self, out_dir='.'):
        """Write the TLS models and amplitudes to files"""

        out_dir = easy_directory(out_dir)

        # Output objects
        amp_table = pandas.DataFrame(index=[mdl.tag for mdl in self.models])
        tls_json = {}

        for i_g, (group, fit, sel) in self.enumerate_groups():
            # ----------------------------------------------------
            # Add TLS models to json
            # -----------i-----------------------------------------
            tlss = fit.fitted_tls_models()
            for i_tls in range(self._n_tls):
                d = tls_json.setdefault(group, {}).setdefault('TLS model {}'.format(i_tls+1), {})
                d['T'] = tlss[i_tls][00:06].tolist()
                d['L'] = tlss[i_tls][06:12].tolist()
                d['S'] = tlss[i_tls][12:21].tolist()
                d['COM'] = None
            # ----------------------------------------------------
            # TLS ampltidues for all datasets
            # ----------------------------------------------------
            amps = fit.fitted_tls_amplitudes()
            # Write histograms of amplitudes
            x_vals = []; [[x_vals.append(amps[:,i_m,i_a]) for i_a in range(3)] for i_m in range(amps.shape[1])]
            self.histograms(filename=os.path.join(out_dir, 'tls-model-amplitudes-group-{}.png'.format(i_g+1)), x_vals=x_vals,
                            titles=numpy.concatenate(['T (group {a})-L (group {a})-S (group {a})'.format(a=i_m+1).split('-') for i_m in range(amps.shape[1])]),
                            x_labs=['']*numpy.product(amps.shape[1:]), rotate_x_labels=True, shape=amps.shape[1:], n_bins=30)
            # Add to table of amplitudes
            for i_tls in range(self._n_tls):
                amp_table['T'+str(i_tls+1)] = amps[:,i_tls,0]
                amp_table['L'+str(i_tls+1)] = amps[:,i_tls,1]
                amp_table['S'+str(i_tls+1)] = amps[:,i_tls,2]
        # ----------------------------------------------------
        # Write output
        # ----------------------------------------------------
        # Write amplitude CSV
        filename = os.path.join(out_dir, 'tls_amplitudes.csv')
        self.log('Writing: {}'.format(filename))
        amp_table.to_csv(filename)
        # Write model JSON
        filename = os.path.join(out_dir, 'tls_models.json')
        self.log('Writing: {}'.format(filename))
        with open(filename, 'w') as fh: fh.write(json.dumps(tls_json, indent=True))

    def write_parameterisation_models(self, out_dir='.'):
        """Write residual B-factors to master hierarchy."""

        easy_directory(out_dir)

        # ----------------------------------------------------
        # Apply the residual B-factors to the master h
        # ----------------------------------------------------
        self.log.subheading('Writing consensus residual atomic B-factors')
        h = self.blank_master_hierarchy()
        for group, fit, sel in self.iterate_groups():
            uij = fit.fitted_uij_residual()
            h.atoms().select(sel).set_b(flex.double(EIGHT_PI_SQ*numpy.mean(uij[:,0:3],axis=1)))
            h.atoms().select(sel).set_uij(flex.sym_mat3_double(uij))
        h.write_pdb_file(os.path.join(out_dir, 'consensus_uij_residual.pdb'))
        # ----------------------------------------------------
        # Apply the TLS B-factors to the master h
        # ----------------------------------------------------
        self.log.subheading('Writing consensus TLS components of atomic B-factors')
        # Combined TLS + individual hierarchies
        h = self.blank_master_hierarchy()
        t = self.blank_master_hierarchy()
        l = self.blank_master_hierarchy()
        s = self.blank_master_hierarchy()
        # Add components from each group
        for i_g, (group, fit, sel) in self.enumerate_groups():
            xyz = h.atoms().select(sel).extract_xyz()
            # Which atoms are in this group anyway?
            g = self.blank_master_hierarchy()
            g.atoms().select(sel).set_b(flex.double(g.atoms().select(sel).size(), 10))
            g.write_pdb_file(os.path.join(out_dir, 'tls-group-{:02d}-atoms.pdb'.format(i_g+1)))
            # ALL TLS contributions
            uij = fit.average_fitted_uij_tls(xyz=xyz)
            h.atoms().select(sel).set_b(flex.double(EIGHT_PI_SQ*numpy.mean(uij[:,0:3],axis=1)))
            h.atoms().select(sel).set_uij(flex.sym_mat3_double(uij))
            # Only T contributions
            uij = fit.average_fitted_uij_tls(xyz=xyz, t=True, l=False, s=False)
            t.atoms().select(sel).set_b(flex.double(EIGHT_PI_SQ*numpy.mean(uij[:,0:3],axis=1)))
            t.atoms().select(sel).set_uij(flex.sym_mat3_double(uij))
            # Only L contributions
            uij = fit.average_fitted_uij_tls(xyz=xyz, t=False, l=True, s=False)
            l.atoms().select(sel).set_b(flex.double(EIGHT_PI_SQ*numpy.mean(uij[:,0:3],axis=1)))
            l.atoms().select(sel).set_uij(flex.sym_mat3_double(uij))
            # Only S contributions
            uij = fit.average_fitted_uij_tls(xyz=xyz, t=False, l=False, s=True)
            s.atoms().select(sel).set_b(flex.double(EIGHT_PI_SQ*numpy.mean(uij[:,0:3],axis=1)))
            s.atoms().select(sel).set_uij(flex.sym_mat3_double(uij))
        # Write out structures
        h.write_pdb_file(os.path.join(out_dir, 'consensus_uij_tls.pdb'))
        t.write_pdb_file(os.path.join(out_dir, 'consensus_uij_tls_t_only.pdb'))
        l.write_pdb_file(os.path.join(out_dir, 'consensus_uij_tls_l_only.pdb'))
        s.write_pdb_file(os.path.join(out_dir, 'consensus_uij_tls_s_only.pdb'))
        # ----------------------------------------------------
        # Write a chain-by-chain summary of the TLS components
        # ----------------------------------------------------
        # Find the edges of each TLS group
        b = self.blank_master_hierarchy()
        for i_g, (group, fit, sel) in self.enumerate_groups():
            # Where are the boundaries of the group?
            sel_int = numpy.array(sel, dtype=int)
            boundaries = numpy.array(b.atoms()[:-1].extract_b(), dtype=bool) + numpy.array((sel_int[:-1]*(1-sel_int[1:]))+(sel_int[1:]*(1-sel_int[:-1])), dtype=bool)
            b.atoms()[:-1].set_b(flex.double(boundaries.tolist()))
        # Iterate through chains
        for chain_id in [c.id for c in h.chains()]:
            sel = h.atom_selection_cache().selection('chain {}'.format(chain_id))
            sel_h = h.select(sel)
            sel_t = t.select(sel)
            sel_l = l.select(sel)
            sel_s = s.select(sel)
            # Write graphs of individual and cumulative Tls contributions
            filename = os.path.join(out_dir, 'tls-contributions-chain-{}.png'.format(chain_id))
            x_vals   = numpy.array(range(len(list(sel_h.residue_groups()))))+1
            x_labels = ['']+[ShortLabeller.format(rg) for rg in sel_h.residue_groups()]
            l_styles = ['ro-', 'go-', 'bo-', 'ko-']
            l_labels = ['T', 'L', 'S', 'TLS']
            # Make plot
            fig, axes = pyplot.subplots(nrows=2, ncols=1, sharey=True)
            # Individual lines
            tls_vals = [numpy.array([numpy.mean(rg.atoms().extract_b()) for rg in sel_t.residue_groups()]),
                        numpy.array([numpy.mean(rg.atoms().extract_b()) for rg in sel_l.residue_groups()]),
                        numpy.array([numpy.mean(rg.atoms().extract_b()) for rg in sel_s.residue_groups()]),
                        numpy.array([numpy.mean(rg.atoms().extract_b()) for rg in sel_h.residue_groups()])]
            # Plot 1
            axes[0].set_title('Individual contributions')
            axes[0].set_xlabel('Residue')
            axes[0].set_ylabel('Isotropic-ised B')
            for i in range(0,3):
                axes[0].plot(x_vals, tls_vals[i], l_styles[i], label=l_labels[i], markersize=2)
            #axes[0].legend()
            axes[0].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)
            # Plot 2
            axes[1].set_title('Cumulative contributions')
            axes[1].set_xlabel('Residue')
            axes[1].set_ylabel('Isotropic-ised B')
            for i in range(3,4):
                axes[1].plot(x_vals, tls_vals[i], l_styles[i], label=l_labels[i], markersize=2)
            axes[1].set_xticklabels([x_labels[int(i)] if (i<len(x_labels)) and (float(int(i))==i) else '' for i in axes[1].get_xticks()])
            pyplot.setp(axes[1].get_xticklabels(), rotation=90)
            #axes[1].legend()
            axes[1].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=1, mode="expand", borderaxespad=0.)
            # Plot boundaries
            v_lines = numpy.where(numpy.array([max(rg.atoms().extract_b()) for rg in b.select(sel).residue_groups()], dtype=bool))[0] + 1.5
            for val in v_lines:
                axes[0].axvline(x=val, ls='dotted')
                axes[1].axvline(x=val, ls='dotted')
            # Format and save
            pyplot.tight_layout()
            self.log('Writing: {}'.format(filename))
            pyplot.savefig(filename)
            pyplot.close(fig)

    def write_parameterisation_analysis(self, out_dir='.'):
        """Write atom-by-atom and dataset-by-dataset graphs"""

        easy_directory(out_dir)
        atm_dir = easy_directory(os.path.join(out_dir, 'atom_by_atom'))
        dst_dir = easy_directory(os.path.join(out_dir, 'dataset_by_dataset'))
        cor_dir = easy_directory(os.path.join(out_dir, 'error_correlations'))

        dst_labels = numpy.array([m.tag for m in self.models])

        # ----------------------------------------------------
        # Write the rmsds from the refined uijs
        # ----------------------------------------------------
        self.log.subheading('Calculating atom-by-atom rmsds to refined B-factors')
        h = self.blank_master_hierarchy()
        for i_g, (group, fit, sel) in self.enumerate_groups():
            # Select the atoms in this group, as list to allow for indexing
            sel_h = h.select(sel)
            sel_a = list(sel_h.atoms())
            # Extract the atom-by-atom rmsds fit/obs for each dataset
            rmsds = fit.uij_fit_obs_all_rmsds()
            # ------------------------
            # Write boxplots for all atoms in each dataset
            # ------------------------
            for i_d in range(0, len(self.models), 50):
                self.boxplot(filename=os.path.join(dst_dir, 'dataset-by-dataset-rmsds-datasets-group-{:02d}-{:04d}-{:04d}.png'.format(i_g+1,i_d+1,i_d+50)),
                             y_vals=rmsds.T[:,i_d:i_d+50],
                             x_labels=dst_labels[i_d:i_d+50].tolist(),
                             title='rmsd of fitted and refined B-factors (by dataset)',
                             x_lab='dataset', y_lab='rmsd', rotate_x_labels=True,
                             x_lim=(0,51), y_lim=(numpy.min(rmsds), numpy.max(rmsds)))
            # ------------------------
            # Write boxplots for all of the atoms in this group
            # ------------------------
            # Breaks between residues
            i_brk = numpy.array([float(i)+1.5 for i, (a1,a2) in enumerate(zip(sel_a,sel_a[1:])) if a1.parent().parent().resid()!=a2.parent().parent().resid()])
            for i_a in range(0, len(sel_a), 50):
                self.boxplot(filename=os.path.join(atm_dir, 'rmsds-group-{:02d}-atoms-{:06d}-{:06d}.png'.format(i_g+1,i_a+1,i_a+50)),
                             y_vals=rmsds[:,i_a:i_a+50],
                             x_labels=[ShortLabeller.format(a) for a in sel_a[i_a:i_a+50]],
                             title='rmsd of fitted and refined B-factors (by atom)',
                             x_lab='atom', y_lab='rmsd', rotate_x_labels=True,
                             x_lim=(0,51), y_lim=(numpy.min(rmsds), numpy.max(rmsds)), vlines=i_brk-i_a)
            # ------------------------
            # Write boxplots for atoms in each residue group separately
            # ------------------------
            for rg in sel_h.residue_groups():
                self.boxplot(filename=os.path.join(atm_dir, 'rmsds-residue-{}-group-{:02d}.png'.format(ShortLabeller.format(rg),i_g+1)),
                             y_vals=rmsds[:, [sel_a.index(a) for a in rg.atoms()]],
                             x_labels=[ShortLabeller.format(a) for a in rg.atoms()],
                             title='average rmsd of fitted and refined B-factors: {}'.format(ShortLabeller.format(rg)),
                             x_lab='atom', y_lab='rmsd', rotate_x_labels=True,
                             y_lim=(numpy.min(rmsds), numpy.max(rmsds)))
            # ------------------------
            # Append to overall array
            # ------------------------
            if i_g == 0: all_rmsds = rmsds
            else:        all_rmsds = numpy.append(all_rmsds, rmsds, axis=1)
        # ----------------------------------------------------
        # Write distribution of rmsds for each dataset
        # ----------------------------------------------------
        self.log.subheading('Calculating dataset-by-dataset rmsds (over all TLS groups)')
        for i_m in range(0, len(self.models), 50):
            m_idxs = numpy.arange(i_m,min(i_m+50,len(self.models)))
            self.boxplot(filename=os.path.join(dst_dir, 'dataset-by-dataset-rmsds-datasets-{:04d}-{:04d}.png'.format(i_m+1, i_m+50)),
                         y_vals=[all_rmsds[:,i].flatten() for i in m_idxs],
                         x_labels=dst_labels[m_idxs].tolist(),
                         title='rmsds for each dataset of fitted and refined B-factors',
                         x_lab='dataset', y_lab='rmsd', rotate_x_labels=True,
                         x_lim=(0,51), y_lim=(numpy.min(all_rmsds), numpy.max(all_rmsds)))

        # ----------------------------------------------------
        # Write averaged rmsds for all atoms
        # ----------------------------------------------------
        # Average values for all atoms + Interquartile range for all atoms
        h_avg = self.blank_master_hierarchy()
        h_iqr = self.blank_master_hierarchy()
        for i_g, (group, fit, sel) in self.enumerate_groups():
            # Extract rmsds for this group
            rmsds = fit.uij_fit_obs_all_rmsds()
            # Average over the datasets to obtain average per atom (also calculate the average fit/obs difference, should be ~0)
            h_avg.select(sel).atoms().set_b(flex.double(numpy.mean(rmsds, axis=0)))
            h_avg.select(sel).atoms().set_uij(flex.sym_mat3_double(fit.uij_fit_obs_atom_averaged_differences()))
            # Calculate IQRs of rmsds for each atom
            h_iqr.select(sel).atoms().set_b(flex.double(numpy.subtract(*numpy.percentile(rmsds, [75, 25],axis=0))))
        # Write out structures
        h_avg.write_pdb_file(os.path.join(out_dir, 'uij_fit_rmsds_averages.pdb'))
        h_iqr.write_pdb_file(os.path.join(out_dir, 'uij_fit_rmsds_iqranges.pdb'))
        # ----------------------------------------------------
        # Select the fitted atoms of this structure for some more graphs!
        # ----------------------------------------------------
        sel_h = h_avg.select(self.atom_selection_all)
        # ----------------------------------------------------
        # Write distribution of average rmsds for each residue residue
        # ----------------------------------------------------
        self.log.subheading('Calculating residue-by-residue rmsds to refined B-factors')
        all_bs = numpy.concatenate([list(rg.atoms().extract_b()) for rg in sel_h.residue_groups()])
        min_b, max_b = numpy.min(all_bs), numpy.max(all_bs)
        for sel_c in sel_h.chains():
            all_rgs = list(sel_c.residue_groups())
            for i_rg in range(0, len(all_rgs), 50):
                this_rgs = all_rgs[i_rg:i_rg+50]
                self.violinplot(filename=os.path.join(atm_dir, 'residue-by-residue-rmsds-chain-{}-residues-{:04d}-{:04d}.png'.format(sel_c.id,i_rg+1,i_rg+50)),
                                y_vals=[numpy.array(rg.atoms().extract_b()) for rg in this_rgs],
                                x_labels=[ShortLabeller.format(rg) for rg in this_rgs],
                                title='averaged rmsds of fitted and refined B-factors (for each residue)',
                                x_lab='residue', y_lab='rmsd', rotate_x_labels=True,
                                x_lim=(0,51), y_lim=(numpy.min(all_bs), numpy.max(all_bs)))

        # ----------------------------------------------------
        # Calculate correlations between the residual uij and the fitting error
        # ----------------------------------------------------
        self.log.subheading('Calculating correlations between residual uij and fitting errors')
        all_corr_table = pandas.DataFrame(index=[mdl.tag for mdl in self.models])
        all_atom_labels = []
        for i_g, (group, fit, sel) in self.enumerate_groups():
            uij_err = fit.uij_fit_obs_differences()
            uij_res = fit.fitted_uij_residual()
            # -------------------
            # Calculate the correlations between fitting error and residual uij
            # -------------------
            corr_table = pandas.DataFrame(index=[mdl.tag for mdl in self.models])
            atom_labels = [ShortLabeller.format(a) for a in self.master_h.atoms().select(sel)]
            all_atom_labels.extend(atom_labels)
            for i_atm, atm in enumerate(atom_labels):
                # Correlations between this atom uij_res and the fitting error across the datasets
                corr = numpy.corrcoef(uij_err[:,i_atm,:], uij_res[i_atm])[-1,:-1]
                assert corr.shape == (uij_err.shape[0],)
                if numpy.isnan(corr).all(): corr = numpy.random.randn(len(corr)) * 0.0000001
                corr_table[atm] = corr
            # Append this group's atoms to the main table
            all_corr_table = all_corr_table.join(corr_table, how="outer")
            # -------------------
            # Make violin plot of the correlations
            # -------------------
            sel_h = h.select(sel)
            sel_a = list(sel_h.atoms())
            i_brk = numpy.array([float(i)+1.5 for i, (a1,a2) in enumerate(zip(sel_a,sel_a[1:])) if a1.parent().parent().resid()!=a2.parent().parent().resid()])
            for i_a in range(0, len(corr_table.columns), 20):
                self.violinplot(filename=os.path.join(cor_dir, 'residual-correlation-group-{:02d}-atoms-{:06d}-{:06d}.png'.format(i_g+1,i_a+1,i_a+20)),
                                y_vals=corr_table.values[:,i_a:i_a+20].T.tolist(),
                                x_labels=[ShortLabeller.format(a) for a in sel_a[i_a:i_a+20]],
                                title='correlation of residuals (by atom)',
                                x_lab='atom', y_lab='correlation', rotate_x_labels=True,
                                x_lim=(0,21), y_lim=(-1,1), vlines=i_brk-i_a)
        # Write out csv
        filename = os.path.join(cor_dir, 'correlation_uij_residual_to_fitting_residual.csv')
        self.log('\nWriting: {}'.format(filename))
        all_corr_table.to_csv(filename)

    def write_tls_subtracted_models(self, out_dir):
        """Write models for each dataset with the TLS component subtracted from refined uij"""

        out_dir = easy_directory(out_dir)

        # ----------------------------------------------------
        # Write models for each dataset with uij(TLS) subtracted
        # ----------------------------------------------------
        for i_mdl, mdl in enumerate(self.models):
            self.log('Writing tls-subtracted uij model for {}'.format(mdl.tag))
            # Model with fitted TLS subtracted from input
            h_sub = mdl.hierarchy.deep_copy()
            for group, fit, sel in self.iterate_groups():
                # Model with fitted TLS subtracted from input
                tls = fit.extract_fitted_uij_tls(datasets=[i_mdl])[0]
                h_sub.atoms().select(sel).set_b(h_sub.atoms().select(sel).extract_b() - flex.double(EIGHT_PI_SQ*numpy.mean(tls[:,0:3],axis=1)))
                h_sub.atoms().select(sel).set_uij(h_sub.atoms().select(sel).extract_uij() - flex.sym_mat3_double(tls))
            # Residual (tls-subtracted)
            f_sub = os.path.join(out_dir, mdl.tag+'.tls-subtracted.pdb')
            self.log('\t> {}'.format(f_sub))
            h_sub.write_pdb_file(f_sub)

    def write_fitted_dataset_models(self, out_dir='./'):
        """Write fitted B-factors to output structures."""

        easy_directory(out_dir)

        # Extract the fitted output for each dataset
        for i_mdl, mdl in enumerate(self.models):
            self.log('Writing structure for model: {}'.format(mdl.filename))
            # Model with full set of fitted B-factors
            h_fit = mdl.hierarchy.deep_copy()
            for group, fit, sel in self.iterate_groups():
                # Model with full set of fitted B-factors
                uij = fit.extract_fitted_uij(datasets=[i_mdl])[0]
                h_fit.atoms().select(sel).set_b(flex.double(h_fit.atoms().select(sel).size(), 0))
                h_fit.atoms().select(sel).set_uij(flex.sym_mat3_double(uij))
            # Create fitted model paths
            mdl.o_pdb = os.path.join(out_dir, mdl.tag+'.mda.pdb')
            mdl.o_mtz = os.path.join(out_dir, mdl.tag+'.mda.mtz')
            if not os.path.exists(mdl.o_mtz):
                rel_symlink(mdl.i_mtz, mdl.o_mtz)
            # Write model
            self.log('\t> {}'.format(mdl.o_pdb))
            h_fit.write_pdb_file(mdl.o_pdb)

    def refine_fitted_dataset_models(self, suffix='-refined'):
        """Refine coordinates of the fitted structures"""

        self.log.heading('Refining coordinates of fitted models')

        refine_phenix.auto = False

        proc_args = []
        for mdl in self.models:
            if not os.path.exists(mdl.o_mtz): rel_symlink(mdl.i_mtz, mdl.o_mtz)
            obj = refine_phenix(pdb_file=mdl.o_pdb, mtz_file=mdl.o_mtz, cif_files=self.cifs,
                                out_prefix=os.path.splitext(mdl.o_pdb)[0]+suffix,
                                strategy='individual_sites+occupancies', n_cycles=1)
            obj.tag = mdl.tag
            proc_args.append(obj)

        # Refine all of the models
        refined = libtbx.easy_mp.pool_map(processes=self._n_cpu, func=wrapper_run, args=proc_args, chunksize=1)

        for mdl, ref in zip(self.models, proc_args):
            assert mdl.tag == ref.tag
            mdl.r_pdb = ref.out_pdb_file
            mdl.r_mtz = ref.out_mtz_file
            assert os.path.exists(mdl.r_pdb), '{} does not exist'.format(mdl.r_pdb)
            assert os.path.exists(mdl.r_mtz), '{} does not exist'.format(mdl.r_mtz)

    def generate_fitted_table_ones(self, out_dir='.'):
        """Write table-ones for each structure before and after fitting"""

        easy_directory(out_dir)

        self.log.subheading('Generating "Table Ones" for input and fitted B-factor models')

        for mdl in self.models:
            if not os.path.exists(mdl.o_mtz): rel_symlink(mdl.i_mtz, mdl.o_mtz)

            assert  os.path.exists(mdl.i_pdb) and \
                    os.path.exists(mdl.i_mtz) and \
                    os.path.exists(mdl.o_pdb) and \
                    os.path.exists(mdl.o_mtz)

        output_eff_orig = os.path.abspath(os.path.join(out_dir, 'table_one_input.eff'))
        output_eff_fitd = os.path.abspath(os.path.join(out_dir, 'table_one_fitted.eff'))
        output_eff_refd = os.path.abspath(os.path.join(out_dir, 'table_one_refined.eff'))

        # Save the names of the csvs (to be created)
        self.table_one_csv_input   = output_eff_orig.replace('.eff', '.csv')
        self.table_one_csv_fitted  = output_eff_fitd.replace('.eff', '.csv')
        self.table_one_csv_refined = output_eff_refd.replace('.eff', '.csv')

        phil = multi_table_ones.master_phil.extract()
        phil.input.dir        = []
        phil.input.labelling  = self.params.input.labelling
        phil.options          = self.params.table_ones_options
        phil.settings.cpus    = self.params.settings.cpus
        phil.settings.verbose = False

        # Run 1
        phil.input.pdb = [mdl.i_pdb for mdl in self.models]
        phil.output.parameter_file = output_eff_orig
        phil.output.output_basename = os.path.splitext(output_eff_orig)[0]
        multi_table_ones.run(params=phil)

        # Run 2
        phil.input.pdb = [mdl.o_pdb for mdl in self.models]
        phil.output.parameter_file = output_eff_fitd
        phil.output.output_basename = os.path.splitext(output_eff_fitd)[0]
        multi_table_ones.run(params=phil)

        # Run 3
        if self.models[0].r_pdb is not None:
            phil.input.pdb = [mdl.r_pdb for mdl in self.models]
            phil.output.parameter_file = output_eff_refd
            phil.output.output_basename = os.path.splitext(output_eff_refd)[0]
            multi_table_ones.run(params=phil)

        self.log.subheading('Running phenix.table_one to calculate R-factors')
        for f in [output_eff_orig,output_eff_fitd,output_eff_refd]:
            if not os.path.exists(f): continue
            cmd = CommandManager('phenix.table_one')
            cmd.add_command_line_arguments([f])
            cmd.print_settings()
            cmd.run()
            cmd.write_output(f.replace('.eff','.log'))

        # Clear all of the symlinks
        for mdl in self.models:
            if os.path.islink(mdl.o_mtz):
                os.remove(mdl.o_mtz)

        assert os.path.exists(self.table_one_csv_input)
        assert os.path.exists(self.table_one_csv_fitted)
        assert os.path.exists(self.table_one_csv_refined)

    def write_combined_csv(self, out_dir='.'):
        """Add data to CSV and write"""

        easy_directory(out_dir)
        self.log.subheading('Writing output csvs')

        # ----------------------------------------------------
        # Group-by-group rmsds CSV
        # ----------------------------------------------------
        for i, g in enumerate(self.groups):
            filename = os.path.join(out_dir, 'all_rmsd_scores_group_{}.csv'.format(i+1))
            g_table = pandas.DataFrame(data    = self.fits[g].uij_fit_obs_all_rmsds(),
                                       index   = [m.tag for m in self.models],
                                       columns = [ShortLabeller.format(a) for a in self.master_h.atoms().select(self.atom_selections[g])])
            self.log('Writing: {}'.format(filename))
            g_table.to_csv(filename)

        # ----------------------------------------------------
        # Main output CSV
        # ----------------------------------------------------
        # Extract dataset-by-dataset RMSDs
        self.table = pandas.DataFrame(index=[m.tag for m in self.models])
        all_dataset_rmsds = numpy.array([numpy.concatenate([self.fits[g].uij_fit_obs_all_rmsds()[i] for g in self.groups]) for i in range(len(self.models))])
        medn_rmsds = numpy.median(all_dataset_rmsds, axis=1)
        self.table['median_rmsds'] = medn_rmsds
        mean_rmsds = numpy.mean(all_dataset_rmsds, axis=1)
        self.table['mean_rmsds'] = mean_rmsds
        # Extract data from the table one CSVs
        self.log.subheading('Looking for table one data')
        for pref, csv in [('old-', self.table_one_csv_input),
                          ('new-', self.table_one_csv_fitted),
                          ('ref-', self.table_one_csv_refined)]:
            self.log('Reading: {}'.format(csv))
            if not os.path.exists(csv):
                raise Exception('Cannot read: {}'.format(csv))
            table_one = pandas.read_csv(csv, index_col=0, dtype=str).transpose()
            table_one['Low Res Limit'], table_one['High Res Limit'] = zip(*table_one['Resolution range'].apply(lambda x: x.split('(')[0].split('-')))
            table_one = table_one[['High Res Limit', 'Low Res Limit', 'Unique reflections','Completeness (%)','Wilson B-factor','R-work','R-free','Average B-factor']]
            for col in table_one:
                self.log('> Formatting col: {}'.format(col))
                table_one[col] = table_one[col].apply(lambda x: x.split('(')[0])
            # Redetermine data types
            table_one = table_one.apply(lambda x: pandas.to_numeric(x, errors='coerce'))
            table_one.columns = pref + table_one.columns
            # Transfer data to other
            self.table = self.table.join(table_one, how="outer")
        # Write output csv
        filename = os.path.join(out_dir, 'dataset_scores.csv')
        self.log('Writing output csv: {}'.format(filename))
        self.table.to_csv(filename)

        # Make graphs for the table
        self.write_r_factor_analysis(table=self.table, out_dir=os.path.join(out_dir,'graphs'))

    def write_r_factor_analysis(self, table, out_dir='./'):
        """Look at pre- and post-refinement graphs"""

        out_dir = easy_directory(out_dir)

        self.log.subheading('Writing final graphs')

        filename = os.path.join(out_dir, 'r_free_change.png')
        ax = table.plot(x='old-R-free', y='new-R-free', kind='scatter')
        pyplot.savefig(filename)
        pyplot.close(ax.get_figure())

        filename = os.path.join(out_dir, 'r_work_change.png')
        ax = table.plot(x='old-R-work', y='new-R-work', kind='scatter')
        pyplot.savefig(filename)
        pyplot.close(ax.get_figure())

        filename = os.path.join(out_dir, 'dataset_mean_rmsds.png')
        ax = table.plot(x='old-High Res Limit', y='mean_rmsds', kind='scatter')
        pyplot.savefig(filename)
        pyplot.close(ax.get_figure())

        return

    def boxplot(self, filename, y_vals, x_labels, title, x_lab='x', y_lab='y', x_lim=None, y_lim=None, rotate_x_labels=True, vlines=None):
        """Generate standard boxplot"""

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
        """Generate standard histogram"""

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

    def violinplot(self, filename, y_vals, x_labels, title, x_lab='x', y_lab='y', x_lim=None, y_lim=None, rotate_x_labels=True, vlines=None):
        """Generate standard violin plot"""

        self.log('Writing: {}'.format(filename))

        fig = pyplot.figure()
        pyplot.rc('font', family='monospace')
        pyplot.title(title)
        pyplot.violinplot(y_vals, showmeans=True)
        pyplot.xticks(range(1,len(x_labels)+1), x_labels)
        #pyplot.violinplot(y_vals, labels=x_labels, showmeans=True)
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

class MultiDatasetTLSFitter(object):

    _tls_weight = 1e6
    _amp_weight = 1e6
    _uij_weight = 1.0
    _ovr_weight = 1.0

    def __init__(self, observed_xyz, observed_uij, tls_params=None, n_tls=None, n_cpu=1, optimisation_datasets=None, log=None):

        if log is None: log = Log(verbose=True)
        self.log = log

        self._test = False
        self._iter = 0

        self._n_tls = n_tls
        self._n_cpu = n_cpu

        self.optimisation_rmsd = numpy.inf
        self.optimisation_penalty = numpy.inf

        # ---------------------------->
        # Input data
        # ---------------------------->
        self.observed_uij = numpy.array(observed_uij)
        self.observed_xyz = numpy.array(observed_xyz)
        # TODO Make this variable over the datasets ? TODO
        self.observed_com = numpy.mean(self.observed_xyz, axis=(0,1))
        self._box_size = (numpy.min(self.observed_xyz, axis=(0,1)), numpy.max(self.observed_xyz, axis=(0,1)))
        self._box_edge = numpy.array([(self._box_size[i][0],self._box_size[j][1],self._box_size[k][2]) for i,j,k in flex.nested_loop((2,2,2))])
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
            #vec_tls = [[0.0]*21]*(num_tls)
            vec_tls = [numpy.mean(self.observed_uij, axis=(0,1)).tolist()+[0.0]*15]+[[0.0]*21]*(num_tls-1)
            # Select the initial Uij to be the Uij with the smallest largest eigenvalue (smallest of the largest eigenvalues)
            #eig_max = numpy.apply_along_axis(func1d=lambda uij: flex.max(self.sym_mat3_eigenvalues(uij)), axis=2, arr=self.observed_uij)
            #uij_min = numpy.round(self.observed_uij[zip(*numpy.where(eig_max == numpy.min(eig_max)))[0]], 3)
            #vec_tls = [uij_min.tolist()+[0.0]*15]+[[0.0]*21]*(num_tls-1)
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
        if (optimisation_datasets is None) or len(optimisation_datasets)==0:
            self._mask_dsets = numpy.ones(self._n_dst, dtype=bool)
        else:
            self._mask_dsets = self._blank_dataset_selection()
            self._mask_dsets.put(optimisation_datasets, 1)
#        self._mask_dsets[:20] = True
        self._mask_atoms = numpy.ones(self._n_atm, dtype=bool)

        # ---------------------------->
        # Minimisation variables
        # ---------------------------->
        # Amount of expected variation for each data type (axes of simplex)
        self._del_tls_mdl = 0.25
        self._del_tls_amp = 0.1
        self._del_uij_res = 0.1
        self._update_del_simplex()
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
        self._update_after_optimisation()

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

    def _normalise_tls_amplitudes(self):
        for i_tls in range(self._n_tls):
            sel_amp = self._sel_tls_amp*self._sel_tls[i_tls]
            sel_mdl = self._sel_tls_mdl*self._sel_tls[i_tls]
            for sel_sub in [self._sel_t, self._sel_l, self._sel_s]:
                amp_mean = numpy.mean(self._var_current[sel_amp*sel_sub])
                if amp_mean < 1e-6: continue
                # Apply normalisation to amplitudes and TLS model
                self._var_current[sel_amp*sel_sub] *= (1.0/amp_mean)
                self._var_current[sel_mdl*sel_sub] *= (1.0*amp_mean)

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

    def _select(self, parameter_selection, datasets=None, atoms=None, max_datasets=None):
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
        # Limit the number of datasets
        if max_datasets is not None:
            self._cur_datasets = self._cur_datasets[:max_datasets]

    def _optimise(self, verbose=False):
        """Run the optimisation"""

        self._verbose=verbose
        self._n_call = 0

        # Initialise the RMSD measure
        self.optimisation_rmsd = 1e6
        self.optimisation_penalty = 0.0
        # Create simplex for these parameters
        cur_simplex = self._get_simplex(self._var_current_sel)
        # Optimise these parameters
        optimised = simplex.simplex_opt(dimension = len(cur_simplex[0]),
                                        matrix    = map(flex.double, cur_simplex),
                                        evaluator = self)
        # Extract and update current values
        self._adopt(optimised.get_solution())
        self._update_after_optimisation()

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
        # Extract parameters (model, amps, residuals)
        parameters = self._extract_parameters()
        # Calculate physical penalties - reject this set if model is not physical
        ppen = self._physical_penalties(parameters=parameters)
        # Print info line if necessary
        if self._verbose:
            if self._n_call%20==0:
                header = '[{}] -> ({:^10}, {:^10})'.format(', '.join(['{:>10}'.format('parameter') for r in sub_vector]), 'fit/rmsd', 'penalty')
                line = '-'*len(header)
                self.log(line)
                self.log(header)
                self.log(line)
            self._n_call += 1
        # Return now if physical penalty if non-zero to save time
        if ppen > 0.0:
            if self._verbose:
                self.log('[{}] -> ({:>10}, {:10.0f})'.format(', '.join(['{:+10.5f}'.format(r) for r in sub_vector]), 'UNPHYSICAL', ppen))
            return ppen
        # Get the fitted and the observed uijs
        uij_fit = self.extract_fitted_uij(datasets=self._cur_datasets, atoms=self._cur_atoms, parameters=parameters)
        uij_obs = self.extract_observed_uij(datasets=self._cur_datasets, atoms=self._cur_atoms)
        # Calculate RMSD
        rmsd = numpy.round(numpy.sqrt(numpy.mean(numpy.power(uij_obs-uij_fit, 2))),3)
        # Calculate fitting penalties (add to rmsd)
        fpen = self._fitting_penalties(uij_fit=uij_fit, uij_obs=uij_obs)
        # Update minima
        if rmsd+fpen < self.optimisation_rmsd+self.optimisation_penalty:
            self.optimisation_rmsd    = rmsd
            self.optimisation_penalty = fpen
        if self._verbose: print '[{}] -> ({:10f}, {:10.0f})'.format(', '.join(['{:+10.5f}'.format(r) for r in sub_vector]), rmsd, fpen)
        return rmsd+fpen

    ################################################################################################
    ###
    ###                             Constraints / Restraints / Search Parameters
    ###
    ################################################################################################

    def _update_del_simplex(self):
        self._del_simplex = numpy.array([self._del_tls_mdl]*self._n_prm_tls_mdl +
                                        [self._del_tls_amp]*self._n_prm_tls_amp +
                                        [self._del_uij_res]*self._n_prm_uij_res )

    def set_simplex_widths(self, tls_model=None, tls_amplitudes=None, uij_residual=None):
        if tls_model is not None:
            self._del_tls_mdl = tls_model
        if tls_amplitudes is not None:
            self._del_tls_amp = tls_amplitudes
        if uij_residual is not None:
            self._del_uij_res = uij_residual
        self._update_del_simplex()
        self.log.subheading('Setting/updating optimisation deltas')
        self.log('TLS Model Step Size:     {}'.format(self._del_tls_mdl))
        self.log('TLS Amplitude Step Size: {}'.format(self._del_tls_amp))
        self.log('Uij Residual Step Size:  {}'.format(self._del_uij_res))

    def _get_simplex(self, selection):
        starting_values = self._var_current[selection]
        starting_deltas = self._del_simplex[selection]
        starting_simplex = numpy.repeat([starting_values], len(starting_values)+1, axis=0)
        for i in range(len(starting_values)):
            starting_simplex[i+1][i] += starting_deltas[i]
        return starting_simplex

    def set_penalty_weights(self, tls_weight=None, amp_weight=None, uij_weight=None, ovr_weight=None):
        """Set penalties for parameters to be invalid"""
        if tls_weight is not None: self._tls_weight = tls_weight
        if amp_weight is not None: self._amp_weight = amp_weight
        if uij_weight is not None: self._uij_weight = uij_weight
        if ovr_weight is not None: self._ovr_weight = ovr_weight
        self.log.subheading('Setting/updating optimisation penalties')
        self.log('Invalid TLS Model Penalty:     {}'.format(self._tls_weight))
        self.log('Invalid Amplitude Penalty:     {}'.format(self._amp_weight))
        self.log('Invalid Uij Penalty:           {}'.format(self._uij_weight))
        self.log('Fitted > Observed Uij Penalty: {}'.format(self._ovr_weight))

    def _physical_penalties(self, parameters):
        tls_mdl, tls_amp, uij_res = parameters
        tls_penalties = [self._tls_penalty(values=v) for v in tls_mdl]
        amp_penalties = [self._amp_penalty(values=v) for v in tls_amp]
        uij_penalties = [self._uij_penalty(values=v) for v in uij_res]
        return numpy.sum(tls_penalties+amp_penalties+uij_penalties)

    def _fitting_penalties(self, uij_fit, uij_obs):
        fit_penalties = []; #[fit_penalties.extend([self._uij_penalty(values=vv) for vv in v]) for v in uij_fit]
        ovr_penalties = []; [ovr_penalties.extend([self._ovr_penalty(*vv) for vv in zip(*v)]) for v in zip(uij_fit,uij_obs)]
        return numpy.sum(fit_penalties+ovr_penalties)

    def _tls_penalty(self, values):
        assert len(values) == 21
        t,l,s = self._unpack_tls_parameters(vals=values)
        t_penalty = flex.max((self.sym_mat3_eigenvalues(t)<0.0).as_int())
        l_penalty = flex.max((self.sym_mat3_eigenvalues(l)<0.0).as_int())
        if numpy.sum(numpy.abs(s)) > 0.0:
            s_uij_values = self.uij_for_tls(xyz=self._box_edge, tls_vectors=numpy.array([[0.0]*12+list(s)]), origin=self.observed_com)
            s_penalty = numpy.max([flex.max((self.sym_mat3_eigenvalues(uij)<0.0).as_int()) for uij in s_uij_values])
        else:
            s_penalty = 0
        return self._tls_weight*numpy.sum([t_penalty, l_penalty, s_penalty])

    def _amp_penalty(self, values):
        return self._amp_weight*numpy.sum(values<0.0)

    def _uij_penalty(self, values):
        assert len(values) == 6
        eig_values = self.sym_mat3_eigenvalues(values)
        return self._uij_weight*flex.max((eig_values<0.0).as_int())

    def _ovr_penalty(self, fit, obs):
        """Add penalty for having fitted B-factors greater than observed"""
        eig_values_fit = self.sym_mat3_eigenvalues(fit)
        eig_values_obs = self.sym_mat3_eigenvalues(obs)
        return self._ovr_weight*(flex.max(eig_values_fit)>flex.max(eig_values_obs))

    ################################################################################################
    ###
    ###                             Update functions
    ###
    ################################################################################################

    def _update_after_optimisation(self):
        pass

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

    def fit(self, n_macro_cycles, n_micro_cycles):
        """Run macro-cycles of parameter optimisation"""

        # Cumulative selection for amplitude optimisation
        prev_opt = self._blank_parameter_selection()

        for i_macro in range(n_macro_cycles):

            #########################################
            self.log.heading('Macrocycle {}'.format(i_macro+1), spacer=True)
            #########################################
            if i_macro > 0:
                self.log.subheading('Updating parameters for next iteration')

                self.log('Removing atoms with high residual uij from TLS optimisation')
                self._update_atom_masks()
                self.log('Removing datasets with high fit rmsds from TLS optimisation')
                self._update_dataset_masks()

                self.log('Resetting tls amplitudes')
                self._reset_tls_amplitudes()
                self.log('Resetting uij residuals')
                self._reset_uij_residual()

                #########################################
                self.log.subheading('Optimising all TLS amplitudes')
                self.set_penalty_weights(ovr_weight=1.0)
                proc_args = []
                for i_dst in range(self._n_dst):
                    n = self.copy()
                    n._select(parameter_selection = self._sel_dst[i_dst],
                              datasets = [i_dst],
                              atoms    = None)
                    proc_args.append(n._prep_for_mp())
                self.log.subheading('Running optimisation')
                self._adopt_from_others(libtbx.easy_mp.pool_map(processes=self._n_cpu, func=wrapper_optimise, args=proc_args, chunksize=1))
                self.optimisation_summary(False)
                #########################################

            #########################################
            # Start-of-cycle summary
            #########################################
            self.settings_summary()

            #########################################
            for i_tls in range(self._n_tls):
                for c_name, c_sel in [('T',self._sel_t), ('L',self._sel_l), ('S',self._sel_s)]:
                    #########################################
                    # Add these parameters to the parameters for optimisation
                    this_opt = self._sel_tls[i_tls]*c_sel
                    if i_macro == 0:
                        prev_opt += this_opt
                    #########################################
                    # Run multiple cycles for each model
                    for i_micro in range(n_micro_cycles):
                        cycle_str = 'cycle: {}-{} of {}-{}'.format(i_macro+1, i_micro+1, n_macro_cycles, n_micro_cycles)
                        self.log.heading('Optimising {} parameters for TLS model {} ({})'.format(c_name, i_tls+1, cycle_str))
                        self._select(parameter_selection = self._sel_tls[i_tls]*self._sel_tls_mdl*c_sel,
                                     datasets = None,
                                     atoms    = None)
                        if (c_name=='T'):
                            self.set_simplex_widths(tls_model=0.1)
                        else:
                            self.set_simplex_widths(tls_model=1.0)
                        # Best fit possible (no penalties)
                        self.set_penalty_weights(ovr_weight=0.0)
                        self.log.subheading('Running optimisation')
                        self._optimise(verbose=True)
                        self.optimisation_summary(False)
                        #########################################
                        if (c_name != "T") or (i_macro>0):
                            self.log.heading('Optimising {} amplitudes for TLS model {} ({})'.format(c_name, i_tls+1, cycle_str))
                            self.set_penalty_weights(ovr_weight=0.01)
                            proc_args = []
                            for i_dst in range(self._n_dst):
                                n = self.copy()
                                n._select(parameter_selection = self._sel_dst[i_dst]*this_opt,
                                          datasets = [i_dst],
                                          atoms    = None)
                                proc_args.append(n._prep_for_mp())
                            self.log.subheading('Running optimisation')
                            self._adopt_from_others(libtbx.easy_mp.pool_map(processes=self._n_cpu, func=wrapper_optimise, args=proc_args, chunksize=1))
                            self.optimisation_summary(False)
                        #########################################
                        self.log.heading('Optimising TLS amplitudes for all datasets ({})'.format(cycle_str))
                        # Run one cycles to find optimum fit, then optimise
                        self.set_penalty_weights(ovr_weight=0.1)
                        proc_args = []
                        for i_dst in range(self._n_dst):
                            n = self.copy()
                            n._select(parameter_selection = self._sel_dst[i_dst]*prev_opt,
                                      datasets = [i_dst],
                                      atoms    = None)
                            proc_args.append(n._prep_for_mp())
                        self.log.subheading('Running optimisation')
                        self._adopt_from_others(libtbx.easy_mp.pool_map(processes=self._n_cpu, func=wrapper_optimise, args=proc_args, chunksize=1))
                        self.log.subheading('Normalising all TLS amplitudes to average of one')
                        self._normalise_tls_amplitudes()
                        self.optimisation_summary(False)
            #########################################
            self.log.heading('Optimising residual Uijs')
            self.set_penalty_weights(ovr_weight=0.0)
            for i, uij_del in enumerate([1.0,0.1]):
                self.log.subheading('Optimising residual Uijs (step {})'.format(i+1, 2))
                self.set_simplex_widths(uij_residual=uij_del)
                proc_args = []
                for i in range(self._n_atm):
                    n = self.copy()
                    n._select(parameter_selection = self._sel_atm[i],
                              datasets = None,
                              atoms    = [i])
                    proc_args.append(n._prep_for_mp())
                self.log.subheading('Running optimisation')
                self._adopt_from_others(libtbx.easy_mp.pool_map(processes=self._n_cpu, func=wrapper_optimise, args=proc_args, chunksize=1))
                self.optimisation_summary(False)
            #########################################
            self.log.subheading('End of macrocycle {}'.format(i_macro+1))
            self.optimisation_summary()

        return self

    def extract_observed_xyz(self, datasets=None, atoms=None):
        datasets, atoms = self._selection_filter(datasets=datasets, atoms=atoms)
        return self.observed_xyz[datasets][:,atoms]

    def extract_observed_uij(self, datasets=None, atoms=None):
        datasets, atoms = self._selection_filter(datasets=datasets, atoms=atoms)
        return self.observed_uij[datasets][:,atoms]

    def uij_for_tls(self, xyz, tls_vectors, origin):
        """Convert a set of parameter vectors to a set of uijs"""
        return numpy.sum([uij_from_tls_vector_and_origin(xyz=xyz, tls_vector=v, origin=origin) for v in tls_vectors], axis=0)

    def extract_fitted_uij_tls(self, datasets=None, atoms=None, parameters=None):
        """Extract tls fitted uijs for a subset of datasets or atoms"""

        # Get atoms and datasets if Nones
        datasets, atoms = self._selection_filter(datasets=datasets, atoms=atoms)
        # Extract the optimised values
        if parameters is None:
            tls_p, tls_a, uij_r = self._extract_parameters()
        else:
            tls_p, tls_a, uij_r = parameters

        # Extract coordinates for selected atoms/datasets
        xyz = self.extract_observed_xyz(datasets=datasets, atoms=atoms)
        assert xyz.shape == (len(datasets), len(atoms), 3)
        # Extract only those that we're interested in
        tls_a = tls_a[datasets]
        assert tls_a.shape == (len(datasets), self._n_tls, 3)
        # Multiply tls amplitudes and models
        tls_f = self._expand_tls_amplitudes(tls_amps=tls_a) * tls_p
        assert tls_f.shape == (len(datasets), self._n_tls, 21)
        assert len(xyz) == len(tls_f)
        # Calculate the tls component of uij
        uij_tls = numpy.array([self.uij_for_tls(xyz=xyz[i], tls_vectors=tls_f[i], origin=self.observed_com) for i in range(len(datasets))])
        assert uij_tls.shape == (len(datasets), len(atoms), 6)

        return uij_tls

    def extract_fitted_uij(self, datasets=None, atoms=None, parameters=None):
        """Extract total fitted uijs for a subset of datasets or atoms"""

        # Get atoms and datasets if Nones
        datasets, atoms = self._selection_filter(datasets=datasets, atoms=atoms)
        # Extract the optimised values
        if parameters is None:
            tls_p, tls_a, uij_r = self._extract_parameters()
        else:
            tls_p, tls_a, uij_r = parameters

        # Extract the TLS components
        uij_fit = self.extract_fitted_uij_tls(datasets=datasets, atoms=atoms, parameters=(tls_p, tls_a, uij_r))
        # Extract the residual
        uij_r = uij_r[atoms]
        assert uij_r.shape == (len(atoms), 6)
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

    def average_fitted_uij_tls(self, xyz, models=None, t=True, l=True, s=True):
        tls = self.fitted_tls_models()
        amp = numpy.mean(self.fitted_tls_amplitudes(), axis=0)
        tls = (self._expand_tls_amplitudes(tls_amps=numpy.array([amp]))[0]) * tls
        if models is not None: tls=tls[models]
        if t is not True: tls[:,00:06] = 0.0
        if l is not True: tls[:,06:12] = 0.0
        if s is not True: tls[:,12:21] = 0.0
        return self.uij_for_tls(xyz=xyz, tls_vectors=tls, origin=self.observed_com)

    def sym_mat3_eigenvalues(self, vals):
        assert len(vals) == 6
        return linalg.eigensystem_real_symmetric(vals).values()

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

    def settings_summary(self):
        """Print the number of atoms/datasets selected for modelling the TLS/Uij residual"""
        self.log.subheading('TLS model optimisation/Uij residual optimisation')
        self.log('\tUsing {}/{} atoms'.format(numpy.sum(self._mask_atoms), self._n_atm))
        self.log('\tUsing {}/{} datasets'.format(numpy.sum(self._mask_dsets), self._n_dst))

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
        self.log('Optimisation RMSD:    {}'.format(self.optimisation_rmsd))
        self.log('Optimisation Penalty: {}'.format(self.optimisation_penalty))
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

    easy_directory(params.output.out_dir)

    assert params.table_ones_options.column_labels
    assert params.table_ones_options.r_free_label

    log = Log(os.path.join(params.output.out_dir, '_fitting.log'), verbose=True)

    log('Building model list')
    if params.input.labelling == 'foldername':
        label_func = lambda f: os.path.basename(os.path.dirname(f))
    else:
        label_func = lambda f: os.path.basename(os.path.splitext(f)[0])
    models = [CrystallographicModel.from_file(f).label(tag=label_func(f)) for f in params.input.pdb]#[:10]

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

    p = MultiDatasetBFactorParameterisation(models = models,
                                            groups = params.input.tls_group,
                                            params = params,
                                            n_cpu  = params.settings.cpus,
                                            log = log)
    p.fit()
    p.write_output()

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



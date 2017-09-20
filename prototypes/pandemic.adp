#!/usr/bin/env ccp4-python

import os, sys, copy, traceback
import math, itertools

import scipy.stats, scipy.cluster
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

from giant.manager import Program
from giant.dataset import CrystallographicModel
from giant.structure.select import protein
from giant.structure.tls import uij_from_tls_vector_and_origin, extract_tls_from_pdb
from giant.structure.formatting import ShortLabeller, PhenixSelection
from giant.structure.select import backbone, sidechains
from giant.xray.crystal import CrystalSummary
from giant.xray.refine import refine_phenix
from giant.xray.tls import phenix_find_tls_groups

from giant.jiffies import multi_table_ones

try:
    import matplotlib
    matplotlib.interactive(False)
    from matplotlib import pyplot
    pyplot.switch_backend('agg')
    pyplot.style.use('ggplot')
    pyplot.interactive(0)
except Exception as e:
    print e

numpy.set_printoptions(threshold=numpy.nan)

from IPython import embed

EIGHT_PI_SQ = 8*math.pi*math.pi

############################################################################

PROGRAM = 'pandemic.adp'

DESCRIPTION = """
    Fit a consensus B-factor model to a series of datasets
"""

############################################################################

blank_arg_prepend = {'.pdb':'pdb=', '.cif':'cif='}

master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .help = "input pdb files - with isotropic/anisotropic b-factors"
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
}
fitting {
    auto_levels = *chain *auto_group *residue *backbone *sidechain
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
            .help = "remove atoms with high residual Uijs from this mask?"
            .type = bool
            .multiple = False
    }
    max_datasets_for_optimisation = 20
        .help = 'how many datasets should be used for optimising the TLS parameters?'
        .type = int
    max_resolution_for_optimisation = 3.0
        .help = 'resolution limit for dataset to be used for TLS optimisation'
        .type = float
    number_of_macro_cycles = 2
        .help = 'how many fitting cycles to run (over all levels)'
        .type = int
    number_of_micro_cycles = 2
        .help = 'how many fitting cycles to run (for each level)'
        .type = int
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
    verbose = True
        .type = bool
    dry_run = False
        .type = bool
}
""", process_includes=True)

############################################################################

def get_t_l_s_from_vector(vals):
    return vals[0:6], vals[6:12], vals[12:21]

def rms(vals, axis=None):
    return numpy.sqrt(numpy.mean(numpy.power(vals,2), axis=axis))

def uij_from_multiple_tls(xyz, tls_vectors, origin):
    """Convert a set of parameter vectors to a set of uijs"""
    return numpy.sum([uij_from_tls_vector_and_origin(xyz=xyz, tls_vector=v, origin=origin) for v in tls_vectors], axis=0)

def uij_to_b(uij):
    return EIGHT_PI_SQ*numpy.mean(uij[:,0:3],axis=1)

############################################################################

#def wrapper_optimise(arg):
#    arg._optimise(verbose=False)
#    return arg

def wrapper_run(arg):
    return arg.run()

def wrapper_fit(args):
    fitter, kw_args = args
    try:
        fitter.optimise(**kw_args)
        return fitter
    except:
        return traceback.format_exc()

############################################################################

class MultiDatasetUijParameterisation(Program):

    master_phil = master_phil

    def __init__(self, models, params, levels, level_labels=None, log=None):
        """Object for fitting a series of TLS models to a set of structures"""

        if log is None: log = Log(verbose=True)
        self.log = log

        self.params = params

        self.out_dir = params.output.out_dir

        self._n_cpu = params.settings.cpus
        self._n_opt = params.fitting.max_datasets_for_optimisation

        self._allow_isotropic = True

        self._opt_datasets_res_limit = params.fitting.max_resolution_for_optimisation
        self._opt_datasets_selection = []

        self.models = models
        self.levels = levels
        self.level_labels = level_labels if level_labels else range(1,len(levels)+1)
        self.fitter = None

        # Misc files
        self.cifs = None

        # Create plot object
        self.plot = MultiDatasetUijPlots(log=self.log)

        # Validate and add output paths, etc.
        self._init_input_models()
        self._init_level_groups()
        self._init_fitter()

        self.table = None
        self.table_one_csv_input   = None
        self.table_one_csv_fitted  = None
        self.table_one_csv_refined = None

        #self.write_running_parameters_to_log(params=params)

    def _init_input_models(self):
        """Prepare the input models"""

        self.log.subheading('Processing input models')

        # Use the first hierarchy as the reference
        self.log('Using {} as reference structure'.format(self.models[0].tag))
        self.master_h = self.models[0].hierarchy.deep_copy()

        errors = []

        for i_m, m in enumerate(self.models):
            # Check that all of the structures are the same
            if not self.master_h.is_similar_hierarchy(m.hierarchy):
                errors.append(Failure("Structures are not all the same. Model {}. File: {}".format(i_m, m.filename)))
                continue

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

        # Check for errors - Structures not all the same?
        if errors:
            for e in errors:
                print str(e)
            raise Failure("Not all structures are the same atomic model of the protein.")
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
            self.log('All atoms are missing anisotropic displacement parameters -- using the isotropic parameters instead')
            observed_uij = numpy.zeros_like(observed_uij)
            # Extract isotropic component from atoms (again only for the atoms we're interested in)
            observed_uij[:,:,0] = numpy.array([a.extract_b()/EIGHT_PI_SQ for a in atoms])[:,self.atom_mask]
            observed_uij[:,:,2] = observed_uij[:,:,1] = observed_uij[:,:,0]
        elif (observed_uij==-1).any():
            raise Failure('Some atoms for fitting (but not all) do not have anisotropically-refined B-factors -- either all atoms for fitting must have anistropic atoms or none')

        # Create the fitting object
        self.fitter = MultiDatasetHierarchicalUijFitter(observed_uij = observed_uij,
                                                        observed_xyz = observed_xyz,
                                                        level_array  = self.level_array,
                                                        level_labels = self.level_labels,
                                                        log = self.log)
        # Select the datasets for be used for TLS optimisation
        self.fitter.set_optimisation_datasets(self._opt_datasets_selection)

        # Write summary of the fitted model (groups & levels)
        model_dir = easy_directory(os.path.join(self.out_dir, 'model'))
        self.log.heading('Writing summary of the hierarchical model')
        self.hierarchy_summary(out_dir=model_dir)

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
            m_a.set_b(flex.double(iso))
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
        structure_dir = easy_directory(os.path.join(self.out_dir, 'structures'))
        # Fully parameterised structures
        self.log.bar()
        self.log('Writing fitted structures (full uij models)')
        self.log.bar()
        pdbs = self.output_structures(uij = uij_all,
                                      iso = map(uij_to_b, uij_all),
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
                               out_dir = structure_dir,
                               model_suffix = '.tls.pdb')
        # Level by level TLS-parameterised structures
        for i_level in range(len(self.levels)):
            self.log.bar()
            self.log('Writing fitted structures (TLS components only - level {})'.format(i_level+1))
            self.log.bar()
            self.output_structures(uij = uij_lvl[i_level],
                                   iso = map(uij_to_b, uij_lvl[i_level]),
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

        model_dir = easy_directory(os.path.join(self.out_dir, 'model'))

        # Write average uijs for each level
        self.log.heading('Summaries of the level-by-level TLS fittings')
        self.tls_level_summary(out_dir=model_dir)

        # Write residuals
        self.log.heading('Summaries of the residual atomic components')
        self.residuals_summary(out_dir=model_dir)

        #------------------------------------------------------------------------------#
        #---#      Plot distribution of parameterised uijs for all structures      #---#
        #------------------------------------------------------------------------------#

        uij_dir = easy_directory(os.path.join(self.out_dir, 'fit_uijs'))

        # Distributions of the uijs for groups
        self.log.heading('Calculating distributions of uijs over the model')
        self.fit_uij_distributions(uij_fit=uij_all,
                                   uij_inp=uij_inp,
                                   out_dir=uij_dir)

        #------------------------------------------------------------------------------#
        #---#        Compare the fittend and input uijs for all structures         #---#
        #------------------------------------------------------------------------------#

        fit_dir = easy_directory(os.path.join(self.out_dir, 'fit_metrics'))

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
        self.generate_fitted_table_ones(out_dir=os.path.join(self.out_dir, 'table_ones'))

        # Write output CSV of... everything
        self.log.heading('Writing output csvs')
        self.write_combined_csv(uij_fit=uij_all,
                                uij_inp=uij_inp,
                                out_dir=self.out_dir)

        #------------------------------------------------------------------------------#
        #                                                                              #
        #------------------------------------------------------------------------------#

    def output_structures(self,  uij, iso=None, out_dir='./', model_suffix='.pdb'):
        """Write sets of uij to models."""

        # Make sure ouput directory exists
        easy_directory(out_dir)

        # Validate AnoUijs
        uij = numpy.array(uij)
        assert uij.shape == (len(self.models), sum(self.atom_mask), 6)
        # Validation IsoBs
        if iso is not None:
            iso = numpy.array(iso)
            assert iso.shape == (len(self.models), sum(self.atom_mask))

        # Mask to allow us to apply uijs back to structure
        sel = flex.bool(self.atom_mask.tolist())
        # List of output pdbs
        pdbs = []

        # Apply to each model and write out
        for i, mdl in enumerate(self.models):
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
            # Create model paths and write model
            mdl_d = easy_directory(os.path.join(out_dir, mdl.tag))
            mdl_f = os.path.join(mdl_d, mdl.tag+model_suffix)
            self.log('{} > {}'.format(mdl.tag, mdl_f))
            h.write_pdb_file(mdl_f)
            pdbs.append(mdl_f)

        return pdbs

    def hierarchy_summary(self, out_dir='./'):
        """Write out the composition of the hierarchical model"""

        out_dir = easy_directory(out_dir)

        # Extract the global mask and convert to flex
        global_sel = flex.bool(self.atom_mask.tolist())

        # Write out the groups for each level
        for i_level, level in enumerate(self.fitter.levels):
            self.log('Writing partition groups for level {}'.format(i_level+1))
            for i_group, sel, fitter in level:
                sel = flex.bool(sel.tolist())
                g = self.blank_master_hierarchy()
                g.atoms().select(global_sel).select(sel).set_b(flex.double(sum(sel), 1))
                filename = os.path.join(out_dir, 'level-{:04d}_group-{:04d}-atoms.pdb'.format(i_level+1,i_group))
                self.log('\t> {}'.format(filename))
                g.write_pdb_file(filename)

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
            filename = os.path.join(out_dir, 'chain-{}-level-array.png'.format(c.id))
            self.log('\t> {}'.format(filename))
            self.plot.level_plots(filename=filename, hierarchies=hierarchies, title='chain {}'.format(c.id))

    def tls_level_summary(self, out_dir='./'):
        """Write the various TLS uijs to the master hierarchy structure"""

        out_dir = easy_directory(out_dir)

        self.log.subheading('Writing TLS models and amplitudes for each level')
        # Iterate through the levels
        for level in self.fitter.levels:
            self.log('Level {}'.format(level.index))
            # Table for TLS model components
            mdl_filename = os.path.join(out_dir, 'tls_models_level_{:04d}.csv'.format(level.index))
            mdl_table = pandas.DataFrame(columns=["group", "i_tls",
                                                  "T11","T22","T33","T12","T13","T23",
                                                  "L11","L22","L33","L12","L13","L23",
                                                  "S11","S12","S13","S21","S22","S23","S31","S32","S33"])
            # Iterate through the groups in this level
            for i_group, sel, fitter in level:
                tls_model, tls_amps = fitter.result()
                for i_tls in range(tls_model.shape[0]):
                    # Add model values to last row of table
                    mdl_table.loc[len(mdl_table.index)] = numpy.concatenate([[i_group, i_tls],tls_model[i_tls]])
                    # Create ampltidue table
                    amp_table = pandas.DataFrame(index=[mdl.tag for mdl in self.models],
                                                 columns=list("TLS"),
                                                 data=tls_amps[:,i_tls,:])
                    amp_filename = os.path.join(out_dir, 'tls_amplitudes_level_{:04d}_group_{:04d}_mode_{:04d}.csv'.format(level.index, i_group, i_tls+1))
                    amp_table.to_csv(amp_filename)
                    # Write histograms of amplitudes
                    x_vals = []; [[x_vals.append(tls_amps[:,i_t,i_a]) for i_a in range(3)] for i_t in range(tls_amps.shape[1])]
                    filename = os.path.join(out_dir, 'tls-model-amplitudes-level-{}-group-{}.png'.format(level.index, i_group))
                    self.log('\t> {}'.format(filename))
                    self.plot.histograms(filename = filename,
                                         x_vals   = x_vals,
                                         titles   = numpy.concatenate(['T (mode {a})-L (mode {a})-S (mode {a})'.format(a=i_t+1).split('-') for i_t in range(tls_amps.shape[1])]),
                                         x_labs   = ['']*numpy.product(tls_amps.shape[1:]), rotate_x_labels=True, shape=tls_amps.shape[1:], n_bins=30)
            # Write model table
            mdl_table.to_csv(mdl_filename)

        self.log.subheading('Writing T-L-S Uij components for each level')
        # Extract the atom mask to apply b-factors
        sel = flex.bool(self.atom_mask.tolist())
        # Write out separated average T-L-S-TLS components for the master hierarchy
        for tls_comp in ['T','L','S','TLS']:
            # Boolean selections
            t = 'T' in tls_comp
            l = 'L' in tls_comp
            s = 'S' in tls_comp
            self.log.bar(blank_before=True)
            self.log('{}{}{} components of fitted model'.format('T' if t else '_','L' if l else '_','S' if s else '_'))
            self.log.bar(blank_after=True)
            # Cumulative uijs (all)
            uij_all = None
            # Iterate through the levels
            for i_level, level in enumerate(self.fitter.levels):
                self.log('Level {}'.format(level.index))
                # Cumulative uijs (level)
                uij_lvl = None
                # Boundaries for this level
                boundaries = self.partition_boundaries_for_level(i_level=i_level)
                # Iterate through the different tls models - TODO FIXME allow n_tls > 0 FIXME TODO
                for i_tls in range(1):
                    # Create copy for resetting T-L-S components
                    l_copy = copy.deepcopy(level)
                    l_copy.reset_model(t = (not t),
                                       l = (not l),
                                       s = (not s))
                    # TODO FIXME allow n_tls > 0 FIXME TODO
                    #l_copy.reset_model(t = True,
                    #                   l = True,
                    #                   s = True,
                    #                   i_tls = [i for i in ...)
                    # Extract uijs
                    uij = l_copy.extract(average=True)
                    # Write out structure & graph of uijs (this level and this group)
                    prefix = os.path.join(out_dir, 'level_{}-mode_{}-{}'.format(level.index, i_tls+1, tls_comp))
                    m_h = self.custom_master_hierarchy(uij=uij, iso=uij_to_b(uij), mask=sel)
                    m_f = prefix+'.pdb'
                    self.log('\t> {}'.format(m_f))
                    m_h.write_pdb_file(m_f)
                    self.log('\t> {}*.png'.format(prefix))
                    self.plot.residue_by_residue(hierarchy=m_h,
                                                 prefix=prefix,
                                                 v_line_hierarchy=boundaries)

                    # Add to cumulative uij (this level)
                    uij_lvl = uij if uij_lvl is None else uij_lvl+uij
                # Write out structure & graph of uijs (this level)
                prefix = os.path.join(out_dir, 'level_{}-{}'.format(level.index, tls_comp))
                m_h = self.custom_master_hierarchy(uij=uij_lvl, iso=uij_to_b(uij_lvl), mask=sel)
                m_f = prefix+'.pdb'
                self.log('\t> {}'.format(m_f))
                m_h.write_pdb_file(m_f)
                #self.log('\t> {}*.png'.format(prefix))
                #self.plot.residue_by_residue(hierarchy=m_h,
                #                             prefix=prefix,
                #                             v_line_hierarchy=boundaries)

                # Add to cumulative uij (all levels)
                uij_all = uij_lvl if uij_all is None else uij_all+uij_lvl
            # Write out structure & graph of uijs (this level)
            prefix = os.path.join(out_dir, 'all-{}'.format(tls_comp))
            m_h = self.custom_master_hierarchy(uij=uij_all, iso=uij_to_b(uij_all), mask=sel)
            m_f = prefix+'.pdb'
            self.log('\t> {}'.format(m_f))
            m_h.write_pdb_file(m_f)
            self.log('\t> {}*.png'.format(prefix))
            self.plot.residue_by_residue(hierarchy=m_h,
                                         prefix=prefix)

    def residuals_summary(self, out_dir='./'):
        """Write the residual uijs to the master hierarchy structure"""

        out_dir = easy_directory(out_dir)

        uij = self.fitter.residual.extract()

        # Write out structure & graph of uijs (this level)
        self.log('Writing structure and plot of residual ADPs')
        prefix = os.path.join(out_dir, 'all-residual')
        m_h = self.custom_master_hierarchy(uij=uij, iso=uij_to_b(uij), mask=flex.bool(self.atom_mask.tolist()))
        m_f = prefix+'.pdb'
        self.log('\t> {}'.format(m_f))
        m_h.write_pdb_file(m_f)
        self.log('\t> {}*.png'.format(prefix))
        self.plot.residue_by_residue(hierarchy=m_h,
                                     prefix=prefix)

    def fit_uij_distributions(self, uij_fit, uij_inp, out_dir='./'):

        self.log('Not currently implemented')

        # Plot the distribution of B-factors for residues
        # Residues within groups
        # distribution of average uijs for members of each group for each leve
        #   > chain level > distribution of average TLS components for each dataset
        pass

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
                for i_d in range(0, len(self.models), max_x_width):
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
            for i_a in range(0, sum(sel), max_x_width):
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
            for i_r in range(0, len(g_res_sels), max_x_width):
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
        self.plot.residue_by_residue(hierarchy=m_h,
                                     prefix=prefix)

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
        self.plot.residue_by_residue(hierarchy=m_h,
                                     prefix=prefix)

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
        for i in range(res.shape[0]):
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
            for i_a in range(0, sum(sel), max_x_width):
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

        self.log.heading('Refining coordinates of fitted models')

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
        for mdl in self.models:
            if os.path.islink(mdl.o_mtz):
                os.remove(mdl.o_mtz)

        assert os.path.exists(self.table_one_csv_input)
        assert os.path.exists(self.table_one_csv_fitted)
        if self.models[0].r_pdb is not None:
            assert os.path.exists(self.table_one_csv_refined)

    def write_combined_csv(self, uij_fit, uij_inp, out_dir='./'):
        """Add data to CSV and write"""

        out_dir = easy_directory(out_dir)

        # Create output table
        self.table = pandas.DataFrame(index=[m.tag for m in self.models])

        self.log.bar()
        self.log('Calculating mean and median fitting rmsds by dataset')
        self.log.bar()
        # Calculate rmsd between input and fitted uijs
        uij_rmsd = rms(uij_inp-uij_fit, axis=2)
        # Extract mean/median dataset-by-dataset RMSDs
        dset_medn_rmsds = numpy.median(uij_rmsd, axis=1)
        dset_mean_rmsds = numpy.mean(uij_rmsd, axis=1)
        self.table['mean_rmsds']   = dset_mean_rmsds
        self.table['median_rmsds'] = dset_medn_rmsds
        for i in range(0, min(10, len(self.models))):
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
        self.table['mean_input_iso_b']  = dset_mean_inp_iso
        self.table['mean_fitted_iso_b'] = dset_mean_fit_iso
        for i in range(0, min(10, len(self.models))):
            self.log('Model {:10}: {:6.3f} (input) -> {:6.3f} (fitted)'.format(self.models[i].tag, dset_mean_inp_iso[i], dset_mean_fit_iso[i]))

        # Extract data from the table one CSVs
        self.log.subheading('Looking for table one data')
        for pref, csv in [('old-', self.table_one_csv_input),
                          ('new-', self.table_one_csv_fitted),
                          ('ref-', self.table_one_csv_refined)]:
            self.log('Reading: {}'.format(csv))
            if not os.path.exists(csv):
                if pref == 'ref-': continue
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
            self.log('')
        # Write output csv
        filename = os.path.join(out_dir, 'dataset_scores.csv')
        self.log('Writing output csv: {}'.format(filename))
        self.table.to_csv(filename)

        # Make graphs for the table
        self.write_r_factor_analysis(table=self.table, out_dir=os.path.join(out_dir,'graphs'))

    def write_r_factor_analysis(self, table, out_dir='./'):
        """Look at pre- and post-refinement graphs"""

        out_dir = easy_directory(out_dir)

        self.log.subheading('Model improvement summary')

        # Extract Column averages (means and medians)
        table_means = self.table.mean().round(3)
        out_str = '{:>30s} | {:>15} | {:>15} | {:>15}'

        self.log.bar()
        self.log('Dataset Averages:')
        self.log.bar()
        # Columns without old/new prefix
        for lab in table_means.index:
            if lab.startswith('new') or lab.startswith('old'):
                continue
            self.log(out_str.format(lab, table_means[lab], '', ''))
        self.log.bar()
        # Columns with old/new prefix
        self.log(out_str.format('', 'Single-dataset', 'Multi-dataset', 'Difference'))
        for new_lab in table_means.index:
            if not new_lab.startswith('new'):
                continue
            lab = new_lab[4:]
            old_lab = 'old-'+lab
            self.log(out_str.format(lab, table_means[old_lab], table_means[new_lab], table_means[new_lab]-table_means[old_lab]))
        self.log.bar()

#        filename = os.path.join(out_dir, 'r_free_change.png')
#        ax = table.plot(x='old-R-free', y='new-R-free', kind='scatter')
#        pyplot.savefig(filename)
#        pyplot.close(ax.get_figure())
#
#        filename = os.path.join(out_dir, 'r_work_change.png')
#        ax = table.plot(x='old-R-work', y='new-R-work', kind='scatter')
#        pyplot.savefig(filename)
#        pyplot.close(ax.get_figure())
#
#        filename = os.path.join(out_dir, 'dataset_mean_rmsds.png')
#        ax = table.plot(x='old-High Res Limit', y='mean_rmsds', kind='scatter')
#        pyplot.savefig(filename)
#        pyplot.close(ax.get_figure())

        return

class MultiDatasetUijPlots(object):

    def __init__(self, log=None):

        if log is None: log = Log(verbose=True)
        self.log = log

    def histograms(self, filename, x_vals, titles, x_labs, rotate_x_labels=True, shape=None, n_bins=30):
        """Generate standard histogram"""

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

    def boxplot(self, filename, y_vals, x_labels, title, x_lab='x', y_lab='y', x_lim=None, y_lim=None, rotate_x_labels=True, vlines=None):
        """Generate standard boxplot"""

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

    def violinplot(self, filename, y_vals, x_labels, title, x_lab='x', y_lab='y', x_lim=None, y_lim=None, rotate_x_labels=True, vlines=None):
        """Generate standard violin plot"""

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

    def level_plots(self, filename, hierarchies, title, rotate_x_labels=True):

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
        axis.set_yticklabels(['Level {}'.format(i) for i in range(1, len(hierarchies)+1)])
        pyplot.xlim((0, len(a_labels)+1))
        pyplot.setp(axis.get_xticklabels(), rotation=90)
        # Invert y-axis
        pyplot.gca().invert_yaxis()
        # Format and save
        pyplot.tight_layout()
        pyplot.savefig(filename, dpi=300)
        pyplot.close(fig)

    def residue_by_residue(self, hierarchy, prefix, title=None, v_line_hierarchy=None):
        """Write out residue-by-residue b-factor graphs"""

        h = hierarchy
        for chain_id in [c.id for c in h.chains()]:
            sel = h.atom_selection_cache().selection('chain {}'.format(chain_id))
            sel_h = h.select(sel)
            y_vals = numpy.array([numpy.mean(rg.atoms().extract_b()) for rg in sel_h.residue_groups()])
            if not y_vals.any(): continue # Skip chains with no Bs
            x_vals = numpy.array(range(len(list(sel_h.residue_groups()))))+1
            x_labels = ['']+[ShortLabeller.format(rg) for rg in sel_h.residue_groups()]
            filename = prefix + '-chain_{}.png'.format(chain_id)
            fig, axis = pyplot.subplots(nrows=1, ncols=1, sharey=True)
            if title is not None: axis.set_title(label=str(title))
            axis.set_xlabel('Residue')
            axis.set_ylabel('Isotropic B')
            axis.plot(x_vals, y_vals, '-ko', markersize=2)
            axis.set_xticklabels([x_labels[int(i)] if (i<len(x_labels)) and (float(int(i))==i) else '' for i in axis.get_xticks()])
            pyplot.setp(axis.get_xticklabels(), rotation=90)
            # Plot boundaries
            if v_line_hierarchy:
                v_lines = numpy.where(numpy.array([max(rg.atoms().extract_b()) for rg in v_line_hierarchy.select(sel).residue_groups()], dtype=bool))[0] + 1.5
                for val in v_lines:
                    axis.axvline(x=val, ls='dotted')
            # Format and save
            pyplot.tight_layout()
            pyplot.savefig(filename)
            pyplot.close(fig)

class MultiDatasetHierarchicalUijFitter(object):

    def __init__(self, observed_uij, observed_xyz, level_array, level_labels=None,  log=None):

        if log is None: log = Log(verbose=True)
        self.log = log

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

        # Level labels, and grouping for each level
        self.level_labels = level_labels if level_labels else ['Level {}'.format(i) for i in range(1, len(level_array)+1)]
        self.level_array = level_array

        # Series of objects to fit the Uij TLS groups
        self.levels = []
        for idx, (lab, group_idxs) in enumerate(zip(self.level_labels, self.level_array)):
            self.levels.append(MultiDatasetUijTLSGroupPartition(observed_uij = observed_uij,
                                                                observed_xyz = observed_xyz,
                                                                group_idxs   = group_idxs,
                                                                index = idx+1,
                                                                label = lab,
                                                                log=self.log))

        # One object to fit all the Uij residuals
        self.residual = MultiDatasetUijResidualFitter(observed_uij = observed_uij)

        assert len(self.level_labels) == len(self.level_array)
        assert len(self.level_labels) == len(self.levels)

        self.apply_masks()

        #self.summary(show=True)

    def _target_uij(self, fitted_uij_by_level, i_level):
        arr = fitted_uij_by_level.copy()
        arr[i_level] = 0.0
        arr_sum = numpy.sum(arr, axis=0)
        assert arr_sum.shape == self.observed_uij.shape
        return self.observed_uij - arr_sum

    def set_optimisation_datasets(self, dataset_indices):
        self.dataset_mask[:] = False
        self.dataset_mask[dataset_indices] = True

    def apply_masks(self):
        for level in self.levels+[self.residual]:
            level.set_dataset_mask(self.dataset_mask)
        for level in self.levels:
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

        self.log.heading('Fitting TLS models for {} levels (+ residual)'.format(len(self.levels)))

        # Update logfiles if more that one cpu
        #logfile = os.path.splitext(log.log_file)[0]+'level-{:04}.log'.format(idx)

        # Cumulative fitted Uij from the different levels (+1 level for the residuals!)
        # Indexing: [level, dataset, atom, uij]
        fitted_uij_by_level = numpy.zeros((len(self.levels)+1,)+self.observed_uij.shape)

        # Run macro cycles of optimisation
        for i_macro in range(n_macro_cycles):
            self.log.heading('Macrocycle {} of {}'.format(i_macro+1, n_macro_cycles), spacer=True)

            # Update masks at beginning of each cycles
            if i_macro > 0:
                self.log.subheading('Updating parameters for next iteration')
                self.log('Removing atoms with high residual uij from TLS optimisation')
                self.update_atomic_mask(percentile=90, atomic_uij=self.residual.extract())
                self.log('Removing datasets with high fit rmsds from TLS optimisation')
                self.update_dataset_mask(percentile=95, observed_fitted_differences=self.observed_uij-fitted_uij_by_level.sum(axis=0))
                self.log('Removing uij residuals from target function')
                fitted_uij_by_level[-1] = 0.0

            # Ensure the masks are up-to-date
            self.log('Setting atom and datasets masks')
            self.apply_masks()

            # Iterate through the TLS levels of the fitting
            for i_level, fitter in enumerate(self.levels):
                self.log.subheading('Fitting TLS Groups (level {} - {})'.format(fitter.index, fitter.label))
                # Update settings
                #self.log('Resetting TLS amplitudes')
                #fitter.reset_amplitudes()
                # Update the target uij by subtracting contributions from other levels
                self.log('Updating target Uijs for optimisation')
                fitter.set_target_uij(target_uij=self._target_uij(fitted_uij_by_level=fitted_uij_by_level, i_level=i_level))
                # Optimise
                self.log('Running optimisation')
                fitted_uij_by_level[i_level] = fitter.run(n_cpus=n_cpus, n_cycles=n_micro_cycles)

            # Fit the residuals
            self.log.subheading('Fitting residual atomic Uijs')
            # Update setttings
            self.log('Updating datasets masks')
            self.residual.set_dataset_mask(mask=self.dataset_mask)
            # Update the target uij by subtracting contributions from other levels
            self.residual.set_target_uij(target_uij=self._target_uij(fitted_uij_by_level=fitted_uij_by_level, i_level=-1))
            # Update fitters and optimise -- always run two cycles of this
            fitted_uij_by_level[-1] = self.residual.run(n_cpus=n_cpus, n_cycles=2)

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
        for level in self.levels+[self.residual]:
            s += '\n\n'
            s += level.summary(show=False).strip()
        if show: self.log(s)
        return s

class _MultiDatasetUijFitter(object):

    _uij_shape = None

    def run(self, n_cpus=1, n_cycles=1):
        jobs = [(fitter, {'n_cycles':n_cycles}) for (i, sel, fitter) in self]
        self.log.heading('Running {} jobs using {} cpus'.format(len(jobs), n_cpus))
        # Run jobs in parallel
        finished_jobs = libtbx.easy_mp.pool_map(processes=n_cpus, func=wrapper_fit, args=jobs, chunksize=1)
        # Record list of errors and raise all at end
        errors = []
        for i_iter, (i, sel, fitter) in enumerate(self):
            ret_job = finished_jobs[i_iter]
            if isinstance(ret_job, str) or (ret_job is None):
                errors.append((fitter,ret_job))
                continue
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
                self.log('Level "{}", Group "{}": Attempting to apply mask of less than one atom -- not applying mask'.format(self.label, fitter.label))

    def set_dataset_mask(self, mask):
        for i, sel, fitter in self:
            fitter.set_dataset_mask(mask)

class MultiDatasetUijTLSGroupPartition(_MultiDatasetUijFitter):

    def __init__(self, observed_uij, observed_xyz, group_idxs, index=0, label=None, log=None):

        if log is None: log = Log(verbose=True)
        self.log = log

        assert observed_uij.shape[:2] == observed_xyz.shape[:2]
        assert observed_uij.shape[2]  == 6
        assert observed_xyz.shape[2]  == 3
        assert observed_uij.shape[1]  == len(group_idxs)

        self.index = index
        self.label = label if label else 'Level {}'.format(index)

        self.group_idxs = group_idxs

        self._n_groups = sum(numpy.unique(self.group_idxs)!=-1)

        self._uij_shape = observed_uij.shape

        self.fitters = {}
        for i, sel, f in self:
            assert f is None
            self.fitters[i] = MultiDatasetUijTLSOptimiser(target_uij = observed_uij[:,sel],
                                                          atomic_xyz = observed_xyz[:,sel],
                                                          label = '{:4d} of {:4d}'.format(i, self._n_groups),
                                                          log = self.log)

    def __iter__(self):
        for i in numpy.unique(self.group_idxs):
            if i==-1: continue
            yield (i, (self.group_idxs==i), self.fitters.get(i, None))

    def extract(self, average=False):
        fitted_uij = numpy.zeros(self._uij_shape)
        for i, sel, fitter in self:
            fitted_uij[:,sel] = fitter.extract()
        if average:
            fitted_uij = fitted_uij.mean(axis=0)
        return fitted_uij

    def reset_amplitudes(self, t=False, l=False, s=False):
        for i, sel, fitter in self:
            fitter.reset_amplitudes(t=t,l=l,s=s)

    def reset_model(self, t=False, l=False, s=False):
        for i, sel, fitter in self:
            fitter.reset_model(t=t,l=l,s=s)

    def summary(self, show=True):
        s = self.log._subheading('TLS Partitions Summary (index {}; label {})'.format(self.index, self.label), blank=True)
        for i, sel, f in self:
            s += '\n'
            s += f.summary(show=False).strip()
        if show: self.log(s)
        return s

class MultiDatasetUijResidualFitter(_MultiDatasetUijFitter):

    def __init__(self, observed_uij, log=None):

        if log is None: log = Log(verbose=True)
        self.log = log

        self._uij_shape = observed_uij.shape

        self._n_atm = self._uij_shape[1]

        self.fitters = {}
        for i, sel, f in self:
            assert f is None
            self.fitters[i] = MultiDatasetUijAtomOptimiser(target_uij = observed_uij[:,sel],
                                                           label = 'atom {:5d} of {:5d}'.format(i,self._n_atm),
                                                           log = self.log)

    def __iter__(self):
        for i in range(self._n_atm):
            yield (i+1, i, self.fitters.get(i+1, None))

    def extract(self, expanded=False):
        fitted_uij = numpy.zeros(self._uij_shape[1:])
        for i, sel, fitter in self:
            fitted_uij[sel] = fitter.extract()
        if expanded:
            fitted_uij = fitted_uij.reshape((1,)+fitted_uij.shape).repeat(self._uij_shape[0], axis=0)
        return fitted_uij

    def summary(self, show=True):
        s = self.log._subheading('Uij Residuals Summary', blank=True)
        for i, sel, f in self:
            s += '\n'
            s += f.summary(show=False).strip()
        if show: self.log(s)
        return s

class _UijPenalties(object):
    _mdl_weight = 1e6
    _amp_weight = 1e6
    _uij_weight = 1.0
    _ovr_weight = 1.0

    def __init__(self, log=None):
        if log is None: log = Log(verbose=True)
        self.log = log

    def set_test_xyz(self, xyz, com):
        self._tst_xyz = xyz
        self._tst_com = com

    def set_weights(self, mdl_weight=None, amp_weight=None, uij_weight=None, ovr_weight=None, log=None):
        """Set penalties for parameters to be invalid"""
        if mdl_weight is not None: self._mdl_weight = mdl_weight
        if amp_weight is not None: self._amp_weight = amp_weight
        if uij_weight is not None: self._uij_weight = uij_weight
        if ovr_weight is not None: self._ovr_weight = ovr_weight
        return self.summary(show=False)

    def summary(self, show=False):
        s = 'Optimisation penalties'
        s += '\nInvalid TLS Model Penalty:     {}'.format(self._mdl_weight)
        s += '\nInvalid Amplitude Penalty:     {}'.format(self._amp_weight)
        s += '\nInvalid Uij Penalty:           {}'.format(self._uij_weight)
        s += '\nFitted > Observed Uij Penalty: {}'.format(self._ovr_weight)
        if show: self.log(s)
        return s

    def amplitudes(self, values):
        return self._amp_weight*numpy.sum(values<0.0)

    def tls_params(self, values):
        assert len(values) == 21
        t,l,s = get_t_l_s_from_vector(vals=values)
        t_penalty = flex.max((self._sym_mat3_eigenvalues(t)<0.0).as_int())
        l_penalty = flex.max((self._sym_mat3_eigenvalues(l)<0.0).as_int())
        if numpy.sum(numpy.abs(s)) > 0.0:
            s_uij_values = uij_from_multiple_tls(xyz=self._tst_xyz, tls_vectors=numpy.array([[0.0]*12+list(s)]), origin=self._tst_com)
            s_penalty = numpy.max([flex.max((self._sym_mat3_eigenvalues(uij)<0.0).as_int()) for uij in s_uij_values])
        else:
            s_penalty = 0
        return self._mdl_weight*numpy.sum([t_penalty, l_penalty, s_penalty])

    def uij_size(self, fitted, target):
        """Add penalty for having fitted B-factors greater than observed"""
        if self._ovr_weight == 0.0: return 0.0
        eig_values = self._sym_mat3_eigenvalues([a-b for a,b in zip(target,fitted)])
        return self._ovr_weight*sum(eig_values<0)
        #eig_values_fit = self._sym_mat3_eigenvalues(fitted)
        #eig_values_tar = self._sym_mat3_eigenvalues(target)
        #return self._ovr_weight*(flex.max(eig_values_fit)>flex.max(eig_values_tar))

    def uij_valid(self, values):
        assert len(values) == 6
        eig_values = self._sym_mat3_eigenvalues(values)
        return self._uij_weight*flex.max((eig_values<0.0).as_int())

    def _sym_mat3_eigenvalues(self, vals):
        assert len(vals) == 6
        return linalg.eigensystem_real_symmetric(vals).values()

class _Simplex(object):
    def get_simplex(self, start, selection=None):
        del_simplex = self._del_simplex
        if selection is not None:
            start = start[selection]
            del_simplex = del_simplex[selection]
        assert len(start) == len(del_simplex)
        start_simplex = numpy.repeat([start], len(start)+1, axis=0)
        for i in range(len(start)):
            start_simplex[i+1][i] += del_simplex[i]
        return start_simplex

class TLSSimplex(_Simplex):
    _del_mdl = 0.25
    _del_amp = 0.1

    def __init__(self, n_mdl_prm, n_amp_prm):
        """Initialise TLS Simplex"""
        self._n_mdl_prm = n_mdl_prm
        self._n_amp_prm = n_amp_prm
        self.set_deltas()

    def set_deltas(self, mdl=None, amp=None):
        if mdl is not None: self._del_mdl = mdl
        if amp is not None: self._del_amp = amp
        self._del_simplex = numpy.array((self._del_mdl,)*self._n_mdl_prm + (self._del_amp,)*self._n_amp_prm)

    def summary(self):
        s = 'Simplex optimisation deltas'
        s += '\nTLS Model Step Size:     {}'.format(self._del_mdl)
        s += '\nTLS Amplitude Step Size: {}'.format(self._del_amp)
        return s

class UijSimplex(_Simplex):
    _del_uij = 0.1

    def __init__(self):
        """Initialise AtomicUij Simplex"""
        self._n_prm = 6
        self.set_deltas()

    def set_deltas(self, uij=None):
        if uij is not None: self._del_uij = uij
        self._del_simplex = numpy.array((self._del_uij,)*self._n_prm)

    def summary(self):
        s = 'Simplex optimisation deltas'
        s += '\nUij Residual Step Size:  {}'.format(self._del_uij)
        return s

class _UijOptimiser(object):

    def __init__(self, target_uij, atomic_xyz=None, label='', log=None):

        if log is None: log = Log(verbose=True)
        self.log = log

        self._n_prm = 0
        self._n_dst = 0
        self._n_atm = 0

        self.target_uij = target_uij
        self.atomic_xyz = atomic_xyz

        self.label = label

        self._var_current = None
        self._var_current_sel = None

        self.optimisation_rmsd = numpy.inf
        self.optimisation_penalty = numpy.inf

        self.penalty = _UijPenalties()

    #===========================================+>
    # Private Functions
    #===========================================+>

    def _blank_atom_selection(self):
        return numpy.zeros(self._n_atm, dtype=bool)

    def _blank_dataset_selection(self):
        return numpy.zeros(self._n_dst, dtype=bool)

    def _blank_parameter_selection(self):
        return numpy.zeros(self._n_prm, dtype=bool)

    def _adopt(self, sub_vector, selection=None):
        """Insert a set of parameters into the complete parameter set"""
        if selection is None: selection = self._var_current_sel
        assert len(sub_vector) == numpy.sum(selection)
        self._var_current[selection] = sub_vector

    def _adopt_from_others(self, others):
        self._var_current_sel = self._blank_parameter_selection()
        for n in others:
            self._var_current_sel += n._var_current_sel
            self._adopt(sub_vector=n._var_current[n._var_current_sel], selection=n._var_current_sel)
        self._update_after_optimisation()

    def _optimise(self, verbose=False):
        """Run the optimisation"""

        # Prep variables for target function
        self._verbose = verbose
        self._n_call = 0
        # Apply the masks to the target uij
        self._update_target()
        # Initialise the RMSD measure
        self.optimisation_rmsd = 1e6
        self.optimisation_penalty = 0.0
        # Create simplex for these parameters
        opt_simplex = self.simplex.get_simplex(start     = self._var_current,
                                               selection = self._var_current_sel)
        # Optimise these parameters
        optimised = simplex.simplex_opt(dimension = len(opt_simplex[0]),
                                        matrix    = map(flex.double, opt_simplex),
                                        evaluator = self,
                                        tolerance = 1e-04)
        # Extract and update current values
        self._adopt(optimised.get_solution())
        self._update_after_optimisation()

    def _select(self, selection):
        """Select variables for optimisation"""
        self._var_current_sel = selection

    def _update_after_optimisation(self):
        pass

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

    def target(self, sub_vector):
        """Target function for the simplex optimisation"""

        # Increment counter
        self._n_call += 1
        # Combine the optimising parameters in the complete parameter set
        self._adopt(sub_vector)
        # Calculate physical penalties - reject this set if model is not physical
        ppen = self._parameter_penalties()
        # Print info line if necessary
        if self._verbose:
            if self._n_call%20==1:
                header = '[{}] -> ({:^10}, {:^10})'.format(', '.join(['{:>10}'.format('parameter') for r in sub_vector]), 'fit/rmsd', 'penalty')
                line = '-'*len(header)
                self.log(line)
                self.log(header)
                self.log(line)
        # Return now if physical penalty if non-zero to save time
        if ppen > 0.0:
            if self._verbose:
                self.log('[{}] -> ({:>10}, {:10.0f})'.format(', '.join(['{:+10.5f}'.format(r) for r in sub_vector]), 'UNPHYSICAL', ppen))
            return ppen
        # Get the fitted uijs (including masked atoms)
        self._update_fitted()
        # Calculate RMSD
        rmsd = numpy.sqrt(numpy.mean(numpy.power(self._target_uij-self._fitted_uij, 2)))
        # Calculate fitting penalties (add to rmsd)
        fpen = self._fitting_penalties(uij_fit=self._fitted_uij, uij_obs=self._target_uij)
        # Update minima
        if rmsd+fpen < self.optimisation_rmsd+self.optimisation_penalty:
            self.optimisation_rmsd    = rmsd
            self.optimisation_penalty = fpen
        if self._verbose:
            self.log('[{}] -> ({:10f}, {:10.0f})'.format(', '.join(['{:+10.5f}'.format(r) for r in sub_vector]), rmsd, fpen))
        return rmsd+fpen

class MultiDatasetUijAtomOptimiser(_UijOptimiser):

    def __init__(self, target_uij, label='', log=None):
        super(MultiDatasetUijAtomOptimiser, self).__init__(target_uij=target_uij, label=label, log=log)

        # Should be n_dataset observations of 6 parameters
        assert len(self.target_uij.shape) == 2
        assert self.target_uij.shape[1] == 6

        self._n_dst = self.target_uij.shape[0]

        # Number of parameters
        self._n_prm = 6

        # Initialise simplex generator
        self.simplex = UijSimplex()
        # Initialse penalty weights
        self.penalty.set_weights(uij_weight=1.0,
                                 ovr_weight=0.0)

        # Initialise loop variables for use during optimisation
        self._var_current = numpy.zeros(self._n_prm)
        self._select(self._blank_parameter_selection()+True)

        # Initialise the masks
        self.set_dataset_mask(range(self._n_dst))

    #===========================================+>
    # Private Functions - common to parent class
    #===========================================+>

    def _fitting_penalties(self, uij_fit, uij_obs):
        return 0.0

    def _parameter_penalties(self):
        return self.penalty.uij_valid(values=self.result())

    def _update_fitted(self):
        self._fitted_uij = self.extract()

    def _update_target(self):
        self._target_uij = self.target_uij[self._mask_dset]

    #===========================================+>
    # Public Functions - common to parent class
    #===========================================+>

    def extract(self):
        """Return the fitted uijs - for all atoms"""
        return tuple(self._var_current)

    def optimise(self, n_cycles=1):
        """Optimise the residual for a set of atoms"""
        uij_del_steps = itertools.cycle([0.1])
        for i_cycle in range(n_cycles):
            self.simplex.set_deltas(uij=uij_del_steps.next())
            self._optimise(verbose=False)
            self.summary(show=True)

    def result(self):
        """Return the fitted parameters (same as extract for this class)"""
        return tuple(self._var_current)

    def summary(self, show=True):
        """Print the number of parameters/input data"""
        uij = self.result()
        s = 'Uij ({}): '.format(self.label)+', '.join(['{:8.3f}'.format(v) for v in uij])
        if show: self.log(s)
        return s

class MultiDatasetUijTLSOptimiser(_UijOptimiser):

    def __init__(self, target_uij, atomic_xyz, n_tls=1, tls_params=None, label='', log=None):
        super(MultiDatasetUijTLSOptimiser, self).__init__(target_uij=target_uij, atomic_xyz=atomic_xyz, label=label, log=log)

        # Should be n_dataset observations of n_atm with 6 parameters
        assert len(self.target_uij.shape) == 3
        assert self.target_uij.shape[2] == 6
        assert self.atomic_xyz is not None

        self._n_dst = self.target_uij.shape[0]
        self._n_atm = self.target_uij.shape[1]

        assert self.target_uij.shape == (self._n_dst, self._n_atm, 6)
        assert self.atomic_xyz.shape == (self._n_dst, self._n_atm, 3)

        # Calculate the centre of mass of the atoms (for the rotation/screw components)
        self.atomic_com = numpy.mean(self.atomic_xyz, axis=1)

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
        self._n_prm = self._n_prm_mdl + self._n_prm_amp

        # Initialise parameter selection vectors
        self._sel_init()

        # Initialise simplex generator
        self.simplex = TLSSimplex(n_mdl_prm=self._n_prm_mdl, n_amp_prm=self._n_prm_amp)
        self.simplex.set_deltas(mdl=0.1, amp=0.1)
        # Initialse penalty set of test points (for identifying physically-valid TLS matrices)
        box_size = (numpy.min(self.atomic_xyz, axis=(0,1)),
                    numpy.max(self.atomic_xyz, axis=(0,1)))
        box_edge = numpy.array([(box_size[i][0],box_size[j][1],box_size[k][2]) for i,j,k in flex.nested_loop((2,2,2))])
        self.penalty.set_test_xyz(xyz=box_edge, com=self.atomic_com.mean(axis=0))

        # Initialise loop variables for use during optimisation
        self._var_current = numpy.zeros(self._n_prm)
        self._select(self._blank_parameter_selection()+True)
        self._adopt(sub_vector=numpy.concatenate(vec_tls),  selection=self._sel_mdl)
        self._adopt(sub_vector=numpy.ones(self._n_prm_amp), selection=self._sel_amp)

        # Initialise the masks
        self.set_atomic_mask(range(self._n_atm))
        self.set_dataset_mask(range(self._n_dst))

    #===========================================+>
    # Private Functions - common to parent class
    #===========================================+>

    def _fitting_penalties(self, uij_fit, uij_obs):
        fit_penalties = []; #[fit_penalties.extend([self._uij_penalty(values=vv) for vv in v]) for v in uij_fit]
        ovr_penalties = []; [ovr_penalties.extend([self.penalty.uij_size(*vv) for vv in zip(*v)]) for v in zip(uij_fit,uij_obs)]
        return numpy.sum(fit_penalties+ovr_penalties)

    def _parameter_penalties(self):
        tls_mdl, tls_amp = self.result()
        tls_penalties = [self.penalty.tls_params(values=v) for v in tls_mdl]
        amp_penalties = [self.penalty.amplitudes(values=v) for v in tls_amp]
        return numpy.sum(tls_penalties+amp_penalties)

    def _update_fitted(self):
        self._fitted_uij = self._extract(mask_datasets=self._mask_dset, mask_atoms=self._mask_atom)

    def _update_target(self):
        self._target_uij = self.target_uij[self._mask_dset][:,self._mask_atom]

    #===========================================+>
    # Private Functions - custom for this class
    #===========================================+>

    def _expand_tls_amplitudes(self, tls_amps):
        """Convert 3-element vector into 21 element vector for TLS multiplication"""
        n_dst = len(tls_amps)
        assert tls_amps.shape == (n_dst,self._n_tls,3)
        t_amps = numpy.repeat(tls_amps[:,:,0], 6, axis=1).reshape((n_dst, self._n_tls, 6))
        l_amps = numpy.repeat(tls_amps[:,:,1], 6, axis=1).reshape((n_dst, self._n_tls, 6))
        s_amps = numpy.repeat(tls_amps[:,:,2], 9, axis=1).reshape((n_dst, self._n_tls, 9))
        exp_tls_amps = numpy.concatenate([t_amps, l_amps, s_amps], axis=2)
        assert exp_tls_amps.shape == (n_dst,self._n_tls,21)
        return exp_tls_amps

    def _extract(self, mask_datasets=None, mask_atoms=None):
        """Calculate the TLS components for a selection of atoms and datasets"""
        # Extract parameters and coordinates
        mdls, amps = self.result()
        xyzs = self.atomic_xyz
        # Apply masks - datasets
        if mask_datasets is not None:
            assert isinstance(mask_datasets, list)
            n_dst = len(mask_datasets)
            xyzs = xyzs[mask_datasets]
            amps = amps[mask_datasets]
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
        assert mdls.shape == (self._n_tls, 21)
        assert amps.shape == (n_dst, self._n_tls, 3)
        assert xyzs.shape == (n_dst, n_atm, 3)
        # Expand the TLS amplitudes to 3s -> 21s
        tls = self._expand_tls_amplitudes(tls_amps=amps) * mdls
        assert tls.shape == (n_dst, self._n_tls, 21)
        # Calculate the TLS components
        uij = numpy.array([uij_from_multiple_tls(xyz=xyzs[i], tls_vectors=tls[i], origin=self.atomic_com[i]) for i in range(len(xyzs))])
        assert uij.shape == (n_dst, n_atm, 6)
        return uij

    def _normalise_tls_amplitudes(self):
        for i_tls in range(self._n_tls):
            # Extract selections for the model and the amplitudes
            sel_i_tls = self._sel_tls(i_tls)
            sel_mdl = self._sel_mdl*sel_i_tls
            sel_amp = self._sel_amp*sel_i_tls
            # Iterate through each T, L, S and normalise model by amplitudes
            for sel_sub in [self._sel_t, self._sel_l, self._sel_s]:
                amp_mean = numpy.mean(self._var_current[sel_amp*sel_sub])
                if amp_mean < 1e-6: continue
                # Apply normalisation to amplitudes and TLS model
                self._var_current[sel_amp*sel_sub] *= (1.0/amp_mean)
                self._var_current[sel_mdl*sel_sub] *= (1.0*amp_mean)

    def _sel_init(self):
        """Initialise a number of selection vectors for parameter groups"""

        # All TLS parameters
        self._sel_mdl = self._blank_parameter_selection()
        self._sel_mdl[0:self._n_prm_mdl] = 1.0
        assert sum(self._sel_mdl) == self._n_prm_mdl
        # All amplitudes
        self._sel_amp = self._blank_parameter_selection()
        self._sel_amp[-self._n_prm_amp:] = 1.0
        assert sum(self._sel_amp) == self._n_prm_amp
        # Any "T" component of a TLS parameter set
        self._sel_t = self._blank_parameter_selection()
        self._sel_t[self._sel_mdl] = ([1]*6 + [0]*6 + [0]*9)*self._n_tls
        self._sel_t[self._sel_amp] = [1,0,0]*self._n_tls*self._n_dst
        assert sum(self._sel_t) == 6*self._n_tls + self._n_tls*self._n_dst
        # Any "L" component of a TLS parameter set
        self._sel_l = self._blank_parameter_selection()
        self._sel_l[self._sel_mdl] = ([0]*6 + [1]*6 + [0]*9)*self._n_tls
        self._sel_l[self._sel_amp] = [0,1,0]*self._n_tls*self._n_dst
        assert sum(self._sel_l) == 6*self._n_tls + self._n_tls*self._n_dst
        # Any "S" component of a TLS parameter set
        self._sel_s = self._blank_parameter_selection()
        self._sel_s[self._sel_mdl] = ([0]*6 + [0]*6 + [1]*9)*self._n_tls
        self._sel_s[self._sel_amp] = [0,0,1]*self._n_tls*self._n_dst
        assert sum(self._sel_s) == 9*self._n_tls + self._n_tls*self._n_dst

        # Check selections do not overlap
        assert sum(self._sel_t*self._sel_l) == 0
        assert sum(self._sel_l*self._sel_s) == 0
        assert sum(self._sel_s*self._sel_t) == 0

    def _sel_dst(self, i):
        """Select all parameters for a particular dataset (ampltidues)"""
        sel = self._blank_parameter_selection()
        sel[self._sel_amp] = [0]*3*self._n_tls*i + [1]*3*self._n_tls + [0]*3*self._n_tls*(self._n_dst-i-1)
        return sel

    def _sel_tls(self, i):
        """Select all parameters for a particular tls model (tls params and amplitdues)"""
        sel = self._blank_parameter_selection()
        sel[self._sel_mdl] = [0]*21*i + [1]*21 + [0]*(self._n_tls-1-i)
        sel[self._sel_amp] = ([0]*3*i + [1]*3  + [0]*(self._n_tls-1-i))*self._n_dst
        return sel

    #===========================================+>
    # Public Functions
    #===========================================+>

    def extract(self):
        """Extract the fitted uij"""
        return self._extract(mask_datasets=None, mask_atoms=None)

    def optimise(self, n_cycles=1):
        """Optimise a (series of) TLS model(s) against the target data"""

        # Extract the masks (so that can be reapplied if changed)
        opt_dset_mask = self.get_dataset_mask()
        opt_atom_mask = self.get_atomic_mask()

        # Cumulative parameter selection
        s_cuml = self._blank_parameter_selection()

        for i_cycle in range(n_cycles):
            self.log.heading('Group {} - Optimisation cycle {} of {}'.format(self.label, i_cycle+1, n_cycles))
            for i_tls in range(self._n_tls):
                self.log.subheading('Optimising TLS model {} of {}'.format(i_tls+1, self._n_tls))
                self.log('Optimising using {} atoms'.format(len(opt_atom_mask)))
                self.log('Optimising using {} datasets'.format(len(opt_dset_mask)))
                # Extract selection for this TLS model's parameters
                s_tls = self._sel_tls(i_tls)
                # Optimise each T, L, S component separately
                for n_cmp, s_cmp in [('T',self._sel_t), ('L',self._sel_l), ('S',self._sel_s)]:
                    # Combine selections and add to cumulative
                    s_this = s_tls*s_cmp
                    s_cuml += s_this
                    # Optimise TLS model parameters
                    self.log.subheading('Optimising {} parameters for TLS model {}'.format(n_cmp, i_tls+1))
                    if n_cmp == 'T':
                        self.penalty.set_weights(ovr_weight=0.1)
                    else:
                        self.penalty.set_weights(ovr_weight=0.0)
                    self._select(self._sel_mdl*s_this)
                    self._optimise(verbose=False)
                    self.model_summary(show=True)
                    # Optimise TLS amplitude parameters (all ampltidues!)
                    self.log.subheading('Optimising TLS amplitudes for all datasets')
                    self.penalty.set_weights(ovr_weight=1.0)
                    for i_dst in range(self._n_dst):
                        self.set_dataset_mask([i_dst])
                        self._select(selection=self._sel_dst(i_dst)*s_cuml)
                        self._optimise(verbose=False)
                        self.log('> dataset {} of {} (rmsd {})'.format(i_dst+1,self._n_dst,self.optimisation_rmsd))
                    self.log.bar(blank_before=True, blank_after=True)
                    self.amplitude_summary(show=True)
                    # Reapply dataset mask
                    self.set_dataset_mask(opt_dset_mask)
            # End of cycle house-keeping and summary
            self._normalise_tls_amplitudes()
            self.log.subheading('End-of-cycle summary')
            self.summary(show=True)
        self.log.subheading('Optimisation complete')

    def reset_amplitudes(self, t=False, l=False, s=False, i_tls=None):
        selection = self._blank_parameter_selection()
        if t is True: selection += self._sel_t
        if l is True: selection += self._sel_l
        if s is True: selection += self._sel_s
        if i_tls is not None:
            tls_sel = self._blank_parameter_selection()
            i_tls = [i_tls] if isinstance(i_tls,int) else i_tls
            for i in i_tls:
                tls_sel += self._sel_tls(i_tls)
            selection *= tls_sel
        selection *= self._sel_amp
        self._var_current[selection] = 1.0

    def reset_model(self, t=False, l=False, s=False, i_tls=None):
        selection = self._blank_parameter_selection()
        if t is True: selection += self._sel_t
        if l is True: selection += self._sel_l
        if s is True: selection += self._sel_s
        if i_tls is not None:
            tls_sel = self._blank_parameter_selection()
            i_tls = [i_tls] if isinstance(i_tls,int) else i_tls
            for i in i_tls:
                tls_sel += self._sel_tls(i_tls)
            selection *= tls_sel
        selection *= self._sel_mdl
        self._var_current[selection] = 0.0

    def result(self):
        """Extract the fitted parameters"""
        tls_mdls = self._var_current[self._sel_mdl].reshape((self._n_tls, 21))
        tls_amps = self._var_current[self._sel_amp].reshape((self._n_dst, self._n_tls, 3) )
        return (tls_mdls,tls_amps)

    def summary(self, show=True):
        """Print the number of parameters/input data"""
        s = self.log._bar()+'\nTLS Group Fit Summary: {}\n'.format(self.label)+self.log._bar()
        s += '\n> Input summary'
        s += '\nNumber of datasets:   {}'.format(self._n_dst)
        s += '\nNumber of atoms:      {}'.format(self._n_atm)
        s += '\ninput uij parameters: {}'.format(self.target_uij.shape)
        s += '\ninput xyz parameters: {}'.format(self.atomic_xyz.shape)
        s += '\nCentre of mass: {}'.format(tuple(self.atomic_com.mean(axis=0).round(2).tolist()))
        s += '\n> Parameterisation summary'
        s += '\nNumber of TLS models: {}'.format(self._n_tls)
        s += '\nNumber of parameters for TLS fitting: {}'.format(self._n_prm)
        s += '\nNumber of observations (all): {}'.format(numpy.product(self.target_uij.shape))
        s += '\nData/parameter ratio (all) is {:.3f}'.format(numpy.product(self.target_uij.shape)*1.0/self._n_prm)
        if hasattr(self,'_target_uij'):
            n_obs_used = numpy.product(self._target_uij.shape)
            s += '\nNumber of observations (used): {}'.format(n_obs_used)
            s += '\nData/parameter ratio (used) is {:.3f}'.format(n_obs_used*1.0/self._n_prm)
        s += '\n> Atoms/Datasets for TLS model optimisation'
        s += '\n\tUsing {}/{} atoms'.format(len(self.get_atomic_mask()), self._n_atm)
        s += '\n\tUsing {}/{} datasets'.format(len(self.get_dataset_mask()), self._n_dst)
        if self.optimisation_rmsd is not numpy.inf:
            s += '\n> Optimisation Summary'
            s += '\nOptimisation RMSD:    {}'.format(self.optimisation_rmsd)
            s += '\nOptimisation Penalty: {}'.format(self.optimisation_penalty)
        s += '\n' + self.model_summary(show=False)
        s += '\n' + self.amplitude_summary(show=False)
        if show: self.log(s)
        return s

    def model_summary(self, show=False):
        """Print summary of the TLS components"""
        tls_prms, tls_amps = self.result()
        s = '> TLS parameters'
        for i, vals in enumerate(tls_prms):
            tls_vars = get_t_l_s_from_vector(vals=vals)
            s += '\nMode {:2}'.format(i+1)
            s += '\n\tT: '+', '.join(['{:8.3f}'.format(v) for v in tls_vars[0]])
            s += '\n\tL: '+', '.join(['{:8.3f}'.format(v) for v in tls_vars[1]])
            s += '\n\tS: '+', '.join(['{:8.3f}'.format(v) for v in tls_vars[2]])
        if show: self.log(s)
        return s

    def amplitude_summary(self, show=False):
        """Print summary of the TLS amplitudes"""
        tls_prms, tls_amps = self.result()
        s = '> TLS amplitudes'
        for i, vals in enumerate(tls_amps):
            s += '\nDataset {:4}'.format(i+1)
            for i, a in enumerate(vals):
                s += '\n\tMode {:2}:'.format(i+1)+' {:8.3f} (T) {:8.3f} (L) {:8.3f} (S)'.format(*a)
        if show: self.log(s)
        return s

############################################################################

def run(params):

    easy_directory(params.output.out_dir)

    assert params.table_ones_options.column_labels
    assert params.table_ones_options.r_free_label

    log = Log(os.path.join(params.output.out_dir, '_fitting.log'), verbose=params.settings.verbose)

    # Report parameters
    log.heading('Processed parameters')
    log(master_phil.format(params).as_str())

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

    # Build the levels for the hierarchical fitting
    log.subheading('Building hierarchy selections')
    levels = []; labels=[];
    # TODO IF statement for if custom levels are defined TODO
    # TODO allow ability to insert
    # FIXME Only run on the protein for the moment FIXME
    filter_h = protein(models[0].hierarchy)
    if 'chain' in params.fitting.auto_levels:
        log('Level {}: Creating level with groups for each chain'.format(len(levels)+1))
        levels.append([PhenixSelection.format(c) for c in filter_h.chains()])
        labels.append('chain')
    if 'auto_group' in params.fitting.auto_levels:
        log('Level {}: Creating level with groups determined by phenix.find_tls_groups'.format(len(levels)+1))
        levels.append([s.strip('"') for s in phenix_find_tls_groups(models[0].filename)])
        labels.append('groups')
    if 'residue' in params.fitting.auto_levels:
        log('Level {}: Creating level with groups for each residue'.format(len(levels)+1))
        levels.append([PhenixSelection.format(r) for r in filter_h.residue_groups()])
        labels.append('residue')
    if 'backbone' in params.fitting.auto_levels:
        log('Level {}: Creating level with groups for each residue backbone'.format(len(levels)+1))
        levels.append([PhenixSelection.format(r)+' and (name C or name CA or name N or name O)'     for r in backbone(filter_h).atom_groups() if (r.resname not in ['ALA','GLY','PRO'])])
        labels.append('backbone')
    if 'sidechain' in params.fitting.auto_levels:
        log('Level {}: Creating level with groups for each residue sidechain'.format(len(levels)+1))
        levels.append([PhenixSelection.format(r)+' and not (name C or name CA or name N or name O)' for r in sidechains(filter_h).atom_groups() if (r.resname not in ['ALA','GLY','PRO'])])
        labels.append('sidechain')
    # TODO Take any custom groups and insert them here TODO (levels.insert)
    # Report
    log('\n> {} levels created\n'.format(len(levels)))
    for i_l, level in enumerate(levels):
        log.bar()
        log('Level {}'.format(i_l+1))
        log.bar()
        for l in level: log('\t'+l)

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

    p.fit_hierarchical_uij_model()
    p.process_results()
    log.heading('Parameterisation complete')

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



import os, sys, glob, copy, shutil, traceback
import math, itertools, operator

from IPython import embed
import time, tqdm, signal

import scipy.cluster
import numpy, pandas

import libtbx.phil, libtbx.easy_mp, libtbx.easy_pickle
import iotbx.pdb
import mmtbx.tls.tools

from libtbx import adopt_init_args
from libtbx.utils import Sorry, Failure
from scitbx.array_family import flex

from bamboo.common import Meta, ListStream, dict_from_class
from bamboo.common.logs import Log, LogStream, ScreenLogger
from bamboo.common.path import easy_directory, rel_symlink, foldername, filename
from bamboo.common.command import CommandManager
from bamboo.maths.functions import rms, Sigmoid

from giant.manager import Program
from giant.dataset import CrystallographicModel, AtomicModel
from giant.structure.uij import uij_positive_are_semi_definite, \
        sym_mat3_eigenvalues, uij_to_b, calculate_uij_anisotropy_ratio, \
        scale_uij_to_target_by_selection
from giant.structure.formatting import ShortLabeller, PhenixSelection, PymolSelection
from giant.structure.select import protein, backbone, sidechains, default_secondary_structure_selections_filled
from giant.structure.pymol import auto_residue_images, auto_chain_images, selection_images
from giant.xray.crystal import CrystalSummary
from giant.xray.refine import refine_phenix, refine_refmac
from giant.xray.tls import phenix_find_tls_groups

from giant.jiffies import multi_table_ones

from pandemic.adp import html as pandemic_html
from pandemic.adp.optimise import MultiDatasetUijFitter, MultiDatasetTLSFitter, InterLevelAmplitudeOptimiser
from mmtbx.tls.utils import uij_eigenvalues, TLSMatrices, TLSAmplitudes

try:
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.interactive(False)
    from matplotlib import pyplot, patches
    pyplot.switch_backend('agg') # yes I know this done twice -- for safety!
    pyplot.style.use('ggplot')
    pyplot.interactive(0)
    pyplot.rc('font', family='monospace')
except Exception as e:
    print e

numpy.set_printoptions(linewidth=numpy.inf, threshold=numpy.nan)

EIGHT_PI_SQ = 8*math.pi*math.pi

tls_amplitudes_hash = {
        'simple_multiplier': TLSAmplitudes,
        }

############################################################################

import multiprocessing

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
    Fit a hierarchical Atomic Displacement Parameter (B-factor) model to a (series of) molecular structure(s)
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
    input_uij_model = isotropic anisotropic *mixed
        .type = choice(multi=False)
        .multiple = False
        .help = "Control whether an error is raised if selected atoms in structures have anisotropic/isotropic/mixed Bs."
    model_type = *crystallographic noncrystallographic
        .help = "is this a crystal structure? or a cryo-EM model? or any other non-crystallographic model?"
        .type = choice(multi=False)
    look_for_reflection_data = True
        .type = bool
    reference_r_values = None
        .type = path
        .help = "Reference R-values against which to compare the fitted and refined models. Labelling the CSV must be the same as the dataset labelling. If not provided will use the R-values of the input structures."
    pickle = None
        .type = path
        .multiple = False
}
output {
    out_dir = pandemic-adp
        .help = "output directory"
        .type = str
    pickle = False
        .type = bool
        .multiple = False
    html = True
        .help = "Write all output data in an html file"
        .type = bool
    images {
        all = False
            .help = "Generate all graphical output - will significantly increase runtime"
            .type = bool
        pymol = none *chain all
            .help = "Write residue-by-residue images of the output B-factors"
            .type = choice(multi=False)
        distributions = False
            .help = "Write distribution graphs for each TLS group"
            .type = bool
    }
    copy_reflection_data_to_output_folder = True
        .type = bool
    clean_up_files = *compress_logs *delete_mtzs
        .help = "Delete unnecessary output files (MTZs)"
        .type = choice(multi=True)
}
levels {
    overall_selection = "(not water) and (not element H)"
        .type = str
    auto_levels = *chain auto_group ss *secondary_structure *residue *backbone *sidechain atom
        .type = choice(multi=True)
    cbeta_in_backbone = True
        .help = "Flag to control whether the c-beta atom is considered part of the backbone or the sidechain"
        .type = bool
    custom_level
        .multiple = True
    {
        depth = None
            .help = "Where to insert this level into the hierarchy? (after auto_levels have been determined). Inserting as '1' will make this the first level and '2' will add as the second level, etc. Any auto-generated levels will be shifted down appropriately."
            .type = int
            .multiple = False
        label = None
            .help = "What do I call the level?"
            .type = str
            .multiple = False
        selection = None
            .help = "list of selections that define groups for this level"
            .type = str
            .multiple = True
    }
}
fitting {
    n_tls_modes_per_tls_group = 1
        .help = 'how many TLS models to fit per group of atoms?'
        .type = int
    tls_amplitude_model = *simple_multiplier
        .help = 'how do amplitudes transform each TLS group in each dataset? \n\tsimple_multiplier = T->aT, L->aL, S->aS.'
        .type = choice(multi=False)
    number_of_macro_cycles = 20
        .help = 'how many fitting cycles to run (over all levels) -- should be more than 0'
        .type = int
    number_of_micro_cycles = 3
        .help = 'how many fitting cycles to run (for each level) -- should be more than 1'
        .type = int
    fit_residual = True
        .help = 'should the residual level be optimised? (it may be preferable to turn this off for small numbers of low-resolution structures)'
        .type = bool
    optimisation {
        dataset_weights = one inverse_resolution inverse_resolution_squared *inverse_resolution_cubed
            .help = 'control how datasets are weighted during optimisation?'
            .type = choice(multi=False)
        atom_weights = one *inverse_mod_U inverse_mod_U_squared inverse_mod_U_cubed
            .help = 'control how atoms are weighted during optimisation?'
            .type = choice(multi=False)
        renormalise_atom_weights_by_dataset = True
            .help = "Should atom weights be normalised to an average of 1 for each dataset?"
            .type = bool
        sort_datasets_by = *resolution random name
            .help = 'how should datasets be ordered when selecting a subset of datasets for optimisation?'
            .type = choice(multi=False)
        random_seed = 0
            .help = 'random number seed for dataset ordering (and reference dataset selection)'
            .type = int
        max_resolution = None
            .help = 'resolution limit for dataset to be used for TLS optimisation'
            .type = float
        max_datasets = None
            .help = 'takes up to this number of datasets for TLS parameter optimisation'
            .type = int
        simplex
            .help = 'set the various step sizes taken during simplex optimisation'
        {
            amplitude_delta_frac = 0.01
                .help = 'Amplitude step size (fractional)'
                .type = float
            amplitude_delta_min = 1e-6
                .help = 'Minimum amplitude step size (absolute)'
                .type = float
            vibration_delta = 1.0
                .help = 'Vibration step size for zero-value matrices (relative scale)'
                .type = float
            libration_delta = 10.0
                .help = 'Libration step size for zero-value matrices (relative scale)'
                .type = float
            uij_delta = 1e-3
                .help = 'Uij step size (absolute)'
                .type = float
        }
        penalties
            .help = 'penalties during optimisation.'
        {
            over_target_values
                .help = 'penalties when the fitted Uij is greater than the target Uij. penalty p is number of atoms with a Uij(fitted) > Uij(target). calculated as number of negative eigenvalues of tensor Uij(fitted)-Uij(target).'
            {
                barrier_height = 0.0005
                    .type = float
                    .help = 'manual multiplier for the sigmoid penalty function'
                barrier_width = 0.02
                    .type = float
                    .help = 'width of the buffer zone (width of buffer approx 6.9 time this value; default 0.02 -> 0.02*6.9*8*pi*pi = ~10A B-factor)'
                barrier_offset = 0.0
                    .type = float
                    .help = "Offset of the form function (e.g. inversion point of the sigmoid function)"
            }
        }
        simplex_convergence = 1e-10
            .help = "cutoff for which the least-squares simplex is considered converged"
            .type = float
        gradient_convergence = 1e-10
            .help = "cutoff for which the least-squares gradient is considered converged"
            .type = float
        amplitude_sum_weight = 1.0
            .help = "weight on the sum of the amplitudes during amplitude optimisation"
            .type = float
        normalised_tls_scale = 1e+0
            .help = "Target scale for normalisation of TLS matrices (angstroms^2)"
            .type = float
        normalised_tls_eps = 1e-3
            .help = "Cutoff for determining when a TLS-matrix has refined to zero-value (this is applied to the matrices after normalisation)."
            .type = float
        normalised_amplitude_eps = 1e-2
            .help = "Cutoff for determining when a TLS-amplitude has refined to zero-value (this is applied to the amplitudes after normalisation)."
            .type = float
    }
    precision
        .help = 'set various levels of precision for calculations'
    {
        tls_tolerance = 1e-6
            .help = "tolerance for validating TLS matrices (cutoff for defining zero when calculating negative eigenvalues, etc)"
            .type = float
        uij_tolerance = 1e-6
            .help = "tolerance for validating Uij values (maximum allowed negative eigenvalues of Uijs)"
            .type = float
        tls_matrices_decimals = 12
            .help = "how many decimals of precision to be used for TLS matrices (PDB maximum is three, but more needed here due to amplitudes less than one)"
            .type = int
        tls_amplitude_decimals = 12
            .help = "how many decimals of precision to be used for TLS amplitudes"
            .type = int
    }
}
analysis {
    refine_output_structures = True
        .help = "Refine the structures after fitting (coordinates and occupancies)"
        .type = bool
    calculate_r_factors = True
        .help = "Recalculate r-factors for the input, output (and refined) models"
        .type = bool
    calculate_electron_density_metrics = False
        .expert_level = 3
        .type = bool
    table_ones_options {
        include scope giant.jiffies.multi_table_ones.options_phil
    }
    generate_diagnostic_graphs = False
        .help = "Write diagnostic graphs that show the fit rmsd as a function of dataset/residue/atom -- may add a significant amount of runtime"
        .type = bool
}
refinement {
    program = refmac *phenix
        .help = "Should refinement be performed on the output models (coordinates and occupancies)"
        .type = choice(multi=False)
        .optional = True
    cif = None
        .help = "Cif files required for refinement"
        .type = str
        .multiple = True
}
settings {
    cpus = 1
        .type = int
    dry_run = False
        .type = bool
        .multiple = False
    verbose = False
        .type = bool
    debug = False
        .type = bool
}
""", process_includes=True)

############################################################################
#                             Utility Functions                            #
############################################################################

def format_no_fail(f_str, val):
    try:
        s = ('{'+f_str+'}').format(val)
    except:
        s = str(val)
    return s

def uij_modulus(uij_array):
    sh = uij_array.shape
    assert sh[-1] == 6
    sh_ = sh[:-1]
    uijs_flex = flex.sym_mat3_double(uij_array.reshape((numpy.product(sh_),6)))
    uijs_eigs = numpy.array(uij_eigenvalues(uijs_flex)).reshape(sh_+(3,))
    uijs_mods = numpy.mean(uijs_eigs, axis=len(sh_))
    assert uijs_mods.shape == sh_
    return uijs_mods

def translate_phenix_selections_to_pymol_selections_simple(selections, verbose=False):
    """Convert simplex phenix selections into simplex pymol selections"""

    easy_regexes = [
            ("chain '?[A-Za-z]{1,2}'?",
                lambda x: x
                ),
            ("chain '?[A-Za-z]{1,2}'? and resid '?[0-9]*?'? through '?[0-9]*'?",
                lambda x: x.replace("resid", "resi").replace(" through ",":")
                ),
                    ]

    log = ScreenLogger(stdout=verbose)

    import re
    output_selections = []
    for s in selections:
        # Start
        o = None
        # Remove double spaces and padding
        while '  ' in s:
            s = s.replace('  ',' ')
        s = s.strip(' ')
        log('Phenix Selection String: {}'.format(s))
        for rgx_s, trans_func in easy_regexes:
            rgx = re.compile(rgx_s)
            mat = rgx.findall(s)
            # No matches or too many matches
            if len(mat) != 1:
                log('> No matches or multiple matches to {}'.format(rgx_s))
                continue
            # If the match is the same as input string, then process
            if mat[0] == s:
                log('Exact match for {}'.format(rgx_s))
                o = trans_func(s)
                log('Translating to pymol: \n\t{} -> {}'.format(s, o))
                break
        # Append processed selection or none
        output_selections.append(o)

    assert len(output_selections) == len(selections)
    return output_selections

def dendrogram(fname, link_mat, labels=None, ylab=None, xlab=None, ylim=None, annotate_y_min=0.25, num_nodes=20):
    from matplotlib import pyplot
    fig = pyplot.figure()
    ax1 = pyplot.subplot(1,1,1)
    dend = scipy.cluster.hierarchy.dendrogram(link_mat, p=num_nodes, truncate_mode='lastp', labels=labels)
    # Change labels if requested
#    if labels: ax1.set_xticklabels([labels[i] for i in dend['leaves']])
    if xlab:   ax1.set_xlabel(xlab)
    if ylab:   ax1.set_ylabel(ylab)
    if ylim:   ax1.set_ylim(ylim)
    # Make sure the labels are rotated
    xlocs, xlabels = pyplot.xticks()
    pyplot.setp(xlabels, rotation=90)
    for i, d in zip(dend['icoord'], dend['dcoord']):
        x = 0.5 * sum(i[1:3])
        y = d[1]
        if y < annotate_y_min: continue
        pyplot.plot(x, y, 'ro')
        pyplot.annotate("%.3g" % y, (x, y), xytext=(0, -8), textcoords='offset points', va='top', ha='center')
    pyplot.tight_layout()
    fig.savefig(fname)
    return fig

############################################################################
#                        Multi-processing functions                        #
############################################################################

def _wrapper_run(arg):
    try:
        ret = arg.run()
    except:
        tr = traceback.format_exc()
        return tr
    return arg

def _wrapper_fit(args):
    fitter, kw_args = args
    fitter.log.toggle(0)
    try:
        msg = 'Error'
        fitter.optimise(**kw_args)
        msg = 'Finished'
        return fitter
    except Exception as e:
        tr = traceback.format_exc()
        fitter.log(tr)
        fitter.log.log_file().write_to_log(clear_data=True)
        return tr

############################################################################

def res_from_model_pdb_input(m):
    return CrystalSummary.from_pdb(pdb_input=m.input).high_res
def res_from_model_mtz(m):
    return CrystalSummary.from_mtz(m.i_mtz).high_res

class MultiDatasetUijParameterisation(Program):

    master_phil = master_phil

    def __init__(self,
                 models,
                 params,
                 levels,
                 level_labels=None,
                 log=None):
        """Object for fitting a series of TLS models to a set of structures"""

        #if log is None: log = Log()
        self.log = log

        # List of non-fatal errors from the program (to be reported at the end)
        self.warnings = []

        self.params = params

        self._n_cpu = params.settings.cpus

        # Flag to control whether input/reference is compared
        self.compare = "Input"

        self.disorder_model = None
        self.has_reflection_data = False
        self._resolution_from_model = res_from_model_pdb_input

        self.optimisation_datasets = None
        self.dataset_weights = None
        self.atom_weights = None

        self.models = models
        self.levels = levels
        self.level_labels = level_labels if level_labels else ['Level {}'.format(i) for i in range(1,len(levels)+1)]
        self.fitter = None

        # Create plot object
        self.plot = MultiDatasetUijPlots(n_levels=len(self.levels)+1)

        # Validate and add output paths, etc.
        self._init_checks()
        self._init_file_system()
        self._init_input_models()
        self._init_level_groups()
        self._init_tables()
        self._init_fitter()

        #self.write_running_parameters_to_log(params=params)

    def _init_checks(self):
        try:
            if self.params.analysis.refine_output_structures is True:
                if self.params.refinement.program == 'phenix':
                    message = 'phenix.refine is required when analysis.refine_output_structures=True and refinement.program=phenix'
                    self.check_programs_are_available(['phenix.refine'])
                elif self.params.refinement.program == 'refmac':
                    message = 'refmac5 is required when analysis.refine_output_structures=True and refinement.program=refmac'
                    self.check_programs_are_available(['refmac5'])
                else:
                    raise Sorry('Must select refinement.program when analysis.refine_output_structures=True')
            if self.params.analysis.calculate_r_factors is True:
                message = 'phenix.table_one is required when analysis.calculate_r_factors is True'
                self.check_programs_are_available(['phenix.table_one'])
            if self.params.analysis.calculate_electron_density_metrics:
                message = 'edstats (script name "edstats.pl") is required when analysis.calculate_electron_density_metrics is True'
                self.check_programs_are_available(['edstats.pl'])
            if self.params.output.images.pymol is not None:
                message = 'pymol is required when output.images.pymol is not set to "none"'
                self.check_programs_are_available(['pymol'])
        except Exception as e:
            self.log(message)
            raise

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
        table_one_cols = set(
                self.params.analysis.table_ones_options.column_labels.split(',') + \
                [self.params.analysis.table_ones_options.r_free_label]
            )

        # Flag for whether we need resolution information
        need_res_info = (self.params.fitting.optimisation.max_resolution is not None) or (self.params.fitting.optimisation.dataset_weights != 'one')

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
            m.hierarchy.write_pdb_file(
                    m.i_pdb,
                    crystal_symmetry=(m.crystal_symmetry if hasattr(m, 'crystal_symmetry') else None),
                    )
            assert os.path.exists(m.i_pdb), 'PDB does not exist: {}'.format(m.i_pdb)

            # Check MTZ exists if requested
            if self.params.input.look_for_reflection_data:
                # Look for same filename as the pdb
                mtz_file = m.filename.replace('.pdb','.mtz')
                # Check all or none of the datasets have data
                if (os.path.exists(mtz_file)):
                    # We have found reflection data -- update flags
                    if self.has_reflection_data is False:
                        # Update flag
                        self.has_reflection_data = True
                        # Update function to extract resolution from model
                        self._resolution_from_model = res_from_model_mtz
                    # Link the input mtz to the output folder
                    if self.params.output.copy_reflection_data_to_output_folder:
                        shutil.copy(mtz_file, m.i_mtz)
                    else:
                        rel_symlink(mtz_file, m.i_mtz)
                    assert os.path.exists(m.i_mtz), 'MTZ does not exist: {}'.format(m.i_mtz)
                elif (self.has_reflection_data is True):
                    errors.append(Sorry('MTZ files have been found for some datasets but not others.\n'\
                                        'To disable searching for MTZs, set input.look_for_reflection_data=False.\n'\
                                        'Missing File: {}'.format(mtz_file)))
                    continue

            # Check we have resolution information, and/or that correct columns are present in the MTZ file
            if self.has_reflection_data:
                # Extract crystal information from the MTZ
                cs = CrystalSummary.from_mtz(m.i_mtz)
                # Check for columns
                if self.params.analysis.calculate_r_factors and table_one_cols.difference(cs.column_labels):
                    errors.append(Failure("MTZ {} does not contain the correct columns.".format(m.filename) + \
                                          "\n\tLooking for: {}".format(','.join(table_one_cols)) + \
                                          "\n\tMTZ contains: {}".format(','.join(cs.column_labels)) + \
                                          "\n\tCan't find: {}".format(','.join(table_one_cols.difference(cs.column_labels))) + \
                                          "\nChange required columns with {}".format("table_ones_options.column_labels or table_ones_options.r_free_label")))

                # Check for high res information
                if (self._resolution_from_model(m) is None) and need_res_info:
                    errors.append(Failure("MTZ does not contain resolution information: {}\n".format(m.i_mtz)))
                    continue
            else:
                # Check for high res information
                if (self._resolution_from_model(m) is None) and need_res_info:
                    errors.append(Failure("PDB does not contain resolution information in the REMARK 3 records: {}\n".format(m.i_pdb) + \
                                          "This is normally present in structures that have come from refinement.\n" + \
                                          "You can also remove the high-resolution cutoff for dataset masking (optimisation.max_resolution=None),\n" + \
                                          "or use equal weighting for each dataset (optimisation.dataset_weights=one)"))
                    continue

        # Check for errors
        if errors:
            self.log.heading('Errors found while processing input files')
            for e in errors:
                self.log(str(e))
                self.log('')
            self.log.bar(False, True)
            raise Sorry("Some structures are invalid or missing information.")

        # Report dataset resolutions
        self.log('\nDataset resolutions (extracted from {} files):'.format('MTZ' if self.has_reflection_data else 'PDB'))
        for m in self.models:
            self.log('> {}: {}'.format(m.tag, self._resolution_from_model(m)))

        # Select optimisation datasets and run checks
        self.log('Selecting datasets to be used for TLS and residual optimisation')
        opt_datasets = self.select_optimisation_datasets(self.models)
        if len(opt_datasets) == 0:
            raise Sorry('No datasets selected for optimisation (e.g. above resolution cutoff: {})'.format(self.params.fitting.optimisation.max_resolution))

        # Order the datasets
        self.log('Ordering optimisation datasets by {}'.format(self.params.fitting.optimisation.sort_datasets_by))
        if self.params.fitting.optimisation.sort_datasets_by == 'resolution':
            opt_datasets = sorted(opt_datasets, key=lambda i: self._resolution_from_model(self.models[i]))
        elif self.params.fitting.optimisation.sort_datasets_by == 'name':
            opt_datasets = sorted(opt_datasets, key=lambda i: self.models[i].tag)
        elif self.params.fitting.optimisation.sort_datasets_by == 'random':
            self.log('Setting random seed: {}'.format(self.params.fitting.optimisation.random_seed))
            numpy.random.seed(self.params.fitting.optimisation.random_seed)
            opt_datasets = numpy.random.permutation(opt_datasets).tolist()
        self.log('After reordering:')
        for i_m in opt_datasets:
            self.log('\t{}: {}'.format(i_m, self.models[i_m].tag))

        # Limit the number of datasets for optimisation
        n_opt = self.params.fitting.optimisation.max_datasets
        if (n_opt is not None) and (len(opt_datasets) > n_opt):
            self.log('\nLimiting list of datasets for TLS optimisation to {} datasets'.format(n_opt))
            opt_datasets = opt_datasets[:n_opt]

        # Store & report
        self.optimisation_datasets = opt_datasets
        # Print optimisation datasets
        self.log('\nUsing {} datasets for TLS and residual parameterisation'.format(len(self.optimisation_datasets)))
        for i_m in self.optimisation_datasets:
            self.log('\t'+self.models[i_m].tag)

        # Get and print optimisation weights
        self.dataset_weights = self.calculate_dataset_weights(self.models)
        self.log('\nOptimisation weights for TLS and residual parameterisation')
        for i_m, m in enumerate(self.models):
            self.log('\t{} (model {}) -> {:+10.6f} (resolution {})'.format(m.tag, i_m, self.dataset_weights[i_m], m.input.get_r_rfree_sigma().high))

        # Update parameter flags
        if not self.has_reflection_data:
            self.log.bar(True, False)
            self.log('Reflection data has not been found in the input folder(s) -- updating settings.')
            self.log('(Added reflection data must have the same filename as input structures')
            self.log('but with different extension, e.g. <filename>.pdb & <filename>.mtz')
            self.log('Without reflection data, it is not possible to refine the fitted model or recalculate R-factors.')
            self.log('Setting analysis.refine_output_structures = False')
            self.params.analysis.refine_output_structures = False
            self.log('Setting analysis.calculate_r_factors = False')
            self.params.analysis.calculate_r_factors = False
            self.log('Setting analysis.calculate_electron_density_metrics = False')
            self.params.analysis.calculate_electron_density_metrics = False
            self.log.bar()

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

        # Create table for tracking progress over cycles
        self.tables.tracking = pandas.DataFrame(columns=['cycle', 'step', 'level',
                                                         'rmsd',
                                                         'u_iso (level)',
                                                         'b_iso (level)',
                                                         'b_min (level)',
                                                         'b_max (level)',
                                                         'u_iso (total)',
                                                         'b_iso (total)'])

        # Create output table for all model statistics
        self.tables.statistics = pandas.DataFrame(
                index   = [m.tag for m in self.models],
                columns = [
                    'High Resolution Limit',
                    'Low Resolution Limit',
                    'R-free Change (Fitted-Reference)',
                    'R-work Change (Fitted-Reference)',
                    'R-gap Change (Fitted-Reference)',
                    'R-free Change (Refined-Reference)',
                    'R-work Change (Refined-Reference)',
                    'R-gap Change (Refined-Reference)',
                    'R-free Change (Fitted-Input)',
                    'R-work Change (Fitted-Input)',
                    'R-gap Change (Fitted-Input)',
                    'R-free Change (Refined-Input)',
                    'R-work Change (Refined-Input)',
                    'R-gap Change (Refined-Input)',
                    'Average B-factor Change (Fitted-Input)',
                    'Average B-factor (fitted atoms) Change (Fitted-Input)',
                          ]
                )

        if self.params.input.reference_r_values is not None:

            self.log.subheading('Reading Reference R-values from {}'.format(self.params.input.reference_r_values))

            # Flag to control whether input/reference is compared
            self.compare = "Reference"

            # Prepare an error string
            err_str = 'A table of reference R-values has been provided (input.reference_r_values={}) '.format(self.params.input.reference_r_values)

            if not self.params.input.reference_r_values.endswith('.csv'):
                raise Sorry(err_str + 'but it must have a .csv extension.')

            # Read reference table
            ref_table = pandas.read_csv(self.params.input.reference_r_values, index_col=0)

            # Validate row names
            missing_ids = set(self.tables.statistics.index).difference(ref_table.index)
            if len(missing_ids) == len(self.tables.statistics.index):
                raise Sorry(err_str + 'but all of the dataset labels are missing. \nThe dataset labels must be present as the first column of the csv file.' + \
                        '\nCurrent datasets labels (should be present in first column in {}): \n\t{}'.format(self.params.input.reference_r_values, '\n\t'.join(self.tables.statistics.index)) + \
                        '\nCurrent first column of {}: \n\t{}'.format(self.params.input.reference_r_values, '\n\t'.join(ref_table.index)))
            elif missing_ids:
                raise Sorry(err_str + 'but some of the dataset labels are missing from the index column.' + \
                        '\nMissing labels: \n\t{}'.format('\n\t'.join(missing_ids)))

            # Validate columns
            ref_cols = ['R-free', 'R-work']
            missing_cols = set(ref_cols).difference(ref_table.columns)
            if missing_cols:
                raise Sorry(err_str + ' but columns are missing. \nMissing Columns: {}'.format(', '.join(missing_cols)))

            # Transfer R-values
            if 'R-gap' not in ref_table.columns:
                self.log('Calculating "R-gap" column for reference R-free & R-work values')
                ref_table['R-gap'] = table_one['R-free'] - table_one['R-work']
            else:
                self.log('"R-gap" column already present in reference table -- using this column.')

            for col in ['R-free', 'R-work', 'R-gap']:
                new_col = col+' (Reference)'
                if new_col in self.tables.statistics.index:
                    raise Failure('Column should not exist: {} in self.tables.statistics.index'.format(new_col))
                new_vals = ref_table[col][self.tables.statistics.index]

                assert list(new_vals) == [ref_table[col][i] for i in self.tables.statistics.index]
                self.tables.statistics[new_col] = new_vals

            self.log('\n> Processed Reference Values\n')
            for i, v in self.tables.statistics[['R-free (Reference)', 'R-work (Reference)', 'R-gap (Reference)']].iterrows():
                self.log(i)
                self.log('\t'+v.to_string().replace('\n','\n\t'))

    def _init_fitter(self):
        """Create fitting object"""

        # Extract all atoms from all datasets
        atoms = [m.hierarchy.atoms() for m in self.models]
        assert len(atoms) > 0, 'No models have been used?!'
        assert len(atoms[0]) > 0, 'No atoms have been extracted from models'

        # Extract uij and xyz from all datasets (and only select atoms we're interested in)
        observed_uij = numpy.array([a.extract_uij() for a in atoms])[:,self.atom_mask]
        observed_xyz = numpy.array([a.extract_xyz() for a in atoms])[:,self.atom_mask]

        # Variables to store details about anisotropic/isotropic atoms
        disorder_model = self.params.input.input_uij_model
        isotropic_mask = None

        self.log.subheading('Validating input uijs')

        # Check all uijs are present
        if (observed_uij==-1).all():
            # meaning: ALL uij are isotropic
            if (self.params.input.input_uij_model == 'isotropic'):
                # Extract isotropic component from atoms (again only for the atoms we're interested in)
                self.log('All atoms are missing anisotropic displacement parameters -- using the isotropic parameters instead')
                iso_b_vals_all = numpy.array([a.extract_b()/EIGHT_PI_SQ for a in atoms])[:,self.atom_mask]
                # Set Uij values of isotropic atoms
                observed_uij = numpy.zeros_like(observed_uij)
                observed_uij[:,:,0] = iso_b_vals_all
                observed_uij[:,:,1] = iso_b_vals_all
                observed_uij[:,:,2] = iso_b_vals_all
                # Set isotropic mask
                isotropic_mask = numpy.ones(observed_uij.shape[1], dtype=bool)
            else:
                raise Sorry('All atoms selected for fitting have isotropic ADPs and input.input_uij_model is currently set to "{}".\nOptions:'.format(disorder_model) + \
                        '\n   1) Re-refine structures with anisotropic ADPs.' + \
                        '\n   2) Set input.input_uij_model=isotropic to allow use of isotropic ADPs.')
        elif (observed_uij==-1).any():
            # meaning: SOME uij are isotropic
            if (self.params.input.input_uij_model == 'mixed'):
                # Find atoms that are isotropic
                isotropic_mask = (observed_uij == -1).any(axis=(0,2)).astype(bool) # Could also make 2d-mask
                # Set -1 values to 0.0
                assert (observed_uij[:,isotropic_mask,:] == -1).all()
                observed_uij[:,isotropic_mask,:] = 0.0
                # Extract B-values and convert to Uijs
                iso_b_vals_all = numpy.array([a.extract_b()/EIGHT_PI_SQ for a in atoms])[:,self.atom_mask]
                iso_b_vals_sel = iso_b_vals_all[:,isotropic_mask]
                # Set Uij values of isotropic atoms
                observed_uij[:,isotropic_mask,0] = iso_b_vals_sel
                observed_uij[:,isotropic_mask,1] = iso_b_vals_sel
                observed_uij[:,isotropic_mask,2] = iso_b_vals_sel
            else:
                raise Sorry('Some atoms for fitting (but not all) have anisotropically-refined ADPs and input.input_uij_model is currently set to "{}".\nOptions:'.format(disorder_model) + \
                        '\n   1) Re-refine structures so that all atoms selected have {} ADPs and set input.input_uij_model={}.'.format(disorder_model, disorder_model) + \
                        '\n   2) Change levels.overall_selection so that atoms with anisotropic/isotropic ADPs are not selected.' + \
                        '\n   3) Change input.input_uij_model to "mixed". This may have result in misleading results that are difficult to interpret.'
                        '\n   (Current overall selection: levels.overall_selection="{}")'.format(self.params.levels.overall_selection))
        elif self.params.input.input_uij_model not in ['anisotropic','mixed']:
            raise Sorry('input.input_uij_model is set to "{}" but all atoms selected have anisotropic uij.'.format(disorder_model))
        else:
            # meaning: All atoms are anisotropic
            assert not (observed_uij==-1).any()
            assert self.params.input.input_uij_model in ['anisotropic','mixed']
            if self.params.input.input_uij_model == 'mixed':
                self.log('Input disorder model is set to "{}", but all input atoms have anisotropic Uij.'.format(self.params.input.input_uij_model))
                disorder_model = 'anisotropic'
                self.log('Updated disorder_model to "{}"'.format(disorder_model))
            # Set isotropic mask to all False
            #self.isotropic_mask = numpy.zeros(observed_uij.shape[1], dtype=bool)

        # All -1 values should now have been removeed from input array
        assert not (observed_uij==-1).any()

        self.log('Disorder model is {}'.format(disorder_model))
        self.log('')
        if isotropic_mask is None:
            self.log('Anisotropic atoms: {}'.format(observed_uij.shape[1]))
        else:
            self.log('Anisotropic atoms: {}'.format(isotropic_mask.size-isotropic_mask.sum()))
            self.log('Isotropic atoms:   {}'.format(isotropic_mask.sum()))

        # Calculate Uij weights
        self.atom_weights = self.calculate_atom_weights(uij_values=observed_uij)

        # Renormalise the atom weights if required
        if self.params.fitting.optimisation.renormalise_atom_weights_by_dataset:
            for i,c in enumerate(self.atom_weights):
                self.atom_weights[i] = c / c.mean()

        assert abs(self.dataset_weights.mean()-1.0) < 1e-3, 'average weight should be approximately 1!'
        assert abs(self.atom_weights.mean()-1.0) < 1e-3,    'average weight should be approximately 1!'

        # Calculate combined weighting
        total_weights = self.dataset_weights.reshape((len(self.dataset_weights),1)) * self.atom_weights

        # If atom weights are pre-normalised by dataset, total weight should already be normalised
        if self.params.fitting.optimisation.renormalise_atom_weights_by_dataset:
            assert abs(total_weights.mean()-1.0) < 1e-3, 'average weight should be approximately 1!'

        # Normalise the total weights
        total_weights /= total_weights.mean()

        # Create the fitting object
        self.fitter = MultiDatasetHierarchicalUijFitter(observed_uij = observed_uij,
                                                        observed_xyz = observed_xyz,
                                                        level_array  = self.level_array,
                                                        level_labels = self.level_labels,
                                                        uij_weights = total_weights,
                                                        isotropic_mask = isotropic_mask,
                                                        params = self.params.fitting,
                                                        verbose = self.params.settings.verbose,
                                                        log = self.log)
        # Select the datasets for be used for TLS optimisation
        self.fitter.set_optimisation_datasets(self.optimisation_datasets)

        # Table to be populated during optimisation
        self.fitter.set_tracking(table    = self.tables.tracking,
                                 csv_path = self.file_manager.add_file(file_name='tracking_data.csv', file_tag='tracking_csv', dir_tag='results'),
                                 trk_path = self.file_manager.add_file(file_name='tracking_data.png', file_tag='tracking_png', dir_tag='results'),
                                 cvg_path = self.file_manager.add_file(file_name='convergence.png',   file_tag='convergence_png', dir_tag='results'),
                                 )

        # Write summary of the fitted model (groups & levels)
        self.log.subheading('Writing summary of the hierarchical model')
        self.file_manager.add_dir(dir_name='partitions', dir_tag='partitions', top_dir_tag='model')
        self.hierarchy_summary(out_dir_tag='partitions')

        # Write summary of the current settings
        self.log.subheading('Writing graphs of the different penalty functions used in optimisation')
        self.file_manager.add_dir(dir_name='penalties', dir_tag='penalties', top_dir_tag='model')
        self.penalty_summary(out_dir_tag='penalties')

        self.log.subheading('Fitter summary')
        self.log(self.fitter.summary(show=False))

        # Calculate Parameter-observed data ratio
        self.log.subheading('Data-Parameter Ratios')
        n_o = self.fitter.n_input_values()
        n_p = self.fitter.n_params(non_zero=False)
        self.log('Number of model parameters: {}'.format(n_p))
        self.log('Number of input uij values: {}'.format(n_o))
        self.log('Parameter ratio Gain (should be <1): {}'.format(float(n_p)/float(n_o)))

    def calculate_dataset_weights(self, models):
        """Calculate weightings for each model"""

        # Dataset resolution weighting function
        if self.params.fitting.optimisation.dataset_weights == 'one':
            return numpy.ones_like(models)
        elif self.params.fitting.optimisation.dataset_weights == 'inverse_resolution':
            dataset_weight = lambda r: r**(-1.0)
        elif self.params.fitting.optimisation.dataset_weights == 'inverse_resolution_squared':
            dataset_weight = lambda r: r**(-2.0)
        elif self.params.fitting.optimisation.dataset_weights == 'inverse_resolution_cubed':
            dataset_weight = lambda r: r**(-3.0)
        else:
            raise Sorry('Invalid parameter provided for fitting.optimisation.dataset_weights: {}'.format(self.params.fitting.optimisation.dataset_weights))

        # Create weights
        weights = numpy.array([dataset_weight(self._resolution_from_model(m)) for m in models])
        return weights / weights.mean()

    def calculate_atom_weights(self, uij_values):
        """Calculate weighting for each atom (in each model)"""

        # Dataset resolution weighting function
        if self.params.fitting.optimisation.atom_weights == 'one':
            atom_weight = lambda r: r**(0.0)
        elif self.params.fitting.optimisation.atom_weights == 'inverse_mod_U':
            atom_weight = lambda r: r**(-1.0)
        elif self.params.fitting.optimisation.atom_weights == 'inverse_mod_U_squared':
            atom_weight = lambda r: r**(-2.0)
        elif self.params.fitting.optimisation.atom_weights == 'inverse_mod_U_cubed':
            atom_weight = lambda r: r**(-3.0)
        else:
            raise Sorry('Invalid parameter provided for fitting.optimisation.atom_weights: {}'.format(self.params.fitting.optimisation.atom_weights))

        # Calculate the average eigenvalue for each atom
        uij_mod = uij_modulus(uij_values)
        # Apply the transformation to get the weighting
        weights = atom_weight(uij_mod)
        return weights / weights.mean()

    def select_optimisation_datasets(self, models):
        """Choose which datasets to use for e.g. shape optimisation"""
        if (self.params.fitting.optimisation.max_resolution is None):
            return range(len(models))
        selection = []
        for i_m, m in enumerate(models):
            if (self._resolution_from_model(m) < self.params.fitting.optimisation.max_resolution):
                selection.append(i_m)
        assert len(selection) > 0, 'no datasets selected for optimisation with resolution cutoff: {}'.format(self.params.fitting.optimisation.max_resolution)
        return selection

    def show_warnings(self):
        """Print non-fatal errors accumulated during running"""
        if len(self.warnings) == 0:
            self.log('Program completed with 0 errors')
            return
        banner = '{} non-fatal errors/warnings occurred during the program (see below)'.format(len(self.warnings))
        self.log.subheading(banner)
        for i, e in enumerate(self.warnings):
            self.log.bar()
            self.log('Error {} of {}'.format(i+1, len(self.warnings)))
            self.log.bar()
            self.log(e)
        self.log.bar()
        self.log.subheading(banner.replace('below','above'))

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

    def extract_results(self, sum_levels=False, sum_modes=True, average_datasets=False):
        """Return the fitted uij for each level (+residual)"""
        uij_lvl = self.fitter.extract_tls(sum_levels=sum_levels, sum_modes=sum_modes, average_datasets=average_datasets)
        uij_res = self.fitter.residual.extract(expanded=False)
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
            if os.path.exists(mdl.i_mtz) and (not os.path.exists(mdl.o_mtz)):
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
        for i_level, j_level in flex.nested_loop((len(self.levels), len(self.levels))):
            if i_level >= j_level: continue
            self.log.bar()
            self.log('Writing fitted structures (cumulative TLS components: level {:4} -> level {:4})'.format(i_level+1,j_level+1))
            self.log.bar()
            cuml_uij = uij_lvl[i_level:j_level+1].sum(axis=0)
            self.output_structures(uij = cuml_uij,
                                   iso = map(uij_to_b, cuml_uij),
                                   headers = self.make_tls_headers(datasets=None, levels=range(i_level, j_level+1)),
                                   out_dir = structure_dir,
                                   model_suffix = '.tls-level-{:04}-to-{:04}.pdb'.format(i_level+1,j_level+1))
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
                h.write_pdb_file(
                        mdl_f,
                        open_append=open_append,
                        crystal_symmetry=(mdl.crystal_symmetry if hasattr(mdl, 'crystal_symmetry') else None),
                        )
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
        if self.params.output.images.pymol is not None:
            for i_level in xrange(len(self.level_array)):
                # Images for each chain (of the partitions) - coloured by group
                self.log('\t> {}'.format(lvl_pml_prt_template.format(i_level+1,'*')))
                s_f = self.file_manager.get_file('level-partition-template').format(i_level+1)
                # Choose the style based on whether interested in atoms or larger regions
                styles = 'cartoon' if self.level_labels[i_level] in ['chain','groups','sec. struct.','residue'] else 'lines+spheres'
                # Create selections for each group in each level
                pymol_selections = translate_phenix_selections_to_pymol_selections_simple(self.levels[i_level])
                # If returned selection is none, create atom-by-atom pymol selection
                blank_h = self.blank_master_hierarchy()
                cache_h = blank_h.atom_selection_cache()
                selections = [s1 if (s1 is not None) else PymolSelection.join_or([PymolSelection.format(a) for a in blank_h.atoms().select(cache_h.selection(s2))]) for s1,s2 in zip(pymol_selections, self.levels[i_level])]
                # Create image of different selections
                auto_chain_images(structure_filename = s_f,
                                  output_prefix = lvl_pml_prt_prefix.format(i_level+1),
                                  style = styles,
                                  het_style = 'lines+spheres',
                                  colours = ['red','green','blue'],
                                  colour_selections = selections,
                                  settings = [('cartoon_oval_length', 0.5),
                                              ('cartoon_discrete_colors', 'on'),
                                              ('sphere_scale', 0.25)],
                                  width=1000, height=750,
                                  delete_script = (not self.params.settings.debug))
                if not glob.glob(lvl_pml_prt_template.format(i_level+1,'*')):
                    self.warnings.append('no partition images have been generated by pymol! ({})'.format(lvl_pml_prt_template.format(i_level+1,'*')))

    def penalty_summary(self, out_dir_tag):
        """Write out the composition of the penalty functions"""

        p = self.params.fitting.optimisation.penalties

        penalties_meta = [
                ('over_target_values',
                    Sigmoid(
                        y_scale=p.over_target_values.barrier_height,
                        x_width=p.over_target_values.barrier_width,
                        x_offset=p.over_target_values.barrier_offset),
                    'Penalty functions for eigenvalues\nof $\Delta$U=U$_{model}$-U$_{target}$',
                    'Eigenvalue of 8*$\pi^{2}$*$\Delta$U',
                    'Penalty Value',
                    (   p.over_target_values.barrier_offset-10.0*p.over_target_values.barrier_width,
                        p.over_target_values.barrier_offset+10.0*p.over_target_values.barrier_width),
                    8.*math.pi*math.pi,
                    ),
                ]

        for key, func, title, x_lab, y_lab, x_lim, x_scale in penalties_meta:
            self.log('Plotting penalty function for {}'.format(key))
            x_vals = numpy.linspace(x_lim[0], x_lim[1], 101)
            y_vals = func(x_vals)
            f_name = self.file_manager.add_file(file_name=key+'.png',
                                                file_tag=key+'_png',
                                                dir_tag=out_dir_tag)
            self.plot.lineplot(
                    filename=f_name,
                    x=x_scale*x_vals, y=y_vals,
                    title=title, x_lab=x_lab, y_lab=y_lab,
                    )

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
            tls_matrices = zip(*[[reduce(operator.add, ms) for ms in zip(*[maa.expand() for maa in f.parameters()])] for f in fitters if (f is not None)])
            #tls_matrices_1d = [reduce(operator.add, ms) for ms in zip(*tls_matrices_2d)]
            # Filter to the selection
            if datasets is not None:
                tls_matrices = [mat for i, mat in enumerate(tls_matrices) if i in datasets]
            assert len(tls_matrices) == n_dst
            assert len(tls_matrices[0]) == len(fitters)-fitters.count(None)

            # Get the origin for the datasets
            group_origins = [f for f in fitters if (f is not None)][0].atomic_com

            # Iterate through the datasets and create TLSOs for each dataset
            for i_dst, tls_list in enumerate(tls_matrices):
                # Add all the levels for each dataset
                total_mat = reduce(operator.add, tls_list)
                # Extract TLS remark lines from the
                tls_obj = mmtbx.tls.tools.tlso(t = tuple(total_mat.get('T')),
                                               l = tuple(total_mat.get('L')),
                                               s = tuple(total_mat.get('S')),
                                               origin = tuple(group_origins[i_dst]))
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

        n_tls = self.params.fitting.n_tls_modes_per_tls_group

        # ------------------------
        self.log.subheading('Writing TLS models and amplitudes for each level')
        # ------------------------
        # Output file names -- standard templates to be used for many things...
        # ------------------------
        mdl_template = fm.add_file(file_name='tls_matrices_level_{:04d}.csv',              file_tag='csv-tls-mat-template', dir_tag='tls-csvs')
        amp_template = fm.add_file(file_name='tls_amplitudes_level_{:04d}.csv',            file_tag='csv-tls-amp-template', dir_tag='tls-csvs')
        dst_template = fm.add_file(file_name='tls-model-amplitudes-level-{}-group-{}.png', file_tag='png-tls-amp-dist-template', dir_tag='tls-graphs')
        # ------------------------
        # Iterate through the levels
        for i_level, level in enumerate(self.fitter.levels):
            self.log('Level {}'.format(level.index))
            # Table for TLS model components
            mdl_filename = mdl_template.format(level.index)
            mdl_table = pandas.DataFrame(columns=["group", "mode", "label",
                                                  "T11","T22","T33","T12","T13","T23",
                                                  "L11","L22","L33","L12","L13","L23",
                                                  "S11","S12","S13","S21","S22","S23","S31","S32","S33"])
            # Create amplitude table
            amp_filename = amp_template.format(level.index)
            amp_table = pandas.DataFrame(columns=["group", "mode", "label"]+[mdl.tag for mdl in self.models], dtype=object)
            # Iterate through the groups in this level
            for i_group, (n_group, sel, fitter) in enumerate(level):
                tls_mats, tls_amps = fitter.result()
                assert tls_mats.shape == (n_tls, 21)
                assert tls_amps.shape == (n_tls, len(self.models))

                # Add to model and amplitudes tables
                for i_tls in xrange(n_tls):
                    # Add model values to last row of table
                    mdl_table.loc[len(mdl_table.index)] = numpy.concatenate([numpy.array([n_group, i_tls+1, self.levels[i_level][i_group]], dtype=object), tls_mats[i_tls]])
                    # Add amplitudes to last row of table
                    amp_table.loc[len(amp_table.index)] = numpy.concatenate([numpy.array([n_group, i_tls+1, self.levels[i_level][i_group]], dtype=object), tls_amps[i_tls,:]])

                # Write histograms of amplitudes -- only for non-zero models
                if (tls_mats.sum() > 0.0) and self.params.output.images.distributions:
                    filename = dst_template.format(level.index, i_group)
                    self.log('\t> {}'.format(filename))
                    titles = ['Mode {}:'.format(i_t+1) for i_t in xrange(n_tls)]
                    x_vals = [tls_amps[i_t,:]          for i_t in xrange(n_tls)]
                    self.plot.multi_histogram(filename  = filename,
                                              x_vals    = x_vals,
                                              titles    = titles,
                                              x_labs    = ['']*n_tls,
                                              rotate_x_labels = True,
                                              shape     = (tls_amps.shape[0], 1),
                                              n_bins    = 30, x_lim=[0, None])
            # Write tables
            self.log('\t> {}'.format(mdl_filename))
            mdl_table.to_csv(mdl_filename)
            self.log('\t> {}'.format(amp_filename))
            amp_table.to_csv(amp_filename)

            self.log.subheading('Amplitude statistics: Level {}'.format(level.index))
            # Plot histogram of amplitudes for each dataset for this level
            for mdl in self.models:
                self.log('> Model {}'.format(mdl.tag))
                amps = amp_table[mdl.tag].values.astype(float)
                mn = amps.mean()
                sd = amps.std()
                self.log('\tMean: {}'.format(mn))
                self.log('\tStd:  {}'.format(sd))
                self.log.bar()
                if len(amps) > 5:
                    try:
                        from ascii_graph import Pyasciigraph
                        g=Pyasciigraph(float_format='{0:d}')
                        counts, bounds = numpy.histogram(a=amps, bins=10)
                        graph_data = [('{:.3f}-{:.3f}'.format(bounds[i],bounds[i+1]), v) for i,v in enumerate(counts)]
                        for l in g.graph(label='Histogram of amplitudes for dataset {}'.format(mdl.tag), data=graph_data):
                            if l.startswith('#######'): continue
                            self.log(l.replace(u"\u2588", '=').replace('= ','> '))
                        self.log.bar()
                    except ImportError:
                        pass
                    except:
                        raise

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
            for i_tls in xrange(n_tls):
                self.log('- Extracting TLS contributions of model {} of level {}'.format(i_tls+1, i_level+1))
                # Create copy of the level for resetting T-L-S components
                l_copy = copy.deepcopy(level)
                # Reset the TLS for all other models
                if n_tls > 1:
                    other_models = range(n_tls)
                    other_models.remove(i_tls)
                    l_copy.zero_amplitudes(models=other_models)
                # Extract uijs
                uij_this = l_copy.extract(sum_modes=True, average_datasets=True) # TODO change this and preceding code to take advantage of new variable "sum_modes"
                # Append to total list
                uij_all.append((i_tls, uij_this))
                # Output structure for MODEL - XXX set to True for developing XXX
                if True or (n_tls > 1):
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
                                  legends       = ['TLS (Mode {})'.format(i+1) for i in mdls],
                                  title         = 'TLS contributions - Level {} ({})'.format(i_level+1, level.label),
                                  v_line_hierarchy = boundaries if (sum(boundaries.atoms().extract_b()) < 0.5*len(list(boundaries.residue_groups()))) else None,
                                  colour_indices = [float(i_level)+(float(i_mode)/float(n_tls)) for i_mode in range(n_tls)])
            if not glob.glob(lvl_plt_template.format(i_level+1,'*')):
                self.warnings.append('no plots have been generated! ({})'.format(lvl_plt_template.format(i_level+1,'*')))

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
            if self.params.output.images.pymol is not None:
                # Generate pymol images
                self.log('- Generating pymol images')
                # Images for each chain
                self.log('\t> {}'.format(lvl_pml_chn_template.format(i_level+1,'*')))
                auto_chain_images(structure_filename = m_f,
                                  output_prefix = lvl_pml_chn_prefix.format(i_level+1),
                                  style = 'lines+ellipsoids',
                                  colours = 'bfactor',
                                  width=1000, height=750)
                if not glob.glob(lvl_pml_chn_template.format(i_level+1,'*')):
                    self.warnings.append('no pymol images have been generated! ({})'.format(lvl_pml_chn_template.format(i_level+1,'*')))
            # Write pymol images of each group
            if self.params.output.images.pymol == 'all':
                # Output structure for LEVEL (scaled!) - for html only
                m_h_scl = scale_uij_to_target_by_selection(hierarchy=m_h,
                                                           selections=self.levels[i_level],
                                                           target=0.1)
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
                if not glob.glob(lvl_pml_grp_template.format(i_level+1,'*')):
                    self.warnings.append('no pymol images have been generated! ({})'.format(lvl_pml_grp_template.format(i_level+1,'*')))
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
                if not glob.glob(lvl_pml_scl_template.format(i_level+1,'*')):
                    self.warnings.append('no pymol images have been generated! ({})'.format(lvl_pml_scl_template.format(i_level+1,'*')))
            # ------------------------
            # Write out structure & graph of anisotropy
            # ------------------------
            self.log('- Calculating anisotropy of Uijs')
            ani = 1.0 - calculate_uij_anisotropy_ratio(uij=uij)
            m_h = self.custom_master_hierarchy(uij=None, iso=ani, mask=atom_sel)
            # Make plot of the uij anisotropy over the chain
            self.log('\t> {}'.format(lvl_aniso_plt_template.format(i_level+1,'*')))
            self.plot.stacked_bar(prefix=lvl_aniso_plt_prefix.format(i_level+1),
                                  hierarchies=[m_h.select(atom_sel)],
                                  legends=['Anisotropy  '],
                                  title='Anisotropy of Level {} ({})\nfully isotropic -> 0 (spheres)\nfully anisotropic -> 1 (lines/disks)'.format(i_level+1, level.label),
                                  y_lim=(0.0,1.05),
                                  y_lab='Anisotropy of Uij ($1 - \\frac{E_{min}}{E_{max}}$)',
                                  colours=['grey'])
            if not glob.glob(lvl_aniso_plt_template.format(i_level+1,'*')):
                self.warnings.append('no plots have been generated! ({})'.format(lvl_aniso_plt_template.format(i_level+1,'*')))

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
                              colour_indices=[len(self.levels)])
        if not glob.glob(res_plt_template.format('*')):
            self.warnings.append('no files have been generated! ({})'.format(res_plt_template.format('*')))
        # Write pymol images for each chain
        if self.params.output.images.pymol is not None:
            self.log('- Generating pymol images')
            self.log('\t> {}'.format(res_pml_chn_template.format('*')))
            auto_chain_images(structure_filename = m_f,
                              output_prefix = res_pml_chn_prefix,
                              style = 'lines+ellipsoids',
                              colours = 'bfactor',
                              width=1000, height=750)
            if not glob.glob(res_pml_chn_template.format('*')):
                self.warnings.append('no files have been generated! ({})'.format(res_pml_chn_template.format('*')))
        # Write pymol images for each group
        if self.params.output.images.pymol == 'all':
            self.log('\t> {}'.format(res_pml_grp_template.format('*')))
            auto_residue_images(structure_filename = m_f,
                                output_prefix = res_pml_grp_prefix,
                                style = 'lines+ellipsoids',
                                width=250, height=250)
            if not glob.glob(res_pml_grp_template.format('*')):
                self.warnings.append('no files have been generated! ({})'.format(res_pml_grp_template.format('*')))
        # ------------------------
        # Write out structure & graph of anisotropy
        # ------------------------
        self.log('- Calculating anisotropy of Uijs')
        ani = 1.0 - calculate_uij_anisotropy_ratio(uij=uij)
        m_h = self.custom_master_hierarchy(uij=None, iso=ani, mask=atom_sel)
        # Make plot of the uij anisotropy over the chain
        self.log('\t> {}'.format(res_aniso_plt_template.format('*')))
        self.plot.stacked_bar(prefix=res_aniso_plt_prefix,
                              hierarchies=[m_h.select(atom_sel)],
                              legends=['Anisotropy  '],
                              title='Anisotropy of Residual Level\nfully isotropic -> 0 (spheres)\nfully anisotropic -> 1 (lines/disks)',
                              y_lim=(0.0,1.05),
                              y_lab='Anisotropy of Uij ($1 - \\frac{E_{min}}{E_{max}}$)',
                              colours=['grey'])
        if not glob.glob(res_aniso_plt_template.format('*')):
            self.warnings.append('no files have been generated! ({})'.format(res_aniso_plt_template.format('*')))

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
        # One hierarchy for the average input structure
        h_inp = self.custom_master_hierarchy(iso=uij_to_b(self.fitter.observed_uij.mean(axis=0)), mask=atom_sel)
        # Create stacked plot
        self.log('\t> {}'.format(stacked_template.format('*')))
        self.plot.stacked_bar(prefix=stacked_prefix,
                              hierarchies=[h.select(atom_sel) for h in uij_hierarchies],
                              legends=[l.capitalize() for l in self.level_labels]+['Residual'],
                              reference_hierarchy=h_inp.select(atom_sel),
                              reference_legend='Input',
                              title='TLS and residual contributions',
                              reverse_legend_order=True,
                              colour_indices=range(len(self.levels)+1))
        if not glob.glob(stacked_template.format('*')):
            self.warnings.append('no files have been generated! ({})'.format(stacked_template.format('*')))

    def calculate_electron_density_metrics(self, out_dir_tag):
        pass

    def calculate_amplitudes_dendrogram(self, out_dir_tag):
        """Cluster the amplitudes across the datasets"""

        fm = self.file_manager
        out_dir = fm.get_dir(out_dir_tag)
        dendro_template = fm.add_file(file_name='amplitude_dendrogram-level_{}.png', file_tag='png-amplitude-dendrogram-template', dir_tag=out_dir_tag)

        self.log.heading('Generating dendrograms of TLS amplitudes between datasets')

        # List of linkage-distances for each level
        level_distance_matrices = []

        labels = [m.tag for m in self.models]

        all_amplitudes = []
        for level in self.fitter.levels:
            # List of amplitudes for each group in this level
            group_amplitudes = []
            for i_group, sel, fitter in level:
                tls_mats, tls_amps = fitter.result()
                group_amplitudes.append(tls_amps)
            group_amplitudes = numpy.array(group_amplitudes)

            n_grp, n_tls, n_dst = group_amplitudes.shape
            self.log('For level {}:'.format(level.index))
            self.log('\tNumber of groups: {}'.format(n_grp))
            self.log('\tNumber of TLS: {}'.format(n_tls))
            self.log('\tNumber of datasets: {}'.format(n_dst))

            amp_mat = group_amplitudes.reshape((n_grp*n_tls, n_dst)).T
            all_amplitudes.append(amp_mat)

            # Plot dendrogram
            fname = dendro_template.format(level.index)
            self.log('> {}'.format(fname))
            link_mat = scipy.cluster.hierarchy.linkage(amp_mat, method='average', metric='euclidean')
            dendrogram(fname=fname, link_mat=link_mat, labels=labels)

        # Plot dendrogram
        fname=dendro_template.format('all')
        self.log('> {}'.format(fname))
        link_mat = scipy.cluster.hierarchy.linkage(numpy.concatenate(all_amplitudes, axis=1), method='average', metric='euclidean')
        dendrogram(fname=fname, link_mat=link_mat, labels=labels)

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

    def analyse_fit_rms_distributions(self, uij_fit, uij_inp, out_dir='./', max_x_width=40):
        """Analyse the dataset-by-dataset and residue-by-residue and atom-by-atom fit qualities"""

        out_dir = easy_directory(out_dir)

        #------------------------------------------------------------------------------#
        # Create output directories
        #------------------------------------------------------------------------------#
        res_dir = easy_directory(os.path.join(out_dir, 'residue_by_residue'))
        dst_dir = easy_directory(os.path.join(out_dir, 'dataset_by_dataset'))
        atm_dir = easy_directory(os.path.join(out_dir, 'atom_by_atom'))

        #------------------------------------------------------------------------------#
        # Extract series of labels for datasets, atoms ... extract rmsd values
        #------------------------------------------------------------------------------#

        # ... dataset labels
        dst_labels = numpy.array([m.tag for m in self.models])
        # ... atom labels
        fit_sel = flex.bool(self.atom_mask.tolist())
        m_h = self.blank_master_hierarchy().select(fit_sel, copy_atoms=True)
        atm_labels = numpy.array([ShortLabeller.format(a) for a in m_h.atoms()])
        res_labels = numpy.array([ShortLabeller.format(a.parent().parent()) for a in m_h.atoms()])

        # Average over uij sym_mat3
        uij_diff = uij_inp-uij_fit
        rmsd_all = rms(uij_diff, axis=2)

        y_lim = (0.0, 1.01*numpy.max(rmsd_all))

        #------------------------------------------------------------------------------#
        # Plot rmsds over datasets (compare datasets)
        #------------------------------------------------------------------------------#

        for i_dst in xrange(0, len(self.models), max_x_width):
            filename = os.path.join(dst_dir, 'rmsds_datasets-{:04d}-{:04d}.png'.format(i_dst+1, i_dst+max_x_width))
            self.log('\t> {}'.format(filename))
            self.plot.violinplot(filename=filename,
                                 y_vals=rmsd_all.T[:,i_dst:i_dst+max_x_width],
                                 x_labels=dst_labels[i_dst:i_dst+max_x_width].tolist(),
                                 title='rmsd of fitted and refined B-factors (by dataset)',
                                 x_lab='dataset', y_lab='rmsd', rotate_x_labels=True,
                                 x_lim=(0, min(len(self.models), max_x_width)+1), y_lim=y_lim)

        #------------------------------------------------------------------------------#
        # Plot rmsds by residue over datasets (compare residues)
        #------------------------------------------------------------------------------#

        self.log.subheading('Calculating dataset-by-dataset rmsds to input Uijs')

        # Format rmsds as hierarchy
        rmsd_hs = [self.custom_master_hierarchy(uij=None, iso=rmsd_all[i].tolist(), mask=fit_sel) for i in xrange(len(self.models))]

        #############################
        # Plot all datasets
        #############################

        # Plot all
        prefix = os.path.join(res_dir, 'rmsds_all_datasets')
        self.plot.alpha_scatter(
                prefix=prefix,
                hierarchies=[h.select(fit_sel, copy_atoms=True) for h in rmsd_hs],
                title='All RMSDs',
                y_lab='RMSD ($\AA^2$)',
                y_lim=y_lim,
                rotate_x_labels=True,
                plot_mean=True,
                )
        # Write as structure
        filename = prefix + '.pdb'
        all_h = self.custom_master_hierarchy(uij=None, iso=rmsd_all.mean(axis=0).tolist(), mask=fit_sel)
        all_h.write_pdb_file(
                filename,
                )

        #############################
        # Plot dataset by dataset
        #############################

        for i_dst, mdl in enumerate(self.models):
            # Plot scatter (with average)
            prefix = os.path.join(res_dir, 'rmsds_dataset_{}'.format(mdl.tag))
            self.plot.alpha_scatter(
                    prefix=prefix,
                    hierarchies=[rmsd_hs[i_dst].select(fit_sel, copy_atoms=True)],
                    title='Dataset {} RMSDs'.format(mdl.tag),
                    y_lab='RMSD ($\AA^2$)',
                    y_lim=y_lim,
                    rotate_x_labels=True,
                    plot_mean=True,
                    )
            # Write as structure
            filename = prefix + '.pdb'
            rmsd_hs[i_dst].write_pdb_file(
                    filename,
                    crystal_symmetry=(mdl.crystal_symmetry if hasattr(mdl, 'crystal_symmetry') else None),
                    )

        #------------------------------------------------------------------------------#
        # Atom-by-atom and residue-by-residue "error" boxplots
        #------------------------------------------------------------------------------#

        self.log.subheading('Calculating atom-by-atom and residue-by-residue rmsds to input Uijs')
        # Use this for the atom-by-atom graphs with atom_label=res_labels
        residue_boundaries = self.partition_boundaries_custom(atom_labels=res_labels, mask=fit_sel).select(fit_sel)
        # Create boxplots of the rmsds by atom (for the first level only!)
        for chain_id in sorted(set([c.id for c in m_h.chains()])):
            self.log.bar()
            self.log('Chain {}'.format(chain_id))
            self.log.bar()
            # Create a selection for each chain and select it in each hierarchy (m_h is already trimmed to fitted atoms)
            sel = numpy.array(m_h.atom_selection_cache().selection('chain {}'.format(chain_id)), dtype=bool)
            # Extract the uij rmsds for this group
            rmsd_sel = rmsd_all[:,sel]
            # Create boxplots of the rmsds by atom
            self.log('Atom-by-atom plots')
            for i_a in xrange(0, sum(sel), max_x_width):
                filename = os.path.join(atm_dir, 'rmsds-chain_{}-atoms-{:06d}-{:06d}.png'.format(chain_id, i_a+1, i_a+max_x_width))
                self.log('\t> {}'.format(filename))
                self.plot.violinplot(filename=filename,
                                     y_vals=rmsd_sel[:,i_a:i_a+max_x_width],
                                     x_labels=atm_labels[sel].tolist()[i_a:i_a+max_x_width],
                                     title='rmsd of fitted and refined B-factors (by atom)',
                                     x_lab='atom', y_lab='rmsd', rotate_x_labels=True,
                                     x_lim=(0,max_x_width+1), y_lim=y_lim,
                                     vlines=(numpy.where(numpy.array(list(residue_boundaries.select(flex.bool(sel.tolist())).atoms().extract_b()[i_a:i_a+max_x_width]), dtype=bool))[0] + 1.5))

    def analyse_fit_residual_correlations(self, uij_diff, out_dir='./', max_x_width=25):
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
        fit_sel = flex.bool(self.atom_mask.tolist())
        m_h = self.blank_master_hierarchy().select(fit_sel, copy_atoms=True)
        atm_labels = numpy.array([ShortLabeller.format(a) for a in m_h.atoms()])
        res_labels = numpy.array([ShortLabeller.format(a.parent().parent()) for a in m_h.atoms()])
        mdl_labels = numpy.array([mdl.tag for mdl in self.models])

        # Extract residue boundaries
        residue_boundaries = self.partition_boundaries_custom(atom_labels=res_labels, mask=fit_sel).select(fit_sel)
        # Make violin plot of the correlations
        for chain_id in sorted(set([c.id for c in m_h.chains()])):
            self.log.bar()
            self.log('Chain {}'.format(chain_id))
            self.log.bar()
            # Create a selection for each chain and select it in each hierarchy (m_h is already trimmed to fitted atoms)
            sel = numpy.array(m_h.atom_selection_cache().selection('chain {}'.format(chain_id)), dtype=bool)
            # Extract the uij rmsds for this group
            corr_sel = corr[:,sel]
            # Create boxplots of the rmsds by atom
            self.log('Atom-by-atom plots')
            for i_a in xrange(0, sum(sel), max_x_width):
                filename = os.path.join(cor_dir, 'residual-correction_chain-{}_atoms-{:06d}-{:06d}.png'.format(chain_id, i_a+1, i_a+max_x_width))
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

        if self.params.refinement.program == 'refmac':
            refine_ = refine_refmac
        elif self.params.refinement.program == 'phenix':
            refine_ = refine_phenix
        else:
            raise Failure('Invalid refinement program selected: {}'.format(self.params.refinement.program))

        self.log('Preparing refinements -- will use {} for refinement'.format(refine_.program))

        proc_args = []
        for mdl in self.models:
            if not os.path.exists(mdl.i_mtz):
                raise Failure('Something has gone wrong -- Trying to refine output structures but mtz has not been provided/has been lost.')
            if not os.path.exists(mdl.o_mtz):
                rel_symlink(mdl.i_mtz, mdl.o_mtz)
            # Create refinement object
            obj = refine_(
                    pdb_file=mdl.o_pdb,
                    mtz_file=mdl.o_mtz,
                    cif_files=(self.params.refinement.cif if self.params.refinement.cif else None),
                    out_prefix=os.path.splitext(mdl.o_pdb)[0]+suffix,
                    strategy='individual_sites+occupancies',
                    n_cycles=5,
                    log=self.log).setup()
            obj.log('')
            obj.print_settings()
            obj.tag = mdl.tag
            if os.path.exists(obj.out_pdb_file):
                raise Failure('Refined PDB already exists! (model {})'.format(mdl.tag))
            proc_args.append(obj)

        # Refine all of the models
        self.log.subheading('Running refinements')
        refined = libtbx.easy_mp.pool_map(processes=self._n_cpu, func=_wrapper_run, args=proc_args, chunksize=1)

        for mdl, ref in zip(self.models, refined):
            if isinstance(ref, str):
                self.log('Failed in refinement -- continuing')
                continue
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
            if not os.path.exists(mdl.i_mtz):
                raise Failure('Something has gone wrong -- Trying to generate table ones but mtz has not been provided/has been lost.')
            if not os.path.exists(mdl.o_mtz):
                rel_symlink(mdl.i_mtz, mdl.o_mtz)

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
        phil.options          = self.params.analysis.table_ones_options
        phil.settings.cpus    = self.params.settings.cpus
        phil.settings.verbose = False

        # Run 1
        phil.input.pdb = [mdl.i_pdb for mdl in self.models]
        phil.input.labelling  = 'foldername'
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
        if self.params.analysis.refine_output_structures:
            phil.input.pdb = [mdl.r_pdb for mdl in self.models if mdl.r_pdb is not None]
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

        try:
            assert os.path.exists(table_one_csv_orig)
            assert os.path.exists(table_one_csv_fitd)
            if self.params.analysis.refine_output_structures:
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

        # Extract output table
        main_table = self.tables.statistics

        # columns to be extracted from each file
        input_cols = ['High Resolution Limit', 'Low Resolution Limit', 'Unique reflections','Completeness (%)','Wilson B-factor']
        other_cols = ['R-work','R-free','Average B-factor']
        other_cols_sort = ['R-free','R-work','R-gap','Average B-factor']

        # Extract data from the table one CSVs
        self.log.subheading('Looking for table one data')
        for suff, csv_lab in [(' (Input)',   'table_one_input'),
                              (' (Fitted)',  'table_one_fitted'),
                              (' (Refined)', 'table_one_refined')]:
            if not self.file_manager.has_file(csv_lab):
                self.log('CSV not generated: {} -- skipping'.format(csv_lab))
                continue
            csv = self.file_manager.get_file(csv_lab)
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
                assert (main_table.index == table_one.index).all()
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
            for suff1, suff2 in [
                    (' (Input)', ' (Fitted)'),
                    (' (Input)', ' (Refined)'),
                    (' (Reference)', ' (Fitted)'),
                    (' (Reference)', ' (Refined)'),
                                ]:
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

        # Store output table (object is changed in the join steps)
        self.tables.statistics = main_table

        # Remove unpopulated columns and write
        self.log('')
        self.log('Writing output csv: {}'.format(out_csv))
        main_table.dropna(axis='columns', how='all').to_csv(out_csv)

        if self.params.analysis.calculate_r_factors:
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
            if lab.lower().endswith('-reference)'): continue
            if lab.lower().endswith('(input)'): continue
            if lab.lower().endswith('(fitted)'): continue
            if lab.lower().endswith('(refined)'): continue
            if lab.lower().endswith('(reference)'): continue
            self.log(shrt_str.format(lab, table_means[lab]))
        self.log.bar()
        # Columns with old/new prefix
        for lab1, lab2 in itertools.product(
                                    ['Input', 'Reference'],
                                    ['Fitted','Refined']
                                            ):
            lab1_suff = ' ({})'.format(lab1)
            lab2_suff = ' ({})'.format(lab2)
            # Check whether the columns are here - check the R-free as always should be present
            if ('R-free ({})'.format(lab1) not in table.columns) or \
                    ('R-free ({})'.format(lab2) not in table.columns):
                continue
            self.log('Comparing {} and {} Uijs'.format(lab1, lab2))
            self.log(long_str.format('', lab1, lab2, 'Difference'))
            for lab in table_means.index:
                # Only select labs for the first label
                if not lab.endswith(lab1_suff): continue
                # Extract common base component of label
                lab_base = lab[:lab.index(lab1_suff)]
                f_lab1 = lab_base+lab1_suff
                f_lab2 = lab_base+lab2_suff
                c_lab = lab_base+' Change ({}-{})'.format(lab2, lab1)
                # Check other label is also in table
                if f_lab1 not in table_means.index: continue
                if f_lab2 not in table_means.index: continue
                # If this is the case then this should be true
                assert c_lab in table_means.index
                # Print
                self.log(long_str.format(lab_base,
                                         '{:.3f}'.format(table_means[f_lab1]),
                                         '{:.3f}'.format(table_means[f_lab2]),
                                         '{:.3f}'.format(table_means[c_lab])))
            self.log.bar()

        return

    def write_table_one_graphs(self, table, results_dir_tag, analysis_dir_tag):
        """Look at pre- and post-fitting graphs"""

        fm = self.file_manager
        result_dir = fm.get_dir(results_dir_tag)
        analys_dir = fm.get_dir(analysis_dir_tag)

        self.log.subheading('Model improvement graphs')

        labs = ['Input']     * ('R-free (Input)'     in table.columns) + \
               ['Fitted']    * ('R-free (Fitted)'    in table.columns) + \
               ['Refined']   * ('R-free (Refined)'   in table.columns) + \
               ['Reference'] * ('R-free (Reference)' in table.columns)

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
                                 title = 'Structure R-free values',
                                 x_lab = 'Resolution ($\AA$)',
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
                                 title = 'Structure R-work values',
                                 x_lab = 'Resolution ($\AA$)',
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
                                 title = 'Structure R-gap values',
                                 x_lab = 'Resolution ($\AA$)',
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
                                 title = 'R-value change from input to fitted models',
                                 x_lab = 'Resolution ($\AA$)',
                                 y_lab = 'R-value change (%)',
                                 rotate_x_labels = True,
                                 min_bin_width = 0.1,
                                 hlines = [0.0])
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
            main_log = os.path.abspath(self.log.log_file().path)
            for log in glob.glob(os.path.join(self.file_manager.get_dir('logs'),'*.log')):
                # Do not compress main log file
                if os.path.abspath(log) == main_log:
                    continue
                if self.params.settings.verbose:
                    self.log('Compressing: {}'.format(log))
                zip_file = compress_file(log, delete_original=True)

class MultiDatasetUijPlots(object):

    def __init__(self, n_levels, cmap_name='rainbow'):
        self.n_levels = n_levels
        self.cmap_name = cmap_name

    def get_colour_map(self):
        return matplotlib.cm.get_cmap(self.cmap_name)

    def get_level_colours_arbitrary(self, indices):
        cm = self.get_colour_map()
        return cm(numpy.array(indices)/float(self.n_levels-1))

    def get_level_colours(self):
        cm = self.get_colour_map()
        return cm(numpy.linspace(0., 1., self.n_levels))

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
    def lineplot(cls, filename,
                 x, y,
                 title='', x_lab='x', y_lab='y'):

        fig, axis = pyplot.subplots(nrows=1, ncols=1)
        axis.set_title(title)
        axis.plot(x, y, 'b-')
        axis.set_xlabel(x_lab)
        axis.set_ylabel(y_lab)
        fig.tight_layout()
        fig.savefig(filename)
        pyplot.close(fig)

        return

    @classmethod
    def binned_boxplot(cls, filename,
                       x, y=None, y_vals=None, legends=None,
                       title='', x_lab='x', y_lab='y',
                       rotate_x_labels=True,
                       plot_scatter=True,
                       max_bins=10, min_bin_width=None,
                       hlines=[], vlines=[],
                       plot_type='violinplot',
                       ):
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
        bar_width = 2./float(1+3*n_y) # half-bar padding between bars
        # Colours of each of y_vals
        colours = pyplot.cm.rainbow(numpy.linspace(0,1,n_y))

        # Create figures
        fig, axis = pyplot.subplots(nrows=1, ncols=1)
        axis.set_title(title)
        # Draw horizontal/vertical lines (do first so they're at the bottom)
        for v in hlines: axis.axhline(y=v, linewidth=1, zorder=1)
        for v in vlines: axis.axvline(x=v, linewidth=1, zorder=1)
        # Store plot objects
        plot_dicts = []
        for i_y, y in enumerate(binned_y):
            # Offset for this bar set relative to (n-1)
            x_offset = 0.5 + (1+1.5*i_y)*bar_width
            positions = numpy.arange(n_bins) + x_offset
            # Filter on values that are actually present
            y_idx = [i for i in xrange(len(y)) if len(y[i])>0]
            # Plot
            if plot_type == 'boxplot':
                plt = axis.boxplot(
                        y,
                        positions=positions,
                        widths=bar_width,
                        showmeans=False,
                        patch_artist=True)
            else:
                plt = axis.violinplot(
                        [y[i] for i in y_idx],
                        positions=[positions[i] for i in y_idx],
                        widths=bar_width,
                        showmeans=True,
                        showextrema=True)
            plot_dicts.append(plt)
            # Color the boxplots
            c = colours[i_y]
            for obj in ['boxes','bodies']:
                for patch in plt.get(obj, []):
                    patch.set_facecolor(c)
                    patch.set_edgecolor('k')
            for obj in ['means','medians']:
                for line in plt.get(obj, []):
                    line.set_facecolor(c)
                    line.set_edgecolor('k')
            for obj in ['cbars','cmeans','cmedians','cmins','cmaxes']:
                lc = plt.get(obj, None)
                if lc is not None:
                    lc.set_facecolor(c)
                    lc.set_edgecolor('k')
            # Plot a superposed scatter plot
            if plot_scatter is True:
                for x, ys in zip(positions, y):
                    if len(ys) == 0: continue
                    xs = [x]*len(ys)
                    jitter = numpy.random.randn(len(ys)) * 0.5 * bar_width / 3.0 # Normalise to be approx width of bar
                    jitter = jitter - jitter.mean()
                    axis.scatter(x=xs+jitter, y=ys, s=10, facecolor=c, edgecolor='k', zorder=10)

        # Make labels
        bin_centres = numpy.arange(n_bins)+1
        axis.set_xticks(bin_centres)
        axis.set_xticklabels(bin_labels)
        # Axis labels
        axis.set_xlabel(x_lab)
        axis.set_ylabel(y_lab)
        # Axis limits
        axis.set_xlim((0.25, n_bins+0.75))
        # Plot v_lines
        if n_y > 1:
            for i in xrange(n_bins+1):
                axis.axvline(i+0.5, c='grey', ls="dashed", lw=1.0)#ls="solid", lw=0.5)
        # X-axis rotations
        if rotate_x_labels:
            pyplot.setp(axis.get_xticklabels(), rotation=45)
        # Make legend
        if legends is not None:
            #handles = [patches.Patch(edgecolor='k', facecolor=colours[i_y], label=legends[i_y]) for i_y in xrange(n_y)]
            #lgd = axis.legend(handles=handles,
            #                  bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            #extra_artists.append(lgd)
            # Plot additional text on the graph
            y_min, y_max = axis.get_ylim()
            fig_scale = (axis.transData.transform((1.0,0.)) - axis.transData.transform((0.0,0.)))[0]
            fig_width = 30.0 * n_y # final width of the fonts in the scale of the image
            x_width = fig_width / fig_scale # Width on the x-axis scale
            t_width = x_width / float(n_y)
            fontsize = 10
            text_artists = []
            max_l_length = max(map(len,legends))
            for i_y, l in enumerate(legends):
                y_pos = (0.5*y_min+0.5*y_max) #(0.05*(y_max-y_min) if (y_max>0.0>y_min) else (0.5*y_min+0.5*y_max))
                t = axis.text(
                        x=(n_bins+0.5)+(0.5+i_y)*t_width, y=y_pos, s=l.center(max_l_length),
                        bbox=dict(boxstyle='round', facecolor=colours[i_y], edgecolor='k', linewidth=0.5, alpha=0.75, pad=0.3),
                        fontsize=fontsize, rotation=90, ha='center', va='center', zorder=1)
                text_artists.append(t)
            axis.set_xlim(right=n_bins+0.5+x_width)

        fig.tight_layout()
        fig.savefig(filename,
                    bbox_extra_artists=extra_artists+text_artists,
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

        # Colours of each of y_vals
        if n_plot == n_line:
            colours = ['b']*n_line
        elif (n_line % n_plot) == 0:
            colours = sorted(numpy.concatenate([pyplot.cm.rainbow(numpy.linspace(0,1,n_line//n_plot))]*n_plot).tolist())
        else:
            colours = pyplot.cm.rainbow(numpy.linspace(0,1,n_line))
        assert len(colours)==n_line

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
                               c=colours[i])
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
        #if rotate_x_labels:
        #    pyplot.setp(axis.get_xticklabels(), rotation=45)
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

    def level_plots(self, filename, hierarchies, title, rotate_x_labels=True):
        """Plot a schematic representation of the hierarchical partitioning"""

        colours = self.get_level_colours()

        fig, axis = pyplot.subplots()
        axis.set_facecolor('w')

        assert len(hierarchies) == (self.n_levels-1)

        for i_h, h in enumerate(hierarchies):
            # Extract B-factors and atom labels
            b_vals = numpy.array(h.atoms().extract_b())
            # Extract colour for this level
            col = colours[i_h]
            # Plot a bar around the atoms (below other bars)
            axis.broken_barh([(-0.25, len(b_vals)-1.0+0.50)], (i_h+0.5, 1.0), edgecolor='k', facecolor='w')
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
                    plot_vals.append((min_i-0.25, max_i-min_i+0.50))
                # Plot
                axis.broken_barh(plot_vals, (i_h+0.5, 1.0), edgecolor='k', facecolor=col)

        # Set axes, etc.
        if title is not None: axis.set_title(label=str(title))
        axis.set_xlabel('Atom')
        axis.set_ylabel('')
        a_labels = numpy.array([ShortLabeller.format(a.parent().parent()) for a in hierarchies[0].atoms()])
        #axis.set_xticklabels([a_labels[int(i)] if ((i>0) and (i<len(a_labels)) and (float(int(i))==i)) else '' for i in axis.get_xticks()])
        axis.set_xticks(range(0, len(a_labels), max(1, len(a_labels)//40)))
        axis.set_xticklabels([a_labels[i] for i in axis.get_xticks()])
        axis.set_yticks(range(1, len(hierarchies)+1))
        axis.set_yticklabels(['Level {}'.format(i) for i in xrange(1, len(hierarchies)+1)])
        axis.set_xlim((-0.5, len(a_labels)-0.5))
        pyplot.setp(axis.get_xticklabels(), rotation=90)
        # Invert y-axis
        #axis.invert_yaxis()
        # Format and save
        fig.tight_layout()
        fig.savefig(filename, dpi=500)
        pyplot.close(fig)

    def stacked_bar(self,
                    prefix,
                    hierarchies,
                    legends,
                    reference_hierarchy=None,
                    reference_legend=None,
                    title=None,
                    y_lab='Isotropic B ($\AA^2$)',
                    y_lim=None,
                    v_line_hierarchy=None,
                    rotate_x_labels=True,
                    reverse_legend_order=False,
                    colour_space=(0,1),
                    colour_indices=None,
                    colours=None):
        """Plot stacked bar plots for a series of hierarchies (plotted values are the average B-factors of each residue of the hierarchies)"""

        # TODO MERGE stacked_bar and multi_bar with mode='stacked' mode='grid'

        legends = list(legends)
        assert len(hierarchies) == len(legends)

        m_h = hierarchies[0]

        if colour_indices is not None:
            assert len(hierarchies) == len(colour_indices)
            colours = self.get_level_colours_arbitrary(colour_indices)
        elif colours is not None:
            assert len(hierarchies) == len(colours)
        else:
            cm = self.get_colour_map()
            colours = cm(numpy.linspace(colour_space[0],colour_space[1],len(hierarchies)))

        # Check all hierarchies are the same
        for h in hierarchies:
            assert m_h.is_similar_hierarchy(h)
        if reference_hierarchy is not None:
            assert m_h.is_similar_hierarchy(reference_hierarchy)

        # Create a plot for each chain
        for chain_id in sorted(set([c.id for c in m_h.chains()])):

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

            # Create colours + hatches
            #hatchs = itertools.cycle(['//', 'x', '\\'])
            hatchs = itertools.cycle([None])

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
                hdl = axis.bar(
                        left=x_vals,
                        height=y_vals,
                        bottom=cuml_y,
                        width=1.0,
                        align='center',
                        color=colours[i_h],
                        label=legends[i_h],
                        hatch=hatchs.next(),
                        linewidth=0,
                        edgecolor='black',
                        zorder=5,
                        )
                handles.append(hdl)

                # Append to cumulative y
                cuml_y += y_vals

            # Plot lines at top of graph for prettiness/as reference hierarchy
            if (reference_hierarchy is not None) or (len(sel_hs)==1):
                if reference_hierarchy is not None:
                    line_h_sel = reference_hierarchy.select(sel)
                else:
                    line_h_sel = sel_hs[0]
                # Extract y-values
                y_vals = numpy.array([numpy.mean(rg.atoms().extract_b()) for rg in line_h_sel.residue_groups()])
                # Duplicate the x-vals and y-vals and shift x-vals apart by half a unit
                x_vals_dup = numpy.concatenate([(x_vals-0.5),(x_vals+0.5)]).reshape((2,x_vals.size)).T.reshape(2*x_vals.size).tolist()
                y_vals_dup = numpy.concatenate([y_vals,y_vals]).reshape((2,y_vals.size)).T.reshape(2*y_vals.size).tolist()
                # Add another point at the beginning and end so starts and ends on the baseline
                x_vals_dup = numpy.concatenate([[x_vals_dup[0]], x_vals_dup, [x_vals_dup[-1]]])
                y_vals_dup = numpy.concatenate([[0.0], y_vals_dup, [0.0]])
                hdl = axis.plot(x_vals_dup, y_vals_dup, 'k-', label=reference_legend, lw=0.5, zorder=5)
                if reference_hierarchy is not None:
                    handles.extend(hdl)

            # Legends (reverse the plot legends!)
            if reverse_legend_order is True:
                handles.reverse()

            # No lines if cuml_y is None
            if cuml_y is None:
                continue

            # Plot boundaries
            if v_line_hierarchy is not None:
                v_lines = numpy.where(numpy.array([max(rg.atoms().extract_b()) for rg in v_line_hierarchy.select(sel).residue_groups()], dtype=bool))[0] + 1.5
                for val in v_lines:
                    axis.axvline(x=val, c='grey', ls='solid', label='boundaries', lw=0.5, zorder=1)
                # Add a line at zero
                h = axis.axvline(x=0.5, c='grey', ls='solid', label='boundaries', lw=0.5, zorder=1)
                handles.append(h)

            # Plot legend
            lgd = axis.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            # Axis ticks & labels
            x_ticks = numpy.arange(1, len(x_labels)+1, int(max(1.0, numpy.floor(len(x_labels)/20))))
            #x_ticks = numpy.unique(numpy.floor(numpy.linspace(1,len(x_labels)+1,20)))
            axis.set_xticks(x_ticks)
            axis.set_xticklabels([x_labels[int(i)] if (i<len(x_labels)) and (float(int(i))==i) else '' for i in axis.get_xticks()])
            # Rotate axis labels
            if rotate_x_labels: pyplot.setp(axis.get_xticklabels(), rotation=90)

            # Format and save
            fig.tight_layout()
            fig.savefig(filename,
                        bbox_extra_artists=[lgd],
                        bbox_inches='tight',
                        dpi=300)
            pyplot.close(fig)

    def alpha_scatter(self,
                      prefix,
                      hierarchies,
                      title=None,
                      y_lab='Isotropic B ($\AA^2$)',
                      y_lim=None,
                      v_line_hierarchy=None,
                      rotate_x_labels=True,
                      plot_mean=True,
                      add_jitter=True,
                      ):
        """Plot scatter plot over atoms (by residue) with alpha to prevent overplotting for a series of hierarchies (plotted values are the average B-factors of each residue of the hierarchies)"""

        m_h = hierarchies[0]

        # Check all hierarchies are the same
        for h in hierarchies:
            assert m_h.is_similar_hierarchy(h)

        colours = pyplot.cm.rainbow(numpy.linspace(0,1,len(hierarchies)))

        # Create a plot for each chain
        for chain_id in sorted(set([c.id for c in m_h.chains()])):

            # Create a selection for each chain and select it in each hierarchy
            sel = m_h.atom_selection_cache().selection('chain {}'.format(chain_id))
            sel_hs = [h.select(sel) for h in hierarchies]
            sel_mh = sel_hs[0]

            # Filename!
            filename = prefix + '-chain_{}.png'.format(chain_id)

            # Create x-values for each residue starting from 1
            x_vals = numpy.array(range(len(list(sel_mh.residue_groups()))))+1
            x_labels = ['']+[ShortLabeller.format(rg) for rg in sel_mh.residue_groups()]

            # Create the output figure
            fig, axis = pyplot.subplots(nrows=1, ncols=1)
            if title is not None: axis.set_title(label=str(title))
            axis.set_xlabel('Residue')
            axis.set_ylabel(y_lab)
            if y_lim:
                axis.set_ylim(y_lim)

            # Alpha value for scatter -- TODO optimise value?
            alpha = 1.0 / (1.0 + float(len(hierarchies))**0.5)

            # Mean values to be plotted as line
            mean_arr = numpy.zeros(x_vals.shape)

            # Iterative though hierarchies
            for i_h, h in enumerate(sel_hs):
                # Extract b-factors from this hierarchy
                y_vals = [numpy.array(rg.atoms().extract_b()) for rg in h.residue_groups()]
                # Append to mean arr
                y_mean = numpy.array([y.mean() for y in y_vals])
                assert len(y_mean) == len(mean_arr)
                mean_arr = mean_arr + y_mean
                # Reformat x,y for plotting
                y = numpy.concatenate(y_vals)
                x = numpy.concatenate([[x]*len(yvs) for x, yvs in zip(x_vals, y_vals)])
                if add_jitter is True:
                    jitter = numpy.random.randn(len(x)) * 0.1 / 3.0 # Normalise to be approx width of bar
                    jitter = jitter - jitter.mean()
                else:
                    jitter = 0.0
                # Plot the scatter
                axis.scatter(x=x+jitter, y=y, s=5, alpha=alpha, facecolor=colours[i_h], lw=0.0)
                # Plot line
                axis.plot(x_vals, [v.mean() for v in y_vals], '-', color=colours[i_h], alpha=0.75, zorder=5, lw=0.75)

            # Skip if no values plotted
            if (mean_arr == 0.0).all():
                continue

            if plot_mean is True:
                mean_arr = mean_arr / float(len(hierarchies))
                axis.plot(x_vals, mean_arr, 'k-', lw=1.0, zorder=6)

            if v_line_hierarchy is not None:
                v_lines = numpy.where(numpy.array([max(rg.atoms().extract_b()) for rg in v_line_hierarchy.select(sel).residue_groups()], dtype=bool))[0] + 1.5
                for val in v_lines:
                    axis.axvline(x=val, c='grey', ls='solid', label='boundaries', lw=0.5, zorder=1)
                # Add a line at zero
                h = axis.axvline(x=0.5, c='grey', ls='solid', label='boundaries', lw=0.5, zorder=1)

            # Axis ticks & labels
            x_ticks = numpy.arange(1, len(x_labels)+1, int(max(1.0, numpy.floor(len(x_labels)/20))))
            #x_ticks = numpy.unique(numpy.floor(numpy.linspace(1,len(x_labels)+1,20)))
            axis.set_xticks(x_ticks)
            axis.set_xticklabels([x_labels[int(i)] if (i<len(x_labels)) and (float(int(i))==i) else '' for i in axis.get_xticks()])
            # Rotate axis labels
            if rotate_x_labels: pyplot.setp(axis.get_xticklabels(), rotation=90)

            # Format and save
            fig.tight_layout()
            fig.savefig(filename,
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
        for chain_id in sorted(set([c.id for c in m_h.chains()])):

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
                hdl = axis.bar(left=x_vals,
                               height=y_vals[i],
                               width=1.0,
                               align='center',
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

    def tracking_plots(self, table, filename, number_to_plot=5):

        # trim the table to certain rows
        n = max(table['cycle'])
        cycles_to_plot = range(0, n, int(n/number_to_plot)+1) + [n]
        cycles_to_plot_bool = table['cycle'].isin(cycles_to_plot)
        table = table[cycles_to_plot_bool]

        fig, axes = pyplot.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
        # Create list if only one plot
        axes = numpy.array(axes).flatten()

        # Group by cycle & step to allow stacking
        grouped = table.groupby(['cycle','step'], sort=False, as_index=False)

        n_total = len(grouped)

        grouped_reduced = grouped.max()
        grouped_reduced['x'] = range(len(grouped_reduced))

        prev_x = prev_r = prev_b = None

        for n_cycle, cycle_info in grouped_reduced.groupby('cycle', sort=False):

            x_vals = cycle_info['x'].values
            r_vals = cycle_info['rmsd'].values
            b_vals = cycle_info['b_iso (total)'].values

            # Create RMSD plot
            hdl1 = axes[0].plot(
                     x_vals, r_vals,
                     'bo-', label='model',
                     lw=1, ms=max(1, min(3, 5-0.1*grouped.ngroups)),
                     )
            if prev_x is not None:
                axes[0].plot(
                         [prev_x, x_vals[0]], [prev_r, r_vals[0]],
                         'b:', label='model',
                         lw=1, ms=max(1, min(3, 5-0.1*grouped.ngroups)),
                         )
            # Create an overall B-iso TOTAL line
            hdl2 = axes[1].plot(
                     x_vals, b_vals,
                     'ko-', label='total',
                     lw=1, ms=max(1, min(3, 5-0.1*grouped.ngroups)),
                     )
            if prev_x is not None:
                axes[1].plot(
                         [prev_x, x_vals[0]], [prev_b, b_vals[0]],
                         'k:', label='total',
                         lw=1, ms=max(1, min(3, 5-0.1*grouped.ngroups)),
                         )

            prev_x = x_vals[-1]
            prev_r = r_vals[-1]
            prev_b = b_vals[-1]

        # Other plot -- done as one

        x_vals = numpy.arange(0, grouped.ngroups)
        x_max = max(x_vals)                 # remove
        x_keys = [v[0] for v in grouped]    # remove
        x_labs = [v[1] for v in x_keys]     # remove

        # Create B-iso lines for each LEVEL
        colours = self.get_level_colours()
        # Bottoms of bar where stacking occurs
        y_cuml = numpy.zeros(len(x_vals))
        # handles for legend at bottom of image
        handles = []
        for lvl in sorted(set(table['level'])):
            if lvl is None: continue
            assert isinstance(lvl, int)
            # Get values
            sel = (table['level'] == lvl)
            sel_t = table[sel]
            # The indices of the x-axis positions
            i_x = [x_keys.index(v) for v in map(tuple,sel_t[['cycle','step']].values.tolist())]
            # Extract y_vals
            y_vals = sel_t['b_iso (level)'].values

            # Plot
            hdl = axes[1].bar(left   = x_vals[i_x],
                         height = y_vals,
                         bottom = y_cuml[i_x],
                         width  = 0.8,
                         color  = colours[lvl-1],
                         edgecolor = 'k',
                         linewidth = 0.5,
                         align  = 'center',
                         label  = 'Level {}'.format(lvl))
            handles.append(hdl)
            # Add to cuml
            y_cuml[i_x] = y_cuml[i_x] + y_vals

        # Axis stuff
        ax = axes[0]
        ax.set_title('Hierarchical Model Fit')
        ax.set_xticks([])
        ax.xaxis.set_ticks_position('bottom')
        ax.set_ylabel('model fit \n$8\pi^2rms(\Delta U)$ ($\AA^2$)')
        ax.set_ylim(bottom=0.0)

        # Axis stuff
        ax = axes[1]
        ax.xaxis.set_ticks_position('bottom')
        ax.set_title('B-factors of Hierarchical Model')
        ax.set_xlabel('Optimisation Stage/Cycle')
        ax.set_ylabel('Isotropic B\n($\AA^2$)')
        ax.set_xticks(x_vals)
        ax.set_xticklabels(x_labs, rotation=90, ha='center')
        ax.tick_params('x', labelsize=max(2, min(10, 14-0.15*len(grouped))))
        ax.set_xlim(left=-0.5, right=max(x_vals)+0.5)
        ax.set_ylim(bottom=0.0)

        # Create legend for axis
        ncol = 3

        # Add legend to first graph for both lines
        lgd0 = axes[0].legend(handles=hdl1+hdl2, bbox_to_anchor=(1.02, 0.95), loc=2, borderaxespad=0.)

        # Other legends
        flip_h = []; [flip_h.extend(handles[i::ncol]) for i in range(ncol)]
        lgd1 = axes[1].legend(
                handles=flip_h, ncol=ncol,
                bbox_to_anchor=(0.5, 0.0),
                bbox_transform=fig.transFigure,
                loc=9, borderaxespad=0.,
                )

        # BOTH AXES -- Add vertical lines between macro-cycles
        start_x = x_vals[[x_keys.index(v) for v in map(tuple,table[['cycle','step']].values.tolist()) if v[1]=='start']]
        n_cycles = len(start_x)
        last_v = None
        delta = None
        for i, v in enumerate(start_x - 0.5):
            # Dashed lines to separate cycles
            if (v > 0) and (v < x_max):
                for ax in axes:
                    ax.axvline(x=v, linewidth=1, zorder=1, color='k', linestyle='--')
            # Text to label each cycle
            if last_v is not None:
                delta = v - last_v
                axes[0].text(x=last_v+delta/2.0,
                             y=0.05*axes[0].get_ylim()[0] + 0.95*axes[0].get_ylim()[1],
                             s='cycle '*(n_cycles<6) +str(cycles_to_plot[i-1]), # This is plotting the previous point so need -1
                             horizontalalignment='center',
                             verticalalignment='top',
                            )
            last_v = v
        # Plot the last point (or do nothing for 1 cycle)
        if delta is not None:
            axes[0].text(x=min(v+delta/2.0, axes[0].get_xlim()[1]),
                         y=0.05*axes[0].get_ylim()[0] + 0.95*axes[0].get_ylim()[1],
                         s='cycle '*(n_cycles<6) +str(cycles_to_plot[i]), # This does not need a -1
                         horizontalalignment='center',
                         verticalalignment='top',
                        )

        fig.tight_layout()
        fig.savefig(filename,
                    bbox_extra_artists=[lgd0,lgd1],
                    bbox_inches='tight',
                    dpi=200)
        pyplot.close(fig)

    def convergence_plots(self, table, filename):

        fig, axes = pyplot.subplots(nrows=1, ncols=2, sharex=True, sharey=False)
        axes = numpy.array(axes).flatten()

        # Extract rows with non-zero values
        #table = table[table['b_iso (level)']!=0.0]
        # Extract only inter-level optimisation values (last step of each cycle)
        table = table[(table['step']=='inter-level')]
        # Labels for each of the series to plot
        m_cyc = 0 if (len(table) == 0) else min(table['cycle'])
        labels = table[table['cycle']==m_cyc]['level'].values

        # Colours for each level
        colours = self.get_level_colours()

        ########################
        # FIRST AXIS
        ########################
        ax = axes[0]
        handles = []
        # Extract common list of x-values
        x_keys = sorted(set(table['cycle'].values))
        x_vals = numpy.array(x_keys)
        # Select labels to plot (up to maximum)
        max_labels = 1 + 10
        plot_every = max(1, 1+((len(x_vals)-1)//max_labels))
        x_labs = [x_vals[i] if (i%plot_every)==0 else '' for i in xrange(len(x_vals))]
        # Cumulative y-values for stacking
        y_cuml = numpy.zeros(len(x_vals))
        # Plot same values as stacked bars
        for l in labels:
            assert isinstance(l, int)
            # Extract relevant rows from table
            l_table = table[table['level']==l]
            # Extract plot vals
            i_x = [x_keys.index(v) for v in l_table['cycle'].values]
            y_vals = l_table['b_iso (level)'].values
            # Plot stacked bar
            hdl = ax.bar(left   = x_vals[i_x],
                         height = y_vals,
                         bottom = y_cuml[i_x],
                         width  = 0.8,
                         color  = colours[l-1],
                         edgecolor = 'k',
                         linewidth = 0.5,
                         align  = 'center',
                         label  = 'Level {}'.format(l))
            handles.append(hdl)
            # Add to cuml
            y_cuml[i_x] = y_cuml[i_x] + y_vals

        # Create legend for axis
        ncol = 3
        flip_h = []; [flip_h.extend(handles[i::ncol]) for i in range(ncol)]
        lgd0 = ax.legend(
                handles=flip_h, ncol=ncol,
                bbox_to_anchor=(0.5, 0.0),
                bbox_transform=fig.transFigure,
                loc=9, borderaxespad=0.,
                )

        ax.xaxis.set_ticks_position('bottom')
        ax.set_xticks(x_vals)
        ax.set_xticklabels(x_labs)
        ax.tick_params('x', labelsize=max(6, min(10, 14-0.15*len(x_keys))))
        ax.set_xlabel('Optimisation Cycle')
        ax.set_ylabel('Average B-factor of Level')
        ax.set_ylim(bottom=0.0)

        ########################
        # SECOND AXIS
        ########################
        ax = axes[1]
        handles = []
        for l in labels:
            assert isinstance(l, int)
            # Extract relevant rows from table
            l_table = table[table['level']==l]
            # Extract plot vals
            x_vals = l_table['cycle'].values
            y_vals = l_table['b_iso (level)'].values

            hd_ = ax.plot(x_vals, y_vals,
                    'ko-', lw=2, ms=max(3, min(5, 7-0.1*len(x_vals))))
            hdl = ax.plot(x_vals, y_vals,
                    'o-', lw=1, ms=max(1, min(3, 5-0.1*len(x_vals))),
                    color=colours[l-1], label='Level {}'.format(l))
            handles.extend(hdl)

        # Axis stuff
        ax.xaxis.set_ticks_position('bottom')
        ax.tick_params('x', labelsize=max(6, min(10, 14-0.15*len(x_keys))))
        ax.set_xlabel('Optimisation Cycle')
        ax.set_ylabel('Average B-factor of Level')
        ax.set_ylim(bottom=0.0)

        # Create legend for axis
        #lgd1 = ax.legend(handles=handles, ncol=3, bbox_to_anchor=(0.00, -0.15), loc=9, borderaxespad=0.)

        t = fig.suptitle('Convergence of level b-factors',
                y = 1.00, verticalalignment='bottom')

        fig.tight_layout()
        fig.savefig(filename,
                    bbox_extra_artists=[t, lgd0],
                    bbox_inches='tight',
                    dpi=200)
        pyplot.close(fig)

class MultiDatasetHierarchicalUijFitter(object):

    def __init__(self,
                 observed_uij,
                 observed_xyz,
                 level_array,
                 level_labels=None,
                 uij_weights=None,
                 isotropic_mask=None,
                 params=None,
                 verbose=False,
                 log=None):

        #if log is None: log = Log()
        self.log = log
        self.params = params
        self.verbose = verbose

        self.warnings = []

        assert observed_uij.shape[1]  == level_array.shape[1]
        assert observed_uij.shape[:-1] == observed_xyz.shape[:-1]
        assert observed_uij.shape[-1]  == 6
        assert observed_xyz.shape[-1]  == 3

        # Store observed values (needed for later)
        self.observed_uij = observed_uij
        self.observed_xyz = observed_xyz

        self.n_datasets = observed_uij.shape[0]
        self.n_atoms = observed_uij.shape[1]

        # Masks to exclude datasets or atoms from optimisation
        self.dataset_mask = numpy.ones(self.n_datasets, dtype=bool)
        self.atomic_mask  = numpy.ones(self.n_atoms, dtype=bool)

        # Weights for each uij
        self.uij_weights = uij_weights

        self.set_isotropic_mask(isotropic_mask=isotropic_mask)

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
            self.levels.append(MultiDatasetTLSGroupLevel(observed_uij = observed_uij,
                                                         observed_xyz = observed_xyz,
                                                         tls_origins      = self.tls_origins,
                                                         tls_origins_hash = self.tls_origins_hash,
                                                         group_idxs = group_idxs,
                                                         n_tls = self.params.n_tls_modes_per_tls_group,
                                                         index = idx+1,
                                                         label = lab,
                                                         weights = self.uij_weights,
                                                         isotropic_mask = self.isotropic_mask,
                                                         params = self.params.optimisation,
                                                         verbose = verbose,
                                                         log = self.log))

        # One object to fit all the Uij residuals
        self.residual = MultiDatasetResidualLevel(observed_uij = observed_uij,
                                                  index = len(self.levels)+1,
                                                  label = 'residual',
                                                  weights = self.uij_weights,
                                                  params = self.params.optimisation,
                                                  verbose = verbose,
                                                  log = self.log)

        assert len(self.level_labels) == len(self.levels)

        self.apply_masks()
        self.make_hierarchy_as_tree()

        self.plot = MultiDatasetUijPlots(n_levels=len(self.levels)+1)

    def __iter__(self):
        for i_level, level in enumerate(self.levels):
            yield (i_level+1, self.level_labels[i_level], level)

    def _calculate_target_uij(self, fitted_uij_by_level, i_level):
        arr = fitted_uij_by_level.copy()
        # Zero the current level (as should not be included in target function!)
        arr[i_level] = 0.0
        # Sum over all other levels to find total currently modelled
        arr_sum = numpy.sum(arr, axis=0)
        # If input data is isotropic, make output isotropic for target function
        if self.disorder_model == "isotropic":
            iso_arr_sum = arr_sum[:,:,0:3].mean(axis=2)
            arr_sum = numpy.zeros_like(arr_sum)
            for i in range(3):
                arr_sum[:,:,i] = iso_arr_sum
        elif self.disorder_model == "mixed":
            # Calculate the isotropic values for the isotropic atoms
            iso_arr_sum = arr_sum[:,self.isotropic_mask,0:3].mean(axis=2)
            # Zero out the values in-place and replace with isotropic values
            arr_sum[:,self.isotropic_mask,:] = 0.0
            for i in range(3):
                arr_sum[:,self.isotropic_mask,i] = iso_arr_sum
        assert arr_sum.shape == self.observed_uij.shape
        # Subtract currently modelled from input data
        return self.observed_uij - arr_sum

    def disable_fitting_penalties(self):
        """Remove flag for calculating fitting penalties during optimisation (unconstrained optimisation)"""
        for l in self.levels:
            for g_n, g_s, g_f in l:
                g_f.use_fitting_penalties = False

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
                    overlap = set(g_this[i_l_test]).intersection(set(g_other[i_l_test])).difference({-1})
                    if overlap:
                        self.log('      ...there is a group on level {} whose atoms are split between this group and another'.format(i_l_test+1))
                        self.log('      ...(groups: {})'.format(', '.join(sorted(overlap))))
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

    def make_hierarchy_as_tree(self):
        """Identify the tree of groups that form each atom"""

        # Dictionary pointing to "owned" groups in the next level(s)
        # Goes from lower levels to higher levels: graph[chain_level] -> groups on ss_level
        graph = {}

        residual_level = numpy.arange(1, self.level_array.shape[1]+1)

        # Search downwards
        for i_above in range(1,self.level_array.shape[0] + 1): # +1 causes inclusion of the residual level
            if i_above == len(self.level_array):
                l_label = 'X'
                l_above = residual_level
            else:
                l_label = i_above+1
                l_above = self.level_array[i_above]

            # Populate
            for v_above in sorted(numpy.unique(l_above)):
                # Skip in not in a group
                if (v_above == -1): continue
                # Select atoms in this group
                v_sel = (l_above == v_above)

                # Find the nearest group "below" (may not be on the adjacent level)
                i_below = i_above-1
                while i_below >= 0:
                    l_below = self.level_array[i_below]
                    if (l_below[v_sel] != -1).sum() > 0:
                        break
                    i_below -= 1
                if i_below < 0:
                    raise Exception('Broken! Couldn\'t find a group in a lower level')

                # Create dictionary for the lower level, pointing to higher levels
                below_dict = graph.setdefault(i_below+1, {})
                # Iterate through the groups in the lower level
                for v_below in sorted(numpy.unique(l_below[v_sel])):
                    if (v_below == -1): continue
                    # Initialise a dictionary for this group
                    below_group_dict = below_dict.setdefault(v_below, {})
                    # Create list of groups in this level that correspond to higher levels
                    above_list = below_group_dict.setdefault(l_label, [])
                    # Add the group
                    above_list.append(v_above)

        if self.verbose or True:
            self.log.subheading('Tree summary of hierarchy:')
            for i in graph.keys():
                self.log.bar()
                self.log('Ownership of level {} (Related groups that will be co-optimised)'.format(i))
                self.log.bar()
                for g, vals in graph[i].items():
                    groups = ['(Level {}, Groups {})'.format(l, ', '.join(map(str,vs))) for l, vs in vals.items()]
                    self.log('Group {} -> {}'.format(g, ', '.join(groups)))
                    #for l, vs in vals.items():
                    #    self.log('\t-> Level {}, Groups {}'.format(l, ', '.join(map(str,vs))))

        self.group_tree = graph

        return graph

    def n_levels(self):
        return len(self.levels)

    def n_params(self, non_zero=False):
        return sum([l.n_params(non_zero=non_zero) for l in self.levels+[self.residual]])

    def n_params_per_atom(self, non_zero=False):
        return self.n_params(non_zero=non_zero) / (self.observed_uij.shape[0] * self.observed_uij.shape[1])

    def n_input_params_per_atom(self):
        return self._input_params_per_atom

    def n_input_values(self):
        if self.n_input_params_per_atom() is None:
            return 'unknown'
        return int(self.observed_uij.shape[0] * self.observed_uij.shape[1] * self.n_input_params_per_atom())

    def optimise_level_amplitudes(self, n_cpus=1, max_recursions=None, include_residual=True, last_cycle=False):
        """Optimise amplitudes for pairs of adjacent levels"""

        self.log.bar(True, False)
        self.log('Running inter-level amplitude optimisation')
        self.log.bar(False, True)

        self.log('CPUs: {}'.format(n_cpus))
        self.log('Max recursions (number of levels to be co-optimised): {}'.format(max_recursions))

        # Create optimisation object -- edits levels in-place so shouldn't need to do anything to unpack results
        ao = InterLevelAmplitudeOptimiser(
                target_uij = self.observed_uij,
                levels     = self.levels,
                residual   = self.residual,
                group_tree = self.group_tree,
                convergence_tolerance = self.params.optimisation.gradient_convergence,
                amplitude_sum_weight = 0.0 if last_cycle is True else self.params.optimisation.amplitude_sum_weight,
                weights    = self.uij_weights,
                isotropic_mask = self.isotropic_mask,
                params     = self.params,
                verbose    = self.verbose,
                log        = self.log)
        ao.set_residual_mask(flex.bool(self.dataset_mask))
        ao.optimise(n_cpus=n_cpus, max_recursions=max_recursions)
        ao.apply_multipliers(self.levels, self.residual)
        self.log(ao.summary(show=True))

        # Update output
        fitted_uij = numpy.zeros((len(self.levels)+1,)+self.observed_uij.shape)
        for i, l in enumerate(self.levels):
            fitted_uij[i] = l.extract()
        fitted_uij[-1] = self.residual.extract()
        assert fitted_uij.shape == (len(self.levels)+1,) + self.observed_uij.shape
        return fitted_uij

    def parameter_ratio_gain(self, non_zero=False):
        """Total number of input values divided by number of model parameters (used or total)"""
        if self.n_input_params_per_atom() is None:
            return 'unknown'
        return self.n_params(non_zero=non_zero) / self.n_input_values()

    def set_isotropic_mask(self, isotropic_mask):
        if isotropic_mask is None:
            disorder_model = 'anisotropic'
            self._input_params_per_atom = 6
        elif isotropic_mask.sum() == isotropic_mask.size:
            disorder_model = 'isotropic'
            self._input_params_per_atom = 1
        else:
            assert isotropic_mask.sum() > 0, 'isotropic mask has been provided that selects no atoms'
            disorder_model = 'mixed'
            self._input_params_per_atom = 6.0 - 5.0*float(isotropic_mask.sum())/float(isotropic_mask.size)

        self.disorder_model = disorder_model
        self.isotropic_mask = isotropic_mask

    def set_optimisation_datasets(self, dataset_indices):
        self.dataset_mask[:] = False
        self.dataset_mask[dataset_indices] = True

    def set_tracking(self, table, csv_path, trk_path, cvg_path):
        self.tracking_data = table
        self.tracking_csv = csv_path
        self.tracking_png = trk_path
        self.convergence_png = cvg_path

    def update_tracking(self, uij_lvl, step, i_cycle, i_level=None, write_graphs=False):
        """Update the tracking table"""

        self.log.subheading('Updating tracking...')

        # Extract uijs for all of the levels for all datasets
        uij_tot = uij_lvl.sum(axis=0)
        # Calculate the rms between fitted and input
        rmsd = EIGHT_PI_SQ*rms(self.observed_uij-uij_tot, axis=None)

        # Average over all datasets
        uij_tot = uij_tot.mean(axis=0)
        uij_lvl = uij_lvl.mean(axis=1)

        assert uij_lvl.shape == (len(self.levels)+1,) + self.observed_uij.shape[1:]
        assert uij_tot.shape == self.observed_uij.shape[1:]

        if not isinstance(i_level, list):
            i_level = [i_level]

        # Iterate through levels to be dumped into table
        for i_l in i_level:

            # Extract the Uij for the selected level(s)
            if isinstance(i_l, int):
                uij_sel = uij_lvl[i_l]
                level = i_l+1
            else:
                uij_sel = None
                level = None

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
            cycle_lab = i_cycle

            # Add to tracking table
            self.tracking_data.loc[len(self.tracking_data.index)] = {'cycle' : cycle_lab,
                                                                     'step'  : step,
                                                                     'level' : level,
                                                                     'rmsd'  : round(rmsd,3),
                                                                     'u_iso (level)' : round(u_iso_sel,3),
                                                                     'b_iso (level)' : round(b_iso_sel,3),
                                                                     'b_min (level)' : round(b_min_sel,3),
                                                                     'b_max (level)' : round(b_max_sel,3),
                                                                     'u_iso (total)' : round(u_iso_tot,3),
                                                                     'b_iso (total)' : round(b_iso_tot,3)}

        self.log(self.tracking_data.loc[len(self.tracking_data)-len(i_level):].to_string())

        if write_graphs:
            # Dump to csv
            self.tracking_data.to_csv(self.tracking_csv)
            # Make plots
            self.plot.tracking_plots(
                    table = self.tracking_data,
                    filename = self.tracking_png
                    )
            self.plot.convergence_plots(
                    table = self.tracking_data,
                    filename = self.convergence_png
                    )

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

    def fit(self, n_cpus=1, n_macro_cycles=1, n_micro_cycles=1):
        """Run macro-cycles of parameter optimisation"""

        self.log('Fitting TLS models for {} levels (+ residual)'.format(len(self.levels)))

        # Cumulative fitted Uij from the different levels (+1 level for the residuals!)
        # Indexing: [level, dataset, atom, uij]
        fitted_uij_by_level = numpy.zeros((len(self.levels)+1,)+self.observed_uij.shape)

        # Run macro cycles of optimisation
        for i_macro in xrange(n_macro_cycles+1):
            self.log.heading('Macrocycle {} of {}'.format(i_macro, n_macro_cycles), spacer=True)

            # Things that are different for the first cycle -- optimise against one dataset only
            if (i_macro == 0):
                # Set dataset mask to the first optimisation dataset to give some reasonable starting values
                optimisation_datasets = numpy.where(self.dataset_mask)[0].tolist()
                self.set_optimisation_datasets([optimisation_datasets[0]])

            # Re-initialise tracking at the beginning of each subsequent cycle
            self.update_tracking(uij_lvl=fitted_uij_by_level, step='start', i_cycle=i_macro)

            # Ensure the masks are up-to-date
            self.apply_masks()

            # Iterate through the TLS levels of the fitting
            for i_level, fitter in enumerate(self.levels):
                self.log.subheading('Macrocycle {} of {}: '.format(i_macro, n_macro_cycles)+'Fitting TLS Groups (level {} - {})'.format(fitter.index, fitter.label))
                # Update the target uij by subtracting contributions from other levels
                self.log('Updating target Uijs for optimisation')
                fitter.set_target_uij(target_uij=self._calculate_target_uij(fitted_uij_by_level=fitted_uij_by_level, i_level=i_level))
                # Optimise
                fitted_uij_by_level[i_level] = fitter.run(n_cpus=n_cpus, n_cycles=(max(5,n_micro_cycles) if i_macro==0 else n_micro_cycles))
                # Optimise the amplitudes between levels
                if i_macro > 0:
                    self.log.subheading('Macrocycle {} of {}: '.format(i_macro, n_macro_cycles)+'Optimising inter-level amplitudes')
                    fitted_uij_by_level = self.optimise_level_amplitudes(n_cpus=1, max_recursions=1, include_residual=True, last_cycle=(i_macro==n_macro_cycles))
                    fitted_uij_by_level = self.optimise_level_amplitudes(n_cpus=1, max_recursions=None, include_residual=True, last_cycle=(i_macro==n_macro_cycles))
                # Update tracking
                self.update_tracking(
                        uij_lvl=fitted_uij_by_level,
                        step='level {}'.format(i_level+1),
                        i_cycle=i_macro,
                        i_level=i_level,
                        write_graphs=(i_macro==0),
                        )

            if self.params.fit_residual:
                # Fit the residuals
                self.log.subheading('Macrocycle {} of {}: '.format(i_macro, n_macro_cycles)+'Optimising residual atomic Uijs')
                # Update the target uij by subtracting contributions from other levels
                self.residual.set_target_uij(target_uij=self._calculate_target_uij(fitted_uij_by_level=fitted_uij_by_level, i_level=-1))
                # Update fitters and optimise -- always run two cycles of this
                fitted_uij_by_level[-1] = self.residual.run(n_cpus=n_cpus, n_cycles=1)
                # Optimise the amplitudes between levels
                if i_macro > 0:
                    self.log.subheading('Macrocycle {} of {}: '.format(i_macro, n_macro_cycles)+'Optimising inter-level amplitudes')
                    fitted_uij_by_level = self.optimise_level_amplitudes(n_cpus=1, max_recursions=1, include_residual=True, last_cycle=(i_macro==n_macro_cycles))
                    fitted_uij_by_level = self.optimise_level_amplitudes(n_cpus=1, max_recursions=None, include_residual=True, last_cycle=(i_macro==n_macro_cycles))
                # Update tracking
                self.update_tracking(
                        uij_lvl=fitted_uij_by_level,
                        step='residual',
                        i_cycle=i_macro,
                        i_level=len(fitted_uij_by_level)-1,
                        write_graphs=(i_macro==0),
                        )

            # Update tracking
            self.update_tracking(
                    uij_lvl=fitted_uij_by_level,
                    step='inter-level',
                    i_cycle=i_macro,
                    i_level=range(len(fitted_uij_by_level)),
                    write_graphs=True,
                    )

            if (i_macro == 0):
                # Restore the original optimisation datasets
                self.set_optimisation_datasets(optimisation_datasets)
                # Turn off fitting penalties for future optimisations
                self.disable_fitting_penalties()

        # Check validity of output model and generate warnings...
        self.validate_fitted_uijs()

        return self.extract()

    def validate_fitted_uijs(self):
        """Check for negative Uij eigenvalues of output models, etc"""

        self.log.subheading('Validating output Uijs')

        num_to_print = 3
        uij_fmt_str = ', '.join(6*['{{:.{:d}f}}'.format(int(-1*numpy.log10(self.params.precision.uij_tolerance)+1.0))])
        eig_fmt_str = ', '.join(3*['{{:.{:d}f}}'.format(int(-1*numpy.log10(self.params.precision.uij_tolerance)+1.0))])

        warnings = []
        for i_level, level in enumerate(self.levels):
            fitted_uij = level.extract()
            for i_dst, dataset_uij in enumerate(fitted_uij):
                uij_valid = uij_positive_are_semi_definite(uij=dataset_uij, tol=0.01) #self.params.precision.uij_tolerance)
                if uij_valid.sum() > 0.0:
                    fmt_list = [(i+1, level.get_atom_group(index=i).label, i_level+1, uij_fmt_str.format(*tuple(dataset_uij[i])), eig_fmt_str.format(*tuple(sym_mat3_eigenvalues(dataset_uij[i])))) for i in numpy.where(uij_valid)[0]]
                    err_msgs = ['Atom {:>5d} (group {}, level {}): \n\t\tUij -> ({})\n\t\tEigenvalues -> ({})'.format(*d) for d in fmt_list]
                    if (len(err_msgs) > num_to_print) and (not self.verbose):
                        n_hidden = len(err_msgs) - num_to_print
                        err_msgs = err_msgs[:num_to_print]
                        err_msgs.append('[...] ({} more similar warning{} not shown)'.format(n_hidden, 's' if n_hidden>1 else ''))
                    w = 'Level {}: Uijs for dataset {} are not positive-semi-definite ({} atoms)\n\t{}'.format(i_level+1, i_dst, uij_valid.sum(), '\n\t'.join(err_msgs))
                    warnings.append(w)
        if warnings:
            self.log('{} Warnings!'.format(len(warnings)))
            for w in warnings: self.log(w)
        else:
            self.log('No warnings during validation.')
        self.warnings.extend(warnings)

    def extract(self, average_datasets=False):
        """Extract the fitted Uij for all structures or the average structure (averaged TLS)"""
        return self.extract_tls(sum_levels=True, sum_modes=True, average_datasets=average_datasets) + self.residual.extract()

    def extract_tls(self, sum_levels=False, sum_modes=True, average_datasets=False):
        """Extract the fitted TLS Uij for all structure or the average structure (averaged TLS)"""
        uij = numpy.array([l.extract(sum_modes=sum_modes, average_datasets=average_datasets) for l in self.levels])
        if sum_levels is True:
            uij = uij.sum(axis=0)
        return uij

    def summary(self, show=False):
        s = ''
        s += '\n> Input Data:'
        s += '\nNumber of datasets: {}'.format(self.n_datasets)
        s += '\nNumber of atoms (in each structure): {}'.format(len(self.atomic_mask))
        s += '\nInput Disorder Type: {}'.format(self.disorder_model)
        if self.disorder_model == 'mixed':
            n_iso = self.isotropic_mask.sum()
            n_ani = self.isotropic_mask.size - n_iso
            s += '\n\t{:d} anisotropic atoms ({:5.2%})'.format(n_ani, n_ani/float(self.isotropic_mask.size))
            s += '\n\t{:d} isotropic atoms   ({:5.2%})'.format(n_iso, n_iso/float(self.isotropic_mask.size))
        s += '\nNumber of uij values per atom (average): {:.2f}'.format(self.n_input_params_per_atom())
        s += '\nNumber of uij values (total): {}'.format(self.n_input_values())
        s += '\n'
        s += '\n> Hierarchical Model:'
        s += '\nNumber of levels (incl. residual): {}'.format(self.n_levels()+1)
        s += '\nNumber of tls groups per level:'
        for l in self.levels:
            s += '\n\tLevel {} ({}): {}'.format(l.index, l.label, str(l.n_groups()))
        s += '\nNumber of TLS modes per group: {}'.format(self.params.n_tls_modes_per_tls_group)
        s += '\nNumber of model parameters per group: {} x (20 + {}) = {}'.format(
                self.params.n_tls_modes_per_tls_group,
                self.n_datasets,
                self.params.n_tls_modes_per_tls_group*(20+self.n_datasets),
                )
        s += '\nNumber of atoms in residual level: {}'.format(self.residual.n_atoms())
        s += '\nNumber of model parameters (total): {}'.format(self.n_params(non_zero=False))
        #s += '\nNumber of model parameters per level:'
        #for l in self.levels+[self.residual]:
        #    s += '\n\tLevel {}: {}'.format(l.index, str(l.n_params(non_zero=False)))
        s += '\n'
        s += '\n> Optimisation:'
        s += '\nDatasets used for TLS and residual optimisation: {} of {}'.format(sum(self.dataset_mask), len(self.dataset_mask))
        s += '\nAtoms used for TLS optimisation: {} of {}'.format(sum(self.atomic_mask), len(self.atomic_mask))
        s += '\nDataset Weighting: {}'.format(self.params.optimisation.dataset_weights)
        s += '\nAtomic Weighting: {}'.format(self.params.optimisation.atom_weights)
        s += '\nAtomic weights renormalised by dataset: {}'.format("Yes" if self.params.optimisation.renormalise_atom_weights_by_dataset else "No")
        s += '\n'
        s += '\n> Parameterisation:'
        s += '\nNumber of parameters per atom:'
        s += '\n... in input data (average over all atoms): {}'.format(format_no_fail(':2.2f', self.n_input_params_per_atom()))
        s += '\n... in fitted model: {:5.2f}'.format(self.n_params_per_atom(non_zero=False))
        s += '\n... in fitted model (non-zero): {:5.2f}'.format(self.n_params_per_atom(non_zero=True))
        s += '\nNon-zero model parameters / number of input data: {}'.format(format_no_fail(':5.2%', self.parameter_ratio_gain(non_zero=True)))
        s += '\n'
        if show: self.log(s)
        return s


class _MultiDatasetUijLevel(object):

    _uij_shape = None
    _chunksize = 1

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
            self.log('')
            sigint = signal.signal(signal.SIGINT, signal.SIG_IGN) # set the signal handler to ignore
            pool = NonDaemonicPool(n_cpus)
            signal.signal(signal.SIGINT, sigint) # Replace the original signal handler
            finished_jobs = []
            chunksize = min(self._chunksize, 1+int(float(len(jobs)-1)/float(n_cpus)))
            try:
                pbar = tqdm.tqdm(total=len(jobs), ncols=100)
                for result in pool.imap_unordered(func=_wrapper_fit, iterable=jobs, chunksize=chunksize):
                    pbar.update(1)
                    finished_jobs.append(result)
                pbar.close()
            except KeyboardInterrupt:
                pbar.close()
                pool.terminate()
                pool.join()
                raise
            # Close pool
            pool.close()
            pool.join()
            # Resort the output so results presented in a readable fashion
            def get_label(x):
                if isinstance(x, str): return x
                return x.label
            finished_jobs = sorted(finished_jobs, key=lambda x: get_label(x))
            #finished_jobs = workers.map(func=_wrapper_fit, iterable=jobs, chunksize=self._chunksize)
        else:
            # Run jobs in serial
            self.log('Running {} job(s) [with {} cpu(s)]'.format(len(jobs), n_cpus_per_job))
            self.log('')
            finished_jobs = []
            pbar = tqdm.tqdm(total=len(jobs), ncols=100)
            for result in (_wrapper_fit(j) for j in jobs):
                pbar.update(1)
                finished_jobs.append(result)
            pbar.close()
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
                self.log('Log file written to {}'.format(str(fitter.log.log_file())))
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

class MultiDatasetTLSGroupLevel(_MultiDatasetUijLevel):

    _chunksize = 20

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
                 isotropic_mask=None,
                 params=None,
                 verbose=False,
                 log=None):

        #if log is None: log = Log()
        self.log = log

        assert observed_uij.shape[:2] == observed_xyz.shape[:2]
        assert observed_xyz.shape[2]  == 3
        assert observed_uij.shape[1]  == len(group_idxs)
        assert observed_uij.shape[2]  == 6
        assert tls_origins.shape[1] == observed_uij.shape[0]
        assert tls_origins.shape[2] == 3
        assert tls_origins_hash.shape[0] == observed_xyz.shape[1]

        if weights is not None:
            assert weights.shape == observed_uij.shape[:-1]

        self.index = index
        self.label = label if label else 'Level {}'.format(index)

        self.group_idxs = group_idxs

        self._n_tls = n_tls
        self._n_groups = sum(numpy.unique(self.group_idxs)!=-1)
        self._n_obj = self._n_groups

        self._uij_shape = observed_uij.shape

        self.fitters = {}
        for i, sel, f in self:
            assert f is None
            # Create a separate log for each fitter, that gets stored in memory rather than being written immediately
            ls = LogStream(path=os.path.splitext(self.log.log_file().path)[0]+'-level{:04d}-group{:06d}.log'.format(self.index, i))
            log = Log(stdout=True).add_output(ls).toggle(status=(self._n_groups==1))
            # Decide which tls_origins to use fo this
            i_origin = numpy.unique(tls_origins_hash[sel])
            assert len(i_origin) == 1
            observed_com = tls_origins[i_origin[0]]
            assert observed_com.shape == (observed_xyz.shape[0], 3)
            # Select isotropic_mask
            if isotropic_mask is not None:
                sel_iso_mask = isotropic_mask[sel]
            else:
                sel_iso_mask = None
            # Create fitter object
            self.fitters[i] = MultiDatasetTLSFitter(target_uij = observed_uij[:,sel],
                                                          atomic_xyz = observed_xyz[:,sel],
                                                          atomic_com = observed_com,
                                                          n_tls = n_tls,
                                                          label = '{:4d} of {:4d}'.format(i, self._n_groups),
                                                          weights = weights[:,sel],
                                                          isotropic_mask = sel_iso_mask,
                                                          params = params,
                                                          verbose = verbose,
                                                          log = log)

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

    def extract(self, sum_modes=True, average_datasets=False):
        shape = self._uij_shape
        dst_axis = 0
        if sum_modes is False:
            shape = (self._n_tls,) + shape
            dst_axis += 1
        fitted_uij = numpy.zeros(shape)
        for i, sel, fitter in self:
            f_uij = fitter.extract(sum_modes=sum_modes)
            if sum_modes is False:
                fitted_uij[:,:,sel] = f_uij
            else:
                fitted_uij[:,sel] = f_uij
        if average_datasets is True:
            fitted_uij = fitted_uij.mean(axis=dst_axis)
        return fitted_uij

    def extract_core(self):
        shape = (self._n_tls,) + self._uij_shape
        core_uijs = numpy.zeros(shape)
        for i, sel, fitter in self:
            f_uij = fitter.extract_core()
            core_uijs[:,:,sel] = f_uij
        return core_uijs

    def summary(self, show=True):
        num_model_params = self.n_params(non_zero=False)
        num_used_params = self.n_params(non_zero=True)
        n_dst, n_atm, _ = self._uij_shape
        s = '> Level {} - {}'.format(self.index, self.label)
        s += '\n\tNumber of Groups: {}'.format(self.n_groups())
        s += '\n\tModel parameters (total): {}'.format(num_model_params)
        s += '\n\tModel parameters (used): {}'.format(num_used_params)
        s += '\n\tModel parameters per atom (total): {:5.3f}'.format(num_model_params/(n_dst*n_atm))
        s += '\n\tModel parameters per atom (used):  {:5.3f}'.format(num_used_params/(n_dst*n_atm))
        s += '\n'
        if show: self.log(s)
        return s

    def zero_amplitudes(self, models):
        # Iterate through each group/partition
        for i, sel, fitter in self:
            fitter.parameters().zero_amplitudes(selection=models)

    #def zero_matrices(self, models, components=None):
    #    # Iterate through each group/partition
    #    for i, sel, fitter in self:
    #        fitter.parameters().zero_matrices(models=models, components=components)

class MultiDatasetResidualLevel(_MultiDatasetUijLevel):

    _chunksize = 250

    def __init__(self,
                 observed_uij,
                 index=-1,
                 label=None,
                 weights=None,
                 params=None,
                 verbose=False,
                 log=None):

        #if log is None: log = Log()
        self.log = log

        self.index = index
        self.label = label if label else 'Level {} (residual)'.format(index)

        self._uij_shape = observed_uij.shape

        self._n_atm = self._uij_shape[1]
        self._n_obj = self._n_atm

        # Create one log for all fitters
        ls = LogStream(path=os.path.splitext(self.log.log_file().path)[0]+'-level-residual.log')
        log = Log(stdout=True).add_output(ls).toggle(0)
        # Create a fitter for each atom
        self.fitters = {}
        for i, sel, f in self:
            assert f is None
            # Create fitter object
            self.fitters[i] = MultiDatasetUijFitter(target_uij = observed_uij[:,sel],
                                                           label = 'atom {:5d} of {:5d}'.format(i,self._n_atm),
                                                           weights = weights[:,sel],
                                                           params = params,
                                                           verbose = verbose,
                                                           log = log)

    def __iter__(self):
        for i in xrange(self._n_atm):
            yield (i+1, i, self.fitters.get(i+1, None))

    def n_atoms(self):
        return self._n_atm

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

############################################################################

def build_levels(model, params, log=None):
    """Build the levels for the hierarchical fitting"""

    log.subheading('Building hierarchy selections')
    levels = []; labels=[];

    filter_h = model.hierarchy
    log('Input hierarchy contains {} atoms'.format(filter_h.atoms().size()))

    # Filter the hierarchy by the overall selection
    if params.levels.overall_selection:
        cache = filter_h.atom_selection_cache()
        choice = cache.selection(params.levels.overall_selection)
        filter_h = filter_h.select(choice, copy_atoms=True)
        log('Overall selection selects {} atoms'.format(filter_h.atoms().size()))

    # Whether cbeta is in backbone
    cbeta_flag = params.levels.cbeta_in_backbone
    backbone_atoms = ['C','CA','N','O']+(['CB']*cbeta_flag)
    backbone_atoms_sel = '({})'.format(' or '.join(['name {}'.format(a) for a in backbone_atoms]))
    back_sel = ' and '+backbone_atoms_sel
    side_sel = ' and not '+backbone_atoms_sel

    log.bar(True, False)
    log('Creating automatic levels:')

    if 'chain' in params.levels.auto_levels:
        log('Level {}: Creating level with groups for each chain'.format(len(levels)+1))
        groups = [PhenixSelection.format(c) for c in filter_h.chains()]
        levels.append(sorted(set(groups))) # Chains can be present multiple times
        labels.append('chain')
    if 'auto_group' in params.levels.auto_levels:
        log('Level {}: Creating level with groups determined by phenix.find_tls_groups'.format(len(levels)+1))
        groups = [s.strip('"') for s in phenix_find_tls_groups(model.filename)]
        levels.append([g for g in groups if sum(cache.selection(g)>0)])
        labels.append('groups')
    if ('secondary_structure' in params.levels.auto_levels) or ('ss' in params.levels.auto_levels):
        log('Level {}: Creating level with groups based on secondary structure'.format(len(levels)+1))
        groups = [s.strip('"') for s in default_secondary_structure_selections_filled(hierarchy=filter_h, verbose=params.settings.verbose)]#model.hierarchy)]
        levels.append([g for g in groups if sum(cache.selection(g))>0])
        labels.append('sec. struct.')
    if 'residue' in params.levels.auto_levels:
        log('Level {}: Creating level with groups for each residue'.format(len(levels)+1))
        levels.append([PhenixSelection.format(r) for r in filter_h.residue_groups()])
        labels.append('residue')
    if 'backbone' in params.levels.auto_levels:
        log('Level {}: Creating level with groups for each residue backbone'.format(len(levels)+1))
        gps = backbone(filter_h, cbeta=cbeta_flag).atom_groups()
        levels.append([PhenixSelection.format(r)+back_sel for r in gps if (r.resname not in ['ALA','GLY','PRO'])])
        labels.append('backbone')
    if 'sidechain' in params.levels.auto_levels:
        log('Level {}: Creating level with groups for each residue sidechain'.format(len(levels)+1))
        gps = sidechains(filter_h, cbeta=(not cbeta_flag)).atom_groups()
        levels.append([PhenixSelection.format(r)+side_sel for r in gps if (r.resname not in ['ALA','GLY','PRO'])])
        labels.append('sidechain')
    if 'atom' in params.levels.auto_levels:
        log('Level {}: Creating level with groups for each atom'.format(len(levels)+1))
        levels.append([PhenixSelection.format(a) for a in filter_h.atoms()])
        labels.append('atom')
    log.bar()

    # Insert custom levels
    if params.levels.custom_level:
        # Print auto levels
        log('> {} automatic levels created:'.format(len(levels)))
        for i_l, level in enumerate(levels):
            log('\tLevel {} ({})'.format(i_l+1, labels[i_l]))
        log.bar()
        # Insert custom levels
        log.bar(True, False)
        log('Inserting custom levels:')
        for l_params in params.levels.custom_level:
            # Skip blank levels that might be inserted
            if (l_params.depth is None) and (l_params.label is None) and (l_params.selection == []):
                continue
            # List index to insert level
            idx = l_params.depth - 1
            log('Inserting level: \n\tLabel: {}\n\tPosition: {}\n\tGroups: {}'.format(l_params.label, l_params.depth, len(l_params.selection)))
            assert len(l_params.selection) > 0, 'No selections provided for this group!'
            levels.insert(idx, l_params.selection)
            labels.insert(idx, l_params.label)
        log.bar()

    # Report
    log.subheading('Hierarchy summary: {} levels created'.format(len(levels)))
    for i_l, level in enumerate(levels):
        log.bar()
        log('Level {} ({})'.format(i_l+1, labels[i_l]))
        log.bar()
        for l in level: log('\t'+l)

    return levels, labels

def run(params, args=None):

    # Load input pickle object
    if (params.input.pickle is None):

        if os.path.exists(params.output.out_dir):
            raise Sorry('Output directory already exists: {}\nPlease delete the directory or provide a new output directory'.format(params.output.out_dir))

        assert params.analysis.table_ones_options.column_labels
        assert params.analysis.table_ones_options.r_free_label

        easy_directory(params.output.out_dir)
        log_dir = easy_directory(os.path.join(params.output.out_dir, 'logs'))
        log = Log(log_file=os.path.join(log_dir, 'fitting.log'))

        # Report parameters
        log.heading('Input command')
        input_command = ' \\\n\t'.join(['pandemic.adp'] + args)
        log(input_command)

        log.heading('Non-default parameters')
        log(master_phil.fetch_diff(source=master_phil.format(params)).as_str())

        log.heading('Input parameters')
        log(master_phil.format(params).as_str())

        with open(os.path.join(params.output.out_dir, 'params.eff'), 'w') as fh:
            fh.write(master_phil.format(params).as_str())

        log.heading('Processing input parameters')

        if len(params.input.pdb) == 0:
            raise Sorry('No structures have been provided for analysis (input.pdb=[...])')

        if params.settings.debug:
            log('DEBUG is turned on -- setting verbose=True')
            params.settings.verbose = True

        # Process/report on input parameters
        log.bar(True, False)
        log('Selected TLS amplitude model: {}'.format(params.fitting.tls_amplitude_model))
        log(tls_amplitudes_hash[params.fitting.tls_amplitude_model](10).description)
        log.bar(False, True)

        # Set precisions, etc
        log.bar()
        log('Setting TLS Matrices and Amplitude Precision')
        # Model precision
        log('TLS Matrices precision  -> {:10} / {} decimals'.format(10.**-params.fitting.precision.tls_matrices_decimals, params.fitting.precision.tls_matrices_decimals))
        TLSMatrices.set_precision(params.fitting.precision.tls_matrices_decimals)
        # Amplitude precision
        log('TLS Amplitude precision -> {:10} / {} decimals'.format(10.**-params.fitting.precision.tls_amplitude_decimals, params.fitting.precision.tls_amplitude_decimals))
        TLSAmplitudes.set_precision(params.fitting.precision.tls_amplitude_decimals)
        log.bar()

        # Set tolerances, etc
        log('Setting model tolerances')
        # TLS Model tolerance
        log('TLS Matrices tolerance  -> {:10}'.format(params.fitting.precision.tls_tolerance))
        TLSMatrices.set_tolerance(params.fitting.precision.tls_tolerance)
        # Uij Positive-semi-definiteness tolerance
        log('Uij tolerance           -> {:10}'.format(params.fitting.precision.uij_tolerance))
        MultiDatasetUijFitter.set_uij_tolerance(params.fitting.precision.uij_tolerance)
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
        # Special settings if only one model is loaded

        if len(params.input.pdb) == 1:
            log.bar(True, False)
            log('One file provided for analysis -- updating settings')
            log('For one model, it is not neccessary to refine the fitted model or recalculate R-factors (or even look for reflection data)')
            log('Setting analysis.refine_output_structures = False')
            params.analysis.refine_output_structures = False
            log('Setting analysis.calculate_r_factors = False')
            params.analysis.calculate_r_factors = False
            log('Setting analysis.calculate_electron_density_metrics = False')
            params.analysis.calculate_electron_density_metrics = False
            log('Setting input.look_for_reflection_data = False')
            params.input.look_for_reflection_data = False
            log('Setting fitting.optimisation.dataset_weights = one')
            params.fitting.optimisation.dataset_weights = 'one'
            log.bar()

        log.heading('Running setup')

        # Load input structures
        log.subheading('Building model list -- {} files'.format(len(params.input.pdb)))
        if params.input.labelling == 'foldername':
            label_func = foldername
        elif params.input.labelling == 'filename':
            label_func = filename
        else:
            raise Exception('Invalid labelling function: {}'.format(params.input.labelling))
        models = []
        for f in params.input.pdb:
            if params.input.model_type == "crystallographic":
                m = CrystallographicModel.from_file(f)
            else:
                m = AtomicModel.from_file(f)
            l = label_func(f)
            if not l:
                if len(params.input.pdb) == 1:
                    log('No label created for label function: {}'.format(params.input.labelling))
                    log('Trying to label by filename instead')
                    label_func = filename
                    l = label_func(f)
                if not l:
                    raise Sorry('No label generated using labelling function "{}"\n\tLabel {}\n\tFile {}'.format(params.input.labelling, l, f))
            m.label(tag=l)
            models.append(m)

        # Check for duplicate labels
        all_labels = [m.tag for m in models]
        unq_labels = sorted(set(all_labels))
        if len(unq_labels) != len(models):
            counts = [(l, all_labels.count(l)) for l in unq_labels]
            dups = ['{} (found {} times)'.format(l,c) for l,c in counts if c>1]
            raise Sorry('Duplicate labels generated for models: \n\t{}'.format('\n\t'.join(dups)))

        # Sort models for convenience
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
        p.write_combined_summary_graphs(out_dir_tag='model')

        if p.params.analysis.calculate_electron_density_metrics:
            p.calculate_electron_density_metrics(out_dir_tag='analysis')

        if len(models) > 1:
            p.calculate_amplitudes_dendrogram(out_dir_tag='analysis')

        if p.params.analysis.generate_diagnostic_graphs is True:
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
        # TODO only transfer the non-default params!
        #p.params = params

    # Print all errors
    p.log.heading('Errors/Warnings')
    p.warnings.extend(p.fitter.warnings)
    p.show_warnings()

    # Write HTML
    p.log.heading('Writing HTML output')
    pandemic_html.write_adp_summary(parameterisation=p, input_command=input_command)

############################################################################

if __name__=='__main__':

    from giant.jiffies import run_default
    from pandemic import module_info
    from functools import partial
    run_default._module_info = module_info
    from bamboo.common.profile import profile_code
    #a = profile_code()
    try:
        run_default(
            run                 = partial(run, args=sys.argv[1:]),
            master_phil         = master_phil,
            args                = sys.argv[1:],
            blank_arg_prepend   = blank_arg_prepend,
            program             = PROGRAM,
            description         = DESCRIPTION)
    except KeyboardInterrupt:
        print '\nProgram terminated by user'
    #a.stop()

import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys

import pathlib as pl
from libtbx import group_args
from libtbx.utils import Sorry, Failure

import numpy
numpy.set_printoptions(linewidth=numpy.inf, threshold=sys.maxsize)

############################################################################

PROGRAM = 'pandemic.adp'

DESCRIPTION = """
    Fit a hierarchical Atomic Displacement Parameter (ADP; B-factor) model to a (series of related) molecular structure(s)
"""

############################################################################

blank_arg_prepend = {'.pdb':'pdb=', '.cif':'cif=', '.pickle':'input.pickle=', '.json':'input.json='}

import libtbx.phil
master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .help = "input pdb files - with isotropic/anisotropic b-factors"
        .multiple = True
        .type = str
    labelling = *filename foldername
        .type = choice(multi=False)
        .multiple = False
    input_adp_model = isotropic anisotropic *mixed
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
    json = None
        .type = path
        .multiple = False
}
output {
    out_dir = pandemic-adp
        .help = "output directory"
        .type = str
    copy_reflection_data_to_output_folder = True
        .type = bool
    clean_up_files = compress_pdbs delete_mtzs
        .help = "Delete unnecessary output files (MTZs)"
        .type = choice(multi=True)
    json = True
        .type = bool
        .multiple = False
    pickle = False
        .type = bool
        .multiple = False
    write_isotropic_output_for_isotropic_atoms = True
        .type = bool
        .help = "Will make the output ADPs isotropic for atoms where the input atoms contained isotropic B-factors."
    html {
        make_html = True
            .help = "Write all output data in an html file"
            .type = bool
        embed_images = False
            .help = "Embed images directly into HTML file?"
            .type = bool
    }
    images {
        all = False
            .help = "Generate all graphical output - will significantly increase runtime"
            .type = bool
        pymol = none *chain
            .help = "Write residue-by-residue images of the output B-factors"
            .type = choice(multi=False)
        colour_map = 'rainbow'
            .help = "Colour map to use for making graphs. Can be any valid matplotlib colour map. Try: plasma, Pastel1, gist_earth, ocean."
            .type = str
        plot_style = 'ggplot'
            .help = "Plot style to use for making graphs. Can be any valid matplotlib plot style."
            .type = str
        font_family = None
            .help = "Font family to use for making graphs. Can be any valid matplotlib font family."
            .type = str
        font_name = None
            .help = "Specify a specific font to use for making graphs (must be available in the font family defined by font_family)."
            .type = str
        dpi = 100
            .help = "dpi of output images"
            .type = int
    }
}
model {
    b_factor_model = *echt
        .type = choice(multi=False)
    overall_selection = "(not water) and (not element H)"
        .type = str
    auto_levels = *chain phenix_find_tls_groups ss *secondary_structure *residue *backbone_sidechain atom
        .type = choice(multi=True)
    assign_het_residues_to_nearest_ss_groups = True
        .help = "Include het residues in secondary structure groups. A residue will be placed in the nearest SS group."
        .type = bool
    assign_het_residues_to_nearest_custom_groups = True
        .help = "Include het residues in nearby custom level groups. A residue will be placed in the nearest group. Useful for defining domains with bound ligands / sugars without having to add them individually."
        .type = bool
    cbeta_in_backbone = True
        .help = "Flag to control whether the c-beta atom is considered part of the backbone or the sidechain"
        .type = bool
    remove_duplicate_groups = none keep_highest_group *keep_lowest_group
        .help = "Flag to control whether identical groups that are present in multiple levels are removed (by default, the group in the lowest level will be kept)."
        .type = choice(multi=False)
    atomic_adp_level = True
        .help = 'should the atomic adp level be optimised? (it may be preferable to turn this off for small numbers of low-resolution structures)'
        .type = bool
    custom_level
        .multiple = True
    {
        insert_before = None
            .help = "Insert this level before another level? (lower in hierarchy; larger scale; smaller level numbers)"
            .type = str
            .multiple = False
        insert_after = None
            .help = "Insert this level before another level? (higher in hierarchy; smaller scale; larger level numbers)"
            .type = str
            .multiple = False
        depth = None
            .help = "Where to insert this level into the hierarchy? (after auto_levels have been determined). Inserting as '1' will make this the first level and '2' will add as the second level, etc. Any auto-generated levels will be shifted down appropriately. Suggest running with dry_run=True to check that you did it right!"
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
    echt {
        n_tls_modes_per_tls_group = 1
            .help = 'how many TLS modes to fit per group of atoms? (only worth doing for more than one dataset)'
            .type = int
        tls_amplitude_model = *simple_multiplier
            .help = 'how do amplitudes transform each TLS group in each dataset? \n\tsimple_multiplier = T->aT, L->aL, S->aS.'
            .type = choice(multi=False)
        precision
            .help = 'set various levels of precision for calculations'
        {
            tls_matrix_decimals = 16
                .help = "how many decimals of precision to be used for TLS matrices (PDB maximum is three, but more needed here due to amplitudes less than one)"
                .type = int
            tls_amplitude_decimals = 16
                .help = "how many decimals of precision to be used for TLS amplitudes"
                .type = int
        }
        eps {
            tls_matrix_eps = 1e-6
                .help = "define non-zero values for matrix elements"
                .type = float
            tls_amplitude_eps = 1e-4
                .help = "define non-zero values for amplitudes (0.01A^2 ~ 1A^2 B-factor, approximately)"
                .type = float
        }
        tolerances {
            tls_matrix_tolerance = 1e-6
                .help = "optimisation tolerance for TLS matrices (maximum allowed negative eigenvalues)"
                .type = float
            tls_amplitude_tolerance = 1e-6
                .help = "optimisation tolerance for TLS amplitudes (maximum allowed negative amplitudes)"
                .type = float
        }
    }
    adp_values {
        uij_eps = 1e-6
            .help = "define non-zero values for uij elements"
            .type = float
        uij_tolerance = 1e-6
            .help = "optimisation tolerance for Uij values (maximum allowed negative eigenvalues of Uijs)"
            .type = float
    }
}
optimisation {
    min_macro_cycles = 30
        .help = 'minimum number of fitting cycles to run (over all levels) -- must be at least 1'
        .type = int
    max_macro_cycles = 200
        .help = 'maximum number of fitting cycles to run (over all levels) -- must be at least 1'
        .type = int
    max_micro_cycles = 50
        .help = 'how many fitting cycles to run (for each level) -- must be at least 1'
        .type = int
    fit_isotropic_atoms_by_magnitude_only = True
        .type = bool
        .help = "Only use the magnitude of isotropic B-factors (True) or treat them as spherical anisotropic ADPs (False)?"
    elastic_net {
        weights {
            sum_of_amplitudes = 0.9
                .help = "weight for sum(amplitudes). minimises the number of TLS components in the optimisation. equivalent to a lasso weighting term."
                .type = float
            sum_of_squared_amplitudes = 0.1
                .help = "weight for sum(amplitudes^2). minimises the variance in the sizes of the TLS components in the optimisation. equivalent to a ridge regression term."
                .type = float
            sum_of_amplitudes_squared = 0.0
                .help = "weight for sum(amplitudes)^2. restrains the total amplitudes sum."
                .type = float
        }
        weight_scale = 1e2
            .help = "global scale applied to all optimisation_weights."
            .type = float
        weight_decay {
            decay_factor = 0.8
                .help = "amount by which optimisation_weights is scaled every cycle. must be less than 1."
                .type = float
            minimum_weight = None
                .help = "minimum weight value. weight decay will not continue once the total elastic-net weight reaches this value. The proportions between the different weights will be kept."
                .type = float
            weights_to_decay = *sum_of_amplitudes *sum_of_amplitudes_squared *sum_of_squared_amplitudes
                .help = "which optimisation_weights to decay every cycle."
                .type = choice(multi=True)
        }
    }
    intermediate_output {
        write_model_every = 10
            .help = 'output summary of hierarchical model every <n> cycles'
            .type = int
        stop_writing_after = 50
            .help = 'stop outputting graphs past this cycle'
            .type = int
        remove_previous = True
            .help = 'remove previous intermediate files at each update or at end of program.'
            .type = bool
    }
    dataset_selection {
        max_resolution = None
            .help = 'resolution limit for dataset to be used for TLS optimisation'
            .type = float
        max_datasets = None
            .help = 'takes up to this number of datasets for TLS parameter optimisation'
            .type = int
        sort_datasets_by = *resolution random name
            .help = 'how should datasets be ordered when selecting a subset of datasets for optimisation?'
            .type = choice(multi=False)
        random_seed = 0
            .help = 'random number seed for dataset ordering (and reference dataset selection)'
            .type = int
    }
    weights {
        dataset_weights = one inverse_resolution inverse_resolution_squared *inverse_resolution_cubed
            .help = 'control how datasets are weighted during optimisation?'
            .type = choice(multi=False)
        atom_weights = one *inverse_mod_U inverse_mod_U_squared inverse_mod_U_cubed
            .help = 'control how atoms are weighted during optimisation?'
            .type = choice(multi=False)
        renormalise_atom_weights_by_dataset = True
            .help = "Should atom_weights be normalised to an average of 1 for each dataset? (dataset_weights are then applied)"
            .type = bool
    }
    simplex_optimisation
        .help = 'set the various step sizes taken during simplex optimisation'
    {
        simplex_convergence = 1e-8
            .help = "cutoff for which the least-squares simplex is considered converged"
            .type = float
        vibration_delta = 1e-2
            .help = 'Vibration step size for zero-value matrices (relative scale)'
            .type = float
        libration_delta = 1e-2
            .help = 'Libration step size for zero-value matrices (relative scale)'
            .type = float
        uij_delta = 1e-6
            .help = 'Uij step size (absolute)'
            .type = float
    }
    gradient_optimisation {
        gradient_convergence = 1e-12
            .help = "cutoff for which the least-squares gradient optimisation is considered converged."
            .type = float
    }
    termination {
        max_b_rmsd = 1e-2
            .help = "Stop optimisation when the rmsd between the input and the output is changing less than this."
            .type = float
        max_b_change = 1.0
            .help = "Stop optimisation when all atoms are changing less that this value over recent cycles (determined by max_b_change_window_frac)."
            .type = float
        max_b_change_window_frac = 0.05
            .help = "Fraction of cycles over which these changes are measured."
            .type = float
        max_b_change_window_min = 3
            .type = float
    }
}
analysis {
    refine_output_structures = False
        .help = "Refine the structures after fitting (coordinates and occupancies)"
        .type = bool
    calculate_r_factors = False
        .help = "Recalculate r-factors for the input, output (and refined) models"
        .type = bool
    table_one_options {
        include scope giant.jiffies.multi_table_ones.options_phil
    }
    tls_group_clustering {
        xyz_cutoff_distance = 4.0
            .type = float
        metric = *overlap_mass
            .type = choice(multi=False)
        overlap_mass {
            comparison = simple_average *weighted_average
                .type = choice(multi=False)
        }
    }
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

def run(params, args=None):

    # Start program timer
    from pandemic.timer import Timer
    timer = Timer().start()

    if os.path.exists(params.output.out_dir):
        raise Sorry('Output directory already exists: {}\nPlease delete the directory or provide a new output directory'.format(params.output.out_dir))

    from pandemic.adp import file_system
    file_system = file_system.PandemicAdpFileSystem(output_directory=params.output.out_dir)

    from giant.logs import setup_logging
    logger = setup_logging(
        name = __name__,
        log_file = os.path.join(file_system.output_directory, 'pandemic.log'),
        warning_handler_name = 'warnings',
        debug = False,
    )

    # Extract warning handler from logging objects
    from giant.logs import get_handler_recursive
    warning_handler = get_handler_recursive(logger=logger, handler_name='warnings')
    assert warning_handler is not None

    # Report program start time
    logger('Program started at {}'.format(timer.format_time(timer.get_start_time())))

    #
    # Report parameters
    #
    if args is not None:
        logger.heading('Input command')
        input_command = ' \\\n\t'.join([PROGRAM] + args)
        logger(input_command)
    else:
        input_command = None

    from giant.phil import log_running_parameters
    log_running_parameters(
        params = params,
        master_phil = master_phil,
        logger = logger,
    )

    with open(os.path.join(file_system.output_directory, 'params-input.eff'), 'w') as fh:
        fh.write(master_phil.format(params).as_str())

    #
    # Imports
    #
    # import plots first to set matplotlib variables
    #
    from pandemic.adp import \
        plots, \
        preprocess, postprocess, \
        weights, tracking, results, \
        elastic_net, \
        hierarchy, json_manager

    #
    # Validate input parameters
    #
    preprocess.validate.validate_parameters(params=params)

    ##################################################
    #                                                #
    #         Initialise generic input tasks         #
    #                                                #
    ##################################################

    load_models = preprocess.load_models.ModelLoader(
        model_type = params.input.model_type,
        labelling = params.input.labelling,
        )

    process_input_models_task = preprocess.process_models.ProcessInputModelsTask(
        output_directory = file_system.structure_directory,
        dataset_selection_params = params.optimisation.dataset_selection,
        weights_params = params.optimisation.weights,
        analysis_params = params.analysis,
        table_one_options = params.analysis.table_one_options,
        look_for_reflection_data = params.input.look_for_reflection_data,
        copy_reflection_data_to_output_folder = params.output.copy_reflection_data_to_output_folder,
        check_column_labels = params.analysis.calculate_r_factors,
        )

    select_optimisation_datasets = preprocess.select_datasets.SelectOptimisationDatasetsTask(
        max_resolution = params.optimisation.dataset_selection.max_resolution,
        max_datasets = params.optimisation.dataset_selection.max_datasets,
        sort_datasets_by = params.optimisation.dataset_selection.sort_datasets_by,
        random_seed = params.optimisation.dataset_selection.random_seed,
        )

    extract_uijs_task = preprocess.extract_uijs.ExtractAndProcessModelUijsTask(
        expected_disorder_model = params.input.input_adp_model,
        )

    uij_weights_task = weights.UijArrayWeightsTask(
        dataset_weighting = params.optimisation.weights.dataset_weights,
        atom_weighting = params.optimisation.weights.atom_weights,
        renormalise_by_dataset = params.optimisation.weights.renormalise_atom_weights_by_dataset,
        )

    create_hierarchy_task = hierarchy.CreateHierarchicalModelTask(
        auto_levels = params.model.auto_levels,
        custom_levels = params.model.custom_level,
        overall_selection = params.model.overall_selection,
        cbeta_in_backbone = params.model.cbeta_in_backbone,
        assign_het_residues_to_nearest_ss_groups = params.model.assign_het_residues_to_nearest_ss_groups,
        assign_het_residues_to_nearest_custom_groups = params.model.assign_het_residues_to_nearest_custom_groups,
        remove_duplicate_groups = params.model.remove_duplicate_groups,
        n_cpus = params.settings.cpus,
        )

    level_optimisation_weights = weights.scale_weights(
        weights = params.optimisation.elastic_net.weights,
        scale = params.optimisation.elastic_net.weight_scale,
        )

    update_optimisation_object = elastic_net.UpdateOptimisationFunction(
        initial_weights = level_optimisation_weights,
        weight_decay_factor = params.optimisation.elastic_net.weight_decay.decay_factor,
        weights_to_update = params.optimisation.elastic_net.weight_decay.weights_to_decay,
        minimum_weight = params.optimisation.elastic_net.weight_decay.minimum_weight,
        output_directory = file_system.optimisation_directory,
        plotting_object = plots.PandemicAdpPlotter(),
        )

    ##################################################
    #                                                #
    #            Graphical/image settings            #
    #                                                #
    ##################################################

    plot_helper = plots.PlotHelper(
        colour_map_name = params.output.images.colour_map,
        plot_style = params.output.images.plot_style,
        font_family = params.output.images.font_family,
        font_name = params.output.images.font_name,
        dpi = params.output.images.dpi,
        )
    # Override defaults
    plots.PandemicAdpPlotter.helper = plot_helper

    ##################################################
    #                                                #
    #            Modify core HTML objects            #
    #                                                #
    ##################################################

    # Update class attributes - needed for loading relative image paths, etc.
    import pandemic.adp.html
    pandemic.adp.html.HtmlSummary.embed_images = params.output.html.embed_images
    pandemic.adp.html.HtmlSummary.output_dir = params.output.out_dir

    ##################################################
    #                                                #
    #         Create model-specific objects          #
    #                                                #
    ##################################################
    #
    # Create tasks (for input model) - varying code in here
    #
    if params.model.b_factor_model == 'echt':
        from pandemic.adp import echt

        echt.validate_parameters(params)

        # Create model task

        create_model_object_task = echt.hierarchy.CreateEchtModelTask(
            n_tls_modes = params.model.echt.n_tls_modes_per_tls_group,
            matrix_decimals = params.model.echt.precision.tls_matrix_decimals,
            amplitude_decimals = params.model.echt.precision.tls_amplitude_decimals,
            matrix_tolerance = params.model.echt.tolerances.tls_matrix_tolerance,
            amplitude_tolerance = params.model.echt.tolerances.tls_amplitude_tolerance,
            )

        # Optimisation functions

        optimise_tls_function = echt.optimise.tls.OptimiseTLSGroup(
            simplex_params = group_args(
                vibration_delta = params.optimisation.simplex_optimisation.vibration_delta,
                libration_delta = params.optimisation.simplex_optimisation.libration_delta,
                ),
            convergence_tolerance = params.optimisation.simplex_optimisation.simplex_convergence,
            eps_values = group_args(
                matrix_eps = params.model.echt.eps.tls_matrix_eps,
                amplitude_eps = params.model.echt.eps.tls_amplitude_eps,
                ),
            )

        if params.model.atomic_adp_level is True:
            optimise_adp_function = echt.optimise.uij.OptimiseUijValue(
                simplex_params = group_args(
                    uij_delta = params.optimisation.simplex_optimisation.uij_delta,
                    ),
                convergence_tolerance = params.optimisation.simplex_optimisation.simplex_convergence,
                uij_eps = params.model.adp_values.uij_eps,
                uij_tolerance = params.model.adp_values.uij_tolerance,
                )
        else:
            optimise_adp_function = None

        optimise_level_amplitudes_function = echt.optimise.inter_level.OptimiseInterLevelAmplitudes(
            convergence_tolerance = params.optimisation.gradient_optimisation.gradient_convergence,
            optimisation_weights =  level_optimisation_weights,
            )

        optimise_model_main = echt.optimise.OptimiseEchtModel(
            optimise_tls_function = optimise_tls_function,
            optimise_adp_function = optimise_adp_function,
            optimise_level_amplitudes_function = optimise_level_amplitudes_function,
            max_n_cycles = params.optimisation.max_micro_cycles,
            n_cpus = params.settings.cpus,
            )

        validate_model = echt.validate.ValidateEchtModel(
            uij_tolerance = params.model.adp_values.uij_tolerance,
            )

        model_analysis_task = echt.analysis.AnalyseEchtModelTask(
            output_directory = file_system.analysis_directory,
            master_phil = master_phil,
            analysis_parameters = params.analysis,
            )

        # Define / override summary classes (these take standard inputs wherever possible)

        ModelTrackingClass = echt.tracking.EchtTracking

        WriteModelSummary = echt.output.WriteEchtModelSummary
        WriteStructures = echt.output.WriteEchtDatasetStructures

        # Define / override Html classes

        import pandemic.adp.echt.html
        ModelHtmlSummary = echt.html.EchtModelHtmlSummary

    else:
        raise Sorry('Unknown Hierarchical B-factor model selected: {}'.format(params.model.b_factor_model))

    ##################################################
    #                                                #
    #            Run preprocessing tasks             #
    #                                                #
    ##################################################

    # Load and process input structures
    models = load_models(pdb_files=params.input.pdb)
    process_input_models_task.run(models)

    # Create hierarchical partitioning
    hierarchy_info = create_hierarchy_task.run(
        hierarchy = models[0].hierarchy,
        )

    # Create model objects
    model_object = create_model_object_task.run(
        models              = models,
        level_group_array   = hierarchy_info.level_group_array,
        level_group_labels  = hierarchy_info.level_group_selection_strings,
        level_labels        = hierarchy_info.level_labels,
        overall_atom_mask   = hierarchy_info.overall_atom_mask,
        )

    # Initialise model_object using json if provided
    if params.input.json is not None:
        input_json_manager = json_manager.JsonDataManager.from_json_file(
            filename = params.input.json,
            )
        input_json_manager.apply_to_model_object(model_object=model_object)

    #
    # Write output json for the starting model
    #
    if params.output.json is True:
        json_manager.write_model_as_json(
            model_object = model_object,
            output_file = os.path.join(file_system.output_directory, 'starting_model.json'),
        )

    # Extract fitted uijs
    extract_uijs_task.run(
        models              = models,
        overall_atom_mask   = hierarchy_info.overall_atom_mask,
        )
    uij_weights_task.run(
        resolutions = [m.resolution for m in models],
        uij_values  = extract_uijs_task.result.model_uij,
        dataset_labels = model_object.dataset_labels,
        )

    # Write parameters being used (may have been updated by program)
    logger.heading('Updated non-default parameters')
    logger(master_phil.fetch_diff(source=master_phil.format(params)).as_str())
    params_file = os.path.join(file_system.output_directory, 'params-running.eff')
    logger('Writing all running parameters to output folder: {}'.format(os.path.basename(params_file)))
    with open(params_file, 'w') as fh:
        fh.write(master_phil.format(params).as_str())

    # Report warnings
    warning_handler.report_new(logger_name=__name__)

    #############################################################
    #                                                           #
    #               Create output objects / tasks               #
    #                                                           #
    #   ** create here after parameters have been modified **   #
    #                                                           #
    #############################################################

    #
    # Plotting object (ensures standard colours / styles between the different output tasks)
    #
    plotting_object = plots.PandemicAdpPlotter(n_levels=model_object.n_levels)

    #
    # Create progress-tracking object (gets a special plotting object)
    #
    main_tracking_object = tracking.PandemicTrackingObject(
        output_directory = file_system.optimisation_directory,
        plotting_object = plotting_object,
        structure_factory = hierarchy_info.masked_structure_factory,
        dataset_names = model_object.dataset_labels,
        level_names = model_object.all_level_names,
        convergence_args = dict(
            max_rmsd_b = params.optimisation.termination.max_b_rmsd,
            max_delta_b = params.optimisation.termination.max_b_change,
            delta_b_window_frac = params.optimisation.termination.max_b_change_window_frac,
            delta_b_window_min = params.optimisation.termination.max_b_change_window_min,
            ),
        )

    model_tracking_object = ModelTrackingClass(
        output_directory = file_system.optimisation_directory,
        plotting_object = plotting_object,
        model_object = model_object,
        )

    #
    # Output results object
    #
    results_object = results.PandemicResultsObject(
        filename = os.path.join(file_system.output_directory, 'output_data.csv'),
        plotting_object = plotting_object,
        models = models,
        )

    # Add reference values to results object
    if params.input.reference_r_values is not None:
        results_object.add_reference_r_values(
            filename = params.input.reference_r_values,
            input_column_labels = {
                'r_free' : 'R-free',
                'r_work' : 'R-work',
                'r_gap'  : 'R-gap',
                },
            )

    #
    # Initialise generic analysis/processing tasks
    #

    # Post-processing
    post_process_task = postprocess.PostProcessTask(
        output_directory = file_system.analysis_directory,
        structure_directory = file_system.structure_directory,
        refine_structures = params.analysis.refine_output_structures,
        calculate_r_factors = params.analysis.calculate_r_factors,
        refinement_program = params.refinement.program, # TODO - refinement_params
        table_one_options = params.analysis.table_one_options,
        plotting_object = plotting_object,
        n_cpus = params.settings.cpus,
        )

    #
    # Initialise output tasks (may have been overridden above)
    #

    # Summary of the input/constructed B-factor hierarchy
    write_hierarchy_summary_task = hierarchy.WriteHierarchicalModelSummaryTask(
        output_directory = file_system.partition_directory,
        master_phil = master_phil,
        pymol_images = params.output.images.pymol,
        )

    # Summary of the model fit (error distribution, etc)
    import pandemic.adp.analysis
    generic_analysis_task = pandemic.adp.analysis.HierarchicalModelAnalysisTask(
        output_directory = file_system.analysis_directory,
        plotting_object = plotting_object,
        master_phil = master_phil,
        )

    # Write output graphs (distributions, profiles, etc) across datasets
    write_fitted_model_summary = WriteModelSummary(
        output_directory = file_system.hierarchy_directory,
        pymol_images = params.output.images.pymol,
        )

    # Copy of this task for writing out during optimisation
    from pandemic.adp.utils import TaskMultipleRunWrapper
    write_fitted_model_summary_intermediate = TaskMultipleRunWrapper(
        output_directory = file_system.optimisation_directory,
        task = write_fitted_model_summary,
        task_output_directory_prefix = 'optimised_model_cycle_',
        )
    # Turn off pymol images
    write_fitted_model_summary_intermediate.task.pymol_images = None

    # Write output structures for each dataset
    write_output_structures = WriteStructures(
        output_directory = file_system.structure_directory,
        )

    ################################
    #                              #
    #  Create output HTML objects  #
    #                              #
    ################################

    # Writes html from input summaries
    write_html_task = pandemic.adp.html.WriteHtmlSummaryTask(
        output_directory = file_system.output_directory,
        )

    # Html objects that are added to iteratively and called at the end to generate HTML
    html_tracking_summary = pandemic.adp.html.HtmlSummaryConcatenator(
        title = 'Model Optimisation',
        summaries = [
            main_tracking_object.as_html_summary(),
            model_tracking_object.as_html_summary(),
            update_optimisation_object.as_html_summary(),
        ],
    )

    # Main model summary -- defined as custom class with standard inputs for now -- could be changed
    html_model_summary = ModelHtmlSummary(
        hierarchy_summary_task = write_hierarchy_summary_task,
        model_summary_task = write_fitted_model_summary,
        model_object = model_object,
        isotropic_mask = extract_uijs_task.result.isotropic_mask,
        parameters = params,
    )

    import pandemic.adp.postprocess.html as postprocesshtml
    html_postprocess_summary = postprocesshtml.PostProcessingHtmlSummary(
        results_object = results_object,
        analysis_task = post_process_task,
    )

    html_analysis_summary = pandemic.adp.html.HtmlSummaryCollator(
        title = 'Model Analysis',
        alt_title = 'Analysis',
        summaries = [
            generic_analysis_task.as_html_summary(),
            model_analysis_task.as_html_summary(),
        ],
    )

    import pandemic.adp.html.parameters as parametershtml
    html_parameter_summary = parametershtml.ParameterHtmlSummary(
        input_command = input_command,
        master_phil = master_phil,
        running_params = params,
    )

    ################################
    #                              #
    #   Preprocessing summaries    #
    #                              #
    ################################

    # Write summary of hierarchical partitioning
    write_hierarchy_summary_task.run(
        reference_hierarchy           = models[0].hierarchy,
        level_group_array             = hierarchy_info.level_group_array,
        level_group_selection_strings = hierarchy_info.level_group_selection_strings,
        level_labels                  = hierarchy_info.level_labels,
        overall_atom_mask             = hierarchy_info.overall_atom_mask,
        plotting_object               = plotting_object,
        )

    # Report warnings
    warning_handler.report_new(logger_name=__name__)

    ################################
    #                              #
    #         End of setup         #
    #                              #
    ################################

    #
    # Exit if requested
    #
    if params.settings.dry_run:
        logger.heading('Exiting after initialisation: dry_run=True')
        raise SystemExit('Program exited normally')

    ####################################
    #                                  #
    #  Prep different isotropic masks  #
    #                                  #
    ####################################

    if (params.optimisation.fit_isotropic_atoms_by_magnitude_only is True):
        optimise_isotropic_mask = extract_uijs_task.result.isotropic_mask
    else:
        optimise_isotropic_mask = extract_uijs_task.result.isotropic_mask.as_fully_anisotropic()

    if params.output.write_isotropic_output_for_isotropic_atoms is True:
        output_isotropic_mask = extract_uijs_task.result.isotropic_mask
    else:
        output_isotropic_mask = extract_uijs_task.result.isotropic_mask.as_fully_anisotropic()

    ################################
    #                              #
    #  Run main optimisation task  #
    #                              #
    ################################

    logger.heading('Optimising Hierarchical Disorder Model', spacer=True)

    main_tracking_object.set_target(
        uij_target = extract_uijs_task.result.model_uij,
        uij_target_weights = uij_weights_task.result.total_weight_array,
        )

    main_tracking_object.n_cycle = 0

    try:

        for _ in range(params.optimisation.max_macro_cycles):

            # Increment before each cycle -- ensures new counter for this cycle.
            main_tracking_object.n_cycle += 1

            logger.subheading('Macrocycle {}'.format(main_tracking_object.n_cycle), spacer=True)

            update_optimisation_object.update(
                model_optimisation_function = optimise_model_main,
                n_cycle = main_tracking_object.n_cycle,
                )

            model_object = optimise_model_main(
                model_object        = model_object,
                level_group_tree    = hierarchy_info.level_group_tree,
                uij_target          = extract_uijs_task.result.model_uij,
                uij_target_weights  = uij_weights_task.result.total_weight_array,
                uij_isotropic_mask  = optimise_isotropic_mask,
                tracking_object     = main_tracking_object,
                )

            model_tracking_object.update(
                model_object = model_object, # TODO REMOVE
                n_cycle = main_tracking_object.n_cycle,
                )

            # Write graphs
            logger.heading('Writing tracking output')
            main_tracking_object.write_output()
            model_tracking_object.write_output()
            update_optimisation_object.write_output()

            if main_tracking_object.is_converged():
                if (main_tracking_object.n_cycle >= params.optimisation.min_macro_cycles):
                    logger.heading('Terminating optimisation -- model has converged', spacer=True)
                    break
            else:
                logger.subheading('Model has not yet reached convergence criteria')

            # Output intermediate model
            if params.output.json is True:
                json_manager.write_model_as_json(
                    model_object = model_object,
                    output_file = str(
                        pl.Path(file_system.optimisation_directory) / 'optimised_model_cycle_{n_cycle:03d}.json'.format(
                            n_cycle = main_tracking_object.n_cycle,
                            )
                        ),
                )

            if params.optimisation.intermediate_output.write_model_every:
                if (params.optimisation.intermediate_output.stop_writing_after is not None) and (main_tracking_object.n_cycle > params.optimisation.intermediate_output.stop_writing_after):
                    pass
                elif (main_tracking_object.n_cycle % params.optimisation.intermediate_output.write_model_every) == 0:
                    logger.subheading('Macrocycle {}'.format(main_tracking_object.n_cycle) + ' - Writing intermediate model summary')
                    # Write out the model during optimisation
                    model_summary_output = write_fitted_model_summary_intermediate(
                        output_directory_suffix = '{:03d}'.format(main_tracking_object.n_cycle),
                        remove_previous = params.optimisation.intermediate_output.remove_previous,
                        # Normal arguments (passed to task)
                        overall_atom_mask = hierarchy_info.overall_atom_mask,
                        level_group_array = hierarchy_info.level_group_array,
                        model_object = model_object,
                        isotropic_mask = output_isotropic_mask,
                        reference_model = models[0],
                        uij_target = extract_uijs_task.result.model_uij,
                        plotting_object = plotting_object,
                        )

    finally:
        # Ensure processes used for optimisation are closed
        optimise_model_main.close_processes()

    logger.heading('Optimisation finished', spacer=True)

    #
    # Validate output model
    #
    logger.heading('Validating output ADPs')
    validate_model(
        model_object = model_object,
        reference_hierarchy = models[0].hierarchy,
        overall_atom_mask = hierarchy_info.overall_atom_mask,
        )

    # Report warnings
    warning_handler.report_new(logger_name=__name__)

    ################################
    #                              #
    # Write model output summaries #
    #                              #
    ################################

    #
    # Output graphs/csvs
    #
    logger.heading('Writing disorder model summary')
    model_summary_output = write_fitted_model_summary(
        overall_atom_mask = hierarchy_info.overall_atom_mask,
        level_group_array = hierarchy_info.level_group_array,
        model_object = model_object,
        isotropic_mask = output_isotropic_mask,
        reference_model = models[0],
        uij_target = extract_uijs_task.result.model_uij,
        plotting_object = plotting_object,
        )
    # Tidy up output from the intermediates...
    if params.optimisation.intermediate_output.remove_previous:
        write_fitted_model_summary_intermediate.delete_output()

    #
    # Write fitted structures for each of the datasets
    #
    logger.heading('Writing output structures for each dataset')
    fitted_structures = write_output_structures(
        level_group_array = hierarchy_info.level_group_array,
        model_object = model_object,
        isotropic_mask = output_isotropic_mask,
        models = models,
        overall_selection = params.model.overall_selection,
        overall_selection_bool = hierarchy_info.overall_atom_mask,
        )

    #
    # Pickle model object
    #
    if params.output.pickle is True:
        from libtbx import easy_pickle
        easy_pickle.dump(os.path.join(file_system.output_directory, 'model.pickle'), model_object)

    #
    # Write model JSON
    #
    if params.output.json is True:
        json_manager.write_model_as_json(
            model_object = model_object,
            output_file = os.path.join(file_system.output_directory, 'optimised_model.json'),
        )

    ################################
    #                              #
    #     Analyse fitted model     #
    #                              #
    ################################

    #
    # Generic model analysis task
    #
    generic_analysis_task.run(
        uij_target =  extract_uijs_task.result.model_uij,
        uij_target_weights = uij_weights_task.result.total_weight_array,
        uij_isotropic_mask = output_isotropic_mask,
        model_object = model_object,
        model_hierarchy_info = hierarchy_info,
        model_hierarchy_files = write_hierarchy_summary_task.result.output_files,
        reference_hierarchy = models[0].hierarchy,
        )

    #
    # Artisinal model analysis task
    #
    model_analysis_task.run(
        model_object = model_object,
        model_summary_output = model_summary_output,
        )

    ################################
    #                              #
    #    Post-processing (slow)    #
    #                              #
    ################################

    #
    # Post-process: Refine fitted structures / Calculate R-factors / etc.
    #
    logger.heading('Post-processing fitted structures')
    post_process_task.run(
        dataset_labels = model_object.dataset_labels,
        input_structures_dict = process_input_models_task.input_structures,
        fitted_structures_dict = fitted_structures['complete_structures'],
        input_reflection_data_dict = process_input_models_task.input_reflection_data,
        cif_files = params.refinement.cif,
        results_object = results_object,
        )

    ################################
    #                              #
    #     Write results + html     #
    #                              #
    ################################

    #
    # Write output csvs/graphs
    #
    results_object.write()
    model_tracking_object.write_output()

    #
    # Write HTML
    #

    if params.output.html.make_html is True:

        write_html_task.run(
            html_objects = [
                html_tracking_summary,
                html_model_summary,
                html_postprocess_summary,
                html_analysis_summary,
                html_parameter_summary,
            ],
        )

    if params.output.clean_up_files:

        from pandemic.adp import file_system
        tidy = file_system.TidyOutputFolder(
            compress_pdbs = ('compress_pdbs' in params.output.clean_up_files),
            delete_mtzs   = ('delete_mtzs'   in params.output.clean_up_files),
            )
        tidy(output_directory=params.output.out_dir)

    # Report all errors accumulated over the lifetime of the program
    warning_handler.report_all(logger_name=__name__)

    logger(timer.report_str())

    return None

############################################################################

def run_pandemic_adp(args):

    from functools import partial
    from giant.jiffies import run_default
    from pandemic import module_info
    run_default._module_info = module_info

    run_default(
        run                 = partial(run, args=args),
        master_phil         = master_phil,
        args                = args,
        blank_arg_prepend   = blank_arg_prepend,
        program             = PROGRAM,
        description         = DESCRIPTION,
        )

if __name__=='__main__':

    from giant.profile import profile_code
    #a = profile_code()
    try:
        run_pandemic_adp(args=sys.argv[1:])
    except Exception as e:
        import traceback
        logger(traceback.format_exc())
        raise
    except KeyboardInterrupt:
        logger('\nProgram terminated by user')
    #finally:
    #    a.stop()


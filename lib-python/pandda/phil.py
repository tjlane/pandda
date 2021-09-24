import libtbx.phil

pandda_phil = libtbx.phil.parse("""
pandda
{
    input
        .help = "File names"
    {
        data_dirs = None
            .type = str
        pdb_style = final.pdb
            .type = str
        mtz_style = None
            .type = str
        lig_style = *.cif
            .type = str
        regex
            .help = "Advanced dataset labelling regexs"
        {
            dir_regex = None
                .type = str
            pdb_regex = None
                .type = str
            mtz_regex = None
                .type = str
        }
        filter
            .help = "Provide a dataset to filter input datasets against"
        {
            pdb = None
                .type = path
        }
        reference
            .help = "Manually define reference dataset to align+scale all other datasets to"
        {
            pdb = None
                .type = path
            mtz = None
                .type = path
        }
        flags
            .help = "Flags for individual datasets"
        {
            ignore_datasets = None
                .help = "Reject these datasets, don't even load them - comma separated list of dataset tags"
                .type = str
            only_datasets = None
                .help = "Only load these datasets, don't load any others - comma separated list of dataset tags"
                .type = str
            ground_state_datasets = None
                .help = "Define a set of known 'ground-state' (e.g. unbound) datasets. ONLY these datasets will be used for characterising the ground-state electron density."
                .type = str
            exclude_from_z_map_analysis = None
                .help = "Don't analyse these datasets, only use them to build the distributions - comma separated list of dataset tags"
                .type = str
            exclude_from_characterisation = None
                .help = "Don't use these datasets to characterise the ground-state electron density, only analyse them - comma separated list of dataset tags"
                .type = str
            test = None
                .help = "Only load these datasets, don't load any others - comma separated list of dataset tags"
                .type = str
            train = None
                .help = "Only load these datasets, don't load any others - comma separated list of dataset tags"
                .type = str
            not_test = None
                .help = "Only load these datasets, don't load any others - comma separated list of dataset tags"
                .type = str
            not_train = None
                .help = "Only load these datasets, don't load any others - comma separated list of dataset tags"
                .type = str
        }
    }

    output
        .help = "Output directory"
    {
        out_dir = './pandda'
            .type = path
        overwrite = False
            .type = bool
        dataset_prefix = ''
            .help = "Prefix to be attached to each dataset name"
            .type = str
        output_maps_for = all *events first_dataset_only
            .help = "Which datasets to write output maps for" 
            .type = choice(multi=False)
    }

    params
        .help = "Algorithm Parameters"
    {
        analysis
            .help = "Settings to control the selection of datasets"
        {
            high_res_increment = 0.05
                .help = 'Increment of resolution shell for map analysis'
                .type = float
            high_res_upper_limit = 0.0
                .help = 'Highest resolution limit (maps are never calulcated above this limit)'
                .type = float
            high_res_lower_limit = 4.0
                .help = 'Lowest resolution limit (datasets below this are ignored)'
                .type = float
        }
        diffraction_data
        {
            structure_factors = 2FOFCWT_iso-fill,PH2FOFCWT_iso-fill FWT,PHWT 2FOFCWT_fill,PH2FOFCWT_fill 2FOFCWT,PH2FOFCWT
                .type = strings
            checks
                .help = "checks on the mtz file data provided for each dataset"
            {
                all_data_are_valid_values = true
                    .help = "check that all reflections in the diffraction data have valid values (are not zero or n/a)"
                    .type = bool
                low_resolution_completeness = 4.0
                    .help = "check that diffraction data is 100% complete up to this resolution cutoff. missing reflections at low resolution may seriously degrade analysis quality. set to none to turn off this check."
                    .type = float
            }
            scaling 
            {
                apply_b_factor_scaling = True
                    .help = "Apply b-factor scaling to reflections? (reciprocal space)"
                    .type = bool
            }
        }
        alignment
            .help = "Settings to control the alignment of the structures"
        {
            method = global *local
                .help = "How should the structures be aligned? 'global' is fast, but requires high structural conservation, whereas local is slower and accounts for structural variation"
                .type = choice
        }
        filtering
            .help = "Settings to control when datasets are rejected from the analysis"
        {
            max_rmsd_to_reference = 1.5
                .help = "Reject datasets that have a calpha rmsd of greater than this to the reference (after global alignment)"
                .type = float
                .multiple = False
            max_rfree = 0.4
                .help = 'Maximum allowed rfree for a structure (datasets above this are rejected)'
                .type = float
                .multiple = False
            same_space_group_only = True
                .help = "Reject datasets that are a different spacegroup to the reference/filter dataset - NOT YET IMPLEMENTED - MUST BE SET TO TRUE"
                .type = bool
            similar_models_only = False
                .help = "Reject datasets that have a different model composition to the reference/filter dataset. All models must have the same number and identity of atoms."
                .type = bool
        }
        excluding
            .help = "Parameters to control when datasets are automatically excluded from characterisation"
        {
            max_wilson_plot_z_score = 5.0
                .help = "Maximum Z-score for RMSD of dataset diffraction wilson curve compared to the reference dataset. Z-scores calculated relative to the rest of the datasets."
                .type = float
                .multiple = False
        }
        maps
            .help = "Settings to control how maps are generated and analysed"
        {
            resolution_factor = 0.25
                .help = 'Sampling factor for fft-ing the maps'
                .type = float
            grid_spacing = 0.5
                .help = 'Spacing of the grid points in the sampled maps (A) - fixed across resolutions'
                .type = float
            density_scaling = none *sigma volume
                .help = "Apply scaling to electron density? (real-space)"
                .type = choice
        }
        masks
            .help = "Parameters to control the masking of grid points around the protein"
        {
            mask_pdb = None
                .help = "A PDB to mask the grid against (if none provided, use reference dataset)"
                .type = str
                .multiple = False
            mask_selection_string = None
                .help = "Define a custom region of the protein to load electron density for (selection + outer_mask - inner_mask). If not defined, uses all protein atoms."
                .type = str
            align_mask_to_reference = False
                .help = "If masks.pdb is supplied, does it require alignment to the reference structure?
                            If selecting a fragment of the structure, masks.pdb must already be aligned prior to running pandda (can't align fragments)."
                .type = bool
            outer_mask = 6.0
                .help = 'include region within outer_mask of the atoms defined by selection_string (or protein if selection_string undefined)'
                .type = float
            inner_mask = 1.8
                .help = 'exclude region within inner_mask of the atoms defined by selection_string (or protein if selection_string undefined)'
                .type = float
            inner_mask_symmetry = 3.0
                .help = 'exclude region within inner_mask of the (symmetry) atoms defined by selection_string (or protein if selection_string undefined)'
                .type = float
        }
        statistical_maps
            .help = "Settings to control the calculation of z-maps"
        {
            min_build_datasets = 30
                .help = 'Minimum number of datasets needed to build distributions'
                .type = int
            max_build_datasets = 60
                .help = 'Maximum number of datasets used to build distributions'
                .type = int
            average_map_type = *mean median
                .help = 'Which statistical map should the uncertainties and Z-map be calculated from?
                    mean: mean map      (simple averaging)
                    medn: median map    (less sensitive to outliers)'
                .type = choice(multi=False)
            fit_mu = True
                .type = bool
            fit_sigma_adjusted = True
                .type = bool
            z_map_type = uncertainty *adjusted+uncertainty quantile
                .help = 'Type of Z-map to calculate'
                .type = choice
        }
        z_map_analysis
            .help = "Settings to control the finding of blobs"
            .expert_level = 1
        {
            clustering_method = *agglomerative_hierarchical 
                .help = "Method of clustering to use"
                .type = choice(multi=False)
            contour_level = 2.5
                .help = 'Contour level when looking for blobs'
                .type = float
            negative_values = False
                .help = 'Look for large negative values as well?'
                .type = bool
            min_blob_volume = 10.0
                .help = 'Blob volume filter for detecting blobs'
                .type = float
            min_blob_z_peak = 3.0
                .help = 'Blob Z-peak filter for detecting blobs'
                .type = float
            masks {
                selection_string = None
                    .help = 'Define a further mask to define a region for z-map analysis (selection + outer_mask - inner_mask). If not defined, uses the same region as defined in params.masks.'
                    .type = str
                outer_mask = None
                    .help = 'include region within outer_mask of the atoms defined by selection_string (or protein if selection_string undefined)'
                    .type = float
                inner_mask = None
                    .help = 'exclude region within inner_mask of the atoms defined by selection_string'
                    .type = float
            }
            agglomerative_hierarchical {
                clustering_cutoff = 2.0
                    .help = "cutoff range for agglomerative clustering"
                    .type = float
            }
        }
        background_correction
            .help = "Parameters to control the estimation of feature background corrections"
            .expert_level = 3
        {
            min_bdc = 0.0
                .type = float
                .help = 'Minimum background correction estimate'
            max_bdc = 1.0
                .type = float
                .help = 'Maximum background correction estimate'
            increment = 0.01
                .type = float
                .help = 'Resolution of background correction estimation'
            output_multiplier = 1.0
                .type = float
                .help = 'Empirical multiplier to be applied to the contrast-estimated value of 1-BDC (truncated after multiplication to range 0-1)'
        }
    }

    results
        .help = "Change ordering/filtering of the output data"
    {
        events {
            order_by = z_peak z_mean *cluster_size
                .help = "How should events be ordered within each site?"
                .type = choice
        }
        sites {
            order_by = *num_events
                .help = "How should sites be ordered?"
                .type = choice
        }
    }

    exit_flags
        .help = "Flags for terminating the program early"
    {
        dry_run = False
            .help = "Find input files on disk, but exit before loading any datasets."
            .type = bool
        setup_only = False
            .help = "Load datasets and create grid partition, but exit before analysis."
            .type = bool
        calculate_first_average_map_only = False
           .help = "Will calculate the highest resolution average map and then exit - used for initial reference modelling."
            .type = bool
    }

    processing
        .help = "Change how the PanDDA is calculated"
    {
        cpus = 1
            .help = "Number of local cpus to use"
            .type = int
        backend = serial *parallel_joblib 
            .help = "What multiprocessing backend to use"
            .type = choice
        pandda_backend = *numpy
            .help = "What technology to use for calculating pandda statistics"
        process_shells = *serial  
            .help = "How to process shells"
            .type = choice
        remote_nodes = None
            .help = "Number of remote nodes to use (e.g. with luigi)"
        cpus_per_remote_node = None
            .help = "Number of cpus to use on a remote machine (e.g. with luigi)"
        h_vmem = 100
            .help = "How to process dicts"
            .type = int
        m_mem_free = 5
            .help = "How to process dicts"
            .type = int
    }

}
settings
    .help = "General Settings"
{
    verbose = False
        .type = bool
    plotting {
        backend = 'agg'
            .help = "Backend to use in matplotlib"
            .type = str
    }
}
""", process_includes=True)


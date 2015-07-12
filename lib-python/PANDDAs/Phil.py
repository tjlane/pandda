
pandda_phil_def = """
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
            reference
                .help = "Manually define reference structure"
            {
                pdb = None
                    .type = path
                mtz = None
                    .type = path
                ampl_label = None
                    .type = str
                phas_label = None
                    .type = str
            }
        }
        output
            .help = "Output directory"
        {
            outdir = './pandda'
                .type = path
            dataset_prefix = ''
                .type = str
        }
        method
            .help = "High Level control of algorithm"
        {
            dry_run = False
                .help = "Setup pandda, print working arguments, and exit"
                .type = bool
            reload_existing_datasets = False
                .help = "Reload existing datasets? (Time-consuming) - if False, will only load new datasets (unprocessed datasets)"
                .type = bool
            recalculate_statistical_maps = False
                .help = "Recalculate all of the statistical maps? (Time-consuming) - if False, looks for existing statistical maps and uses those (reverts to True if none are found)"
                .type = bool
            reprocess_existing_datasets = False
                .help = "Reprocess existing datasets? (Time-consuming) - if False, will only calculate z-maps for new datasets (unprocessed datasets)"
                .type = bool
        }
        params
            .help = "Algorithm Parameters"
        {
            alignment
                .help = "Settings to control the alignment of the structures"
            {
                method = global *local
                    .type = choice
                rmsd_cutoff = 1.0
                    .help = "Reject datasets that have a calpha rmsd of greater than this to the reference"
                    .type = float
            }
            maps
                .help = "Settings to control the generation of the dataset maps"
            {
                maps_to_analyse = *2FOFC FOFC
                    .type = choice
                ampl_label = FWT *2FOFCWT DELFWT FOFCWT
                    .type = choice
                phas_label = PHWT *PH2FOFCWT PHDELWT PHFOFCWT
                    .type = choice
                scaling = none *sigma volume
                    .type = choice
                resolution_factor = 0.33
                    .help = 'Sampling factor for fft-ing the maps'
                    .type = float
                grid_spacing = 0.5
                    .help = 'Spacing of the grid points (fixed across resolutions)'
                    .type = float
                padding = 3
                    .help = "Padding around the edge of the maps"
                    .type = float
            }
            masks
                .help = "Parameters to control the masking of grid points around the protein"
            {
                inner_mask = 1.8
                    .help = "Points are masked within this distance of protein atoms"
                    .type = float
                outer_mask = 6
                    .help = "Points are masked outside this distance of protein atoms"
                    .type = float
            }
            analysis
                .help = "Settings to control the selection of datasets"
            {
                min_datasets = 40
                    .help = 'Minimum number of datasets needed to build distributions'
                    .type = int
                max_datasets = 75
                    .help = 'Maximum number of datasets used to build distributions'
                    .type = int
                high_res_increment = 0.05
                    .help = 'Increment of resolution shell for map analysis'
                    .type = float
                high_res_lower_limit = 4.0
                    .help = 'Lowest resolution limit (datasets below this are rejected)'
                    .type = float
                high_res_upper_limit = 0.0
                    .help = 'Highest resolution limit (maps are never calulcated above this limit)'
                    .type = float
                dynamic_res_limits = True
                    .help = 'Allow the limits to change depending on the dataset resolution ranges'
                    .type = bool
                max_rfree = 0.4
                    .help = 'Maximum allowed rfree for a structure (datasets above this are rejected)'
                    .type = float
                no_analyse = None
                    .help = 'Don\'t analyse these datasets, only use them to build the distributions - comma separated list of dataset tags'
                    .type = str
                no_build = None
                    .help = 'Don\'t use these datasets to build the distributions, only analyse them - comma separated list of dataset tags'
                    .type = str
                ignore_datasets = None
                    .help = 'Reject these datasets, don\'t even load them - comma separated list of dataset tags'
                    .type = str
            }
            blob_search
                .help = "Settings to control the finding of blobs"
            {
                contour_level = 2.5
                    .help = 'Contour level when looking for blobs'
                    .type = float
                min_blob_volume = 10.0
                    .help = 'Blob volume filter for detecting blobs'
                    .type = float
                min_blob_z_peak = 3.0
                    .help = 'Blob Z-peak filter for detecting blobs'
                    .type = float
                blobs_to_image = 2
                    .help = 'Maximum number of blobs to image per dataset'
                    .type = int
                clustering
                    .help = "Settings to control the clustering of blob points"
                {
                    criterion = *distance
                        .type = choice
                        .optional = False
                    metric    = *euclidean
                        .type = choice
                        .optional = False
                    linkage   = *single
                        .type = choice
                        .optional = False
                }
            }
            occupancy_estimation
                .help = "Parameters to control the estimation of feature occupancies"
            {
                min_occ = 0.0
                    .type = float
                    .help = 'Minimum occupancy estimate'
                max_occ = 1.0
                    .type = float
                    .help = 'Maximum occupancy estimate'
                occ_increment = 0.01
                    .type = float
                    .help = 'Resolution of occupancy estimation'
            }
        }
        settings
            .help = "General Settings"
        {
            verbose = True
                .type = bool
            cpus = 1
                .type = int
            max_memory = 25.0
                .type = float
            testing = False
                .help = "Validate hits against a known list of ligands"
                .type = bool
            plot_graphs = True
                .help = "Output graphs using matplotlib"
                .type = bool
            developer
                .help = "Development Settings (Not needed by most users)"
            {
                write_all_z_maps = False
                    .help = "Output all possible types of Z-maps (Takes a lot of space on disk)"
                    .type = bool
            }
            pickling
                .help = "Pickle Settings"
            {
                full_pickle = False
                    .help = "Pickle the entire PANDDA object for reloading - takes a lot of space"
                    .type = bool
            }
        }
    }
"""

pandda_twiddle_def = """
    twiddle
    {
        direction = toref *fromref
            .type = choice
        file
        {
            input  = None
                .type = path
            output = None
                .type = path
        }
        dataset_pickle = ./pickles/dataset.pickle
            .type = path
    }

"""


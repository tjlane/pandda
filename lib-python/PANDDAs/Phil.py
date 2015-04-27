
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
        }
        method
            .help = "High Level control of algorithm"
        {
            dry_run = False
                .help = "Setup pandda, print working arguments, and exit"
                .type = bool
            recalculate_statistical_maps = False
                .help = "Recalculate all of the statistical maps? (Time-consuming)"
                .type = bool
            reprocess_all_datasets = False
                .help = "Recalculate and process all dataset z-maps (For different contour levels, etc)"
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
                not_used_obs_ampl_label = FWT *2FOFCWT
                    .type = choice
                not_used_obs_phas_label = PHWT *PH2FOFCWT
                    .type = choice
                not_used_diff_ampl_label = DELFWT *FOFCWT
                    .type = choice
                not_used_diff_phas_label = PHDELWT *PHFOFCWT
                    .type = choice
                scaling = none *sigma volume
                    .type = choice
                resolution_factor = 0.33
                    .type = float
                border_padding = 5
                    .type = float
            }
            masks
                .help = "Parameters to control the masking of grid points around the protein"
            {
                inner_mask = 1.8
                    .help = "Points are masked with this distance of protein atoms"
                    .type = float
                outer_mask = 5
                    .help = "Points are masked outside this distance of protein atoms"
                    .type = float
            }
            analysis
                .help = "Settings to control the selection of datasets"
            {
                max_rfree = 0.4
                    .help = 'Maximum allowed rfree for a structure (datasets above this are rejected)'
                    .type = float
                high_res_lower_limit = 3.0
                    .help = 'Lowest resolution limit (datasets below this are rejected)'
                    .type = float
                high_res_upper_limit = 0.0
                    .help = 'Highest resolution limit (maps are never calulcated above this limit)'
                    .type = float
                high_res_increment = 0.1
                    .help = 'Increment of resolution shell for map analysis'
                    .type = float
                max_datasets = 1000
                    .help = 'Maximum number of datasets to process at once'
                    .type = int
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


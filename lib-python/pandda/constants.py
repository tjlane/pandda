
STRUCTURE_MASK_NAMES = [    'bad structure - non-identical structures',
                            'bad structure - different space group'    ]
CRYSTAL_MASK_NAMES   = [    'bad crystal - isomorphous crystal',
                            'bad crystal - isomorphous structure',
                            'bad crystal - space group',
                            'bad crystal - data quality',
                            'bad crystal - data correlation',
                            'bad crystal - rfree'  ]
REJECT_MASK_NAMES    = [    'rejected - total',
                            'rejected - crystal',
                            'rejected - structure',
                            'rejected - unknown'  ]
FLAG_MASK_NAMES      = [    'noisy zmap',
                            'analysed'  ]

DATASET_INFO_FIELDS  = [    'high_resolution',
                            'low_resolution',
                            'r_work',
                            'r_free',
                            'rmsd_to_reference',
                            'space_group',
                            'uc_a',
                            'uc_b',
                            'uc_c',
                            'uc_alpha',
                            'uc_beta',
                            'uc_gamma',
                            'uc_vol',
                            'rejection_reason'                ]
DATASET_MAP_FIELDS   = [    'analysed_resolution',
                            'map_uncertainty',
                            'obs_map_mean',
                            'obs_map_rms',
                            'z_map_mean',
                            'z_map_std',
                            'z_map_skew',
                            'z_map_kurt'            ]
DATASET_EVENT_FIELDS = [    'site_idx',
                            '1-BDC',
                            'z_peak',
                            'z_mean',
                            'cluster_size',
                            'x','y','z',
                            'refx','refy','refz'    ]
SITE_TABLE_FIELDS    = [    'centroid',
                            'num_events',
                            'nearest_residue 1',
                            'nearest_residue 2',
                            'nearest_residue 3',
                            'native_centroid',
                            'near_crystal_contacts'     ]

class PanddaMaskNames:
    structure_mask_names = STRUCTURE_MASK_NAMES
    crystal_mask_names   = CRYSTAL_MASK_NAMES
    reject_mask_names    = REJECT_MASK_NAMES
    flag_mask_names      = FLAG_MASK_NAMES
    custom_mask_names    = ['no_analyse', 'no_build']
    all_mask_names       = structure_mask_names + crystal_mask_names + reject_mask_names + flag_mask_names + custom_mask_names

class PanddaTableFields:
    all_dataset_fields     = DATASET_INFO_FIELDS
    all_dataset_map_fields = DATASET_MAP_FIELDS
    all_event_fields       = DATASET_EVENT_FIELDS
    all_site_fields        = SITE_TABLE_FIELDS

########################################################################

class PanddaDatasetFilenames:
    # Input data...
    input_structure        = '{!s}-pandda-input.pdb'
    input_data             = '{!s}-pandda-input.mtz'
    # Dataset Information
    dataset_info           = '{!s}-info.csv'
    dataset_log            = '{!s}.log'
    dataset_pickle         = 'dataset.pickle'
    # Structure Files...
    aligned_structure      = '{!s}-aligned-structure.pdb'
    symmetry_copies        = '{!s}-aligned-sym-contacts.pdb'
    # Modelled Structures...
    modelled_structure     = '{!s}-pandda-model.pdb'
    ensemble_structure     = '{!s}-ensemble-model.pdb'
    # Dataset Maps
    sampled_map            = '{!s}-aligned-map.ccp4'
    # Z-Maps...
    mean_diff_map          = '{!s}-difference-from-mean.ccp4'
    z_map                  = '{!s}-z_map.ccp4'
    z_map_naive            = '{!s}-z_map_naive.ccp4'
    z_map_naive_norm       = '{!s}-z_map_naive_normalised.ccp4'
    z_map_uncertainty      = '{!s}-z_map_uncertainty.ccp4'
    z_map_uncertainty_norm = '{!s}-z_map_uncertainty_normalised.ccp4'
    z_map_corrected        = '{!s}-z_map_adjusted.ccp4'
    z_map_corrected_norm   = '{!s}-z_map_adjusted_normalised.ccp4'
    # Events...
    event_map              = '{!s}-event_{!s}_1-BDC_{!s}_map.ccp4'
    # Masks...
    high_z_mask            = '{!s}-high_z_mask.ccp4'
    grid_mask              = '{!s}-masked_grid.ccp4'
    # Ligands...
    ligand_coordinates     = '{!s}-ligand.pdb'
    ligand_restraints      = '{!s}-ligand.cif'
    ligand_image           = '{!s}-ligand.png'
    # Native (Output) Maps
    native_obs_map         = '{!s}-observed.native.ccp4'
    native_z_map           = '{!s}-z_map.native.ccp4'
    native_event_map       = '{!s}-event_{!s}_1-BDC_{!s}_map.native.ccp4'
    # Misc...
    z_peaks_csv            = '{!s}-z_map_peaks.csv'
    pymol_script           = 'load_maps.pml'
    ccp4mg_script          = 'ccp4mg_{!s}_{!s}.py'
    ccp4mg_png             = 'blob_{!s}_{!s}.png'

class PanddaDatasetPNGFilenames:
    s_map_png                   = '{!s}-obsv_map_dist.png'
    d_mean_map_png              = '{!s}-diff_mean_map_dist.png'
    z_map_naive_png             = '{!s}-z_map_dist_naive.png'
    z_map_naive_norm_png        = '{!s}-z_map_dist_naive_normalised.png'
    z_map_uncertainty_png       = '{!s}-z_map_dist_uncertainty.png'
    z_map_uncertainty_norm_png  = '{!s}-z_map_dist_uncertainty_normalised.png'
    z_map_corrected_png         = '{!s}-z_map_dist_adjusted.png'
    z_map_corrected_norm_png    = '{!s}-z_map_dist_adjusted_normalised.png'
    z_map_qq_plot_png           = '{!s}-z_map_dist_qq_plot.png'
    bdc_est_png                 = '{!s}-event_{!s}_BDC_estimation.png'
    unc_qqplot_png              = '{!s}-uncertainty-qqplot.png'
    obs_qqplot_sorted_png       = '{!s}-mean-v-obs-sorted-qqplot.png'
    obs_qqplot_unsorted_png     = '{!s}-mean-v-obs-unsorted-plot.png'

class PanddaAnalyserFoldernames:
    pass

class PanddaAnalyserFilenames:
    initial_html            = 'pandda_initial.html'
    analyse_html            = 'pandda_analyse.html'
    analyse_site_graph      = 'pandda_analyse_sites_graph.png'
    analyse_site_graph_mult = 'pandda_analyse_sites_graph_{!s}.png'
    # Output csvs
    event_info              = 'pandda_analyse_events.csv'
    site_info               = 'pandda_analyse_sites.csv'
    dataset_info            = 'all_datasets_info.csv'
    dataset_map_info        = 'all_datasets_info_maps.csv'
    dataset_combined_info   = 'all_datasets_info_combined.csv'
    # Statistical Maps
    mean_map                = '{!s}A-mean_map.ccp4'
    medn_map                = '{!s}A-medn_map.ccp4'
    stds_map                = '{!s}A-stds_map.ccp4'
    sadj_map                = '{!s}A-sadj_map.ccp4'
    skew_map                = '{!s}A-skew_map.ccp4'
    kurt_map                = '{!s}A-kurt_map.ccp4'
    bimo_map                = '{!s}A-bimo_map.ccp4'
    # Reference Files
    reference_structure     = 'reference.pdb'
    reference_dataset       = 'reference.mtz'
    reference_on_origin     = 'reference.shifted.pdb'
    reference_symmetry      = 'reference.symmetry.pdb'
    # Pymol images
    pymol_sites_png_1       = 'pandda_analyse_sites_pymol_1.png'
    pymol_sites_png_2       = 'pandda_analyse_sites_pymol_2.png'
    pymol_sites_py          = 'pandda_analyse_sites_pymol.py'
    pymol_sites_pml         = 'pandda_analyse_sites_pymol.pml'

class PanddaInspectorFilenames:
    inspect_html            = 'pandda_inspect.html'
    inspect_site_graph      = 'pandda_inspect_sites_graph.png'
    inspect_site_graph_mult = 'pandda_inspect_sites_graph_{!s}.png'
    # output csvs
    event_info              = 'pandda_inspect_events.csv'
    site_info               = 'pandda_inspect_sites.csv'

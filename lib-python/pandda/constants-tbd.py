
DATASET_INFO_FIELDS  = [    'high_resolution',
                            'low_resolution',
                            'r_work',
                            'r_free',
                            'rmsd_to_reference',
                            'space_group',
                            'uc_a',     'uc_b',    'uc_c',
                            'uc_alpha', 'uc_beta', 'uc_gamma',
                            'uc_vol',
                            'unscaled_wilson_B',   'scaled_wilson_B',
                            'applied_b_factor_scaling',
                            'unscaled_wilson_rmsd_all',   'unscaled_wilson_rmsd_>4A',   'unscaled_wilson_rmsd_<4A',
                            'unscaled_wilson_rmsd_all_z', 'unscaled_wilson_rmsd_<4A_z', 'unscaled_wilson_rmsd_>4A_z',
                            'unscaled_wilson_ln_rmsd',    'unscaled_wilson_ln_rmsd_z',
                            'unscaled_wilson_ln_dev',     'unscaled_wilson_ln_dev_z',
                            'scaled_wilson_rmsd_all',     'scaled_wilson_rmsd_>4A',     'scaled_wilson_rmsd_<4A',
                            'scaled_wilson_rmsd_all_z',   'scaled_wilson_rmsd_<4A_z',   'scaled_wilson_rmsd_>4A_z',
                            'scaled_wilson_ln_rmsd',      'scaled_wilson_ln_rmsd_z',
                            'scaled_wilson_ln_dev',       'scaled_wilson_ln_dev_z',
                            'rejection_reason'                ]

########################################################################


class PanddaHtmlFilenames:

    initial_html            = 'pandda_initial.html'

    map_html                = 'pandda_map_{:.2f}A.html'

    analyse_html            = 'pandda_analyse.html'
    analyse_site_graph      = 'pandda_analyse_sites_graph.png'
    analyse_site_graph_mult = 'pandda_analyse_sites_graph_{!s}.png'
    pymol_sites_png_1       = 'pandda_analyse_sites_pymol_1.png'
    pymol_sites_png_2       = 'pandda_analyse_sites_pymol_2.png'
    pymol_sites_py          = 'pandda_analyse_sites_pymol.py'
    pymol_sites_pml         = 'pandda_analyse_sites_pymol.pml'

    inspect_html            = 'pandda_inspect.html'
    inspect_site_graph      = 'pandda_inspect_sites_graph.png'
    inspect_site_graph_mult = 'pandda_inspect_sites_graph_{!s}.png'

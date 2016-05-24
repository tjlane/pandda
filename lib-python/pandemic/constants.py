
DATASET_FIELDS       = [    'high_resolution',
                            'low_resolution'
                            'centre_of_mass'        ]

RESIDUE_LABEL_FIELDS = [    'resname',
                            'conformers',
                            'num_atoms'             ]

RESIDUE_PARAM_FIELDS = [    # Coordinates of the calpha
                            'calpha_x',
                            'calpha_y',
                            'calpha_z',
                            # Mean & Stdv of B-factors
                            'mean_b_all',
                            'stdv_b_all',
                            'mean_b_back',
                            'stdv_b_back',
                            'mean_b_side',
                            'stdv_b_side',
                            # Mean & Stdv of B-factors normalised internally to z-scores
                            'mean_b_z_all',
                            'stdv_b_z_all',
                            'mean_b_z_back',
                            'stdv_b_z_back',
                            'mean_b_z_side',
                            'stdv_b_z_side',
                            # Variation of the alignments at this residue
                            'var_r_mx',
                            'var_t_mx',
                            # Dihedral angle of this residue
                            'mean_dihedral',
                            'stdv_dihedral'         ]

# Individual parameters of each residue in each dataset
DST_RES_PARAM_FIELDS = [    # Mean & Stdv of B-factors
                            'mean_b_all',
                            'stdv_b_all',
                            'mean_b_back',
                            'stdv_b_back',
                            'mean_b_side',
                            'stdv_b_side',
                            # Mean & Stdv of B-factors normalised internally to z-scores
                            'mean_b_z_all',
                            'stdv_b_z_all',
                            'mean_b_z_back',
                            'stdv_b_z_back',
                            'mean_b_z_side',
                            'stdv_b_z_side',
                            # Dihedral angle of this residue
                            'dihedral'              ]

# Individual parameters of each residue in each dataset, compared to the ensemble as Z-scores
DST_RES_Z_VAL_FIELDS = [    'mean_b_all_z',
                            'stdv_b_all_z',
                            'mean_b_back_z',
                            'stdv_b_back_z',
                            'mean_b_side_z',
                            'stdv_b_side_z',
                            'mean_b_z_all_z',
                            'stdv_b_z_all_z',
                            'mean_b_z_back_z',
                            'stdv_b_z_back_z',
                            'mean_b_z_side_z',
                            'stdv_b_z_side_z',
                            'dihedral_z'            ]

class PandemicTableFields:
    all_dataset_fields = DATASET_FIELDS
    all_residue_fields = RESIDUE_LABEL_FIELDS + RESIDUE_PARAM_FIELDS
    all_dst_res_fields = DST_RES_PARAM_FIELDS + DST_RES_Z_VAL_FIELDS

########################################################################

class PandemicDatasetFilenames:
    pass

class PandemicResidueFilenames:
    pass

class PandemicAnalyserFilenames:
    pass

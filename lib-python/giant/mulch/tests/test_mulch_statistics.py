import pytest

def test_ExtractBasicXrayStatistics(five_baz2b_test_datasets_mcd):

    from giant.mulch.statistics.xray import ExtractBasicXrayStatistics

    extract = ExtractBasicXrayStatistics()

    dataset_dict = extract(
        dataset = five_baz2b_test_datasets_mcd.datasets["BAZ2BA-x434.dimple.pdb"],
        )

    print(dataset_dict)

    assert dataset_dict["r_free"] == pytest.approx(0.20442)
    assert dataset_dict["r_work"] == pytest.approx(0.17842)
    assert dataset_dict["space_group"] == "C 2 2 21"
    assert dataset_dict["unit_cell"] == (81.403, 96.831, 57.784, 90, 90, 90)

def test_ExtractWilsonStatistics(five_baz2b_test_datasets_mcd):

    ##

    reference_dataset = five_baz2b_test_datasets_mcd.datasets["BAZ2BA-x434.dimple.pdb"]
    target_dataset = five_baz2b_test_datasets_mcd.datasets["BAZ2BA-x430.dimple.pdb"]

    ##

    from giant.mulch.xray.get_data import LoadFirstValidMillerArrayFromListOfColumnOptions

    get_miller_array = LoadFirstValidMillerArrayFromListOfColumnOptions(
        structure_factor_pairs = [("2FOFCWT", "PH2FOFCWT")],
        )

    from giant.mulch.xray.scale_data import GetIsotropicMillerArrayScaler

    get_scaler = GetIsotropicMillerArrayScaler(
        reference_miller_array = get_miller_array(reference_dataset),
        )

    scaling_object = get_scaler(
        miller_array = get_miller_array(target_dataset),
        )

    ##

    from giant.mulch.statistics.xray import ExtractWilsonStatistics

    extract = ExtractWilsonStatistics()

    wilson_statistics = extract(
        dataset = target_dataset, 
        get_miller_array = get_miller_array,
        get_scaling_object = (lambda d: scaling_object),
        )

    print(wilson_statistics)

    assert wilson_statistics['wilson_b_input'] == pytest.approx(10.323167)
    assert wilson_statistics['wilson_b_scaled'] == pytest.approx(10.321833)

def test_ExtractDatasetStatistics(five_baz2b_test_datasets_mcd):

    from giant.mulch.statistics import ExtractDatasetStatistics
    from giant.mulch.statistics.xray import ExtractBasicXrayStatistics

    extract = ExtractDatasetStatistics(
        extracters = [ExtractBasicXrayStatistics()],
        )

    dataset_statistics = extract(
        mcd = five_baz2b_test_datasets_mcd,
        )

    dataframe = dataset_statistics.as_dataframe()

    print(dataframe)

    assert dataframe.loc["BAZ2BA-x434.dimple.pdb"]["r_free"] == pytest.approx(0.20442)
    assert dataframe.loc["BAZ2BA-x434.dimple.pdb"]["r_work"] == pytest.approx(0.17842)
    assert dataframe.loc["BAZ2BA-x434.dimple.pdb"]["space_group"] == "C 2 2 21"
    assert dataframe.loc["BAZ2BA-x434.dimple.pdb"]["unit_cell"] == (81.403, 96.831, 57.784, 90, 90, 90)

def test_ClassifyDatasetStatistics(five_baz2b_test_datasets_mcd):

    from giant.mulch.statistics import ExtractDatasetStatistics, ClassifyDatasetStatistics
    from giant.mulch.statistics.xray import ExtractBasicXrayStatistics
    from giant.mulch.statistics.classifiers import ZScoreClassifier

    extract = ExtractDatasetStatistics(
        extracters = [ExtractBasicXrayStatistics()],
        )

    dataset_statistics = extract(
        mcd = five_baz2b_test_datasets_mcd,
        )

    print(dataset_statistics.as_dataframe())

    ##

    classify = ClassifyDatasetStatistics(
        classifiers = [
            ZScoreClassifier(
                columns = ["r_free"],
                outlier_partition = "test_partition",
                column_suffix = "_test",
                max_z_score = 1.0,
                use_absolute_z_score = True,
                z_score_type = "normal",
                ),
            ],
        )

    dataframe = dataset_statistics.as_dataframe()
    classifications = classify(dataframe)

    print(dataframe)
    print(classifications)

    assert dataframe.loc["BAZ2BA-x434.dimple.pdb"]["r_free"] == pytest.approx(0.20442)
    assert dataframe.loc["BAZ2BA-x434.dimple.pdb"]["r_free_test"] == pytest.approx(-1.321304)

    assert sorted(classifications["test_partition"]) == ['BAZ2BA-x432.dimple.pdb', 'BAZ2BA-x434.dimple.pdb']

    ##

    classify = ClassifyDatasetStatistics(
        classifiers = [
            ZScoreClassifier(
                columns = ["r_free"],
                outlier_partition = "test_partition",
                column_suffix = "_test",
                max_z_score = 1.0,
                use_absolute_z_score = False, # one-sided
                z_score_type = "normal",
                ),
            ],
        )

    dataframe = dataset_statistics.as_dataframe()
    classifications = classify(dataframe)

    print(dataframe)
    print(classifications)
    
    assert dataframe.loc["BAZ2BA-x434.dimple.pdb"]["r_free"] == pytest.approx(0.20442)
    assert dataframe.loc["BAZ2BA-x434.dimple.pdb"]["r_free_test"] == pytest.approx(-1.321304)

    assert sorted(classifications["test_partition"]) == ['BAZ2BA-x432.dimple.pdb']

    ##

    classify = ClassifyDatasetStatistics(
        classifiers = [
            ZScoreClassifier(
                columns = ["r_free"],
                outlier_partition = "test_partition",
                column_suffix = "_test",
                max_z_score = 1.0,
                use_absolute_z_score = True,
                z_score_type = "modified", # uses median-based rejection useful for heavy tails and skew
                ),
            ],
        )

    dataframe = dataset_statistics.as_dataframe()
    classifications = classify(dataframe)

    print(dataframe)
    print(classifications)

    assert dataframe.loc["BAZ2BA-x434.dimple.pdb"]["r_free"] == pytest.approx(0.20442)
    assert dataframe.loc["BAZ2BA-x434.dimple.pdb"]["r_free_test"] == pytest.approx(-2.914804)

    assert sorted(classifications["test_partition"]) == ['BAZ2BA-x432.dimple.pdb', 'BAZ2BA-x434.dimple.pdb']

    ##

    classify = ClassifyDatasetStatistics(
        classifiers = [
            ZScoreClassifier(
                columns = ["r_free"],
                outlier_partition = "test_partition",
                column_suffix = "_test",
                max_z_score = 3.0, # should have no outliers now
                use_absolute_z_score = True,
                z_score_type = "normal",
                ),
            ],
        )

    dataframe = dataset_statistics.as_dataframe()
    classifications = classify(dataframe)

    print(dataframe)
    print(classifications)
    
    assert dataframe.loc["BAZ2BA-x434.dimple.pdb"]["r_free"] == pytest.approx(0.20442)
    assert dataframe.loc["BAZ2BA-x434.dimple.pdb"]["r_free_test"] == pytest.approx(-1.321304)

    assert sorted(classifications["test_partition"]) == []

    ##

    # missing columns

    classify = ClassifyDatasetStatistics(
        classifiers = [
            ZScoreClassifier(
                columns = ["missing_column1", "missing_column2"],
                outlier_partition = "test_partition",
                column_suffix = "_test",
                max_z_score = 3.0, # should have no outliers now
                use_absolute_z_score = True,
                z_score_type = "normal",
                ),
            ],
        )

    dataframe = dataset_statistics.as_dataframe()

    with pytest.raises(ValueError) as e:
        classifications = classify(dataframe)

    assert "'missing_column1'" in str(e.value)
    assert "'missing_column2'" in str(e.value)

    ##

def test_ExtractAndClassifyDatasetStatistics(five_baz2b_test_datasets_mcd):

    from giant.mulch.statistics import ExtractAndClassifyDatasetStatistics
    from giant.mulch.statistics.xray import ExtractBasicXrayStatistics
    from giant.mulch.statistics.classifiers import ZScoreClassifier

    extract_and_classify = ExtractAndClassifyDatasetStatistics(
        extracters = [
            ExtractBasicXrayStatistics(),
            ],
        classifiers = [
            ZScoreClassifier(
                columns = ["r_free"],
                outlier_partition = "test_partition",
                column_suffix = "_test",
                max_z_score = 1.0, 
                use_absolute_z_score = True,
                z_score_type = "normal",
                ),            
            ],
        )

    statistics_info = extract_and_classify(
        mcd = five_baz2b_test_datasets_mcd,
        )
    
    print(str(statistics_info))

    assert statistics_info.statistics.data["BAZ2BA-x434.dimple.pdb"]["r_free"] == pytest.approx(0.20442)

    assert statistics_info.dataframe.loc["BAZ2BA-x434.dimple.pdb"]["r_free"] == pytest.approx(0.20442)
    assert statistics_info.dataframe.loc["BAZ2BA-x434.dimple.pdb"]["r_free_test"] == pytest.approx(-1.321304)

    assert statistics_info.classifications["test_partition"] == ['BAZ2BA-x432.dimple.pdb', 'BAZ2BA-x434.dimple.pdb']

    ## 

def test_ExtractAndClassifyScalingStatistics(five_baz2b_test_datasets_mcd_labelled):

    ##

    reference_dataset = five_baz2b_test_datasets_mcd_labelled.datasets["BAZ2BA-x434.dimple.pdb"]

    ##

    # Prep scaling objects

    from giant.mulch.xray.get_data import LoadFirstValidMillerArrayFromListOfColumnOptions

    get_miller_array = LoadFirstValidMillerArrayFromListOfColumnOptions(
        structure_factor_pairs = [("2FOFCWT", "PH2FOFCWT")],
        )

    from giant.mulch.xray.scale_data import GetIsotropicMillerArrayScaler

    get_scaler = GetIsotropicMillerArrayScaler(
        reference_miller_array = get_miller_array(reference_dataset),
        )

    from giant.mulch.xray.scale_data import GetMultiDatasetMillerArrayScalers

    get_scalings = GetMultiDatasetMillerArrayScalers(
        get_miller_array = get_miller_array,
        get_scaler = get_scaler,
        )

    scaling_objects_dict = get_scalings(
        mcd = five_baz2b_test_datasets_mcd_labelled,
        )

    ##

    from giant.mulch.statistics import ExtractAndClassifyDatasetStatistics
    from giant.mulch.statistics.scaling import ExtractScalingStatistics, ClassifyScalingStatistics

    extract_and_classify = ExtractAndClassifyDatasetStatistics(
        extracters = [
            ExtractScalingStatistics(),
            ],
        classifiers = [
            ClassifyScalingStatistics(
                columns = None, # forces use of default columns
                outlier_partition = "test_partition",
                column_suffix = "_test",
                max_z_score = 1.0, 
                use_absolute_z_score = False, # only want positive outliers
                z_score_type = "normal",
                ),
            ],
        )

    statistics_info = extract_and_classify(
        mcd = five_baz2b_test_datasets_mcd_labelled,
        get_scaling_object = lambda d: scaling_objects_dict[d.tag],
        )

    print(statistics_info)

    # Output from extract scaling statistics 

    assert statistics_info.dataframe.loc["BAZ2BA-x430.dimple.pdb"]["b_factor_scaling"] == pytest.approx(-0.001)
    assert statistics_info.dataframe.loc["BAZ2BA-x431.dimple.pdb"]["b_factor_scaling"] == pytest.approx(-0.293)
    assert statistics_info.dataframe.loc["BAZ2BA-x432.dimple.pdb"]["b_factor_scaling"] == pytest.approx( 0.244)
    assert statistics_info.dataframe.loc["BAZ2BA-x433.dimple.pdb"]["b_factor_scaling"] == pytest.approx(-0.197)
    assert statistics_info.dataframe.loc["BAZ2BA-x434.dimple.pdb"]["b_factor_scaling"] == pytest.approx(-0.000)

    assert statistics_info.dataframe.loc["BAZ2BA-x430.dimple.pdb"]["scaled_wilson_rmsd"] == pytest.approx(4015.332)
    assert statistics_info.dataframe.loc["BAZ2BA-x431.dimple.pdb"]["scaled_wilson_rmsd"] == pytest.approx(6025.885)
    assert statistics_info.dataframe.loc["BAZ2BA-x432.dimple.pdb"]["scaled_wilson_rmsd"] == pytest.approx(4889.924)
    assert statistics_info.dataframe.loc["BAZ2BA-x433.dimple.pdb"]["scaled_wilson_rmsd"] == pytest.approx(2953.761)
    assert statistics_info.dataframe.loc["BAZ2BA-x434.dimple.pdb"]["scaled_wilson_rmsd"] == pytest.approx(0.000)

    assert statistics_info.dataframe.loc["BAZ2BA-x430.dimple.pdb"]["unscaled_wilson_rmsd"] == pytest.approx(4012.384)
    assert statistics_info.dataframe.loc["BAZ2BA-x431.dimple.pdb"]["unscaled_wilson_rmsd"] == pytest.approx(4546.437)
    assert statistics_info.dataframe.loc["BAZ2BA-x432.dimple.pdb"]["unscaled_wilson_rmsd"] == pytest.approx(7667.431)
    assert statistics_info.dataframe.loc["BAZ2BA-x433.dimple.pdb"]["unscaled_wilson_rmsd"] == pytest.approx(3695.816)
    assert statistics_info.dataframe.loc["BAZ2BA-x434.dimple.pdb"]["unscaled_wilson_rmsd"] == pytest.approx(0.000)

    # Output from classify scaling statistics 

    assert statistics_info.classifications['test_partition'] == ["BAZ2BA-x431.dimple.pdb"]

    assert statistics_info.dataframe.loc["BAZ2BA-x430.dimple.pdb"]["scaled_wilson_rmsd_test"] == pytest.approx(0.213405, abs=1e-3)
    assert statistics_info.dataframe.loc["BAZ2BA-x431.dimple.pdb"]["scaled_wilson_rmsd_test"] == pytest.approx(1.192213, abs=1e-3)
    assert statistics_info.dataframe.loc["BAZ2BA-x432.dimple.pdb"]["scaled_wilson_rmsd_test"] == pytest.approx(0.639187, abs=1e-3)
    assert statistics_info.dataframe.loc["BAZ2BA-x433.dimple.pdb"]["scaled_wilson_rmsd_test"] == pytest.approx(-0.303405, abs=1e-3)
    assert statistics_info.dataframe.loc["BAZ2BA-x434.dimple.pdb"]["scaled_wilson_rmsd_test"] == pytest.approx(-1.741400, abs=1e-3)

    assert statistics_info.dataframe.loc["BAZ2BA-x430.dimple.pdb"]["scaled_wilson_rmsd_high_test"] == pytest.approx(0.205981, abs=1e-3)
    assert statistics_info.dataframe.loc["BAZ2BA-x431.dimple.pdb"]["scaled_wilson_rmsd_high_test"] == pytest.approx(1.199517, abs=1e-3)
    assert statistics_info.dataframe.loc["BAZ2BA-x432.dimple.pdb"]["scaled_wilson_rmsd_high_test"] == pytest.approx(0.639640, abs=1e-3)
    assert statistics_info.dataframe.loc["BAZ2BA-x433.dimple.pdb"]["scaled_wilson_rmsd_high_test"] == pytest.approx(-0.309021, abs=1e-3)
    assert statistics_info.dataframe.loc["BAZ2BA-x434.dimple.pdb"]["scaled_wilson_rmsd_high_test"] == pytest.approx(-1.736116, abs=1e-3)

    assert statistics_info.dataframe.loc["BAZ2BA-x430.dimple.pdb"]["scaled_wilson_rmsd_low_test"] == pytest.approx(0.962035, abs=1e-3)
    assert statistics_info.dataframe.loc["BAZ2BA-x431.dimple.pdb"]["scaled_wilson_rmsd_low_test"] == pytest.approx(0.039969, abs=1e-3)
    assert statistics_info.dataframe.loc["BAZ2BA-x432.dimple.pdb"]["scaled_wilson_rmsd_low_test"] == pytest.approx(0.578210, abs=1e-3)
    assert statistics_info.dataframe.loc["BAZ2BA-x433.dimple.pdb"]["scaled_wilson_rmsd_low_test"] == pytest.approx(0.325694, abs=1e-3)
    assert statistics_info.dataframe.loc["BAZ2BA-x434.dimple.pdb"]["scaled_wilson_rmsd_low_test"] == pytest.approx(-1.905908, abs=1e-3)

    assert statistics_info.dataframe.loc["BAZ2BA-x430.dimple.pdb"]["scaled_wilson_ln_rmsd_test"] == pytest.approx(0.747113, abs=1e-3)
    assert statistics_info.dataframe.loc["BAZ2BA-x431.dimple.pdb"]["scaled_wilson_ln_rmsd_test"] == pytest.approx(0.747113, abs=1e-3)
    assert statistics_info.dataframe.loc["BAZ2BA-x432.dimple.pdb"]["scaled_wilson_ln_rmsd_test"] == pytest.approx(0.556523, abs=1e-3)
    assert statistics_info.dataframe.loc["BAZ2BA-x433.dimple.pdb"]["scaled_wilson_ln_rmsd_test"] == pytest.approx(-0.167719, abs=1e-3)
    assert statistics_info.dataframe.loc["BAZ2BA-x434.dimple.pdb"]["scaled_wilson_ln_rmsd_test"] == pytest.approx(-1.883029, abs=1e-3)

    assert statistics_info.dataframe.loc["BAZ2BA-x430.dimple.pdb"]["scaled_wilson_ln_dev_test"] == pytest.approx(0.559273, abs=1e-3)
    assert statistics_info.dataframe.loc["BAZ2BA-x431.dimple.pdb"]["scaled_wilson_ln_dev_test"] == pytest.approx(0.744385, abs=1e-3)
    assert statistics_info.dataframe.loc["BAZ2BA-x432.dimple.pdb"]["scaled_wilson_ln_dev_test"] == pytest.approx(0.670450, abs=1e-3)
    assert statistics_info.dataframe.loc["BAZ2BA-x433.dimple.pdb"]["scaled_wilson_ln_dev_test"] == pytest.approx(-0.055644, abs=1e-3)
    assert statistics_info.dataframe.loc["BAZ2BA-x434.dimple.pdb"]["scaled_wilson_ln_dev_test"] == pytest.approx(-1.918464, abs=1e-3)


import pytest 
import numpy as np

def test_LoadFirstValidMillerArrayFromListOfColumnOptions(five_baz2b_test_datasets):

    from giant.mulch.xray.get_data import LoadFirstValidMillerArrayFromListOfColumnOptions

    loader = LoadFirstValidMillerArrayFromListOfColumnOptions(
        structure_factor_pairs = [
            ("FWT", "PHWT"),
            ("2FOFCWT", "PH2FOFCWT"),
            ],
        )

    dataset = five_baz2b_test_datasets[0]

    ma = loader(dataset)

    assert ma.is_complex_array()

    ##

    loader = LoadFirstValidMillerArrayFromListOfColumnOptions(
        structure_factor_pairs = [
            ("NOTF", "NOTPHF"),
            ("NOTG", "NOTPHG"),
            ],
        )

    with pytest.raises(ValueError) as e:
        loader(dataset)

    assert str(e.value).startswith("No matching structure factors")

    ##

def test_GetIsotropicMillerArrayScaler(five_baz2b_test_datasets_mcd):

    reference_dataset = five_baz2b_test_datasets_mcd.datasets["BAZ2BA-x434.dimple.pdb"]
    target_dataset = five_baz2b_test_datasets_mcd.datasets["BAZ2BA-x430.dimple.pdb"]

    from giant.mulch.xray.get_data import LoadFirstValidMillerArrayFromListOfColumnOptions

    get_miller_array = LoadFirstValidMillerArrayFromListOfColumnOptions(
        structure_factor_pairs = [("2FOFCWT", "PH2FOFCWT")],
        )

    reference_miller_array = get_miller_array(reference_dataset)
    target_miller_array = get_miller_array(target_dataset)

    from giant.mulch.xray.scale_data import GetIsotropicMillerArrayScaler

    get_scaler = GetIsotropicMillerArrayScaler(
        reference_miller_array = reference_miller_array,
        )

    scale_data = get_scaler(
        miller_array = target_miller_array,
        )

    target_miller_array_scaled = scale_data(
        miller_array = target_miller_array,
        )

    assert scale_data.info()["scaling_b_factor"] == pytest.approx(-0.0013033300638198853)
    assert scale_data.info()["scaled_ln_dev"] == pytest.approx(89.94759742092617)
    assert scale_data.info()["unscaled_ln_dev"] == pytest.approx(90.15607943723163)
    assert scale_data.info()["scaled_ln_rmsd"] == pytest.approx(0.06946983579455338)
    assert scale_data.info()["unscaled_ln_rmsd"] == pytest.approx(0.06970549973307504)

    from giant.xray.data import estimate_wilson_b_factor
    wilson_b_reference = estimate_wilson_b_factor(reference_miller_array)
    wilson_b_target_unscaled = estimate_wilson_b_factor(target_miller_array)
    wilson_b_target_scaled = estimate_wilson_b_factor(target_miller_array_scaled)

    assert (wilson_b_target_scaled - wilson_b_target_unscaled) == pytest.approx(-0.0013331343517322125, abs=1e-6)

def test_GetMultiDatasetMillerArrayScalers(five_baz2b_test_datasets_mcd):
    
    reference_dataset = five_baz2b_test_datasets_mcd.datasets["BAZ2BA-x434.dimple.pdb"]

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

    scaling_objects_dict = get_scalings(five_baz2b_test_datasets_mcd)

    assert scaling_objects_dict["BAZ2BA-x434.dimple.pdb"].info()["scaling_b_factor"] == pytest.approx(0.)
    assert scaling_objects_dict["BAZ2BA-x434.dimple.pdb"].info()["scaled_ln_dev"] == pytest.approx(0.)
    assert scaling_objects_dict["BAZ2BA-x434.dimple.pdb"].info()["unscaled_ln_dev"] == pytest.approx(0.)
    assert scaling_objects_dict["BAZ2BA-x434.dimple.pdb"].info()["scaled_ln_rmsd"] == pytest.approx(0.)
    assert scaling_objects_dict["BAZ2BA-x434.dimple.pdb"].info()["unscaled_ln_rmsd"] == pytest.approx(0.)

    assert scaling_objects_dict["BAZ2BA-x432.dimple.pdb"].info()["scaling_b_factor"] == pytest.approx(0.24365846067667007)
    assert scaling_objects_dict["BAZ2BA-x432.dimple.pdb"].info()["scaled_ln_dev"] == pytest.approx(93.9842415431974)
    assert scaling_objects_dict["BAZ2BA-x432.dimple.pdb"].info()["unscaled_ln_dev"] == pytest.approx(174.4268864939686)
    assert scaling_objects_dict["BAZ2BA-x432.dimple.pdb"].info()["scaled_ln_rmsd"] == pytest.approx(0.06364252415606744)
    assert scaling_objects_dict["BAZ2BA-x432.dimple.pdb"].info()["unscaled_ln_rmsd"] == pytest.approx(0.1015309431155802)

def test_DataGetterPrepper(five_baz2b_test_datasets):

    reference_dataset = five_baz2b_test_datasets[0]
    target_dataset = five_baz2b_test_datasets[1]

    from giant.mulch.xray.get_data import LoadFirstValidMillerArrayFromListOfColumnOptions

    get_miller_array = LoadFirstValidMillerArrayFromListOfColumnOptions(
        structure_factor_pairs = [("2FOFCWT", "PH2FOFCWT")],
        )

    from giant.mulch.xray.truncate_data import CommonSetMillerArrayTruncator

    truncate_data = CommonSetMillerArrayTruncator(
        miller_arrays = list(map(get_miller_array, [reference_dataset, target_dataset])),
        )

    from giant.mulch.xray.scale_data import GetIsotropicMillerArrayScaler

    get_scaler = GetIsotropicMillerArrayScaler(
        reference_miller_array = get_miller_array(reference_dataset),
        )

    scale_data = get_scaler(
        miller_array = get_miller_array(target_dataset),
        )

    from giant.mulch.xray import DataGetterPrepper

    get_data = DataGetterPrepper(
        get_data = get_miller_array,
        scale_data = scale_data,
        truncate_data = truncate_data,
        )

    target_miller_array = get_data(target_dataset)

    from giant.xray.data import estimate_wilson_b_factor
    wilson_b_out = estimate_wilson_b_factor(target_miller_array)

    assert wilson_b_out == pytest.approx(10.366082766331772)

def test_MapGetterPrepper(five_baz2b_test_datasets_mcd):

    from giant.mulch.xray.get_data import LoadFirstValidMillerArrayFromListOfColumnOptions

    get_miller_array = LoadFirstValidMillerArrayFromListOfColumnOptions(
        structure_factor_pairs = [("2FOFCWT", "PH2FOFCWT")],
        )

    from giant.mulch.xray import MapGetterPrepper

    if True: 

        get_map = MapGetterPrepper(
            get_miller_array = get_miller_array,
            map_scaling = "none",
            map_resolution_factor = 0.5,
            )

        ed_map = get_map(
            dataset = five_baz2b_test_datasets_mcd.datasets["BAZ2BA-x434.dimple.pdb"],
            resolution = 2.0,
            )

        assert ed_map.n_real() == (90, 100, 60)

        ed_data = ed_map.real_map_unpadded()

        assert ed_data[(1,4,10)] == pytest.approx(3153.602242736506)
        assert ed_data[(1,4,10)] / np.std(ed_data) == pytest.approx(0.0315718853509) # rmsd-scaled value
        assert ed_data[(1,4,10)] / ed_map.unit_cell().volume() == pytest.approx(0.00692379942548) # volume-scaled value

    if True: 

        get_map = MapGetterPrepper(
            get_miller_array = get_miller_array,
            map_scaling = "rmsd",
            map_resolution_factor = 0.5,
            )

        ed_map = get_map(
            dataset = five_baz2b_test_datasets_mcd.datasets["BAZ2BA-x434.dimple.pdb"],
            resolution = 2.0,
            )

        assert ed_map.n_real() == (90, 100, 60)

        ed_data = ed_map.real_map_unpadded()

        assert ed_data[(1,4,10)] == pytest.approx(0.0315718853509)

    if True: 

        get_map = MapGetterPrepper(
            get_miller_array = get_miller_array,
            map_scaling = "volume",
            map_resolution_factor = 0.5,
            )

        ed_map = get_map(
            dataset = five_baz2b_test_datasets_mcd.datasets["BAZ2BA-x434.dimple.pdb"],
            resolution = 2.0,
            )

        assert ed_map.n_real() == (90, 100, 60)

        ed_data = ed_map.real_map_unpadded()

        assert ed_data[(1,4,10)] == pytest.approx(0.00692379942548)

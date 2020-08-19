import copy
import pytest
import numpy as np
from giant import processors

def test_MaskedNativeMapMaker(five_baz2b_test_datasets_mcd_labelled, tmp_path):

    from giant.mulch.reference import DefaultReferenceDataset

    input_dataset = five_baz2b_test_datasets_mcd_labelled.datasets["BAZ2BA-x434.dimple.pdb"]
    reference_dataset = DefaultReferenceDataset(
        model = copy.deepcopy(input_dataset.model),
        data = copy.deepcopy(input_dataset.data),
        )

    #

    from giant.mulch.transform.maps import GetWarpedMapLoader, GetWarpedMapGrid

    get_map_grid = GetWarpedMapGrid(
        map_grid_spacing = 2.0,
        outer_mask_radius = 5.0,
        inner_mask_radius = 1.8,
        symmetry_mask_radius = 3.0,
        mask_pdb = None,
        align_mask_to_reference = False, # Not applicable if mask_pdb not provided
        create_grid_selection_string = "resname EDO", # define analysis region
        mask_grid_selection_string = "not hetero", # define regions
        partition_grid_selection_string = "pepnames and name CA", # define partitions
        processor = None,
        )

    map_grid = get_map_grid(
        reference_dataset = reference_dataset,
        )

    #

    from giant.mulch.transform.align import AlignDatasets

    align_datasets = AlignDatasets(
        method = "local",
        alignment_kwargs = dict(
            require_hierarchies_identical = False,
            ),
        processor = None,
        )

    alignments = align_datasets(
        mcd = five_baz2b_test_datasets_mcd_labelled,
        reference_dataset = reference_dataset,
        )

    #

    # Create some map data to embed

    map_data = np.ones(
        map_grid.masks["inner"].mask_size
        )

    from giant.mulch.transform.maps.output import MaskedNativeMapMaker

#    print(str(map_grid))

    make_native_map = MaskedNativeMapMaker(
        map_grid = map_grid,
        map_mask = map_grid.masks["inner"],
        )

    #

    output_map_data = make_native_map(
        map_data = map_data,
        dataset = five_baz2b_test_datasets_mcd_labelled.datasets["BAZ2BA-x434.dimple.pdb"],
        filename = str(tmp_path / "BAZ2BA-x434.ccp4"),
        )

    assert output_map_data.mean() == pytest.approx(0.0160, abs=1e-4)
    assert output_map_data.std() == pytest.approx(0.0752, abs=1e-4)

    large_select = (output_map_data > 0.5)
    large_values = zip(*list(np.where(large_select)))

    assert len(large_values) == 510
    assert large_values[:16] == [
        (4, 18, 17),(4, 18, 18),(4, 18, 19),(4, 19, 17),
        (4, 19, 18),(4, 23, 13),(4, 24, 13),(4, 26, 17),
        (4, 27, 17),(4, 31, 12),(4, 31, 13),(4, 32, 11),
        (4, 32, 12),(4, 32, 13),(5, 18, 17),(5, 18, 18),
        ]

    assert np.array(output_map_data[large_select]).mean() == pytest.approx(0.5762, abs=1e-4)
    assert np.array(output_map_data[large_select]).std() == pytest.approx(0.0495, abs=1e-4)

    #

    output_map_data = make_native_map(
        map_data = map_data,
        dataset = five_baz2b_test_datasets_mcd_labelled.datasets["BAZ2BA-x432.dimple.pdb"],
        filename = str(tmp_path / "BAZ2BA-x432.ccp4"),
        )

    assert output_map_data.mean() == pytest.approx(0.0178, abs=1e-4)
    assert output_map_data.std() == pytest.approx(0.0842, abs=1e-4)

    large_select = (output_map_data > 0.5)
    large_values = zip(*list(np.where(large_select)))

    assert len(large_values) == 552
    assert large_values[:16] == [
        (4, 18, 17),(4, 18, 18),(4, 18, 19),(4, 19, 17),
        (4, 19, 18),(4, 23, 13),(4, 24, 13),(4, 26, 17),
        (4, 27, 17),(4, 31, 12),(4, 31, 13),(4, 32, 11),
        (4, 32, 12),(4, 32, 13),(5, 18, 17),(5, 18, 18),
        ]

    assert np.array(output_map_data[large_select]).mean() == pytest.approx(0.6409, abs=1e-4)
    assert np.array(output_map_data[large_select]).std() == pytest.approx(0.0924, abs=1e-4)

    #

def test_NativeMapMaker(five_baz2b_test_datasets_mcd_labelled, tmp_path):

    from giant.mulch.reference import DefaultReferenceDataset

    input_dataset = five_baz2b_test_datasets_mcd_labelled.datasets["BAZ2BA-x434.dimple.pdb"]
    reference_dataset = DefaultReferenceDataset(
        model = copy.deepcopy(input_dataset.model),
        data = copy.deepcopy(input_dataset.data),
        )

    #

    from giant.mulch.transform.maps import GetWarpedMapLoader, GetWarpedMapGrid

    get_map_grid = GetWarpedMapGrid(
        map_grid_spacing = 2.0,
        outer_mask_radius = 5.0,
        inner_mask_radius = 1.8,
        symmetry_mask_radius = 3.0,
        mask_pdb = None,
        align_mask_to_reference = False, # Not applicable if mask_pdb not provided
        create_grid_selection_string = "resname EDO", # define analysis region
        mask_grid_selection_string = "not hetero", # define regions
        partition_grid_selection_string = "pepnames and name CA", # define partitions
        processor = None,
        )

    map_grid = get_map_grid(
        reference_dataset = reference_dataset,
        )

    #

    from giant.mulch.transform.align import AlignDatasets

    align_datasets = AlignDatasets(
        method = "local",
        alignment_kwargs = dict(
            require_hierarchies_identical = False,
            ),
        processor = None,
        )

    alignments = align_datasets(
        mcd = five_baz2b_test_datasets_mcd_labelled,
        reference_dataset = reference_dataset,
        )

    #

    # Create some map data to embed

    map_data = np.ones(
        map_grid.grid_size_1d()
        ).reshape(
        map_grid.grid_size()
        )

    from giant.mulch.transform.maps.output import NativeMapMaker

    make_native_map = NativeMapMaker(
        map_grid = map_grid,
        sites_mask = [(10., 35., 25.),(10., 35., 30.)],
        sites_mask_radius = 2.0,
        )

    #

    output_map_data = make_native_map(
        map_data = map_data,
        dataset = five_baz2b_test_datasets_mcd_labelled.datasets["BAZ2BA-x434.dimple.pdb"],
        filename = str(tmp_path / "BAZ2BA-x434.ccp4"),
        )

    assert output_map_data.mean() == pytest.approx(0.0010, abs=1e-4)
    assert output_map_data.std() == pytest.approx(0.0282, abs=1e-4)

    large_select = (output_map_data > 0.5)
    large_values = zip(*list(np.where(large_select)))

    assert len(large_values) == 64
    assert large_values[:16] == [
        (5, 18, 12), (5, 18, 13), (5, 18, 15), (5, 18, 16),
        (5, 32, 14), (5, 32, 15), (5, 32, 17), (5, 32, 18),
        (6, 18, 12), (6, 18, 13), (6, 18, 15), (6, 18, 16),
        (6, 32, 14), (6, 32, 15), (6, 32, 17), (6, 32, 18),
        ]

    assert np.array(output_map_data[large_select]).mean() == pytest.approx(0.8502, abs=1e-4)
    assert np.array(output_map_data[large_select]).std() == pytest.approx(0.2402, abs=1e-4)

    #

    output_map_data = make_native_map(
        map_data = map_data,
        dataset = five_baz2b_test_datasets_mcd_labelled.datasets["BAZ2BA-x432.dimple.pdb"],
        filename = str(tmp_path / "BAZ2BA-x432.ccp4"),
        )

    assert output_map_data.mean() == pytest.approx(0.0010, abs=1e-4)
    assert output_map_data.std() == pytest.approx(0.0291, abs=1e-4)

    large_select = (output_map_data > 0.5)
    large_values = zip(*list(np.where(large_select)))

    assert len(large_values) == 56
    assert large_values[:16] == [
        (5, 18, 12), (5, 18, 13), (5, 18, 15), (5, 18, 16),
        (5, 32, 14), (5, 32, 15), (5, 32, 17), (5, 32, 18),
        (6, 18, 12), (6, 18, 13), (6, 18, 15), (6, 18, 16),
        (6, 32, 14), (6, 32, 15), (6, 32, 17), (6, 32, 18),
        ]

    assert np.array(output_map_data[large_select]).mean() == pytest.approx(0.9411, abs=1e-4)
    assert np.array(output_map_data[large_select]).std() == pytest.approx(0.1887, abs=1e-4)

    #

import copy
import pytest
import numpy as np
from giant import processors

def test_GetWarpedMapLoader(five_baz2b_test_datasets_mcd_labelled):

    from giant.mulch.reference import DefaultReferenceDataset

    input_dataset = five_baz2b_test_datasets_mcd_labelled.datasets["BAZ2BA-x434.dimple.pdb"]
    reference_dataset = DefaultReferenceDataset(
        model = copy.deepcopy(input_dataset.model),
        data = copy.deepcopy(input_dataset.data),
        )

    #

    from giant.mulch.xray.get_data import LoadFirstValidMillerArrayFromListOfColumnOptions

    get_miller_array = LoadFirstValidMillerArrayFromListOfColumnOptions(
        structure_factor_pairs = [("2FOFCWT", "PH2FOFCWT")],
        )
    
    from giant.mulch.xray import MapGetterPrepper

    reference_map_getter = MapGetterPrepper(
        get_miller_array = get_miller_array, 
        map_scaling = "rmsd",
        map_resolution_factor = 0.5,
        )

    dataset_map_getters = {
        dkey: MapGetterPrepper(
            get_miller_array = get_miller_array, 
            map_scaling = "none", # to test other map scaling
            map_resolution_factor = 0.5,
            )
        for dkey, d in five_baz2b_test_datasets_mcd_labelled.datasets.items()
    }

    #

    from giant.mulch.transform.maps import GetWarpedMapLoader, GetWarpedMapGrid
    get_warped_loader = GetWarpedMapLoader(
        get_map_grid = GetWarpedMapGrid(
            map_grid_spacing = 2.0,
            outer_mask_radius = 10.0, 
            inner_mask_radius = 1.8,
            symmetry_mask_radius = 3.0,
            mask_pdb = None,
            align_mask_to_reference = False, # Not applicable if mask_pdb not provided
            create_grid_selection_string = "resname EDO", # define analysis region
            mask_grid_selection_string = "not hetero", # define regions
            partition_grid_selection_string = "pepnames and name CA", # define partitions
            processor = None,
            ),
        processor = processors.ProcessorJoblib(n_cpus=1),
        )

    load_maps = get_warped_loader(
        reference_dataset = reference_dataset,
        reference_map_getter = reference_map_getter,
        dataset_map_getters = dataset_map_getters,
        )    

    #

    # Loading maps requires aligned datasets! 

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

    map_manager = load_maps(
        datasets = five_baz2b_test_datasets_mcd_labelled.datasets,
        map_resolution = 2.0,
        )

    # 

    assert hasattr(map_manager, 'map_data')

    map_data = map_manager.as_dict(
        delete_internal_references = False,
        )

    assert hasattr(map_manager, 'map_data') # should still have internal reference

    map_data = map_manager.as_dict(
        delete_internal_references = True,
        )

    assert not hasattr(map_manager, 'map_data') # should have deleted internal reference

    #

    # The maps should come out as approximately rmsd-scaled 
    # because the reference map is rmsd-scaled! 

    assert map_data["BAZ2BA-x430.dimple.pdb"].mean() == pytest.approx(0.06928193304074237)
    assert map_data["BAZ2BA-x430.dimple.pdb"].std() == pytest.approx(1.1109177550460645)

    assert map_data["BAZ2BA-x431.dimple.pdb"].mean() == pytest.approx(0.04372304816627063)
    assert map_data["BAZ2BA-x431.dimple.pdb"].std() == pytest.approx(1.1218759358124852)

    assert map_data["BAZ2BA-x432.dimple.pdb"].mean() == pytest.approx(0.04845977984781397)
    assert map_data["BAZ2BA-x432.dimple.pdb"].std() == pytest.approx(1.1079459493273707)

    assert map_data["BAZ2BA-x433.dimple.pdb"].mean() == pytest.approx(0.05004279909365189)
    assert map_data["BAZ2BA-x433.dimple.pdb"].std() == pytest.approx(1.1157386526104653)

    assert map_data["BAZ2BA-x434.dimple.pdb"].mean() == pytest.approx(0.0786052458163574)
    assert map_data["BAZ2BA-x434.dimple.pdb"].std() == pytest.approx(1.1305115800134509)

    #

    reference_map = load_maps.get_reference_map(
        dataset = load_maps.reference_dataset,
        resolution = 2.0,
        )

    warped_reference_map = load_maps.warp_reference_map(
        dataset = load_maps.reference_dataset,
        fft_map = reference_map,
        )

    correlation_coeffs = [
        np.corrcoef([warped_reference_map, map_data[dkey]])[0,1]
        for dkey in [
            "BAZ2BA-x430.dimple.pdb",
            "BAZ2BA-x431.dimple.pdb",
            "BAZ2BA-x432.dimple.pdb",
            "BAZ2BA-x433.dimple.pdb",
            "BAZ2BA-x434.dimple.pdb",
        ]
    ]

    assert correlation_coeffs == pytest.approx([
        0.9693018, 0.9673743, 0.95408347, 0.96689988, 1.0,
        ])

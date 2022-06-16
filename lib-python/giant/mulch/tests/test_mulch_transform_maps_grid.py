import copy
import pytest
import numpy as np

def test_MapGrid():
    
    from giant.mulch.transform.maps.grid import MapGrid

    ##

    map_grid = MapGrid(
        grid_spacing = 1.2, 
        cart_origin = (12., -5., 17.),
        cart_approx_max = (14., 0., 20.),
        )

    assert map_grid.grid_spacing() == pytest.approx(1.2)
    assert map_grid.grid_size() == (3, 6, 4)
    assert map_grid.grid_size_1d() == 3*6*4
    assert map_grid.cart_size() == pytest.approx((3.6, 7.2, 4.8))
    assert map_grid.cart_origin() == (12., -5., 17.)

    grid_points = map_grid.grid_points()

    assert grid_points.shape == (map_grid.grid_size_1d(), 3)
    assert list(map(tuple, grid_points[:10])) == [
        (0,0,0),(0,0,1),(0,0,2),(0,0,3),(0,1,0),
        (0,1,1),(0,1,2),(0,1,3),(0,2,0),(0,2,1),
        ]
    assert tuple(grid_points[10]) == (0, 2, 2)
    assert tuple(grid_points[42]) == (1, 4, 2)

    # 

    cart_points = map_grid.cart_points(origin_shift=False)
    assert tuple(cart_points[10]) == pytest.approx((0., 2.4, 2.4))
    assert tuple(cart_points[42]) == pytest.approx((1.2, 4.8, 2.4))

    mapped_points = map_grid.grid2cart([(0, 2, 2),(1, 4, 2)], origin_shift=False)
    assert tuple(mapped_points[0]) == pytest.approx((0., 2.4, 2.4))
    assert tuple(mapped_points[1]) == pytest.approx((1.2, 4.8, 2.4))

    mapped_points = map_grid.cart2grid([(0., 2.4, 2.4),(1.2, 4.8, 2.4)], origin_shift=False)
    assert tuple(mapped_points[0]) == pytest.approx((0, 2, 2))
    assert tuple(mapped_points[1]) == pytest.approx((1, 4, 2))

    #

    cart_points = map_grid.cart_points(origin_shift=True)
    assert tuple(cart_points[10]) == pytest.approx((12., -2.6, 19.4))
    assert tuple(cart_points[42]) == pytest.approx((13.2, -0.2, 19.4))

    mapped_points = map_grid.grid2cart([(0, 2, 2),(1, 4, 2)], origin_shift=True)
    assert tuple(mapped_points[0]) == pytest.approx((12., -2.6, 19.4))
    assert tuple(mapped_points[1]) == pytest.approx((13.2, -0.2, 19.4))

    mapped_points = map_grid.cart2grid([(12., -2.6, 19.4),(13.2, -0.2, 19.4)], origin_shift=True)
    assert tuple(mapped_points[0]) == pytest.approx((0, 2, 2))
    assert tuple(mapped_points[1]) == pytest.approx((1, 4, 2))

    #

    unit_cell = map_grid.unit_cell()
    assert unit_cell.parameters() == pytest.approx((3.6, 7.2, 4.8, 90., 90., 90.))

    ##

    with pytest.raises(ValueError) as e: 

        map_grid = MapGrid(
            grid_spacing = 1.2, 
            cart_origin = (12., -5., 17.),
            cart_approx_max = (14., -6., 20.), # negative grid size
            )

    #

def test_GetWarpedMapGrid(five_baz2b_test_datasets_mcd):

    from giant.mulch.reference import DefaultReferenceDataset

    input_dataset = five_baz2b_test_datasets_mcd.datasets["BAZ2BA-x434.dimple.pdb"]
    reference_dataset = DefaultReferenceDataset(
        model = copy.deepcopy(input_dataset.model),
        data = copy.deepcopy(input_dataset.data),
        )

    from giant.mulch.transform.maps.grid import GetWarpedMapGrid

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
        )

    map_grid = get_map_grid(
        reference_dataset = reference_dataset,
        )

    assert map_grid.grid_size_1d() == 2028
    assert map_grid.grid_size() == (13, 13, 12)
    assert map_grid.cart_origin() == pytest.approx((3.811, 29.831, 19.88), abs=1e-3)
    assert map_grid.masks["outer"].mask_size == 1986
    assert map_grid.masks["inner"].mask_size == 639
    assert map_grid.masks["symmetry"].mask_size == 500
    assert map_grid.masks["total"].mask_size == 903
    assert len(map_grid.partition.partition_sites) == 116
    assert list(map_grid.partition.query_by_grid_indices(list(range(10)))) == [
        35, 35, 35, 34, 31, 31, 31, 31, 28, 28,
        ]
    print(list(np.where(map_grid.partition.embed().flatten() == -1)[0]))
    assert list(np.where(map_grid.partition.embed().flatten() == -1)[0]) == [
        103, 104, 115, 116, 117, 118, 119, 127, 128, 129, 130, 131, 139, 
        140, 141, 142, 143, 151, 152, 153, 154, 155, 272, 273, 284, 285, 
        286, 287, 296, 297, 298, 299, 308, 309, 310, 311, 464, 465, 466, 
        467, 2023, 2024,
        ]

    ### Test without mask string

    get_map_grid = GetWarpedMapGrid(
        map_grid_spacing = 2.0,
        outer_mask_radius = 10.0, 
        inner_mask_radius = 1.8,
        symmetry_mask_radius = 3.0,
        mask_pdb = None,
        align_mask_to_reference = False, # Not applicable if mask_pdb not provided
        create_grid_selection_string = "not hetero", # define analysis region
        mask_grid_selection_string = "not hetero", # define regions
        partition_grid_selection_string = "pepnames and name CA", # define partitions
        processor = None,
        )

    map_grid = get_map_grid(
        reference_dataset = reference_dataset,
        )

def test_GridMask():

    from giant.mulch.transform.maps.grid.mask import GridMask

    grid_size = (3, 6, 4)
    bool_mask = np.zeros(shape=grid_size, dtype=bool)

    assert np.product(grid_size) == 72
    assert bool_mask.shape == grid_size

    #

    grid_mask = GridMask(
        grid_size = grid_size,
        mask_indices = list(range(5)),
        mask_binary = None,
        )

    assert grid_mask.mask_size == 5
    assert grid_mask.get_mask_binary().shape == (72,) # should be flattened
    assert grid_mask.get_mask_indices().shape == (5,)

    assert list(grid_mask.get_mask_binary()) == [True]*5 + [False]*67
    assert list(grid_mask.get_mask_indices()) == list(range(5))
    assert list(map(tuple, grid_mask.get_mask_points())) == [
        (0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 0, 3), (0, 1, 0),
        ]

    #

    grid_mask = GridMask(
        grid_size = grid_size,
        mask_indices = [2, 20, 65],
        mask_binary = None,
        )

    assert grid_mask.mask_size == 3
    assert list(grid_mask.get_mask_indices()) == [2, 20, 65]
    assert list(np.where(grid_mask.get_mask_binary())[0]) == [2, 20, 65]

    assert list(map(tuple, grid_mask.get_mask_points())) == [
        (0, 0, 2), (0, 5, 0), (2, 4, 1),
        ]

    embedded_data = grid_mask.embed_data(
        data = list(range(1,4)),
        )
    
    assert embedded_data.shape == grid_size

    assert list(zip(*list(map(tuple, np.where(embedded_data))))) == [
        (0, 0, 2), (0, 5, 0), (2, 4, 1),
        ]
    assert embedded_data[(0, 0, 2)] == 1
    assert embedded_data[(0, 5, 0)] == 2
    assert embedded_data[(2, 4, 1)] == 3

    assert list(np.where(embedded_data.flatten())[0]) == [2, 20, 65]
    assert list((embedded_data != 0).flatten()) == list(grid_mask.get_mask_binary())

    assert embedded_data.flatten()[2] == 1
    assert embedded_data.flatten()[20] == 2
    assert embedded_data.flatten()[65] == 3

    with pytest.raises(ValueError) as e: 

        grid_mask.embed_data(
            data = list(range(4)),
            )

    ##

    if True: 

        with pytest.raises(ValueError) as e:

            grid_mask = GridMask(
                grid_size = grid_size,
                mask_indices = None,
                mask_binary = None,
                )

        with pytest.raises(ValueError) as e:

            grid_mask = GridMask(
                grid_size = grid_size,
                mask_indices = [72],
                mask_binary = np.ones(72, dtype=bool),
                )

        with pytest.raises(ValueError) as e:

            grid_mask = GridMask(
                grid_size = grid_size,
                mask_indices = [72],
                mask_binary = None,
                )

        with pytest.raises(ValueError) as e:

            grid_mask = GridMask(
                grid_size = grid_size,
                mask_indices = [-1],
                mask_binary = None,
                )

        with pytest.raises(ValueError) as e:

            grid_mask = GridMask(
                grid_size = grid_size,
                mask_indices = None,
                mask_binary = np.ones(73, dtype=bool),
                )

    ##

def test_GridIndexers():

    from giant.mulch.transform.maps.grid import \
        GridIndexer, GridInverter, MultiPointGridIndexer, MultiPointGridInverter

    grid_gp2idx = GridIndexer(
        grid_size = (3, 6, 4),
        )

    assert grid_gp2idx((0,0,0)) == 0
    assert grid_gp2idx((2,5,3)) == 71
    assert grid_gp2idx((0,0,3)) == 3
    assert grid_gp2idx((0,5,0)) == 20
    assert grid_gp2idx((2,0,0)) == 48

    grid_idx2gp = GridInverter(
        grid_size = (3, 6, 4),
        )

    assert grid_idx2gp(0) == (0,0,0)
    assert grid_idx2gp(71) == (2,5,3)
    assert grid_idx2gp(3) == (0,0,3)
    assert grid_idx2gp(20) == (0,5,0)
    assert grid_idx2gp(48) == (2,0,0)

    grid_gps2idx = MultiPointGridIndexer(
        grid_size = (3, 6, 4),
        )

    assert list(
        grid_gps2idx([(0,0,0),(2,5,3),(0,0,3),(0,5,0),(2,0,0)])
        ) == [0, 71, 3, 20, 48]

    grid_idx2gps = MultiPointGridInverter(
        grid_size = (3, 6, 4),
        )

    assert list(map(tuple, grid_idx2gps([0, 71, 3, 20, 48]))) == [
        (0,0,0),(2,5,3),(0,0,3),(0,5,0),(2,0,0),
        ]

def test_GetSitesMask():

    from giant.mulch.transform.maps.grid import MapGrid
    from giant.mulch.transform.maps.grid.mask import GetSitesMask

    map_grid = MapGrid(
        grid_spacing = 1.2, 
        cart_origin = (12., -5., 17.),
        cart_approx_max = (14., 0., 20.),
        )
    assert map_grid.grid_size() == (3, 6, 4)

    get_map_mask = GetSitesMask.from_map_grid(
        map_grid = map_grid,
        mask_dist = 1.5, # should 3x3 square around a grid point
        )

    ##
    
    cart_point = map_grid.grid2cart([(1,1,1)])[0]

    assert tuple(cart_point) == pytest.approx((13.2, -3.8, 18.2))

    map_mask = get_map_mask(
        sites_cart = [cart_point], 
        )

    assert tuple(map_mask.get_mask_indices()) == (5, 25, 28, 29, 30, 33, 53)
    assert list(map(tuple, map_mask.get_mask_points())) == [
        (0, 1, 1), (1, 0, 1), (1, 1, 0), 
        (1, 1, 1), 
        (1, 1, 2), (1, 2, 1), (2, 1, 1),
        ]

    #

    map_mask = get_map_mask(
        sites_cart = [(13.8, -3.2, 18.8)], 
        )

    assert tuple(map_mask.get_mask_indices()) == (29, 30, 33, 34, 53, 54, 57, 58)
    assert list(map(tuple, map_mask.get_mask_points())) == [
        (1, 1, 1), (1, 1, 2), (1, 2, 1), (1, 2, 2), 
        (2, 1, 1), (2, 1, 2), (2, 2, 1), (2, 2, 2),
        ]

    #

    map_mask = get_map_mask(
        sites_cart = [
            (13.2, -3.8, 18.2),
            (13.8, -3.2, 18.8),
            ], 
        )

    assert tuple(map_mask.get_mask_indices()) == (5, 25, 28, 29, 30, 33, 34, 53, 54, 57, 58)
    assert list(map(tuple, map_mask.get_mask_points())) == [
        (0, 1, 1), (1, 0, 1), (1, 1, 0), 
        (1, 1, 1), 
        (1, 1, 2), (1, 2, 1), (1, 2, 2), 
        (2, 1, 1), (2, 1, 2), (2, 2, 1), 
        (2, 2, 2),
        ]

    ##

    # test periodicity

    cart_point = map_grid.grid2cart([(0,0,0)])[0]

    assert tuple(cart_point) == pytest.approx((12., -5., 17.))

    map_mask = get_map_mask(
        sites_cart = [cart_point], 
        )

    assert tuple(map_mask.get_mask_indices()) == (0, 1, 3, 4, 20, 24, 48)
    assert list(map(tuple, map_mask.get_mask_points())) == [
        (0, 0, 0), (0, 0, 1), 
        (0, 0, 3), 
        (0, 1, 0), 
        (0, 5, 0), 
        (1, 0, 0), (2, 0, 0),
        ]

    ##

def test_GetNonPeriodicSitesMask():

    from giant.mulch.transform.maps.grid import MapGrid
    from giant.mulch.transform.maps.grid.mask import GetNonPeriodicSitesMask

    map_grid = MapGrid(
        grid_spacing = 1.2, 
        cart_origin = (12., -5., 17.),
        cart_approx_max = (14., 0., 20.),
        )
    assert map_grid.grid_size() == (3, 6, 4)

    get_map_mask = GetNonPeriodicSitesMask.from_map_grid(
        map_grid = map_grid,
        mask_dist = 1.5, # should 3x3 square around a grid point
        )

    #####

    cart_point = map_grid.grid2cart([(0,0,0)])[0]

    assert tuple(cart_point) == pytest.approx((12., -5., 17.))

    map_mask = get_map_mask(
        sites_cart = [cart_point], 
        )

    assert tuple(map_mask.get_mask_indices()) == (0, 1, 4, 24)
    assert list(map(tuple, map_mask.get_mask_points())) == [
        (0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0),
        ]

    # 

    cart_point = map_grid.grid2cart([(0,0,-1)])[0]

    map_mask = get_map_mask(
        sites_cart = [cart_point], 
        )

    assert tuple(map_mask.get_mask_indices()) == (0,)
    assert list(map(tuple, map_mask.get_mask_points())) == [
        (0, 0, 0),
        ]

    #

    cart_point = map_grid.grid2cart([(3,6,4)])[0]

    map_mask = get_map_mask(
        sites_cart = [cart_point], 
        )

    assert tuple(map_mask.get_mask_indices()) == tuple()
    assert list(map(tuple, map_mask.get_mask_points())) == []

    #

    cart_point = map_grid.grid2cart([(3,5,3)])[0]

    map_mask = get_map_mask(
        sites_cart = [cart_point], 
        )

    assert tuple(map_mask.get_mask_indices()) == (71,)
    assert list(map(tuple, map_mask.get_mask_points())) == [
        (2, 5, 3),
        ]

    #

    cart_points = map_grid.grid2cart(
        [(-1, 3, 3), (2, 5, 1)]
        )

    map_mask = get_map_mask(
        sites_cart = cart_points, 
        )

    assert tuple(map_mask.get_mask_indices()) == (15, 45, 65, 68, 69, 70)
    assert list(map(tuple, map_mask.get_mask_points())) == [
        (0, 3, 3), (1, 5, 1), (2, 4, 1), 
        (2, 5, 0), (2, 5, 1), (2, 5, 2),
        ]
    
    #####

    # test all edges

    cart_points = map_grid.grid2cart([
        (-1, 3, 2), (+3, 3, 2), 
        (2, -1, 2), (2, +6, 2), 
        (2, 3, -1), (2, 3, +4)
        ])

    map_mask = get_map_mask(
        sites_cart = cart_points, 
        )

    assert tuple(map_mask.get_mask_indices()) == (14, 50, 60, 62, 63, 70)
    assert list(map(tuple, map_mask.get_mask_points())) == [
        (0, 3, 2), (2, 0, 2), (2, 3, 0), 
        (2, 3, 2), (2, 3, 3), (2, 5, 2)
        ]

    # edges and 1 central point

    cart_points = map_grid.grid2cart([
        (1, 1, 1), 
        (-1, 3, 2), (+3, 3, 2), 
        (2, -1, 2), (2, +6, 2), 
        (2, 3, -1), (2, 3, +4),
        ])

    map_mask = get_map_mask(
        sites_cart = cart_points, 
        )

    assert tuple(map_mask.get_mask_indices()) == (5, 14, 25, 28, 29, 30, 33, 50, 53, 60, 62, 63, 70)
    assert list(map(tuple, map_mask.get_mask_points())) == [
        (0, 1, 1), (0, 3, 2), (1, 0, 1), (1, 1, 0), (1, 1, 1), 
        (1, 1, 2), (1, 2, 1), (2, 0, 2), (2, 1, 1), (2, 3, 0), 
        (2, 3, 2), (2, 3, 3), (2, 5, 2),
        ]

    # all edges just out of range 

    cart_points = map_grid.grid2cart([
        (-2, 3, 2), (+4, 3, 2), 
        (2, -2, 2), (2, +7, 2), 
        (2, 3, -2), (2, 3, +5),
        ])

    map_mask = get_map_mask(
        sites_cart = cart_points, 
        )

    assert tuple(map_mask.get_mask_indices()) == tuple()
    assert list(map(tuple, map_mask.get_mask_points())) == []

    # all edges out of range and one central point

    cart_points = map_grid.grid2cart([
        (1, 1, 1), 
        (-2, 3, 2), (+4, 3, 2), 
        (2, -2, 2), (2, +7, 2), 
        (2, 3, -2), (2, 3, +5),
        ])

    map_mask = get_map_mask(
        sites_cart = cart_points, 
        )

    print(tuple(map_mask.get_mask_indices()))
    print(list(map(tuple, map_mask.get_mask_points())))

    assert tuple(map_mask.get_mask_indices()) == (5, 25, 28, 29, 30, 33, 53)
    assert list(map(tuple, map_mask.get_mask_points())) == [
        (0, 1, 1), (1, 0, 1), (1, 1, 0), 
        (1, 1, 1), 
        (1, 1, 2), (1, 2, 1), (2, 1, 1),
        ]

    #####

    # now select points that're way off the grid

    cart_point = map_grid.grid2cart([(-4,-12,-5)])[0]

    map_mask = get_map_mask(
        sites_cart = [cart_point], 
        )

    assert tuple(map_mask.get_mask_indices()) == tuple()
    assert list(map(tuple, map_mask.get_mask_points())) == []

    #

    cart_point = map_grid.grid2cart([(4,12,-5)])[0]

    map_mask = get_map_mask(
        sites_cart = [cart_point], 
        )

    assert tuple(map_mask.get_mask_indices()) == tuple()
    assert list(map(tuple, map_mask.get_mask_points())) == []

    ##

def test_compound_grid_masks():

    from giant.mulch.transform.maps.grid import MapGrid
    from giant.mulch.transform.maps.grid.mask import GetSitesMask, compound_grid_masks

    map_grid = MapGrid(
        grid_spacing = 1.2, 
        cart_origin = (0., 0., 0.),
        cart_approx_max = (4.8, 6.0, 7.2),
        )
    assert map_grid.grid_size() == (5, 6, 7)

    get_map_mask = GetSitesMask.from_map_grid(
        map_grid = map_grid,
        mask_dist = 1.5, # should 3x3 square around a grid point
        )

    mask_1 = get_map_mask(map_grid.grid2cart([(1,1,0)]))
    mask_2 = get_map_mask(map_grid.grid2cart([(1,1,1)]))
    mask_3 = get_map_mask(map_grid.grid2cart([(1,1,2)]))

    #

    combined = compound_grid_masks(
        positive_masks = [mask_1, mask_2, mask_3],
        negative_masks = [],
        )

    assert list(map(tuple, combined.get_mask_points())) == [
        (0, 1, 0), (0, 1, 1), (0, 1, 2), 
        (1, 0, 0), (1, 0, 1), (1, 0, 2), 
        (1, 1, 0), 
        (1, 1, 1), (1, 1, 2), (1, 1, 3), 
        (1, 1, 6), 
        (1, 2, 0), (1, 2, 1), (1, 2, 2), 
        (2, 1, 0), (2, 1, 1), (2, 1, 2),
        ]

    #

    with pytest.raises(ValueError) as e:
        combined = compound_grid_masks(
            positive_masks = [],
            negative_masks = [mask_1, mask_2, mask_3],
        )

    # 

    combined = compound_grid_masks(
        positive_masks = [mask_1, mask_2],
        negative_masks = [mask_3],
        )

    assert list(map(tuple, combined.get_mask_points())) == [
        (0, 1, 0), (0, 1, 1), (1, 0, 0), 
        (1, 0, 1), (1, 1, 0), 
        (1, 1, 6), 
        (1, 2, 0), (1, 2, 1), (2, 1, 0), 
        (2, 1, 1),
        ]

    # 

    combined = compound_grid_masks(
        positive_masks = [mask_1, mask_3],
        negative_masks = [mask_2],
        )

    print(list(map(tuple, combined.get_mask_points())))

    assert list(map(tuple, combined.get_mask_points())) == [
        (0, 1, 0), (0, 1, 2), 
        (1, 0, 0), (1, 0, 2), 
        (1, 1, 3), 
        (1, 1, 6), 
        (1, 2, 0), (1, 2, 2), 
        (2, 1, 0), (2, 1, 2),
        ]

    # 

def test_MakeVoronoiGridPartition(five_baz2b_test_datasets_mcd):

    reference_dataset = five_baz2b_test_datasets_mcd.datasets["BAZ2BA-x434.dimple.pdb"]
     
    from giant.mulch.transform.maps.grid import MapGrid

    map_grid = MapGrid(
        grid_spacing = 2.0, 
        cart_origin = (12., 40., 28.), 
        cart_approx_max = (16., 44., 32.),
        )

    assert map_grid.grid_size_1d() == 27

    from giant.mulch.transform.maps.grid.mask import GridMask

    map_grid.masks["test"] = GridMask(
        grid_size = map_grid.grid_size(),
        mask_indices = list(range(10)),
        mask_binary = None,
        )

    from giant.mulch.transform.maps.grid.partition import MakeVoronoiGridPartition

    make_voronoi_partition = MakeVoronoiGridPartition()

    ##

    # Test independent class functions 

    # points should be assigned to themselves! 
    mappings = make_voronoi_partition.get_mappings(
        ref_sites = map_grid.grid_points(),
        query_sites = map_grid.grid_points(),
        )
    assert list(mappings) == pytest.approx(list(range(map_grid.grid_size_1d())))

    ##

    grid_partition_1 = make_voronoi_partition(
        map_grid = map_grid, 
        dataset = reference_dataset,
        selection_string = "pepnames and name CA",
        mask_name = None,
        )

    sites, sites_counts = np.unique(grid_partition_1.partition_mappings, return_counts=True)

    assert list(sites) == [32, 37, 42, 45, 85, 88, 89, 95]
    assert list(sites_counts) == [2, 4, 8, 2, 2, 3, 3, 3]

    # Check that all of the points are assigned to the right place

    this_partition = grid_partition_1
    for i in [32, 37, 42, 45, 85, 88, 89, 95]:
        # Get a site
        site_xyz = this_partition.partition_sites[i]
        # Get all points assigned to this site
        partition_xyzs_grid = list(zip(*np.where(map_grid.embed_data(this_partition.partition_mappings) == i)))
        partition_xyzs_cart = map_grid.grid2cart(partition_xyzs_grid)
        # Check that the closest site is the assigned one
        for xyz in partition_xyzs_cart:
            diff_xyzs_cart = (this_partition.partition_sites - xyz)
            dists = np.sqrt(np.power(diff_xyzs_cart, 2).sum(axis=1)) # TODO: update to use giant.common.maths function
            assert np.where(dists == min(dists))[0] == i

    #

    grid_partition_2 = make_voronoi_partition(
        map_grid = map_grid, 
        dataset = reference_dataset,
        selection_string = "pepnames and name CA",
        mask_name = "test",
        )
    
    sites, sites_counts = np.unique(grid_partition_2.partition_mappings, return_counts=True)

    assert list(sites) == [32, 37, 42, 95]
    assert list(sites_counts) == [1, 2, 5, 2]

    # Check that all of the points are assigned to the right place

    this_partition = grid_partition_2
    for i in [32, 37, 42, 95]:
        # Get a site
        site_xyz = this_partition.partition_sites[i]
        # Get all points assigned to this site
        partition_xyzs_grid = list(zip(*np.where(map_grid.masks["test"].embed_data(this_partition.partition_mappings) == i)))
        partition_xyzs_cart = map_grid.grid2cart(partition_xyzs_grid)
        # Check that the closest site is the assigned one
        for xyz in partition_xyzs_cart:
            diff_xyzs_cart = (this_partition.partition_sites - xyz)
            dists = np.sqrt(np.power(diff_xyzs_cart, 2).sum(axis=1)) # TODO: update to use giant.common.maths function
            assert np.where(dists == min(dists))[0] == i

    #

    assert list(grid_partition_2.partition_mappings) == list(grid_partition_1.partition_mappings[:10])

    ##


    

from pytest import approx, raises

import numpy

comparison_matrix = numpy.array([
    [  6.975,   4.839,   4.482,   5.798,   5.992,   4.286,   1.836,   3.056,   1.925,   5.494,   3.616],
    [  3.96 ,   8.603,   5.964,  12.111,  10.036,   8.298,   2.801,   4.776,   2.39 ,   9.5  ,   6.306],
    [  4.817,   6.118,   7.741,  10.55 ,   8.281,   6.028,   3.222,   6.589,   2.443,  10.898,   6.035],
    [  4.149,   3.637,   6.091,  14.135,   9.968,   6.665,   2.775,   3.769,   2.365,   7.958,   4.686],
    [  5.897,   6.074,   5.27 ,   8.859,  11.313,   7.234,   2.222,   4.789,   1.993,   7.986,   5.191],
    [  3.344,   6.091,   4.372,   9.341,   6.47 ,   9.856,   3.106,   6.489,   2.245,   7.474,   4.699],
    [  2.933,   3.395,   3.426,   3.582,   3.366,   2.609,   3.239,   3.183,   1.808,   4.84 ,   3.668],
    [  3.117,   3.879,   4.93 ,   8.546,   6.898,   4.94 ,   2.457,   8.437,   2.031,   8.117,   5.174],
    [  1.823,   2.93 ,   2.637,   4.44 ,   4.185,   2.422,   1.377,   1.845,   2.489,   3.788,   2.381],
    [  4.874,   5.889,   7.363,  10.054,   7.654,   4.988,   3.231,   6.782,   2.217,  13.373,   7.022],
    [  3.767,   5.759,   5.404,   9.96 ,   8.593,   3.699,   2.374,   4.901,   2.44 ,  10.045,   8.315],
    ])

connectivity_matrix = numpy.array([
    [1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1],
    [1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1],
    [0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1],
    [0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0],
    [0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0],
    [0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1],
    [1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1],
    ], dtype=bool)


def test_identify_modules_exceptions():

    from pandemic.adp.echt.analysis.cluster_tls_groups.identify_modules import IdentifyModules

    for metric in ['similarity', 'distance']:
        with raises(AssertionError) as e:    
            im = IdentifyModules(
                threshold_start = -1.0,
                threshold_delta = 1.0,
                comparison_metric = metric,
                )    
        assert str(e.value) == "start_threshold must be positive: {} thresholds must be positive".format(metric)

    with raises(AssertionError) as e:    
        im = IdentifyModules(
            threshold_start = 5.0,
            threshold_delta = -1.0,
            comparison_metric = 'similarity',
            )    
    assert str(e.value) == "delta_threshold must be positive: search is from small thresholds to large thresholds"

    with raises(AssertionError) as e:    
        im = IdentifyModules(
            threshold_start = 5.0,
            threshold_delta = 1.0,
            comparison_metric = 'distance',
            )    
    assert str(e.value) == "delta_threshold must be negative: search is from large thresholds to small thresholds"

    with raises(AssertionError) as e:    
        im = IdentifyModules(
            threshold_start = -1.0,
            threshold_delta = 1.0,
            comparison_metric = 'x',
            )    
    assert str(e.value) == "Invalid comparison_metric (x) provided"        

def test_identify_modules():

    from pandemic.adp.echt.analysis.cluster_tls_groups.identify_modules import IdentifyModules

    im = IdentifyModules(
        threshold_start = 1.0,
        threshold_delta = 1.0,
        comparison_metric = 'similarity',
        )

    result = im(
        connectivity = connectivity_matrix,
        comparison_matrix = comparison_matrix,
        filter_complete_sets = False,
        )

    assert result.cumulative_hierarchy == [
        [(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)], 
        [(1, 2, 3, 6, 7, 9, 10)], 
        [(1, 2, 9, 10), (3, 4, 8)], 
        [(2, 9, 10)], 
        [(1, 9, 10)], 
        [(3, 4, 5), (2, 10)], 
        [(2, 3, 4)], 
        [(3, 4)], 
        [(9, 10), (2, 3)],
        ]

    assert result.reduced_hierarchy == [
        [(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)], 
        [(1, 2, 3, 6, 7, 9, 10)], 
        [(1, 2, 9, 10), (3, 4, 8)], 
        [(2, 9, 10), (3, 4, 5)], 
        [(1, 9, 10), (2, 3, 4)], 
        [(2, 10), (3, 4)], 
        [(9, 10), (2, 3)],
        ]

    um = result.threshold_unique_modules
    assert list(um.keys()) == [10.0, 9.0, 8.0, 7.0, 6.0, 4.0, 3.0, 2.0]
    assert um[10.0] == [(9, 10), (2, 3)]
    assert um[9.0] == [(3, 4)]
    assert um[8.0] == [(2, 3, 4)]
    assert um[7.0] == [(3, 4, 5), (2, 10)]
    assert um[6.0] == [(2, 9, 10), (1, 9, 10)]
    assert um[2.0] == [(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)]

    uh = result.threshold_unique_hierarchies
    assert list(uh.keys()) == [10.0, 9.0, 8.0, 7.0, 6.0, 4.0, 3.0, 2.0]
    assert uh[10.0] == [[(9, 10), (2, 3)]]
    assert uh[9.0] == [[(3, 4)]]
    assert uh[8.0] == [[(2, 3, 4)]]
    assert uh[7.0] == [[(3, 4, 5), (2, 10)]]
    assert uh[6.0] == [[(2, 9, 10)], [(1, 9, 10)]]
    assert uh[2.0] == [[(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)]]

    ch = result.threshold_cumulative_hierarchies
    assert list(ch.keys()) == [10.0, 9.0, 8.0, 7.0, 6.0, 4.0, 3.0, 2.0]
    assert ch[10.0] == [
        [(9, 10), (2, 3)],
        ]
    assert ch[9.0] == [
        [(3, 4)], 
        [(9, 10), (2, 3)],
        ]
    assert ch[8.0] == [
        [(2, 3, 4)], 
        [(3, 4)], 
        [(9, 10), (2, 3)],
        ]
    assert ch[7.0] == [
        [(3, 4, 5), (2, 10)], 
        [(2, 3, 4)], 
        [(3, 4)], 
        [(9, 10), (2, 3)],
        ]
    assert ch[6.0] == [
        [(2, 9, 10)], 
        [(1, 9, 10)], 
        [(3, 4, 5), (2, 10)], 
        [(2, 3, 4)], 
        [(3, 4)], 
        [(9, 10), (2, 3)],
        ]
    assert ch[2.0] == [
        [(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)], 
        [(1, 2, 3, 6, 7, 9, 10)], 
        [(1, 2, 9, 10), (3, 4, 8)], 
        [(2, 9, 10)], 
        [(1, 9, 10)], 
        [(3, 4, 5), (2, 10)], 
        [(2, 3, 4)], 
        [(3, 4)], 
        [(9, 10), (2, 3)],
        ]

    rh = result.threshold_reduced_hierarchies
    assert list(rh.keys()) == [10.0, 9.0, 8.0, 7.0, 6.0, 4.0, 3.0, 2.0]
    assert rh[10.0] == [
        [(9, 10), (2, 3)],
        ]
    assert rh[9.0] == [
        [(3, 4), (9, 10)], 
        [(2, 3)],
        ]
    assert rh[8.0] == [
        [(2, 3, 4), (9, 10)], 
        [(3, 4)], [(2, 3)],
        ]
    assert rh[7.0] == [
        [(3, 4, 5), (2, 10)], 
        [(2, 3, 4), (9, 10)], 
        [(3, 4)], [(2, 3)],
        ]
    assert rh[6.0] == [
        [(2, 9, 10), (3, 4, 5)], 
        [(1, 9, 10), (2, 3, 4)], 
        [(2, 10), (3, 4)], 
        [(9, 10), (2, 3)],
        ]
    assert rh[2.0] == [
        [(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)], 
        [(1, 2, 3, 6, 7, 9, 10)], 
        [(1, 2, 9, 10), (3, 4, 8)], 
        [(2, 9, 10), (3, 4, 5)], 
        [(1, 9, 10), (2, 3, 4)], 
        [(2, 10), (3, 4)], 
        [(9, 10), (2, 3)],
        ]

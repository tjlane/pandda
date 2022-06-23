from pytest import approx, raises

import numpy

uijs = [
    [
        [ 0.26348568,  0.08803938,  0.92739069,  0.81716141,  0.40285524,  0.65731412],
        [ 0.28450075,  0.45968414,  0.88807511,  0.85285115,  0.03640805,  0.15434743],
        [ 0.42020582,  0.24677701,  0.96400143,  0.42411038,  0.44711358,  0.22477824],
        [ 0.75970411,  0.9087415 ,  0.40797275,  0.03218416,  0.57227292,  0.34488995],
        ],
    [
        [ 0.85047116,  0.66417643,  0.88500985,  0.54441866,  0.58203759,  0.43958623],
        [ 0.99508565,  0.15226961,  0.78888134,  0.51321502,  0.26370695,  0.88747190],
        [ 0.6920631 ,  0.76229571,  0.36675838,  0.62061572,  0.17007319,  0.79687877],
        [ 0.43784481,  0.93796117,  0.3786025 ,  0.84214932,  0.01153883,  0.76836497],
        ],
    ]

resolutions = [1.4, 3.2]

n_datasets = 2    
n_atoms = 4

def test_atom_weight_calculator():
    from pandemic.adp.weights import AtomWeightCalculator

    uij_array = numpy.array(uijs)
    assert uij_array.shape == (n_datasets, n_atoms, 6)

    uij_mods = uij_array[:,:,0:3].mean(axis=2)
    assert uij_mods.shape == (n_datasets, n_atoms)

    for w_str, ref_wgts in [
            (
                'one', 
                numpy.ones((n_datasets, n_atoms)),
                ),
            (
                'inverse_mod_U',
                1. / uij_mods,
                ),
            (
                'inverse_mod_U_squared',
                1. / (uij_mods*uij_mods),
                ),
            (
                'inverse_mod_U_cubed',
                1. / (uij_mods*uij_mods*uij_mods),
                ),
            ]:

        wc = AtomWeightCalculator(
            weighting = w_str,
            renormalise_by_dataset = False,
            )
        wgts = wc(uij_array)
        assert wgts == approx(ref_wgts / ref_wgts.mean())

        #############################

        wc = AtomWeightCalculator(
            weighting = w_str,
            renormalise_by_dataset = True,
            )
        wgts = wc(uij_array)
        ref_wgts_norm = numpy.array([w / w.mean() for w in ref_wgts])
        assert wgts == approx(ref_wgts_norm / ref_wgts_norm.mean())

def test_dataset_weight_calculator():
    from pandemic.adp.weights import DatasetWeightCalculator

    res_array = numpy.array(resolutions)
    assert res_array.shape == (n_datasets,)

    for w_str, ref_wgts in [
            (
                'one', 
                numpy.ones(n_datasets),
                ),
            (
                'inverse_resolution',
                1. / res_array,
                ),
            (
                'inverse_resolution_squared',
                1. / (res_array*res_array),
                ),
            (
                'inverse_resolution_cubed',
                1. / (res_array*res_array*res_array),
                ),
            ]:

        wc = DatasetWeightCalculator(
            weighting = w_str,
            )
        wgts = wc(
            resolutions = resolutions,
            dataset_labels = None,
            )
        assert wgts == approx(ref_wgts / ref_wgts.mean())

        #############################

        if w_str == 'one':
            wgts = wc(
                resolutions = None,
                dataset_labels = list(range(2)),
                )
            assert wgts == approx(ref_wgts / ref_wgts.mean())
        else:
            with raises(Exception):
                wgts = wc(
                    resolutions = None,
                    dataset_labels = list(range(2)),
                    )
            with raises(Exception):
                wgts = wc(
                    resolutions = None,
                    dataset_labels = None,
                    )

def test_uij_array_weight_task():
    from pandemic.adp.weights import UijArrayWeightsTask
    
    uij_array = numpy.array(uijs)
    assert uij_array.shape == (n_datasets, n_atoms, 6)

    uij_mods = uij_array[:,:,0:3].mean(axis=2)
    assert uij_mods.shape == (n_datasets, n_atoms)

    res_array = numpy.array(resolutions)
    assert res_array.shape == (n_datasets,)

    for a_w_str, a_ref_wgts in [
            (
                'one', 
                numpy.ones((n_datasets, n_atoms)),
                ),
            (
                'inverse_mod_U',
                1. / uij_mods,
                ),
            (
                'inverse_mod_U_squared',
                1. / (uij_mods*uij_mods),
                ),
            (
                'inverse_mod_U_cubed',
                1. / (uij_mods*uij_mods*uij_mods),
                ),
            ]:

        for d_w_str, d_ref_wgts in [
                (
                    'one', 
                    numpy.ones(n_datasets),
                    ),
                (
                    'inverse_resolution',
                    1. / res_array,
                    ),
                (
                    'inverse_resolution_squared',
                    1. / (res_array*res_array),
                    ),
                (
                    'inverse_resolution_cubed',
                    1. / (res_array*res_array*res_array),
                    ),
                ]:

            wc = UijArrayWeightsTask(
                dataset_weighting = d_w_str,
                atom_weighting = a_w_str,
                renormalise_by_dataset = False,
                )
            wgts = wc.run(
                uij_values = uij_array,
                resolutions = list(resolutions),
                dataset_labels = None,
                )

            ref_wgts = numpy.array([d_ref_wgts[i] * a_ref_wgts[i] for i in range(n_datasets)])
            assert wgts.atom_weight_array    == approx(a_ref_wgts / a_ref_wgts.mean())
            assert wgts.dataset_weight_array == approx(d_ref_wgts / d_ref_wgts.mean())
            assert wgts.total_weight_array   == approx(ref_wgts / ref_wgts.mean())

            #############################
        
            wc = UijArrayWeightsTask(
                dataset_weighting = d_w_str,
                atom_weighting = a_w_str,
                renormalise_by_dataset = True,
                )
            wgts = wc.run(
                uij_values = uij_array,
                resolutions = list(resolutions),
                dataset_labels = None,
                )

            a_ref_wgts_norm = numpy.array([a / a.mean() for a in a_ref_wgts])
            ref_wgts = numpy.array([d_ref_wgts[i] * a_ref_wgts_norm[i] for i in range(n_datasets)])
            assert wgts.atom_weight_array    == approx(a_ref_wgts_norm / a_ref_wgts_norm.mean())
            assert wgts.dataset_weight_array == approx(d_ref_wgts / d_ref_wgts.mean())
            assert wgts.total_weight_array   == approx(ref_wgts / ref_wgts.mean())
        
            #############################

            wc = UijArrayWeightsTask(
                dataset_weighting = d_w_str,
                atom_weighting = a_w_str,
                renormalise_by_dataset = True,
                )

            if d_w_str == 'one':
                wgts = wc.run(
                    uij_values = uij_array,
                    resolutions = None,
                    dataset_labels = None,
                    )
                assert wgts.total_weight_array == approx(ref_wgts / ref_wgts.mean())
            else:
                with raises(Exception):
                    wgts = wc.run(
                        uij_values = uij_array,
                        resolutions = None,
                        dataset_labels = list(range(n_datasets)),
                        )
                with raises(Exception):
                    wgts = wc.run(
                        uij_values = uij_array,
                        resolutions = None,
                        dataset_labels = None,
                        )
























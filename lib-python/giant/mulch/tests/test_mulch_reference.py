import pytest


def test_DefaultReferenceDataset(five_baz2b_test_datasets_mcd):

    from giant.mulch.reference import DefaultReferenceDataset

    input_dataset = five_baz2b_test_datasets_mcd.datasets["BAZ2BA-x434.dimple.pdb"]

    reference_dataset = DefaultReferenceDataset(
        model = input_dataset.model,
        data = input_dataset.data,
        )

    assert reference_dataset.model is input_dataset.model
    assert reference_dataset.data is input_dataset.data

def test_DefaultReferenceSelector(five_baz2b_test_datasets_mcd):

    from giant.mulch.reference import DefaultReferenceSelector

    ##

    selector = DefaultReferenceSelector(
        selection_method = "resolution",
        )

    selected_dataset = selector(
        mcd = five_baz2b_test_datasets_mcd,
        )

    get_rhigh = lambda d: d.data.mtz_object().max_min_resolution()[1]

    assert get_rhigh(selected_dataset) == min([
        get_rhigh(d) 
        for d in five_baz2b_test_datasets_mcd.datasets.values()
        ])
    assert get_rhigh(selected_dataset) == pytest.approx(1.6545, abs=1e-4)

    ##

    selector = DefaultReferenceSelector(
        selection_method = "rfree",
        )

    selected_dataset = selector(
        mcd = five_baz2b_test_datasets_mcd,
        )

    get_rfree = lambda d: d.model.input.get_r_rfree_sigma().r_free

    assert get_rfree(selected_dataset) == min([
        get_rfree(d) 
        for d in five_baz2b_test_datasets_mcd.datasets.values()
        ])
    assert get_rfree(selected_dataset) == pytest.approx(0.2044, abs=1e-4)


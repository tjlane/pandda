import os, copy
import pytest

#####

def test_ManualDatasetFilter(five_baz2b_test_datasets_mcd):

    from giant.mulch.filters import ManualDatasetFilter
    
    # Initialise datasets to be rejected (with reasons)
    filter = ManualDatasetFilter(
        rejections = {
            "BAZ2BA-x432.dimple.pdb" : "rejected",
            },
        )

    mcd = filter(five_baz2b_test_datasets_mcd)

    assert sorted(mcd.datasets.keys()) == [
        "BAZ2BA-x430.dimple.pdb",
        "BAZ2BA-x431.dimple.pdb",
        "BAZ2BA-x433.dimple.pdb",
        "BAZ2BA-x434.dimple.pdb",
        ]

    ##

    # Should also work from list
    filter = ManualDatasetFilter(
        rejections = [
            "BAZ2BA-x432.dimple.pdb",
            ],
        )

    mcd = filter(five_baz2b_test_datasets_mcd)

    assert sorted(mcd.datasets.keys()) == [
        "BAZ2BA-x430.dimple.pdb",
        "BAZ2BA-x431.dimple.pdb",
        "BAZ2BA-x433.dimple.pdb",
        "BAZ2BA-x434.dimple.pdb",
        ]

def test_InclusiveDatasetTagFilter(five_baz2b_test_datasets_mcd_labelled):

    from giant.mulch.filters import InclusiveDatasetTagFilter
    
    filter = InclusiveDatasetTagFilter(
        dataset_tags = [
            "BAZ2BA-x432.dimple.pdb",
            ],
        )

    mcd = filter(five_baz2b_test_datasets_mcd_labelled)

    assert sorted(mcd.datasets.keys()) == [
        "BAZ2BA-x432.dimple.pdb",
        ]

def test_ExclusiveDatasetTagFilter(five_baz2b_test_datasets_mcd_labelled):

    from giant.mulch.filters import ExclusiveDatasetTagFilter
    
    filter = ExclusiveDatasetTagFilter(
        dataset_tags = [
            "BAZ2BA-x432.dimple.pdb",
            ]
        )

    mcd = filter(five_baz2b_test_datasets_mcd_labelled)

    assert sorted(mcd.datasets.keys()) == [
        "BAZ2BA-x430.dimple.pdb",
        "BAZ2BA-x431.dimple.pdb",
        "BAZ2BA-x433.dimple.pdb",
        "BAZ2BA-x434.dimple.pdb",
        ]

def test_HighResolutionFilter(five_baz2b_test_datasets_mcd):

    from giant.mulch.filters import HighResolutionFilter
    
    filter = HighResolutionFilter(high_resolution_cutoff=1.7)

    mcd = filter(five_baz2b_test_datasets_mcd)

    assert sorted(mcd.datasets.keys()) == [
        "BAZ2BA-x433.dimple.pdb",
        "BAZ2BA-x434.dimple.pdb",
        ]

def test_RValueFilter(five_baz2b_test_datasets_mcd):

    from giant.mulch.filters import RValueFilter
    
    ##

    filter = RValueFilter(max_rfree=0.210)

    mcd = filter(five_baz2b_test_datasets_mcd)

    assert sorted(mcd.datasets.keys()) == [
        "BAZ2BA-x433.dimple.pdb",
        "BAZ2BA-x434.dimple.pdb",
        ]

    ##

    filter = RValueFilter(max_rwork=0.182)

    mcd = filter(five_baz2b_test_datasets_mcd)

    assert sorted(mcd.datasets.keys()) == [
        "BAZ2BA-x430.dimple.pdb",
        "BAZ2BA-x434.dimple.pdb",
        ]

    ##

    filter = RValueFilter(max_rfree=0.210, max_rwork=0.182)

    mcd = filter(five_baz2b_test_datasets_mcd)

    assert sorted(mcd.datasets.keys()) == [
        "BAZ2BA-x434.dimple.pdb",
        ]

def test_SpaceGroupFilter(five_baz2b_test_datasets_mcd):

    from giant.mulch.filters import SpaceGroupFilter    

    filter = SpaceGroupFilter(
        reference_dataset = list(five_baz2b_test_datasets_mcd.datasets.values())[0],
        )

    mcd = filter(five_baz2b_test_datasets_mcd)

    assert sorted(mcd.datasets.keys()) == [
        "BAZ2BA-x430.dimple.pdb",
        "BAZ2BA-x431.dimple.pdb",
        "BAZ2BA-x432.dimple.pdb",
        "BAZ2BA-x433.dimple.pdb",
        "BAZ2BA-x434.dimple.pdb",
        ]

    ## Now edit a couple of datasets to have different SGs

    from cctbx.sgtbx import space_group
    five_baz2b_test_datasets_mcd_copy = copy.deepcopy(five_baz2b_test_datasets_mcd)
    assert five_baz2b_test_datasets_mcd_copy.datasets["BAZ2BA-x432.dimple.pdb"].model.crystal.space_group is not None
    assert five_baz2b_test_datasets_mcd_copy.datasets["BAZ2BA-x434.dimple.pdb"].model.crystal.space_group is not None
    five_baz2b_test_datasets_mcd_copy.datasets["BAZ2BA-x432.dimple.pdb"].model.crystal.space_group = space_group("P1")
    five_baz2b_test_datasets_mcd_copy.datasets["BAZ2BA-x434.dimple.pdb"].model.crystal.space_group = space_group("P1")

    mcd = filter(five_baz2b_test_datasets_mcd_copy)

    assert sorted(mcd.datasets.keys()) == [
        "BAZ2BA-x430.dimple.pdb",
        "BAZ2BA-x431.dimple.pdb",
        "BAZ2BA-x433.dimple.pdb",
        ]

def test_IdenticalHierarchyFilter(five_baz2b_test_datasets_mcd):

    from giant.mulch.filters import IdenticalHierarchyFilter

    filter = IdenticalHierarchyFilter(
        reference_dataset = list(five_baz2b_test_datasets_mcd.datasets.values())[0],
        atom_selection_string = "pepnames",
        )

    mcd = filter(five_baz2b_test_datasets_mcd)

    assert sorted(mcd.datasets.keys()) == [
        "BAZ2BA-x430.dimple.pdb",
        "BAZ2BA-x431.dimple.pdb",
        "BAZ2BA-x432.dimple.pdb",
        "BAZ2BA-x433.dimple.pdb",
        "BAZ2BA-x434.dimple.pdb",
        ]

##

def test_DefaultDatasetFilter(five_baz2b_test_datasets_mcd):

    from giant.mulch.filters import DefaultDatasetFilter

    filter = DefaultDatasetFilter(
        same_space_group_only = True,
        similar_models_only = False,
        max_rfree = 0.210,
        max_rwork = 0.182,
        reference_dataset = list(five_baz2b_test_datasets_mcd.datasets.values())[0],
        )

    mcd = filter(five_baz2b_test_datasets_mcd)
    
    assert sorted(mcd.datasets.keys()) == [
        "BAZ2BA-x434.dimple.pdb",
        ]

import os
import pytest

#####

@pytest.fixture
def five_baz2b_test_datasets_mcd(five_baz2b_test_datasets):

    from giant.mulch.collection import MultiCrystalDataset

    mcd = MultiCrystalDataset({
        os.path.basename(d.model.filename) : d 
        for d in five_baz2b_test_datasets
        })

    return mcd

def test_DummyDatasetFinder(five_baz2b_test_datasets_mcd):

    from giant.mulch.finders import DummyDatasetFinder

    finder = DummyDatasetFinder()

    dkey = finder(five_baz2b_test_datasets_mcd)

    assert dkey == list(five_baz2b_test_datasets_mcd.datasets.keys())[0] # dummy picks the first key

def test_HighestResolutionFinder(five_baz2b_test_datasets_mcd):

    from giant.mulch.finders import HighestResolutionFinder

    finder = HighestResolutionFinder()

    dkey = finder(five_baz2b_test_datasets_mcd) 

    assert dkey == "BAZ2BA-x434.dimple.pdb" # has highest resolution

def test_LowestRfreeFinder(five_baz2b_test_datasets_mcd):

    from giant.mulch.finders import LowestRfreeFinder

    finder = LowestRfreeFinder()

    dkey = finder(five_baz2b_test_datasets_mcd)

    assert dkey == "BAZ2BA-x434.dimple.pdb" # has lowest rfree
import pytest

#####

def test_MultiCrystalDataset(five_baz2b_test_datasets):

    from giant.mulch.collection import MultiCrystalDataset

    datasets = {
        i:d 
        for i, d in enumerate(five_baz2b_test_datasets)
    }

    mcd = MultiCrystalDataset(datasets=datasets)

    ##

    datasets_list = list(mcd)
    expected_list = list(datasets.values())

    assert [
        (d1 is d2) 
        for d1, d2 in zip(datasets_list, expected_list)
        ].count(True) == 5

    ##

    mcd.dataset_keys.sort()

    datasets_list = list(mcd)
    expected_list = [datasets[i] for i in range(5)]

    assert [
        (d1 is d2) 
        for d1, d2 in zip(datasets_list, expected_list)
        ].count(True) == 5

    ##

    assert mcd.resolution_high == pytest.approx(1.65)
    assert mcd.resolution_low == pytest.approx(1.84)

    ##

    assert mcd.n_datasets() == 5

    ##

    new_datasets = {
        i:d 
        for i, d in enumerate(five_baz2b_test_datasets[:3])
        }

    mcd_2 = mcd.new_from_datasets(datasets=new_datasets)

    assert mcd_2.dataset_keys == list(new_datasets.keys())
    assert list(mcd_2.datasets.keys()) == list(new_datasets.keys())

    for i in range(3):
        assert datasets[i] is mcd_2.datasets[i]

    ##

    from giant.mulch.partitioners import DatasetKeyPartitioner

    partitioner = DatasetKeyPartitioner(
        partition_1 = [0,1,2],
        partition_2 = [1,4],
        )

    partitioned = mcd.partition(partitioner)

    assert sorted(partitioned.keys()) == ["partition_1","partition_2"]

    assert list(partitioned["partition_1"].datasets.keys()) == [0,1,2]
    assert list(partitioned["partition_2"].datasets.keys()) == [1,4]

    assert partitioned["partition_1"].datasets[1] is datasets[1]
    assert partitioned["partition_2"].datasets[4] is datasets[4]

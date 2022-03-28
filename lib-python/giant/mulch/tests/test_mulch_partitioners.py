import os
import pytest

#####

def test_DatasetKeyPartitioner(five_baz2b_test_datasets_mcd):

    from giant.mulch.partitioners import DatasetKeyPartitioner

    partitioner = DatasetKeyPartitioner(
        partition_1 = [
            "BAZ2BA-x430.dimple.pdb",
            "BAZ2BA-x434.dimple.pdb",
            ],
        partition_2 = [
            "BAZ2BA-x430.dimple.pdb",
            "BAZ2BA-x431.dimple.pdb",
            "BAZ2BA-x432.dimple.pdb",
            "nonexistent",
            ],
        partition_3 = [
            "BAZ2BA-x430.dimple.pdb",
            "BAZ2BA-x431.dimple.pdb",
            "BAZ2BA-x432.dimple.pdb",
            "BAZ2BA-x433.dimple.pdb",
            "BAZ2BA-x434.dimple.pdb",
            "nonexistent",
            ],
        partition_4 = [
            "nonexistent1",
            "nonexistent2",
            ],
        )

    partitioned = partitioner(five_baz2b_test_datasets_mcd.datasets)

    assert sorted(partitioned.keys()) == [
        "partition_1",
        "partition_2",
        "partition_3",
        "partition_4",
        ]

    assert sorted(partitioned["partition_1"].keys()) == [
            "BAZ2BA-x430.dimple.pdb",
            "BAZ2BA-x434.dimple.pdb",
            ]
    assert sorted(partitioned["partition_2"].keys()) == [
            "BAZ2BA-x430.dimple.pdb",
            "BAZ2BA-x431.dimple.pdb",
            "BAZ2BA-x432.dimple.pdb",
            ]
    assert sorted(partitioned["partition_3"].keys()) == [
            "BAZ2BA-x430.dimple.pdb",
            "BAZ2BA-x431.dimple.pdb",
            "BAZ2BA-x432.dimple.pdb",
            "BAZ2BA-x433.dimple.pdb",
            "BAZ2BA-x434.dimple.pdb",
            ]
    assert sorted(partitioned["partition_4"].keys()) == []

    assert partitioned["partition_1"]["BAZ2BA-x430.dimple.pdb"] is five_baz2b_test_datasets_mcd.datasets["BAZ2BA-x430.dimple.pdb"]
    assert partitioned["partition_2"]["BAZ2BA-x432.dimple.pdb"] is five_baz2b_test_datasets_mcd.datasets["BAZ2BA-x432.dimple.pdb"]
    assert partitioned["partition_3"]["BAZ2BA-x434.dimple.pdb"] is five_baz2b_test_datasets_mcd.datasets["BAZ2BA-x434.dimple.pdb"]

def test_TestTrainPartitioner(five_baz2b_test_datasets_mcd):

    from giant.mulch.partitioners import TestTrainPartitioner

    ##

    if True: 
        partitioner = TestTrainPartitioner(
            test = ["BAZ2BA-x430.dimple.pdb","BAZ2BA-x431.dimple.pdb"],
            train = ["BAZ2BA-x431.dimple.pdb","BAZ2BA-x432.dimple.pdb","BAZ2BA-x433.dimple.pdb"],
            not_test = None,
            not_train = None,
            test_selector = None,
            train_selector = None,
            )

        partitioned = partitioner(five_baz2b_test_datasets_mcd.datasets)

        assert sorted(partitioned.keys()) == ["test","train"]
        assert sorted(partitioned["test"].keys()) == [
            "BAZ2BA-x430.dimple.pdb",
            "BAZ2BA-x431.dimple.pdb",
            ]
        assert sorted(partitioned["train"].keys()) == [
            "BAZ2BA-x431.dimple.pdb",
            "BAZ2BA-x432.dimple.pdb",
            "BAZ2BA-x433.dimple.pdb",
            ]

    ##

    if True:
        partitioner = TestTrainPartitioner(
            test = ["BAZ2BA-x430.dimple.pdb","BAZ2BA-x431.dimple.pdb"],
            train = None,
            not_test = None,
            not_train = None,
            test_selector = None,
            train_selector = None,
            )

        partitioned = partitioner(five_baz2b_test_datasets_mcd.datasets)

        assert sorted(partitioned.keys()) == ["test","train"]
        assert sorted(partitioned["test"].keys()) == [
            "BAZ2BA-x430.dimple.pdb",
            "BAZ2BA-x431.dimple.pdb",
            ]
        assert sorted(partitioned["train"].keys()) == [
            "BAZ2BA-x430.dimple.pdb",
            "BAZ2BA-x431.dimple.pdb",
            "BAZ2BA-x432.dimple.pdb",
            "BAZ2BA-x433.dimple.pdb",
            "BAZ2BA-x434.dimple.pdb",
            ]

    ##

    if True: 
        partitioner = TestTrainPartitioner(
            test = None,
            train = ["BAZ2BA-x431.dimple.pdb","BAZ2BA-x432.dimple.pdb","BAZ2BA-x433.dimple.pdb"],
            not_test = None,
            not_train = None,
            test_selector = None,
            train_selector = None,
            )

        partitioned = partitioner(five_baz2b_test_datasets_mcd.datasets)

        assert sorted(partitioned.keys()) == ["test","train"]
        assert sorted(partitioned["test"].keys()) == [
            "BAZ2BA-x430.dimple.pdb",
            "BAZ2BA-x431.dimple.pdb",
            "BAZ2BA-x432.dimple.pdb",
            "BAZ2BA-x433.dimple.pdb",
            "BAZ2BA-x434.dimple.pdb",
            ]
        assert sorted(partitioned["train"].keys()) == [
            "BAZ2BA-x431.dimple.pdb",
            "BAZ2BA-x432.dimple.pdb",
            "BAZ2BA-x433.dimple.pdb",
            ]

    ##
    
    if True: 
        partitioner = TestTrainPartitioner(
            test = None,
            train = None,
            not_test = ["BAZ2BA-x430.dimple.pdb","BAZ2BA-x431.dimple.pdb"],
            not_train = ["BAZ2BA-x431.dimple.pdb","BAZ2BA-x432.dimple.pdb","BAZ2BA-x433.dimple.pdb"],
            test_selector = None,
            train_selector = None,
            )

        partitioned = partitioner(five_baz2b_test_datasets_mcd.datasets)

        assert sorted(partitioned.keys()) == ["test","train"]
        assert sorted(partitioned["test"].keys()) == [
            "BAZ2BA-x432.dimple.pdb",
            "BAZ2BA-x433.dimple.pdb",
            "BAZ2BA-x434.dimple.pdb",
            ]
        assert sorted(partitioned["train"].keys()) == [
            "BAZ2BA-x430.dimple.pdb",
            "BAZ2BA-x434.dimple.pdb",
            ]

    ##

    if True: 
        partitioner = TestTrainPartitioner(
            test = None,
            train = None,
            not_test = None,
            not_train = ["BAZ2BA-x431.dimple.pdb","BAZ2BA-x432.dimple.pdb","BAZ2BA-x433.dimple.pdb"],
            test_selector = None,
            train_selector = None,
            )

        partitioned = partitioner(five_baz2b_test_datasets_mcd.datasets)

        assert sorted(partitioned.keys()) == ["test","train"]
        assert sorted(partitioned["test"].keys()) == [
            "BAZ2BA-x430.dimple.pdb",
            "BAZ2BA-x431.dimple.pdb",
            "BAZ2BA-x432.dimple.pdb",
            "BAZ2BA-x433.dimple.pdb",
            "BAZ2BA-x434.dimple.pdb",
            ]
        assert sorted(partitioned["train"].keys()) == [
            "BAZ2BA-x430.dimple.pdb",
            "BAZ2BA-x434.dimple.pdb",
            ]

    ##

    if True: 
        partitioner = TestTrainPartitioner(
            test = None,
            train = None,
            not_test = ["BAZ2BA-x430.dimple.pdb","BAZ2BA-x431.dimple.pdb"],
            not_train = None,
            test_selector = None,
            train_selector = None,
            )

        partitioned = partitioner(five_baz2b_test_datasets_mcd.datasets)

        assert sorted(partitioned.keys()) == ["test","train"]
        assert sorted(partitioned["test"].keys()) == [
            "BAZ2BA-x432.dimple.pdb",
            "BAZ2BA-x433.dimple.pdb",
            "BAZ2BA-x434.dimple.pdb",
            ]
        assert sorted(partitioned["train"].keys()) == [
            "BAZ2BA-x430.dimple.pdb",
            "BAZ2BA-x431.dimple.pdb",
            "BAZ2BA-x432.dimple.pdb",
            "BAZ2BA-x433.dimple.pdb",
            "BAZ2BA-x434.dimple.pdb",
            ]

    ##

    # illogical selections should raise errors

    with pytest.raises(ValueError) as e:
        partitioner = TestTrainPartitioner(
            test = ["BAZ2BA-x430.dimple.pdb","BAZ2BA-x431.dimple.pdb"],
            train = None,
            not_test = ["BAZ2BA-x430.dimple.pdb"], 
            not_train = None,
            test_selector = None,
            train_selector = None,
            )

    with pytest.raises(ValueError) as e:
        partitioner = TestTrainPartitioner(
            test = None,
            train = ["BAZ2BA-x430.dimple.pdb","BAZ2BA-x431.dimple.pdb"],
            not_test = None, 
            not_train = ["BAZ2BA-x430.dimple.pdb"],
            test_selector = None,
            train_selector = None,
            )

    with pytest.raises(ValueError) as e:
        partitioner = TestTrainPartitioner(
            test = ["BAZ2BA-x433.dimple.pdb","BAZ2BA-x434.dimple.pdb"],
            train = ["BAZ2BA-x430.dimple.pdb","BAZ2BA-x431.dimple.pdb"],
            not_test = ["BAZ2BA-x433.dimple.pdb"], 
            not_train = ["BAZ2BA-x430.dimple.pdb"],
            test_selector = None,
            train_selector = None,
            )

    ##

    from giant.mulch.selectors import SortedDatasetSelector
    from giant.mulch.sorters import ArbitraryDatasetSorter

    if True:

        partitioner = TestTrainPartitioner(
            test = [
                "BAZ2BA-x431.dimple.pdb",
                "BAZ2BA-x432.dimple.pdb",
                "BAZ2BA-x433.dimple.pdb",
                "BAZ2BA-x434.dimple.pdb",
                ],
            train = [
                "BAZ2BA-x430.dimple.pdb",
                "BAZ2BA-x434.dimple.pdb",
                ],
            not_test = None,
            not_train = None,
            test_selector = SortedDatasetSelector(
                max_datasets = 3,
                sort_datasets = ArbitraryDatasetSorter(
                    sort_function = lambda d: d.model.filename,
                    ),
                ),
            train_selector = None,
            )

        partitioned = partitioner(five_baz2b_test_datasets_mcd.datasets)

        assert sorted(partitioned.keys()) == ["test","train"]
        assert sorted(partitioned["test"].keys()) == [
            "BAZ2BA-x431.dimple.pdb",
            "BAZ2BA-x432.dimple.pdb",
            "BAZ2BA-x433.dimple.pdb",
            ]
        assert sorted(partitioned["train"].keys()) == [
            "BAZ2BA-x430.dimple.pdb",
            "BAZ2BA-x434.dimple.pdb",
            ]

    if True: 

        partitioner = TestTrainPartitioner(
            test = [
                "BAZ2BA-x430.dimple.pdb",
                "BAZ2BA-x434.dimple.pdb",
                ],
            train = [
                "BAZ2BA-x431.dimple.pdb",
                "BAZ2BA-x432.dimple.pdb",
                "BAZ2BA-x433.dimple.pdb",
                "BAZ2BA-x434.dimple.pdb",
                ],
            not_test = None,
            not_train = None,
            test_selector = None,
            train_selector = SortedDatasetSelector(
                max_datasets = 3,
                sort_datasets = ArbitraryDatasetSorter(
                    sort_function = lambda d: d.model.filename,
                    ),
                ),
            )

        partitioned = partitioner(five_baz2b_test_datasets_mcd.datasets)

        assert sorted(partitioned.keys()) == ["test","train"]
        assert sorted(partitioned["test"].keys()) == [
            "BAZ2BA-x430.dimple.pdb",
            "BAZ2BA-x434.dimple.pdb",
            ]
        assert sorted(partitioned["train"].keys()) == [
            "BAZ2BA-x431.dimple.pdb",
            "BAZ2BA-x432.dimple.pdb",
            "BAZ2BA-x433.dimple.pdb",
            ]

    ##

    if True: 

        partitioner = TestTrainPartitioner(
            test = [
                "BAZ2BA-x430.dimple.pdb",
                "BAZ2BA-x434.dimple.pdb",
                ],
            train = [
                "BAZ2BA-x431.dimple.pdb",
                "BAZ2BA-x432.dimple.pdb",
                "BAZ2BA-x433.dimple.pdb",
                "BAZ2BA-x434.dimple.pdb",
                ],
            not_test = None,
            not_train = None,
            test_selector = SortedDatasetSelector(
                max_datasets = 1,
                sort_datasets = ArbitraryDatasetSorter(
                    sort_function = lambda d: d.model.filename,
                    ),
                ),
            train_selector = SortedDatasetSelector(
                max_datasets = 3,
                sort_datasets = ArbitraryDatasetSorter(
                    sort_function = lambda d: d.model.filename,
                    ),
                ),
            )

        partitioned = partitioner(five_baz2b_test_datasets_mcd.datasets)

        assert sorted(partitioned.keys()) == ["test","train"]
        assert sorted(partitioned["test"].keys()) == [
            "BAZ2BA-x430.dimple.pdb",
            ]
        assert sorted(partitioned["train"].keys()) == [
            "BAZ2BA-x431.dimple.pdb",
            "BAZ2BA-x432.dimple.pdb",
            "BAZ2BA-x433.dimple.pdb",
            ]



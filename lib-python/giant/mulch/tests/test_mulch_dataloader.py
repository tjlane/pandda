import pytest

def test_MultiDatasetDataloader(tmp_path):

    from giant.resources.test_data import TEST_DATA

    TEST_DATA.extract_to_temporary_directory(
        dataset = "BAZ2B",
        tmp_dir = str(tmp_path),
        n = 5,
        )

    from giant.mulch.dataloader import MultiDatasetDataloader
    
    ##

    if True: 

        loader = MultiDatasetDataloader(
            data_dirs = str(tmp_path / "data" / "*"),
            pdb_style = "*.dimple.pdb",
            mtz_style = "*.dimple.mtz",
            pdb_regex = None,
            mtz_regex = None,
            dir_regex = None,
            only_datasets = None,
            ignore_datasets = None,
            dataset_prefix = None,
            )

        mcd = loader()
        datasets = mcd.datasets

        assert len(datasets) == 5

        assert sorted(datasets.keys()) == [
            "BAZ2BA-x430","BAZ2BA-x431","BAZ2BA-x432","BAZ2BA-x433","BAZ2BA-x434",
            ]

        # Check string function does not error
        _ = str(loader)

        # investigate a couple of the datasets
        d = datasets["BAZ2BA-x433"]
        assert d.model.filename == str(tmp_path / "data" / "BAZ2BA-x433" / "BAZ2BA-x433.dimple.pdb")
        assert d.data.filename == str(tmp_path / "data" / "BAZ2BA-x433" / "BAZ2BA-x433.dimple.mtz")

        ##

        loader = MultiDatasetDataloader(
            data_dirs = str(tmp_path / "data" / "*"),
            pdb_style = "*.dimple.pdb",
            mtz_style = None,
            pdb_regex = None,
            mtz_regex = None,
            dir_regex = None,
            only_datasets = None,
            ignore_datasets = None,
            dataset_prefix = None,
            )

        mcd = loader()
        datasets = mcd.datasets

        assert len(datasets) == 5

        assert sorted(datasets.keys()) == [
            "BAZ2BA-x430","BAZ2BA-x431","BAZ2BA-x432","BAZ2BA-x433","BAZ2BA-x434",
            ]

        ##

        loader = MultiDatasetDataloader(
            data_dirs = str(tmp_path / "data" / "*"),
            pdb_style = None,
            mtz_style = "*.dimple.mtz",
            pdb_regex = None,
            mtz_regex = None,
            dir_regex = None,
            only_datasets = None,
            ignore_datasets = None,
            dataset_prefix = None,
            )

        mcd = loader()
        datasets = mcd.datasets

        assert len(datasets) == 5

        assert sorted(datasets.keys()) == [
            "BAZ2BA-x430","BAZ2BA-x431","BAZ2BA-x432","BAZ2BA-x433","BAZ2BA-x434",
            ]

    ##

    if True: 

        loader = MultiDatasetDataloader(
            data_dirs = str(tmp_path / "data" / "*"),
            pdb_style = "*.dimple.pdb",
            mtz_style = "*.dimple.mtz",
            pdb_regex = None,
            mtz_regex = None,
            dir_regex = None,
            only_datasets = ["BAZ2BA-x431","BAZ2BA-x433"],
            ignore_datasets = None,
            dataset_prefix = None,
            )

        mcd = loader()
        datasets = mcd.datasets

        assert len(datasets) == 2

        assert sorted(datasets.keys()) == [
            "BAZ2BA-x431","BAZ2BA-x433",
            ]

        ##

        loader = MultiDatasetDataloader(
            data_dirs = str(tmp_path / "data" / "*"),
            pdb_style = "*.dimple.pdb",
            mtz_style = "*.dimple.mtz",
            pdb_regex = None,
            mtz_regex = None,
            dir_regex = None,
            only_datasets = None,
            ignore_datasets = ["BAZ2BA-x431","BAZ2BA-x433"],
            dataset_prefix = None,
            )

        mcd = loader()
        datasets = mcd.datasets

        assert len(datasets) == 3

        assert sorted(datasets.keys()) == [
            "BAZ2BA-x430","BAZ2BA-x432","BAZ2BA-x434",
            ]

        ##

        loader = MultiDatasetDataloader(
            data_dirs = str(tmp_path / "data" / "*"),
            pdb_style = "*.dimple.pdb",
            mtz_style = "*.dimple.mtz",
            pdb_regex = None,
            mtz_regex = None,
            dir_regex = None,
            only_datasets = ["BAZ2BA-x431","BAZ2BA-x432"],
            ignore_datasets = ["BAZ2BA-x431","BAZ2BA-x433"],
            dataset_prefix = None,
            )

        mcd = loader()
        datasets = mcd.datasets

        assert len(datasets) == 1

        assert sorted(datasets.keys()) == [
            "BAZ2BA-x432",
            ]

    ##
    
    for regexes in [
        (
            "BAZ2BA-(.*).dimple.pdb", 
            None, 
            None,
            ),
        (
            None, 
            "BAZ2BA-(.*).dimple.mtz", 
            None,
            ),
        (
            None, 
            None, 
            "BAZ2BA-(.*)",
            ),
        (
            "(.*).dimple.pdb", 
            "(.*).dimple.mtz",    
            "BAZ2BA-(.*)"), # this should be preferred
        (
            "BAZ2BA-(.*).dimple.pdb", # this should be preferred
            "(.*).dimple.mtz", 
            None,
            ),
        ]:

        loader = MultiDatasetDataloader(
            data_dirs = str(tmp_path / "data" / "*"),
            pdb_style = "*.dimple.pdb",
            mtz_style = "*.dimple.mtz",
            pdb_regex = regexes[0],
            mtz_regex = regexes[1],
            dir_regex = regexes[2],
            only_datasets = None,
            ignore_datasets = None,
            dataset_prefix = None,
            )

        mcd = loader()
        datasets = mcd.datasets

        assert len(datasets) == 5

        assert sorted(datasets.keys()) == [
            "x430","x431","x432","x433","x434",
            ]

    ##

    if True: 

        loader = MultiDatasetDataloader(
            data_dirs = str(tmp_path / "data" / "*"),
            pdb_style = "*.dimple.pdb",
            mtz_style = "*.dimple.mtz",
            pdb_regex = None,
            mtz_regex = None,
            dir_regex = None,
            only_datasets = None,
            ignore_datasets = None,
            dataset_prefix = "prefix",
            )

        mcd = loader()
        datasets = mcd.datasets

        assert len(datasets) == 5

        assert sorted(datasets.keys()) == [
            "prefixBAZ2BA-x430",
            "prefixBAZ2BA-x431",
            "prefixBAZ2BA-x432",
            "prefixBAZ2BA-x433",
            "prefixBAZ2BA-x434",
            ]

    
    

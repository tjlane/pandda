import os, copy
import pytest

def get_test_dataset_files(output_directory, n=1, dataset='BAZ2B'):

    from giant.resources.test_data import TEST_DATA

    TEST_DATA.extract_to_temporary_directory(
        dataset = dataset,
        tmp_dir = str(output_directory),
        n = n,
        )

    from giant.paths import resolve_glob
    pdb_files = resolve_glob(
        os.path.join(output_directory, 'data', '*', '*.dimple.pdb'),
        n = n,
    )

    mtz_files = [os.path.splitext(f)[0]+'.mtz' for f in pdb_files]

    return pdb_files, mtz_files

@pytest.fixture(scope="module")
def five_baz2b_test_datasets(tmp_path_factory):

    tmp_dir = tmp_path_factory.mktemp('datasets')

    pdb_files, mtz_files = get_test_dataset_files(str(tmp_dir), n=5, dataset='BAZ2B')

    from giant.mulch.dataset import CrystallographicDataset

    datasets = [
        CrystallographicDataset.from_file(
            model_filename = p,
            data_filename = m,
            )
        for p,m in zip(pdb_files, mtz_files)
        ]

    return datasets

@pytest.fixture
def five_baz2b_test_datasets_mcd(five_baz2b_test_datasets):

    from giant.mulch.collection import MultiCrystalDataset

    label = lambda d: os.path.basename(d.model.filename)

    mcd = MultiCrystalDataset({
        label(d) : copy.deepcopy(d) # copy to avoid contamination 
        for d in five_baz2b_test_datasets
        })

    return mcd

@pytest.fixture
def five_baz2b_test_datasets_mcd_labelled(five_baz2b_test_datasets):

    from giant.mulch.collection import MultiCrystalDataset

    label = lambda d: os.path.basename(d.model.filename)

    mcd = MultiCrystalDataset({
        label(d) : copy.deepcopy(d).label(tag=label(d)) # datasets are not automatically labelled! 
        for d in five_baz2b_test_datasets
        })

    return mcd

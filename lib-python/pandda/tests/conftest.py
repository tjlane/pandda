
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

# run once per module, then check output in different tests
@pytest.fixture(scope="module") 
def pandda_baz2b_test_data(tmp_path_factory):

    tmp_dir = tmp_path_factory.mktemp('input_data')

    pdb_files, mtz_files = get_test_dataset_files(str(tmp_dir), n=10, dataset='BAZ2B')

    pandda_data_dir = (tmp_dir / 'data')
    assert pandda_data_dir.exists()

    pandda_data_dir_star = (pandda_data_dir / '*')
    
    pandda_output_dir = (tmp_dir / 'pandda')
    assert not pandda_output_dir.exists()

    return {
        'work_dir' : str(tmp_dir),
        'out_dir' : str(pandda_output_dir),
        'data_dirs' : str(pandda_data_dir_star), 
        'pdb_style' : '*.dimple.pdb',
        'mtz_style' : '*.dimple.mtz',
        'pdb_files' : pdb_files,
        'mtz_files' : mtz_files,
    }

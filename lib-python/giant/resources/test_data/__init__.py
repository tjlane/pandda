import giant.logs as lg
logger = lg.getLogger(__name__)

import os, glob, shutil

from giant.paths import decompress_all_files_in_directory

from libtbx import group_args

_here = os.path.abspath(os.path.dirname(__file__))

class _TestData(object):
    datasets = {
        'BAZ2B' : os.path.join(_here, 'baz2b_test_data'),
    }

    def __init__(self):
        pass

    def extract_to_temporary_directory(self, dataset, tmp_dir=None, n=None):
        """Copy a test dataset to a temporary directory -- return path of temporary directory"""
        if (n is not None):
            assert n > 0
        data_path = self.get_path(dataset)
        if tmp_dir is None:
            import tempfile
            tmp_dir = tempfile.mkdtemp(prefix='test-data')
        assert os.path.exists(tmp_dir)
        target_dir = os.path.join(tmp_dir, 'data')
        os.mkdir(target_dir)
        for i, in_d in enumerate(
                sorted(
                    glob.glob(
                        os.path.join(data_path,'*/')
                    )
                )
            ):
            if (n is not None) and (i >= n):
                break
            out_d = os.path.join(target_dir, os.path.basename(in_d.strip('/')))
            shutil.copytree(in_d, out_d)
        # Unzip directory
        decompress_all_files_in_directory(path=target_dir)
        return tmp_dir

    def get_path(self, dataset):
        return self.datasets[dataset]

TEST_DATA = _TestData()

def get_test_data(output_directory, n=1, dataset='BAZ2B'):

    TEST_DATA.extract_to_temporary_directory(
        dataset = dataset,
        tmp_dir = str(output_directory),
        n = n,
        )

    from giant.paths import resolve_glob
    pdb_files = resolve_glob(
        os.path.join(str(output_directory), 'data', '*', '*.dimple.pdb'),
        n = n,
    )

    return pdb_files

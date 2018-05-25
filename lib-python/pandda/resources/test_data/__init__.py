import os, tempfile, shutil

from bamboo.common.logs import Log
from bamboo.common.file import decompress_all_files_in_directory

from libtbx import group_args

_here = os.path.abspath(os.path.dirname(__file__))

class _TestData:
    __hash__ = {'BAZ2B_TEST_DATA' : os.path.join(_here, 'baz2b_test_data')}

    def __init__(self):
        keys = self.__hash__.keys()
        self.keys = group_args(**dict(zip(map(str,keys), keys)))

    def extract_to_temporary_directory(self, choice, verbose=True):
        """Copy a test dataset to a temporary directory -- return path of temporary directory"""
        data_path = self.get_path(choice)
        tmp_dir = tempfile.mkdtemp(prefix='pandda-test-')
        assert os.path.exists(tmp_dir)
        target_dir = os.path.join(tmp_dir, 'data')
        shutil.copytree(data_path, target_dir)
        # Unzip directory
        decompress_all_files_in_directory(path=target_dir)
        return tmp_dir

    def get_path(self, choice):
        assert self.__hash__.has_key(choice), 'invalid choice'
        return self.__hash__[choice]

TEST_DATA = _TestData()

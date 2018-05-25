import os, shutil

import pytest

from libtbx.utils import Failure, Sorry

from pandda.analyse import pandda_phil, run
from pandda.resources.test_data import TEST_DATA

class TestBadInputData(object):

    def edit(self, filename, old_line, new_line):
        assert os.path.exists(filename)
        contents = open(filename, 'r').read()
        assert old_line in contents
        contents = contents.replace(old_line, new_line)
        with open(filename, 'w') as fh:
            fh.write(contents)

    def run(self, tmp_dir):
        os.chdir(tmp_dir)
        # Input options
        options_phil = """
        data_dirs='./data/*'
        pdb_style='*.dimple.pdb'
        high_res_upper_limit=2.0
        high_res_lower_limit=2.0
        dynamic_res_limits=False
        min_build_datasets=1
        """
        # Parse input args
        cmd_interpr = pandda_phil.command_line_argument_interpreter()
        working_phil = pandda_phil.fetch(sources=[cmd_interpr.process(options_phil)])
        # Make pandda object and run setup
        pandda = run(params=working_phil.extract())

    def test_no_cryst_line(self, capsys):
        cwd = os.getcwd()
        tmp_dir = TEST_DATA.extract_to_temporary_directory(choice=TEST_DATA.keys.BAZ2B_TEST_DATA)
        self.edit(filename=os.path.join(tmp_dir, 'data', 'BAZ2BA-x430', 'BAZ2BA-x430.dimple.pdb'),
                  old_line="CRYST1   82.128   96.768   57.955  90.00  90.00  90.00 C 2 2 21                 \n",
                  new_line="")
        with pytest.raises(SystemExit) as err_info:
            self.run(tmp_dir)
        captured = capsys.readouterr()
        assert "Dataset BAZ2BA-x430: There is no crystal symmetry for this structure -- does the input PDB file have a valid CRYST line?" in captured.out
        assert "Failure: 1 datasets had problems during loading. Error messages printed above." in captured.out
        os.chdir(cwd)
        shutil.rmtree(tmp_dir)

    def test_no_unit_cell(self, capsys):
        cwd = os.getcwd()
        tmp_dir = TEST_DATA.extract_to_temporary_directory(choice=TEST_DATA.keys.BAZ2B_TEST_DATA)
        self.edit(filename=os.path.join(tmp_dir, 'data', 'BAZ2BA-x430', 'BAZ2BA-x430.dimple.pdb'),
                  old_line="CRYST1   82.128   96.768   57.955  90.00  90.00  90.00 C 2 2 21                 \n",
                  new_line="CRYST1                                                 C 2 2 21                 \n")
        with pytest.raises(SystemExit) as err_info:
            self.run(tmp_dir)
        captured = capsys.readouterr()
        assert "Dataset BAZ2BA-x430: There is no unit cell information for this structure -- does the input PDB file have a valid CRYST line?" in captured.out
        assert "Failure: 1 datasets had problems during loading. Error messages printed above." in captured.out
        os.chdir(cwd)
        shutil.rmtree(tmp_dir)

    def test_no_space_group(self, capsys):
        cwd = os.getcwd()
        tmp_dir = TEST_DATA.extract_to_temporary_directory(choice=TEST_DATA.keys.BAZ2B_TEST_DATA)
        self.edit(filename=os.path.join(tmp_dir, 'data', 'BAZ2BA-x430', 'BAZ2BA-x430.dimple.pdb'),
                  old_line="CRYST1   82.128   96.768   57.955  90.00  90.00  90.00 C 2 2 21                 \n",
                  new_line="CRYST1   82.128   96.768   57.955  90.00  90.00  90.00                          \n")
        with pytest.raises(SystemExit):
            self.run(tmp_dir)
        captured = capsys.readouterr()
        assert "Dataset BAZ2BA-x430: There is no space group information for this structure -- does the input PDB file have a valid CRYST line?" in captured.out
        assert "Failure: 1 datasets had problems during loading. Error messages printed above." in captured.out
        os.chdir(cwd)
        shutil.rmtree(tmp_dir)


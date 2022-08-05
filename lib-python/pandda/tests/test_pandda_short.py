import os 
import pytest 

options_phil_template = """

data_dirs={data_dirs}
pdb_style={pdb_style}
out_dir={out_dir}
cpus={cpus}

high_res_upper_limit=2.0
high_res_lower_limit=2.0

grid_spacing=1.0
min_build_datasets=1

mask_selection_string="chain B and resname EDO and resseq 1"
outer_mask=10

ignore_datasets='BAZ2BA-x432'
not_train='BAZ2BA-x434'
not_test='BAZ2BA-x430,BAZ2BA-x431'
train='BAZ2BA-x430,BAZ2BA-x433,BAZ2BA-x435,BAZ2BA-x436,BAZ2BA-x437,BAZ2BA-x439'

"""
 
def run_standard_pandda(phil_str):

    from pandda.phil import pandda_phil
    from pandda.config import Config
    from pandda.programs.analyse import standard_pandda

    # Parse input args
    cmd_interpr = pandda_phil.command_line_argument_interpreter()
    working_phil = pandda_phil.fetch(
        sources=[cmd_interpr.process(phil_str)]
        )

    # Make pandda args object
    pandda_config = Config(working_phil.extract())

    # Run
    standard_pandda(pandda_config)


@pytest.mark.slow
class TestStandardPandda(object):

    @pytest.mark.dependency(name="test_standard_pandda")
    def test_standard_pandda(self, pandda_baz2b_test_data):

        import multiprocessing 

        formatted_phil = options_phil_template.format(
            cpus = multiprocessing.cpu_count(),
            **pandda_baz2b_test_data
            )

        run_standard_pandda(
            phil_str = formatted_phil,
            )

    @pytest.mark.dependency(depends=["test_standard_pandda"])
    def test_validate_event_csv(self, pandda_baz2b_test_data):
        pass

    @pytest.mark.dependency(depends=["test_standard_pandda"])
    def test_new(self, pandda_baz2b_test_data):
        pass

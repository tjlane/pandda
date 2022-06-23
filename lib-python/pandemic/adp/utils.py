import giant.logs as lg
logger = lg.getLogger(__name__)

import os, copy, shutil, collections
import pathlib as pl
from libtbx import adopt_init_args
from giant.paths import easy_directory


class TaskMultipleRunWrapper(object):
    """Class for calling a task multiple times with different output directories"""

    def __init__(self,
        output_directory,
        task,
        task_output_directory_prefix,
        ):

        output_directory = easy_directory(output_directory)
        task = copy.deepcopy(task)

        adopt_init_args(self, locals())

        self.previous_runs = collections.OrderedDict()
        self.counter = 1

    def __call__(self,
        output_directory_suffix = '',
        remove_previous = False,
        **task_args
        ):

        # Remove previous output
        if remove_previous is True:
            self.delete_output()

        # Create new output directory
        new_o_dir = os.path.join(self.output_directory, self.task_output_directory_prefix+output_directory_suffix)
        assert not os.path.exists(new_o_dir), 'output directory already exists: {}'.format(new_o_dir)

        # Create copy of the task and change the output directory
        task_copy = copy.deepcopy(self.task)
        task_copy.setup(output_directory = pl.Path(new_o_dir))

        # Run the task
        task_copy(**task_args)

        self.previous_runs[self.counter] = new_o_dir

        self.counter += 1

        return task_copy

    def delete_output(self):
        for o_dir in self.previous_runs.values():
            if os.path.exists(o_dir):
                shutil.rmtree(o_dir)


def show_file_dict(self, file_dict, indent=0):
    s = '  '
    for k, v in file_dict.items():
        if isinstance(v, dict):
            logger(s*indent + '> {}'.format(k))
            show_file_dict(self, v, indent+1)
        elif isinstance(v, str):
            logger(s*indent + '> {}: {}'.format(k, v))
        else:
            logger(s*indent + '> {}'.format(k))
            try:
                for vv in v:
                    logger(s*(indent+1)+str(vv))
            except:
                logger(s*(indent+1)+str(v))

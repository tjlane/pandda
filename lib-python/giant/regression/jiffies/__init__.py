import os

def check_files_exist(file_list):
    for f in file_list:
        assert os.path.exists(f), f

def run_module(args, module):

    from giant.jiffies import run_default
    run_default(
        run                 = module.run,
        master_phil         = module.master_phil,
        args                = args,
        blank_arg_prepend   = module.blank_arg_prepend,
        program             = module.PROGRAM,
        description         = module.DESCRIPTION,
    )

def run_program(args, program):

    from giant.dispatcher import Dispatcher

    prog = Dispatcher(program)
    prog.extend_args(args)
    prog.run()

def run_jiffy(args, module):

    prog_name = module.__name__.replace('.jiffies', '')

    run_program(
        args = args,
        program = prog_name,
    )

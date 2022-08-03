import os

from giant.jiffies import (
    parse_phil_args,
    )

# def get_params_from_phil_and_args(master_phil, args):

#     cmd_interpr = master_phil.command_line_argument_interpreter()

#     arg_sources = [
#         cmd_interpr.process(arg=a) for a in args
#         ]

#     working_phil = master_phil.fetch(
#         sources = (
#             arg_sources
#             ),
#         )

#     params = working_phil.extract()

#     return params

def run_module(module, args):

    params = parse_phil_args(
        master_phil = module.master_phil,
        args = args,
        blank_arg_prepend = (
            module.blank_arg_prepend
            if hasattr(module, 'blank_arg_prepend')
            else None
            ),
        ).extract()

    return module.run(params)

def run_program(program, args):
    """
    Run program from the command line with args
    """

    from giant.dispatcher import Dispatcher

    prog = Dispatcher(program)
    prog.extend_args(args)
    prog.run()

    return prog

def run_jiffy(args, module):
    """
    Run a module as command line argument, assuming that
    the name of the program is derivable from the module
    name.
    """

    prog_name = (
        module.__name__.replace('.jiffies', '')
        )

    run_program(
        program = prog_name,
        args = args,
    )

###

def check_files(files_list, check_list, error_remaining=True):

    # files_list = map(str, files_list)
    # check_list = map(str, check_list)

    for f in check_list:
        try:
            files_list.remove(
                str(f)
                )
        except ValueError:
            raise Exception(
                'expected file not present: {f}'.format(
                    f = str(f),
                    )
                )

    if error_remaining is True:

        if len(files_list) > 0:

            raise Exception(
                'Unexpected output files: {}'.format(
                    '\n'.join(files_list),
                    ),
                )

###

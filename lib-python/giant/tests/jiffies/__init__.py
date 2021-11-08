#

def get_params_from_phil_and_args(master_phil, args):

    cmd_interpr = master_phil.command_line_argument_interpreter()

    arg_sources = [
        cmd_interpr.process(arg=a) for a in args
        ]

    working_phil = master_phil.fetch(
        sources = (
            arg_sources
            ),
        )

    params = working_phil.extract()

    return params

def check_files(files_list, check_list, error_remaining=True):

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


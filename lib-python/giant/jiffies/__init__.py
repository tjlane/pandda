import os, sys, copy

from libtbx.utils import Sorry

def parse_phil_args(master_phil, args, home_scope=None, blank_arg_prepend=None):

    if blank_arg_prepend is None:
        pass
    elif isinstance(blank_arg_prepend, dict):
        for item in blank_arg_prepend.values():
            assert '=' in item
    elif isinstance(blank_arg_prepend, str):
        assert '=' in blank_arg_prepend
    else:
        raise Exception('blank_arg_prepend must be str or dict')

    # Copy the args so that we can remove items from the list without affecting args etc
    args = copy.copy(args)
    # Construct interpreter
    cmd_interpr = master_phil.command_line_argument_interpreter(home_scope=home_scope)
    # Process input arguments
    phil_objects = []
    for arg in args:
        try:
            # Prepend if blank
            if '=' not in arg:
                if isinstance(blank_arg_prepend, dict):
                    found_key = False
                    for key in blank_arg_prepend.keys():
                        if key is None:
                            continue
                        if arg.endswith(key):
                            arg = blank_arg_prepend[key]+arg
                            found_key = True
                            break
                    if (found_key == False) and (None in blank_arg_prepend.keys()):
                        arg = blank_arg_prepend[None]+arg
                elif isinstance(blank_arg_prepend, str):
                    arg = blank_arg_prepend+arg
            # Attempt to process arg
            cmd_line_args = cmd_interpr.process(arg=arg)
        except KeyboardInterrupt:
            raise
        except Exception:
            raise Sorry("Unknown file or keyword: %s" % arg)
        else:
            phil_objects.append(cmd_line_args)
    # Extract Scope object
    working_phil = master_phil.fetch(sources=phil_objects)

    return working_phil

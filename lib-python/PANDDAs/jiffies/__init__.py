import os, sys, copy

from libtbx.utils import Sorry

def parse_phil_args(master_phil, args, home_scope=None):

    # Copy the args so that we can remove items from the list without affecting args etc
    args = copy.copy(args)
    # Construct interpreter
    cmd_interpr = master_phil.command_line_argument_interpreter(home_scope=home_scope)
    # Process input arguments
    phil_objects = []
    for arg in args:
        try:
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

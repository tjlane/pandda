import os, sys, copy

import libtbx.phil
from libtbx.utils import Sorry

def show_defaults_and_exit_maybe(master_phil, args):
    """Show master_phil and exit if requested"""

    # Show Defaults (just values)
    if '--show-defaults' in args:
        print('=========>\nShowing Default Parameters\n=========>')
        master_phil.show(expert_level=0, attributes_level=0)
    # Show Defaults (including help)
    elif '--help' in args:
        print('=========>\nShowing Default Parameters & Help\n=========>')
        master_phil.show(expert_level=0, attributes_level=1)
    # Show Defaults (including help and types)
    elif '--help-and-types' in args:
        print('=========>\nShowing Default Parameters & Help+Types\n=========>')
        master_phil.show(expert_level=0, attributes_level=2)
    # Show Defaults (everything)
    elif '--expert' in args:
        print('=========>\nShowing Default Parameters (Expert Level)\n=========>')
        master_phil.show(expert_level=3, attributes_level=3)
    else:
        return

    raise SystemExit('=========>\nNow Exiting.\n=========>')

def parse_phil_args(master_phil, args, blank_arg_prepend=None, home_scope=None):

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

    # Process any args that are eff files
    eff_files = [f for f in args if ((f.endswith('.eff') or f.endswith('.def')) and (not f.count('=')) and os.path.isfile(f))]
    # Remove them from the original lists
    [args.remove(f) for f in eff_files]
    # Parse the 'eff' files - these should contain phils
    eff_sources = [libtbx.phil.parse(open(f, 'r').read()) for f in eff_files]

    # Process input arguments
    arg_sources = []
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
            arg_sources.append(cmd_line_args)
    # Extract Scope object (putting eff sources first so that they're overridden if double-defined)
    working_phil = master_phil.fetch(sources=eff_sources+arg_sources)

    return working_phil

def extract_params_default(master_phil, args, blank_arg_prepend=None, home_scope=None):
    """Extract the parameters by a default script"""
    show_defaults_and_exit_maybe(master_phil, args)
    working_phil = parse_phil_args(master_phil, args, blank_arg_prepend=blank_arg_prepend, home_scope=home_scope)
    return working_phil

def run_default(run, master_phil, args, blank_arg_prepend=None):
    """Run a program via a standard setup of functions and objects"""
    working_phil = extract_params_default(master_phil=master_phil, args=args, blank_arg_prepend=blank_arg_prepend)
    out = run(params=working_phil.extract())
    return 0

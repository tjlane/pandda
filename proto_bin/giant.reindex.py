#!/usr/bin/env pandda.python

import os, sys, copy, re

import libtbx.phil
from Bamboo.Common.command import commandManager
from Giant.jiffies import parse_phil_args
from Giant.Xray.Data import crystalSummary

blank_arg_prepend = 'mtz='

master_phil = libtbx.phil.parse("""
input {
    mtz = None
        .type = path
        .multiple = True
    reference_mtz = None
        .type = path
}
options {
    free_r_flag = FreeR_flag
        .type = str
}
""")

def run(params):

    assert params.input.reference_mtz is not None

    for mtz_file in params.input.mtz:

        reind_output_mtz = os.path.join(os.path.dirname(mtz_file), 'reind_'+os.path.basename(mtz_file))
        cad_output_mtz   = os.path.join(os.path.dirname(mtz_file), 'cad_'+os.path.basename(mtz_file))
        free_output_mtz  = os.path.join(os.path.dirname(mtz_file), 'free_'+os.path.basename(mtz_file))

        cm = commandManager('pointless')
        cm.add_command_line_arguments([ 'hklin', mtz_file,
                                        'hklref', params.input.reference_mtz,
                                        'hklout', reind_output_mtz   ])
        cm.add_standard_input([ 'tolerance 5' ])
        cm.print_settings()
        cm.run()
        print '============================>'
        print cm.output
        print '============================>'
        print cm.error
        print '============================>'

        cm = commandManager('cad')
        cm.add_command_line_arguments([ 'hklin1', reind_output_mtz,
                                        'hklin2', params.input.reference_mtz,
                                        'hklout', cad_output_mtz   ])
        cm.add_standard_input([ 'labin file_number 1 E1=IMEAN E2=SIGIMEAN',
                                'labin file_number 2 E1={}'.format(params.options.free_r_flag),
                                'labout file_number 2 E1=FreeR_flag',
                                'END'     ])
        cm.print_settings()
        cm.run()
        print '============================>'
        print cm.output
        print '============================>'
        print cm.error
        print '============================>'

        cm = commandManager('freerflag')
        cm.add_command_line_arguments([ 'hklin', cad_output_mtz,
                                        'hklout', free_output_mtz   ])
        cm.add_standard_input([ 'COMPLETE FREE=FreeR_flag',
                                'END'     ])
        cm.print_settings()
        cm.run()
        print '============================>'
        print cm.output
        print '============================>'
        print cm.error
        print '============================>'



if __name__ == '__main__':

    # Show Defaults (just values)
    if '--show-defaults' in sys.argv:
        master_phil.show(attributes_level=0)
    # Show Defaults (including information)
    elif '--help' in sys.argv:
        master_phil.show(attributes_level=2)
    # ... or just run ...
    elif '--expert' in sys.argv:
        master_phil.show(attributes_level=4)
    # ... or just run ...
    else:
        working_phil = parse_phil_args(master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)
        out = run(params=working_phil.extract())
    # Exit (unnecessary, but eh)
    sys.exit()


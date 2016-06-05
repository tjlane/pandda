import os, sys, copy, re

import libtbx.phil
from bamboo.common.command import CommandManager
from giant.xray.data import CrystalSummary

#######################################

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

#######################################

def run(params):

    assert params.input.reference_mtz is not None

    for mtz_file in params.input.mtz:

        reind_output_mtz = os.path.join(os.path.dirname(mtz_file), 'reind_'+os.path.basename(mtz_file))
        cad_output_mtz   = os.path.join(os.path.dirname(mtz_file), 'cad_'+os.path.basename(mtz_file))
        free_output_mtz  = os.path.join(os.path.dirname(mtz_file), 'free_'+os.path.basename(mtz_file))

        ######################################
        # REINDEX
        ######################################
        cm = CommandManager('pointless')
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

        ######################################
        # Transfer Free-R Flags
        ######################################
        cm = CommandManager('cad')
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

        ######################################
        # Fill in missing reflections (nans)
        ######################################
        cm = CommandManager('uniqueify')
        cm.add_command_line_arguments([ '-f', 'FreeR_flag',
                                        cad_output_mtz,
                                        free_output_mtz   ])
        cm.print_settings()
        cm.run()
        print '============================>'
        print cm.output
        print '============================>'
        print cm.error
        print '============================>'

#        ######################################
#        # Complete Free-R Flags
#        ######################################
#        cm = CommandManager('freerflag')
#        cm.add_command_line_arguments([ 'hklin', cad_output_mtz,
#                                        'hklout', free_output_mtz   ])
#        cm.add_standard_input([ 'COMPLETE FREE=FreeR_flag',
#                                'END'     ])
#        cm.print_settings()
#        cm.run()
#        print '============================>'
#        print cm.output
#        print '============================>'
#        print cm.error
#        print '============================>'

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)

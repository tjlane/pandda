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
    columns = F,SIGF
        .type = str
}
output {
    column = IMEAN
        .type = str
    file_prefix = intensities_
        .type = str
}
""")

def run(params):

    for mtz_file in params.input.mtz:

        output_mtz = os.path.join(os.path.dirname(mtz_file), params.output.file_prefix+os.path.basename(mtz_file))

        cs = crystalSummary.from_mtz(mtz_file=mtz_file)

        if params.output.column in cs.column_labels:
            print 'OUTPUT COLUMN ALREADY PRESENT IN MTZ_FILE: {}'.format(mtz_file)
            os.symlink(mtz_file, output_mtz)
            continue

        cm = commandManager('phenix.reflection_file_converter')
        cm.add_command_line_arguments([ mtz_file,
                                        '--mtz-root-label='+params.output.column,
                                        '--label='+params.input.columns,
                                        '--write-mtz-intensities',
                                        '--mtz='+output_mtz
                                    ])
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


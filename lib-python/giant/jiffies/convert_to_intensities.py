import os, sys

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
    columns = F,SIGF
        .type = str
}
method {
    package = *phenix ccp4
        .type = choice
}
output {
    column = IMEAN
        .type = str
    file_prefix = intensities_
        .type = str
}
""")

#######################################

def run(params):

    input_cols = params.input.columns.split(',')

    for mtz_file in params.input.mtz:

        output_mtz = os.path.join(os.path.dirname(mtz_file), params.output.file_prefix+os.path.basename(mtz_file))

        cs = CrystalSummary.from_mtz(mtz_file=mtz_file)

        if params.output.column in cs.column_labels:
            print 'OUTPUT COLUMN ALREADY PRESENT IN MTZ_FILE: {}'.format(mtz_file)
            if not os.path.exists(output_mtz):
                os.symlink(os.path.basename(mtz_file), output_mtz)
            continue

        cols_correct  = True
        for col in input_cols:
            if col not in cs.column_labels:
                cols_correct = False
                break
        if cols_correct == False:
            print 'Column not present in mtz_file:'
            print cs.column_labels
            print 'SKIPPING'
            continue

        if params.method.package == 'phenix':
            cm = CommandManager('phenix.reflection_file_converter')
            cm.add_command_line_arguments([ mtz_file,
                                            '--mtz-root-label='+params.output.column,
                                            '--label='+params.input.columns,
                                            '--non-anomalous',
                                            '--write-mtz-intensities',
                                            '--mtz='+output_mtz     ])
        elif params.method.package == 'ccp4':
            raise Exception('NOT IMPLEMENTED')

        cm.print_settings()
        cm.run()
        print '============================>'
        print cm.output
        print '============================>'
        print cm.error
        print '============================>'

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)

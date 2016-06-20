import os, sys, copy, re, shutil

import libtbx.phil
import libtbx.easy_mp
from bamboo.common.command import CommandManager
from giant.xray.data import CrystalSummary

#######################################

blank_arg_prepend = {'.pdb':'input.pdb=', '.mtz':'input.mtz='}

master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .type = path
    mtz = None
        .type = path
}
options {
    n_cycles = 0
        .type = int
    add_hydrogens = True
        .type = bool
    output_hydrogens = True
        .type = bool
    build_absent_atoms = True
        .type = bool
}
output {
    pdb = None
        .type = path
    mtz = None
        .type = path
    log = None
        .type = path
    suffix = '.refmac'
        .type = str
}
""")

#######################################

def run(params):

    in_pdb = params.input.pdb
    assert os.path.exists(in_pdb), 'PDB does not exist: {}'.format(in_pdb)
    in_mtz = params.input.mtz
    assert os.path.exists(in_mtz), 'MTZ does not exist: {}'.format(in_mtz)

    if params.output.pdb: out_pdb = params.output.pdb
    else:                 out_pdb = os.path.splitext(in_pdb)[0] + params.output.suffix + '.pdb'

    if params.output.mtz: out_mtz = params.output.mtz
    else:                 out_mtz = os.path.splitext(in_mtz)[0] + params.output.suffix + '.mtz'

    if params.output.log: out_log = params.output.log
    else:                 out_log = os.path.splitext(in_pdb)[0] + params.output.suffix + '.log'

    ###########################################
    # COMMAND LINE COMMANDS
    ###########################################
    cm = CommandManager('refmac5')
    cm.add_command_line_arguments(['xyzin',in_pdb,'xyzout',out_pdb,'hklin',in_mtz,'hklout',out_mtz])

    ###########################################
    # MAIN PARAMETERS
    ###########################################
    cm.add_standard_input('MAKE NCYC {}'.format(params.options.n_cycles))

    ###########################################
    # HYDROGENS
    ###########################################
    if params.options.add_hydrogens:
        cm.add_standard_input('MAKE HYDR A')
    else:
        cm.add_standard_input('MAKE HYDR Y')

    if params.options.output_hydrogens:
        cm.add_standard_input('MAKE HOUT Y')
    else:
        cm.add_standard_input('MAKE HOUT N')

    ###########################################
    # MISSING ATOMS
    ###########################################
    if params.options.build_absent_atoms:
        cm.add_standard_input('MAKE BUIL Y')
    else:
        cm.add_standard_input('MAKE BUIL N')

    ###########################################
    # LASTLY
    ###########################################
    cm.add_standard_input('END')

    ###########################################
    # RUN
    ###########################################
    print 'Running Refmac:'
    cm.print_settings()
    ret_code = cm.run()
    if ret_code != 0:
        print '============================>'
        print 'Refmac returned with an error'
        print '============================>'
        print cm.output
        print '============================>'
        print cm.error
        print '============================>'

    with open(out_log, 'a') as fh:
        fh.write('============================>'+'\n')
        fh.write('Refmac Log File'+'\n')
        fh.write('============================>'+'\n')
        fh.write('> STDOUT'+'\n')
        fh.write('============================>'+'\n')
        fh.write(cm.output+'\n')
        fh.write('============================>'+'\n')
        fh.write('> STDERR'+'\n')
        fh.write('============================>'+'\n')
        fh.write(cm.error+'\n')
        fh.write('============================>'+'\n')

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)
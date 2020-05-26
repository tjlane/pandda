import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys, copy, re, shutil

from giant.exceptions import Sorry, Failure

#######################################

blank_arg_prepend = {
    '.pdb':'input.pdb=',
    '.mtz':'input.mtz=',
}

import libtbx.phil
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
    in_mtz = params.input.mtz

    if (in_pdb is None):
        raise IOError('No PDB file has been provided')

    if not os.path.exists(in_pdb):
        raise IOError('PDB does not exist: {}'.format(in_pdb))

    if (in_mtz is not None) and not os.path.exists(in_mtz):
        raise IOError('MTZ does not exist: {}'.format(in_mtz))

    if (params.output.pdb is not None):
        out_pdb = params.output.pdb
    else:
        out_pdb = os.path.splitext(in_pdb)[0] + params.output.suffix + '.pdb'

    if (params.output.mtz is not None):
        out_mtz = params.output.mtz
    elif in_mtz:
        out_mtz = os.path.splitext(in_mtz)[0] + params.output.suffix + '.mtz'
    else:
        out_mtz = None

    if (params.output.log is not None):
        out_log = params.output.log
    else:
        out_log = os.path.splitext(in_pdb)[0] + params.output.suffix + '.log'

    ###########################################
    # COMMAND LINE COMMANDS
    ###########################################
    from giant.dispatcher import Dispatcher
    prog = Dispatcher('refmac5')
    prog.extend_args(['xyzin',in_pdb])
    prog.extend_args(['xyzout',out_pdb])
    if (in_mtz is not None):
        prog.extend_args(['hklin',in_mtz])
        prog.extend_args(['hklout',out_mtz])

    ###########################################
    # MAIN PARAMETERS
    ###########################################
    prog.append_stdin('NCYC {}'.format(params.options.n_cycles))

    ###########################################
    # OPTIONS
    ###########################################
    if (params.options.add_hydrogens is True):
        prog.append_stdin('MAKE HYDR A')
    else:
        prog.append_stdin('MAKE HYDR Y')

    if (params.options.output_hydrogens is True):
        prog.append_stdin('MAKE HOUT Y')
    else:
        prog.append_stdin('MAKE HOUT N')

    if params.options.build_absent_atoms:
        prog.append_stdin('MAKE BUIL Y')
    else:
        prog.append_stdin('MAKE BUIL N')

    ###########################################
    # LASTLY
    ###########################################
    prog.append_stdin('END')

    ###########################################
    # RUN
    ###########################################
    logger.subheading('Running REFMAC')
    logger(prog.as_string())

    result = prog.run()

    prog.write_output(log_file=out_log)

    if result['exitcode'] != 0:
        raise Sorry('REFMAC returned an error -- see log file: {}'.format(out_log))

    logger.subheading('finished normally')

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
    )

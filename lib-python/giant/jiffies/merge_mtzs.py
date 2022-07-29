import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys

from giant.exceptions import Sorry

#######################################

blank_arg_prepend = {
    '.mtz':'input.mtz=',
    None:'column.label=',
}

import libtbx.phil
master_phil = libtbx.phil.parse("""
input {
    mtz = None
        .multiple = True
        .type = path
}
output {
    mtz = merged.mtz
        .type = path
    log = merged.log
        .type = path
}
options {
    label_suffix = *incremental filename foldername
        .help = 'suffix added to the column name in the output file'
        .type = choice
        .multiple = False
    column {
        label = None
            .help = 'columns to be extracted from the input mtz files'
            .type = str
            .multiple = True
    }
}
""")

#######################################

def run(params):

    # unpack params
    in_mtzs = params.input.mtz
    out_mtz = params.output.mtz
    out_log = params.output.log
    label_suffix = params.options.label_suffix

    # Setup logging
    logger = lg.setup_logging(
        name = __name__,
        log_file = out_log,
    )

    if len(in_mtzs) < 2:
        raise IOError('Need to provide at least two mtz files')

    # command line caller
    from giant.dispatcher import Dispatcher
    cm = Dispatcher('cad')

    # iterate through input mtzs
    for i_mtz, mtz in enumerate(in_mtzs):
        # Use numbering from 1
        n_mtz = i_mtz + 1
        # Create an mtz suffix
        if   label_suffix == 'incremental':
            suffix = '-{}'.format(n_mtz)
        elif label_suffix == 'filename':
            suffix = '-{}'.format(os.path.splitext(os.path.basename(mtz))[0])
        elif label_suffix == 'foldername':
            suffix = '-{}'.format(os.path.basename(os.path.dirname(mtz)))
        else:
            raise NotImplementedError('Not yet implemented, sorry')
        # Add to the command line
        cm.extend_args(['hklin{}'.format(n_mtz), mtz])
        # Select column labels
        cm.append_stdin(
            'labin file_number {0} {1}'.format(
                n_mtz,
                ' '.join([
                    'E{0}={1}'.format(i+1,c)
                    for i,c in enumerate(params.options.column.label)
                ])
            )
        )
        cm.append_stdin(
            'labout file_number {0} {1}'.format(
                n_mtz,
                ' '.join([
                    'E{0}={1}'.format(i+1,c+suffix)
                    for i,c in enumerate(params.options.column.label)
                ])
            )
        )

    # output files
    cm.extend_args(['hklout', out_mtz])

    # run
    logger.subheading('Running CAD')
    logger(cm.as_string())
    result = cm.run()
    cm.write_output(log_file=out_log)
    if result['exitcode'] != 0:
        logger(str(cm.result.stdout))
        logger(str(cm.result.stderr))
        raise Sorry('{} returned an error'.format(cm.progam))

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
    )

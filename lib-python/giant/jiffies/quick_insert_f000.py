import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys, copy, re, shutil

from giant.exceptions import Sorry, Failure

#######################################

blank_arg_prepend = {
    '.mtz':'input.mtz=',
    '.pdb':'input.pdb=',
    None:'column.label=',
}

import libtbx.phil
master_phil = libtbx.phil.parse("""
input {
    mtz = None
        .multiple = False
        .type = path
    pdb = None
        .multiple = False
        .type = path
}
options {
    column {
        label = 2FOFCWT,PH2FOFCWT
            .help = 'columns to be extracted from the input mtz files'
            .type = str
            .multiple = False
    }
}
output {
    suffix = '.with-F000.mtz'
        .type = str
}
""")

#######################################

def run(params):

    assert params.input.pdb and os.path.exists(params.input.pdb), 'Need to provide an pdb file'
    assert params.input.mtz and os.path.exists(params.input.mtz), 'Need to provide an mtz file'

    out_mtz = os.path.splitext(params.input.mtz)[0]+params.output.suffix
    assert not os.path.exists(out_mtz)
    out_log = out_mtz.replace('.mtz', '.log')

    ###########################################
    # Load diffraction data
    ###########################################
    import iotbx.mtz
    diff_data = iotbx.mtz.object(params.input.mtz)

    ###########################################
    # COMMAND LINE COMMANDS
    ###########################################
    prog = Dispatcher('phenix.f000')
    prog.add_command_line_arguments(params.input.pdb)

    ###########################################
    # RUN
    ###########################################
    logger(prog.as_string())
    result = prog.run()
    prog.write_output(log_file=out_log)
    if result['exitcode'] != 0:
        logger.subheading('phenix.f000 returned with an error')
        logger(prog.stdout)
        logger(prog.stderr)

    f000_regex = re.compile('Estimate of F\(0,0,0\)=(.*?) given mean bulk-solvent density (.*?) and fraction (.*?)\n')
    f000_details = f000_regex.findall(cm.output)

    assert len(f000_details) == 1, 'Error extracting f000 from output of phenix.f000'
    f000, solvent_density, solvent_fraction = list(map(float,f000_details[0]))
    logger('F000 is {}'.format(f000))

    ###########################################
    # Append to diffraction data
    ###########################################

    col_labs = params.options.column.label.split(',')
    logger.subheading('Extracting column labels: {}'.format(col_labs))

    miller_dict = diff_data.as_miller_arrays_dict()
    ampl_data = [d for l,d in miller_dict.items() if l[2] in col_labs]
    assert ampl_data, 'Matching AMPLITUDE DATA not in mtz data: {}'.format(ampl)
    assert len(ampl_data) == 2
    ampl_data = ampl_data[0]
    assert ampl_data.info().labels == col_labs

    # Create new objects
    new_data = ampl_data.customized_copy()
    # Append the F000 term
    new_data.indices().append((0,0,0))
    new_data.data().append(f000)
    # Sort for kicks
    new_data = new_data.sort('packed_indices')

    # Create new "dataset" from the miller array of the amplitudes
    mtz_dataset = new_data.as_mtz_dataset(column_root_label="2FOFCWT-ABS")
    # Convert "dataset" to an "object"
    mtz_object = mtz_dataset.mtz_object()
    mtz_object.write(file_name = out_mtz)

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
    )

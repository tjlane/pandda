import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys, shutil

from giant.exceptions import Sorry, Failure

from giant.paths import splice_ext, easy_directory
from giant.dispatcher import Dispatcher

############################################################################

PROGRAM = 'giant.datasets.prepare'

DESCRIPTION = """
    Take input MTZs, reindex, transfer R-free flags from a reference, fill missing reflections, and optionally run refinement.
"""

############################################################################

blank_arg_prepend = {
    '.mtz':'mtz=',
    '.pdb':'reference_pdb=',
}

import libtbx.phil
master_phil = libtbx.phil.parse("""
input {
    mtz = None
        .type = path
        .multiple = True
    reference_pdb = None
        .type = path
    reference_mtz = None
        .type = path
}
output {
    output_directory = None
        .type = str
    labelling = none *filename foldername
        .type = choice
    keep_intermediate_files = False
        .type = bool
}
actions = *reindex_to_reference \
          *fill_missing_reflections \
          *transfer_rfree_flags \
          *refinement_pipeline
    .type = choice(multi=True)
options {
    reindex {
        tolerance = 5
            .type  = float
    }
    rfree_flag {
        input = FreeR_flag
           .type = str
        output = None
            .type = str
    }
    missing_reflections {
        fill_high = 4.0
            .type = float
        fill_low = 9999.0
            .type = float
    }
    refinement {
        pipeline = *dimple
            .type = choice(multi=True)
        col_label = IMEAN
            .type = str
    }
}
settings {
    cpus = 1
        .help = "number of cpus to use for processing"
        .type = int
}
""")

############################################################################

def run_program(prog):

    logger(prog.as_string())

    result = prog.run()

    if result['exitcode'] != 0:

        logger.heading('{} exited with an error'.format(prog.program))
        logger.subheading('Stdout')
        logger(result.stdout)
        logger.subheading('Stderr')
        logger(result.stderr)

        raise Sorry('{} exited with an error'.format(prog.program))

############################################################################

def reindex_mtz_to_reference(in_mtz, out_mtz, reference_mtz, tolerance):
    """Reindex the data in one mtz to a reference mtz"""

    logger.subheading('Running reindexing')

    prog = Dispatcher('pointless')
    prog.extend_args([
        'hklin',  in_mtz,
        'hklref', reference_mtz,
        'hklout', out_mtz,
    ])
    prog.extend_args([
        'tolerance {}'.format(tolerance),
    ])

    run_program(prog)

    if not os.path.exists(out_mtz):
        logger(prog.result.stdout)
        logger(prog.result.stderr)
        raise Failure('reindexing has failed -- {} does not exist'.format(out_mtz))

    return prog.result

def fill_missing_reflections(in_mtz, out_mtz, fill_resolution_low, fill_resolution_high, delete_tmp_files=True):
    """Complete the set of miller indices in an MTZ file"""

    logger.subheading('Filling missing reflections')

    tmp_mtz_1 = splice_ext(path=out_mtz, new='step1-truncate')
    tmp_mtz_2 = splice_ext(path=out_mtz, new='step2-uniquify')
    tmp_mtz_3 = splice_ext(path=out_mtz, new='step3-remerged')

    # Stage 1 - truncate dataset, fill missing reflections, change column name
    logger('Running CAD...')
    prog = Dispatcher('cad')
    prog.extend_args([
        'hklin1', in_mtz,
        'hklout', tmp_mtz_1,
    ])
    prog.extend_stdin([
        'monitor BRIEF',
        'labin file_number 1 ALL',
        'resolution file 1 {} {}'.format(fill_resolution_low, fill_resolution_high),
    ])

    run_program(prog)

    if not os.path.exists(tmp_mtz_1):
        logger(prog.result.stdout)
        logger(prog.result.stderr)
        raise Failure('filling of missing reflections has failed -- {} does not exist'.format(tmp_mtz_1))

    # Stage 2 - Uniqueify the file
    logger('Running uniqueify...')
    prog = Dispatcher('uniqueify')
    prog.extend_args([
        '-p', '0.05',
        tmp_mtz_1,
        tmp_mtz_2,
    ])

    run_program(prog)

    if not os.path.exists(tmp_mtz_2):
        logger(prog.result.stdout)
        logger(prog.result.stderr)
        raise Failure('filling of missing reflections has failed -- {} does not exist'.format(tmp_mtz_2))

    # Stage 3 - remerge the two files
    logger('Running CAD...')
    prog = Dispatcher('cad')
    prog.extend_args([
        'hklin1', in_mtz,
        'hklin2', tmp_mtz_2,
        'hklout', tmp_mtz_3,
    ])
    prog.extend_stdin([
        'monitor BRIEF',
        'labin file_number 1 ALL',
        'labin file_number 2 E1=FreeR_flag',
        'labout file_number 2 E1=dummy',
    ])

    run_program(prog)

    if not os.path.exists(tmp_mtz_3):
        logger(prog.result.stdout)
        logger(prog.result.stderr)
        raise Failure('filling of missing reflections has failed -- {} does not exist'.format(tmp_mtz_3))

    # Stage 4 - remove the dummy column
    prog = Dispatcher('mtzutils')
    prog.extend_args([
        'hklin1', tmp_mtz_3,
        'hklout', out_mtz,
    ])
    prog.extend_stdin([
        'HEADER BRIEF',
        'EXCLUDE 1 dummy',
        'ONEFILE',
        'END',
    ])

    run_program(prog)

    if not os.path.exists(tmp_mtz_3):
        logger(prog.result.stdout)
        logger(prog.result.stderr)
        raise Failure('filling of missing reflections has failed -- {} does not exist'.format(out_mtz))

    if (delete_tmp_files is True):
        os.remove(tmp_mtz_1)
        os.remove(tmp_mtz_2)
        os.remove(tmp_mtz_3)

def transfer_rfree_flags(in_mtz, out_mtz, reference_mtz, input_free_r_flag, output_free_r_flag, delete_tmp_files=True):
    """Copy R-free flags from reference mtz"""

    tmp_mtz = splice_ext(path=out_mtz, new='step1-transfer')

    # Stage 1 - transfer R-free from reference
    logger.subheading('Transferring R-free flags')
    prog = Dispatcher('cad')
    prog.extend_args([
        'hklin1', in_mtz,
        'hklin2', reference_mtz,
        'hklout', tmp_mtz,
    ])
    prog.extend_stdin([
        'labin file_number 1 ALL',
        'labin file_number 2 E1={}'.format(input_free_r_flag),
        'labout file_number 2 E1={}'.format(output_free_r_flag),
        'END',
    ])

    run_program(prog)

    if not os.path.exists(tmp_mtz):
        logger(prog.result.stdout)
        logger(prog.result.stderr)
        raise Failure('transfer of R-free flags has failed -- {} does not exist'.format(tmp_mtz))

    # Stage 2 - populate missing R-free values
    logger.subheading('Completing R-free flags')
    prog = Dispatcher('freerflag')
    prog.extend_args([
        'hklin',  tmp_mtz,
        'hklout', out_mtz,
    ])
    prog.extend_stdin([
        'COMPLETE FREE={}'.format(output_free_r_flag),
        'END',
    ])

    run_program(prog)

    if not os.path.exists(out_mtz):
        raise Failure('expanding of R-free flags has failed -- {} does not exist'.format(out_mtz))

    if (delete_tmp_files is True):
        os.remove(tmp_mtz)

def run_refinement_pipeline(in_mtz, ref_pdb, out_dir, program='dimple'):
    """Run refinement of the input MTZ file against a reference PDB file"""

    logger.subheading('Running refinement pipelines')

    if program == 'dimple':
        # Define output files
        out_pdb = os.path.join(out_dir, 'final.pdb')
        out_mtz = os.path.join(out_dir, 'final.mtz')
        # Create command manager for dimple
        prog = Dispatcher('dimple')
        progcmd.extend_args([
            '--jelly', '5',
            in_mtz,
            ref_pdb,
            out_dir,
        ])
    else:
        raise ValueError("no stop that. you're doing it wrong.")

    run_program(prog)

    if not os.path.exists(out_pdb):
        raise Failure('running refinement with {} has failed -- {} does not exist'.format(program, out_pdb))
    if not os.path.exists(out_mtz):
        raise Failure('running refinement with {} has failed -- {} does not exist'.format(program, out_mtz))

    return out_pdb, out_mtz

############################################################################

def reprocess_mtz(in_mtz, params):
    """Prepare mtz by various standard sanitisations"""

    from giant.paths import foldername, filename

    # Determine naming for this dataset
    if params.output.labelling == 'foldername':
        label = foldername(in_mtz)
    elif params.output.labelling == 'filename':
        label = filename(in_mtz)
    else:
        label = None

    if label is not None:
        logger.heading('Processing dataset {}'.format(label))

    # Create output directory or identify it
    if params.output.output_directory is not None:
        out_dir = easy_directory(os.path.join(params.output.output_directory, label))
        # Copy input mtz to output directory and update variable
        new_in_mtz = os.path.join(out_dir, os.path.basename(in_mtz))
        shutil.copy(in_mtz, new_in_mtz)
        # Update variable to point to new file
        in_mtz = new_in_mtz
    else:
        # Simply take current directory -- likely not used
        out_dir = None

    # Save the original mtz as variable
    orig_mtz = in_mtz

    # List of intermediate files to be deleted at the end
    intermediate_files = []

    stage = 1

    if 'reindex_to_reference' in params.actions:

        out_mtz = splice_ext(path=orig_mtz, new='{}.reindexed'.format(stage))
        reindex_mtz_to_reference(
            in_mtz  = in_mtz,
            out_mtz = out_mtz,
            reference_mtz = params.input.reference_mtz,
            tolerance = params.options.reindex.tolerance,
        )
        intermediate_files.append(out_mtz)
        # Update new input mtz
        in_mtz = out_mtz
        stage += 1

    if 'fill_missing_reflections' in params.actions:

        out_mtz = splice_ext(path=orig_mtz, new='{}.filled'.format(stage))
        fill_missing_reflections(
            in_mtz  = in_mtz,
            out_mtz = out_mtz,
            fill_resolution_low  = params.options.missing_reflections.fill_low,
            fill_resolution_high = params.options.missing_reflections.fill_high,
            delete_tmp_files = (not params.output.keep_intermediate_files),
        )
        intermediate_files.append(out_mtz)
        # Update new input mtz
        in_mtz = out_mtz
        stage += 1

    if 'transfer_rfree_flags' in params.actions:

        out_mtz = splice_ext(path=orig_mtz, new='{}.with-free'.format(stage))
        transfer_rfree_flags(
            in_mtz  = in_mtz,
            out_mtz = out_mtz,
            reference_mtz = params.input.reference_mtz,
            input_free_r_flag  = params.options.rfree_flag.input,
            output_free_r_flag = params.options.rfree_flag.output,
            delete_tmp_files = (not params.output.keep_intermediate_files),
        )
        intermediate_files.append(out_mtz)
        # Update new input mtz
        in_mtz = out_mtz
        stage += 1

    # Copy to final output mtz
    out_mtz = splice_ext(path=orig_mtz, new='prepared')
    shutil.copy(in_mtz, out_mtz)
    # Update new input mtz
    in_mtz = out_mtz

    if 'refinement_pipeline' in params.actions:

        for pipeline in params.options.refinement.pipeline:

            if out_dir is None:
                ref_dir = os.path.splitext(in_mtz)[0] + '_' + pipeline
            else:
                ref_dir = os.path.join(out_dir, pipeline)

            o_pdb, o_mtz = run_refinement_pipeline(
                in_mtz  = in_mtz,
                ref_pdb = params.input.reference_pdb,
                out_dir = ref_dir,
                program = pipeline,
            )

    # Tidy up
    if (not params.output.keep_intermediate_files):
        for f in intermediate_files:
            os.remove(f)

    return out_mtz

############################################################################

SUCCESS_MESSAGE = 'processed successfully'

def wrapper_reprocess(args):
    params, mtz = args
    try:
        reprocess_mtz(in_mtz=mtz, params=params)
    except Exception as e:
        return (mtz, e)
    return (mtz, SUCCESS_MESSAGE)

############################################################################

def run(params):

    if len(params.input.mtz) == 0:
        raise IOError('no mtz files supplied (input.mtz=XXX)')

    if (params.output.output_directory is not None):
        params.output.output_directory = easy_directory(params.output.output_directory)
        if (params.output.labelling == 'none') or (not params.output.labelling):
            raise ValueError('Must provide labelling choice when using output directory (output.labelling=XXX)')

    if 'reindex_to_reference' in params.actions:
        if (params.input.reference_mtz is None):
            raise IOError('no reference mtz supplied (reference_mtz=XXX)')

    if 'fill_missing_reflections' in params.actions:
        if (params.options.missing_reflections.fill_high > params.options.missing_reflections.fill_low):
            raise ValueError('missing_reflections.fill_high must be less than missing_reflections.fill_low')

    if 'transfer_rfree_flags' in params.actions:
        if (params.input.reference_mtz is None):
            raise IOError('no reference mtz supplied (reference_mtz=XXX)')
        if (params.options.rfree_flag.input is None):
            raise ValueError('R-free flag must be specified (rfree_flag.input=XXX)')
        if (params.options.rfree_flag.output is None):
            params.options.rfree_flag.output = params.options.rfree_flag.input

    if 'refinement_pipeline' in params.actions:
        if params.input.reference_pdb is None:
            raise IOError('must provide reference_pdb to run refinement_pipeline (reference_pdb=XXX)')

    #################

    # Build argument list for multi-processing
    arg_list = [(params, m) for m in params.input.mtz]

    logger.heading('Re-processing {} mtz files'.format(len(arg_list)))
    logger.setLevel(lg.WARNING) # Turn off logging for multi-processing
    import libtbx.easy_mp
    returned = libtbx.easy_mp.pool_map(
        fixed_func   = wrapper_reprocess,
        args         = arg_list,
        processes    = params.settings.cpus,
        chunksize    = 1,
    )
    logger.setLevel(lg.INFO) # re-enable logging

    logger.subheading('Re-processing results')
    n_success = 0; errors = []
    for mtz, message in returned:
        if (message == success_message):
            logger('\t{}'.format(mtz))
            n_success += 1
        else:
            errors.append((mtz, message))

    if len(errors) > 0:
        logger.subheading('Errored datasets')
        for mtz, err in errors:
            logger('\nFile: {}\n\tError: {}'.format(mtz, type(message), str(message)))

    logger.bar()
    logger('Successfully processed:     {}'.format(n_success))
    logger('Errors during reprocessing: {}'.format(len(errors)))
    logger.bar()

############################################################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
        program             = PROGRAM,
        description         = DESCRIPTION,
    )

import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys, shutil
import traceback

import pathlib as pl

from giant.paths import (
    rel_symlink,
    )

from giant.processors import (
    ProcessorJoblib,
    )

from giant.dispatcher import (
    Dispatcher,
    )

from giant.mulch.labelling import (
    PathLabeller,
    )

from giant.mulch.dataset import (
    CrystallographicData,
    )

############################################################################

PROGRAM = 'giant.prepare_datasets'

DESCRIPTION = """
    Take input MTZs, reindex, transfer R-free flags from a reference, 
    fill missing reflections, and optionally run refinement pipelines.
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
    output_log = prepared_datasets.log
        .type = str
    output_directory = prepared_datasets
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
        pipeline = *dimple *pipedream
            .type = choice(multi=True)
    }
}
settings {
    cpus = 1
        .help = "number of cpus to use for processing"
        .type = int
    verbose = False
        .type = bool
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

        raise Exception(
            '{} exited with an error'.format(
                prog.program
                )
            )

def raise_missing(self, filepath, result):

    logger.subheading('Stdout')
    logger(result.stdout)
    
    logger.subheading('Stderr')
    logger(result.stderr)

    raise Exception(
        'Failed: {} does not exist'.format(
            str(filepath)
            )
        )

############################################################################


class ReindexMtzToReference(object):
    """
    Reindex the data in one mtz to a reference mtz
    """

    def __init__(self, 
        reference_mtz,
        out_suffix = ".reindex.mtz",
        tolerance = 5,
        ):

        self.reference_mtz = reference_mtz
        self.out_suffix = out_suffix
        self.tolerance = tolerance

        assert self.reference_mtz is not None

    def __call__(self, 
        in_mtz, 
        ):

        logger.subheading('Running reindexing')

        out_mtz = in_mtz.with_suffix(self.out_suffix)

        prog = Dispatcher('pointless')

        prog.extend_args([
            'hklin',  str(in_mtz),
            'hklref', str(self.reference_mtz),
            'hklout', str(out_mtz),
        ])

        prog.extend_stdin([
            'tolerance {}'.format(self.tolerance),
        ])

        run_program(prog)

        if not out_mtz.exists():
            raise_missing(
                filepath = out_mtz,
                result = prog.result,
                )

        return out_mtz


class FillMissingReflections(object): 

    def __init__(self, 
        fill_resolution_low = None, 
        fill_resolution_high = None,
        out_suffix = ".filled.mtz",
        delete_tmp_files = True,
        ):

        self.fill_resolution_low = fill_resolution_low
        self.fill_resolution_high = fill_resolution_high
        self.out_suffix = out_suffix
        self.delete_tmp_files = delete_tmp_files

    def __call__(self,
        in_mtz, 
        ):
        """Complete the set of miller indices in an MTZ file"""

        logger.subheading('Filling missing reflections')

        out_mtz = in_mtz.with_suffix(self.out_suffix)

        tmp_mtz_1 = out_mtz.with_suffix('.step1-truncate.mtz')
        tmp_mtz_2 = out_mtz.with_suffix('.step2-uniquify.mtz')
        tmp_mtz_3 = out_mtz.with_suffix('.step3-remerged.mtz')

        fill_low, fill_high = self.get_fill_limits(in_mtz)

        # Stage 1 - truncate dataset
        self.truncate(
            in_mtz = in_mtz, 
            out_mtz = tmp_mtz_1,
            resolution_low = fill_low, 
            resolution_high = fill_high,
            )

        # Stage 2 - Uniqueify the file
        self.uniqueify(
            in_mtz = tmp_mtz_1, 
            out_mtz = tmp_mtz_2,
            )

        # Stage 3 - remerge the two files
        self.remerge_with_dummy_column(
            in_mtz_1 = in_mtz,
            in_mtz_2 = tmp_mtz_2,
            out_mtz = tmp_mtz_3,
            )

        # Stage 4 - remove the dummy column
        self.remove_dummy_column(
            in_mtz = tmp_mtz_3,
            out_mtz = out_mtz,
            )

        self.cleanup(
            tmp_files = [
                tmp_mtz_1, 
                tmp_mtz_2, 
                tmp_mtz_3, 
                ],
            )

        return out_mtz

    def get_fill_limits(self, 
        in_mtz,
        ):

        fill_low = self.fill_resolution_low
        fill_high = self.fill_resolution_high

        m = CrystallographicData.from_file(
            str(in_mtz)
            )

        if fill_low is None: 
            fill_low = m.crystal.resolution_low
        
        if fill_high is None: 
            fill_high = m.crystal.resolution_high

        if fill_low is None: 
            raise Exception('No fill_low supplied or extracted from mtz')

        if fill_high is None: 
            raise Exception('No fill_high supplied or extracted from mtz')

        return (fill_low, fill_high)

    def truncate(self, 
        in_mtz, 
        out_mtz, 
        resolution_low, 
        resolution_high,
        ):

        assert resolution_low is not None
        assert resolution_high is not None

        prog = Dispatcher('cad')

        prog.extend_args([
            'hklin1', str(in_mtz),
            'hklout', str(out_mtz),
        ])

        prog.extend_stdin([
            'monitor BRIEF',
            'labin file_number 1 ALL',
            'resolution file 1 {low} {high}'.format(
                low = resolution_low, 
                high = resolution_high,
                ),
        ])

        run_program(prog)

        if not out_mtz.exists():
            raise_missing(
                filepath = out_mtz, 
                result = prog.result,
                )

    def uniqueify(self, 
        in_mtz, 
        out_mtz,
        ):

        prog = Dispatcher('uniqueify')

        prog.extend_args([
            '-p', '0.05',
            str(in_mtz),
            str(out_mtz),
        ])

        run_program(prog)

        if not out_mtz.exists():
            raise_missing(
                filepath = out_mtz, 
                result = prog.result,
                )

        # Annoying uncontrollable log file! (in the cwd)
        possible_log = pl.Path(out_mtz.name).with_suffix('.log')
        if possible_log.exists():
            os.remove(str(possible_log))

    def remerge_with_dummy_column(self, in_mtz_1, in_mtz_2, out_mtz):

        prog = Dispatcher('cad')

        prog.extend_args([
            'hklin1', str(in_mtz_1),
            'hklin2', str(in_mtz_2),
            'hklout', str(out_mtz),
        ])

        prog.extend_stdin([
            'monitor BRIEF',
            'labin file_number 1 ALL',
            'labin file_number 2 E1=FreeR_flag',
            'labout file_number 2 E1=dummy',
        ])

        run_program(prog)

        if not out_mtz.exists():
            raise_missing(
                filepath = out_mtz, 
                result = prog.result,
                )

    def remove_dummy_column(self, in_mtz, out_mtz):

        prog = Dispatcher('mtzutils')

        prog.extend_args([
            'hklin1', str(in_mtz),
            'hklout', str(out_mtz),
        ])

        prog.extend_stdin([
            'HEADER BRIEF',
            'EXCLUDE 1 dummy',
            'ONEFILE',
            'END',
        ])

        run_program(prog)

        if not out_mtz.exists():
            raise_missing(
                filepath = out_mtz, 
                result = prog.result,
                )

    def cleanup(self, tmp_files):

        if (self.delete_tmp_files is True):
            for f in tmp_files:
                os.remove(str(f))


class TransferRFreeFlags(object):

    def __init__(self, 
        reference_mtz, 
        input_free_r_flag, 
        output_free_r_flag = None,
        out_suffix = ".free.mtz",
        delete_tmp_files = True,
        ):
        
        if output_free_r_flag is None:
            output_free_r_flag = input_free_r_flag

        self.reference_mtz = reference_mtz
        self.input_free_r_flag = input_free_r_flag 
        self.output_free_r_flag = output_free_r_flag
        self.out_suffix = out_suffix
        self.delete_tmp_files = delete_tmp_files

    def __call__(self,
        in_mtz, 
        ):
        """Copy R-free flags from reference mtz"""

        logger.subheading('Transferring and completing R-free flags')

        out_mtz = in_mtz.with_suffix(self.out_suffix)

        tmp_mtz = out_mtz.with_suffix('.step1-transfer.mtz')

        # Stage 1 - transfer R-free from reference
        self.transfer_rfree_flags(
            in_mtz = in_mtz,
            reference_mtz = self.reference_mtz,
            out_mtz = tmp_mtz,
            input_free_r_flag = self.input_free_r_flag,
            output_free_r_flag = self.output_free_r_flag,
            )

        # Stage 2 - populate missing R-free values
        self.fill_missing_flags(
            in_mtz = tmp_mtz,
            out_mtz = out_mtz,
            free_r_flag = self.output_free_r_flag,
            )

        if (self.delete_tmp_files is True):
            os.remove(str(tmp_mtz))

        return out_mtz

    def transfer_rfree_flags(self, 
        in_mtz,
        reference_mtz,
        out_mtz,
        input_free_r_flag,
        output_free_r_flag,
        ):

        prog = Dispatcher('cad')

        prog.extend_args([
            'hklin1', str(in_mtz),
            'hklin2', str(reference_mtz),
            'hklout', str(out_mtz),
        ])

        prog.extend_stdin([
            'labin file_number 1 ALL',
            'labin file_number 2 E1={}'.format(
                input_free_r_flag
                ),
            'labout file_number 2 E1={}'.format(
                output_free_r_flag
                ),
            'END',
        ])

        run_program(prog)

        if not out_mtz.exists():
            raise_missing(
                filepath = out_mtz, 
                result = prog.result,
                )

    def fill_missing_flags(self, 
        in_mtz, 
        out_mtz, 
        free_r_flag,
        ):

        prog = Dispatcher('freerflag')

        prog.extend_args([
            'hklin',  str(in_mtz),
            'hklout', str(out_mtz),
        ])

        prog.extend_stdin([
            'COMPLETE FREE={}'.format(free_r_flag),
            'END',
        ])

        run_program(prog)

        if not out_mtz.exists():
            raise_missing(
                filepath = out_mtz, 
                result = prog.result,
                )


class RefinementPipelineResults(object):

    def __init__(self, pdb, mtz):

        self.pdb = pdb
        self.mtz = mtz


class RunPipeline(object):

    def __init__(self, **kwargs):

        self.setup(**kwargs)

        self.validate()


class RunDimplePipeline(RunPipeline):

    def setup(self,
        reference_mtz = None,
        n_jelly_body = 5,
        ):

        self.reference_mtz = reference_mtz
        self.n_jelly_body = n_jelly_body

    def validate(self):

        assert isinstance(self.n_jelly_body, int)

    def make_args(self, prog):

        if self.n_jelly_body is not None:
            prog.extend_args([
                '--jelly', str(self.n_jelly_body),
            ])

    def __call__(self, 
        in_pdb,
        in_mtz,
        out_dir,
        ):

        self.validate()

        prog = Dispatcher('dimple')

        self.make_args(prog)

        prog.extend_args([
            str(in_pdb),
            str(in_mtz),
            str(out_dir),
        ])

        run_program(prog)

        # Define output files
        out_pdb = (out_dir / 'final.pdb')
        out_mtz = (out_dir / 'final.mtz')

        if not out_pdb.exists():
            raise_missing(
                filepath = out_pdb, 
                result = prog.result,
                )

        if not out_mtz.exists():
            raise_missing(
                filepath = out_mtz, 
                result = prog.result,
                )

        return RefinementPipelineResults(
            pdb = out_pdb,
            mtz = out_mtz,
            )


class RunPipedreamPipeline(RunPipeline):

    def setup(self,
        protocol = "default",
        keep_waters = False,
        reference_mtz = None,
        restrain_to_input = False,
        n_threads = 1,
        ):

        self.protocol = protocol
        self.keep_waters = keep_waters
        self.reference_mtz = reference_mtz
        self.restrain_to_input = restrain_to_input
        self.n_threads = n_threads

    def validate(self):

        assert self.protocol in ['default','thorough','quick']

        assert isinstance(self.keep_waters, bool)

        if self.reference_mtz is not None:
            assert self.reference_mtz.exists()

        assert isinstance(self.restrain_to_input, bool)

    def make_args(self, prog):

        if self.protocol == 'thorough':
            prog.append_arg('-thorough')
        elif self.protocol == 'quick': 
            prog.append_arg('-quick')
        else: 
            pass

        if self.keep_waters is True:
            prog.append_arg('-keepwater')

        if self.reference_mtz is not None:
            prog.extend_args([
                '-hklref', str(self.reference_mtz),
                ])
        else:
            prog.append_arg('-nofreeref')

        if self.restrain_to_input is True:
            prog.append_arg('-target')

        if self.n_threads is not None: 
            prog.extend_args([
                '-nthreads', str(self.n_threads),
                ])

    def __call__(self, 
        in_pdb,
        in_mtz,
        out_dir,
        ):

        self.validate()

        prog = Dispatcher('pipedream')

        # pipedream \
        #     -imagedir <directory> \
        #     -d <output directory> \
        #     -xyzin input.pdb \
        #     -hklref input.mtz

        prog.extend_args([
            '-xyzin', str(in_pdb),
            '-hklin', str(in_mtz),
            '-d', str(out_dir),
        ])

        self.make_args(prog)

        run_program(prog)

        # Define output files
        out_pdb = (out_dir / 'refine' / 'refine.pdb')
        out_mtz = (out_dir / 'refine' / 'refine.mtz')

        if not out_pdb.exists():
            raise_missing(
                filepath = out_pdb, 
                result = prog.result,
                )

        if not out_mtz.exists():
            raise_missing(
                filepath = out_mtz, 
                result = prog.result,
                )

        return RefinementPipelineResults(
            pdb = out_pdb,
            mtz = out_mtz,
            )


class RunRefinementPipeline(object):
    """Run refinement of the input MTZ file against a reference PDB file"""

    def __init__(self, program='dimple', **program_kwargs):

        if program == 'dimple':
            run_pipeline = RunDimplePipeline(**program_kwargs)
        elif program == 'pipedream':
            run_pipeline = RunPipedreamPipeline(**program_kwargs)
        else: 
            raise NotImplemented()

        self.program = str(program)
        self.run_pipeline = run_pipeline

    def __call__(self, in_pdb, in_mtz, out_dir):

        logger.subheading(
            'Running {} refinement pipeline'.format(
                self.program
                )
            )

        pipeline_dir = (
            out_dir / self.program
            )

        result = self.run_pipeline(
            in_pdb = in_pdb,
            in_mtz = in_mtz,
            out_dir = pipeline_dir,
            )

        out_pdb = (
            out_dir / (self.program + '-final.pdb')
            )

        out_mtz = (
            out_pdb.with_suffix('.mtz')
            )

        rel_symlink(
            str(result.pdb), str(out_pdb),
            )

        rel_symlink(
            str(result.mtz), str(out_mtz),
            )

        return RefinementPipelineResults(
            pdb = out_pdb, 
            mtz = out_mtz,
            )


############################################################################


class PrepareMTZWrapper(object):

    def __init__(self, func):

        self.func = func

    def __call__(self, in_mtz, *args, **kwargs):

        try:

            self.func(in_mtz=in_mtz, *args, **kwargs)

        except Exception as e:
            return (in_mtz, e)

        return (in_mtz, None)


class PrepareMTZ(object):

    def __init__(self, params):

        self.output_directory = (
            params.output.output_directory
            if (params.output.output_directory is not None)
            else None
            )

        self.delete_tmp_files = (
            not params.output.keep_intermediate_files
            )

        self.labeller = PathLabeller(
            method = params.output.labelling,
            )
        
        self.reindex_to_reference = (
            ReindexMtzToReference(
                reference_mtz = params.input.reference_mtz,
                tolerance = params.options.reindex.tolerance,
                )
            if ('reindex_to_reference' in params.actions)
            else
            None
            )

        self.fill_missing_reflections = (
            FillMissingReflections(
                fill_resolution_low = params.options.missing_reflections.fill_low,
                fill_resolution_high = params.options.missing_reflections.fill_high,
                delete_tmp_files = (not params.output.keep_intermediate_files),
                )
            if ('fill_missing_reflections' in params.actions)
            else 
            None
            )

        self.transfer_rfree_flags = (
            TransferRFreeFlags(
                reference_mtz = params.input.reference_mtz,
                input_free_r_flag  = params.options.rfree_flag.input,
                output_free_r_flag = params.options.rfree_flag.output,
                delete_tmp_files = (not params.output.keep_intermediate_files),
                )
            if ('transfer_rfree_flags' in params.actions)
            else
            None
            )

        self.refinement_pipelines = (
            [
                RunRefinementPipeline(
                    program = program,
                    reference_mtz = params.input.reference_mtz,
                    )
                for program in params.options.refinement.pipeline
                ]
            if ('refinement_pipeline' in params.actions)
            else 
            None
            )

    def __call__(self, in_mtz, in_pdb=None):
        """Prepare mtz by various standard sanitisations"""

        in_mtz = pl.Path(in_mtz)

        label = self.labeller(str(in_mtz))

        if self.output_directory is not None:

            out_dir = (
                self.output_directory / label
                )

            if not out_dir.is_dir():
                out_dir.mkdir(parents=True)

            # Make a copy of the file in the new directory
            in_mtz = self.copy_to_directory(
                filepath = in_mtz,
                directory = out_dir,
                )

        else: 

            # All done in place
            out_dir = None 

        handler = self.start_logging(
            logpath = in_mtz.with_suffix('.prepared.log'),
            )

        # logging starts
        logger.heading(
            'Processing dataset: {}'.format(
                label
                )
            )

        try: 
            out_mtz = self.run(
                in_mtz = in_mtz, 
                in_pdb = in_pdb, 
                out_dir = out_dir,
                )
        except Exception as e: 
            logger(traceback.format_exc())
            raise

        self.stop_logging(handler)

        return out_mtz

    def run(self, in_mtz, in_pdb, out_dir):

        # Save the original mtz as variable
        orig_mtz = in_mtz

        # List of intermediate files to be deleted at the end
        intermediate_files = []

        if self.reindex_to_reference is not None:
            
            if str(orig_mtz) != str(in_mtz):
                intermediate_files.append(in_mtz)

            in_mtz = self.reindex_to_reference(
                in_mtz  = in_mtz,
            )

        if self.fill_missing_reflections is not None:
            
            if str(orig_mtz) != str(in_mtz):
                intermediate_files.append(in_mtz)

            in_mtz = self.fill_missing_reflections(
                in_mtz  = in_mtz,
                )

        if self.transfer_rfree_flags is not None: 

            if str(orig_mtz) != str(in_mtz):
                intermediate_files.append(in_mtz)

            in_mtz = self.transfer_rfree_flags(
                in_mtz  = in_mtz,
                )

        ###

        if str(orig_mtz) != str(in_mtz):
            intermediate_files.append(in_mtz)

        out_mtz = orig_mtz.with_suffix('.prepared.mtz')

        shutil.copy(
            str(in_mtz), 
            str(out_mtz),
            )

        ###

        if self.refinement_pipelines is not None: 

            if out_dir is not None:
                pipelines_dir = (
                    out_dir
                    )
            else: 
                pipelines_dir = (
                    out_mtz.parent / (out_mtz.stem + "_refinement_pipelines")
                    )

            if not pipelines_dir.is_dir():
                pipelines_dir.mkdir()

            for run_refinement_pipeline in self.refinement_pipelines:

                result = run_refinement_pipeline(
                    in_pdb = in_pdb,
                    in_mtz = out_mtz,
                    out_dir = pipelines_dir,
                    )

        # Tidy up
        if self.delete_tmp_files:
            for f in intermediate_files:
                os.remove(str(f))

        return out_mtz

    def start_logging(self, logpath):

        handler = lg.lg.FileHandler(
            filename = str(logpath),
            )

        logger.addHandler(handler)

        return handler

    def stop_logging(self, handler):

        logger.removeHandler(
            handler
            )

    def copy_to_directory(self, filepath, directory):

        out_path = (
            directory / filepath.name
            )

        shutil.copy(
            str(filepath), 
            str(out_path),
            )

        return out_path


############################################################################

def validate_params(params):

    if len(params.input.mtz) == 0:
        raise IOError(
            'no mtz files supplied (input.mtz=XXX)'
            )

    if (params.input.reference_mtz is not None):

        params.input.reference_mtz = pl.Path(params.input.reference_mtz)

        assert params.input.reference_mtz.exists()

    if (params.output.output_directory is not None):

        params.output.output_directory = pl.Path(params.output.output_directory)

        if not params.output.output_directory.is_dir():
            params.output.output_directory.mkdir(parents=True)

        if (params.output.labelling == 'none') or (not params.output.labelling):
            raise ValueError(
                'Must provide labelling choice when using output directory (output.labelling=XXX)'
                )

    if 'reindex_to_reference' in params.actions:

        if (params.input.reference_mtz is None):
            raise IOError(
                'no reference mtz supplied (reference_mtz=XXX)'
                )

    if 'fill_missing_reflections' in params.actions:

        if (
            (params.options.missing_reflections.fill_high is not None) and 
            (params.options.missing_reflections.fill_low is not None) and
            (params.options.missing_reflections.fill_high > params.options.missing_reflections.fill_low)
            ):
            raise ValueError(
                'missing_reflections.fill_high must be less than missing_reflections.fill_low'
                )

    if 'transfer_rfree_flags' in params.actions:

        if (params.input.reference_mtz is None):
            raise IOError(
                'no reference mtz supplied (reference_mtz=XXX)'
                )

        if (params.options.rfree_flag.input is None):
            raise ValueError(
                'R-free flag must be specified (rfree_flag.input=XXX)'
                )

        if (params.options.rfree_flag.output is None):
            params.options.rfree_flag.output = params.options.rfree_flag.input

    if 'refinement_pipeline' in params.actions:

        if params.input.reference_pdb is None:
            raise IOError(
                'must provide reference_pdb to run refinement_pipeline (reference_pdb=XXX)'
                )

        if params.input.reference_mtz is None:
            raise IOError(
                'must provide reference_mtz to run refinement_pipeline (reference_mtz=XXX)'
                )

def show_results(results):

    logger.subheading('Results')
    
    n_success = 0; errors = []

    logger('Processing messages')

    for mtz, message in results:

        if (message is None):

            logger('\t{} - {}'.format(mtz, 'processed successfully'))
            n_success += 1

        else:

            logger('\t{} - {}'.format(mtz, 'errored'))
            errors.append((mtz, message))

    if len(errors) > 0:

        logger.subheading('Errored datasets')

        logger('> See log files in output directory\n')
        
        for mtz, err in errors:

            logger(
                '\nFile: {}\n\tError: {}'.format(
                    mtz, str(message),
                    )
                )

    logger.bar()
    logger('Successfully processed:     {}'.format(n_success))
    logger('Errors during reprocessing: {}'.format(len(errors)))
    logger.bar()

def run(params):

    logpath = pl.Path(
        params.output.output_log
        )

    logger = lg.setup_logging(
        name = __name__,
        log_file = str(logpath),
    )

    validate_params(params)
    
    #################

    processor = ProcessorJoblib(
        n_cpus = params.settings.cpus,
        )

    prepare_mtz = PrepareMTZ(params)

    # Build argument list for multi-processing
    arg_list = [
        processor.make_wrapper(
            func = PrepareMTZWrapper(prepare_mtz), 
            in_mtz = m,
            in_pdb = params.input.reference_pdb,
            ) 
        for m in params.input.mtz
        ]
    
    logger.heading(
        'Re-processing {} mtz files'.format(
            len(arg_list)
            )
        )

    prev_propagate_state = logger.propagate

    if not params.settings.verbose:
        logger.propagate = False

    try: 
        results = processor(arg_list)
    except: 
        logger.propagate = prev_propagate_state
        raise
    finally:
        logger.propagate = prev_propagate_state

    show_results(
        results, 
        )


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

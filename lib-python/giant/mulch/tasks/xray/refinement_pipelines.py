import giant.logs as lg
logger = lg.getLogger(__name__)

import traceback

import pathlib as pl

from giant.mulch.tasks.utils import (
    run_program,
    raise_missing,
    )

from giant.mulch.tasks.io import (
    TaskReturn,
    TaskReturnStatus,
    ModelDataInputOutput,
    )

from giant.mulch.labelling import (
    PathLabeller,
    )

from giant.dispatcher import (
    Dispatcher,
    )

from giant.paths import (
    rel_symlink,
    )


class _Pipeline(object):

    def __init__(self, **kwargs):

        # Function to store args for the specific pipeline
        self.setup(**kwargs)

        # Validate status with provided args
        self.validate()


class DimplePipeline(_Pipeline):

    def setup(self,
        reference_mtz = None,
        n_jelly_body = None,
        ):

        self.reference_mtz = reference_mtz
        self.n_jelly_body = n_jelly_body

    def validate(self):

        if self.n_jelly_body is not None: 
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

        return ModelDataInputOutput(
            input_pdb = in_pdb,
            input_mtz = in_mtz,
            output_pdb = out_pdb,
            output_mtz = out_mtz,
            )


class PipedreamPipeline(_Pipeline):

    def setup(self,
        protocol = None,
        reference_mtz = None,
        keep_waters = None,
        restrain_to_input = None,
        mr_step = None,
        n_threads = None,
        ):

        self.protocol = protocol
        self.reference_mtz = reference_mtz
        self.keep_waters = keep_waters
        self.restrain_to_input = restrain_to_input
        self.mr_step = mr_step
        self.n_threads = n_threads

    def validate(self):

        if self.protocol is not None: 
            assert self.protocol in ['default','thorough','quick']

        if self.reference_mtz is not None:
            assert self.reference_mtz.exists()

        if self.keep_waters is not None: 
            assert isinstance(self.keep_waters, bool)

        if self.restrain_to_input is not None: 
            if self.restrain_to_input is not True: 
                assert isinstance(self.restrain_to_input, str)

        if self.mr_step is not None:
            assert isinstance(self.mr_step, bool)

        if self.n_threads is not None: 
            assert isinstance(self.n_threads, int)

    def make_args(self, prog):

        if self.protocol == 'thorough':
            prog.append_arg('-thorough')
        elif self.protocol == 'quick': 
            prog.append_arg('-quick')
        else: 
            pass

        if self.reference_mtz is False:
            prog.append_arg(
                '-nofreeref'
                )
        elif self.reference_mtz is not None:
            prog.extend_args([
                '-hklref', str(self.reference_mtz),
                ])
        else: 
            raise InputError('reference_mtz must be provided or set to False')

        if self.keep_waters is True:
            prog.append_arg('-keepwater')
        else:
            pass

        if str(self.restrain_to_input) == 'True':
            prog.append_arg(
                '-target'
                )
        elif str(self.restrain_to_input) == 'False':
            pass
        elif self.restrain_to_input is not None: 
            prog.extend_args([
                '-target',
                str(self.restrain_to_input),
                ])
        else: 
            pass

        if self.mr_step is False:
            prog.append_arg(
                '-nolmr'
                )
        else:
            pass

        if self.n_threads is not None: 
            prog.extend_args([
                '-nthreads', str(self.n_threads),
                ])
        else: 
            pass

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

        return ModelDataInputOutput(
            input_pdb = in_pdb,
            input_mtz = in_mtz,
            output_pdb = out_pdb,
            output_mtz = out_mtz,
            )


class RefinementPipelineTask(object):
    """Run refinement of the input MTZ file against a reference PDB file"""

    def __init__(self, 
        program = 'dimple', 
        program_options = None,
        output_directory = None,
        path_labeller = None,
        ):

        if program == 'dimple':
            run_pipeline = (
                DimplePipeline(
                    n_jelly_body = program_options.n_jelly_body,
                    )
                if 
                (program_options is not None)
                else 
                DimplePipeline()
                )
        elif program == 'pipedream':
            run_pipeline = (
                PipedreamPipeline(
                    protocol = program_options.protocol,
                    reference_mtz = program_options.reference_mtz,
                    keep_waters = program_options.keep_waters,
                    restrain_to_input = program_options.restrain_to_input,
                    mr_step = program_options.mr_step,
                    n_threads = program_options.n_threads,
                    )
                if 
                (program_options is not None)
                else
                PipedreamPipeline()
                )
        else: 
            raise NotImplemented()

        self.program = (
            str(program)
            )

        self.run_pipeline = (
            run_pipeline
            )

        self.output_directory = (
            output_directory
            )

        self.path_labeller = (
            path_labeller
            if 
            (path_labeller is not None)
            else 
            PathLabeller()
            )

    def __call__(self, in_pdb, in_mtz):

        in_pdb = pl.Path(in_pdb)
        in_mtz = pl.Path(in_mtz)

        label = self.path_labeller(str(in_mtz))

        if self.output_directory is not None:

            out_dir = (
                self.output_directory / label 
                )

        else: 

            out_dir = (
                in_mtz.parent / (in_mtz.stem + "_pipelines")
                )

        if not out_dir.is_dir():
            out_dir.mkdir(parents=True)

        file_handler, warning_handler = self.start_logging(
            logpath = (out_dir / (self.program + '.log')),
            )

        logger.subheading(
            'Running {} pipeline for {}'.format(
                self.program,
                label,
                )
            )

        try: 

            out_pdb, out_mtz = self.run(
                in_mtz = in_mtz,
                in_pdb = in_pdb,
                out_dir = out_dir,
                )

        except Exception as e: 

            # Ensure usable information is logged in the local logfile
            logger.warning(traceback.format_exc())

            self.stop_logging(
                handlers = [file_handler, warning_handler],
                )

            return TaskReturn(
                output = ModelDataInputOutput(
                    input_pdb = in_pdb,
                    input_mtz = in_mtz,
                    ),
                status = TaskReturnStatus(
                    success = False,
                    errors = [e],
                    warnings = list(warning_handler.list()),
                    ),
                )

        self.stop_logging(
            handlers = [file_handler, warning_handler],
            )

        return TaskReturn(
            output = ModelDataInputOutput(
                input_pdb = in_pdb,
                input_mtz = in_mtz,
                output_pdb = out_pdb,
                output_mtz = out_mtz,
                ),
            status = TaskReturnStatus(
                success = True,
                errors = None,
                warnings = list(warning_handler.list()),
                ),
            )


    def run(self, in_mtz, in_pdb, out_dir):

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
            str(result.output_pdb), str(out_pdb),
            )

        rel_symlink(
            str(result.output_mtz), str(out_mtz),
            )

        return out_pdb, out_mtz

    def start_logging(self, logpath):

        logger = lg.getLogger('giant')

        f_handler = lg.lg.FileHandler(
            filename = str(logpath),
            )

        logger.addHandler(f_handler)

        w_handler = lg.add_warning_handler(
            logger = logger,
            )

        return f_handler, w_handler

    def stop_logging(self, handlers):

        logger = lg.getLogger('giant')

        for h in handlers:
            logger.removeHandler(h)



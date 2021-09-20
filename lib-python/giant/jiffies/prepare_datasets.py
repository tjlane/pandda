import giant.logs as lg
logger = lg.getLogger(__name__) 

import sys

import pathlib as pl

from giant.phil import (
    log_running_parameters,
    )

from giant.processors import (
    ProcessorJoblib,
    )

from giant.mulch.labelling import (
    PathLabeller,
    )

from giant.mulch.tasks.xray.mtz_utils import (
    PrepareMTZTask,
    RemoveColumns,
    ReindexMtzToReference,
    FillMissingReflections,
    TransferRFreeFlags,
    )

from giant.mulch.tasks.xray.refinement_pipelines import (
    RefinementPipelineTask,
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
    labelling = none basename *filename foldername
        .type = choice
    keep_intermediate_files = False
        .type = bool
}
actions = remove_columns \
          *reindex_to_reference \
          *fill_missing_reflections \
          *transfer_rfree_flags \
          *refinement_pipeline
    .type = choice(multi=True)
options {
    remove_columns {
        columns = None
            .type = str
            .multiple = True
    }
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
        fill_high = None
            .type = float
        fill_low = 9999.0
            .type = float
    }
    refinement {
        pipelines = *dimple *pipedream
            .type = choice(multi=True)
        dimple {
            n_jelly_body = None
                .type = int
        }
        pipedream {
            protocol = quick *default thorough
                .type = choice(multi=False)
            reference_mtz = None
                .type = path
            keep_waters = None
                .type = bool
            restrain_to_input = None
                .type = str
            mr_step = None
                .type = bool
            n_threads = 1
                .type = int
        }
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

def show_results(results):

    logger.subheading('Results')
    
    n_success = 0; errored = []

    logger('Processing messages')

    for r in results:

        if r.status.success is True: 

            logger(
                '\t{} - {}'.format(
                    r.output.input_mtz, 
                    'processed successfully'
                    )
                )
            n_success += 1

        else:

            logger(
                '\t{} - {}'.format(
                    r.output.input_mtz, 
                    'errored'
                    )
                )
            errored.append(r)

    if len(errored) > 0:

        logger.subheading('Errored datasets')

        logger('> See log files in output directory\n')
        
        for r in errored:

            logger(
                '\nFile: {}\n\tError(s): \n\t\t{}'.format(
                    r.output.input_mtz, 
                    '\n\t\t'.join([
                        str(e) 
                        for e in r.status.errors
                        ]),
                    )
                )

    logger.bar()
    logger('Successfully processed:     {}'.format(n_success))
    logger('Errors during processing: {}'.format(len(errored)))
    logger.bar()


class RunRefinementPipelineJobs(object):

    def __init__(self, 
        pipeline_functions,
        processor,
        verbose = False,
        ):

        self.pipeline_functions = pipeline_functions
        self.processor = processor
        self.verbose = verbose

    def __call__(self, 
        mtz_list,
        reference_pdb,
        ):

        arg_list = []

        for pf in self.pipeline_functions:

            arg_list.extend([
                self.processor.make_wrapper(
                    func = pf,
                    in_mtz = m,
                    in_pdb = reference_pdb,
                    )
                for m in mtz_list
                ])

        # Control module logging
        glg = lg.getLogger('giant')

        prev_propagate_state = glg.propagate

        if not self.verbose:
            glg.propagate = False

        try: 

            results = self.processor(arg_list)

        except Exception as e: 

            glg.propagate = prev_propagate_state
            raise

        finally:

            glg.propagate = prev_propagate_state

        show_results(results)

        return results

    
class RunPrepareMTZJobs(object):

    def __init__(self, 
        prepare_mtz,
        processor, 
        verbose = False,
        ):

        self.prepare_mtz = prepare_mtz
        self.processor = processor
        self.verbose = verbose

    def __call__(self,
        mtz_list,
        ):

        arg_list = [
            self.processor.make_wrapper(
                func = self.prepare_mtz, 
                in_mtz = m,
                ) 
            for m in mtz_list
            ]

        # Control module logging
        glg = lg.getLogger('giant')

        prev_propagate_state = glg.propagate

        if not self.verbose:
            glg.propagate = False

        try: 

            results = self.processor(arg_list)

        except Exception as e: 

            glg.propagate = prev_propagate_state
            raise

        finally:

            glg.propagate = prev_propagate_state

        show_results(results)

        return results


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

        # hack alert
        if params.options.refinement.pipedream.reference_mtz is None:
            params.options.refinement.pipedream.reference_mtz = (
                params.input.reference_mtz
                )

def run(params):

    logpath = pl.Path(
        params.output.output_log
        )

    logger = lg.setup_logging(
        name = __name__,
        log_file = str(logpath),
    )

    log_running_parameters(
        params = params,
        master_phil = master_phil,
        logger = logger,
    )

    validate_params(params)
    
    pipeline_options = {
        'dimple' : params.options.refinement.dimple,
        'pipedream' : params.options.refinement.pipedream,
    }

    #################

    # setup

    path_labeller = PathLabeller(
        method = params.output.labelling,
        )

    processor = ProcessorJoblib(
        n_cpus = params.settings.cpus,
        )

    prepare_mtzs = RunPrepareMTZJobs(
        prepare_mtz = PrepareMTZTask(
            remove_columns = (
                RemoveColumns(
                    columns_to_remove = params.options.remove_columns.columns,
                    )
                if ('remove_columns' in params.actions)
                else 
                None
                ),
            reindex_to_reference = (
                ReindexMtzToReference(
                    reference_mtz = params.input.reference_mtz,
                    tolerance = params.options.reindex.tolerance,
                    )
                if ('reindex_to_reference' in params.actions)
                else
                None
                ),
            fill_missing_reflections = (
                FillMissingReflections(
                    fill_resolution_low = params.options.missing_reflections.fill_low,
                    fill_resolution_high = params.options.missing_reflections.fill_high,
                    delete_tmp_files = (not params.output.keep_intermediate_files),
                    )
                if ('fill_missing_reflections' in params.actions)
                else 
                None
                ),
            transfer_rfree_flags = (
                TransferRFreeFlags(
                    reference_mtz = params.input.reference_mtz,
                    input_free_r_flag  = params.options.rfree_flag.input,
                    output_free_r_flag = params.options.rfree_flag.output,
                    delete_tmp_files = (not params.output.keep_intermediate_files),
                    )
                if ('transfer_rfree_flags' in params.actions)
                else
                None
                ),
            output_directory = (
                params.output.output_directory
                ),
            path_labeller = (
                path_labeller
                ),
            delete_tmp_files = (
                not params.output.keep_intermediate_files
                ),
            ),
        processor = processor, 
        verbose = params.settings.verbose,
        )

    run_refinement_pipelines = (
        RunRefinementPipelineJobs(
            pipeline_functions = [
                RefinementPipelineTask(
                    program = p,
                    program_options = pipeline_options[p],
                    output_directory = (
                        None #params.output.output_directory
                        ),
                    path_labeller = (
                        path_labeller
                        ),
                    )
                for p in params.options.refinement.pipelines
                ],
            processor = processor,
            verbose = params.settings.verbose,
            )
        if 
        ('refinement_pipeline' in params.actions)
        else 
        None
        )

    # run

    logger.heading('Preparing MTZs')

    prepared_mtzs = prepare_mtzs(
        mtz_list = params.input.mtz,
        )

    if run_refinement_pipelines is not None: 

        logger.heading('Running refinement pipelines')

        run_refinement_pipelines(
            reference_pdb = params.input.reference_pdb,
            mtz_list = [
                p.output.output_mtz 
                for p in prepared_mtzs
                if 
                (p.status.success is True)
                ],
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

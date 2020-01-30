from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure

import pandemic.adp.preprocess.load_models
import pandemic.adp.preprocess.process_models
import pandemic.adp.preprocess.select_datasets
import pandemic.adp.preprocess.extract_uijs

from bamboo.common.command import check_programs_are_available

def validate_parameters(params, log=None):

    log.heading('Processing input parameters')

    if len(params.input.pdb) == 0:
        raise Sorry('No structures have been provided for analysis (input.pdb=[...])')

    if params.settings.debug:
        log('DEBUG is turned on -- setting verbose=True')
        params.settings.verbose = True

    # Output images
    if params.output.images.all:
        log('params.output.images.all = {}'.format(params.output.images.all))
        # ----------------
        new = 'chain'
        log('Setting params.output.images.pymol to {}'.format(new))
        params.output.images.pymol = new
        # ----------------
        new = True
        log('Setting params.output.images.distributions to {}'.format(new))
        params.output.images.distributions = new
        # ----------------
        log.bar()

    # Special settings if only one model is loaded
    if (len(params.input.pdb) == 1):
        log.bar(True, False)
        log('One file provided for analysis -- updating settings')
        log('Setting fitting.optimisation.dataset_weights = one')
        params.optimisation.weights.dataset_weights = 'one'
        log.bar()

    try:
        # Initalise variable to blank string (overwritten at each check)
        message = ''
        # Standard string for missing programs
        missing_program_info_string = """
        The program is not currently available in any of the PATH directories.
        [ It must be available as an executable script/link -- not as an alias.  ]
        """

        if params.analysis.refine_output_structures is True:
            if params.refinement.program == 'phenix':
                message = """
                phenix.refine is required when analysis.refine_output_structures=True and refinement.program=phenix.
                {missing_program_info_string}
                """.format(missing_program_info_string=missing_program_info_string)
                check_programs_are_available(['phenix.refine'])
            elif params.refinement.program == 'refmac':
                message = """
                refmac5 is required when analysis.refine_output_structures=True and refinement.program=refmac.
                {missing_program_info_string}
                """.format(missing_program_info_string=missing_program_info_string)
                check_programs_are_available(['refmac5'])
            else:
                raise Sorry('Must select refinement.program when analysis.refine_output_structures=True')

        if params.analysis.calculate_r_factors is True:
            message = """
            phenix.table_one is required when analysis.calculate_r_factors=True.
            {missing_program_info_string}
            To turn off this function add analysis.calculate_r_factors=False to command line options.
            """.format(missing_program_info_string=missing_program_info_string)
            check_programs_are_available(['phenix.table_one'])

        #if params.analysis.calculate_electron_density_metrics:
        #    message = 'edstats (script name "edstats.pl") is required when analysis.calculate_electron_density_metrics is True'
        #    check_programs_are_available(['edstats.pl'])

        if params.output.images.pymol is not None:
            message = """
            pymol is required when output.images.pymol is not set to "none".
            {missing_program_info_string}
            To turn off pymol add output.images.pymol=none to the command line options.
            """.format(missing_program_info_string=missing_program_info_string)
            check_programs_are_available(['pymol'])
    except Exception as e:
        log(message)
        raise

    # Simple checks
    if params.optimisation.min_macro_cycles > params.optimisation.max_macro_cycles:
        raise Sorry("max_macro_cycles must be greater than or equal to min_macro_cycles! ({} < {})".format(
            params.optimisation.min_macro_cycles,
            params.optimisation.max_macro_cycles,
            ))


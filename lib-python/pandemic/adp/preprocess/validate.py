import textwrap

from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure

from bamboo.common.command import check_programs_are_available

def phil_path(params, attribute):
    return params.__phil_path__(attribute)

def phil_value(params, attribute):
    return '{!s}={!s}'.format(*params.__phil_path_and_value__(attribute))

def phil_new_value(params, attribute, value):
    return '{!s}={!s}'.format(phil_path(params, attribute), value)

def update_parameter(params, attribute, value):
    update_string = 'setting {!s}'.format(phil_new_value(params, attribute, value))
    setattr(params, attribute, value)
    return update_string

def validate_parameters(params, log=None):

    log.heading('Processing input parameters')

    validate_input_phil(params, log)
    validate_settings_phil(params, log)
    validate_optimisation_phil(params, log)
    validate_output_phil(params, log)
    validate_available_programs(params, log)

def validate_input_phil(params, log=None):

    if len(params.input.pdb) == 0:
        raise Sorry('No structures have been provided for analysis (input.pdb=[...])')

    # Special settings if only one model is loaded
    if (len(params.input.pdb) == 1):
        update_string = update_parameter(params.optimisation.weights, 'dataset_weights', 'one')
        log("\nOne file provided for analysis: {update_string}.".format(update_string=update_string))

def validate_settings_phil(params, log=None):

    if params.settings.debug:
        update_string = update_parameter(params.settings, 'verbose', True)
        log('\nDEBUG is turned on: {update_string}.'.format(update_string=update_string))

def validate_optimisation_phil(params, log=None):

    # Simple checks
    if params.optimisation.min_macro_cycles > params.optimisation.max_macro_cycles:
        current_string = textwrap.dedent("""
            min_macro_cycles is greater than max_macro_cycles ({min_macro_cycles} > {max_macro_cycles}).
            """).strip().format(
                min_macro_cycles=params.optimisation.min_macro_cycles, 
                max_macro_cycles=params.optimisation.max_macro_cycles,
                )
        update_string = update_parameter(params.optimisation, 'min_macro_cycles', params.optimisation.max_macro_cycles)
        log('\n'+textwrap.dedent("""
            {current_string}.
            {update_string}.
            """).strip().format(
                current_string=current_string,
                update_string=update_string,
                )
            )

    if params.optimisation.min_macro_cycles < 1: 
        raise Sorry('min_macro_cycles must be greater than 0 ({current_params})'.format(phil_value(params.optimisation, 'min_macro_cycles')))

    if params.optimisation.max_macro_cycles < 1: 
        raise Sorry('max_macro_cycles must be greater than 0 ({current_params})'.format(phil_value(params.optimisation, 'max_macro_cycles')))

def validate_output_phil(params, log=None):

    # Output images
    if params.output.images.all:
        
        log('\n'+phil_value(params.output.images, 'all'))

        update_string = update_parameter(params.output.images, 'pymol', 'chain')
        log(update_string.capitalize())

def validate_available_programs(params, log=None):

    try:
        # Initalise variable to blank string (overwritten at each check)
        message = ''
        # Standard string for missing programs
        missing_program_info_string = textwrap.dedent("""
        The program is not currently available in any of the PATH directories.
        [ It must be available as an executable script/link -- not as an alias.  ]
        """)

        if params.analysis.refine_output_structures is True:

            # Format string showing the relevant current parameters
            current_params_string = ' and '.format([
                phil_value(params.analysis, 'refine_output_structures'),
                phil_value(params.refinement, 'program'),
                ])

            if params.refinement.program == 'phenix':
                message = """
                phenix.refine is required when {current_params}.
                {missing_program_info_string}
                """.format(
                    current_params=current_params_string,
                    missing_program_info_string=missing_program_info_string,
                    )
                check_programs_are_available(['phenix.refine'])
            elif params.refinement.program == 'refmac':
                message = """
                refmac5 is required when {current_params}.
                {missing_program_info_string}
                """.format(
                    current_params=current_params_string,
                    missing_program_info_string=missing_program_info_string,
                    )
                check_programs_are_available(['refmac5'])
            else:
                raise Sorry("""
                    Must select refinement.program when {requirement}.
                    """.format(
                        requirement=phil_value(params.analysis, 'refine_output_structures'),
                        )
                    )

        if params.analysis.calculate_r_factors is True:
            message = textwrap.dedent("""
            phenix.table_one is required when {current_params}.
            {missing_program_info_string}
            To turn off this function add {fix_params} to command line options.
            """).format(
                current_params=' and '.format([
                    phil_value(params.analysis, 'calculate_r_factors'),
                    ]),
                missing_program_info_string=missing_program_info_string,
                fix_params=phil_new_value(params.analysis, 'calculate_r_factors', False),
                )
            check_programs_are_available(['phenix.table_one'])

        #if params.analysis.calculate_electron_density_metrics:
        #    message = 'edstats (script name "edstats.pl") is required when analysis.calculate_electron_density_metrics is True'
        #    check_programs_are_available(['edstats.pl'])

        if params.output.images.pymol is not None:
            message = textwrap.dedent("""
            pymol is required when {current_params}.
            {missing_program_info_string}
            To turn off pymol add {fix_params} to the command line options.
            """).format(
                current_params=' and '.format([
                    phil_value(params.output.images, 'pymol'),
                    ]),
                missing_program_info_string=missing_program_info_string,
                fix_params=phil_new_value(params.output.images, 'pymol', 'none'),
                )
            check_programs_are_available(['pymol'])
    except Exception as e:
        log(message)
        raise

    return True

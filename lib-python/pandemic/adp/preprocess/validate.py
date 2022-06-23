import giant.logs as lg
logger = lg.getLogger(__name__)

import textwrap
import multiprocessing

from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure

def elastic_net_weight_summary(params):

    summary_lines = []
    summary_lines.append('Overall weight scale: {}'.format(params.weight_scale))
    for k, w in params.weights.__dict__.items():
        if k.startswith('_'):
            continue
        summary_lines.append('> {}: {}'.format(k, w))
    summary_lines.append(
        'sum_of_amplitudes + sum_of_squared_amplitudes = {}'.format(
            params.weights.sum_of_amplitudes + params.weights.sum_of_squared_amplitudes
        )
    )

    return '\n'.join(summary_lines)

####

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

####

def validate_parameters(params):

    logger.heading('Processing input parameters')

    validate_input_phil(params)
    validate_settings_phil(params)
    validate_optimisation_phil(params)
    validate_output_phil(params)
    validate_available_programs(params)

def validate_input_phil(params):

    if len(params.input.pdb) == 0:
        raise Sorry('No structures have been provided for analysis (input.pdb=[...])')

    # Special settings if only one model is loaded
    if (len(params.input.pdb) == 1):
        update_string = update_parameter(params.optimisation.weights, 'dataset_weights', 'one')
        logger("\nOne file provided for analysis: {update_string}.".format(update_string=update_string))

def validate_settings_phil(params):

    if params.settings.debug:
        update_string = update_parameter(params.settings, 'verbose', True)
        logger('\nDEBUG is turned on: {update_string}.'.format(update_string=update_string))

    if params.settings.cpus is None: 
        update_string = update_parameter(params.settings, 'cpus', multiprocessing.cpu_count())
        logger('\nDEBUG is turned on: {update_string}.'.format(update_string=update_string))


def validate_optimisation_phil(params):

    # Simple checks
    if params.optimisation.min_macro_cycles > params.optimisation.max_macro_cycles:
        current_string = textwrap.dedent("""
            min_macro_cycles is greater than max_macro_cycles ({min_macro_cycles} > {max_macro_cycles}).
            """).strip().format(
                min_macro_cycles=params.optimisation.min_macro_cycles,
                max_macro_cycles=params.optimisation.max_macro_cycles,
                )
        update_string = update_parameter(params.optimisation, 'min_macro_cycles', params.optimisation.max_macro_cycles)
        logger('\n'+textwrap.dedent("""
            {current_string}.
            {update_string}.
            """).strip().format(
                current_string=current_string,
                update_string=update_string,
                )
            )

    # Cycle checks
    if params.optimisation.min_macro_cycles < 1:
        raise Sorry(
            '{param_name} must be greater than 0 ({current_params})'.format(
                param_name = 'min_macro_cycles',
                current_params = phil_value(params.optimisation, 'min_macro_cycles'),
            )
        )
    if params.optimisation.max_macro_cycles < 1:
        raise Sorry(
            '{param_name} must be greater than 0 ({current_params})'.format(
                param_name = 'max_macro_cycles',
                current_params = phil_value(params.optimisation, 'max_macro_cycles'),
            )
        )

    ##########################
    # Elastic net parameters #
    ##########################

    e_net_params = params.optimisation.elastic_net

    # Check weights are positive or zero
    if e_net_params.weights.sum_of_amplitudes < 0.0:
        raise Sorry(
            '{param_name} must be greater than 0 ({current_params})'.format(
                param_name = 'sum_of_amplitudes',
                current_params = phil_value(e_net_params.weights, 'sum_of_amplitudes')
            )
        )
    if e_net_params.weights.sum_of_squared_amplitudes < 0.0:
        raise Sorry(
            '{param_name} must be greater than 0 ({current_params})'.format(
                param_name = 'sum_of_squared_amplitudes',
                current_params = phil_value(e_net_params.weights, 'sum_of_squared_amplitudes')
            )
        )
    if e_net_params.weights.sum_of_amplitudes_squared < 0.0:
        raise Sorry(
            '{param_name} must be greater than 0 ({current_params})'.format(
                param_name = 'sum_of_amplitudes_squared',
                current_params = phil_value(e_net_params.weights, 'sum_of_amplitudes_squared')
            )
        )

    # Ensure that the sum of sum_of_amplitudes and sum_of_square_amplitudes is 1
    wgt_sum = (e_net_params.weights.sum_of_amplitudes + e_net_params.weights.sum_of_squared_amplitudes)
    if (wgt_sum > 0.0) and (wgt_sum != 1.0):

        # Scale to apply to the weights
        wgt_scale = 1.0 / wgt_sum

        # Log current values
        logger('\nElastic net weight normalisation')
        logger('\nInput elastic net weights')
        logger(elastic_net_weight_summary(e_net_params))

        logger('\nScaling input weights sum_of_amplitudes and sum_of_squared_amplitudes to sum to 1.')
        logger('Scaling weights internally by {}'.format(round(wgt_scale,6)))

        # Scale weights
        from pandemic.adp.weights import scale_weights
        e_net_params.weights = scale_weights(weights=e_net_params.weights, scale=wgt_scale)
        e_net_params.weight_scale = e_net_params.weight_scale / wgt_scale

        # Log new values
        logger('\nUpdated weights')
        logger(elastic_net_weight_summary(e_net_params))

def validate_output_phil(params):

    # Output images
    if params.output.images.all:

        logger('\n'+phil_value(params.output.images, 'all'))

        update_string = update_parameter(params.output.images, 'pymol', 'chain')
        logger(update_string.capitalize())

def validate_available_programs(params):

    from giant.paths import is_available

    # Initalise variable to blank string (overwritten at each check)
    message = ''
    # Standard string for missing programs
    missing_program_info_string = textwrap.dedent("""
    The program is not currently available in any of the PATH directories.
    [ It must be available as an executable script/link -- not as an alias.  ]
    """)

    if params.analysis.refine_output_structures is True:

        # Format string showing the relevant current parameters
        current_params_string = ' and '.join([
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
            if not is_available('phenix.refine'):
                raise Sorry(message)
        elif params.refinement.program == 'refmac':
            message = """
            refmac5 is required when {current_params}.
            {missing_program_info_string}
            """.format(
                current_params=current_params_string,
                missing_program_info_string=missing_program_info_string,
                )
            if not is_available('refmac5'):
                raise Sorry(message)
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
            current_params=' and '.join([
                phil_value(params.analysis, 'calculate_r_factors'),
                ]),
            missing_program_info_string=missing_program_info_string,
            fix_params=phil_new_value(params.analysis, 'calculate_r_factors', False),
            )
        if not is_available('phenix.table_one'):
            raise Sorry(message)

    #if params.analysis.calculate_electron_density_metrics:
    #    message = 'edstats (script name "edstats.pl") is required when analysis.calculate_electron_density_metrics is True'
    #    check_programs_are_available(['edstats.pl'])

    if params.output.images.pymol is not None:
        message = textwrap.dedent("""
        pymol is required when {current_params}.
        {missing_program_info_string}
        To turn off pymol add {fix_params} to the command line options.
        """).format(
            current_params=' and '.join([
                phil_value(params.output.images, 'pymol'),
                ]),
            missing_program_info_string=missing_program_info_string,
            fix_params=phil_new_value(params.output.images, 'pymol', 'none'),
            )
        if not is_available('pymol'):
            raise Sorry(message)

    return True

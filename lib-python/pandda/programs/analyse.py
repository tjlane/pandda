import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys, time, copy

import pathlib as pl

from pandda import (
    ModuleBanner,
    phil,
    config,
    options,
    output,
    )

from pandda.utils import (
    DataCollator,
    merge_dicts,
    )

DESCRIPTION = """
pandda.analyse is a program for identifying "changed states" in crystallographic datasets.
It requires a (set of) "ground state" datasets to use as a reference to identify "events" in other datasets.
These "events" can be binding ligands, or any other change. PanDDA constructs a statistical model based on
the ground-state datasets and then compares each "test" datasets to see if there is a significant deviation. 
Once events have been identified, the occupancy of the event is estimated and a background correction is applied
to remove the confounding superposed ground state density. By doing this, PanDDA can very clearly reveal partial-
occupancy ligands (or other changes) that are not clear in the original electron density.
"""

program_banner = ModuleBanner(
    program = "pandda.analyse",
    description = DESCRIPTION,
    )

def standard_pandda(pandda_config):

    # Create output dictionary
    data_collator = DataCollator()
    
    # Options: maps a config to code abstraction
    pandda_options = options.Options(pandda_config)

    # Write parameters as computer-readable json
    pandda_options.dump_config_to_json(pandda_config)

    # This function is doing a lot of heavy lifting behind the scenes
    mcd = pandda_options.dataset_initialiser()

    data_collator({
        'output_files' : pandda_options.dataset_initialiser.output_files
        })

    # Make the main page 
    html_files = pandda_options.make_html_output()

    data_collator({
        'output_files' : {'html': html_files},
        })

    # Prepare to load maps (make grid, etc) -- also lots of heavy lifting
    map_loader = pandda_options.get_warped_map_loader(
        reference_dataset = mcd.reference_dataset,
        reference_map_getter = pandda_options.dataset_initialiser.reference_map_getter,
        dataset_map_getters = pandda_options.dataset_initialiser.map_getters,
        )

    # Main processing loop
    pandda_results = pandda_options.run_pandda_model(
        datasets = mcd.datasets,
        load_maps = map_loader,
        )

    data_collator(pandda_results)

    # Write output csvs, html, etc
    results_files = pandda_options.write_results(
        mcd = mcd,
        pandda_results = data_collator.data, # give all output so far
        )

    data_collator({
        'output_files' : results_files,
        })

    # Output results as a json so computer-readable
    pandda_options.dump_results_to_json(
        records_dict = data_collator.get_sorted(),
        )

    return data_collator.data

def run(args):

    ####################
    # Program Setup
    ####################

    # Get start time
    pandda_start_time = time.time()

    # Parse Config files and command line arguments
    from giant import jiffies
    working_phil = jiffies.extract_params_default(
        master_phil = phil.pandda_phil,
        args = args,
        blank_arg_prepend = None,
        home_scope = None,
    ).extract()

    # Maps options to code abstraction: Phil -> Config
    pandda_config = config.Config(working_phil)

    # Initialise logging
    logger = lg.setup_logging(
        name = __name__,
        log_file = str(pandda_config.output.out_dir / "pandda.log"),
        warning_handler_name = 'warnings',
        debug = False,
        )

    # Extract warning handler from logging objects
    from giant.logs import get_handler_recursive
    warning_handler = get_handler_recursive(
        logger = logger, 
        handler_name = 'warnings',
        )
    assert warning_handler is not None

    # Log banner
    logger(str(program_banner))

    # Write starup info
    from giant.phil import startup_parameters_logging
    startup_parameters_logging(
        output_directory = pandda_config.output.out_dir,
        master_phil = phil.pandda_phil,
        working_phil = working_phil,
        )

    try:

        standard_pandda(pandda_config)

    except Exception as e: 

        raise

    # Report warnings at the end
    warning_handler.report_all(logger_name=__name__)

    # We done it! 
    logger.heading("Pandda finished normally!", spacer=True)

def main():

    try:

        lg.setup_logging_basic(__name__)

        run(args=sys.argv[1:])

    except Exception as e:

        import traceback
        logger.subheading("Error - stack trace")
        logger(traceback.format_exc())#limit=10))
        logger.subheading("PanDDA exited with an error")
        logger("Error: %s", str(e))
        logger.subheading("PanDDA exited with an error")
        sys.exit(1)

if __name__ == "__main__":

    main()

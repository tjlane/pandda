import giant.logs as lg
logger = lg.getLogger(__name__)

import pathlib as pl

settings_phil = """
settings
    .help = "General Settings"
{
    cpus = 1
        .type = int
    verbose = False
        .type = bool
    plot_graphs = True
        .help = "Output graphs using matplotlib"
        .type = bool
    plotting {
        backend = 'agg'
            .help = "Backend to use in matplotlib"
            .type = str
    }
}
"""

image_phil = """
image
    .help = "Settings for image generation"
{
    size = 300x200
        .help = "image size of width x height (inches)"
        .type = str
    format = *png jpg svg
        .help = "select output image format
        .type = choice
}
"""

def log_running_parameters(
    params, 
    master_phil, 
    logger = None,
    ):

    if logger is None:
        logger = lg.getLogger(__name__)

    logger.heading('Processed parameters')
    logger(
        master_phil.format(
            python_object = params,
        ).as_str()
    )

    logger.heading('Non-default parameters')
    logger(
        master_phil.fetch_diff(
            source = master_phil.format(
                python_object = params,
            )
        ).as_str()
    )

def dump_config_to_json(
    config,
    output_path,
    ):

    record = {
        "data_dirs": str(config.input.data_dirs),
        "out_dir": str(config.output.out_dir),
    }

    import json
    json_string = json.dumps(record)

    with open(str(output_path), "w") as f:
        f.write(json_string)

def dump_params_to_eff(
    master_phil,
    working_phil,
    output_path,
    ):

    fmt_phil = master_phil.format(working_phil)

    with open(str(output_path), 'w') as fh:
        fh.write(
            fmt_phil.as_str()
        )

def startup_parameters_logging(
    output_directory,
    master_phil,
    working_phil,
    working_config,
    ):

    # Show input objects
    log_running_parameters(
        params = working_phil,
        master_phil = master_phil,
        logger = None,
    )

    # Write input params
    dump_config_to_json(
        config = working_config,
        output_path = str( pl.Path(output_directory) / "params.json"),
    )
    
    dump_params_to_eff(
        master_phil = master_phil,
        working_phil = working_phil,
        output_path = str( pl.Path(output_directory) / "params.eff"),
    )

import giant.logs as lg
logger = lg.getLogger(__name__)

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

def log_running_parameters(params, master_phil, logger=None):

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

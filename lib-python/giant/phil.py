
settings_phil = """
settings
    .help = "General Settings"
{
    cpus = 1
        .type = int
    verbose = True
        .type = bool
    plot_graphs = True
        .help = "Output graphs using matplotlib"
        .type = bool
}
"""

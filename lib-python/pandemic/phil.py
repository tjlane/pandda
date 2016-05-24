import libtbx.phil

pandemic_phil = libtbx.phil.parse("""
pandemic
    .help = "parameters to control the pandemic analysis"
{
    input
    {
        pandda_dir = './pandda'
            .type = path
            .multiple = False
    }
    output
    {
        out_dir = './pandemic'
            .type = path
            .multiple = False
    }
    params
    {
        resolution = 1.7
            .type = float
            .multiple = False
    }
}
include scope pandda.phil.pandda_phil
""", process_includes=True)

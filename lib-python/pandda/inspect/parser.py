import argparse

parser = argparse.ArgumentParser(
    prog = "pandda.inspect",
    )

parser.add_argument(
    "-m", 
    "--mode", 
    help = "mode to run in", 
    choices = [
        "events", "datasets",
        ],
    default = "events",
    type = str,
    )

parser.add_argument(
    "-p",
    "--pandda_directory",
    help = "input/output pandda directory", 
    default = ".",
    type = str,
    )

# parser.parse_args()
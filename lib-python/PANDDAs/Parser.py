from Bamboo.Common.Parser import argParser

def build_pandda_parser():
    """Parser for PANDDA Fragment Screen Analysis Program"""

    des = """Parser for use with the Pan-Dataset Density Analysis (PANDDA) Program"""

    # Initialise the main parser object
    p = argParser(prog='pandda', description=des, version='0.1')

    # Output Options
    p.add_argument('--outdir', '-o', metavar='DIR', default='./pandda', help='Output Directory')

    # Data input arguments
    p.add_argument('--ref-pdb', metavar='FILE', default='./reference.pdb', help='PDB File for the APO Structure of the Protein')
    p.add_argument('--ref-mtz', metavar='FILE', default='./reference.mtz', help='Associated MTZ File for the Reference Structure')
    p.add_argument('--data-dirs', metavar='DIR', required=True, help='Directorys Containing Diffraction Data, PDB Files and Smile Strings')
    p.add_argument('--pdb-style', metavar='STRING', default='final.pdb', help='Naming style of the pdb files in --data-dir. Can contain wildcards. i.e. "*-refine.pdb"')
    p.add_argument('--mtz-style', metavar='STRING', default=None, help='Naming style of the mtz files in --data-dir. Can contain wildcards. i.e. "*-refine.mtz". If None, looks for files named similarly to pdb_style')
    p.add_argument('--lig-style', metavar='STRING', default='*.cif', help='Naming style of the ligand files in --data-dir. Can contain wildcards. i.e. "lig-*.cif". If None, looks for files of the style "*.cif"')

    # Processing Arguments
    p.add_argument('--res-low', metavar='FLOAT', default=None, help='Highest resolution to perform analysis at')
    p.add_argument('--res-high', metavar='FLOAT', default=None, help='Lowest resolution to perform analysis at')
    p.add_argument('--res-jump', metavar='FLOAT', default=0.1, help='Increment in resolution cutoff between analyses')

    p.add_argument('--z-cut', metavar='FLOAT', default=2, help='Z-cutoff when looking for bound Z-blobs')
    p.add_argument('--min-blob-size', metavar='FLOAT', default=10, help='Minimum Z-blob volume (A)')

    # Add command arguments
    p.add_argument('--verbose', action='store_true', default=True, help='Set Verboseness')

    # Add Computer Settings
    p.add_argument('--max-memory', '-rmax', metavar='N', default=25, type=int, help='Maximum Amount of Memory to Use (GB) - Will fail if uses above this')
    p.add_argument('--cpus', '--num-cpus', metavar='N', default=1, type=int, help='Maximum Number of Processors to use')

    return p

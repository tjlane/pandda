from Bamboo.Common.Parser import argParser

def panddaParser():
    """Parser for PANDDA Fragment Screen Analysis Program"""

    des = """Parser for use with the Pan-Dataset Density Analysis (PANDDA) Program"""
    # Initialise the main parser object
    p = argParser(prog='pandda', description=des, version='1.0')
    # Add command arguments
    p.add_argument('--verbose', action='store_true', default=False, help='Set Verboseness')
    # Output Options
    p.add_argument('--outdir', '-o', metavar='DIR', default='./pandda', help='Output Directory')
    # Add Computer Settings
    p.add_argument('--maxmemory', '-RMAX', metavar='N', default=25, type=int, help='Maximum Amount of Memory to Use (GB) - Will fail if uses above this')
    p.add_argument('--keep-in-memory', '-K', action='store_true', default=False, help='Set whether the maps are held in memory - this will increase speed by use more memory')
    # Data input arguments
    p.add_argument('--ref-pdb', metavar='FILE', default='./reference.pdb', help='PDB File for the APO Structure of the Protein')
    p.add_argument('--ref-mtz', metavar='FILE', default='./reference.mtz', help='Associated MTZ File for the Reference Structure')
    p.add_argument('--data-dir', metavar='DIR', default='./data', help='Directory Containing Diffraction Data and Smile Strings')

    return p

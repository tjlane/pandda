from Bamboo.Common.Parser import argParser

def build_pandda_parser():
    """Parser for PANDDA Fragment Screen Analysis Program"""

    des = """Parser for use with the Pan-Dataset Density Analysis (PANDDA) Program"""
    # Initialise the main parser object
    p = argParser(prog='pandda', description=des, version='0.1')
    # Add command arguments
    p.add_argument('--verbose', action='store_true', default=False, help='Set Verboseness')
    # Output Options
    p.add_argument('--outdir', '-o', metavar='DIR', default='./pandda', help='Output Directory')
    # Add Computer Settings
    p.add_argument('--max-memory', '-rmax', metavar='N', default=25, type=int, help='Maximum Amount of Memory to Use (GB) - Will fail if uses above this')
    p.add_argument('--keep-maps-in-memory', '-K', action='store_true', default=False, help='Set whether the maps are held in memory - this will increase speed by use more memory')
    p.add_argument('--cpus', '--num-cpus', metavar='N', default=1, type=int, help='Maximum Number of Processors to use')

    # Add Algorithm Settings
    p.add_argument('--run_mol_subst', action='store_true', default=False, help='No PDBs supplied. Run refinement of reference model against all diffraction data')
    # Data input arguments
    p.add_argument('--ref-pdb', metavar='FILE', default='./reference.pdb', help='PDB File for the APO Structure of the Protein')
    p.add_argument('--ref-mtz', metavar='FILE', default='./reference.mtz', help='Associated MTZ File for the Reference Structure')
    p.add_argument('--data-dirs', metavar='DIR', default='./Processing/*', help='Directorys Containing Diffraction Data, PDB Files and Smile Strings')
    p.add_argument('--pdb-style', metavar='STRING', default='*/refine.pdb', help='Naming style of the pdb files in --data-dir. Can contain wildcards. i.e. "*/refine.pdb"')
    p.add_argument('--mtz-style', metavar='STRING', default=None, help='Naming style of the mtz files in --data-dir. Can contain wildcards. i.e. "*/refine.mtz". If None, looks for files named similarly to pdb_style')

    return p

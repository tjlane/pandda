#! /usr/local/python/python2.7.3-64bit/bin/python

from Bamboo.Common.Parser import argParser

def SoakParser():
    """Parser for Modelling Fragment-Soaking Experiments"""

    des = """Parser for use with Soaking Experiments"""
    # Initialise main parser object
    p = argParser(description=des)
    # Add arguments
    p.add_argument('--verbose','-v', action='store_true', default=False, help='Set Verboseness')
    # Reference PDB
    p.add_argument('--refpdb','--pdb', metavar='FILE', required=True, help='Reference Structure of the Protein (previously solved crystal form).')
    # Cif files for input protein
    p.add_argument('--cif', metavar='FILE', help='Other cif files for the structure')
    # CSV File
    p.add_argument('--csv', metavar='FILE', required=True, help='CSV Containing: {"DATASET LABEL", "SMILE STRING", "MTZ FILE LOCATION"}')
    # File Arguments
    p.add_argument('--outdir', '-o', metavar='DIR', default='./ProcessedFragmentSoak', help='Output Directory')
    # Meta
    p.add_argument('--experiment', '--exp', metavar='NAME', default='FragmentSoak', help ='Name describing the experiment i.e. PROTEIN_MODEL_1')

    return p








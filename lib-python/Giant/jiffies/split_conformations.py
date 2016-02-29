import os, sys, copy

import libtbx.phil

import numpy

import iotbx.pdb

from scitbx.array_family import flex

from Giant.Utils.pdb import strip_pdb_to_input

############################################################################

master_phil = libtbx.phil.parse("""
pdb = None
    .help = 'The ensemble of protein conformations (normally the unbound and bound structures)'
    .type = str
    .multiple = True

res = LIG,UNL
    .help = 'Residues to be define selected conformations (comma separated list of residue names, i.e. res=LIG or res=LIG,UNL)'
    .type = str
conf = None
    .help = 'Define selected conformations explicitly (comma separated list of conformer IDs, i.e. conf=A or conf=A,B)'
    .type = str

suffix = '.split.pdb'
    .help = 'output suffix for input files'
    .type = str

overwrite = False
    .type = bool
verbose = False
    .type = bool
""")

############################################################################

def run(params):

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'VALIDATING INPUT PARAMETERS'
    ######################################################################

    assert params.pdb, 'No PDB files given'
    assert params.res or params.conf, 'Both residue IDs and conformer IDs given: can only supply one'

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'SPLITTING ENSEMBLES''
    ######################################################################

    # Iterate through the input structures and extract the conformation
    for ensemble_file in params.pdb:

        # Create output file path
        output_file = os.path.split_ext(ensemble_file)[0] + params.suffix

        # Read in the ligand file and set each residue to the requested conformer
        ens_obj = strip_pdb_to_input(ensemble_file, remove_ter=True)

        # Check that ... something
        ens_obj.hierarchy.only_model()

        # Create a new copy of the structures
        new_ens = ens_obj.hierarchy.deep_copy()

        if params.res:
            confs_to_select =

        ######################################################################
        if params.verbose: print '===========================================>>>'
        if params.verbose: print 'PRUNING OUTPUT STRUCTURE'
        ######################################################################



        ######################################################################
        if params.verbose: print '===========================================>>>'
        if params.verbose: print 'WRITING OUTPUT STRUCTURE'
        ######################################################################

        # Update the atoms numbering
#        new_ens.atoms_reset_serial()
        # Write output file
        new_ens.write_pdb_file(output_file)

    return


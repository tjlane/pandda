import os, sys, copy

import libtbx.phil

import numpy

import iotbx.pdb

from scitbx.array_family import flex

from Giant.Structure.Utils import resolve_residue_id_clashes, normalise_occupancies, merge_hierarchies

############################################################################

systematic_letters = iotbx.pdb.systematic_chain_ids()

############################################################################

master_phil = libtbx.phil.parse("""
major = None
    .help = 'The major conformation of the protein (normally the unbound or reference structure)'
    .type = str

minor = None
    .help = 'The minor conformation of the protein (normally the bound or "interesting" structure)'
    .type = str

output = None
    .help = 'output pdb file'
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
    if params.verbose: print 'VALIDATING GIVEN PARAMETERS'
    ######################################################################

    assert params.major
    assert params.minor
    if not params.output: params.output='./output.pdb'
    if os.path.exists(params.output):
        if params.overwrite: os.remove(params.output)
        else: raise Exception('File already exists: {}'.format(params.output))

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'READING INPUT FILES'
    ######################################################################

    # Read in the ligand file and set each residue to the requested conformer
    maj_obj = iotbx.pdb.hierarchy.input(params.major)
    min_obj = iotbx.pdb.hierarchy.input(params.minor)

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'VALIDATING INPUT MODELS'
    ######################################################################

    # Check that ... something
    maj_obj.hierarchy.only_model()
    min_obj.hierarchy.only_model()

    # Create a new copy of the structures
    new_major = maj_obj.hierarchy.deep_copy()
    new_minor = min_obj.hierarchy.deep_copy()

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'RESOLVING RESIDUE ID CLASHES (BY ADDING MINOR RESIDUES TO NEW CHAINS)'
    if params.verbose: print '===========================================>>>'
    ######################################################################

    new_minor = resolve_residue_id_clashes(ref_hierarchy=new_major, mov_hierarchy=new_minor)

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'MERGE THE STRUCTURES'
    ######################################################################

    final_struct = merge_hierarchies(ref_hierarchy=new_major, mov_hierarchy=new_minor)

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'NORMALISING OUTPUT STRUCTURE'
    ######################################################################

    final_struct = normalise_occupancies(hierarchy=final_struct)

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'WRITING OUTPUT STRUCTURE'
    ######################################################################

    # Update the atoms numbering
    final_struct.atoms_reset_serial()
    # Write output file
    final_struct.write_pdb_file(params.output)

    return


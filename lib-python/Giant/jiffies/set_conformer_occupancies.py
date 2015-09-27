import os, sys, copy

import libtbx.phil

import numpy

import iotbx.pdb

from scitbx.array_family import flex

from Giant.Structure.Utils import normalise_occupancies, set_conformer_occupancy

############################################################################

master_phil = libtbx.phil.parse("""
conformer
    .multiple = True
    {
    id = None
        .type = str
    occ = 0.9
        .type = float
    }

input = None
    .type = path
output = None
    .type = path

minimum_occ = 0.0
    .type = float
maximum_occ = 1.0
    .type = float

overwrite = False
    .type = bool
verbose = False
    .type = bool
""")

############################################################################

def run(params):

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'SETTING OCCUPANCIES'
    ######################################################################

    assert params.conformer, 'No Conformers given!'
    assert params.input, 'No Input file given!'

    if not params.output: params.output='./occ-edited.pdb'
    if os.path.exists(params.output):
        if params.overwrite: os.remove(params.output)
        else: raise Exception('File already exists: {}'.format(params.output))

    for conf in params.conformer:
        print '{} -> {}'.format(conf.id,conf.occ)
        assert conf.id
        assert conf.occ

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'READING INPUT FILES'
    ######################################################################
    # Read in the ligand file and set each residue to the requested conformer
    prot_h = iotbx.pdb.hierarchy.input(params.input).hierarchy

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'VALIDATING INPUT MODELS'
    ######################################################################
    # Check that ... something
    prot_h.only_model()

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'SETTING OCCUPANCIES'
    ######################################################################

    for conf in params.conformer:
        if params.verbose: print 'Setting Conformer {} to {}'.format(conf.id, conf.occ)
        set_conformer_occupancy(hierarchy=prot_h, conf_id=conf.id, conf_occ=conf.occ, in_place=True)

    prot_h = normalise_occupancies( hierarchy=prot_h,
                                    exclude_conformers=[c.id for c in params.conformer],
                                    max_occ=params.maximum_occ,
                                    min_occ=params.minimum_occ  )

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'WRITING OUTPUT STRUCTURE'
    ######################################################################

    # Write output file
    prot_h.write_pdb_file(params.output)

    return


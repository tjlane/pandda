import os, sys, copy

import scipy

import libtbx.phil
import libtbx.easy_pickle

import iotbx.pdb

master_phil = libtbx.phil.parse("""
direction = toref *fromref
    .type = choice
file
{
    input  = None
        .type = path
    output = None
        .type = path
}
dataset_pickle = ./pickles/dataset.pickle
    .type = path
verbose = False
""")

def run(params):

    if params.dataset_pickle:
        assert os.path.exists(params.dataset_pickle), 'PICKLE DOES NOT EXIST: {!s}'.format(params.dataset_pickle)
    else:
        # Find the dataset pickle
        raise Exception('PLEASE PROVIDE DATASET PICKLE')

    if params.file.input:
        assert os.path.exists(params.file.input), 'FILE DOES NOT EXIST: {!s}'.format(params.file.input)
    else:
        # Find the dataset pickle
        raise Exception('PLEASE PROVIDE INPUT FILE')

    # Create output file if needed
    if not params.file.output:
        params.file.output = os.path.splitext(os.path.basename(params.file.input))[0] + '.{!s}.pdb'
        if params.direction == 'toref':   params.file.output = params.file.output.format('ref')
        else:                           params.file.output = params.file.output.format('native')
    assert not os.path.exists(params.file.output), 'FILE ALREADY EXISTS: {!s}'.format(params.file.output)

    if params.verbose: print 'TRANSFORMING: {!s}'.format(params.file.input)
    if params.verbose: print 'USING DATASET FROM: {!s}'.format(params.dataset_pickle)
    if params.verbose: print 'WRITING FILE TO: {!s}'.format(params.file.output)

    # Load the dataset handler from the pickle
    d_handler = libtbx.easy_pickle.load(params.dataset_pickle)

    # Load the ligand to be twiddled
    hier = iotbx.pdb.hierarchy.input(params.file.input).hierarchy

    if params.direction == 'toref':
        trans_points = d_handler.transform_to_reference(    points=hier.atoms().extract_xyz(),
                                                            method='global')
    elif params.direction == 'fromref':
        trans_points = d_handler.transform_from_reference(  points=hier.atoms().extract_xyz(),
                                                            method='global')

    # Create new hierarchy and write out
    new_hier = hier.deep_copy()
    new_hier.atoms().set_xyz(trans_points)
    new_hier.write_pdb_file(params.file.output)

    return


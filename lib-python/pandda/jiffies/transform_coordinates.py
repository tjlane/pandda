import os, sys, copy

import scipy

import libtbx.phil
import libtbx.easy_pickle

import iotbx.pdb

master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .type = path
    dataset_pickle = ./pickles/dataset.pickle
        .type = path
}
output {
    pdb = None
        .type = path
}
options {
    direction = toref *fromref
        .type = choice
}
settings {
    overwrite = False
        .type = bool
    verbose = False
        .type = bool
}
""")

def run(params):

    assert params.input.pdb is not None, 'No input pdb file provided (params.input.pdb).'
    assert os.path.exists(params.input.pdb), 'Input pdb file does not exist: {}'.format(params.input.pdb)
    assert params.input.dataset_pickle is not None, 'No dataset pickle file provided (params.input.dataset_pickle).'
    assert os.path.exists(params.input.dataset_pickle), 'Dataset pickle file does not exist: {!s}'.format(params.input.dataset_pickle)
    assert params.output.pdb is not None, 'No output pdb file provided (params.output.pdb).'
    assert not os.path.exists(params.output.pdb), 'Output pdb file does not exist: {}'.format(params.output.pdb)

    if params.verbose:
        print 'Input file: {!s}'.format(params.input.pdb)
        print 'Output file: {!s}'.format(params.output.pdb)
        print 'Transforming using the dataset from: {!s}'.format(params.input.dataset_pickle)

    # Load the dataset handler from the pickle
    dataset = libtbx.easy_pickle.load(params.input.dataset_pickle)
    # Load the ligand to be twiddled
    inpt = iotbx.pdb.hierarchy.input(params.input.pdb)
    hier = inpt.hierarchy

    # Transform the coordinates
    if params.options.direction == 'toref':
        hier.atoms().set_xyz(dataset.alignment.nat2ref(coordinates=hier.atoms().extract_xyz()))
    else:
        hier.atoms().set_xyz(dataset.alignment.ref2nat(coordinates=hier.atoms().extract_xyz()))

    # Try to add symmetry to the output file (if sending to native frame)
    crystal_symmetry = None
    if params.options.direction == 'fromref':
        try:    crystal_symmetry = dataset.data.summary.symmetry
        except: pass

    # Create new hierarchy and write out
    hier.write_pdb_file(    file_name        = params.output.pdb,
                            crystal_symmetry = crystal_symmetry     )

    return

#######################################

if __name__ == '__main__':
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)

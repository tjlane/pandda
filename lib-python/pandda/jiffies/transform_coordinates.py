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
overwrite = False
    .type = bool
verbose = False
    .type = bool
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
    if os.path.exists(params.file.output):
        if params.overwrite: os.remove(params.file.output)
        else: raise Exception('FILE ALREADY EXISTS: {!s}'.format(params.file.output))

    if params.verbose: print 'TRANSFORMING: {!s}'.format(params.file.input)
    if params.verbose: print 'USING DATASET FROM: {!s}'.format(params.dataset_pickle)
    if params.verbose: print 'WRITING FILE TO: {!s}'.format(params.file.output)

    # Load the dataset handler from the pickle
    d_handler = libtbx.easy_pickle.load(params.dataset_pickle)

    # Load the ligand to be twiddled
    hier = iotbx.pdb.hierarchy.input(params.file.input).hierarchy

    if d_handler.local_alignment_transforms:
        method = 'local'
        if params.direction == 'toref':
            map_hier = d_handler.new_structure().hierarchy
        else:
            map_hier = d_handler.new_structure().hierarchy
            map_hier.atoms().set_xyz(d_handler.transform_to_reference(points=map_hier.atoms().extract_xyz(), method=method))

        mappings = d_handler.find_nearest_calpha(   points = hier.atoms().extract_xyz(),
                                                    hierarchy = map_hier    )
    else:
        method = 'global'
        mappings = None

    if params.direction == 'toref':
        trans_points = d_handler.transform_to_reference(    points = hier.atoms().extract_xyz(),
                                                            method = method,
                                                            point_mappings = mappings   )
    elif params.direction == 'fromref':
        trans_points = d_handler.transform_from_reference(  points = hier.atoms().extract_xyz(),
                                                            method = method,
                                                            point_mappings = mappings   )

    # Try to add symmetry to the output file
    try:    crystal_symmetry = d_handler.mtz_summary.symmetry
    except: crystal_symmetry = None

    # Create new hierarchy and write out
    new_hier = hier.deep_copy()
    new_hier.atoms().set_xyz(trans_points)
    new_hier.write_pdb_file(    file_name        = params.file.output,
                                crystal_symmetry = crystal_symmetry     )

    return

#######################################

if __name__ == '__main__':
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)

import os, sys, copy, glob

import libtbx.phil

import numpy

from PANDDAs.jiffies import transform_coordinates, export_files
from Giant.jiffies import merge_conformations

############################################################################

master_phil = libtbx.phil.parse("""
process_and_export {
    merge_and_export_fitted_structures = True
        .help = 'export fitted ligand structures'
        .type = bool

    transform_and_export_defaults = True
        .type = bool

    pdbs_to_transform_and_export = None
        .type = path
        .multiple = True

    maps_to_transform_and_export = None
        .type = path
        .multiple = True
}

include scope PANDDAs.jiffies.export_files.master_phil

""",
process_includes = True)

############################################################################

def prepend_prefix_to_basename(prefix, path):
    return os.path.join(os.path.dirname(path),prefix+os.path.basename(path))

def run(params):

    # Prefix to append to the transformed files
    temp_prefix = 'temp-'
    temp_files_to_delete = []

    # Default files paths
    aligned_template = '{!s}-aligned.pdb'
    fitted_template = 'modelled_structures/fitted-current.pdb'
    merged_template = '{!s}-pandda-model.pdb'

    # Merge, Transform and Export fitted Structures?
    if params.process_and_export.merge_and_export_fitted_structures:
        # Merged structure (minor + major conformations)
        params.process_and_export.pdbs_to_transform_and_export.append(merged_template)

    # Transform and Export default files?
    if params.process_and_export.transform_and_export_defaults:
        # Un-merged structure (major conformation/reference)
        params.process_and_export.pdbs_to_transform_and_export.append(aligned_template)
        # Fitted structure (minor conformation)
        params.process_and_export.pdbs_to_transform_and_export.append(fitted_template)

    # Add files to be transformed to those to be exported
    [params.export.files_to_export.append(prepend_prefix_to_basename(temp_prefix, f).format('*')) for f in params.process_and_export.pdbs_to_transform_and_export]

    ############################################################################

    # Find the dataset directories to be exported
    export_dirs = sorted(glob.glob(os.path.join(params.input.pandda_dir, params.export.datasets_to_export, '*')))

    # Merge the fitted structures
    for e_dir in export_dirs:

        dir_tag = os.path.basename(e_dir)
        print '==========================================+>'
        print 'Processing Directory:', dir_tag

        ############################################################################

        # Check to see if this folder should be skipped (export fitted folders only)
        if params.export.required_file_for_export and (not os.path.exists(os.path.join(e_dir, params.export.required_file_for_export))):
            print 'NO FITTED MODEL - SKIPPING: {}'.format(e_dir)
            continue

        ############################################################################

        # MERGING WITH ALIGNED STRUCTURE

        if params.process_and_export.merge_and_export_fitted_structures:
            # Extract parameters for the merging and set them
            merging_params = merge_conformations.master_phil.extract()
            merging_params.major =  os.path.join(e_dir, aligned_template.format(dir_tag))
            merging_params.minor =  os.path.join(e_dir, fitted_template)
            merging_params.output = os.path.join(e_dir, merged_template.format(dir_tag))
            merging_params.verbose = params.verbose

            # Merge the structures
            print '=========================+>'
            if os.path.exists(merging_params.minor):
                print 'MERGING FITTED STRUCTURE WITH APO STRUCTURE'
                if params.verbose: print '=========================+>'
                if params.verbose: print 'MAJOR:', merging_params.major
                if params.verbose: print 'MINOR:', merging_params.minor
                if params.verbose: print 'MERGED:', merging_params.output
                merge_conformations.run(merging_params)
            else:
                print 'NOT MERGING - INPUT DOES NOT EXIST:', merging_params.minor

        ############################################################################

        # TRANSFORMING

        # Extract parameters for the transformation and set them
        transform_params = transform_coordinates.master_phil.extract()
        transform_params.dataset_pickle = os.path.join(e_dir, 'pickles', 'dataset.pickle')
        transform_params.verbose = params.verbose

        # Transform the files in the list
        for template_file in params.process_and_export.pdbs_to_transform_and_export:

            # Fill in the name of the file from the template
            dataset_file = template_file.format(dir_tag)

            transform_params.file.input  = os.path.join(e_dir, dataset_file)
            transform_params.file.output = os.path.join(e_dir, prepend_prefix_to_basename(temp_prefix, dataset_file))

            # Transform the structures
            print '=========================+>'
            if os.path.exists(transform_params.file.input):
                print 'MAPPING FILES TO NATIVE FRAME'
                if params.verbose: print '=========================+>'
                if params.verbose: print 'MAPPING INPUT:', transform_params.file.input
                if params.verbose: print 'MAPPING OUTPUT:', transform_params.file.output
                transform_coordinates.run(transform_params)
                # Will want to delete this later
                temp_files_to_delete.append(transform_params.file.output)
            else:
                print 'NOT MAPPING - INPUT DOES NOT EXIST:', transform_params.file.input

    # Export the pandda folder
    export_files.run(params)

    # Go through the files and remove the temp_prefix prefix
    print '=========================+>'
    print 'Renaming Temporary Files:'
    for root, dirs, files in os.walk(params.output.out_dir):
        for f in files:
            if f.startswith(temp_prefix):
                new_f = os.path.join(root, f[len(temp_prefix):])
                if params.verbose: print root, f, '->', new_f
                assert not os.path.exists(new_f), 'New file should not exist: {!s}'.format(new_f)
                os.rename(os.path.join(root, f), new_f)

    # Go through the created files and remove them
    print '=========================+>'
    print 'Deleting Temporary Files:'
    for del_f in temp_files_to_delete:
        assert os.path.isfile(del_f), 'This should be a file: {!s}'.format(del_f)
        if params.verbose: print del_f, '->', '-'
        os.remove(del_f)


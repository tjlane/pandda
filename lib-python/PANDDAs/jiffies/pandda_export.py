import os, sys, copy, glob

import libtbx.phil

import numpy

from PANDDAs.jiffies import transform_coordinates, export_files
from Giant.jiffies import merge_conformations, create_occupancy_params

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

    generate_occupancy_groupings = True
        .type = bool
}

templates {
    temp_prefix = 'temp-'
        .help = 'Prefix to append to temporary files'
        .type = str
    aligned_template = '{!s}-aligned.pdb'
        .help = 'Un-merged structure (major conformation/reference)'
        .type = str
    fitted_template = 'modelled_structures/pandda-model.pdb'
        .help = 'Fitted structure (minor conformation)'
        .type = str
    merged_template = '{!s}-combined-model.pdb'
        .help = 'Merged structure (minor + major conformations)'
        .type = str
    refmac_refinement = 'refmac_refine.params'
        .help = 'refmac occupancy groupings/refinement parameters'
        .type = str
    phenix_refinement = 'phenix_refine.params'
        .help = 'phenix occupancy groupings/refinement parameters'
        .type = str
}

include scope PANDDAs.jiffies.export_files.master_phil

""",
process_includes = True)

############################################################################

def prepend_prefix_to_basename(prefix, path):
    return os.path.join(os.path.dirname(path),prefix+os.path.basename(path))

def set_defaults(params):

    # Transform and Export default files?
    if params.process_and_export.transform_and_export_defaults:
        params.process_and_export.pdbs_to_transform_and_export.append(params.templates.aligned_template)
        params.process_and_export.pdbs_to_transform_and_export.append(params.templates.fitted_template)

    # Merge, Transform and Export fitted Structures?
    if params.process_and_export.merge_and_export_fitted_structures:
        params.process_and_export.pdbs_to_transform_and_export.append(params.templates.merged_template)

    # Add default export files
    if params.export.export_defaults: export_files.set_defaults(params)
    # Export occupancy groupings + refinement parameter files?
    if params.process_and_export.generate_occupancy_groupings:
        params.export.files_to_export.append(params.templates.refmac_refinement)
        params.export.files_to_export.append(params.templates.phenix_refinement)
    # Add files to be transformed to those to be exported
    for f in params.process_and_export.pdbs_to_transform_and_export:
        params.export.files_to_export.append(prepend_prefix_to_basename(params.templates.temp_prefix, f).format('*'))

def process_and_export_folder(dir, params):
    """Merge structures, transform them and export a subset of a folders contents"""

    dir_name = os.path.basename(dir)
    print '\n\n==========================================+>'
    print 'Processing Directory:', dir_name
    print '==========================================+>\n'

    ############################################################################

    # Check to see if this folder should be skipped (export fitted folders only)
    if params.export.required_file_for_export and (not os.path.exists(os.path.join(dir, params.export.required_file_for_export))):
        print 'NO FITTED MODEL - SKIPPING: {}'.format(dir)
        return

    ############################################################################

    # MERGING WITH ALIGNED STRUCTURE

    if params.process_and_export.merge_and_export_fitted_structures:
        # Extract parameters for the merging and set them
        merging_params = merge_conformations.master_phil.extract()
        merging_params.major =  os.path.join(dir, params.templates.aligned_template.format(dir_name))
        merging_params.minor =  os.path.join(dir, params.templates.fitted_template)
        merging_params.output = os.path.join(dir, params.templates.merged_template.format(dir_name))
        merging_params.overwrite = params.overwrite
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

    print '=========================+>'
    print 'TRANSFORMING FILES FROM PANDDAS FRAME TO NATIVE FRAME'

    # Extract parameters for the transformation and set them
    transform_params = transform_coordinates.master_phil.extract()
    transform_params.dataset_pickle = os.path.join(dir, 'pickles', 'dataset.pickle')
    transform_params.overwrite = params.overwrite
    transform_params.verbose = params.verbose

    # Store files to be deleted
    temp_files_to_delete = []

    # Transform the files in the list
    for template_file in params.process_and_export.pdbs_to_transform_and_export:

        # Fill in the name of the file from the template
        dataset_file = template_file.format(dir_name)
        input_files  = glob.glob(os.path.join(dir, dataset_file))

        for f in input_files:
            transform_params.file.input = f
            transform_params.file.output = prepend_prefix_to_basename(params.templates.temp_prefix, f)

            if os.path.exists(transform_params.file.input):
                # Transform the structures
                transform_coordinates.run(transform_params)
                # Will want to delete this later
                temp_files_to_delete.append(transform_params.file.output)
            else:
                print 'NOT MAPPING - INPUT DOES NOT EXIST: {}'.format(transform_params.file.input)

    ############################################################################

    print '=========================+>'
    print 'GENERATING OCCUPANCY REFINEMENT PARAMETERS'

    # Extract parameters for the occupancy parameter generation and set them
    occupancy_params = create_occupancy_params.master_phil.extract()
    occupancy_params.pdb = merging_params.output
    occupancy_params.refmac_occ_out = os.path.join(dir, params.templates.refmac_refinement)
    occupancy_params.phenix_occ_out = os.path.join(dir, params.templates.phenix_refinement)
    occupancy_params.overwrite = params.overwrite
    occupancy_params.verbose = params.verbose

    if params.verbose: print 'ANALYSING FILE: {}'.format(occupancy_params.pdb)
    if params.verbose: print 'OUTPUTTING REFMAC SETTINGS TO: {}'.format(occupancy_params.refmac_occ_out)
    if params.verbose: print 'OUTPUTTING PHENIX SETTINGS TO: {}'.format(occupancy_params.phenix_occ_out)

    try: create_occupancy_params.run(occupancy_params)
    except Exception as e:
        print 'OCCUPANCY PARAMETER GENERATION FAILED: {} - {}'.format(e, e.message)

    ############################################################################

    print '=========================+>'
    print 'EXPORTING FOLDER: {}'.format(dir)

    # Export the pandda folder
    export_files.export_folder(dir=dir, params=params)

    ############################################################################

    # Go through the files and remove the temp_prefix prefix
    print '=========================+>'
    print 'Renaming Temporary Files'
    for root, dirs, files in os.walk(os.path.join(params.output.out_dir, dir_name)):
        for f in files:
            if f.startswith(params.templates.temp_prefix):
                new_f = os.path.join(root, f[len(params.templates.temp_prefix):])
                if params.verbose: print root, f, '->', new_f
                if os.path.exists(new_f):
                    if params.overwrite: os.remove(new_f)
                    else: raise Exception('File already exists: {}'.format(new_f))
                os.rename(os.path.join(root, f), new_f)

    # Go through the temporary files and remove them
    print '=========================+>'
    print 'Deleting Temporary Files'
    for del_f in temp_files_to_delete:
        assert os.path.isfile(del_f), 'This should be a file: {!s}'.format(del_f)
        if params.verbose: print del_f, '->', '-'
        os.remove(del_f)

def run(params):

    # Default files paths
    set_defaults(params)

    ############################################################################

    # Find the dataset directories to be exported
    export_dirs = sorted(glob.glob(os.path.join(params.input.pandda_dir, params.export.datasets_to_export, '*')))
    assert export_dirs, 'No Export Directories Found'

    # Create output directory
    if not os.path.exists(params.output.out_dir): os.mkdir(params.output.out_dir)

    # Merge the fitted structures
    for dir in export_dirs:
        process_and_export_folder(dir=dir, params=params)




import os, sys, copy, glob
import shutil

import libtbx.phil

import numpy

from pandda.jiffies import transform_coordinates
from pandda.constants import PanddaDatasetFilenames
from giant.jiffies import merge_conformations, create_occupancy_params

############################################################################

blank_arg_prepend = {None:'select_dir='}

master_phil = libtbx.phil.parse("""
input {
    pandda_dir = ./pandda
        .help = 'Path to the pandda directory to export files from'
        .type = str
    select_dir = None
        .help = 'Manually specify directories to export'
        .type = str
        .multiple = True
}
output {
    export_dir = ./pandda-export
        .type = str

    dir_prefix = ''
        .help = 'Prefix to be added to each output (sub)directory'
        .type = str

    file_prefix = ''
        .help = 'Prefix to be added to each output file'
        .type = str

    generate_occupancy_groupings = False
        .type = bool
}
export {
    datasets_to_export = *interesting_datasets processed_datasets
        .type = choice

    files_to_export = None
        .type = path
        .multiple = True

    export_defaults = True
        .type = bool

    export_ligands = True
        .type = bool

    required_file_for_export = *model None
        .help = 'only export folder if this file exists (set to None to turn off)'
        .type = choice
}
process_and_export {
    merge_and_export_fitted_structures = True
        .help = 'export fitted ligand structures'
        .type = bool

    transform_and_export_defaults = True
        .type = bool

    pdbs_to_transform_and_export = None
        .type = path
        .multiple = True
}
templates {
    temp_prefix = 'temp-'
        .help = 'Prefix to append to temporary files'
        .type = str
    refmac_refinement = 'occ_refmac.params'
        .help = 'refmac occupancy groupings/refinement parameters'
        .type = str
    phenix_refinement = 'occ_phenix.params'
        .help = 'phenix occupancy groupings/refinement parameters'
        .type = str
}

overwrite = True
    .type = bool
verbose = False
    .type = bool

""")

############################################################################

def prepend_prefix_to_basename(prefix, path):
    return os.path.join(os.path.dirname(path),prefix+os.path.basename(path))

def set_defaults(params):

    # Add default export files
    if params.export.export_defaults:
        # Add default files to list
        params.export.files_to_export.append('*-pandda-input.pdb')
        params.export.files_to_export.append('*-pandda-input.mtz')
        params.export.files_to_export.append('*.native.ccp4')
    # Export ligands files?
    if params.export.export_ligands:
        params.export.files_to_export.append('ligand_files/*.cif')
        params.export.files_to_export.append('ligand_files/*.pdb')

    # Transform and Export default files?
    if params.process_and_export.transform_and_export_defaults:
#        params.process_and_export.pdbs_to_transform_and_export.append(PanddaDatasetFilenames.aligned_structure)
        params.process_and_export.pdbs_to_transform_and_export.append(os.path.join('modelled_structures',PanddaDatasetFilenames.modelled_structure))
    # Merge, Transform and Export fitted Structures?
    if params.process_and_export.merge_and_export_fitted_structures:
        params.process_and_export.pdbs_to_transform_and_export.append(PanddaDatasetFilenames.ensemble_structure)

    # Export occupancy groupings + refinement parameter files?
    if params.output.generate_occupancy_groupings:
        params.export.files_to_export.append(params.templates.refmac_refinement)
        params.export.files_to_export.append(params.templates.phenix_refinement)
    # Add files to be transformed to those to be exported
    for f in params.process_and_export.pdbs_to_transform_and_export:
        params.export.files_to_export.append(prepend_prefix_to_basename(params.templates.temp_prefix, f).format('*'))

def export_folder(dir, params):
    """Export a subset of a folders contents"""

    # Extract folder name
    dir_name = os.path.basename(dir)
    if params.verbose: print 'Processing Folder: {!s}'.format(dir_name)
    # Create output dir
    export_dir = os.path.join(params.output.export_dir, params.output.dir_prefix+dir_name)
    if not os.path.exists(export_dir):  os.mkdir(export_dir)
    # Export files
    if params.verbose: print 'Exporting files from {!s} to {!s}'.format(dir, export_dir)
    for proc_file_template in params.export.files_to_export:
        proc_files = glob.glob(os.path.join(dir, proc_file_template))
        for proc_file in proc_files:
            # Check that the file exists
            if not os.path.exists(proc_file):
                if params.verbose: print 'FILE DOES NOT EXIST: {!s}'.format(proc_file)
                continue
            # Exported file path
            export_file = os.path.join(export_dir, params.output.file_prefix+os.path.basename(proc_file))
            if params.verbose: print 'Exporting {!s} to {!s}'.format(proc_file, export_file)
            # Check to see if file already exists and delete if overwrite
            if os.path.exists(export_file):
                if params.overwrite: os.remove(export_file)
                else: raise Exception('File Already Exists: {}'.format(export_file))
            shutil.copy(proc_file, export_file)

def process_and_export_folder(dir, params):
    """Merge structures, transform them and export a subset of a folders contents"""

    dir_name = os.path.basename(dir)
    print '\n\n==========================================+>'
    print 'Processing Directory:', dir_name
    print '==========================================+>\n'

    ############################################################################

    # Check to see if this folder should be skipped (export fitted folders only)
    if params.export.required_file_for_export:
        if params.export.required_file_for_export == 'model':
            if not os.path.exists(os.path.join(dir, 'modelled_structures', PanddaDatasetFilenames.modelled_structure.format(dir_name))):
                print 'NO FITTED MODEL - SKIPPING: {}'.format(dir)
                return

    ############################################################################

    # MERGING WITH ALIGNED STRUCTURE

    if params.process_and_export.merge_and_export_fitted_structures:
        # Extract parameters for the merging and set them
        merging_params = merge_conformations.master_phil.extract()
        merging_params.major =  os.path.join(dir, PanddaDatasetFilenames.aligned_structure.format(dir_name))
        merging_params.minor =  os.path.join(dir, 'modelled_structures', PanddaDatasetFilenames.modelled_structure.format(dir_name))
        merging_params.output = os.path.join(dir, PanddaDatasetFilenames.ensemble_structure.format(dir_name))
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

    if params.output.generate_occupancy_groupings:
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
    export_folder(dir=dir, params=params)

    ############################################################################

    # Go through the files and remove the temp_prefix prefix
    print '=========================+>'
    print 'Renaming Temporary Files'
    for root, dirs, files in os.walk(os.path.join(params.output.export_dir, dir_name)):
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

    # Change paths to absolute paths
    params.input.pandda_dir  = os.path.abspath(params.input.pandda_dir)
    params.output.export_dir = os.path.abspath(params.output.export_dir)

    # Must be in the pandda directory (pandda objects use relative paths)
    os.chdir(params.input.pandda_dir)

    ############################################################################

    # Find the dataset directories to be exported
    if params.input.select_dir:
        export_dirs = sorted([os.path.join(params.input.pandda_dir, params.export.datasets_to_export, p) for p in params.input.select_dir])
        export_dirs = [p for p in export_dirs if os.path.exists(p)]
        print 'Exporting:\n\t', '\n\t'.join(export_dirs)
    else:
        export_dirs = sorted(glob.glob(os.path.join(params.input.pandda_dir, params.export.datasets_to_export, '*')))
    assert export_dirs, 'No Export Directories Found'

    # Create output directory
    if not os.path.exists(params.output.export_dir): os.mkdir(params.output.export_dir)

    # Merge the fitted structures
    for dir in export_dirs:
        process_and_export_folder(dir=dir, params=params)

#######################################

if __name__ == '__main__':
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)

import os, sys, shutil, copy, glob

import libtbx.phil

############################################################################

master_phil = libtbx.phil.parse("""
export {
    datasets_to_export = *interesting_datasets processed_datasets
        .type = choice

    files_to_export = None
        .type = path
        .multiple = True

    export_defaults = True
        .type = bool

    required_file_for_export = 'modelled_structures/pandda-model.pdb'
        .help = 'only export folder if this file exists (set to None to turn off)'
        .type = path
}
input {
    pandda_dir = ./pandda
        .help = 'Path to the pandda directory to export files from'
        .type = str
}
output {
    out_dir = ./exported_pandda
        .type = str

    dir_prefix = ''
        .help = 'Prefix to be added to each output (sub)directory'
        .type = str

    file_prefix = ''
        .help = 'Prefix to be added to each output file'
        .type = str
}

overwrite = True
    .type = bool
verbose = False
    .type = bool
""")

############################################################################

def export_folder(dir, params):
    """Export a subset of a folders contents"""

    # Extract folder name
    dir_name = os.path.basename(dir)
    if params.verbose: print 'Processing Folder: {!s}'.format(dir_name)
    # Check to see if this folder should be skipped (export fitted folders only)
    if params.export.required_file_for_export and (not os.path.exists(os.path.join(dir, params.export.required_file_for_export))):
        print 'REQUIRED FILE MISSING - SKIPPING: {}'.format(dir)
        return
    # Create output dir
    export_dir = os.path.join(params.output.out_dir, params.output.dir_prefix+dir_name)
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

def set_defaults(params):
    # Add default files to list
    params.export.files_to_export.append('*-input.pdb')
    params.export.files_to_export.append('*-input.mtz')
    params.export.files_to_export.append('ligand_files/*.cif')
    params.export.files_to_export.append('ligand_files/*.pdb')
    params.export.files_to_export.append('*.native.ccp4')
    return

def run(params):
    # Create output directory
    if not os.path.exists(params.output.out_dir):   os.mkdir(params.output.out_dir)
    # Check input directory
    assert os.path.exists(params.input.pandda_dir)
    # Add defaults
    if params.export.export_defaults: set_defaults(params)
    # Export files
    print 'Exporting Files:\n\t', '\n\t'.join(params.export.files_to_export)
    for dir in sorted(glob.glob(os.path.join(params.input.pandda_dir, params.export.datasets_to_export, '*'))):
        export_folder(dir=dir, params=params)
    return


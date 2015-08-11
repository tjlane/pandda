import os, sys, shutil, copy, glob

import libtbx.phil

############################################################################

master_phil = libtbx.phil.parse("""
files_to_export = None
    .type = path
    .multiple = True

export_defaults = True
    .type = bool

datasets_to_export = *interesting_datasets processed_datasets
    .type = choice

out_dir = ./exported_pandda
    .type = str

pandda_dir = ./pandda
    .help = 'Path to the pandda directory to export files from'
    .type = str

dir_prefix = ''
    .help = 'Prefix to be added to each output (sub)directory'
    .type = str

file_prefix = ''
    .help = 'Prefix to be added to each output file'
    .type = str

verbose = False
    .type = bool
""")

############################################################################

def run(params):

    # Create output directory
    if not os.path.exists(params.out_dir):   os.mkdir(params.out_dir)
    # Check input directory
    assert os.path.exists(params.pandda_dir)

    # Add defaults
    if params.export_defaults:
        params.files_to_export.append('*-input.pdb')
        params.files_to_export.append('*-input.mtz')
        params.files_to_export.append('ligand_files/*.cif')
        params.files_to_export.append('ligand_files/*.pdb')

    ############################################################################

    print 'Exporting Files:\n\t', '\n\t'.join(params.files_to_export)

    for proc_dir in sorted(glob.glob(os.path.join(params.pandda_dir, params.datasets_to_export, '*'))):

        dir_tag = os.path.basename(proc_dir)
        print 'Processing Dataset: {!s}'.format(dir_tag)

        export_dir = os.path.join(params.out_dir, params.dir_prefix+dir_tag)
        print 'Exporting files from {!s} to {!s}'.format(proc_dir, export_dir)

        if os.path.exists(export_dir):
            print 'DIRECTORY ALREADY EXISTS'
            continue
        else:
            os.mkdir(export_dir)

        for proc_file_template in params.files_to_export:
            proc_files = glob.glob(os.path.join(proc_dir, proc_file_template))
            if proc_files:
                for proc_file in proc_files:
                    # Check that the file exists
                    if not os.path.exists(proc_file):
                        print 'FILE DOES NOT EXIST: {!s}'.format(proc_file)
                        continue
                    # Exported file path
                    export_file = os.path.join(export_dir, params.file_prefix+os.path.basename(proc_file))
                    if params.verbose: print 'Exporting {!s} to {!s}'.format(proc_file, export_file)
                    shutil.copy(proc_file, export_file)

    return


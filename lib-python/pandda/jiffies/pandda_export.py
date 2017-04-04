import os, sys, copy, glob
import shutil

import libtbx.phil

import numpy

from bamboo.common.logs import Log
from pandda.constants import PanddaDatasetFilenames
from giant.jiffies import merge_conformations, make_restraints

############################################################################

PROGRAM = 'pandda.export'

DESCRIPTION = """
Extract the structures and maps from a pandda analysis and prepare them for refinement
"""

############################################################################

blank_arg_prepend = {None:'select_datasets='}

master_phil = libtbx.phil.parse("""
input {
    pandda_dir = None
        .help = 'Path to the pandda directory to export files from'
        .type = str
    select_datasets = None
        .help = 'Manually specify dataset labels to export - multiple may be provided as a comma separated list.'
        .type = str
        .multiple = True
}
output {
    export_dir = pandda-export
        .type = str
    log_file = pandda-export.log
        .type = str
    dir_prefix = ''
        .help = 'Prefix to be added to each output (sub)directory'
        .type = str
    file_prefix = ''
        .help = 'Prefix to be added to each output file'
        .type = str
}
options {
    required_file_for_export = *model none
        .help = 'only export folder if this file exists (otherwise export all folders)'
        .type = choice
    generate_restraints = True
        .type = bool
}
settings {
    overwrite = True
        .type = bool
    verbose = False
        .type = bool
    cpus = 1
        .type = int
}

""")

############################################################################

def prepend_prefix_to_basename(prefix, path):
    return os.path.join(os.path.dirname(path),prefix+os.path.basename(path))

def get_file_list(dir):

    file_list = []
    file_list.extend(glob.glob(os.path.join(dir,PanddaDatasetFilenames.input_model.format('*'))))
    file_list.extend(glob.glob(os.path.join(dir,PanddaDatasetFilenames.input_data.format('*'))))
    file_list.extend(glob.glob(os.path.join(dir,PanddaDatasetFilenames.native_obs_map.format('*'))))
    file_list.extend(glob.glob(os.path.join(dir,PanddaDatasetFilenames.native_z_map.format('*'))))
    file_list.extend(glob.glob(os.path.join(dir,PanddaDatasetFilenames.native_event_map.format('*','*','*'))))
    file_list.extend(glob.glob(os.path.join(dir,PanddaDatasetFilenames.native_mean_map.format('*'))))
    file_list.extend(glob.glob(os.path.join(dir,'ligand_files','*.pdb')))
    file_list.extend(glob.glob(os.path.join(dir,'ligand_files','*.cif')))
    file_list.extend(glob.glob(os.path.join(dir,'modelled_structures',PanddaDatasetFilenames.modelled_structure.format('*'))))
    file_list.extend(glob.glob(os.path.join(dir,'modelled_structures',PanddaDatasetFilenames.ensemble_structure.format('*'))))

    return file_list

def export_folder(dir, params, log=Log()):
    """Export a subset of a folders contents"""

    # Extract folder name and report
    dir_name = os.path.basename(dir)
    # Get the file list for this folder
    file_list = get_file_list(dir=dir)
    # Create output dir
    exp_dir = os.path.join(params.output.export_dir, params.output.dir_prefix+dir_name)
    if not os.path.exists(exp_dir):
        os.mkdir(exp_dir)
    # Report
    log.bar()
    log('Exporting \n\tfrom {!s} \n\t  to {!s}'.format(dir, exp_dir))
    log.bar()
    log('Exporting files:')
    for f in file_list:
        log('\t'+os.path.relpath(f, start=dir))
    log.bar()
    # Export files
    for proc_file in file_list:
        # Check that the file exists
        if not os.path.exists(proc_file):
            log('FILE DOES NOT EXIST: {!s}'.format(proc_file))
            continue
        # Exported file path
        export_file = os.path.join(exp_dir, params.output.file_prefix+os.path.basename(proc_file))
        if params.settings.verbose:
            log('Copying {!s}\n     to {!s}'.format(proc_file, export_file))
        # Check to see if file already exists and delete if overwrite
        if os.path.exists(export_file):
            if params.settings.overwrite:
                os.remove(export_file)
            else:
                raise Exception('File already exists: {}. Need to set overwrite=True to continue.'.format(export_file))
        shutil.copy(proc_file, export_file)

    return exp_dir

def process_and_export_folder(dir, params, log=Log()):
    """Merge structures, transform them and export a subset of a folders contents"""

    dir_name = os.path.basename(dir)
    log.heading('Processing directory: {}'.format(dir_name), spacer=True)

    # Check to see if this folder should be skipped (export fitted folders only)
    if params.options.required_file_for_export == 'model':
        if not os.path.exists(os.path.join(dir, 'modelled_structures', PanddaDatasetFilenames.modelled_structure.format(dir_name))):
            log('No modelled structure in modelled_structures folder.')
            log('SKIPPING: {}'.format(dir))
            return

    ############################################################################
    # Export the pandda folder to output directory
    ############################################################################

    log.subheading('Exporting folder: {}'.format(dir))
    exp_dir = export_folder(dir=dir, params=params, log=log)

    ############################################################################
    # Merge input and pandda-modelled structures
    ############################################################################

    # Extract parameters for the merging and set them
    merging_params = merge_conformations.master_phil.extract()
    merging_params.input.major = os.path.join(exp_dir, params.output.file_prefix+PanddaDatasetFilenames.input_model.format(dir_name))
    merging_params.input.minor = os.path.join(exp_dir, params.output.file_prefix+PanddaDatasetFilenames.modelled_structure.format(dir_name))
    merging_params.output.pdb  = os.path.join(exp_dir, params.output.file_prefix+PanddaDatasetFilenames.ensemble_structure.format(dir_name))
    merging_params.output.log  = os.path.splitext(merging_params.output.pdb)[0]+'.log'
    merging_params.output.make_restraints = True
    # Apply settings
    merging_params.settings.overwrite = params.settings.overwrite
    merging_params.settings.verbose = params.settings.verbose
    # Change the restraints settings
    merging_params.restraints.output.phenix = os.path.splitext(os.path.basename(merging_params.output.pdb))[0]+'.restraints-phenix.params'
    merging_params.restraints.output.refmac = os.path.splitext(os.path.basename(merging_params.output.pdb))[0]+'.restraints-refmac.params'
    merging_params.restraints.output.log    = os.path.splitext(os.path.basename(merging_params.output.pdb))[0]+'.restraints.log'
    # Check files exist
    if not os.path.exists(merging_params.input.minor):
        raise Exception('File does not exist: {}'.format(merging_params.input.minor))
    if not os.path.exists(merging_params.input.major):
        raise Exception('File does not exist: {}'.format(merging_params.input.major))
    # Print and run
    log.subheading('Merging event-map model with input model')
    merge_conformations.run(params=merging_params)

############################################################################

def run(params):

    # Create log object
    log = Log(log_file=os.path.abspath(params.output.log_file), verbose=True)

    # Change paths to absolute paths
    params.input.pandda_dir  = os.path.abspath(params.input.pandda_dir)
    params.output.export_dir = os.path.abspath(params.output.export_dir)
    # Must be in the pandda directory (pandda objects use relative paths)
    os.chdir(params.input.pandda_dir)

    # Report modifed phil
    log.heading('Processed parameters')
    log(master_phil.format(params).as_str())

    ############################################################################

    log.heading('Identifying folders to export')

    # Find the dataset directories to be exported
    if params.input.select_datasets:
        selected_datasets = []; [selected_datasets.extend(s.split(',')) for s in params.input.select_datasets]
        export_dirs = sorted([os.path.join(params.input.pandda_dir, 'processed_datasets', p) for p in selected_datasets])
        # Filter by existence of path
        export_dirs = [p for p in export_dirs if os.path.exists(p)]
    else:
        export_dirs = sorted(glob.glob(os.path.join(params.input.pandda_dir, 'processed_datasets', '*')))
    assert export_dirs, 'No Export Directories Found'

    # Report
    log('Exporting:\n\t'+'\n\t'.join(export_dirs))

    # Create output directory
    if not os.path.exists(params.output.export_dir):
        os.mkdir(params.output.export_dir)

    # Merge the fitted structures
    for dir in export_dirs:
        process_and_export_folder(dir=dir, params=params, log=log)

    log.heading('FINISHED')

#######################################

if __name__ == '__main__':
    from giant.jiffies import run_default
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
        program             = PROGRAM,
        description         = DESCRIPTION)

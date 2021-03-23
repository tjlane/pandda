import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys, glob, shutil

from giant.exceptions import Sorry, Failure

from giant.jiffies import split_conformations

############################################################################

PROGRAM = 'giant.quick_refine'
DESCRIPTION = """
    A tool to simplify the launching of standard refinement jobs in REFMAC or PHENIX.

    1) Simple usage:
        > giant.quick_refine input.pdb input.mtz ligand.cif

    2) Defining custom names for the output folder and the output files
        > giant.quick_refine ... dir_prefix='XXX' out_prefix='XXX'

    3) Specify additonal parameter file (see giant.make_restraints)
        > giant.quick_refine ... params=restraints.params
"""

blank_arg_prepend = {
    '.mtz':'input.mtz=',
    '.pdb':'input.pdb=',
    '.cif':'input.cif=',
    '.params':'input.params=',
}

############################################################################

import libtbx.phil
master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .type = str
    mtz = 'refine.mtz'
        .type = str
    cif = None
        .type = str
        .multiple = True
    params = None
        .help = "File that contains additional parameters for the refinement program."
        .type = str
    args = None
        .help = "Pass any additional arguments to the program? (Command Line arguments). Pass as quoted strings."
        .type = str
        .multiple = True
}
output {
    dir_prefix = 'refine_'
        .type = str
    out_prefix = 'output'
        .type = str
    link_prefix = 'refine'
        .type = str
}
options {
    program = phenix *refmac
        .type = choice
    split_conformations = False
        .help = "Split the output structure into different conformations for modelling."
        .type = bool
}
split_conformations {
    include scope giant.jiffies.split_conformations.master_phil
}
settings {
    verbose = False
        .type = bool
}
""", process_includes=True)

############################################################################

def get_new_output_directory(dir_prefix):

    # Identify any existing output directories
    current_dirs = [d for d in sorted(glob.glob(dir_prefix+'*')) if os.path.isdir(d)]
    if not current_dirs:
        next_int = 1
    else:
        current_nums = [s.replace(dir_prefix, '') for s in current_dirs]
        next_int = sorted(map(int, current_nums))[-1]+1
    # Create output directory name from int
    out_dir = dir_prefix + '{:04}'.format(next_int)
    # Create and return directory
    os.mkdir(out_dir)
    return out_dir, current_dirs

def validate_params(params):

    logger.subheading('Validating input parameters')

    if (params.input.pdb is None):
        raise IOError('No PDB provided as input (input.pdb=...)')
    if (params.input.mtz is None):
        raise IOError('No MTZ provided as input (input.mtz=...)')

    if not os.path.exists(params.input.pdb):
        raise IOError('Input PDB file does not exist: {}'.format(params.input.pdb))
    if not os.path.exists(params.input.mtz):
        raise IOError('Input MTZ file does not exist: {}'.format(params.input.mtz))

    if os.path.islink(params.input.mtz):
        logger('Input mtz is symbolic link: resolving to real path.')
        real_path = os.path.realpath(params.input.mtz)
        logger('Updating input.mtz: \n\t{} -> \n\t{}'.format(params.input.mtz, real_path))
        params.input.mtz = real_path

def copy_input_files_to_output_folder(out_dir, params):

    logger.subheading('Copying/linking files to refinement folder')

    # Output filenames for "input" files
    input_pdb = os.path.abspath(os.path.join(out_dir, 'input.pdb'))
    input_mtz = os.path.abspath(os.path.join(out_dir, 'input.mtz'))
    input_params = os.path.abspath(os.path.join(out_dir, 'input.params'))
    input_args = os.path.abspath(os.path.join(out_dir, 'input.args'))

    # Copy input PDB file
    shutil.copy(params.input.pdb, input_pdb)

    # Only link the input mtz
    from giant.paths import rel_symlink
    rel_symlink(params.input.mtz, input_mtz)

    # Copy parameter file if it exists
    if params.input.params:
        shutil.copy(params.input.params, input_params)

    # Write additional arguments to file
    if params.input.args:
        with open(input_args, 'a') as fh:
            for a in params.input.args:
                fh.write(str(a)+'\n')

    logger('Input PDB file: {}'.format(input_pdb))
    logger('Input MTZ file: {}'.format(input_mtz))
    logger('Parameter file (if any): {}'.format(input_params))
    logger('Other arguments (if any): {}'.format(input_args))

def create_link(real_file, link_file):

    logger('Linking {} -> {}'.format(link_file, real_file))

    if os.path.islink(link_file):
        logger('Removing existing link: {}'.format(link_file))
        os.remove(link_file)

    from giant.paths import rel_symlink
    rel_symlink(real_file, link_file)

def process_output_files(params, refinement):

    # Find output files
    output_pdb = refinement.output_files['pdb']
    output_mtz = refinement.output_files['mtz']

    # Create links to main output files
    create_link(output_pdb, params.output.link_prefix+'.pdb')
    create_link(output_mtz, params.output.link_prefix+'.mtz')

    # Split conformations
    if params.options.split_conformations:
        logger.subheading('Splitting refined structure conformations')
        # Update params
        params.split_conformations.settings.verbose = params.settings.verbose
        # Running split conformations
        out_files = split_conformations.split_conformations(
            filename = output_pdb,
            params = params.split_conformations,
        )
        # Link output files to top
        for real_file in out_files:
            link_file = params.output.link_prefix + os.path.basename(real_file.replace(os.path.splitext(output_pdb)[0], ''))
            create_link(real_file, link_file)

############################################################################

def run(params):

    # Get next output directory
    out_dir, previous_dirs = get_new_output_directory(dir_prefix=params.output.dir_prefix)

    # Create log object
    log_file = os.path.join(out_dir, params.output.out_prefix+'.quick-refine.log')
    logger = lg.setup_logging(
        name = __name__,
        log_file = log_file,
    )

    # Report now that log_file is created
    if previous_dirs:
        logger('Found existing refinement directories: \n\t{}\n'.format('\n\t'.join(previous_dirs)))
    logger('Created new output directory: {}'.format(out_dir))

    # Validate input parameters
    validate_params(params)

    copy_input_files_to_output_folder(
        out_dir = out_dir,
        params = params,
    )

    # Create output prefixes
    output_prefix = os.path.join(out_dir, params.output.out_prefix)
    logger('Real output file path prefix: {}'.format(output_prefix))
    logger('Link output file path prefix: {}'.format(params.output.link_prefix))

    # Create command objects
    logger.subheading('Preparing command line input for refinement program')

    from giant.refinement.wrappers import get_refiner
    Refinement = get_refiner( params.options.program )

    refinement = Refinement(
        pdb_file = params.input.pdb,
        mtz_file = params.input.mtz,
        cif_files = params.input.cif,
        out_prefix = output_prefix,
    )

    # Parameter files
    if params.input.params and os.path.exists(params.input.params):
        refinement.add_parameter_file( params.input.params )

    # Pass additional command line arguments?
    if params.input.args:
        refinement.dispatcher.extend_args( params.input.args )

    # Report
    logger.subheading('Refinement parameters')
    logger( refinement.as_string() )

    logger.subheading('Running refinement ({})'.format(refinement.program))
    result = refinement.run()

    logger('Refinement output written to {}'.format(refinement.output_files['log']))

    logger.subheading('Post-processing output files')
    process_output_files(
        params = params,
        refinement = refinement,
    )

    logger.heading('finished normally! huzzah!')

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
        program             = PROGRAM,
        description         = DESCRIPTION,
    )

from __future__ import print_function

import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys, shutil

import numpy as np
import pathlib as pl

from pandda import (
    ModuleBanner,
    )

############################################################################

PROGRAM = 'pandda.export'

DESCRIPTION = """
Extract the structures and maps from a pandda analysis and prepare them for refinement
"""

program_banner = ModuleBanner(
    program = PROGRAM,
    description = DESCRIPTION,
    )

############################################################################

blank_arg_prepend = {
    None : 'select_datasets=',
}

import libtbx.phil
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
    overwrite = False
        .type = bool
    verbose = False
        .type = bool
    cpus = 1
        .type = int
}

""")

############################################################################

def insert_prefix(path, prefix):
    return path.parent / (prefix + path.name)

def get_one(iterator):
    iterator = list(iterator)
    if len(iterator) == 0: 
        raise ValueError('None found')
    elif len(iterator) > 1: 
        raise ValueError('More than one found')
    return iterator[0]

def check_dirs(directories):

    # Raise errors for non-existant
    missing = False

    for p in directories: 
        if not p.exists(): 
            logger.warning(
                'Directory does not exist: {}'.format(
                    str(p)
                    )
                )

    if missing is True: 
        raise ValueError('Some directories do not exist.')

def get_export_dirs(
    pandda_dir,
    selected_datasets = None,
    ):

    logger.heading('Identifying export folders')

    e_dir = pl.Path(pandda_dir) / 'processed_datasets'

    # Find the dataset directories to be exported
    if selected_datasets is not None:

        selected_datasets = sorted(selected_datasets)

        logger(
            'Selected folders: \n\t{}'.format(
                '\n\t'.join(selected_datasets)
                )
            )

        export_dirs = [
            (e_dir / d) for d in selected_datasets
            ]

        check_dirs(export_dirs)

    else:

        export_dirs = [
            p for p in e_dir.glob('*')
            if p.is_dir()
            ]

    if not export_dirs:
        logger('\nNo Export Directories Found!\n')
        sys.exit(1)

    return sorted(export_dirs)


class ExportFolder: 

    def __init__(self, 
        export_dir, 
        export_dir_prefix,
        export_file_prefix,
        overwrite = True,
        ):

        self.export_dir = pl.Path(export_dir)
        self.export_dir_prefix = str(export_dir_prefix)
        self.export_file_prefix = str(export_file_prefix)

        self.overwrite = bool(overwrite)

    def __call__(self, directory):

        i_dir = pl.Path(directory)

        e_dir = self.export_dir / (self.export_dir_prefix + i_dir.name)

        if not e_dir.exists(): 
            e_dir.mkdir(parents=True)

        logger(
            'Exporting\n\tfrom {inp}\n\tto{out}\n'.format(
                inp = str(i_dir),
                out = str(e_dir),
                )
            )
        
        export_list = self.get_export_list(directory=i_dir)


        logger(
            'Export list: \n\t{s}'.format(
                s = '\n\t'.join([
                    str(p[0].relative_to(i_dir)) 
                    for p in export_list
                    ]),
                )
            )

        # Export files
        for p_in, p_out_rel in export_list:

            # Check that the path exists
            if not p_in.exists():
                logger.warning('Input path does not exist: {!s}'.format(str(p_in)))
                continue

            p_out = (e_dir / p_out_rel)

            # Insert file prefix
            if p_out.is_file():
                p_out = insert_prefix(p_out, self.export_file_prefix)

            # Check to see if file already exists and delete if overwrite
            self.resolve_existing(path=p_out)

            self.copy(src_path=p_in, dst_path=p_out)

        return e_dir

    def get_export_list(self, directory):

        export_list = []

        for p in directory.glob('*.pdb'): 
            export_list.append(
                (p, p.name)
                )

        for p in directory.glob('*.mtz'): 
            export_list.append(
                (p, p.name)
                )

        for p in directory.glob('*.ccp4'): 
            export_list.append(
                (p, p.name)
                )

        for p in (directory/'modelled_structures').glob('*modelled.pdb'):
            export_list.append(
                (p, p.name)
                )

        for p in (directory / 'modelled_structures').glob('*ensemble.pdb'):
            export_list.append(
                (p, p.name)
                )

        export_list.append(
            (directory/'ligand_files', 'ligand_files')
            )

        return export_list

    def resolve_existing(self, path):

        if path.exists():
            if self.overwrite is True:
                if path.is_dir():
                    shutil.rmtree(str(path))
                else:
                    os.remove(str(path))
            else:
                raise Exception('Output path already exists: {}. Need to set overwrite=True to continue.'.format(export_file))

    def copy(self, src_path, dst_path):

        assert src_path.exists()

        if src_path.is_file():
            shutil.copytree(str(src_path), str(dst_path))
        else: 
            shutil.copyfile(str(src_path), str(dst_path))


class PostProcessFolder: 

    def __init__(self,
        major_state_glob = "*-pandda-input.pdb",
        minor_state_glob = "*-pandda-model.pdb",
        merged_template = "{dir_name}-ensemble-model.pdb",
        make_restraints = True,
        overwrite = True,
        ):
        
        self.major_state_glob = (
            major_state_glob
            )

        self.minor_state_glob = (
            minor_state_glob
            )

        self.merged_template = (
            merged_template
            )

        self.merging_params = self.get_merging_params(
            make_restraints = make_restraints,
            overwrite = overwrite,
            )

    def __call__(
        dir_path,
        ):

        dir_name = dir_path.name

        try: 

            major_state_path = get_one(dir_path.glob(self.major_state_glob))
            minor_state_path = get_one(dir_path.glob(self.minor_state_glob))

        except ValueError as e: 

            logger(
                "Looking for\n\t{d}\n\tand {m1}\n\tin {m2}".format(
                    d = str(dir_path),
                    m1 = major_state_glob,
                    m2 = minor_state_glob,
                    )
                )

            raise IOError(
                "Error finding files in output folder: {}".format(
                    str(e)
                    )
                )

        merged_path = self.merged_template.format(
            dir_name = dir_name,
            )

        assert major_state_path.exists()
        assert minor_state_path.exists()

        self.run_merging(
            major_state_path = major_state_path,
            minor_state_path = minor_state_path,
            merged_path = merged_path,
            prefix = dir_name,
            )

    def get_merging_params(self,
        make_restraints,
        overwrite = True,
        ):

        from giant.jiffies import merge_conformations as mc

        # Extract parameters for the merging and set them
        merging_params = mc.master_phil.extract()
        merging_params.output.make_restraints = bool(make_restraints)
        merging_params.settings.overwrite = overwrite

        return merging_params

    def run_merging(self, 
        major_state_path, 
        minor_state_path,
        merged_path,
        ):

        m_params = copy.deepcopy(self.merging_params)

        d_name = major_state_path.parent.name

        m_params.input.major = major_state_path
        m_params.input.minor = minor_state_path

        m_params.output.pdb = merged_path
        m_params.output.log = merged_path.with_suffix('.log')

        m_params.restraints.output.phenix = merged_path.with_suffix('.restraints-phenix.params')
        m_params.restraints.output.refmac = merged_path.with_suffix('.restraints-refmac.params')
        m_params.restraints.output.log = merged_path.with_suffix('.restraints.log')
        
        if not os.path.exists(m_params.input.minor):
            raise IOError('Input file does not exist: {}'.format(m_params.input.minor))

        if not os.path.exists(m_params.input.major):
            raise IOError('Input file does not exist: {}'.format(m_params.input.major))

        logger.subheading('Merging event-map model with input model')
        mc.run(params=m_params)


class ValidateFolder: 

    def __init__(self, 
        required_file_for_export,
        required_model_glob = "modelled_structures/*-pandda-model.pdb",
        ):

        self.required_file_for_export = required_file_for_export
        self.model_template_glob = model_template_glob

    def __call__(self,
        dir_path,
        ):

        if self.required_file_for_export == 'model':

            if list(dir_path.glob(self.required_model_glob)):
                return True
            else: 
                return False

        return True

def standard_pandda_export(params):

    export_dirs = get_export_dirs(
        pandda_dir = params.input.pandda_dir,
        selected_datasets = (
            np.concatenate(
                [s.split(',') for s in params.input.select_datasets]
                ).tolist(),
            if (params.input.select_datasets is not None)
            else None
            )
        )

    # Report
    logger(
        'Exporting:\n\t{}'.format(
            '\n\t'.join(map(str, export_dirs))
            )
        )

    # Create output directory
    if not os.path.exists(params.output.export_dir):
        os.mkdir(params.output.export_dir)

    validate_folder = ValidateFolder(
        required_file_for_export = params.output.required_file_for_export,
        )
    
    export_folder = ExportFolder(
        export_dir = params.output.export_dir,
        export_dir_prefix = params.output.dir_prefix,
        export_file_prefix = params.output.file_prefix,
        overwrite = params.settings.overwrite,
        )

    postprocess_folder = PostProcessFolder(
        make_restraints = params.options.generate_restraints,
        overwrite = params.settings.overwrite,
        )

    # Merge the fitted structures
    for d_in in export_dirs:

        logger.heading('Processing directory: {}'.format(d_in.name), spacer=True)

        if not validate_folder(d_in): 
            continue

        d_out = export_folder(d_in)

        postprocess_folder(d_out)

    return

############################################################################

def run(args):

    from giant import jiffies
    params = jiffies.extract_params_default(
        master_phil = master_phil,
        args = args,
        blank_arg_prepend = blank_arg_prepend,
        home_scope = None,
    ).extract()

    logger = lg.setup_logging(
        name = __name__,
        log_file = params.output.log_file,
        warning_handler_name = 'warnings',
        debug = False,
        )

    # Report modifed phil
    logger.heading('Processed parameters')
    logger(master_phil.format(params).as_str())

    ############################################################################

    standard_pandda_export(params)

    logger.heading('FINISHED')

#######################################

def main():

    try: 

        lg.setup_logging_basic(__name__)

        run(args=sys.argv[1:])

    except Exception as e: 

        import traceback
        logger.subheading("Error - stack trace")
        logger(traceback.format_exc())#limit=10))
        logger.subheading("PanDDA exited with an error")
        logger("Error: %s", str(e))
        logger.subheading("PanDDA exited with an error")
        sys.exit(1)

if __name__ == '__main__':
    
    main()


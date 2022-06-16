import giant.logs as lg
logger = lg.getLogger(__name__)

import sys
import pathlib as pl

from giant.phil import (
    log_running_parameters,
    )

import matplotlib as mpl

############################################################################

PROGRAM = 'giant.mass_refine'

DESCRIPTION = """
    A tool to refine a single model against multiple datasets.
"""

############################################################################

blank_arg_prepend = {
    'None' : 'mode=',
    '.mtz' : 'import_data.mtz=',
    '.pdb' : 'refine_model.pdb=',
    '.cif' : 'refine_model.cif_restraints=',
    '.params' : 'refine_model.refinement_parameters='
}

import libtbx.phil
master_phil = libtbx.phil.parse("""
mode = None
    .type = str
work_folder = mass_refine
    .type = path
log = None
    .type = path
configure {
    backend = serial *parallel_joblib
        .type = choice(multi=False)
}
import_data {
    mtz = None
        .type = path
        .multiple = True
    labelling = *filename foldername none
        .type = choice(multi=False)
    meta_mtz_regex = None
        .help = "GI_(?P<temp>[0-9]*)C_APO"
        .type = str
    meta_csv_path = None
        .type = str
    rename_using_meta = True
        .type = bool
}
refine_model {
    pdb = None
        .type = str
        .multiple = False
    labelling = *filename foldername none
        .type = choice(multi=False)
    refinement_parameters = None
        .type = str
        .multiple = True
    cif_restraints = None
        .type = str
        .multiple = True
}
settings {
    verbose = False
        .type = bool
}
""", process_includes=True)

############################################################################

def set_get_log(params):

    if params.log is not None: 
        pass
    elif params.work_folder is not None: 
        params.log = str(
            pl.Path(params.work_folder).with_suffix('.log')
            )

    return params.log

def validate_params(params):

    if params.mode is None: 
        raise Exception('No mode specified (mode=X)')
    elif '.' in params.mode: 
        raise Exception('Invalid model: {}'.format(mode))


class MassRefineProtocol(object):

    update_config = None
    import_data = None
    refine_model = None

    def __init__(self, params):

        self.run_setup(params)

        if params.mode == "configure": 

            self.update_config = self.get_update_config(
                params = params.configure,
                )

        elif params.mode == "import": 

            self.import_data = self.get_import_data(
                params = params.import_data,
                )

            self.import_data_args = {
                'data' : list(params.import_data.mtz),
                }
        
        elif params.mode == "refine": 

            self.import_model  = self.get_import_model(
                params = params.refine_model,
                )

            self.refine_model = self.get_refine_model(
                params = params.refine_model,
                )

            self.import_model_args = {
                'model' : (
                    pl.Path(params.refine_model.pdb)
                    ),
                'cif_restraints' : (
                    params.refine_model.cif_restraints
                    ),
                'refinement_parameters' : (
                    params.refine_model.refinement_parameters
                    ),
                }

        else: 
            raise Exception(
                'mode not recognised: {}'.format(params.mode)
                )

    def run_setup(self, params):

        from giant.mass_refine import (
            MassRefineFileSystem,
            MassRefineConfiguration,
            )

        self.file_system = MassRefineFileSystem(
            str(params.work_folder),
            )

        self.config = MassRefineConfiguration(
            json_path = self.file_system.config_path,
            )

    def run(self):

        self.file_system.update()

        if self.update_config is not None:

            logger.heading('Configuring setup')

            self.update_config(
                self.config,
                )

        if self.import_data is not None: 

            logger.heading('Importing data')
            
            new_datasets = self.import_data.run(
                data_paths = (
                    self.import_data_args['data']
                    ),
                )

            logger(str(new_datasets))

        if self.refine_model is not None: 
            
            # Import the model into the system
            model_folder = self.import_model.run(
                model_path = (
                    self.import_model_args['model']
                    ),
                refinement_parameters_paths = (
                    self.import_model_args['refinement_parameters']
                    ),
                cif_restraints_paths = (
                    self.import_model_args['cif_restraints']
                    ),
                )

            self.refine_model.run(
                model_folder = model_folder,
                )

        return self

    def get_update_config(self,
        params,
        ):

        p = params

        from giant.mass_refine.processing import (
            UpdateConfiguration,
            )

        update_config = (
            UpdateConfiguration(
                backend = p.backend,
                )
            )

        return update_config

    def get_import_data(self, 
        params,
        ): 

        p = params
        fs = self.file_system

        from giant.mass_refine.import_data import (
            ImportDatasets,
            GetMetaFromTable,
            GetMetaFromRegex,
            )

        from giant.mulch.labelling import (
            PathLabeller,
            )

        get_meta = (
            GetMetaFromTable(
                p.meta_csv_path
                )
            if (p.meta_csv_path is not None)
            else 
            GetMetaFromRegex(
                p.meta_mtz_regex
                )
            if (p.meta_mtz_regex is not None)
            else None
            )

        import_data = ( 
            ImportDatasets(
                file_system = fs,
                get_meta = get_meta,
                get_label = PathLabeller(
                    method = p.labelling,
                    ),
                use_meta_as_label = bool(
                    p.rename_using_meta, 
                    ),
                )
            )

        return import_data

    def get_import_model(self,
        params,
        ):

        p = params
        fs = self.file_system

        from giant.mass_refine.import_model import (
            ImportModel,
            )

        from giant.mulch.labelling import (
            PathLabeller,
            )

        import_model = (
            ImportModel(
                file_system = fs,
                get_label = PathLabeller(
                    method = p.labelling,
                    ),
                )
            )

        return import_model

    def get_refine_model(self,
        params,
        ): 

        p = params
        fs = self.file_system
        proc = self.config.processor

        from giant.mass_refine.refine_model import (
            RefineModel,
            )
        
        from giant.refinement.wrappers import (
            get_refiner,
            )

        refine_model = (
            RefineModel(
                file_system = fs,
                refiner = get_refiner("refmac"),
                processor = proc,
                )
            )

        return refine_model


def run(params):

    logger = lg.setup_logging(
        name = __name__,
        log_file = set_get_log(params),
        debug = bool(params.settings.verbose),
        )

    # Report
    logger.subheading(
        'Validating input parameters and input files'
        )

    validate_params(params)

    log_running_parameters(
        params = params,
        master_phil = master_phil,
        logger = logger,
    )

    logger.subheading(
        'Setting up...'
        )

    mass_refine = MassRefineProtocol(params)

    logger(str(mass_refine))

    logger.heading('running')

    mass_refine.run()

    logger.heading('mass_refine finished normally')

############################################################################

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

import giant.logs as lg
logger = lg.getLogger(__name__)


class ImportModel(object):

    def __init__(self, 
        file_system,
        get_label = None,
        get_meta = None,
        ):

        self.file_system = (
            file_system
            )

        self.get_label = (
            get_label
            )

    def run(self, 
        model_path, 
        refinement_parameters_paths = None,
        cif_restraints_paths = None,
        ):

        self.file_system.update()

        model_folder = self.file_system.models.import_model(
            label = self.get_label(model_path),
            model_path = model_path,
            refinement_parameters_paths = refinement_parameters_paths,
            cif_restraints_paths = cif_restraints_paths,
            )

        return model_folder

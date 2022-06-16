import giant.logs as lg
logger = lg.getLogger(__name__)

import shutil 


class RefinementGenerator(object):

    def __init__(self, 
        refiner,
        ):

        self.refiner = refiner

    def __call__(self,
        input_model,
        input_data,
        cif_files,
        parameter_files,
        output_prefix,
        ):

        refine = self.refiner(
            pdb_file = str(input_model),
            mtz_file = str(input_data),
            cif_files = (
                list(map(str,cif_files))
                if cif_files is not None
                else None
                ),
            out_prefix = str(output_prefix),
            )

        if (parameter_files is not None):

            for p in parameter_files: 
                
                refine.add_parameter_file(
                    str(p),
                    )

        return refine


class RefineSingle(object):

    def __init__(self, 
        refinement_folder,
        refinement_generator,
        ):

        # Object to be returned
        self.refinement_folder = (
            refinement_folder
            )

        # Make temporary output folder 
        self.tmp_dir = (
            self.refinement_folder.tmp_dir()
            )

        input_model = self.refinement_folder.input_model
        input_data = self.refinement_folder.input_data

        self.refinement = refinement_generator(
            input_model = input_model.model(), # must exist
            input_data = input_data.data(), # must exist
            cif_files = input_model.cif_restraints(),
            parameter_files = input_model.refinement_parameters(),
            output_prefix = str(self.tmp_dir / "refine"),
            )

    def __call__(self):

        result = self.refinement.run()

        self.process_result(result)

        self.copy_output_files()

        return self.refinement_folder

    def process_result(self, result):
        
        pass

    def copy_output_files(self):

        shutil.copy(
            str(self.refinement.output_files['pdb']),
            str(self.refinement_folder.model_path),
            )

        shutil.copy(
            str(self.refinement.output_files['mtz']),
            str(self.refinement_folder.data_path),
            )


class RefineModel(object):

    def __init__(self, 
        file_system, 
        refiner, 
        processor,
        ):

        self.file_system = file_system

        self.processor = processor

        self.refinement_generator = RefinementGenerator(
            refiner = refiner,
            )

    def run(self, model_folder):

        refine_folders = self.get_refinement_folders(
            model_folder = model_folder,
            )

        self.run_refinements(
            refine_folders = refine_folders,
            )

        return refine_folders

    def get_refinement_folders(self, model_folder):

        model_label = model_folder.get_meta()['label']

        # Create output folders by dataset
        #
        datasets = (
            self.file_system.datasets.get_all_as_dict()
            )
        #
        refine_folders = {}
        #
        for dkey, d_folder in list(datasets.items()):

            r_folder = d_folder.refinements.create_refinement_folder(
                label = model_label,
                model_folder = model_folder,
                meta_dict = None, # TODO
                )

            refine_folders[dkey] = r_folder

        return refine_folders

    def run_refinements(self, refine_folders):

        refinements = {
            dkey : RefineSingle(
                refinement_folder = r_folder,
                refinement_generator = self.refinement_generator,
                )
            for dkey, r_folder in refine_folders.items()
        }

        for dkey, drefine in sorted(refinements.items()):
            logger(str(drefine.refinement)[:2000])

        proc = self.processor.as_dict_processor()

        return proc(refinements)


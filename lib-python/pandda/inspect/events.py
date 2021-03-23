import os 
import pathlib as pl
from giant.paths import rel_symlink


class Event:

    def __init__(self, info):

        self.index = info.name

        self.dtag = info.name[0]
        self.event_num = int(info.name[1])

        self.map_resolution = round(info['analysed_resolution'], 2)
        self.map_uncertainty = round(info['map_uncertainty'], 2)

        self.rwork_rfree = (round(info['r_work'], 3), round(info['r_free'], 3))

        self.site_num = int(info['site_num'])

        self.est_bdc = round(info['bdc'], 2)
        self.est_1_bdc = round(1-info['bdc'], 2)

        self.z_peak = round(info['z_peak'], 1)
        self.z_mean = info['z_mean']
        self.cluster_size = int(info['cluster_size'])

        self.xyz = info[['x', 'y', 'z']]

        self.added_ligand_names = [
            s for s in info['Ligand Names'].split(',') if s
            ]


class GetEventFiles: 

    def __init__(self,
        pandda_directory,
        pandda_files_dict, 
        pandda_path_prefix,
        ):

        self.pandda_directory = pandda_directory
        self.pandda_files_dict = pandda_files_dict
        self.pandda_path_prefix = pandda_path_prefix

    def __call__(self,
        event,
        ):

        files_dict = self.pandda_files_dict

        # Indentify main dataset directory 

        dataset_dir = (
            self.pandda_directory / 'processed_datasets' / event.dtag
            )

        assert dataset_dir.exists(), 'Directory does not exist: {}'.format(event.dtag)

        # Identify dataset subdirectories and files

        ligand_dir = (
            dataset_dir / 'ligand_files'
            )

        if not ligand_dir.exists(): 
            ligand_dir.mkdir(parents=True)

        model_dir = (
            dataset_dir / 'modelled_structures'
            )

        if not model_dir.exists(): 
            model_dir.mkdir(parents=True)

        # The most recent model of the protein in the pandda maps

        output_model_link = (
            model_dir / '{dkey}-pandda-model.pdb'.format(dkey=event.dtag)
            )

        # Files output by pandda

        dataset_files = files_dict['dataset_files'][event.dtag]

        import json
        print(json.dumps(dataset_files, indent=2))

        input_model = (
            self.pandda_directory / pl.Path(
                dataset_files['structure']
                ).relative_to(self.pandda_path_prefix)
            )

        input_data = (
            self.pandda_directory / pl.Path(
                dataset_files['data']
                ).relative_to(self.pandda_path_prefix)
            )

        output_data = (
            self.pandda_directory / pl.Path(
                dataset_files['output_data']
                ).relative_to(self.pandda_path_prefix)
            )

        event_data = (
            self.pandda_directory / pl.Path(
                dataset_files['event_data'][str(event.event_num)]
                ).relative_to(self.pandda_path_prefix)
            )

        return {
            'dataset_dir' : dataset_dir,
            'ligand_dir' : ligand_dir,
            'output_model_dir' : model_dir,
            'output_model_link' : str(output_model_link),
            'input_model' : str(input_model),
            'input_data' : str(input_data),
            'output_data' : str(output_data),
            'event_data' : str(event_data),
            # 'ligand_pdbs' : map(str, ligand_pdbs),
            # 'ligand_cifs' : map(str, ligand_cifs),
        }


class GetNextModelFile:

    output_prefix = 'fitted-v'

    def __init__(self,
        model_directory,
        output_model_link, 
        ):

        self.model_directory = model_directory
        self.output_model_link = str(output_model_link)

    def __call__(self):

        return self.next_file_and_update()

    def next_file_and_update(self):

        next_path = self.next_file()

        self.update_link(next_path)

        return next_path

    def last_file(self):
        """Get the most recent saved model of this protein"""

        fitted_outputs = map(
            str,
            sorted(
                self.model_directory.glob(self.output_prefix+'*')
                )
            )

        if fitted_outputs:
            print 'Current models: \n\t{}'.format('\n\t'.join(fitted_outputs))
            return fitted_outputs[-1]
        else:
            print 'No current models'
            return None

    def next_file(self):

        current = self.last_file()

        print 'Most recent saved model: {!s}'.format(current)

        if current:
            last_idx = int(current[-8:-4]) # matches to {:04d} below
        else:
            last_idx = 0 # will be incremented

        new_fitted = self.output_prefix + '{:04d}.pdb'.format(last_idx + 1)

        return str(self.model_directory / new_fitted)

    def update_link(self, path=None):
        
        if path is None: 
            path = self.last_file()

        # Always want to remove the link
        if os.path.islink(self.output_model_link):
            os.remove(self.output_model_link)

        # Check if it's a file just in case...
        if os.path.exists(self.output_model_link):
            os.remove(self.output_model_link)

        # No files
        if path is None: 
            return

        print 'Linking {!s} -> {!s}'.format(
            os.path.basename(path), 
            os.path.basename(self.output_model_link),
            )

        # Create new link the most recent file
        from giant.paths import rel_symlink
        rel_symlink(
            path, 
            self.output_model_link,
            )

        return self.output_model_link

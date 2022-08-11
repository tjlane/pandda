import giant.logs as lg
logger = lg.getLogger(__name__)

import os

import pathlib as pl
import numpy as np

import gemmi

from pandda.utils import merge_dicts


class ConvertMapsToMtz(object):

    def __init__(self,
        dataset_name = 'pandda',
        delete_maps = False,
        ):

        self.dataset_name = dataset_name
        self.delete_maps = delete_maps

    def __call__(self,
        map_paths,
        map_labels,
        filepath,
        d_min,
        ):

        mtz_obj = self.get_mtz_obj()

        all_mtz_data = []
        first_map = True

        for m_path, m_label in zip(map_paths, map_labels):

            m_sfs = self.get_structure_factors_from_map(
                map_path = m_path,
                d_min = d_min,
                )

            mtz_data = self.prepare_mtz_data(
                mtz_obj = mtz_obj,
                map_sfs = m_sfs,
                map_label = m_label,
                add_base = first_map,
                )

            # Add columns from this map
            all_mtz_data.extend(mtz_data)

            # Prevent adding of further miller arrays
            first_map = False

        mtz_obj.set_data(
            np.concatenate(all_mtz_data, axis=1)
            )

        mtz_obj.write_to_file(
            str(filepath),
            )

        if self.delete_maps is True:
            self.remove_files(map_paths)

        return str(filepath)

    def get_mtz_obj(self):

        mtz = gemmi.Mtz(with_base=True)
        
        mtz.add_dataset(self.dataset_name)

        return mtz

    def get_structure_factors_from_map(self, map_path, d_min):

        map_obj = gemmi.read_ccp4_map(map_path)

        # Extract original spacegroup, and swap for P1 spacegroup
        orig_sg = map_obj.grid.spacegroup
        map_obj.grid.spacegroup = gemmi.SpaceGroup('P1') # this is a hack

        # Now do setup
        map_obj.setup(default_value=0.0)
        
        # Now extract the structure factors and reapply spacegroup
        sf = gemmi.transform_map_to_f_phi(
            map_obj.grid,
            half_l = True,
            )
        sf.spacegroup = orig_sg # end of hack

        # Extract the ASU
        asu_sf = sf.prepare_asu_data(
            dmin = d_min,
            with_000 = True,
            )

        return asu_sf

    def prepare_mtz_data(self,
        mtz_obj,
        map_sfs,
        map_label,
        add_base = False,
        ):

        #####

        mtz_obj.add_column(
            'F{}'.format(map_label),
            'F'
            )
        mtz_obj.add_column(
            'PH{}'.format(map_label),
            'P',
            )

        #####

        out_data = []

        if add_base is True:

            mtz_obj.spacegroup = map_sfs.spacegroup

            mtz_obj.set_cell_for_all(map_sfs.unit_cell)

            out_data.append(map_sfs.miller_array)

        out_data.append(
            np.absolute(
                map_sfs.value_array
                ).reshape(
                (map_sfs.value_array.size, 1)
                )
            )
        out_data.append(
            np.angle(
                map_sfs.value_array, deg=True,
                ).reshape(
                (map_sfs.value_array.size, 1)
                )
            )

        return out_data

    def remove_files(self, filepaths):

        import os
        for f in filepaths:
            os.remove(str(f))


class MakePanddaMainMtz(object):

    output_key = "output_data"
    output_filename_template = "{dkey}-pandda-output.mtz"

    def __init__(self,
        convert_maps_to_mtz,
        output_directory,
        ):

        self.convert_maps_to_mtz = convert_maps_to_mtz
        self.output_directory = output_directory

    def __call__(self,
        dataset_key, 
        dataset_files,
        d_min,
        ):

        filepath = (
            self.output_directory / self.output_filename_template.format(
                dkey = dataset_key,
                )
            )

        map_paths, map_labels = self.get_map_paths_and_labels(
            dataset_files = dataset_files,
            )

        if len(map_paths) == 0:
            return None

        self.convert_maps_to_mtz(
            map_paths = map_paths,
            map_labels = map_labels,
            filepath = str(filepath),
            d_min = d_min,
            )

        return {self.output_key : str(filepath)}

    def get_map_paths_and_labels(self, dataset_files):

        map_paths = []
        map_labels = []

        if ('dataset_map' in dataset_files):

            map_paths.append(
                dataset_files['dataset_map']
                )
            map_labels.append(
                'DATASET'
                )

        if ('ground_map' in dataset_files):

            map_paths.append(
                dataset_files['ground_map']
                )
            map_labels.append(
                'GROUND'
                )

        if ('z_map' in dataset_files):

            map_paths.append(
                dataset_files['z_map']
                )
            map_labels.append(
                'ZVALUES'
                )

        return (map_paths, map_labels)


class MakePanddaEventMtz(object):

    output_key = "event_data"
    output_filename_template = "{dkey}-pandda-output-event-{event_num:03d}.mtz"

    def __init__(self,
        convert_maps_to_mtz,
        output_directory,
        ):

        self.convert_maps_to_mtz = convert_maps_to_mtz
        self.output_directory = output_directory

    def __call__(self,
        event_num,
        event_bdc,
        dataset_key, 
        dataset_files,
        d_min,
        ):

        filepath = (
            self.output_directory / self.output_filename_template.format(
                dkey = dataset_key,
                event_num = event_num,
                )
            )

        map_filepath = filepath.with_suffix(
            '.tmp{event_num}.ccp4'.format(event_num=event_num)
            )
        assert not map_filepath.exists()

        self.make_event_map(
            event_bdc = event_bdc,
            dataset_files = dataset_files,
            filepath = str(map_filepath),
            )

        self.convert_maps_to_mtz(
            map_paths = [str(map_filepath)],
            map_labels = ['EVENT'],
            filepath = str(filepath),
            d_min = d_min,
            )

        os.remove(
            str(map_filepath)
            )

        return {self.output_key : {event_num : str(filepath)}}

    def make_event_map(self,
        event_bdc,
        dataset_files,
        filepath,
        ):

        ground_map = gemmi.read_ccp4_map(
            dataset_files['ground_map'],
            )

        dataset_map = gemmi.read_ccp4_map(
            dataset_files['dataset_map'],
            )

        event_map_values = (
            np.array(dataset_map.grid) - (event_bdc * np.array(ground_map.grid))
            )

        event_map_grid = gemmi.FloatGrid(
            event_map_values,
            cell = dataset_map.grid.unit_cell,
            spacegroup = dataset_map.grid.spacegroup,
            )

        # This is now the event map
        dataset_map.grid = event_map_grid

        dataset_map.write_ccp4_map(
            str(filepath)
            )
        
###


class MakePanddaEvaluationMtzs(object):

    output_key = "dataset_files"

    def __init__(self,
        make_maps,
        output_dir,
        dataset_dir = None,
        dataset_subdir = "",
        processor = None,
        delete_maps = True,
        ):

        if (processor is None):
            from giant.processors import basic_processor
            processor = basic_processor

        if (dataset_dir is None):
            dataset_dir = output_dir

        self.make_maps = make_maps

        self.make_mtz = ConvertMapsToMtz(
            delete_maps = False, # Done at the end
            )

        self.output_dir = pl.Path(output_dir)
        self.dataset_dir = dataset_dir
        self.dataset_subdir = dataset_subdir

        self.delete_maps = delete_maps

        self.processor = processor

    def __call__(self,
        datasets,
        datasets_map_dict,
        datasets_results,
        statistical_model,
        map_resolution,
        *args, **kwargs
        ):

        map_files_dict = self.make_maps(
            datasets = datasets,
            datasets_map_dict = datasets_map_dict,
            datasets_results = datasets_results,
            statistical_model = statistical_model,
            ).get(
            self.make_maps.output_key, None
            )

        if map_files_dict is None: 
            raise Exception('Failed to make maps.')

        ###

        logger('Writing mtzs...')

        job_list = []
        job_keys = []

        for dkey, dataset_files in map_files_dict.items():

            d_jobs = self.get_jobs(
                dataset_key = dkey,
                dataset_files = dataset_files,
                dataset_events = datasets_results[dkey],
                map_resolution = map_resolution,
                )

            if d_jobs is None:
                continue
            
            for j in d_jobs: 
                job_list.append(j)
                job_keys.append(dkey)

        output_files = {}

        if len(job_list) > 0:

            results = self.processor(job_list)

            for dkey, of in zip(job_keys, results):

                if of is None: 
                    continue

                merge_dicts(
                    master_dict = output_files, 
                    merge_dict = {dkey : of}, # every sub-dict "of" is required to have a unique structure
                    )

        if (self.delete_maps is True):
            self.delete_files(
                file_dict = map_files_dict,
                )

        return {self.output_key : output_files}

    def get_jobs(self,
        dataset_key,
        dataset_files,
        dataset_events,
        map_resolution,
        ):

        main_job = self.get_main_job(
            dataset_key = dataset_key,
            dataset_files = dataset_files,
            map_resolution = map_resolution,
            )

        # No maps to create
        if main_job is None: 
            return None

        event_jobs = self.get_event_jobs(
            dataset_key = dataset_key,
            dataset_files = dataset_files,
            dataset_event_bdcs = {e['event_num'] : e['bdc'] for e in dataset_events},
            map_resolution = map_resolution,
            )

        job_list = [main_job]
        job_list.extend(event_jobs)

        return job_list

    def get_main_job(self,
        dataset_key,
        dataset_files,
        map_resolution,
        ):

        make_main_mtz = MakePanddaMainMtz(
            convert_maps_to_mtz = self.make_mtz,
            output_directory = (
                self.dataset_dir / dataset_key / self.dataset_subdir
                ),
            )

        job = self.get_wrapper(
            function = make_main_mtz,
            dataset_key = dataset_key,
            dataset_files = dataset_files,
            d_min = map_resolution,
            )

        return job

    def get_event_jobs(self,
        dataset_key,
        dataset_files,
        dataset_event_bdcs,
        map_resolution,
        ):

        make_event_mtz = MakePanddaEventMtz(
            convert_maps_to_mtz = self.make_mtz,
            output_directory = (
                self.dataset_dir / dataset_key / self.dataset_subdir
                ),
            )

        jobs = [
            self.get_wrapper(
                function = make_event_mtz,
                event_num = event_num,
                event_bdc = event_bdc,
                dataset_key = dataset_key, 
                dataset_files = dataset_files,
                d_min = map_resolution,
                )
            for event_num, event_bdc 
            in dataset_event_bdcs.items()
            ]

        return jobs

    def get_wrapper(self,
        function,
        *args, 
        **kwargs
        ):

        from giant.processors import ProcessWrapper
        return ProcessWrapper(
            func = function,
            *args,
            **kwargs
            )

    def delete_files(self,
        file_dict,
        ):

        for dkey, d_files in file_dict.items():

            for m_key, m_file in d_files.items():

                if str(m_file).endswith('.ccp4'): 

                    os.remove(str(m_file))


class GetMakePanddaEvaluationMtzs(object):

    def __init__(self,
        get_dataset_map_writer,
        output_dir,
        dataset_dir = None,
        dataset_subdir = "",
        processor = None,
        output_requires_events = True,
        ):

        from .maps import (
            GetMakePanddaEvaluationMaps,
            )

        self.get_make_maps = GetMakePanddaEvaluationMaps(
            get_dataset_map_writer = get_dataset_map_writer,
            output_dir = output_dir,
            dataset_dir = dataset_dir,
            dataset_subdir = dataset_subdir,
            processor = processor,
            output_requires_events = output_requires_events,
            )

        self.output_dir = output_dir
        self.dataset_dir = dataset_dir
        self.dataset_subdir = dataset_subdir

        self.processor = processor

    def __call__(self, label=""):

        return MakePanddaEvaluationMtzs(
            make_maps = self.get_make_maps(
                label = label,
                ),
            output_dir = (
                self.output_dir / label
                ),
            dataset_dir = (
                self.dataset_dir
                ),
            dataset_subdir = (
                self.dataset_subdir
                ),
            processor = self.processor,
            )

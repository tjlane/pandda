import giant.logs as lg
logger = lg.getLogger(__name__)

import pathlib as pl

import gemmi
import numpy as np


class ConvertMapsToMtz:

    output_key = "output_maps"

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

            m_grid = self.load_map(
                map_path = m_path,
                )

            mtz_data = self.prepare_mtz_data(
                mtz_obj = mtz_obj,
                map_grid = m_grid,
                map_label = m_label,
                d_min = d_min,
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

        return {self.output_key : str(filepath)}

    def get_mtz_obj(self):

        mtz = gemmi.Mtz(with_base=True)
        mtz.add_dataset(self.dataset_name)

        return mtz

    def load_map(self, map_path):

        map_obj = gemmi.read_ccp4_map(map_path)
        map_obj.setup()

        return map_obj.grid

    def prepare_mtz_data(self,
        mtz_obj,
        map_grid,
        map_label,
        d_min,
        add_base = False,
        ):

        #####

        sf = gemmi.transform_map_to_f_phi(
            map_grid,
            half_l = True,
            )

        asu_sf = sf.prepare_asu_data(
            dmin = d_min,
            with_000 = True,
            )

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

            mtz_obj.spacegroup = sf.spacegroup

            mtz_obj.set_cell_for_all(sf.unit_cell)

            out_data.append(asu_sf.miller_array)

        out_data.append(
            np.absolute(
                asu_sf.value_array
                ).reshape(
                (asu_sf.value_array.size, 1)
                )
            )
        out_data.append(
            np.angle(
                asu_sf.value_array, deg=True,
                ).reshape(
                (asu_sf.value_array.size, 1)
                )
            )

        return out_data

    def remove_files(self, filepaths):

        import os
        for f in filepaths:
            os.remove(str(f))


###


class MakePanddaEvaluationMtzs:

    output_key = "dataset_files"
    output_filename_template = "{dkey}-pandda-output.mtz"

    def __init__(self,
        make_maps,
        output_dir,
        dataset_dir = None,
        dataset_subdir = "",
        # processor = None,
        ):

        # Currently set to serial -- change later?
        processor = None

        if (processor is None):
            from giant.processors import basic_processor
            processor = basic_processor

        if (dataset_dir is None):
            dataset_dir = output_dir

        self.make_maps = make_maps

        self.make_mtz = ConvertMapsToMtz(
            delete_maps = False,
            )

        self.output_dir = pl.Path(output_dir)
        self.dataset_dir = dataset_dir
        self.dataset_subdir = dataset_subdir

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

        ###

        logger('Writing mtzs...')

        jobs = []
        keys = []

        for dkey, dataset_files in map_files_dict.items():

            j = self.get_job(
                dataset_key = dkey,
                dataset_files = dataset_files,
                map_resolution = map_resolution,
                )

            if j is not None:
                jobs.append(j)
                keys.append(dkey)

        if jobs:
            results = self.processor(jobs)
            output_files = dict(zip(keys, results))
        else:
            output_files = {}

        return {self.output_key : output_files}

    def get_job(self,
        dataset_key,
        dataset_files,
        map_resolution,
        ):

        o_dir = (
            self.dataset_dir / dataset_key / self.dataset_subdir
            )

        o_filename = (
            self.output_filename_template.format(dkey=dataset_key)
            )

        o_path = str(
            o_dir / o_filename
            )

        ###

        map_paths = []
        map_labels = []

        if 'dataset_map' in dataset_files:

            map_paths.append(
                dataset_files['dataset_map']
                )
            map_labels.append(
                'DATASET'
                )

        if 'ground_map' in dataset_files:

            map_paths.append(
                dataset_files['ground_map']
                )
            map_labels.append(
                'GROUND'
                )

        if 'z_map' in dataset_files:

            map_paths.append(
                dataset_files['z_map']
                )
            map_labels.append(
                'ZVALUES'
                )

        ###

        if len(map_paths) == 0:
            return None

        wrapper = self.get_wrapper(
            map_paths = map_paths,
            map_labels = map_labels,
            filepath = o_path,
            d_min = map_resolution,
            )

        return wrapper

    def get_wrapper(self,
        map_paths,
        map_labels,
        filepath,
        d_min,
        ):

        from giant.processors import ProcessWrapper
        return ProcessWrapper(
            func = self.make_mtz,
            map_paths = map_paths,
            map_labels = map_labels,
            filepath = filepath,
            d_min = d_min,
            )


class GetMakePanddaEvaluationMtzs:

    def __init__(self,
        get_dataset_map_writer,
        output_dir,
        dataset_dir = None,
        dataset_subdir = "",
        processor = None,
        ):

        from .maps import GetMakePanddaEvaluationMaps

        self.get_make_maps = GetMakePanddaEvaluationMaps(
            get_dataset_map_writer = get_dataset_map_writer,
            output_dir = output_dir,
            dataset_dir = dataset_dir,
            dataset_subdir = dataset_subdir,
            processor = processor,
            )

        if (processor is None):
            from giant.processors import basic_processor
            processor = basic_processor

        self.output_dir = output_dir
        self.dataset_dir = dataset_dir
        self.dataset_subdir = dataset_subdir

        self.processor = processor

    def __call__(self, label=""):

        return MakePanddaEvaluationMtzs(
            make_maps = self.get_make_maps(label=label),
            output_dir = (
                self.output_dir / label
                ),
            dataset_dir = (
                self.dataset_dir
                ),
            dataset_subdir = (
                self.dataset_subdir
                ),
            # processor = self.processor,
            )

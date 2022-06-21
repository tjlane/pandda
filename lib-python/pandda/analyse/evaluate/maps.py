import giant.logs as lg
logger = lg.getLogger(__name__)

import collections
import pathlib as pl

from pandda.analyse.events import background_correction as bc


class GetStandardSparseMapWriter(object):
    """docstring for GetDatasetMapMaker"""

    def __init__(self, 
        get_map_grid, 
        mask_name,
        ):
        
        self.get_map_grid = get_map_grid
        self.mask_name = mask_name

    def __call__(self):

        map_grid = self.get_map_grid()

        from giant.mulch.transform.maps.output import (
            MaskedNativeMapMaker,
            )

        return MaskedNativeMapMaker(
            map_grid = map_grid,
            map_mask = map_grid.masks[self.mask_name],
            )


class MakePanddaMapsJobGenerator(object):

    output_dataset_map = True
    output_ground_map = True
    output_dataset_z_maps = True
    output_dataset_event_maps = False

    def __init__(self, 
        write_map, 
        statistical_model, 
        dataset,
        dataset_map_data,
        bdc_values,
        output_dir,
        dataset_key = None,
        output_requires_events = True,
        ):

        self.write_map = write_map

        self.statistical_model = statistical_model

        self.dataset = dataset
        self.dataset_map_data = dataset_map_data
        self.dataset_key = dataset_key

        self.bdc_values = bdc_values

        self.output_files = None
        self.output_dir = pl.Path(output_dir)

        self.output_requires_events = output_requires_events

    def __call__(self, processor):

        processor(iter(self))

        return self.output_files

    def __iter__(self): 

        # Make this an attribute so that results can be picked up later
        self.output_files = {}

        if bool(self.output_requires_events) and not self.bdc_values:
            return

        if (self.output_dataset_map is True):

            filename = self.output_files.setdefault(
                "dataset_map", str(self.output_dir / "pandda_input_map.ccp4"),
                )

            yield self.get_job(
                map_data = self.dataset_map_data,
                map_path = filename,
                )

        if (self.output_ground_map is True):

            filename = self.output_files.setdefault(
                "ground_map", str(self.output_dir / "pandda_ground_state_map.ccp4"),
                )

            yield self.get_job(
                map_data = self.statistical_model.mu_map_data,
                map_path = filename,
                )

        if (self.output_dataset_z_maps is True): 

            filename = self.output_files.setdefault(
                "z_map", str(self.output_dir / "pandda_z_map.ccp4"),
                )

            yield self.get_job(
                map_data = self.statistical_model.get_z_map(
                    dataset_key = self.dataset_key,
                    map_data = self.dataset_map_data,
                    ),
                map_path = filename,
                )

        if (self.output_dataset_event_maps is True):

            event_maps = self.output_files.setdefault(
                "event_maps", {},
                )

            for i, bdc in enumerate(self.bdc_values):

                filename = event_maps.setdefault(
                    i, str(self.output_dir / "pandda_event_map_{:03d}.ccp4".format(i+1)),
                    )

                yield self.get_job(
                    map_data = bc.make_bdc_map(
                        query_map_data = self.dataset_map_data, 
                        background_map_data = self.statistical_model.mu_map_data, 
                        bdc_value = bdc, 
                        ),
                    map_path = filename,
                    )

    def get_job(self, map_data, map_path):

        from giant.processors import ProcessWrapper
        return ProcessWrapper(
            func = self.write_map,
            dataset = self.dataset,
            map_data = map_data,
            filename = map_path,
            )


###


class MakePanddaEvaluationMaps(object):

    output_key = "dataset_files"

    def __init__(self,
        get_dataset_map_writer,
        output_dir,
        dataset_dir = None,
        dataset_subdir = "",
        processor = None,
        output_requires_events = True,
        ):

        if (processor is None):
            from giant.processors import basic_processor
            processor = basic_processor

        if (dataset_dir is None):
            dataset_dir = output_dir

        self.get_dataset_map_writer = get_dataset_map_writer

        self.output_dir = pl.Path(output_dir)
        self.dataset_dir = dataset_dir
        self.dataset_subdir = dataset_subdir

        self.output_requires_events = output_requires_events
        
        self.processor = processor

    def __call__(self, 
        datasets,
        datasets_map_dict,
        datasets_results,
        statistical_model,
        ):
        
        logger('Writing maps... (this can take a long time!)')

        # Dataset-specific map_maker
        write_map = self.get_dataset_map_writer()

        # Iterators which will produce map writer jobs to be consumed
        dataset_keys = list(datasets.keys())
        dataset_generators = []

        for dkey in dataset_keys:

            events = datasets_results[dkey]

            dataset_generator = MakePanddaMapsJobGenerator(
                write_map = write_map, 
                statistical_model = statistical_model, 
                dataset = datasets[dkey],
                dataset_map_data = datasets_map_dict[dkey],
                bdc_values = [e['bdc'] for e in events],
                output_dir = str(self.dataset_dir / dkey / self.dataset_subdir),
                dataset_key = dkey, # allows use of map uncertainties cache
                output_requires_events = self.output_requires_events,
                )
            dataset_generators.append(dataset_generator)

        from pandda.wrappers import IteratorIterator
        job_iterator = IteratorIterator(
            iterators = dataset_generators,
            )

        self.processor(
            funcs = job_iterator,
            )

        # Extract output files for each dataset
        output_files = {}
        for dkey, d_generator in zip(dataset_keys, dataset_generators):
            # if d_generator.output_files:
            output_files[dkey] = d_generator.output_files

        return {self.output_key : output_files}


class GetMakePanddaEvaluationMaps(object):

    def __init__(self, 
        get_dataset_map_writer,
        output_dir,
        dataset_dir = None,
        dataset_subdir = "",
        processor = None,
        output_requires_events = True,
        ):
    
        if (processor is None):
            from giant.processors import basic_processor
            processor = basic_processor
            
        self.get_dataset_map_writer = get_dataset_map_writer

        self.output_dir = output_dir
        self.dataset_dir = dataset_dir
        self.dataset_subdir = dataset_subdir

        self.output_requires_events = output_requires_events

        self.processor = processor

    def __call__(self, label=""):

        return MakePanddaEvaluationMaps(
            get_dataset_map_writer = self.get_dataset_map_writer,
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
            output_requires_events = self.output_requires_events,
            )

import giant.logs as lg
logger = lg.getLogger(__name__)

import copy, functools, collections

import json
import pathlib as pl

from giant.utils import (
    make_sorted_dict,
    show_dict,
    merge_dicts,
    pretty_format_list,
    )

from .statistical_model import (
    GetPanddaStatisticalModelFitter,
    GetPanddaStatisticalModelOutputter,
    )

from .evaluate import (
    GetPanddaDatasetEvaluator,
    GetPanddaDatasetEvaluatorOutputter,
    )

# # Make this class instead of the chunked test sections in StandardPanddaModelProcessor
# class ProcessPanddaModelDatasets:

#     def __init__(self):
#         pass

#     def __call__(self):
#         pass


class PanddaModelProcessor(object):

    def __init__(self,
        load_maps,
        fit_model,
        evaluate_model,
        write_model_output,
        write_dataset_output,
        test_chunk_size = 50,
        label = None,
        output_json_path = None,
        ):

        self.load_maps = load_maps
        self.fit_model = fit_model
        self.evaluate_model = evaluate_model

        self.write_model_output = write_model_output
        self.write_dataset_output = write_dataset_output

        self.label = label

        self.test_chunk_size = test_chunk_size

        self.output_json_path = (
            pl.Path(output_json_path)
            if (output_json_path is not None) 
            else None
            )

        if (self.output_json_path is not None) and self.output_json_path.exists(): 
            raise IOError(
                'Output json file {path} already exists!'.format(
                    path = str(self.output_json_path),
                    )
                )

    def __call__(self,
        train_datasets,
        test_datasets,
        map_resolution,
        ):

        train_keys = list(
            train_datasets.keys()
            )
        test_keys = list(
            test_datasets.keys()
            )
        all_keys = list(
            set(train_keys).union(test_keys)
            )
        train_only_keys = list(
            set(train_keys).difference(test_keys)
            )
        train_and_test_keys = list(
            set(train_keys).intersection(test_keys)
            )

        all_datasets = {}
        all_datasets.update(train_datasets)
        all_datasets.update(test_datasets)

        # Single output dict
        all_results = {}

        logger.heading(
            'Starting Analysis of "{run}"'.format(
                run = self.label,
                )
            )

        logger(
            pretty_format_list(
                values = sorted(train_datasets.keys()),
                key = "> Training Datasets ({n})".format(
                    n = len(train_datasets),
                    ),
                ).strip().replace('\n','\n\t')
            )

        logger(
            pretty_format_list(
                values = sorted(test_datasets.keys()),
                key = "> Test Datasets ({n})".format(
                    n = len(test_datasets),
                    ),
                ).strip().replace('\n','\n\t')
            )

        ###
        # Prepare map loader (for xray, set data truncators based on keys to be used)
        logger.subheading('Preparing map loader')
        self.load_maps.initialise(
            datasets = all_datasets,
            )

        ###
        # Load maps
        logger.subheading('Loading Training Maps')
        training_map_manager = self.load_maps(
            datasets = train_datasets,
            map_resolution = map_resolution,
            )
        #
        ###

        ###
        # Fit statistical model
        logger.subheading('Fitting Statistical Model')
        fitted_model = self.fit_model(
            map_manager = training_map_manager,
            )
        ###

        ###
        # Reformat training maps for use in evaluation
        training_maps = training_map_manager.as_dict(
            delete_internal_references = True,
            )
        #
        # Remove unnecessary maps (training set only)
        for dkey in train_only_keys:
            training_maps.pop(dkey)
        #
        ###

        ###
        # Evaluate and output for all of the test datasets
        #
        # Evaluate in chunks (make chunks; already loaded maps first)
        test_dataset_chunks = self.chunk_order_datasets(
            datasets = test_datasets,
            first_keys = train_and_test_keys,
            )
        #
        for i_chunk, test_dataset_chunk in enumerate(test_dataset_chunks):

            ###
            # Get maps from dict or load as appropriate
            #
            logger.subheading(
                'Loading Test Maps (Chunk {i}/{n})'.format(
                    i=i_chunk+1, n=len(test_dataset_chunks),
                    )
                )
            #
            test_maps = self.load_maps_or_get_from_existing(
                datasets = test_dataset_chunk,
                existing_maps = training_maps,
                map_resolution = map_resolution,
                )
            #
            ###

            ###
            # Get the uncertainties for each dataset for output info (these are then cached)
            dataset_info = self.extract_dataset_info(
                datasets = test_dataset_chunk,
                datasets_map_dict = test_maps,
                statistical_model = fitted_model,
                map_resolution = map_resolution,
                )
            #
            merge_dicts(
                master_dict = all_results,
                merge_dict = {
                    'dataset_info' : make_sorted_dict(dataset_info),
                    },
                )
            #
            ###

            ###
            # Analyse each dataset against the parameterised statistical model
            #
            logger.subheading(
                'Evaluating Statistical Model (Chunk {i}/{n})'.format(
                    i=i_chunk+1, n=len(test_dataset_chunks),
                    )
                )
            #
            events_dict = self.evaluate_model(
                datasets = test_dataset_chunk,
                datasets_map_dict = test_maps,
                statistical_model = fitted_model,
                )
            #
            self.update_event_info(
                events_dict = events_dict,
                map_resolution = map_resolution,
                )
            #
            event_output_files = self.extract_output_files_from_events(
                events_dict = events_dict,
                )
            #
            show_dict(
                make_sorted_dict(events_dict),
                logger = logger,
                )
            #
            merge_dicts(
                master_dict = all_results,
                merge_dict = {
                    'events' : make_sorted_dict(events_dict),
                    'output_files' : make_sorted_dict(event_output_files),
                    },
                )
            #
            ###

            ###
            # Write output
            #
            logger.subheading(
                'Writing Output Pandda Files (Chunk {i}/{n})'.format(
                    i=i_chunk+1, n=len(test_dataset_chunks),
                    )
                )
            #
            output_files = self.write_dataset_output(
                datasets = test_dataset_chunk,
                datasets_map_dict = test_maps,
                datasets_results = events_dict,
                statistical_model = fitted_model,
                map_resolution = map_resolution,
                )
            #
            show_dict(
                make_sorted_dict(output_files),
                logger = logger,
                )
            #
            merge_dicts(
                master_dict = all_results,
                merge_dict = {
                    'output_files' : output_files,
                    },
                )
            #
            ###

        #####
        # Write model output
        #
        # Extract uncertainties
        train_dataset_map_uncertainties = {
            dkey : fitted_model.get_sigma_uncertainty(
                dataset_key = dkey,
                )
            for dkey in train_datasets.keys()
            }
        test_dataset_map_uncertainties = {
            dkey : (
                all_results['dataset_info']['map_uncertainty'][dkey]
                )
            for dkey in test_datasets.keys()
            }
        #
        output_files = self.write_model_output(
            statistical_model = fitted_model,
            train_dataset_map_uncertainties = train_dataset_map_uncertainties,
            test_dataset_map_uncertainties = test_dataset_map_uncertainties,
            )
        #
        merge_dicts(
            master_dict = all_results,
            merge_dict = {
                'shell_info' : {
                    'map_uncertainty' : train_dataset_map_uncertainties,
                    },
                'output_files' : {
                    'graphs' : {
                        'shells' : {
                            self.label : output_files,
                            },
                        },
                    },
                },
            )
        #
        ####

        return self.return_results(all_results)

    def extract_dataset_info(self,
        datasets,
        datasets_map_dict,
        statistical_model,
        map_resolution,
        ):

        dataset_info = {
            'map_uncertainty' : {
                dkey : (
                    round(
                        statistical_model.get_sigma_uncertainty(
                            map_data = datasets_map_dict[dkey],
                            dataset_key = dkey,
                            ),
                        3,
                        )
                    )
                for dkey in datasets.keys()
            },
            'analysed_resolution' : {
                dkey : (
                    round(map_resolution, 3)
                    )
                for dkey in datasets.keys()
            },
        }

        return dataset_info

    def extract_output_files_from_events(self,
        events_dict,
        ):

        output_files = {}

        for dkey, events in events_dict.items(): 

            for e in events: 

                ekey = e['event_num']

                for k, v in e.pop('output_files', {}).items():

                    # Invert the order of the dictionary keys
                    of = output_files.setdefault(k,{}).setdefault(dkey,{})
                    of[ekey] = v

        return {'graphs' : output_files}

    def update_event_info(self,
        events_dict,
        map_resolution,
        ):

        for dkey, d_events in events_dict.items():

            for e in d_events: 
                
                e['map_resolution'] = (
                    round(map_resolution, 3)
                    )

    def chunk_order_datasets(self, datasets, first_keys=None):
        """Create processing chunks"""

        if (first_keys is None):
            first_keys = []

        chunk_size = self.test_chunk_size

        # Chunk list and first chunk
        output_chunks = []
        current_chunk = None

        # Generate set of all possible keys
        dataset_keys = set(datasets.keys())

        # Overlap between selected keys and all dataset keys
        first_intersection = set(first_keys).intersection(dataset_keys)

        if len(first_intersection) > 0:
            # Put the preferred keys in first
            for dtag in first_intersection:
                # New chunk?
                if (current_chunk is None) or (len(current_chunk) >= chunk_size):
                    current_chunk = {}
                    output_chunks.append(current_chunk)
                # Add to chunk
                current_chunk[dtag] = datasets[dtag]
            # Remove first keys from set
            dataset_keys = dataset_keys.difference(first_keys)
            # Force a new set
            current_chunk = None

        # Add rest of the keys
        for dtag in dataset_keys:
            # New chunk?
            if (current_chunk is None) or (len(current_chunk) >= chunk_size):
                current_chunk = {}
                output_chunks.append(current_chunk)
            # Add to chunk
            current_chunk[dtag] = datasets[dtag]

        return output_chunks

    def load_maps_or_get_from_existing(self,
        datasets,
        existing_maps,
        map_resolution,
        ):

        # Dataset maps for this block
        dataset_maps = {}
        # Datasets to load maps for
        datasets_for_loading = {}
        #
        for dtag, dataset in datasets.items():
            if dtag in existing_maps:
                dataset_maps[dtag] = existing_maps.pop(dtag)
            else:
                datasets_for_loading[dtag] = dataset
        #
        if len(datasets_for_loading) > 0:
            dataset_maps_rest = self.load_maps(
                datasets = datasets_for_loading,
                map_resolution = map_resolution,
                ).as_dict(
                delete_internal_references = True,
                )
            dataset_maps.update(dataset_maps_rest)

        return dataset_maps

    def return_results(self, results_dict):

        if self.output_json_path is not None: 

            assert not self.output_json_path.exists()

            # Write results to json file
            with open(str(self.output_json_path), 'w') as fh:

                fh.write(
                    json.dumps(
                        results_dict,
                        indent = 2,
                        )
                    )

            logger(
                'Output written to {}'.format(
                    str(self.output_json_path)
                    )
                )

            return None

        return results_dict


class GetPanddaModelProcessor(object):

    def __init__(self,
        get_model_fitter,
        get_model_outputter,
        get_dataset_evaluator,
        get_dataset_outputter,
        ):

        self.get_model_fitter = get_model_fitter
        self.get_model_outputter = get_model_outputter
        self.get_dataset_evaluator = get_dataset_evaluator
        self.get_dataset_outputter = get_dataset_outputter

    def __call__(self,
        load_maps,
        run_label,
        output_json_path = None,
        ):

        return PanddaModelProcessor(
            load_maps = load_maps,
            fit_model = self.get_model_fitter(),
            evaluate_model = self.get_dataset_evaluator(),
            write_model_output = self.get_model_outputter(
                label = run_label,
                ),
            write_dataset_output = self.get_dataset_outputter(
                label = run_label,
                ),
            label = run_label,
            output_json_path = output_json_path,
            )


#####


class RunPanddaModel(object):

    def __init__(self,
        get_pandda_model_task,
        partition_shells,
        output_dir,
        processor = None,
        ):

        if (processor is None):
            from giant.processors import basic_processor
            processor = basic_processor

        self.get_pandda_model_task = get_pandda_model_task
        self.partition_shells = partition_shells
        self.output_dir = output_dir
        self.processor = processor

        if not self.output_dir.exists():
            self.output_dir.mkdir(parents=True)

        assert self.output_dir.is_dir()

    def __call__(self,
        datasets,
        load_maps,
        ):

        ###
        #
        # Processing objects
        shell_processors = []
        #
        # Results objects to be returned
        shell_records = []
        #
        ###

        logger.heading('Running Core Pandda Analysis', spacer=True)

        logger.subheading('Partitioning datasets into shells')
        logger(self.partition_shells)

        dataset_shells = self.partition_shells(
            datasets = datasets,
            )

        for (shell_res_high, shell_res_low), shell_datasets in dataset_shells.items():

            run_label = 'shell-{high}-{low}'.format(
                high = str(shell_res_high),
                low = str(shell_res_low),
                )

            output_json_path = (
                self.output_dir / "{}-results.json".format(run_label)
                )

            run_pandda_model = self.get_pandda_model_task(
                load_maps = load_maps,
                run_label = run_label,
                output_json_path = output_json_path,
                )

            shell_processors.append(
                self.processor.make_wrapper(
                    func = run_pandda_model,
                    train_datasets = shell_datasets['train'],
                    test_datasets = shell_datasets['test'],
                    map_resolution = shell_res_low, # get dynamically? resolution of the lowest dataset?
                    )
                )

            # Record which datasets are test/train at each resolution
            shell_records.append(
                collections.OrderedDict([
                    (
                        'label', run_label,
                        ),
                    (
                        'json_path', str(output_json_path),
                        ),
                    (
                        'resolution_high', round(shell_res_high, 3),
                        ),
                    (
                        'resolution_low', round(shell_res_low, 3),
                        ),
                    (
                        'map_resolution', round(shell_res_low, 3),
                        ),
                    (
                        'train', sorted(shell_datasets['train'].keys()),
                        ),
                    (
                        'test', sorted(shell_datasets['test'].keys()),
                        ),
                    ])
                )

        ###
        # Report
        logger.subheading('Test/Train datasets per shell')
        self.log_shells(shell_records)
        #
        ###

        ###
        # Run !
        logger.heading('Running analyses')
        #
        shell_results = self.processor(
            funcs = shell_processors,
            )
        #
        ###

        ###
        #
        return_dict = self.unpack_results(
            #shell_results = shell_results, # shell_results now passed via json file instead
            shell_records = shell_records,
            )
        #
        ###

        logger.subheading('all return information')
        show_dict(return_dict, logger=logger)

        return return_dict

    def unpack_results(self,
        shell_records,
        ):

        # concatenate into one dict/list
        all_events = []
        all_output_files = {}
        all_dataset_records = {}

        for i, s_records in enumerate(shell_records):

            # shell results now passed via json
            s_json = pl.Path(
                s_records['json_path']
                )
            assert s_json.exists()

            with open(str(s_json), 'r') as fh:
                s_results = json.loads(fh.read())

            ###
            # Output files for this run
            #
            merge_dicts(
                master_dict = all_output_files,
                merge_dict = s_results['output_files'],
                )
            #
            ###

            ###
            # Merge shell records
            #
            merge_dicts(
                master_dict = s_records,
                merge_dict = s_results['shell_info'],
                )

            ###
            # Merge dataset info
            #
            merge_dicts(
                master_dict = all_dataset_records,
                merge_dict = s_results['dataset_info'],
                )
            #
            ###

            ###
            # Merge event lists 
            #
            for dkey, d_events in s_results['events'].items():
                all_events.extend(d_events)
            #
            ###

        return_dict = {
            'events' : sorted(
                all_events, 
                key = lambda e: (e['dtag'], e['event_num'])
                ),
            'shell_records' : (
                shell_records
                ),
            'dataset_records' : make_sorted_dict(
                all_dataset_records
                ),
            'output_files' : make_sorted_dict(
                all_output_files, 
                recursive = True,
                ),
            }

        return return_dict

    def log_shells(self, shell_records):

        for i, d in enumerate(shell_records):

            logger('> Processing Shell {i}'.format(i=i+1))

            for k,v in d.items():

                if k in ['test','train']:
                    key = "{k} ({n})".format(k=k, n=len(v))
                    s = pretty_format_list(values=v, key=key)
                else:
                    s = '{k} : {v}'.format(k=k, v=str(v))

                logger(
                    '| {s}'.format(
                        s = s.strip().replace('\n','\n|\t'),
                        )
                    )

            logger('`---->\n')

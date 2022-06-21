import collections

from .output import (
    GetPanddaDatasetEvaluatorOutputter
    )


class PanddaEvaluateSingleDataset(object):

    def __init__(self,
        map_grid,
        map_mask,
        search_mask,
        find_events,
        analyse_event,
        ): 

        self.map_grid = map_grid

        self.map_mask = map_mask
        self.search_mask = search_mask

        self.find_events = find_events
        self.analyse_event = analyse_event

    def __call__(self,
        dataset_key,
        dataset, 
        dataset_map_data,
        statistical_model,
        ):

        events = self.find_events(
            map_grid = self.map_grid,
            map_mask = self.map_mask,
            search_mask = self.search_mask,
            dataset_key = dataset_key,
            dataset = dataset,
            dataset_map_data = dataset_map_data,
            statistical_model = statistical_model,
            )

        dataset_info = self.get_dataset_event_info(
            dataset_key = dataset_key,
            dataset = dataset, 
            dataset_map_data = dataset_map_data,
            statistical_model = statistical_model,
            )

        for e in events: # number the event here! 

            # Add the constant information
            e.update(dataset_info)

            # Analyse event to find BDC, etc.
            e_info = self.analyse_event(
                event = e,
                dataset = dataset,
                dataset_map_data = dataset_map_data,
                reference_map_data = statistical_model.mu_map_data,
                map_mask = self.map_mask,
                map_grid = self.map_grid,
                )
            e.update(e_info)

            # This does not need to be returned
            # Remove to same time with pickling if nothing else
            e.pop('cluster')

        return events

    def get_dataset_event_info(self, 
        dataset_key,
        dataset, 
        dataset_map_data,
        statistical_model,
        ):

        d = collections.OrderedDict([
            (
                'map_uncertainty',
                round(
                    statistical_model.get_sigma_uncertainty(
                        dataset_key = dataset_key, 
                        map_data = dataset_map_data, # only needed if values are not cached! 
                        ), 
                    3,
                    ),
                ),
            ]) 

        return d


class PanddaDatasetEvaluator(object):

    def __init__(self,
        evaluate_dataset,
        processor = None,
        ):

        if (processor is None):
            from giant.processors import basic_processor
            processor = basic_processor

        self.evaluate_dataset = evaluate_dataset
        self.processor = processor

    def __call__(self,
        datasets,
        datasets_map_dict,
        statistical_model,
        ):

        dataset_keys = list(datasets.keys())
        jobs = []

        for dkey in dataset_keys:

            jobs.append(
                self.processor.make_wrapper(
                    func = self.evaluate_dataset,
                    dataset_key = dkey,
                    dataset = datasets[dkey],
                    dataset_map_data = datasets_map_dict[dkey],
                    statistical_model = statistical_model,
                    )
                )

        results = self.processor(jobs)

        return dict(zip(dataset_keys, results))


class GetPanddaDatasetEvaluator(object):

    def __init__(self, 
        get_map_grid,
        find_events, 
        analyse_event,
        map_mask_name,
        search_mask_name,
        processor = None,  
        ):

        self.get_map_grid = get_map_grid
        self.map_mask_name = map_mask_name
        self.search_mask_name = search_mask_name

        self.find_events = find_events
        self.analyse_event = analyse_event

        self.processor = processor
        
    def __call__(self):

        map_grid = self.get_map_grid()

        return PanddaDatasetEvaluator(
            evaluate_dataset = PanddaEvaluateSingleDataset(
                map_grid = map_grid,
                map_mask = map_grid.masks[self.map_mask_name],
                search_mask = map_grid.masks[self.search_mask_name],
                find_events = self.find_events,
                analyse_event = self.analyse_event,
                ),
            processor = self.processor,
            )

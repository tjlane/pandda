from __future__ import absolute_import
import collections

import numpy as np

from . import background_correction as bc

from .graphs import (
    MakeBdcEstimationGraph
    )


class EventAnalyser(object):

    def __init__(self,
        min_bdc = 0.0,
        max_bdc = 1.0, 
        increment = 0.01, 
        output_multiplier = 1.0,
        output_directory = None,
        ):

        self.min_bdc = min_bdc
        self.max_bdc = max_bdc
        self.increment = increment
        self.output_multiplier = output_multiplier

        self.write_image = (
            MakeBdcEstimationGraph(
                output_path_template = str(
                    output_directory / "{label}-event{event_num:03d}.png"
                    ),
                )
            if (output_directory is not None)
            else None
            )

    def __call__(self,
        event,
        dataset,
        dataset_map_data,
        reference_map_data,
        map_mask,
        map_grid,
        output_directory = None,
        ):

        # Get the mask for the area around the event
        event_mask = self.get_event_mask(
            event = event, 
            map_grid = map_grid,
            )

        # Index the event mask on the analysis mask
        event_mask_reindex, event_mask_reindex_sel = event_mask.index_on_other(
            other = map_mask,
            )

        all_event_fractions,\
        all_event_correlations, \
        all_global_correlations = self.calculate_bdc_series_correlations(
            reference_map_data = reference_map_data,
            dataset_map_data = dataset_map_data,
            event_selection = event_mask_reindex_sel,
            )

        event_fraction = self.estimate_event_fraction(
            event_fractions = all_event_fractions,
            event_correlations = all_event_correlations,
            global_correlations = all_global_correlations,
            )

        # Could replace this with find-in-list function
        local_correlation, global_correlation = self.calculate_single_correlation(
            event_fraction = event_fraction,
            reference_map_data = reference_map_data,
            dataset_map_data = dataset_map_data,
            event_selection = event_mask_reindex_sel,
            )

        output_files = (
            self.write_image(
                dataset = dataset,
                event_num = event['event_num'],
                event_fractions = all_event_fractions,
                event_correlations = all_event_correlations,
                global_correlations = all_global_correlations,
                )
            if (self.write_image is not None)
            else {}
            )

        # TODO move this to other function
        # atm = find_nearest_atoms(
        #     atoms=list(protein(dataset.model.hierarchy).atoms_with_labels()),
        #     query=dataset.model.alignment.ref2nat(
        #         grid.grid2cart([event.cluster.centroid])
        #         )[0],
        #     )

        event_info = collections.OrderedDict([
            (
                'event_fraction', round(event_fraction,3),
                ),
            (
                'bdc', round(1.0-event_fraction,3),
                ),
            (
                'global_correlation', round(global_correlation, 3),
                ),
            (
                'local_correlation', round(local_correlation, 3),
                ),
            (
                'output_files', output_files,
                ),
            ])

        return event_info

    def get_event_mask(self, event, map_grid):

        from giant.mulch.transform.maps.grid.mask import GetSitesMask

        get_event_mask = GetSitesMask.from_map_grid(
            map_grid = map_grid,
            mask_dist = 2.0,
            )

        event_mask = get_event_mask(
            sites_cart = map_grid.grid2cart(
                event['cluster'].points,
                ),
            )

        return event_mask

    def calculate_bdc_series_correlations(self,
        reference_map_data,
        dataset_map_data,
        event_selection,
        ):

        # Calculate the correlations for different BDCs
        event_fractions, \
        event_correlations, \
        global_correlations = bc.calculate_feature_fraction_correlations(
            query_map_data = dataset_map_data,
            background_map_data = reference_map_data,
            feature_mask = event_selection,
            min_fraction = (1.0 - self.max_bdc),
            max_fraction = (1.0 - self.min_bdc),
            increment = self.increment,
            )

        return (
            event_fractions, 
            event_correlations, 
            global_correlations,
            )

    def estimate_event_fraction(self,
        event_fractions,
        event_correlations,
        global_correlations,
        ):
        
        # Find the point of maximal descrepancy
        estimated_fraction = bc.find_maximum_series_discrepancy(
            labels = event_fractions,
            series_1 = global_correlations,
            series_2 = event_correlations,
            )

        # Apply output multiplier
        estimated_fraction = (estimated_fraction * self.output_multiplier)

        # Check still within bounds
        estimated_fraction = min( # checking against the upper limit
            estimated_fraction,
            1.0 - self.min_bdc, # the maximum possible 
            )

        estimated_fraction = max( # checking against the lower limit
            estimated_fraction,
            1.0 - self.max_bdc, # the minimum possible 
            )

        return estimated_fraction

    def calculate_single_correlation(self,
        event_fraction,
        reference_map_data,
        dataset_map_data,
        event_selection,
        ):

        event_map_data = bc.make_bdc_map(
            query_map_data = dataset_map_data, 
            background_map_data = reference_map_data, 
            bdc_value = (1.0 - event_fraction),
            )

        global_correlation = np.corrcoef(
            reference_map_data, 
            dataset_map_data,
            )[0, 1]

        local_correlation = np.corrcoef(
            reference_map_data[event_selection],
            dataset_map_data[event_selection],
            )[0, 1]

        return local_correlation, global_correlation

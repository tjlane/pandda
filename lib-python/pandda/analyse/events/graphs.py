import giant.logs as lg
logger = lg.getLogger(__name__)

import collections
import numpy as np

from pandda.graphs import (
    PanddaDatasetPlotter,
    )


class MakeBdcEstimationGraph(PanddaDatasetPlotter):

    output_key = "bdc_estimation"

    def __call__(self, dataset, event_num, *args, **kwargs):

        dataset_label = self.get_label(dataset)

        fig = self.plot(
            dataset = dataset,
            *args, **kwargs
            )

        filename = self.get_path(
            dataset_label = dataset_label,
            event_num = event_num,
            )

        self.save(
            fig = fig,
            filename = filename,
            )

        return {self.output_key : filename}

    def plot(self,
        event_fractions,
        event_correlations,
        global_correlations,
        *args, **kwargs
        ):

        diff_values = tuple(
            np.array(global_correlations) - np.array(event_correlations)
            )

        # Get the peak location
        x_values = event_fractions
        max_x = x_values[diff_values.index(max(diff_values))]

        ###

        fig, axes = self.setup(ncols=1, nrows=2)

        a1, a2 = tuple(axes)

        #

        line_1_1, = a1.plot(
            x_values, global_correlations, 
            'g--', 
            label = 'Global corr.', 
            )
        line_1_2, = a1.plot(
            x_values, event_correlations, 
            'k--', 
            label='Local corr.', 
            )

        a1.set_ylabel('Corr. to\nground state', color='k')
        a1.set_ylim((-1, 1))

        #

        line_2_1, = a2.plot(
            x_values, diff_values, 
            'b-', 
            label='Difference', 
            linewidth = 1,
            )

        a2.set_xlabel(
            'Event fraction (1-BDC)', 
            color='k',
            )
        a2.set_ylabel(
            'Corr. \ndifference', 
            color='k',
            )
        a2.set_ylim(
            (min(diff_values)-0.2, max(diff_values)+0.2)
            )

        #

        # Plot line at the maximum
        _ = a1.plot([max_x,max_x],[-1,1], 'k-', linewidth=1)
        _ = a2.plot([max_x,max_x],[-1,1], 'k-', linewidth=1)

        text_2_1 = a2.text(
            0.02 + max_x, 
            min(diff_values), 
            'BDC={}'.format(1-max_x), 
            )

        #

        a1.legend(
            handles = [line_1_1, line_1_2, line_2_1], 
            loc = 4,
            )

        return fig
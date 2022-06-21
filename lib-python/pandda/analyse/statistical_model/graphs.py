import numpy as np

from pandda.graphs import (
    PanddaPlotter,
    )

from pandda.utils import (
    show_dict, 
    merge_dicts,
    )


class StatisticalMapDistributionsPlotter(PanddaPlotter):

    output_key = "statistical_map_distribution"

    def plot(self,
        statistical_model,
        *args, **kwargs
        ):

        mean_map_values = (
            statistical_model.mu_map_data
            )
        sadj_map_values = (
            statistical_model.sigma_adjusted_map_data
            )

        ###

        fig, axes = self.setup(ncols=2, nrows=1)

        a1, a2 = tuple(axes)

        a1.set_title('Mean Map Distribution')
        a1.set_xlabel('Map Value')
        a1.set_ylabel('Density')

        self.make_histogram(
            axis = a1,
            values = mean_map_values,
            )

        a2.set_title('SAdj Map Distribution')
        a2.set_xlabel('Map Value')
        a2.set_ylabel('Density')

        self.make_histogram(
            axis = a2,
            values = sadj_map_values,
            )

        fig.subplots_adjust(hspace=0.4)

        return fig

    def make_histogram(self,
        axis, 
        values,
        n_bins = 60,
        ):

        axis.hist(x=values, bins=n_bins, density=True)


class StatisticalMapUncertaintyDistributionPlotter(PanddaPlotter):

    output_key = "uncertainty_distribution"

    def plot(self,
        statistical_model,
        train_dataset_map_uncertainties,
        test_dataset_map_uncertainties,
        *args, **kwargs
        ):

        train_uncertainties = list(train_dataset_map_uncertainties.values()) 

        test_uncertainties = list(test_dataset_map_uncertainties.values())

        x_range = (
            0.,
            # min(
            #     train_uncertainties.min(),
            #     test_uncertainties.min(),
            #     ),
            np.max([
                np.max(train_uncertainties),
                np.max(test_uncertainties),
                ]),
            )

        ###

        fig, axes = self.setup(ncols=1, nrows=2, sharex=True)

        a1, a2 = tuple(axes)

        ###
        
        a1.set_title("Training Dataset Uncertainties")
        a1.set_xlabel("Map Uncertainty")
        a1.set_ylabel("Count")

        self.make_histogram(
            axis = a1,
            values = train_uncertainties,
            facecolor = 'gray',
            label = "train datasets",
            range = x_range,
            )

        a1.legend(
            loc = "upper right",
            bbox_to_anchor = (0.9, 0.9),
            )

        ###
        
        a2.set_title("Test Dataset Uncertainties")
        a2.set_xlabel("Map Uncertainty")
        a2.set_ylabel("Count")
        
        self.make_histogram(
            axis = a2,
            values = test_uncertainties,
            facecolor = 'blue',
            label = "test datasets",
            range = x_range,
            )

        a2.legend(
            loc = "upper right",
            bbox_to_anchor = (0.9, 0.9),
            )

        return fig

    def make_histogram(self,
        axis, 
        values,
        n_bins = 30,
        **kwargs
        ):

        axis.hist(x=values, bins=n_bins, density=False, **kwargs)


###


class MakePanddaStatisticalModelGraphs(object):

    def __init__(self, 
        output_dir,
        processor = None,
        ):

        if (processor is None):
            from giant.processors import basic_processor
            processor = basic_processor

        self.plotters = [
            StatisticalMapDistributionsPlotter(
                output_path = str(
                    output_dir / "statistical_map_distributions.png"
                    ),
                ),
            StatisticalMapUncertaintyDistributionPlotter(
                output_path = str(
                    output_dir / "dataset_uncertainty_distribution.png"
                    ),
                ),
            ]

        self.processor = processor

    def __call__(self, 
        statistical_model,
        *args, **kwargs
        ):

        output_files = {}

        for plotter in self.plotters:

            of = plotter(
                statistical_model = statistical_model,
                *args, **kwargs
                )

            merge_dicts(
                master_dict = output_files,
                merge_dict = of,
                )

        return output_files


class GetMakePanddaStatisticalModelGraphs(object): 

    def __init__(self, 
        output_dir,
        processor = None,
        ):

        if (processor is None):
            from giant.processors import basic_processor
            processor = basic_processor

        self.output_dir = output_dir
        self.processor = processor

    def __call__(self, label):

        return MakePanddaStatisticalModelGraphs(
            output_dir = (
                self.output_dir / label
                ),
            processor = self.processor,
            )
        

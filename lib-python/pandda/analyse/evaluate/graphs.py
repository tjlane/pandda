import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy as np
import scipy.stats as st

from pandda.utils import (
    merge_dicts,
    )

from pandda.graphs import (
    PanddaDatasetPlotter,
    PanddaMultiDatasetPlotter,
    )


class MapDistributionPlotter(PanddaDatasetPlotter):

    output_key = "map_distribution"

    def plot(self, 
        dataset_label,
        datasets_map_dict, 
        statistical_model, 
        *args, **kwargs
        ):

        difference_map_values = (
            datasets_map_dict.get(dataset_label) - 
            statistical_model.mu_map_data
            )

        dataset_uncertainty = statistical_model.get_sigma_uncertainty(
            map_data = datasets_map_dict.get(dataset_label),
            dataset_key = dataset_label,
            )

        ###

        fig, axis = self.setup()

        axis.set_title('Distribution of difference map values')
        axis.set_xlabel('Map value')
        axis.set_ylabel('Density')

        self.make_histogram(
            axis = axis,
            values = difference_map_values,
            )

        self.make_normal_distribution(
            axis = axis,
            mean_std = (
                difference_map_values.mean(), 
                difference_map_values.std(),
                ),
            linestyle = 'k--',
            )

        self.make_normal_distribution(
            axis = axis,
            mean_std = (
                0.0, 
                dataset_uncertainty,
                ),
            linestyle = 'b-',
            )

        return fig

    def make_histogram(self,
        axis, 
        values,
        n_bins = 60,
        ):

        axis.hist(x=values, bins=n_bins, density=True)

    def make_normal_distribution(self, 
        axis,
        mean_std,
        z_range = (-5.,5.),
        n_points = 101,
        linestyle = 'k--',
        ):

        # Plot the distribution for the observed distribution
        x_vals = np.linspace(
            z_range[0] * mean_std[1], 
            z_range[1] * mean_std[1], 
            n_points,
            )
        y_vals = st.norm.pdf(
            x_vals, 
            loc = mean_std[0], 
            scale = mean_std[1],
            )
        axis.plot(x_vals, y_vals, linestyle)

        return None     


class MapScatterPlotter(PanddaDatasetPlotter):

    output_key = "map_scatter"

    def plot(self,
        dataset_label,
        datasets_map_dict,
        statistical_model,
        *args, **kwargs
        ):

        dataset_map_values = datasets_map_dict.get(dataset_label)
        mean_map_values = statistical_model.mu_map_data

        ###

        fig, axes = self.setup(ncols=2, nrows=1)

        a1, a2 = tuple(axes)

        a1.set_title('Unsorted Map Values')
        a1.set_xlabel('Dataset map')
        a1.set_ylabel('Mean map')

        self.make_scatter(
            axis = a1,
            map_1 = dataset_map_values,
            map_2 = mean_map_values,
            )

        a2.set_title('Sorted Map Values')
        a2.set_xlabel('Dataset map')
        a2.set_ylabel('Mean map')

        self.make_scatter(
            axis = a2,
            map_1 = np.sort(dataset_map_values),
            map_2 = np.sort(mean_map_values),
            )        

        return fig

    def make_scatter(self, axis, map_1, map_2):

        min_v = min(min(map_1), min(map_2))
        max_v = max(max(map_1), max(map_2))

        axis.plot([min_v, min_v], [max_v, max_v], 'b--')
        axis.plot(map_1, map_2, 'go')


class QQDistributionPlotter(PanddaDatasetPlotter):

    output_key = "qq_plot"

    def plot(self,         
        dataset_label,
        datasets_map_dict, 
        statistical_model, 
        *args, **kwargs
        ):

        difference_map_values = (
            datasets_map_dict.get(dataset_label) - 
            statistical_model.mu_map_data
            )

        map_uncertainty = statistical_model.get_sigma_uncertainty(
            map_data = datasets_map_dict.get(dataset_label),
            dataset_key = dataset_label,
            )

        ###

        fig, axis = self.setup()

        axis.set_title('Quantile-Quantile plot of dataset difference from mean')
        axis.set_xlabel('Observed quantiles')
        axis.set_ylabel('Theoretical quantiles')

        self.make_qq_plot(
            axis = axis,
            difference_map_values = difference_map_values,
            map_uncertainty = map_uncertainty,
            )

        return fig

    def make_qq_plot(self,
        axis,
        difference_map_values,
        map_uncertainty,
        ):
        
        o_values = np.sort(difference_map_values)

        p_values = np.linspace(0., 1., len(difference_map_values)+2)
        e_values = st.norm.ppf(p_values, 0., 1.)[1:-1] # chop ends which are inf

        e_min, e_max = (min(e_values), max(e_values))

        # Line if observed is N(0,1)
        # axis.plot(
        #     [e_min, e_max], 
        #     [e_min, e_max], 
        #     'b--',
        #     )

        # Line if observed is N(0,map_uncertainty)
        axis.plot(
            [e_min*map_uncertainty, e_max*map_uncertainty], 
            [e_min, e_max], 
            'b--',
            )

        # Observed values against N(0,1)
        axis.plot(
            o_values, 
            e_values, 
            'go-',
            )


###


class MakePanddaEvaluationGraphs(object):

    output_key = "graphs"
    
    def __init__(self, 
        output_dir,
        dataset_dir = None,
        dataset_subdir = "{label}",
        processor = None,
        ):

        if (processor is None):
            from giant.processors import basic_processor
            processor = basic_processor

        if (dataset_dir is None):
            dataset_dir = output_dir

        d_dir = (dataset_dir / dataset_subdir)

        self.plotters = [
            MapDistributionPlotter(
                output_path_template = str(
                    d_dir / "{label}-distribution.png"
                    ),
                ),
            MapScatterPlotter(
                output_path_template = str(
                    d_dir / "{label}-scatter.png"
                    ),
                ),
            QQDistributionPlotter(
                output_path_template = str(
                    d_dir / "{label}-qqplot.png"
                    ),
                ),
            ]

        self.processor = processor

    def __call__(self, 
        datasets, 
        datasets_map_dict, 
        statistical_model,
        *args, **kwargs
        ):

        output_files = {}

        for plotter in self.plotters:

            of = plotter(
                datasets = datasets, 
                datasets_map_dict = datasets_map_dict,
                statistical_model = statistical_model,
                *args, **kwargs
                )

            merge_dicts(
                master_dict = output_files,
                merge_dict = of,
                )

        return {self.output_key : output_files}


class GetMakePanddaEvaluationGraphs(object): 

    def __init__(self, 
        output_dir,
        dataset_dir = None,
        dataset_subdir = "dataset_graphs",
        processor = None,
        ):

        if (processor is None):
            from giant.processors import basic_processor
            processor = basic_processor
            
        self.output_dir = output_dir
        self.dataset_dir = dataset_dir
        self.dataset_subdir = dataset_subdir

        self.processor = processor

    def __call__(self, label):

        output_dir = (self.output_dir / label)

        if (self.dataset_dir is not None):
            dataset_dir = self.dataset_dir
        else: 
            dataset_dir = output_dir

        return MakePanddaEvaluationGraphs(
            output_dir = (
                output_dir
                ),
            dataset_dir = (
                dataset_dir
                ),
            dataset_subdir = (
                self.dataset_subdir
                ),
            processor = self.processor,
            )



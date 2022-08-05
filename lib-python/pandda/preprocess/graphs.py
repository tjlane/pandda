import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy as np

from pandda.graphs import (
    PanddaDatasetPlotter,
    PanddaMultiDatasetPlotter,
    )

from pandda.utils import (
    show_dict, 
    merge_dicts,
    )


class WilsonPlotter(PanddaDatasetPlotter):

    output_key = "wilson_plot"

    def plot(self,
        dataset,
        dataset_label,
        data_getters,
        reference_miller_array = None,
        *args, **kwargs
        ):

        get_dataset_data = data_getters.get(dataset_label)

        ###

        fig, axis = self.setup()

        axis.set_title('Dataset wilson plot')
        axis.set_xlabel('Resolution ($\\AA$)')
        axis.set_ylabel('ln(mean amplitude)')

        if (reference_miller_array is not None):
            self.wilson_plot(
                axis = axis,
                miller_array = reference_miller_array,
                linestyle = 'r--',
                )

        self.wilson_plot(
            axis = axis,
            miller_array = get_dataset_data(dataset),
            linestyle = 'k-',
            )

        # Use to avoid warnings
        axis.set_xticks(
            axis.get_xticks()
            )

        axis.set_xticklabels([
            round(t**-0.5, 1) 
            if t >0 
            else "" 
            for t in axis.get_xticks()
            ])

        return fig

    def wilson_plot(self, 
        axis,
        miller_array, 
        linestyle="k-",
        ):

        x, y = self.get_wilson_plot_values(
            miller_array.as_amplitude_array(),
            )

        axis.plot(x, y, linestyle, linewidth=1)

    def get_wilson_plot_values(self, 
        miller_array,
        ):
        
        binner = miller_array.setup_binner(auto_binning=True)
        binned = miller_array.wilson_plot(use_binning=True)
        x_bin_cent = binner.bin_centers(1)
        y_bin_data = binned.data[1:-1]
        assert len(x_bin_cent) == len(y_bin_data)
        return np.power(x_bin_cent,2), np.log(y_bin_data)


class ResolutionDistributionPlotter(PanddaMultiDatasetPlotter):

    output_key = "resolution_distribution"

    def plot(self,
        dataset_labels,
        dataset_list,
        dataset_dicts,
        *args, **kwargs
        ):

        resolution_high_values = list(dataset_dicts['high_resolution'].values())
        
        resolution_low_values = list(dataset_dicts['low_resolution'].values())

        ###

        fig, axes = self.setup(nrows=2)

        fig.suptitle('Dataset resolution histograms')

        #

        a1 = axes[0]
        a2 = axes[1]
        
        #

        a1.set_xlabel('High Resolution Limit ($\\AA$)')
        a1.set_ylabel('Count')

        self.make_histogram(
            axis = a1,
            values = resolution_high_values,
            )

        #

        a2.set_xlabel('Low Resolution Limit ($\\AA$)')
        a2.set_ylabel('Count')

        self.make_histogram(
            axis = a2,
            values = resolution_low_values,
            )

        #

        fig.subplots_adjust(hspace=0.3)
    
        return fig

    def make_histogram(self,
        axis, 
        values,
        n_bins = 30,
        ):

        axis.hist(x=values, bins=n_bins, density=False)


class RFactorDistributionPlotter(PanddaMultiDatasetPlotter):

    output_key = "rfactor_distribution"

    def plot(self,
        dataset_labels,
        dataset_list,
        dataset_dicts,
        *args, **kwargs
        ):

        r_free_values = list(dataset_dicts['r_free'].values())

        r_work_values = list(dataset_dicts['r_work'].values())

        ###

        fig, axes = self.setup(nrows=2)

        fig.suptitle('Dataset R-factor histograms')

        a1 = axes[0]
        a2 = axes[1]

        #

        a1.set_xlabel('R-free Values ($\\AA$)')
        a1.set_ylabel('Count')

        self.make_histogram(
            axis = a1,
            values = r_free_values,
            )

        #

        a2.set_xlabel('R-work Values ($\\AA$)')
        a2.set_ylabel('Count')

        self.make_histogram(
            axis = a2,
            values = r_work_values,
            )

        #

        fig.subplots_adjust(hspace=0.3)
    
        return fig

    def make_histogram(self,
        axis, 
        values,
        n_bins = 30,
        ):

        axis.hist(x=values, bins=n_bins, density=False)


class UnitCellDistributionPlotter(PanddaMultiDatasetPlotter):

    output_key = "unit_cell_distribution"

    def plot(self,
        dataset_labels,
        dataset_list,
        dataset_dicts,
        *args, **kwargs
        ):

        unit_cell_parameters = list(zip(*list(dataset_dicts['unit_cell'].values())))

        ###

        fig, axes = self.setup(nrows=2, ncols=3)

        fig.suptitle('Unit Cell histograms')

        #

        a1 = axes[0][0]
        a2 = axes[0][1]
        a3 = axes[0][2]
        a4 = axes[1][0]
        a5 = axes[1][1]
        a6 = axes[1][2]

        #

        a1.set_xlabel('a ($\\AA$)')
        a1.set_ylabel('Count')

        self.make_histogram(
            axis = a1,
            values = unit_cell_parameters[0],
            )

        #

        a2.set_xlabel('b ($\\AA$)')
        a2.set_ylabel('Count')

        self.make_histogram(
            axis = a2,
            values = unit_cell_parameters[1],
            )

        #

        a3.set_xlabel('c ($\\AA$)')
        a3.set_ylabel('Count')

        self.make_histogram(
            axis = a3,
            values = unit_cell_parameters[2],
            )

        #

        a4.set_xlabel('$alpha ($\\AA$)')
        a4.set_ylabel('Count')

        self.make_histogram(
            axis = a4,
            values = unit_cell_parameters[3],
            )

        #

        a5.set_xlabel('beta ($\\AA$)')
        a5.set_ylabel('Count')

        self.make_histogram(
            axis = a5,
            values = unit_cell_parameters[4],
            )

        #

        a6.set_xlabel('gamma ($\\AA$)')
        a6.set_ylabel('Count')

        self.make_histogram(
            axis = a6,
            values = unit_cell_parameters[5],
            )

        #

        fig.subplots_adjust(hspace=0.3)
    
        return fig

    def make_histogram(self,
        axis, 
        values,
        n_bins = 30,
        ):

        axis.hist(x=values, bins=n_bins, density=False)


###


class MakePanddaDatasetSummaryGraphs(object):

    def __init__(self,
        output_dir,
        dataset_dir = None,
        dataset_subdir = "dataset_graphs",
        processor = None,
        ):

        if (processor is None):
            from giant.processors import basic_processor
            processor = basic_processor
            
        if (dataset_dir is None):
            dataset_dir = output_dir

        dataset_dir = (dataset_dir / dataset_subdir)

        self.plotters = [
            WilsonPlotter(
                output_path_template = str(
                    dataset_dir / "{label}-wilson.png"
                    ),
                ),
            ResolutionDistributionPlotter(
                output_path = str(
                    output_dir / "dataset_resolutions.png"
                    ),
                ),
            RFactorDistributionPlotter(
                output_path = str(
                    output_dir / "dataset_rfactors.png"
                    ),
                ),
            UnitCellDistributionPlotter(
                output_path = str(
                    output_dir / "dataset_unitcells.png"
                    ),
                ),
            ]

        self.processor = processor

    def __call__(self, 
        datasets, 
        dataset_dicts,
        data_getters,
        reference_miller_array,
        ):
        
        output_files = {}

        for plotter in self.plotters:

            of = plotter(
                datasets = datasets, 
                dataset_dicts = dataset_dicts,
                data_getters = data_getters,
                reference_miller_array = reference_miller_array,
                )

            merge_dicts(
                master_dict = output_files,
                merge_dict = of,
                )

        return output_files


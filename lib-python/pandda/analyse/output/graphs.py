import giant.logs as lg
logger = lg.getLogger(__name__)

import collections 
import numpy as np
import pathlib as pl

from pandda.graphs import (
    PanddaPlotter,
    )

from pandda.utils import (
    show_dict, 
    merge_dicts,
    )


class EventSitePlotter(PanddaPlotter):

    output_key = "all_events"

    def __init__(self, output_path_template):

        self.output_path_template = str(output_path_template)

        assert ('{site_num}' in output_path_template)

    def __call__(self, event_dicts, *args, **kwargs):

        event_data = self.unpack_events(event_dicts)

        output_files = collections.OrderedDict()

        for site_num, e_data in event_data.items():

            fig = self.plot(
                site_num = site_num, 
                event_data = e_data, 
                *args, **kwargs
                )

            filename = self.get_path(
                site_num = site_num,
                )

            self.save(
                fig = fig, 
                filename = filename,
                )

            output_files[site_num] = filename

        return {self.output_key : output_files}

    def plot(self,
        site_num,
        event_data,
        *args, **kwargs
        ):

        ###

        fig, axis = self.setup()

        axis.set_title(
            'Events for Site {site_num}'.format(
                site_num = site_num,
                )
            )
        axis.set_xlabel('Event #')
        axis.set_ylabel('Z-Peak Value')

        self.make_bar_plot(
            axis = axis,
            plot_vals = event_data,
            )

        return fig

    def make_bar_plot(self, axis, plot_vals, colour_bool=None, colour_vals=None, min_x=10):
        """Plot set of bar graphs in one figure"""

        assert [colour_bool,colour_vals].count(None) >= 1, 'Provide either colour_bool or colour_vals'

        n = len(plot_vals)

        assert n > 0

        if colour_bool:
            assert len(colour_bool) == n
            bar_colours = ['limegreen' if b else 'red' for b in colour_bool]

        elif colour_vals:
            assert len(colour_vals) == n
            bar_colours = colour_vals

        else: 
            bar_colours = ['slategray'] * n

        axis.bar(
            x = np.arange(n) + 1.0,
            height = plot_vals, 
            width = 0.8, 
            color = bar_colours,
            )

        axis.set_xticks(
            range(1, n+1)
            )

        axis.set_yticks(
            range(0, int(max(plot_vals)+0.5))
            )

        axis.set_xlim(
            [0.5, max(min_x, n) + 2]
            )

    def get_path(self, **kwargs):

        p = pl.Path(
            self.output_path_template.format(
                **kwargs
                )
            )
        
        if not p.parent.exists():
            p.parent.mkdir(parents=True)

        return str(p)

    def unpack_events(self, events):

        site_nums = sorted(
            [e['site_num'] for e in events]
            )

        site_events = {
            i : [
                e['z_peak'] for e in events 
                if (e['site_num'] == i)
                ]
            for i in site_nums
            }

        return site_events
        

class EventResolutionsPlotter(PanddaPlotter):

    output_key = "event_resolutions"

    def plot(self, event_dicts, analysed_resolution, *args, **kwargs):

        event_resolutions = [
            analysed_resolution[e['dtag']]
            for e in event_dicts
            ]

        ###

        fig, axis = self.setup()

        axis.set_title('Event Resolutions')
        axis.set_xlabel('Resolution ($\AA$)')
        axis.set_ylabel('Count')

        self.make_histogram(
            axis = axis,
            values = event_resolutions,
            )

        return fig

    def make_histogram(self,
        axis, 
        values,
        n_bins = 30,
        ):

        axis.hist(x=values, bins=n_bins, density=False)


class EventFractionsPlotter(PanddaPlotter):

    output_key = "event_fractions"

    def plot(self, event_dicts, *args, **kwargs):

        event_fractions = [
            e['event_fraction']
            for e in event_dicts
            ]

        ###

        fig, axis = self.setup()

        axis.set_title('Event Fractions (approx. occupancies)')
        axis.set_xlabel('Fraction')
        axis.set_ylabel('Count')

        self.make_histogram(
            axis = axis,
            values = event_fractions,
            )

        return fig

    def make_histogram(self,
        axis, 
        values,
        n_bins = 30,
        ):

        axis.hist(x=values, bins=n_bins, density=False)


class AnalysedResolutionDistributionPlotter(PanddaPlotter):

    output_key = "analysed_resolution"

    def plot(self,
        analysed_resolution,
        *args, **kwargs
        ):

        x, y = self.unpack(analysed_resolution)

        fig, axis = self.setup()

        axis.set_title('Analysed Resolutions')
        axis.set_xlabel('Resolution ($\AA$)')
        axis.set_ylabel('Count')

        axis.bar(
            x = x,
            height = y,
            width = 0.1,
            )

        axis.set_xlim(
            [0, max(x)+0.2]
            )

        return fig

    def unpack(self, 
        dkey_dict,
        ):

        x, y = np.unique(
            dkey_dict.values(), 
            return_counts = True,
            )

        return (x, y)


class MapUncertaintyDistributionPlotter(PanddaPlotter):

    output_key = "map_uncertainties"

    def plot(self,
        map_uncertainty,
        *args, **kwargs
        ):

        uncertainty_values = map_uncertainty.values()

        ###

        fig, axis = self.setup()

        axis.set_title('Map Uncertainties')
        axis.set_xlabel('Uncertainty (rmsd)')
        axis.set_ylabel('Count')

        self.make_histogram(
            axis = axis,
            values = uncertainty_values,
            )

        return fig

    def make_histogram(self,
        axis, 
        values,
        n_bins = 30,
        ):

        axis.hist(x=values, bins=n_bins, density=False)


###


class MakePanddaResultsGraphs:
    
    def __init__(self, 
        output_directory, 
        ):

        self.plotters = [
            EventSitePlotter(
                output_path_template = str(
                    output_directory / "analyse_events_site_{site_num}.png"
                    ),
                ),
            EventResolutionsPlotter(
                output_path = str(
                    output_directory / "event_resolutions.png"
                    ),
                ),
            EventFractionsPlotter(
                output_path = str(
                    output_directory / "event_fractions.png"
                    ),
                ),
            AnalysedResolutionDistributionPlotter(
                output_path = str(
                    output_directory / "analysed_resolutions.png"
                    ),
                ),
            MapUncertaintyDistributionPlotter(
                output_path = str(
                    output_directory / "map_uncertainties.png"
                    ),
                ),
            ]

    def __call__(self, 
        event_dicts,
        shell_dicts,
        dataset_dicts,
        ):

        output_files = {}

        for plotter in self.plotters:

            of = plotter(
                event_dicts = event_dicts,
                shell_dicts = shell_dicts,
                **dataset_dicts # provide these as columns rather than the dict itself
                )

            merge_dicts(
                master_dict = output_files,
                merge_dict = of,
                )

        return output_files


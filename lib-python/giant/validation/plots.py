import giant.logs as lg 
logger = lg.getLogger(__name__)

import copy, itertools

import numpy as np
import pandas as pd

from matplotlib import pyplot as plt

from giant.utils import (
    merge_dicts,
    )

from giant.plot import (
    Radar,
    )


class ValidationRadarPlot(object):

    PLOT_DEFAULTS = {
        'plot_order' : [
            'rscc',
            'rszd',
            'rszo',
            'b_factor_ratio',
            'rmsd',
            ],
        'input_key_hash' : {
            'rscc' : 'RSCC',
            'rszd' : 'RSZD',
            'rszo' : 'RSZO/OCC',
            'b_factor_ratio' : 'B-factor Ratio',
            'rmsd' : 'Model RSMD',
        },
        'axis_params' : {
            'rscc' : {
                'title' : '\n\nModel\nQuality\n(RSCC)',
                'axis_min' : 0.60,
                'axis_max' : 0.85,
                'axis_invert' : True,
            },
            'rszd' : {
                'title' : 'Model\nAccuracy\n(RSZD)',
                'axis_min' : 1.50,
                'axis_max' : 4.00,
                'axis_invert' : False,
            },
            'rszo' : {
                'title' : 'Density\nPrecision\n(RSZO/OCC)',
                'axis_min' : 0.00,
                'axis_max' : 2.00,
                'axis_invert' : True,
            },
            'b_factor_ratio' : {
                'title' : 'B-Factor\nStability\n(B-factor Ratio)',
                'axis_min' : 1.00,
                'axis_max' : 3.00,
                'axis_invert' : False,
            },
            'rmsd' : {
                'title' : 'Coordinate\nStability\n(RMSD)',
                'axis_min' : 0.00,
                'axis_max' : 1.50,
                'axis_invert' : False,
            },
        },
        'colours' : ['r', 'g', 'b', 'y'],
        'markers' : ['o','^','s','D','*','+'],
        'linestyle' : ['-','--'],
        'markersize' : 5,
    }

    def __init__(self,
        plot_defaults = None,
        simplified_plot_labels = False,
        remove_blank_entries = True,
        ):

        self.plot_parameters = copy.deepcopy(
            self.PLOT_DEFAULTS
            )

        if plot_defaults is not None:

            merge_dicts(
                master_dict = self.plot_parameters,
                merge_dict = plot_defaults,
                dict_class = dict,
                overwrite = True,
                )

        self.simplified_plot_labels = bool(
            simplified_plot_labels
            )

        self.remove_blank_entries = bool(
            remove_blank_entries
            )

    def __call__(self,
        values_dicts,
        out_path,
        ):

        from giant.plot import Radar

        plot_df = self.get_plot_df(
            values_dicts = values_dicts,
            )

        r = self.make_plot(
            plot_df = plot_df,
            )

        r.plot()
        r.legend()
        r.axis_limits()
        r.savefig(str(out_path))
        r.close()

    def get_plot_df(self, values_dicts):

        # Internal naming
        column_keys = self.plot_parameters['plot_order']
        # External/input naming (optional)
        column_names = [
            self.plot_parameters['input_key_hash'].get(k,k)
            for k in column_keys
        ]

        row_labels = []
        row_values = []

        for i, v_dict in enumerate(values_dicts):

            row_labels.append(
                v_dict.get('label', i+1)
                )

            row_values.append([
                v_dict.get(k)
                for k in column_names
                ])

        df = pd.DataFrame(
            data = row_values,
            index = row_labels,
            columns = column_keys, # use internal names
            )

        if (self.remove_blank_entries is True):

            df = df.dropna(axis=1, how='all')

        return df

    def make_plot(self, plot_df):

        pp = copy.deepcopy(
            self.plot_parameters
            )

        column_keys = list(plot_df.columns)

        axis_params = [
            pp['axis_params'][k]
            for k in column_keys
            ]

        colours = itertools.cycle(pp['colours'])
        markers = itertools.cycle(pp['markers'])
        linestyle = itertools.cycle(pp['linestyle'])
        markersize = pp['markersize']

        r = Radar(
            titles = [p['title'] for p in axis_params],
            )

        for row_label, row_values in plot_df.iterrows():

            r.add(
                list(row_values),
                color = next(colours),
                linewidth = 2,
                linestyle = next(linestyle),
                marker = next(markers),
                markersize = markersize,
                markeredgecolor = 'k',
                label = str(row_label),
            )

        r.set_inversion([
            p['axis_invert']
            for p in axis_params
            ])

        r.set_limits([
            (p['axis_min'], p['axis_max'])
            for p in axis_params
            ])

        ###

        if (self.simplified_plot_labels is False) or len(plot_df.index) == 1:

            r.set_ticks(
                labels = [
                    list(map(self.label,plot_df[c].values))
                    for c in column_keys
                    ],
                values = [
                    list(map(self.value,plot_df[c].values))
                    for c in column_keys
                    ],
                )

        return r

    @staticmethod
    def label(value):
        try:
            label = round(value, 2)
        except Exception as e:
            return 'n/a'
        return label

    @staticmethod
    def value(value):
        try:
            label = round(value, 2)
        except Exception as e:
            return None
        return value


class ValidationDistributionPlots(object):

    PLOT_DEFAULTS = {}

    def __init__(self,
        plot_defaults = None,
        simplified_plot_labels = False,
        remove_blank_entries = True,
        ):

        self.plot_parameters = copy.deepcopy(
            self.PLOT_DEFAULTS
            )

        if plot_defaults is not None:

            merge_dicts(
                master_dict = self.plot_parameters,
                merge_dict = plot_defaults,
                dict_class = dict,
                overwrite = True,
                )

    def __call__(self,
        scores_df,
        out_path,
        split_on_column = None,
        ):

        if split_on_column is not None:
            cluster_dfs = self.split_df_on_column(
                scores_df,
                split_on_column = split_on_column,
                )
        else:
            cluster_dfs = [scores_df]

        fig, axes = self.make_plot(
            df_list = cluster_dfs,
            )

        fig.tight_layout()
        fig.savefig(str(out_path))
        plt.close(fig)

        return out_path

    def split_df_on_column(self, df, split_on_column):

        split_values = df[split_on_column]
        split_values_list = list(split_values)

        column_vals = sorted(
            set(split_values),
            key = lambda v: split_values_list.index(v),
            )

        out_dfs = []

        for v in column_vals:

            sel_df = df[split_values == v]

            out_dfs.append(
                sel_df.drop(split_on_column, axis=1)
                )

        return out_dfs

    def make_plot(self, df_list):

        cat_columns = np.concatenate([d.columns for d in df_list]).tolist()

        srt_columns = sorted(
            set(cat_columns),
            key = lambda c: cat_columns.index(c),
            )

        fig, axes = plt.subplots(
            nrows = len(df_list),
            ncols = len(srt_columns),
            squeeze = False,
            )

        for i_df, df in enumerate(df_list):

            for i_col, col in enumerate(srt_columns):

                axis = axes[i_df, i_col]

                if col not in df.columns:
                    continue

                try:
                    col_values = list(map(float, df[col].values))
                except:
                    # invalid column
                    continue

                try:
                    axis.hist(
                        x = col_values,
                        bins = 25,
                        density = False,
                        histtype="stepfilled",
                        alpha = 0.8,
                        )
                except ValueError as e:
                    logger.warning(str(e))
                    continue

        return fig, axes


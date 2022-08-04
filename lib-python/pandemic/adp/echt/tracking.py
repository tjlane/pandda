import os, collections
import numpy, pandas

from libtbx import adopt_init_args, group_args


class EchtTracking(object):

    csv_name = 'tracking_echt.csv'

    def __init__(self,
        output_directory,
        plotting_object,
        model_object,
        ):
        # Create table for tracking progress over cycles
        table = pandas.DataFrame(
            columns=['cycle', 'type', 'average'] + model_object.dataset_labels,
            )

        self.output_files = collections.OrderedDict(
            tracking_csv = os.path.join(output_directory, self.csv_name),
        )

        write_graphs = EchtTrackingPlotter(
            parent = self,
            plotting_object = plotting_object,
            output_directory = output_directory,
            )

        adopt_init_args(self, locals())

    def update(self,
        model_object,
        n_cycle,
        ):

        amplitudes = self.extract_amplitudes(model_object)

        amplitudes_sum = amplitudes.values.sum(axis=0)
        amplitudes_squared = amplitudes ** 2
        amplitudes_squared_sum = amplitudes_squared.sum(axis=0)

        self.table.loc[len(self.table)] = [n_cycle, 'sum of amplitudes',     amplitudes_sum.mean()] + list(amplitudes_sum)
        self.table.loc[len(self.table)] = [n_cycle, 'sum of (amplitudes^2)', amplitudes_squared_sum.mean()] + list(amplitudes_squared_sum)

    def write_output(self):
        self.table.to_csv(self.output_files['tracking_csv'])
        of = self.write_graphs()
        self.output_files.update(of)

    def extract_amplitudes(self,
        model_object,
        ):

        amplitudes = pandas.DataFrame(columns=['level_group']+model_object.dataset_labels)

        for i_l, level_name in enumerate(model_object.all_level_names):

            if level_name not in model_object.tls_level_names:
                continue

            i_l_tls = model_object.tls_level_names.index(level_name)

            tls_objects = model_object.tls_objects[i_l_tls]

            for i_g, tlso in enumerate(tls_objects):
                dataset_amplitudes = numpy.mean([o.amplitudes.values for o in tlso.tls_parameters], axis=0)
                amplitudes.loc[len(amplitudes)] = [(i_l, i_g)] + dataset_amplitudes.tolist()

        amplitudes = amplitudes.set_index('level_group')

        return amplitudes

    def as_html_summary(self):
        from pandemic.adp.echt.html.tracking import EchtTrackingHtmlSummary
        return EchtTrackingHtmlSummary(self)


class EchtTrackingPlotter(object):

    _amplitudes_png = 'tracking_amplitudes.png'

    def __init__(self,
        parent,
        plotting_object,
        output_directory,
        ):

        output_files = collections.OrderedDict(
            amplitudes_lineplot = os.path.join(output_directory, self._amplitudes_png),
            )

        adopt_init_args(self, locals())

    def __call__(self):

        # Extract variables from parent
        self.amplitudes_lineplot(
            table = self.parent.table,
            filename = self.output_files['amplitudes_lineplot'],
            )

        return self.output_files

    def amplitudes_lineplot(self,
        table,
        filename,
        ):

        # Useful functions
        helper = self.plotting_object.helper

        weight_labels = sorted(set(table['type']))

        x_vals_array = []
        y_vals_array = []
        for label in weight_labels:
            l_table = table[table['type']==label]
            x_vals_array.append(l_table['cycle'].values)
            y_vals_array.append(l_table['average'].values)
        nx = max(map(len,x_vals_array))

        x_ticks = list(map(int,sorted(set(table['cycle'].values))))

        self.plotting_object.lineplot(
            x_vals_array = x_vals_array,
            y_vals_array = y_vals_array,
            title = 'Model parameters/penalties across cycles',
            x_label = 'Cycle',
            y_label = 'Amplitudes ($\\AA^2$ or $\\AA^4$)',
            x_ticks = x_ticks,
            legends = weight_labels,
            filename = filename,
            background_line_type = 'chunky',
            )

import os
import numpy, pandas

from libtbx import adopt_init_args, group_args

class EchtTracking:

    def __init__(self,
        output_directory,
        plotting_object,
        model_object,
        verbose = False,
        log = None,
        ):
        if log is None: log = Log()
        # Create table for tracking progress over cycles
        table = pandas.DataFrame(
            columns=['cycle', 'type', 'average'] + model_object.dataset_labels,
            )
        tracking_csv = os.path.join(output_directory, 'tracking_echt.csv')
        adopt_init_args(self, locals())

    def update(self,
        model_object,
        i_cycle,
        ):

        log = self.log
        log.subheading('Updating ECHT tracking...')

        amplitudes = self.extract_amplitudes(model_object)

        amplitudes_sum = amplitudes.values.sum(axis=0)
        amplitudes_squared = amplitudes ** 2
        amplitudes_squared_sum = amplitudes_squared.sum(axis=0)

        self.table.loc[len(self.table)] = [i_cycle, 'sum of amplitudes',     amplitudes_sum.mean()] + list(amplitudes_sum)
        self.table.loc[len(self.table)] = [i_cycle, 'sum of (amplitudes^2)', amplitudes_squared_sum.mean()] + list(amplitudes_squared_sum)

        self.table.to_csv(self.tracking_csv)

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
        from pandemic.adp.echt.html import EchtTrackingHtmlSummary
        return EchtTrackingHtmlSummary(self)

import giant.logs as lg
logger = lg.getLogger(__name__)

import os, copy, collections
import numpy, pandas

from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure
from giant.structure.uij import uij_to_b

from pandemic.adp import constants
from pandemic.adp.tracking.uijs import UijHistory
from pandemic.adp.tracking.convergence import PandemicConvergenceChecker
from pandemic.adp.tracking import plots
from pandemic.functions import rms


class PandemicTrackingObject(object):

    _csv_name1 = 'tracking_levels.csv'
    _csv_name2 = 'tracking_rmsds.csv'

    _ConvergenceCheckerClass = PandemicConvergenceChecker

    def __init__(self,
            output_directory,
            plotting_object,
            structure_factory,
            dataset_names,
            level_names,
            convergence_args = {},
            ):

        # Create table for tracking progress over cycles
        table = pandas.DataFrame(
            columns=[
                'cycle', 'step', 'level', 'level#', 'rmsd',
                'b_min', 'b_avg', 'b_max', 'b_avg (total)',
                ],
            )
        # Table for tracking dataset rmsds
        table_by_dataset = pandas.DataFrame(
            columns=['cycle', 'type', 'overall'] + dataset_names,
            )

        output_files = collections.OrderedDict(
            tracking_csv1 = os.path.join(output_directory, self._csv_name1),
            tracking_csv2 = os.path.join(output_directory, self._csv_name2),
            )

        n_cycle = 0

        # Create plotting object from provided plotting object
        write_graphs = plots.PandemicTrackingPlotter(
            parent = self,
            output_directory = output_directory,
            plotting_object = plotting_object,
            )

        # History of the Uijs for optimisation
        uij_history = UijHistory()

        # Class to check when the model has converged
        convergence_checker = self._ConvergenceCheckerClass(
            parent = self,
            **convergence_args
            )

        adopt_init_args(self, locals(), exclude=('convergence_args',))

    def set_target(self, uij_target, uij_target_weights):
        self.uij_target                  = uij_target
        self.uij_target_weights          = uij_target_weights
        self.uij_target_weights_expanded = uij_target_weights.reshape(uij_target_weights.shape+(1,)).repeat(6, axis=-1)
        self.uij_target_weights_average  = uij_target_weights.mean(axis=0)  # (DATASET, atoms) -> (atoms, )

    def update(self,
        uijs,
        step,
        i_level=None,
        ):
        """Update the tracking table"""

        logger.subheading('Updating tracking...')

        # Must be list, if only one given
        if not isinstance(i_level, list):
            i_level = [i_level]

        self._update_average_table(uijs=uijs, step=step, i_level=i_level)

        logger(self.table.loc[len(self.table)-len(i_level):].to_string())

        if (step == 'end'):

            # Update the dataset-by-dataset (rmsd) statistics
            self._update_dataset_table(uijs=uijs)

            # Update history
            self.uij_history.add(n_cycle=self.n_cycle, uijs=uijs)

            # Update convergence data
            self.convergence_checker.update()

        return

    def _update_average_table(self,
        uijs,
        step,
        i_level,
        ):

        # Sum along various axes
        uij_dst_lvl = uijs                      # (LEVEL, DATASET, atoms, 6)
        uij_dst_tot = uij_dst_lvl.sum(axis=0)   # (x, DATASET, atoms, 6)
        # Average over all datasets
        uij_avg_lvl = uij_dst_lvl.mean(axis=1)  # (LEVEL, x, atoms, 6)
        uij_avg_tot = uij_avg_lvl.sum(axis=0)   # (atoms, 6)

        # Calculate the rms between fitted and input
        b_rmsd = constants.EIGHTPISQ * rms(
            values = self.uij_target - uij_dst_tot,
            weights = self.uij_target_weights_expanded,
            )
        # Calculate B-iso for complete model
        b_iso_tot = numpy.mean(uij_to_b(uij_avg_tot))

        # Iterate through levels to be dumped into table
        for i_l in i_level:

            # Defaults to be overridden
            uij_sel = None
            level_number = None
            level_name = None
            b_iso_sel = 0.0
            b_min_sel = numpy.nan
            b_max_sel = numpy.nan

            # Extract the Uij for the selected level(s)
            if isinstance(i_l, int):
                # Extract level information
                level_number = i_l+1
                level_name = self.level_names[i_l]
                # Calculate min/mean/max B-iso for selected level
                b_iso = uij_to_b(uij_avg_lvl[i_l])
                b_avg = numpy.mean(b_iso)
                b_min = numpy.min(b_iso)
                b_max = numpy.max(b_iso)

            # Add to tracking table
            self.table.loc[len(self.table.index)] = {
                'cycle'  : self.n_cycle,
                'step'   : step,
                'level#' : level_number,
                'level'  : level_name,
                'rmsd'   : round(b_rmsd,6),
                'b_avg' : round(b_avg,3),
                'b_min' : round(b_min,3),
                'b_max' : round(b_max,3),
                'b_avg (total)' : round(b_iso_tot,3),
                }

    def _update_dataset_table(self,
        uijs,
        ):

        uijs_dst = uijs.sum(axis=0) # (x, DATASET, atoms, 6)

        total_rmsd = constants.EIGHTPISQ * rms(
            values = self.uij_target - uijs_dst,
            weights = self.uij_target_weights_expanded,
            )

        # Store by-dataset RMSDs
        dataset_rmsds = [constants.EIGHTPISQ * rms(
            values = dst_diffs,
            weights = self.uij_target_weights_expanded[i_dst],
            ) for i_dst, dst_diffs in enumerate(self.uij_target-uijs_dst)
        ]
        self.table_by_dataset.loc[len(self.table_by_dataset)] = [self.n_cycle, 'rmsd', total_rmsd] + list(dataset_rmsds)

    def write_output(self):
        self.table.to_csv(self.output_files['tracking_csv1'])
        self.table_by_dataset.to_csv(self.output_files['tracking_csv2'])
        of = self.write_graphs()
        self.output_files.update(of)

    def is_converged(self):
        logger.subheading('Checking convergence')
        self.convergence_checker.show()
        cvgd = self.convergence_checker.is_converged()
        return cvgd

    def as_html_summary(self):
        from pandemic.adp.tracking.html import TrackingHtmlSummary
        return TrackingHtmlSummary(self)

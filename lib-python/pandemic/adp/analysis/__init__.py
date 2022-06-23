import giant.logs as lg
logger = lg.getLogger(__name__)

import os, collections

from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure

from scitbx.array_family import flex

import pathlib as pl

class HierarchicalModelAnalysisTask(object):


    def __init__(self,
        output_directory,
        plotting_object,
        master_phil,
        analyse_residuals = True,
        assess_hierarchy_groups = True,
        ):

        output_directory = pl.Path(output_directory)

        if analyse_residuals is True:
            from pandemic.adp.analysis.residuals import AnalyseResidualsTask
            analyse_residuals = AnalyseResidualsTask(
                output_directory = (output_directory / 'residuals'),
                plotting_object = plotting_object,
                )

        if assess_hierarchy_groups is True:
            from pandemic.adp.analysis.hierarchy import AssessHierarchyGroupsTask
            assess_hierarchy_groups = AssessHierarchyGroupsTask(
                output_directory = (output_directory / 'hierarchy_groups'),
                plotting_object = plotting_object,
                )

        # Calculate fit to electron density - TODO
        #calculate_model_density_fit = CalculateModelDensityFitTask()

        adopt_init_args(self, locals())

    def run(self,
        uij_target,
        uij_target_weights,
        uij_isotropic_mask,
        model_object,
        model_hierarchy_info,
        model_hierarchy_files,
        reference_hierarchy,
        ):

        logger.heading('Analysing Fitted Hierarchical Model')

        dataset_labels = model_object.dataset_labels
        level_names = model_object.all_level_names
        uij_fitted = model_object.uijs()

        from pandemic.adp.hierarchy.custom_levels import MakeNewCustomHierarchyEffFilesFromIndices
        write_levels_function = MakeNewCustomHierarchyEffFilesFromIndices.from_model_object(
            model_object = model_object,
            master_phil = self.master_phil,
            custom_level_scope_name = 'model.custom_level',
            )

        # Convert to flex
        overall_atom_selection = flex.bool(model_hierarchy_info.overall_atom_mask)

        if self.analyse_residuals is not False:
            logger.subheading('Analysing Model-Target Residuals', spacer=True)
            residuals_out = self.analyse_residuals.run(
                uij_fitted = uij_isotropic_mask(uij_fitted),
                uij_target = uij_target,
                uij_target_weights = uij_target_weights,
                level_names = level_names,
                dataset_labels = dataset_labels,
                reference_hierarchy = reference_hierarchy.select(overall_atom_selection, copy_atoms=True),
                )
        else:
            residuals_out = None

        if self.assess_hierarchy_groups is not False:
            logger.subheading('Assessing Hierarchical Group Partitions', spacer=True)
            assess_hierarchy_out = self.assess_hierarchy_groups.run(
                model_object = model_object,
                model_hierarchy_info = model_hierarchy_info,
                model_hierarchy_files = model_hierarchy_files,
                reference_hierarchy = reference_hierarchy,
                overall_atom_selection = overall_atom_selection,
                write_levels_function = write_levels_function,
                )
        else:
            assess_hierarchy_out = None

        self.result = group_args(
            residuals = residuals_out,
            assess_hierarchy = assess_hierarchy_out,
            )

        return self.result

    def as_html_summary(self):

        from pandemic.adp.html import (
            HtmlSummaryCollator,
            as_html_summaries_maybe,
            )

        return HtmlSummaryCollator(
            title = 'Analysis of Hierarchical model ADPs',
            alt_title = 'Fitting Analysis',
            summaries = as_html_summaries_maybe([
                self.analyse_residuals,
                self.assess_hierarchy_groups,
                ]),
            )

import giant.logs as lg
logger = lg.getLogger(__name__)

import os, copy, collections
import numpy

from libtbx import adopt_init_args, group_args

from giant.paths import rel_symlink, easy_directory

from pandemic.adp.utils import show_file_dict

from pandemic.adp.postprocess.refine import RefineStructures
from pandemic.adp.postprocess.table_one import CalculateTableOnes


class PostProcessTask(object):


    RefineStructuresClass = RefineStructures
    CalculateTableOneClass = CalculateTableOnes

    refine_structures = None
    calculate_r_factors = None

    input_suffix = ' (Input)'
    fitted_suffix = ' (Fitted)'
    refined_suffix = ' (Refined)'
    difference_suffix = ' Diff.'

    # columns to be extracted from each file
    reflection_data_columns = collections.OrderedDict([
        ('high_res',  'High Resolution Limit'),
        ('low_res',   'Low Resolution Limit'),
        ('uniq_refl', 'Unique reflections'),
        ('complete',  'Completeness (%)'),
        ('wilson_b',  'Wilson B-factor'),
        ])

    model_columns = collections.OrderedDict([
        ('r_free', 'R-free'),
        ('r_work', 'R-work'),
        ('r_gap',  'R-gap'),
        ('av_b',   'Average B-factor'),
        ])

    r_factor_keys = ['r_free', 'r_work', 'r_gap']

    show_file_dict = show_file_dict

    def __init__(self,
        output_directory,
        structure_directory,
        refine_structures = True,
        calculate_r_factors = True,
        refinement_program = 'phenix',
        table_one_options = None,
        plotting_object = None,
        n_cpus = 1,
        ):

        if refine_structures is True:
            refine_structures = self.RefineStructuresClass(
                output_directory = structure_directory,
                refinement_program = refinement_program,
                n_cpus = n_cpus,
                )

        if calculate_r_factors is True:
            calculate_r_factors = self.CalculateTableOneClass(
                    output_directory = easy_directory(os.path.join(output_directory, 'table_ones')),
                    table_one_options = table_one_options,
                    n_cpus = n_cpus,
                    )

        adopt_init_args(self, locals())

    def run(self,
        dataset_labels,
        input_structures_dict,
        input_reflection_data_dict,
        fitted_structures_dict,
        cif_files = None,
        results_object = None,
        ):

        # Extract ordered variables from dicts
        input_structures      = [input_structures_dict[l] for l in dataset_labels]
        input_reflection_data = [input_reflection_data_dict[l] for l in dataset_labels]
        fitted_structures     = [fitted_structures_dict[l] for l in dataset_labels]

        # Do we have reference r-factors?
        has_reference_data = (results_object.reference_r_values is not None)

        # Extract variables to local
        input_suffix = self.input_suffix
        fitted_suffix = self.fitted_suffix
        refined_suffix = self.refined_suffix

        if has_reference_data:
            # Extract suffix for comparisons
            reference_suffix = results_object.reference_suffix
            # SMALL HACK -- add dummy columns to output table if they do not exist to allow comparison!
            for col_base in self.model_columns.values():
                col = (col_base+reference_suffix)
                if col not in results_object.table.columns:
                    results_object.table.loc[:,col] = None

        difference_suffix_template = ' {} '.format(self.difference_suffix.strip()) + '({}-{})'
        # Characters that are removed when formatting difference template
        strip_chars = ' ()[]{}'

        if (self.refine_structures is not False) or (self.calculate_r_factors is not False):

            # Needed because the table one calculation assumes same base name for pdbs and mtzs
            self.link_reflection_data_to_fitted_models(
                input_reflection_data = input_reflection_data,
                fitted_structures = fitted_structures,
                )

        ###########################################

        if self.refine_structures is not False:

            logger.subheading('Refining fitted structures')

            refined_structures_dict = self.refine_structures(
                    input_structures = input_structures,
                    input_reflection_data = input_reflection_data,
                    labels = dataset_labels,
                    cif_files = cif_files,
                    output_suffix = '-refined',
                    )
            # ORDER OF THESE SHOULDN'T BE NEEDED AS SOME REFINEMENTS MIGHT FAIL
            refined_structures_unsorted = [pdb for (pdb, mtz) in refined_structures_dict.values()]
        else:
            refined_structures_unsorted = None

        ###########################################

        # Labels for sets of R-factors
        all_suffixes = []
        if has_reference_data: all_suffixes.append(reference_suffix)

        # Sets for r-factor differences
        comparison_triplets = []

        if self.calculate_r_factors is not False:

            logger.subheading('Running phenix.table_one to calculate R-factors')

            # Input Structures
            inp_csv = self.calculate_r_factors(
                output_prefix = 'table_one_input',
                structures = input_structures,
                )
            # Extract data columns (only needed once)
            results_object.add_table_one(
                csv_file = inp_csv,
                column_labels = list(self.reflection_data_columns.values()),
                output_suffix = None,
                )
            # Extract model columns
            results_object.add_table_one(
                csv_file = inp_csv,
                column_labels = list(self.model_columns.values()),
                output_suffix = self.input_suffix,
                )
            all_suffixes.append(self.input_suffix)

            # Fitted Structures
            fit_csv = self.calculate_r_factors(
                output_prefix = 'table_one_fitted',
                structures = fitted_structures,
                )
            # Extract model columns
            results_object.add_table_one(
                csv_file = fit_csv,
                column_labels = list(self.model_columns.values()),
                output_suffix = self.fitted_suffix,
                )
            all_suffixes.append(self.fitted_suffix)

            # Create suffix for differences
            fitted_input_suffix = difference_suffix_template.format(
                fitted_suffix.strip(strip_chars),
                input_suffix.strip(strip_chars),
                )
            # Calculate differences between columns
            results_object.calculate_column_differences(
                column_labels = list(self.model_columns.values()),
                input_suffixes = [fitted_suffix, input_suffix],
                output_suffix = fitted_input_suffix,
                )
            comparison_triplets.append((fitted_suffix, input_suffix, fitted_input_suffix))

            # Also do comparisons to reference values
            if has_reference_data:

                fitted_reference_suffix = difference_suffix_template.format(
                    fitted_suffix.strip(strip_chars),
                    reference_suffix.strip(strip_chars),
                    )
                results_object.calculate_column_differences(
                    column_labels = list(self.model_columns.values()),
                    input_suffixes = [fitted_suffix, reference_suffix],
                    output_suffix = fitted_reference_suffix,
                    )
                comparison_triplets.append((fitted_suffix, reference_suffix, fitted_reference_suffix))

        if (self.calculate_r_factors is not False) and (refined_structures_unsorted is not None):

            # Refined Structures
            ref_csv = self.calculate_r_factors(
                output_prefix = 'table_one_refined',
                structures = refined_structures_unsorted,
                )
            # Extract model columns
            results_object.add_table_one(
                csv_file = ref_csv,
                column_labels = list(self.model_columns.values()),
                output_suffix = self.refined_suffix,
                )
            all_suffixes.append(self.refined_suffix)

            # Create suffix for differences
            refined_input_suffix = difference_suffix_template.format(
                refined_suffix.strip(strip_chars),
                input_suffix.strip(strip_chars),
                )
            # Calculate differences between columns
            results_object.calculate_column_differences(
                column_labels = list(self.model_columns.values()),
                input_suffixes = [refined_suffix, input_suffix],
                output_suffix = refined_input_suffix,
                )
            comparison_triplets.append((refined_suffix, input_suffix, refined_input_suffix))

            # Also do comparisons to reference values
            if has_reference_data:

                refined_reference_suffix = difference_suffix_template.format(
                    refined_suffix.strip(strip_chars),
                    reference_suffix.strip(strip_chars),
                    )
                results_object.calculate_column_differences(
                    column_labels = list(self.model_columns.values()),
                    input_suffixes = [refined_suffix, reference_suffix],
                    output_suffix = refined_reference_suffix,
                    )
                comparison_triplets.append((refined_suffix, reference_suffix, refined_reference_suffix))

        if all_suffixes:
            logger('Calculated R-factors for these models: \n\t{}'.format('\n\t'.join(all_suffixes)))

        ###########################################

        output_dict = collections.OrderedDict()

        if self.calculate_r_factors is not False:

            logger.subheading('R-factor statistics for structures', spacer=True)

            results_object.show_average_statistics_summary(
                single_column_labels_dict = self.reflection_data_columns,
                comparison_column_prefixes_dict = self.model_columns,
                comparison_column_suffix_triplets = comparison_triplets,
                print_strip_chars = strip_chars,
                label_hash = {s[2]:'Difference' for s in comparison_triplets},
                )

            logger.subheading('Making R-factor summary graphs')

            for key in self.r_factor_keys:

                filename = self.make_r_factor_graph(
                    table = results_object.table,
                    filename_base = 'r_factor_model_comparisons-{}'.format(key),
                    column_prefix = self.model_columns[key],
                    column_suffixes = all_suffixes,
                    column_prefix_label = self.model_columns[key],
                    column_suffixes_labels = [s.strip(strip_chars) for s in all_suffixes],
                    title = '{} values for all models'.format(self.model_columns[key]),
                    )
                output_dict.setdefault('r_values', collections.OrderedDict())[key] = filename

            # Plot difference in R-factors for each model type
            for (col1, col2, del_col) in comparison_triplets:

                # make compatible with filenames
                col1 = col1.strip(strip_chars).replace(' ','').lower()
                col2 = col2.strip(strip_chars).replace(' ','').lower()

                filename = self.make_r_factor_graph(
                    table = results_object.table,
                    filename_base = 'r_factor_changes-{}-{}'.format(col1, col2),
                    column_prefix = '',
                    column_suffixes = [self.model_columns[k]+del_col for k in self.r_factor_keys],
                    column_prefix_label = 'R-value'+self.difference_suffix,
                    column_suffixes_labels = [self.model_columns[k] for k in self.r_factor_keys],
                    title = 'R-value differences between models:\n{}'.format(del_col),
                    horizontal_lines = [0.0],
                    )
                output_dict.setdefault('r_value_differences', collections.OrderedDict())[(col1.strip(strip_chars), col2.strip(strip_chars))] = filename

        self.show_file_dict(output_dict)

        self.result = group_args(
            output_files = output_dict,
            model_suffixes = all_suffixes,
            comparison_triplets = comparison_triplets,
            )

        return self.result

    def link_reflection_data_to_fitted_models(self,
        input_reflection_data,
        fitted_structures,
        ):

        for inp_mtz, fit_pdb in zip(
                input_reflection_data,
                fitted_structures,
                ):
            fit_mtz = fit_pdb.replace('.pdb', '.mtz')
            rel_symlink(inp_mtz, fit_mtz)

    def make_r_factor_graph(self,
        table,
        filename_base,
        column_prefix,
        column_suffixes,
        column_prefix_label = None,
        column_suffixes_labels = None,
        x_column = 'High Resolution Limit',
        x_label = 'Resolution ($\\AA$)',
        title = 'Structure R-values',
        r_factors_as_percentages = True,
        horizontal_lines = [],
        ):
        """Look at pre- and post-fitting graphs"""

        if column_prefix_label is None:
            column_prefix_label = column_prefix
        if column_suffixes_labels is None:
            column_suffixes_labels = column_suffixes

        # Cast table to fix type problems
        table = table.astype(float)

        # Extract x-values from table
        x_values = table[x_column]

        if r_factors_as_percentages is True:
            mult = 100.0
            column_prefix_label = column_prefix_label + ' (%)'
        else:
            mult = 1.0

        # Create filename
        filename = os.path.join(self.output_directory, '{}.png'.format(filename_base))
        # Extract y values
        y_values = [mult*table[column_prefix+suff] for suff in column_suffixes]
        # Plot
        self.plotting_object.binned_boxplot(
            filename = filename,
            x = x_values,
            y_vals = y_values,
            legends = column_suffixes_labels,
            title = title,
            x_lab = x_label,
            y_lab = column_prefix_label,
            rotate_x_labels = True,
            max_bins = 8,
            min_bin_width = 0.1,
            hlines = horizontal_lines,
            )

        return filename



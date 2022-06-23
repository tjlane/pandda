import giant.logs as lg
logger = lg.getLogger(__name__)

import os
import math, collections
import numpy, pandas

from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure

from giant.structure.uij import uij_to_b


class PandemicResultsObject(object):


    initial_columns = [
        'High Resolution Limit',
        'Low Resolution Limit',
        'R-free Change (Fitted-Reference)',
        'R-work Change (Fitted-Reference)',
        'R-gap Change (Fitted-Reference)',
        'R-free Change (Refined-Reference)',
        'R-work Change (Refined-Reference)',
        'R-gap Change (Refined-Reference)',
        'R-free Change (Fitted-Input)',
        'R-work Change (Fitted-Input)',
        'R-gap Change (Fitted-Input)',
        'R-free Change (Refined-Input)',
        'R-work Change (Refined-Input)',
        'R-gap Change (Refined-Input)',
        'Average B-factor Change (Fitted-Input)',
        'Average B-factor (fitted atoms) Change (Fitted-Input)',
        ]

    output_column_labels_base = {
        'r_free' : 'R-free',
        'r_work' : 'R-work',
        'r_gap'  : 'R-gap',
        }

    r_factor_keys = ['r_free', 'r_work', 'r_gap']

    def __init__(self,
            filename,
            plotting_object,
            models = [],
            ):

        reference_r_values = None
        reference_suffix = None

        adopt_init_args(self, locals())

        # Create output table for all model statistics
        self.table = pandas.DataFrame(
                index   = [m.tag for m in self.models],
                columns = self.initial_columns,
                )

    def add_reference_r_values(self,
        filename,
        input_column_labels = {
            'r_free' : 'R-free',
            'r_work' : 'R-work',
            'r_gap'  : 'R-gap',
            },
        output_suffix = ' (Reference)',
        ):

        table = self.table

        # Internal keys for dicts
        column_keys = self.r_factor_keys

        # Store suffix for reference
        self.reference_suffix = output_suffix

        # Form input and output label dictionaries
        i_labs = input_column_labels
        o_labs = {k:v+output_suffix for k,v in self.output_column_labels_base.items()}

        # Check keys
        assert not set(column_keys).difference(i_labs.keys())
        assert not set(column_keys).difference(o_labs.keys())

        # Extract strings from dict
        i_rf_s = i_labs['r_free']
        i_rw_s = i_labs['r_work']
        i_rg_s = i_labs['r_gap']

        o_rf_s = o_labs['r_free']
        o_rw_s = o_labs['r_work']
        o_rg_s = o_labs['r_gap']

        logger.subheading('Adding reference values to results table')

        logger('Reading Reference R-values from {}'.format(filename))

        # Prepare an error string
        err_str = 'A table of reference R-values has been provided ({}) '.format(filename)

        # Check has expected extension
        if not filename.endswith('.csv'):
            raise Sorry(err_str + 'but it must have a .csv extension.')

        # Read reference table
        self.reference_r_values = pandas.read_csv(filename, index_col=0)

        # Validate row names
        missing_ids = set(table.index).difference(self.reference_r_values.index)
        if len(missing_ids) == len(table.index):
            raise Sorry(err_str + 'but all of the dataset labels are missing. \nThe dataset labels must be present as the first column of the csv file.' + \
                    '\nCurrent datasets labels (should be present in first column in {}): \n\t{}'.format(filename, '\n\t'.join(table.index)) + \
                    '\nCurrent first column of {}: \n\t{}'.format(filename, '\n\t'.join(self.reference_r_values.index)))
        elif missing_ids:
            raise Sorry(err_str + 'but some of the dataset labels are missing from the index column.' + \
                    '\nMissing labels: \n\t{}'.format('\n\t'.join(missing_ids)))

        # Transfer R-values
        if (i_rg_s not in self.reference_r_values.columns):
            logger('\nCalculating R-gap column for reference R-free & R-work values')
            self.reference_r_values.loc[:,i_rg_s] = self.reference_r_values.loc[:,i_rf_s] - self.reference_r_values.loc[:,i_rw_s]
        else:
            logger('\nR-gap column already present in reference table -- using this column.')

        # Validate columns
        missing_cols = set(input_column_labels.values()).difference(self.reference_r_values.columns)
        if missing_cols:
            raise Sorry(err_str + 'but columns are missing. \nMissing Columns: {}'.format(', '.join(missing_cols)))

        for key in column_keys:

            i_l = i_labs[key]
            o_l = o_labs[key]

            if o_l in table.index:
                raise Failure('Column already exists in results table: {}'.format(o_l))

            # Extract values and allow for different ordering between tables
            new_vals = self.reference_r_values.loc[:,i_l][table.index]

            assert list(new_vals) == [self.reference_r_values.loc[i,i_l] for i in table.index]

            # Transfer to output table
            table[o_l] = new_vals

        #############################
        #           Report          #
        #############################

        logger('\n> Processed Reference Values\n')
        for i, v in table[[o_labs[k] for k in column_keys]].iterrows():
            logger(i)
            logger('\t'+v.to_string().replace('\n','\n\t'))

    def add_table_one(self,
        csv_file,
        column_labels,
        output_suffix = None,
        ):

        table_one = pandas.read_csv(csv_file, index_col=0, dtype=str).transpose()

        # Separate high- and low-resolution limits
        if 'Resolution range' in table_one.columns:
            table_one['Low Resolution Limit'], table_one['High Resolution Limit'] = list(zip(*table_one['Resolution range'].apply(lambda x: x.split('(')[0].split('-'))))

        # Select columns that exist
        table_one = table_one[table_one.columns.intersection(column_labels)]

        # Remove "high-resolution shell" statistics
        for col in table_one.columns:
            table_one.loc[:,col] = table_one.loc[:,col].apply(lambda x: x.split('(')[0])

        # Redetermine data types
        table_one = table_one.apply(lambda x: pandas.to_numeric(x, errors='coerce'))

        # Calculate additional columns -- only after coersion to numeric!
        if ('R-free' in table_one.columns) and ('R-work' in table_one.columns):
            table_one.loc[:,'R-gap'] = table_one['R-free'] - table_one['R-work']

        # Check all requested columns are present
        missing_cols = set(column_labels).difference(table_one.columns)
        if missing_cols:
            raise Sorry('Requested labels are not present in supplied csv file ({}): \n{}'.format(csv_file, ', '.join(missing_cols)))

        # Add suffix
        if output_suffix is not None:
            table_one.columns = table_one.columns + output_suffix

        # Create blank columns
        for col in table_one.columns.difference(self.table.columns):
            self.table.loc[:,col] = None
        # Transfer data to other
        self.table.update(table_one, overwrite=True)

        return self.table

    def calculate_column_differences(self,
        column_labels,
        input_suffixes,
        output_suffix,
        ):

        assert len(input_suffixes) == 2

        table = self.table

        suff1 = input_suffixes[0]
        suff2 = input_suffixes[1]

        # Create columns for the deltas between variables
        for col in column_labels:
            col1 = col + suff1
            col2 = col + suff2
            col3 = col + output_suffix
            if {col1, col2}.difference(table.columns):
                raise Failure('Columns {} and {} not present in table'.format(col1, col2))
            table.loc[:,col3] = table.loc[:,col1] - table.loc[:,col2]

    def show_average_statistics_summary(self,
        single_column_labels_dict,
        comparison_column_prefixes_dict,
        comparison_column_suffix_triplets,
        print_strip_chars = '',
        label_hash = {},
        ):
        """Display overall statistics of columns and compare sets of related columns"""

        table = self.table

        # Extract Column averages (means and medians)
        table_means = table.mean().round(3)
        shrt_str = '{:>35s} | {:>15} |'
        long_str = '{:>35s} | {:>15} | {:>15} | {:>15}'

        logger.subheading('Average statistics')

        # Columns without old/new prefix
        for col_key, col in single_column_labels_dict.items():
            logger(shrt_str.format(col, table_means[col]))

        logger.subheading('Average difference between different models')

        for triplet in comparison_column_suffix_triplets:

            logger(long_str.format(
                    '',
                    label_hash.get(triplet[0],triplet[0]).strip(print_strip_chars),
                    label_hash.get(triplet[1],triplet[1]).strip(print_strip_chars),
                    label_hash.get(triplet[2],triplet[2]).strip(print_strip_chars))
                )

            for col_prefix_key, col_prefix in comparison_column_prefixes_dict.items():

                if col_prefix_key in self.r_factor_keys:
                    multiplier = 100.0
                    col_print_suffix = ' (%)'
                    str_format = '{:.1f}'
                else:
                    multiplier = 1.0
                    col_print_suffix = ''
                    str_format = '{:.3f}'

                col1, col2, col3 = [col_prefix+suff for suff in triplet]

                logger(long_str.format(col_prefix+col_print_suffix,
                    str_format.format(multiplier*table_means[col1]),
                    str_format.format(multiplier*table_means[col2]),
                    str_format.format(multiplier*table_means[col3])))

            logger.bar()

        return

    def write(self):
        """Add data to CSV and write"""
        self.table.dropna(axis='columns', how='all').to_csv(self.filename)

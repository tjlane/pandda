import os
import math, collections
import numpy, pandas

from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure

from bamboo.maths.functions import rms
from giant.structure.uij import uij_to_b

EIGHT_PI_SQ = 8*math.pi*math.pi


class PandemicTrackingObject:


    def __init__(self,
            output_directory,
            plotting_object,
            verbose = False,
            log = None,
            ):
        if log is None: log = Log()
        # Create table for tracking progress over cycles
        table = pandas.DataFrame(columns=[
            'cycle', 'step', 'level','rmsd',
            'u_iso (level)', 'b_iso (level)',
            'b_min (level)', 'b_max (level)',
            'u_iso (total)', 'b_iso (total)',
            ])
        tracking_csv = os.path.join(output_directory, 'tracking_data.csv')
        tracking_png1 = os.path.join(output_directory, 'tracking_data_1.png')
        tracking_png2 = os.path.join(output_directory, 'tracking_data_2.png')
        i_cycle = 0
        adopt_init_args(self, locals())

    def update(self,
            uij_target,
            uij_lvl,
            step,
            i_level=None,
            write_graphs=False,
            ):
        """Update the tracking table"""

        log = self.log
        log.subheading('Updating tracking...')

        # Extract uijs for all of the levels for all datasets
        uij_tot = uij_lvl.sum(axis=0)
        # Calculate the rms between fitted and input
        rmsd = EIGHT_PI_SQ*rms(uij_target-uij_tot, axis=None)

        # Average over all datasets
        uij_tot = uij_tot.mean(axis=0)
        uij_lvl = uij_lvl.mean(axis=1)

        if not isinstance(i_level, list):
            i_level = [i_level]

        # Iterate through levels to be dumped into table
        for i_l in i_level:

            # Extract the Uij for the selected level(s)
            if isinstance(i_l, int):
                uij_sel = uij_lvl[i_l]
                level = i_l+1
            else:
                uij_sel = None
                level = None

            # Calculate U-iso & B-iso for selected level
            if uij_sel is not None:
                # average values
                b_iso_sel = numpy.mean(uij_to_b(uij_sel))
                u_iso_sel = b_iso_sel / EIGHT_PI_SQ
                # min/max values
                b_min_sel = numpy.min(uij_to_b(uij_sel))
                b_max_sel = numpy.max(uij_to_b(uij_sel))
            else:
                # average values
                b_iso_sel = 0.0
                u_iso_sel = 0.0
                # min/max values
                b_min_sel = numpy.nan
                b_max_sel = numpy.nan

            # Calculate U-iso & B-iso for complete model
            b_iso_tot = numpy.mean(uij_to_b(uij_tot))
            u_iso_tot = b_iso_tot / EIGHT_PI_SQ

            # Create human-readable cycle number
            cycle_lab = self.i_cycle

            # Add to tracking table
            self.table.loc[len(self.table.index)] = {
                'cycle' : cycle_lab,
                'step'  : step,
                'level' : level,
                'rmsd'  : round(rmsd,3),
                'u_iso (level)' : round(u_iso_sel,3),
                'b_iso (level)' : round(b_iso_sel,3),
                'b_min (level)' : round(b_min_sel,3),
                'b_max (level)' : round(b_max_sel,3),
                'u_iso (total)' : round(u_iso_tot,3),
                'b_iso (total)' : round(b_iso_tot,3),
                }

        log(self.table.loc[len(self.table)-len(i_level):].to_string())

        if write_graphs:
            # Dump to csv
            self.table.to_csv(self.tracking_csv)
            # Make plots
            self.plotting_object.tracking_plots(
                    table = self.table,
                    filename = self.tracking_png1,
                    )
            self.plotting_object.convergence_plots(
                    table = self.table,
                    filename = self.tracking_png2,
                    )


class PandemicResultsObject:


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
            verbose = False,
            log = None,
            ):
        if log is None: log = Log()

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
        log = self.log

        # Internal keys for dicts
        column_keys = self.r_factor_keys

        # Store suffix for reference
        self.reference_suffix = output_suffix

        # Form input and output label dictionaries 
        i_labs = input_column_labels
        o_labs = {k:v+output_suffix for k,v in self.output_column_labels_base.iteritems()}
        
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

        log.subheading('Adding reference values to results table')
        
        log('Reading Reference R-values from {}'.format(filename))

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
            log('\nCalculating R-gap column for reference R-free & R-work values')
            self.reference_r_values.loc[:,i_rg_s] = self.reference_r_values.loc[:,i_rf_s] - self.reference_r_values.loc[:,i_rw_s]
        else:
            self.log('\nR-gap column already present in reference table -- using this column.')

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

        self.log('\n> Processed Reference Values\n')
        for i, v in table[[o_labs[k] for k in column_keys]].iterrows():
            log(i)
            log('\t'+v.to_string().replace('\n','\n\t'))

    def add_table_one(self,
        csv_file,
        column_labels,
        output_suffix = None,
        ):

        table_one = pandas.read_csv(csv_file, index_col=0, dtype=str).transpose()

        # Separate high- and low-resolution limits
        if 'Resolution range' in table_one.columns:
            table_one['Low Resolution Limit'], table_one['High Resolution Limit'] = zip(*table_one['Resolution range'].apply(lambda x: x.split('(')[0].split('-')))

        # Select columns that exist
        #table_one = table_one[table_one.intersection(column_labels)]
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
        log = self.log

        # Extract Column averages (means and medians)
        table_means = table.mean().round(3)
        shrt_str = '{:>35s} | {:>15} |'
        long_str = '{:>35s} | {:>15} | {:>15} | {:>15}'

        log.subheading('Average statistics')

        # Columns without old/new prefix
        for col_key, col in single_column_labels_dict.items():
            log(shrt_str.format(col, table_means[col]))

        log.subheading('Average difference between different models')

        for triplet in comparison_column_suffix_triplets:

            log(long_str.format(
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

                log(long_str.format(col_prefix+col_print_suffix,
                    str_format.format(multiplier*table_means[col1]),
                    str_format.format(multiplier*table_means[col2]),
                    str_format.format(multiplier*table_means[col3])))
            
            log.bar()

        return

    # TODO
    def calculate_uij_statistics(self, 
        uij_target,
        uij_levels,
        uij_residual,
        ):

        self.log.bar(True, False)
        self.log('Calculating mean and median fitting rmsds by dataset')
        self.log.bar()
        # Calculate rmsd between input and fitted uijs
        uij_rmsd = rms(uij_inp-uij_fit, axis=2)
        # Extract mean/median dataset-by-dataset RMSDs
        dset_medn_rmsds = numpy.median(uij_rmsd, axis=1)
        dset_mean_rmsds = numpy.mean(uij_rmsd, axis=1)
        main_table['Mean Fitting RMSD']   = dset_mean_rmsds
        main_table['Median Fitting RMSD'] = dset_medn_rmsds
        for i in xrange(0, min(10, len(self.models))):
            self.log('Model {:10}: {:6.3f} (mean), {:6.3f} (median)'.format(self.models[i].tag, dset_mean_rmsds[i], dset_medn_rmsds[i]))

        self.log.bar(True, False)
        self.log('Calculating isotropic ADPs for input and fitted ADPs by dataset')
        self.log.bar()
        # Calculate isotropic ADPs for input and fitted uijs
        uij_inp_iso = numpy.array(map(uij_to_b, uij_inp))
        uij_fit_iso = numpy.array(map(uij_to_b, uij_fit))
        # Calculate mean/median ADPs for each atom
        dset_mean_inp_iso = numpy.mean(uij_inp_iso, axis=1)
        dset_mean_fit_iso = numpy.mean(uij_fit_iso, axis=1)
        main_table['Average B-factor (fitted atoms) (Input)']               = dset_mean_inp_iso
        main_table['Average B-factor (fitted atoms) (Fitted)']              = dset_mean_fit_iso
        main_table['Average B-factor (fitted atoms) Change (Fitted-Input)'] = dset_mean_fit_iso - dset_mean_inp_iso
        for i in xrange(0, min(10, len(self.models))):
            self.log('Model {:10}: {:6.3f} (input) -> {:6.3f} (fitted)'.format(self.models[i].tag, dset_mean_inp_iso[i], dset_mean_fit_iso[i]))

        # Store output table (object is changed in the join steps)
        self.tables.statistics = main_table

    def write(self):
        """Add data to CSV and write"""
        self.table.dropna(axis='columns', how='all').to_csv(self.filename)

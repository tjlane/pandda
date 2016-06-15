from __future__ import print_function

import os, sys, glob, time, re
import copy, warnings, traceback

import numpy, pandas

from bamboo.common import Meta, Info
from bamboo.common.logs import Log
from bamboo.common.file import FileManager
from bamboo.common.path import easy_directory, rel_symlink

from giant.manager import Program
from giant.structure import make_label

from pandda.constants import PanddaMaskNames

from pandemic import PANDEMIC_TOP, PANDEMIC_TEXT, PANDEMIC_VERSION
from pandemic.phil import pandemic_phil
from pandemic.constants import *

class PandemicMultiDatasetAnalyser(Program):
    """Class for the identification and characterisation of heterogeneity and flexibility of a series of isomorphous-ish crystals"""

    _NAME    = 'pandemic.analyse'
    _TEXT    = PANDEMIC_TEXT
    _VERSION = PANDEMIC_VERSION

    def __init__(self, params):
        """Class for the identification and characterisation of heterogeneity and flexibility of a series of isomorphous-ish crystals"""

        # Log init time
        self._init_time = time.time()

        # ===============================================================================>
        # PROCESS INPUT ARGUMENTS
        # ===============================================================================>

        # Read in the master phil
        self.master_phil = pandemic_phil
        # Store the whole params object
        self._input_params = params
        # Pull out the python object of the arguments (at the `pandemic` level)
        self.args = self._input_params.pandemic
        # Most of the variables are contained within the params object - create a shortcut to save typing
        self.params = self._input_params.pandemic.params
        # Program settings (num cpus, etc...)
        self.settings = self._input_params.settings

        # ===============================================================================>
        # OUTPUT FILES STUFF
        # ===============================================================================>

        assert self.args.output.out_dir, 'pandemic.output.out_dir IS NOT DEFINED'
        self.out_dir = easy_directory(os.path.abspath(self.args.output.out_dir))

        # Create a log for the object
        self.log = Log(log_file=os.path.join(self.out_dir, 'pandemic.log'), verbose=self.settings.verbose)

        # ===============================================================================>
        # SETTINGS STUFF
        # ===============================================================================>
        self._high_resolution = None
        self._low_resolution = None

        # ===============================================================================>
        # ANALYSIS OBJECTS
        # ===============================================================================>

        # Create tables object (to hold pandas dataframe objects)
        self.tables = Meta()
        # Record global information about the datasets
        self.tables.dataset_info        = pandas.DataFrame( data    = None,
                                                            index   = pandas.Index(data=[], name='dtag'),
                                                            columns = PandemicTableFields.all_dataset_fields      )
        # Record information about each residue/conformer
        self.tables.residue_info        = pandas.DataFrame( data    = None,
                                                            index   = pandas.MultiIndex(levels=[[],[]], labels=[[],[]], names=['chain','residue']),
                                                            columns = PandemicTableFields.all_residue_fields      )

        # ===============================================================================>
        # DIRECTORY STUFF
        # ===============================================================================>

        # Set up the output folders and files
        self._directory_setup()
        # Set up the pickled object filenames
        self._pickle_setup()

    def _directory_setup(self):
        """Initialise the pandemic directory system"""

        # Create a file and directory organiser
        self.output_handler = FileManager(rootdir=self.out_dir)

        # Filename templates
        f = PandemicAnalyserFilenames

        # ===============================================================================>
        # Global Directories that do not change from run to run
        # ===============================================================================>

        self.output_handler.add_dir(dir_name='analyses',              dir_tag='analyses',             top_dir_tag='root', create=False, exists=False)
        self.output_handler.add_dir(dir_name='pickles',               dir_tag='pickles',              top_dir_tag='root', create=False, exists=False)

        # ================================================>
        # Input + Status parameters
        # ================================================>
        self.output_handler.add_file(file_name='pandemic.eff',        file_tag='settings',         dir_tag='root')
        self.output_handler.add_file(file_name='pandemic.{}',         file_tag='status',           dir_tag='root')

        self.output_handler.add_dir(dir_name='residue_maps', dir_tag='residue_maps', top_dir_tag='root', create=False, exists=False)
        self.output_handler.add_file(file_name='{}_residue_{}.ccp4',  file_tag='residue_map',      dir_tag='residue_maps')
        self.output_handler.add_dir(dir_name='dataset_maps', dir_tag='dataset_maps', top_dir_tag='root', create=False, exists=False)
        self.output_handler.add_file(file_name='dataset_{}.ccp4',     file_tag='dataset_map',      dir_tag='dataset_maps')


    def _pickle_setup(self):
        """Initialise all of the pickle filenames"""

        # Pickle Handler
        self.pickle_handler = FileManager(rootdir=self.output_handler.get_dir('pickles'))

    def run_analysis_init(self):
        """Set up the pandemic for a new analysis (doing this will override links to analyses)"""

        # ===============================================================================>
        # Validate the input parameters
        # ===============================================================================>
        self._validate_parameters()

        # ===============================================================================>
        # Update the FileManager to make sure all directories are now created
        # ===============================================================================>
        self.output_handler.check_and_create_directories()

        # ===============================================================================>
        # Write the header to the log file
        # ===============================================================================>
        # Print logo and parameters to log
        self.log(PANDEMIC_TEXT, True)
        self.write_running_parameters_to_log()
        # Write the used parameters to file
        with open(self.output_handler.get_file('settings'), 'w') as out_file:
            out_file.write( '\n'.join([ '# Command Line Args',
                                        '# ',
                                        # This line assures that strings are quoted for easy copy-pasting
                                        '# pandemic.analyse '+' '.join(sys.argv[1:]).replace(' ','" ').replace('=','="')+'"',
                                        '',
                                        '# Used Settings:',
                                        '',
                                        self.master_phil.format(python_object=self._input_params).as_str() ]))

        # ===============================================================================>
        # PRINT SOME HELPFUL INFORMATION
        # ===============================================================================>
        self.log('===================================>>>', True)
        self.log('RUNNING FROM: {!s}'.format(sys.argv[0]), True)
        self.log('===================================>>>', True)
        self.log('READING INPUT PANDDA FROM: {!s}'.format(self.args.input.pandda_dir))
        self.log('===================================>>>', True)
        self.log('WRITING OUTPUT TO: {!s}'.format(self.out_dir), True)

        # ===============================================================================>
        # LOOK FOR MATPLOTLIB TO SEE IF WE CAN GENERATE GRAPHS
        # ===============================================================================>
        if self.settings.plot_graphs:
            if self.check_for_matplotlib(): pass
            else: self.settings.plot_graphs = False

        # ===============================================================================>
        # LOG THE START TIME
        # ===============================================================================>
        # Store start time and print to log
        self.log('===================================>>>', True)
        self.log('Analysis Started: {!s}'.format(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime(self._init_time))), True)
        # ===============================================================================>
        # Update Status
        # ===============================================================================>
        self.update_status('running')

    def _validate_parameters(self):
        """Validate and preprocess the loaded parameters"""

        p = self.args

    def initialise_dataset_masks_and_tables(self, pandda):
        """Use the datasets in the pandda to populate the fields in the pandemic dataset tables"""

        self.log('===================================>>>', True)
        self.log('Initialising PanDEMIC Dataset Information Tables.', True)

        # Add dataset tags as rows in the tables
        self.tables.dataset_info = self.tables.dataset_info.append(pandas.DataFrame(index=[d.tag for d in pandda.datasets.all()]), verify_integrity=True)

    def initialise_residue_masks_and_tables(self, pandda):
        """Use the residues of the structures to populate the fields in the pandemic residue tables"""

        self.log('===================================>>>', True)
        self.log('Initialising PanDEMIC Residue Information Tables.', True)

        # Extract the residue labels for the reference structure
        residue_labels = [make_label(rg) for rg in pandda.reference_dataset().hierarchy().residue_groups()]
        for res_lab in residue_labels: self.tables.residue_info.loc[res_lab,:]=None

    def select_datasets_and_load_maps(self, pandda, high_res_large_cutoff, high_res_small_cutoff=0):
        """Select datasets for the multi-conformer analysis and return a map analyser object for them"""

        # Mask name - add pandemic so that does not conflict with pandda names
        analysis_mask_name = 'pandemic analysis @ {}A'.format(high_res_large_cutoff)
        # Select datasets for analysis
        analysis_mask = pandda.select_for_analysis(high_res_large_cutoff = high_res_large_cutoff,
                                                   high_res_small_cutoff = high_res_small_cutoff,
                                                   analysis_mask_name = analysis_mask_name)
        if sum(analysis_mask) == 0:
            raise Exception('NO DATASETS TO ANALYSE for {!s}'.format(analysis_mask_name))
        else:
            self.log('{} Datasets selected for analysis at {}A'.format(sum(analysis_mask), high_res_large_cutoff))

        # Truncate data and load maps
        pandda.truncate_scaled_data(dataset_handlers = pandda.datasets.mask(mask_name=analysis_mask_name))
        # Reference map
        ref_map_holder = pandda.load_reference_map( map_resolution = high_res_large_cutoff )
        # Dataset maps
        map_holder_list = pandda.load_and_morph_maps( dataset_handlers = pandda.datasets.mask(mask_name=analysis_mask_name),
                                                      ref_map_holder   = ref_map_holder,
                                                      map_resolution   = high_res_large_cutoff )

        return map_holder_list

    def find_residue_locales(self, pandda):
        """Mask each residue on the grid, and convolve with the grid partition for each residue"""

        # Extract grid objects
        ref_grid = pandda.reference_grid()
        ref_part = ref_grid.partition()
        out_mask = ref_grid.global_mask().outer_mask_binary()

        # Hash for assigning grid points to each residue
        residue_grid_hash = {}

        # Iterate through each residue and find the grid points assigned to each residue
        for i_res, res in enumerate(pandda.reference_dataset().calpha_labels()):
            # Get the indices for this residue from the partition
            res_idxs = numpy.multiply(out_mask, ref_part.nn_groups==i_res).nonzero()
            # Record the indices for this residue (only need the first element since 1d array)
            residue_grid_hash[res] = res_idxs[0]

        return residue_grid_hash

    def pickle_the_pandemic(self, components=None, all=False, datasets=None):
        pass

    def exit(self, error=False):
        """Exit the PANDEMIC, record runtime etc..."""

        self._finish_time = time.time()
        self.log('Runtime: {!s}'.format(time.strftime("%H hours:%M minutes:%S seconds", time.gmtime(self._finish_time - self._init_time))))

        # If error, don't make meta or self pickle
        if error:
            self.update_status('errored')
            self.log('===================================>>>', True)
            self.log('PANDEMIC exited with an error')
            self.log('===================================>>>', True)
            self.log('Error Traceback: ')
            self.log('===================================>>>', True)
            self.log(traceback.format_exc())
            self.log('===================================>>>', True)
            return
        else:
            self.update_status('done')
            self.log('===================================>>>', True)
            self.log('. FINISHED! PANDEMIC EXITED NORMALLY .', True)
            self.log('===================================>>>', True)





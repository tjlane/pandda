from __future__ import print_function

import os, sys, glob, time, re
import copy, warnings, traceback

import numpy, pandas

from libtbx import easy_mp
from libtbx.utils import Sorry, Failure
from libtbx.math_utils import ifloor, iceil

import scitbx.matrix

from scitbx.array_family import flex

from bamboo.common import Meta, Info
from bamboo.common.logs import Log
from bamboo.common.file import FileManager
from bamboo.common.path import easy_directory, rel_symlink, delete_with_glob
from bamboo.common.math import round_no_fail

from bamboo.common.status import status_bar, status_bar_2
from bamboo.common.command import CommandManager

from bamboo.plot import bar, simple_histogram

from giant.common import ModelAndData, ElectronDensityMap
from giant.holders import HolderList
from giant.manager import Program
from giant.grid import Grid
from giant.grid.masks import AtomicMask, GridMask
from giant.structure.align import GlobalAlignment
from giant.structure.select import calphas, protein, sel_altloc
from giant.structure.formatting import Labeller, ShortLabeller

from pandda.phil import pandda_phil
from pandda.misc import write_array_as_map, write_indices_as_map, write_bool_as_map
from pandda.analyse.events import cluster_events
from pandda.analyse.functions import DatasetAligner, MapLoader, DensityStatistics, UncertaintyCalculator, wrapper_run
from pandda.constants import *
from pandda import HEADER_TEXT, VERSION


class PanddaDataset(ModelAndData):


    def __init__(self, model, data):
        """Subclass of ModelAndData used for PanDDA Analysis"""

        super(PanddaDataset, self).__init__(model=model, data=data)

        self.child = None
        self.events = []


class PanddaReferenceDataset(ModelAndData):


    _origin_shift=None

    def set_origin_shift(self, shift):
        """Creates an alignment corresponding to an origin shift"""

        self._origin_shift = shift
        r = scitbx.matrix.rec([1,0,0,0,1,0,0,0,1], (3,3))
        t = scitbx.matrix.rec(shift, (3,1))
        rt = scitbx.matrix.rt((r,t))
        self.model.alignment = GlobalAlignment(alignment_mx=rt, alignment_sites=calphas(self.model.hierarchy).atoms().extract_xyz(),  id='ref')
        return self.model.alignment

    def origin_shift(self):
        return self._origin_shift

    def nat2grid(self, *args, **kwargs):
        return self.model.alignment.nat2ref(*args, **kwargs)
    def grid2nat(self, *args, **kwargs):
        return self.model.alignment.ref2nat(*args, **kwargs)


class PanddaDatasetList(HolderList):


    _holder_class = PanddaDataset
    _reference_class = PanddaReferenceDataset
    _reference = None

    def _get_num(self, item):
        return item.num
    def _get_tag(self, item):
        return item.tag

    def set_reference(self, dataset):
        assert isinstance(dataset, self._reference_class), 'reference not of the right class: {}'.format(self._reference_class)
        self._reference = dataset

    def reference(self):
        return self._reference


class MapHolderList(HolderList):


    _holder_class = ElectronDensityMap

    def _get_num(self, item):
        return item.meta.num
    def _get_tag(self, item):
        return item.meta.tag


class MapList(Info):


    _map_names = []

    def __init__(self, map_names=None):
        if map_names is None: map_names=[]
        assert self._map_names+map_names, 'No Maps defined'
        for m in self._map_names+map_names:
            self.__dict__[m] = None
        self.meta = Meta()
        self._initialized = True

    def __get_item__(self, item):
        return self.__dict__[item]


class PanddaStatMapList(MapList):


    _map_names = ['mean_map','medn_map','stds_map','sadj_map','skew_map','kurt_map','bimo_map']


class PanddaMultipleStatMapList(object):


    def __init__(self):
        """Store lots of statistical maps"""
        self.map_lists = {}

    def get_resolutions(self):
        return sorted(self.map_lists.keys())

    def add(self, stat_map_list, resolution, overwrite=False):
        assert isinstance(resolution, float), 'Resolution of map must be of type float. Type given: {!s}'.format(type(resolution))
        if overwrite is not False: assert resolution not in self.map_lists.keys(), 'MAPS OF THIS RESOLUTION ALREADY ADDED'
        assert isinstance(stat_map_list, PanddaStatMapList), 'stat_map_list must be of type PanddaMultipleStatMapList. Type given: {!s}'.format(type(stat_map_list))
        self.map_lists[resolution] = stat_map_list

    def get(self, resolution):
        return self.map_lists[resolution]


class PanddaMultiDatasetAnalyser(Program):


    _NAME    = 'pandda.analyse'
    _VERSION = VERSION
    _TEXT    = HEADER_TEXT

    def __init__(self, params):
        """Main PanDDA Class"""

        # Log init time
        self._init_time = time.time()

        # ===============================================================================>
        # PROCESS INPUT ARGUMENTS
        # ===============================================================================>

        self.master_phil = pandda_phil

        self._input_params = params
        self.args = self._input_params.pandda
        self.params = self._input_params.pandda.params
        self.settings = self._input_params.settings

        # ===============================================================================>
        # OUTPUT FILES STUFF
        # ===============================================================================>

        assert self.args.output.out_dir, 'pandda.output.out_dir IS NOT DEFINED'
        self.out_dir = easy_directory(os.path.abspath(self.args.output.out_dir))

        # Create a log for the object
        self.log = Log(log_file=os.path.join(self.out_dir, 'pandda-{}.log'.format(time.strftime("%Y-%m-%d-%H%M", time.gmtime()))), verbose=self.settings.verbose)

        # ===============================================================================>
        # SETTINGS STUFF
        # ===============================================================================>
        self._high_resolution = None
        self._low_resolution = None

        # ===============================================================================>
        # DATA AND MAPS STUFF
        # ===============================================================================>

        self._input_files = []
        self.grid = None
        self.datasets = PanddaDatasetList()
        self.pickled_dataset_meta = None
        self.stat_maps = PanddaMultipleStatMapList()

        # ===============================================================================>
        # ANALYSIS OBJECTS
        # ===============================================================================>

        # Create tables object (to hold pandas dataframe objects)
        self.tables = Meta()
        # Record global information about the datasets
        self.tables.dataset_info        = pandas.DataFrame( data    = None,
                                                            index   = pandas.Index(data=[], name='dtag'),
                                                            columns = PanddaTableFields.all_dataset_fields      )
        # Record information about the created maps for each dataset
        self.tables.dataset_map_info    = pandas.DataFrame( data    = None,
                                                            index   = pandas.Index(data=[], name='dtag'),
                                                            columns = PanddaTableFields.all_dataset_map_fields  )
        # Record the events detected in each dataset
        self.tables.event_info          = pandas.DataFrame( data    = None,
                                                            index   = pandas.MultiIndex(levels=[[],[]], labels=[[],[]], names=['dtag','event_idx']),
                                                            columns = PanddaTableFields.all_event_fields        )
        # Record information about the clustered events (cluster of events = site)
        self.tables.site_info           = pandas.DataFrame( data    = None,
                                                            index   = pandas.Index(data=[], name='site_idx'),
                                                            columns = PanddaTableFields.all_site_fields         )

        # ===============================================================================>
        # DIRECTORY STUFF
        # ===============================================================================>

        self._directory_setup()
        self._pickle_setup()

        if os.path.exists(self.pickle_handler.get_file('dataset_meta')):
            self._new_pandda = False
        else:
            self._new_pandda = True

    def _directory_setup(self):
        """Initialise the pandda directory system"""

        # Create a file and directory organiser
        fm = self.initialise_file_manager(rootdir=self.out_dir)

        # Filename templates
        f = PanddaAnalyserFilenames

        # ===============================================================================>
        # Global Directories that do not change from run to run
        # ===============================================================================>
#        fm.add_dir(dir_name='empty_directories',    dir_tag='empty_directories',    top_dir_tag='root', create=False, exists=False)
        fm.add_dir(dir_name='processed_datasets',   dir_tag='processed_datasets',   top_dir_tag='root', create=False, exists=False)
#        fm.add_dir(dir_name='rejected_datasets',    dir_tag='rejected_datasets',    top_dir_tag='root', create=False, exists=False)
#        fm.add_dir(dir_name='interesting_datasets', dir_tag='interesting_datasets', top_dir_tag='root', create=False, exists=False)
        fm.add_dir(dir_name='aligned_structures',   dir_tag='aligned_structures',   top_dir_tag='root', create=False, exists=False)
        fm.add_dir(dir_name='analyses',             dir_tag='analyses',             top_dir_tag='root', create=False, exists=False)

        # ================================================>
        # Input + Status parameters
        # ================================================>
        fm.add_file(file_name='pandda.eff',        file_tag='settings',         dir_tag='root')
        fm.add_file(file_name='pandda.{}',         file_tag='status',           dir_tag='root')

        # Somewhere to store the analysis summaries - for the user
        fm.add_dir(dir_name='html_summaries', dir_tag='output_summaries', top_dir_tag='analyses', create=False, exists=False)
        fm.add_file(file_name=f.initial_html,                      file_tag='initial_html',            dir_tag='output_summaries')
        fm.add_file(file_name=f.analyse_html,                      file_tag='analyse_html',            dir_tag='output_summaries')
        fm.add_file(file_name=f.analyse_site_graph,                file_tag='analyse_site_graph',      dir_tag='output_summaries')
        fm.add_file(file_name=f.analyse_site_graph_mult,           file_tag='analyse_site_graph_mult', dir_tag='output_summaries')
        fm.add_file(file_name=f.pymol_sites_py,                    file_tag='pymol_sites_py',          dir_tag='output_summaries')
        fm.add_file(file_name=f.pymol_sites_pml,                   file_tag='pymol_sites_pml',         dir_tag='output_summaries')
        fm.add_file(file_name=f.pymol_sites_png_1,                 file_tag='pymol_sites_png_1',       dir_tag='output_summaries')
        fm.add_file(file_name=f.pymol_sites_png_2,                 file_tag='pymol_sites_png_2',       dir_tag='output_summaries')
        # Store dataset summary graphs
        fm.add_dir(dir_name='dataset_summary_graphs', dir_tag='d_graphs', top_dir_tag='analyses', create=False, exists=False)
        fm.add_file(file_name='dataset_resolutions.png',           file_tag='d_resolutions',           dir_tag='d_graphs')
        fm.add_file(file_name='dataset_rfactors.png',              file_tag='d_rfactors',              dir_tag='d_graphs')
        fm.add_file(file_name='dataset_global_rmsd_to_ref.png',    file_tag='d_global_rmsd_to_ref',    dir_tag='d_graphs')
        fm.add_file(file_name='dataset_cell_axes.png',             file_tag='d_cell_axes',             dir_tag='d_graphs')
        fm.add_file(file_name='dataset_cell_angles.png',           file_tag='d_cell_angles',           dir_tag='d_graphs')
        fm.add_file(file_name='dataset_cell_volumes.png',          file_tag='d_cell_volumes',          dir_tag='d_graphs')
        # Store map analysis graphs
        fm.add_dir(dir_name='map_analysis_summary_graphs', dir_tag='m_graphs', top_dir_tag='analyses', create=False, exists=False)
        fm.add_file(file_name='{}A-reference_v_mean_unsorted.png',  file_tag='ref_v_mean_map_unsort',  dir_tag='m_graphs')
        fm.add_file(file_name='{}A-reference_v_mean_sorted.png',    file_tag='ref_v_mean_map_sort',    dir_tag='m_graphs')
        fm.add_file(file_name='{}A-reference_map_distribution.png', file_tag='ref_map_dist',           dir_tag='m_graphs')
        # Somewhere to store the dataset information (general values)
        fm.add_file(file_name=f.dataset_info,                      file_tag='dataset_info',            dir_tag='analyses')
        fm.add_file(file_name=f.dataset_map_info,                  file_tag='dataset_map_info',        dir_tag='analyses')
        fm.add_file(file_name=f.dataset_combined_info,             file_tag='dataset_combined_info',   dir_tag='analyses')
        fm.add_file(file_name=f.dataset_masks,                     file_tag='dataset_masks',           dir_tag='analyses')
        # Somewhere to store the dataset information (identified events + sites)
        fm.add_file(file_name=f.event_info,                        file_tag='event_info',              dir_tag='analyses')
        fm.add_file(file_name=f.site_info,                         file_tag='site_info',               dir_tag='analyses')
        fm.add_file(file_name='_point_distributions.csv',          file_tag='point_distributions',     dir_tag='analyses')

        # Somewhere to store the pickled objects
        fm.add_dir(dir_name='pickled_data', dir_tag='pickles', top_dir_tag='root', create=False, exists=False)

        # ===============================================================================>
        # Reference Structure Files (should only be needed once for writing and then only for reloading)
        # ===============================================================================>
        fm.add_dir(dir_name='reference', dir_tag='reference', top_dir_tag='root', create=False, exists=False)
        fm.add_file(file_name=f.reference_structure,               file_tag='reference_structure', dir_tag='reference')
        fm.add_file(file_name=f.reference_dataset,                 file_tag='reference_dataset',   dir_tag='reference')
        fm.add_file(file_name=f.reference_on_origin,               file_tag='reference_on_origin', dir_tag='reference')
        fm.add_file(file_name=f.reference_symmetry,                file_tag='reference_symmetry',  dir_tag='reference')
        fm.add_file(file_name='grid-voronoi-{}.ccp4',              file_tag='grid_voronoi',        dir_tag='reference')

        # ===============================================================================>
        # Standard template files that will be populated when needed
        # ===============================================================================>
        fm.add_dir(dir_name='statistical_maps', dir_tag='statistical_maps', top_dir_tag='reference', create=False, exists=False)
        fm.add_file(file_name=f.mean_map,                          file_tag='mean_map',            dir_tag='statistical_maps')
        fm.add_file(file_name=f.medn_map,                          file_tag='medn_map',            dir_tag='statistical_maps')
        fm.add_file(file_name=f.stds_map,                          file_tag='stds_map',            dir_tag='statistical_maps')
        fm.add_file(file_name=f.sadj_map,                          file_tag='sadj_map',            dir_tag='statistical_maps')
        fm.add_file(file_name=f.skew_map,                          file_tag='skew_map',            dir_tag='statistical_maps')
        fm.add_file(file_name=f.kurt_map,                          file_tag='kurt_map',            dir_tag='statistical_maps')
        fm.add_file(file_name=f.bimo_map,                          file_tag='bimo_map',            dir_tag='statistical_maps')

    def _pickle_setup(self):
        """Initialise all of the pickle filenames"""

        # Pickle Handler
        self.pickle_handler = FileManager(rootdir=self.file_manager.get_dir('pickles'))
        # Pickled Reference Objects
        self.pickle_handler.add_file(file_name='grid.pickle',    file_tag='grid')
        self.pickle_handler.add_file(file_name='reference_dataset.pickle', file_tag='reference_dataset')
        # Pickled Datasets
        self.pickle_handler.add_file(file_name='dataset_masks.pickle',     file_tag='dataset_masks')
        self.pickle_handler.add_file(file_name='dataset_meta.pickle',      file_tag='dataset_meta')
        # Pickled Statistical Maps
        self.pickle_handler.add_file(file_name='statistical_maps.pickle',  file_tag='stat_maps')
        # Map Analysers
        self.pickle_handler.add_file(file_name='map_analyser_{}A.pickle',  file_tag='map_analyser')
        # Pickled Pandda -- (main object)
        self.pickle_handler.add_file(file_name='my_pandda.pickle',         file_tag='my_pandda')

    def run_analysis_init(self):
        """Set up the pandda for a new analysis (doing this will override links to analyses)"""

        # ================================================>
        # Validate the input parameters
        # ================================================>
        self._validate_parameters()

        # ================================================>
        # Create a new analysis directory for analyses/summaries
        # ================================================>
        # New directories will be created for each run (so that data is not overwritten) by the time of the run
        if self.args.output.new_analysis_dir or (not os.path.exists(self.file_manager.get_dir('analyses'))):
            analysis_time_name = 'analyses-{}'.format(time.strftime("%Y-%m-%d-%H%M", time.gmtime(self._init_time)))
            analysis_time_path = easy_directory(os.path.join(self.file_manager.get_dir('root'), analysis_time_name))
            analysis_link_path = self.file_manager.get_dir('analyses')
            # Remove old analysis link if it exists and link in the new analysis directory
            if os.path.exists(analysis_link_path) and os.path.islink(analysis_link_path): os.unlink(analysis_link_path)
            rel_symlink(orig=analysis_time_path, link=analysis_link_path)

        assert os.path.exists(self.file_manager.get_dir('analyses')), 'Output analysis directory does not exist'

        # ===============================================================================>
        # Update the FileManager to make sure all directories are now created
        # ===============================================================================>
        self.file_manager.check_and_create_directories()

        # ===============================================================================>
        # Report
        # ===============================================================================>
        self.log('', True)
        self.log('################################### <~> ###################################', True)
        self.log('              Pre-analysis checks and directory setup complete', True)
        self.log('################################### <~> ###################################', True)
        self.log('', True)

        # ===============================================================================>
        # Write the header to the log file
        # ===============================================================================>
        self.log(HEADER_TEXT, True)
        self.write_running_parameters_to_log()
        # Write the used parameters to file
        with open(self.file_manager.get_file('settings'), 'w') as out_file:
            out_file.write( '\n'.join([ '# Command Line Args',
                                        '# ', # The line below assures that strings are quoted for easy copy-pasting
                                        '# '+self._NAME+' '+' '.join(sys.argv[1:]).replace(' ','" ').replace('=','="')+'"',
                                        '',
                                        '# Used Settings:',
                                        '',
                                        self.master_phil.format(python_object=self._input_params).as_str() ]))

        # ===============================================================================>
        # PRINT SOME HELPFUL INFORMATION
        # ===============================================================================>

        self.log('', True)
        self.log('################################### <~> ###################################', True)
        self.log('                         FILES/MODULE INFORMATION', True)
        self.log('################################### <~> ###################################', True)
        self.log('', True)
        self.log('Running from: {!s}'.format(sys.argv[0]), True)
        self.log('----------------------------------->>>', True)
        self.log('Reading input from : {!s}'.format(self.args.input.data_dirs), True)
        self.log('----------------------------------->>>', True)
        self.log('Writing output to: {!s}'.format(self.out_dir), True)

        # ===============================================================================>
        # LOOK FOR MATPLOTLIB TO SEE IF WE CAN GENERATE GRAPHS
        # ===============================================================================>

        if self.settings.plot_graphs:
            self.log('----------------------------------->>>', True)
            if self.check_for_matplotlib(backend=self.settings.plotting.backend, interactive=False): pass
            else: self.settings.plot_graphs = False

        # ===============================================================================>
        # CHANGE INTO OUTPUT DIRECTORY
        # ===============================================================================>
        self.log('----------------------------------->>>', True)
        self.log('Changing into output directory: {}'.format(self.out_dir), True)
        os.chdir(self.out_dir)

        # ===============================================================================>
        # REPOPULATE PANDDA FROM PREVIOUS RUNS
        # ===============================================================================>
        # Load any objects from previous runs
        self.log('----------------------------------->>>', True)
        self.log('Checking for existing analyses', True)
        self.load_pickled_objects()

        # Reload reference dataset
        if (not self.datasets.reference()) and os.path.exists(self.file_manager.get_file('reference_structure')) and os.path.exists(self.file_manager.get_file('reference_dataset')):
            self.log('----------------------------------->>>', True)
            self.log('Loading Reference Dataset', True)
            self.load_reference_dataset(ref_pdb=self.file_manager.get_file('reference_structure'), ref_mtz=self.file_manager.get_file('reference_dataset'))

        # ===============================================================================>
        # LOG THE START TIME AND UPDATE STATUS
        # ===============================================================================>
        self.log('', True)
        self.log('################################### <~> ###################################', True)
        self.log('                PROGRAM SETUP COMPLETE - starting analysis.')
        self.log('################################### <~> ###################################', True)
        self.log('', True)
        self.log('----------------------------------->>>', True)
        self.log('Analysis Started: {!s}'.format(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime(self._init_time))), True)
        self.log('----------------------------------->>>', True)
        self.update_status('running')

    def _validate_parameters(self):
        """Validate and preprocess the loaded parameters"""

        self.log('')
        self.log('################################### <~> ###################################', True)
        self.log('                 validating and updating input parameters', True)
        self.log('################################### <~> ###################################', True)
        self.log('')

        p = self.args

        # Input
        assert p.input.data_dirs is not None, 'pandda.input.data_dirs IS NOT DEFINED'
        assert p.input.pdb_style, 'pandda.input.pdb_style IS NOT DEFINED'
        if not p.input.mtz_style: self.log('mtz_style not provided - using pdb style: {} -> {}'.format(p.input.pdb_style, p.input.pdb_style.replace('.pdb','.mtz')))
        assert p.input.regex.pdb_regex or p.input.regex.mtz_regex or p.input.regex.dir_regex or \
               (p.input.pdb_style and ('*' in p.input.pdb_style)) or \
               (p.input.mtz_style and ('*' in p.input.mtz_style)) or \
               (p.input.data_dirs and ('*' in p.input.data_dirs)), 'No method has been provided for labelling the datasets'
        if not p.input.lig_style: self.log('pandda.input.lig_style has not been provided - ligand files will not be detected')

        assert p.params.filtering.flags.same_space_group_only, 'PanDDA does not currently support the analysis of different spacegroups'

        # Make fullpath so we can run on the eff file from anywhere and change directories without worrying about relative paths
        p.input.data_dirs = os.path.abspath(p.input.data_dirs)
        p.output.out_dir  = os.path.abspath(p.output.out_dir)
        if p.input.filter.pdb:
            p.input.filter.pdb = os.path.abspath(p.input.filter.pdb)

        # If any datasets are set to be reprocessed, reload all datasets (need to change this to allow for "reload_selected_datasets")
        if p.method.reprocess_existing_datasets or self.args.method.reprocess_selected_datasets:
            self.log('----------------------------------->>>', True)
            self.log('Setting method.reload_existing_datasets = True')
            self.log('Old (previously processed) datasets will be reloaded')
            p.method.reload_existing_datasets = True

        if self.is_new_pandda() or self.args.method.reprocess_existing_datasets:
            self.log('----------------------------------->>>', True)
            self.log('Setting output.new_analysis_dir = True')
            self.log('A new output analysis directory will be created')
            p.output.new_analysis_dir = True

        self.log('----------------------------------->>>', True)

    def load_pickled_objects(self):
        """Loads any pickled objects it finds"""

        self.log('----------------------------------->>>', True)
        self.log('Looking for pickled files from previous runs in: {!s}'.format(os.path.relpath(self.pickle_handler.get_dir('root'))), True)

        # Record whether any pickled objects are loaded
        pickles_found = False

        # ==============================>
        # Load Grid
        if os.path.exists(self.pickle_handler.get_file('grid')):
            pickles_found = True
            self.log('-> Loading reference grid')
            self.grid = self.unpickle(self.pickle_handler.get_file('grid'))

        # ==============================>
        # Load Reference Dataset
        if os.path.exists(self.pickle_handler.get_file('reference_dataset')):
            pickles_found = True
            self.log('-> Loading reference dataset')
            self.datasets.set_reference(dataset=self.unpickle(self.pickle_handler.get_file('reference_dataset')))

        # ==============================>
        # Load the datasets
        if os.path.exists(self.pickle_handler.get_file('dataset_meta')):
            pickles_found = True
            self.log('-> Loading old dataset information (existing datasets)')
            self.pickled_dataset_meta = self.unpickle(self.pickle_handler.get_file('dataset_meta'))
            if self.args.method.reload_existing_datasets:
                pickled_dataset_list = self.pickled_dataset_meta.dataset_pickle_list
                for filename in pickled_dataset_list:
                    assert os.path.isfile(os.path.join(self.out_dir, filename)), 'File does not exist: {!s}'.format(filename)
                self.log('-> Reloading old datasets')
                self.datasets.add([self.unpickle(os.path.join(self.out_dir,f)) for f in pickled_dataset_list])
            else:
                self.log('-> Not reloading old datasets')
        else:
            # No datasets to load - this must be False
            self.args.method.reload_existing_datasets = False
            self.log('-> No old datasets found')

        # ==============================>
        # Load Statistical Maps
        if os.path.exists(self.pickle_handler.get_file('stat_maps')):
            pickles_found = True
            self.log('-> Loading old statistical maps')
            self.stat_maps = self.unpickle(self.pickle_handler.get_file('stat_maps'))

        if not pickles_found:
            self.log('-> No Pickles Found', True)
        self.log('----------------------------------->>>', True)

    def pickle_the_pandda(self, components=None, all=False, datasets=None):
        """Pickles it's major components for quick loading..."""

        self.log('')
        self.log('################################### <~> ###################################', True)
        self.log('                          pickling pandda objects', True)
        self.log('################################### <~> ###################################', True)
        self.log('')

        if all == True:
            self.log('----------------------------------->>>', True)
            self.log('Pickling the Pandda', True)
        elif not components:
            self.log('----------------------------------->>>', True)
            self.log('Pickling NOTHING', True)
            return
        else:
            self.log('----------------------------------->>>', True)
            self.log('Selective Pickling: {!s}'.format(', '.join(components).upper()), True)

        if all or ('grid' in components):
            self.log('----------------------------------->>>')
            if self.grid is not None:
                self.log('Pickling Map Grid')
                self.pickle(pickle_file=self.pickle_handler.get_file('grid'),
                            pickle_object=self.grid, overwrite=False)
            else:
                self.log('No Reference Grid to Pickle')

        if all or ('datasets' in components):
            self.log('----------------------------------->>>')

            if self.datasets.reference():
                self.log('Pickling Reference Dataset')
                self.pickle(pickle_file=self.pickle_handler.get_file('reference_dataset'),
                            pickle_object=self.datasets.reference().get_pickle_copy(), overwrite=True)

            # If no datasets given, pickle them all
            if not datasets:
                datasets = self.datasets.all()

            # Pickle the datasets (individual pickle files)
            if datasets:
                self.log('Pickling Datasets')
                for d in datasets:
                    self.pickle(pickle_file=d.file_manager.get_file('dataset_pickle'),
                                pickle_object=d.get_pickle_copy(), overwrite=True)
            else:
                self.log('No Datasets to Pickle')

        if all or ('stat_maps' in components):
            self.log('----------------------------------->>>')
            if self.stat_maps is not None:
                self.log('Pickling Statistical Maps')
                self.pickle(pickle_file   = self.pickle_handler.get_file('stat_maps'),
                            pickle_object = self.stat_maps,
                            overwrite = True)
            else:
                self.log('No Statistical Maps to Pickle')

    def exit(self, error=False):
        """Exit the PANDDA, record runtime etc..."""

        self._finish_time = time.time()
        self.log('Runtime: {!s}'.format(time.strftime("%H hours:%M minutes:%S seconds", time.gmtime(self._finish_time - self._init_time))))

        # If error, don't make meta or self pickle
        if error:
            self.update_status('errored')
            self.log('', True)
            self.log('################################### <~> ###################################', True)
            self.log('                       PANDDA exited with an error')
            self.log('################################### <~> ###################################', True)
            self.log('Error Traceback: ')
            self.log('----------------------------------->>>', True)
            self.log(traceback.format_exc())
            self.log('################################### <~> ###################################', True)
            return
        else:
            self.update_status('done')
            self.log('', True)
            self.log('################################### <~> ###################################', True)
            self.log('                  .. FINISHED! PANDDA EXITED NORMALLY ..', True)
            self.log('################################### <~> ###################################', True)
            self.log('', True)

        self.log('----------------------------------->>>', True)
        self.log('Writing dataset information for any future runs', True)
        try:
            # Extract meta about the datasets
            if self.pickled_dataset_meta and (not self.args.method.reload_existing_datasets):
                self.log('Combining old dataset meta with new meta for pickle')
                number_of_datasets  = self.pickled_dataset_meta.number_of_datasets  + self.datasets.size()
                dataset_labels      = self.pickled_dataset_meta.dataset_labels      + [d.tag for d in self.datasets.all()]
                dataset_pickle_list = self.pickled_dataset_meta.dataset_pickle_list + [os.path.relpath(d.file_manager.get_file('dataset_pickle'), start=self.out_dir) for d in self.datasets.all()]
            else:
                self.log('Creating new meta for pickle')
                number_of_datasets  = self.datasets.size()
                dataset_labels      = [d.tag for d in self.datasets.all()]
                dataset_pickle_list = [os.path.relpath(d.file_manager.get_file('dataset_pickle'), start=self.out_dir) for d in self.datasets.all()]
            # Create a dictionary to be stored
            dataset_meta = Meta({'number_of_datasets'    : number_of_datasets,
                                 'dataset_labels'        : dataset_labels,
                                 'dataset_pickle_list'   : dataset_pickle_list})
            # Pickle the list of locations of the dataset pickles
            self.pickle(pickle_file=self.pickle_handler.get_file('dataset_meta'), pickle_object=dataset_meta, overwrite=True)
        except:
            self.log('FAILED TO PICKLE META')

        # Lastly, pickle myself (if required)
        try:
            if self.args.output.pickling.pickle_complete_pandda:
                # Pickle myself
                self.log('----------------------------------->>>', True)
                self.log('Pickling the whole PANDDA object (for developer access)')
                self.pickle(pickle_file=self.pickle_handler.get_file('my_pandda'), pickle_object=self, overwrite=True)
        except:
            self.log('FAILED TO PICKLE MYSELF')

        self.log('', True)

    def is_new_pandda(self):
        """Is this the first time the program has been run?"""
        return self._new_pandda

    def set_high_resolution(self, res):
        self._high_resolution = res
    def get_high_resolution(self):
        return self._high_resolution

    def set_low_resolution(self, res):
        self._low_resolution = res
    def get_low_resolution(self):
        return self._low_resolution

    def new_files(self):
        """Get all of the files that were added on this run"""
        return self._input_files

    def initialise_dataset_masks_and_tables(self):
        """Add blank masks to the mask objects, based on how many datasets have been loaded"""

        self.log('----------------------------------->>>', True)
        self.log('Initialising Dataset Masks.', True)

        # Set the dataset ids (dataset tags)
        self.datasets.all_masks().set_index_ids(ids=self.datasets.all_tags())
        # Initialise standard blank masks
        for mask_name in PanddaMaskNames.all_mask_names:
            self.datasets.all_masks().add_mask(name=mask_name, values=False)

        # Initialise masks for datasets that shouldn't be analysed
        if self.args.input.flags.no_analyse:
            no_analyse_tags = self.args.input.flags.no_analyse.split(',')
            self.log('Not analysing {!s} Datasets: {!s}'.format(len(no_analyse_tags), ', '.join(no_analyse_tags)))
            no_analyse_mask = [True if d.tag in no_analyse_tags else False for d in self.datasets.all()]
            self.datasets.all_masks().add_mask(name='no_analyse', values=no_analyse_mask)

        # Initialise mask for datasets that shouldn't be used for building
        if self.args.input.flags.no_build:
            no_build_tags = self.args.input.flags.no_build.split(',')
            self.log('Not building distributions from {!s} Datasets: {!s}'.format(len(no_build_tags), ', '.join(no_build_tags)))
            no_build_mask = [True if d.tag in no_build_tags else False for d in self.datasets.all()]
            self.datasets.all_masks().add_mask(name='no_build', values=no_build_mask)

        # Initialise mask for datasets that have been previously pickled
        self.datasets.all_masks().add_mask(name='old datasets', values=False)
        if self.pickled_dataset_meta and self.args.method.reload_existing_datasets:
            for tag in self.pickled_dataset_meta.dataset_labels:
                self.datasets.all_masks().set_value(name='old datasets', id=tag, value=True)
            self.log('Considering {!s} datasets as "New Datasets"'.format(self.datasets.size(mask_name='old datasets', invert=True)))
            self.log('Considering {!s} datasets as "Old Datasets"'.format(self.datasets.size(mask_name='old datasets')))
        else:
            self.log('Considering all {!s} datasets as "New Datasets"'.format(self.datasets.size(mask_name='old datasets', invert=True)))
            assert self.datasets.size(mask_name='old datasets', invert=True) == self.datasets.size(), 'Total datasets should be same as total new datasets'

        self.log('----------------------------------->>>', True)
        self.log('Initialising dataset data tables.', True)

        # Add dataset tags as rows in the tables
        self.tables.dataset_info     = self.tables.dataset_info.append(pandas.DataFrame(index=[d.tag for d in self.datasets.all()]), verify_integrity=True)
        self.tables.dataset_map_info = self.tables.dataset_map_info.append(pandas.DataFrame(index=[d.tag for d in self.datasets.all()]), verify_integrity=True)

        old_datasets = self.datasets.mask(mask_name='old datasets')
        if (not self.args.method.reprocess_existing_datasets) and old_datasets:
            self.log('Syncing old dataset information to dataset tables.', True)
            self.sync_datasets(datasets=old_datasets)
            self.log('Syncing old dataset events to output tables.', True)
            for dataset in old_datasets:
                if dataset.events:
                    for e in dataset.events:
                        self.add_event_to_event_table(dataset=dataset, event=e)

    # TODO MOVE TO PANDDA DATASET LIST? TODO
    def select_reference_dataset(self, method='resolution', max_rfree=0.4, min_resolution=5):
        """Select dataset to act as the reference - scaling, aligning etc"""

        assert method in ['resolution','rfree'], 'METHOD FOR SELECTING THE REFERENCE DATASET NOT RECOGNISED: {!s}'.format(method)

        # Create a mask of the datasets that can be selected as the reference dataset
        potential_reference_mask = self.datasets.all_masks().combine_masks(names=['no_build', 'rejected - total'], invert_output=True)
        self.datasets.all_masks().add_mask(name='potential reference datasets', values=potential_reference_mask)
        # Get the potential reference datasets
        filtered_datasets = self.datasets.mask(mask_name='potential reference datasets')
        if not filtered_datasets: raise Failure("Can't select a reference dataset - NO SUITABLE (NON-REJECTED) DATASETS REMAINING")

        self.log('---------->>>', True)
        self.log('Selecting Reference Dataset by: {!s}'.format(method), True)
        if method == 'rfree':
            # Get RFrees of datasets (set to dummy value of 999 if resolution is too high so that it is not selected)
            r_frees = [d.model.input.get_r_rfree_sigma().r_free if (d.data.mtz_object().max_min_resolution()[1] < min_resolution) else 999 for d in filtered_datasets]
            if len(r_frees) == 0: raise Exception('NO DATASETS BELOW RESOLUTION CUTOFF {!s}A - CANNOT SELECT REFERENCE DATASET'.format(min_resolution))
            ref_dataset_index = r_frees.index(min(r_frees))
        elif method == 'resolution':
            # Get Resolutions of datasets (set to dummy value of 999 if r-free is too high so that it is not selected)
            resolns = [d.data.mtz_object().max_min_resolution()[1] if (d.model.input.get_r_rfree_sigma().r_free < max_rfree) else 999 for d in filtered_datasets]
            if len(resolns) == 0: raise Exception('NO DATASETS BELOW RFREE CUTOFF {!s} - CANNOT SELECT REFERENCE DATASET'.format(max_rfree))
            ref_dataset_index = resolns.index(min(resolns))

        reference = filtered_datasets[ref_dataset_index]
        self.log('Reference Selected: {!s}'.format(reference.tag), True)
        self.log('Resolution: {!s}, RFree: {!s}'.format(reference.data.mtz_object().max_min_resolution()[1], reference.model.input.get_r_rfree_sigma().r_free), True)

        return reference.model.filename, reference.data.filename

    def load_reference_dataset(self, ref_pdb, ref_mtz):
        """Set the reference dataset, to which all other datasets will be aligned and scaled"""

        self.log('---------->>>', True)
        self.log('Loading Reference Dataset: {!s}'.format(ref_mtz), True)

        link_ref_pdb = self.file_manager.get_file('reference_structure')
        link_ref_mtz = self.file_manager.get_file('reference_dataset')

        # Remove old links?
        if os.path.abspath(ref_pdb) != os.path.abspath(link_ref_pdb):
            if os.path.exists(link_ref_pdb): os.unlink(link_ref_pdb)
            if os.path.exists(link_ref_mtz): os.unlink(link_ref_mtz)
        # Create links to dataset
        if not os.path.exists(link_ref_pdb): rel_symlink(orig=ref_pdb, link=link_ref_pdb)
        if not os.path.exists(link_ref_mtz): rel_symlink(orig=ref_mtz, link=link_ref_mtz)

        # Create and set reference dataset
        ref_dataset = PanddaReferenceDataset.from_file(model_filename=os.path.relpath(link_ref_pdb, start=self.out_dir),
                                                       data_filename=os.path.relpath(link_ref_mtz, start=self.out_dir)).label(num=-1, tag='reference')

        # Calculate the shift required to move the reference structure into the positive quadrant
        buffer = self.params.masks.outer_mask + self.params.maps.padding
        sites_min = protein(ref_dataset.model.hierarchy).atoms().extract_xyz().min()
        ref_dataset.set_origin_shift(shift=-1*flex.double(sites_min)+buffer)
        self.log('Origin Shift for reference structure: {!s}'.format(tuple([round(s,3) for s in ref_dataset.origin_shift()])))

        # Set as the reference dataset for the analysis
        self.datasets.set_reference(dataset=ref_dataset)

        # Write out the structure in the reference frame
        tmp_r_hierarchy = ref_dataset.model.hierarchy.deep_copy()

        if not os.path.exists(self.file_manager.get_file('reference_on_origin')):
            tmp_r_hierarchy.atoms().set_xyz(ref_dataset.nat2grid(tmp_r_hierarchy.atoms().extract_xyz()))
            tmp_r_hierarchy.write_pdb_file(self.file_manager.get_file('reference_on_origin'))

        if not os.path.exists(self.file_manager.get_file('reference_symmetry')):
            ref_sym_copies = ref_dataset.model.crystal_contacts(distance_cutoff=self.args.params.masks.outer_mask+5, combine_copies=True)
            ref_sym_copies.atoms().set_xyz(ref_dataset.nat2grid(ref_sym_copies.atoms().extract_xyz()))
            ref_sym_copies.write_pdb_file(self.file_manager.get_file('reference_symmetry'))

        return self.datasets.reference()

    def create_reference_grid(self, dataset, grid_spacing):
        """Create a grid over the given dataset"""

        self.log('----------------------------------->>>', True)
        self.log('Creating Reference Grid', True)

        sites_cart = protein(dataset.model.hierarchy).atoms().extract_xyz()
        buffer = self.params.masks.outer_mask + self.params.maps.padding

        grid_min = flex.double([s-buffer for s in sites_cart.min()])
        grid_max = flex.double([s+buffer for s in sites_cart.max()])

        # TODO origin -> grid_min, approx_max -> grid_max TODO
        self.grid = Grid(grid_spacing   = grid_spacing,
                         origin         = (0,0,0),
                         approx_max     = tuple(grid_max-grid_min),
                         verbose        = self.settings.verbose)

        self.log(self.grid.summary())

        return self.grid

    def mask_reference_grid(self, dataset):
        """Create masks for the reference grid based on distances from atoms in the reference structure"""

        self.log('----------------------------------->>>', True)
        self.log('Masking Reference Grid', True)

        # ============================================================================>
        # Get main and neighbouring symmetry copies of the reference structures
        # ============================================================================>
        ref_sites_cart = dataset.model.alignment.nat2ref(protein(dataset.model.hierarchy).atoms().extract_xyz())
        sym_copies = dataset.model.crystal_contacts(distance_cutoff=self.args.params.masks.outer_mask+5, combine_copies=True)
        sym_sites_cart = dataset.model.alignment.nat2ref(protein(sym_copies).atoms().extract_xyz())
        # ============================================================================>
        # Global mask used for removing points in the bulk solvent regions
        # ============================================================================>
        if self.grid.global_mask() is None:
            self.log('---------->>>', True)
            self.log('Generating Protein Mask')
            global_mask = AtomicMask(parent=self.grid, sites_cart=ref_sites_cart,
                                     max_dist=self.params.masks.outer_mask,
                                     min_dist=self.params.masks.inner_mask)
            self.grid.set_global_mask(global_mask)
        # ============================================================================>
        # Global mask used for removing points close to symmetry copies of the protein
        # ============================================================================>
        if self.grid.symmetry_mask() is None:
            self.log('---------->>>', True)
            self.log('Generating Symmetry Mask')
            symmetry_mask = GridMask(parent=self.grid, sites_cart=sym_sites_cart,
                                     max_dist=self.params.masks.outer_mask,
                                     min_dist=self.params.masks.inner_mask_symmetry)
            self.grid.set_symmetry_mask(symmetry_mask)

        # ============================================================================>
        # Write masked maps
        # ============================================================================>
        # Write protein masked map
        write_indices_as_map(grid=self.grid, indices=self.grid.global_mask().total_mask_indices(),
                             f_name=self.file_manager.get_file('reference_dataset').replace('.mtz','.totalmask.ccp4'))
        # Write symmetry masked map
        write_indices_as_map(grid=self.grid, indices=self.grid.symmetry_mask().total_mask_indices(),
                             f_name=self.file_manager.get_file('reference_dataset').replace('.mtz','.symmask.ccp4'))

        return self.grid

    def partition_reference_grid(self, dataset, altlocs=['','A']):

        self.log('----------------------------------->>>', True)
        self.log('Partitioning Reference Grid', True)

        # Select the sites for generating the voronoi alignments (calphas)
        partition_h = calphas(sel_altloc(dataset.model.hierarchy, altlocs=altlocs))
        site_cart_ca = dataset.nat2grid(partition_h.atoms().extract_xyz())

        t1 = time.time()
        self.grid.create_grid_partition(sites_cart=site_cart_ca)
        self.grid.partition.partition(mask  = self.grid.global_mask(),
                                      cpus  = self.settings.cpus)
        t2 = time.time()
        self.log('> Grid partitioning complete > Time Taken: {!s} seconds'.format(int(t2-t1)))

        self.log('----------------------------------->>>', True)
        self.log('Partition Summary:', True)
        self.log('----------------------------------->>>', True)
        voronoi_counts = dict(zip(*numpy.unique(self.grid.partition.nn_groups, return_counts=True)))
        # Cell-by-Cell summary of the voronoi cells
        self.log('--------------------------', True)
        self.log('CHN - RES -  RESID  - ATOM - ALT :    VORONOI VOLUME')
        self.log('--------------------------', True)
        for i_atom, atom in enumerate(partition_h.atoms_with_labels()):
            self.log('{:<3} - {:<3} - {:<7} - {:<4} - {:<3} : {:>10} points'.format(atom.chain_id, atom.resname, atom.resid(), atom.name, atom.altloc, voronoi_counts.get(i_atom,0)))
        self.log('--------------------------', True)
        self.log('Unpartitioned space: {} points'.format(voronoi_counts.get(-1,0)))
        self.log('--------------------------', True)
        # Chain-by-chain summary of the voronoi cells
        for c in partition_h.chains():
            self.log('Chain {:1} - {:5} regions - ({:5} residues)'.format(c.id, len(c.atoms()), len(c.residue_groups())))
        self.log('--------------------------', True)
        self.log('Total: {} regions ({} chains, {} residues)'.format(len(partition_h.atoms()), len(list(partition_h.chains())), len(list(partition_h.residue_groups()))), True)

        # Write grid summary for developer purposes
        if self.args.output.developer.write_grid_masks:
            self.log('----------------------------------->>>', True)
            self.log('Writing Voronoi grid masks:', True)
            # Write out the un-partitioned section of the grid
            write_bool_as_map(grid=self.grid, array=(self.grid.partition.nn_groups==-1),
                              f_name=self.file_manager.get_file('grid_voronoi').format('unpartitioned'))
            # Write out the voronoi masks for each atom
            for i_cell in range(0,10):
                write_bool_as_map(grid=self.grid, array=((self.grid.partition.nn_groups%10)==i_cell)*(self.grid.partition.nn_groups>=0).astype(int),
                                  f_name=self.file_manager.get_file('grid_voronoi').format('{:04}'.format(i_cell)))

            # Write out pymol script to allow results to be access easily
            from bamboo.pymol_utils import PythonScript, Sphere
            pml = PythonScript()
            pml.set_normalise_maps(False)
            for i_atom, atom in enumerate(partition_h.atoms_with_labels()):
                ca_xyz = site_cart_ca[i_atom]
                ca_sph = Sphere(centre=tuple(ca_xyz), radius=min(5.0,max(0.1,0.0002*voronoi_counts.get(i_atom,0))))
                ca_nam = 'voronoi_centres'
                pml.add_shape(shape=ca_sph, obj=ca_nam)
                pml.colour(obj=ca_nam, colour='white')
            pdb = pml.load_pdb(f_name=self.file_manager.get_file('reference_on_origin'))
            pml.colour(obj=pdb, colour='grey')
            for grid_file in sorted(glob.glob(self.file_manager.get_file('grid_voronoi').format('*'))):
                if not grid_file.endswith('.ccp4'): continue
                map_name = pml.load_map(f_name=grid_file)
                mes_name = pml.make_mesh(map_name, contour_level=0)
                pml.colour(mes_name)
            pml.write_script(f_name=self.file_manager.get_file('grid_voronoi').format('cell-centres').replace('.ccp4','.py'))

        return self.grid

    def build_input_list(self):
        """Builds a list of input files from the command line arguments passed"""

        self.log('----------------------------------->>>', True)
        self.log('Building List of Datasets')

        dir_style = self.args.input.data_dirs.strip('./')
        pdb_style = self.args.input.pdb_style.lstrip('/')
        if self.args.input.mtz_style:
            mtz_style = self.args.input.mtz_style.lstrip('/')
        else:
            assert pdb_style.endswith('.pdb'), 'pdb_style does not end in .pdb'
            mtz_style = pdb_style.replace('.pdb','.mtz')

        self.log('Looking for folders that match {}'.format(dir_style))
        self.log('...and for pdb files that match "{}" in each folder'.format(pdb_style))
        self.log('...and for mtz files that match "{}" in each folder'.format(mtz_style))

        # Datasets that are already added
        new_files = []
        empty_directories = []

        for dir in sorted(glob.glob(self.args.input.data_dirs)):
            pdb_files = [f for f in glob.glob(os.path.join(dir, pdb_style)) if os.path.exists(f)]
            mtz_files = [f for f in glob.glob(os.path.join(dir, mtz_style)) if os.path.exists(f)]
            if not (pdb_files and mtz_files):
                print('EMPTY DIRECTORY: {!s}'.format(dir))
                empty_directories.append(dir)
            elif not pdb_files:
                print('NO PDB IN DIRECTORY: {!s}'.format(dir))
                empty_directories.append(dir)
            elif not mtz_files:
                print('NO MTZ IN DIRECTORY: {!s}'.format(dir))
                empty_directories.append(dir)
            else:
                assert len(pdb_files) == 1, 'More than one matching PDB file found: {!s}'.format(os.path.join(dir, pdb_style))
                assert len(mtz_files) == 1, 'More than one matching MTZ file found: {!s}'.format(os.path.join(dir, mtz_style))

                new_pdb = pdb_files[0]
                new_mtz = mtz_files[0]
                dataset_tag = [None]

                # Do regex matching on the file pairs
                if '*' in pdb_style:
                    pdb_base = os.path.basename(new_pdb)
                    if self.args.input.regex.pdb_regex:
                        pdb_regex = self.args.input.regex.pdb_regex
                    else:
                        pdb_regex = pdb_style.replace('*', '(.*)')
                    pdb_tag = re.findall(pdb_regex, pdb_base)
                    assert pdb_tag, 'NO PDB TAG FOUND: {!s} -> {!s}'.format(pdb_regex, pdb_base)
                    if isinstance(pdb_tag[0], tuple):
                        self.log('More than one PDB TAG found - choosing the first one of {!s}'.format(pdb_tag[0]))
                        pdb_tag = list(pdb_tag[0])[0:1]
                else: pdb_regex = pdb_tag = None

                if '*' in mtz_style:
                    mtz_base = os.path.basename(new_mtz)
                    if self.args.input.regex.mtz_regex:
                        mtz_regex = self.args.input.regex.mtz_regex
                    else:
                        mtz_regex = mtz_style.replace('*', '(.*)')
                    mtz_tag = re.findall(mtz_regex, mtz_base)
                    assert mtz_tag, 'NO MTZ TAG FOUND: {!s} -> {!s}'.format(mtz_regex, mtz_base)
                    if isinstance(mtz_tag[0], tuple):
                        self.log('More than one MTZ TAG found - choosing the first one of {!s}'.format(mtz_tag[0]))
                        mtz_tag = list(mtz_tag[0])[0:1]
                else: mtz_regex = mtz_tag = None

                if '*' in dir_style:
                    dir_base = os.path.dirname(pdb_files[0])
                    if self.args.input.regex.dir_regex:
                        dir_regex = self.args.input.regex.dir_regex
                    else:
                        dir_regex = dir_style.replace('*', '(.*)')
                    dir_tag = re.findall(dir_regex, dir_base)
                    assert dir_tag, 'NO DIR TAG FOUND: {!s} -> {!s}'.format(dir_regex, dir_base)
                    if isinstance(dir_tag[0], tuple):
                        self.log('More than one DIR TAG found - choosing the first one of {!s}'.format(dir_tag[0]))
                        dir_tag = list(dir_tag[0])[0:1]
                else: dir_regex = dir_tag = None

                if pdb_tag and mtz_tag: assert pdb_tag == mtz_tag, 'PDB-MTZ TAGS ARE NOT IDENTICAL: {} != {}'.format(pdb_tag, mtz_tag)
                if dir_tag and pdb_tag: assert dir_tag == pdb_tag, 'DIR-PDB TAGS ARE NOT IDENTICAL: {} != {}'.format(dir_tag, pdb_tag)
                if dir_tag and mtz_tag: assert dir_tag == mtz_tag, 'DIR-MTZ TAGS ARE NOT IDENTICAL: {} != {}'.format(dir_tag, mtz_tag)

                if   dir_tag: dataset_tag = dir_tag
                elif pdb_tag: dataset_tag = pdb_tag
                elif mtz_tag: dataset_tag = mtz_tag

                # Add prefix
                if isinstance(dataset_tag[0], str): dataset_tag = [self.args.output.dataset_prefix + dataset_tag[0]]
                else:                               assert dataset_tag[0] is None

                new_files.append(pdb_files+mtz_files+dataset_tag)

        # Filter out the already added files
        if self.pickled_dataset_meta:
            filtered_new_files = []
            for i, (pdb, mtz, tag) in enumerate(new_files):
                if tag in self.pickled_dataset_meta.dataset_labels:
                    self.log('Dataset with this tag has already been loaded: {!s} - Not loading'.format(tag))
                else:
                    filtered_new_files.append(new_files[i])
        else:
            filtered_new_files = new_files

        # Filter out manually labelled datasets to ignore
        if self.args.input.flags.ignore_datasets:
            ignore_tags = self.args.input.flags.ignore_datasets.split(',')
            self.log('Ignoring {!s} Datasets: {!s}'.format(len(ignore_tags), ', '.join(ignore_tags)))
            re_filtered_new_files = []
            for i, (pdb, mtz, tag) in enumerate(filtered_new_files):
                if tag in ignore_tags:
                    self.log('Ignoring Dataset: {!s}'.format(tag))
                else:
                    re_filtered_new_files.append(filtered_new_files[i])
            filtered_new_files = re_filtered_new_files

#        # Get the list of already linked empty_directories
#        empty_dir_prefix = 'Dir_'
#        link_old_empty_dirs = glob.glob(os.path.join(self.file_manager.get_dir('empty_directories'),'*'))
#        real_old_empty_dirs = [os.path.realpath(p) for p in link_old_empty_dirs]
#        # Pull out the highest current idx
#        emp_num_i = 0
#        emp_num_offset = max([0]+[int(os.path.basename(v).strip(empty_dir_prefix)) for v in link_old_empty_dirs])
#        assert emp_num_offset == len(link_old_empty_dirs), 'Numbering of empty directories is not consecutive'
#        # Link the empty directories into the same directory
#        for dir in empty_directories:
#            # Already linked
#            if os.path.realpath(dir) in real_old_empty_dirs:
#                continue
#            # Increment counter
#            emp_num_i += 1
#            # Create new dir
#            empty_dir = os.path.join(self.file_manager.get_dir('empty_directories'), empty_dir_prefix+'{:05d}'.format(emp_num_i+emp_num_offset))
#            if not os.path.exists(empty_dir):
#                os.symlink(dir, empty_dir)
#            else:
#                raise Exception('THIS DIRECTORY SHOULD NOT EXIST: {!s}'.format(empty_dir))

        # Record number of empty datasets
        self.log('---------->>>', True)
        self.log('{!s} EMPTY DIRECTORIES FOUND:'.format(len(empty_directories)), True)
        self.log('---------->>>')
        for d in empty_directories:
            self.log('Empty Directory: {}'.format(d))
#        self.log('{!s} EMPTY DIRECTORIES FOUND (TOTAL)'.format(emp_num_i+emp_num_offset), True)
#        self.log('{!s} EMPTY DIRECTORIES FOUND (NEW)'.format(emp_num_i), True)

        # Record total number of datasets, and total number of new datasets
        if self.pickled_dataset_meta: num_old = self.pickled_dataset_meta.number_of_datasets
        else:                         num_old = 0
        self.log('---------->>>', True)
        self.log('{!s} DATASETS FOUND (TOTAL)'.format(len(filtered_new_files)+num_old), True)
        self.log('{!s} DATASETS FOUND (NEW)'.format(len(filtered_new_files)), True)

        return filtered_new_files

    def add_new_files(self, input_files):
        """Add (pdb, mtz) file pairs to the datasets to be processed"""

        self._input_files = input_files

        self.log('----------------------------------->>>', True)
        self.log('{!s} Datasets Added'.format(len(input_files)), True)

    def load_new_datasets(self):
        """Read in maps for the input datasets"""

        if not self.datasets.all() and self.is_new_pandda():
            self.log('Adding First Datasets to Pandda')
        else:
            self.log('Adding more datasets')
            self.log('{!s} already loaded'.format(self.datasets.size()))
            self.log('{!s} not loaded'.format(self.pickled_dataset_meta.number_of_datasets - self.datasets.size()))

        # Counting offset for dataset index
        if self.pickled_dataset_meta: n_offset = self.pickled_dataset_meta.number_of_datasets
        else:                         n_offset = 0


        start = time.time()
        self.log('----------------------------------->>>', True)
        print('Loading Datasets...')
        loaded_datasets = [PanddaDataset.from_file(model_filename=pdb, data_filename=mtz).label(num=num+n_offset, tag=dtag) for num, (pdb, mtz, dtag) in enumerate(self.new_files())]
        finish = time.time()
        self.log('> Adding Datasets > Time Taken: {!s} seconds'.format(int(finish-start)), True)
        self.log('----------------------------------->>>', True)

        lig_style = self.args.input.lig_style.strip('/')

        # Output Path Templates
        f = PanddaDatasetFilenames
        p = PanddaDatasetPNGFilenames

        for dataset in loaded_datasets:

            # Intialise the meta for the dataset
            dataset.meta.analysed = False
            dataset.meta.dataset_info = None
            dataset.meta.dataset_map_info = None

            # Create a file manager object
            dataset.initialise_output_directory(dir=os.path.join(self.file_manager.get_dir('processed_datasets'), dataset.tag))

            # Main input/output files
            dataset.file_manager.add_file(file_name=f.input_structure.format(dataset.tag),                    file_tag='input_structure'              )
            dataset.file_manager.add_file(file_name=f.input_data.format(dataset.tag),                         file_tag='input_data'                   )
            dataset.file_manager.add_file(file_name=f.dataset_info.format(dataset.tag),                       file_tag='dataset_info'                 )
            dataset.file_manager.add_file(file_name=f.dataset_log.format(dataset.tag),                        file_tag='dataset_log'                  )
            dataset.file_manager.add_file(file_name=f.aligned_structure.format(dataset.tag),                  file_tag='aligned_structure'            )
            dataset.file_manager.add_file(file_name=f.symmetry_copies.format(dataset.tag),                    file_tag='symmetry_copies'              )
            dataset.file_manager.add_file(file_name=f.sampled_map.format(dataset.tag),                        file_tag='sampled_map'                  )
            dataset.file_manager.add_file(file_name=f.mean_diff_map.format(dataset.tag),                      file_tag='mean_diff_map'                )
            dataset.file_manager.add_file(file_name=f.z_map.format(dataset.tag),                              file_tag='z_map'                        )
            dataset.file_manager.add_file(file_name=f.z_map_naive.format(dataset.tag),                        file_tag='z_map_naive'                  )
            dataset.file_manager.add_file(file_name=f.z_map_naive_norm.format(dataset.tag),                   file_tag='z_map_naive_normalised'       )
            dataset.file_manager.add_file(file_name=f.z_map_uncertainty.format(dataset.tag),                  file_tag='z_map_uncertainty'            )
            dataset.file_manager.add_file(file_name=f.z_map_uncertainty_norm.format(dataset.tag),             file_tag='z_map_uncertainty_normalised' )
            dataset.file_manager.add_file(file_name=f.z_map_corrected.format(dataset.tag),                    file_tag='z_map_corrected'              )
            dataset.file_manager.add_file(file_name=f.z_map_corrected_norm.format(dataset.tag),               file_tag='z_map_corrected_normalised'   )
            dataset.file_manager.add_file(file_name=f.event_map.format(dataset.tag, '{!s}', '{!s}'),          file_tag='event_map'                    )

            # Miscellaneous files
            dataset.file_manager.add_file(file_name=f.high_z_mask.format(dataset.tag), file_tag='high_z_mask')
            dataset.file_manager.add_file(file_name=f.grid_mask.format(dataset.tag),   file_tag='grid_mask')

            # Links to ligand files (if they've been found)
            dataset.file_manager.add_dir(dir_name='ligand_files', dir_tag='ligand', top_dir_tag='root')
#            dataset.file_manager.add_file(file_name=f.ligand_coordinates.format(dataset.tag),                 file_tag='ligand_coordinates'   )
#            dataset.file_manager.add_file(file_name=f.ligand_restraints.format(dataset.tag),                  file_tag='ligand_restraints'    )
#            dataset.file_manager.add_file(file_name=f.ligand_image.format(dataset.tag),                       file_tag='ligand_image'         )

            # Native (back-rotated/transformed) maps
            dataset.file_manager.add_file(file_name=f.native_obs_map.format(dataset.tag),                     file_tag='native_obs_map'       )
            dataset.file_manager.add_file(file_name=f.native_z_map.format(dataset.tag),                       file_tag='native_z_map'         )
            dataset.file_manager.add_file(file_name=f.native_event_map.format(dataset.tag,'{!s}','{!s}'),     file_tag='native_event_map'     )
            dataset.file_manager.add_file(file_name=f.native_mean_map.format(dataset.tag),                    file_tag='native_mean_map'      )

            # Fitted structures when modelled with pandda.inspect
            dataset.file_manager.add_dir(dir_name='modelled_structures', dir_tag='models', top_dir_tag='root')

            # Output images
            dataset.file_manager.add_dir(dir_name='output_images', dir_tag='images', top_dir_tag='root')
            # Smapled map
            dataset.file_manager.add_file(file_name=p.s_map_png.format(dataset.tag),                          file_tag='s_map_png',                        dir_tag='images')
            dataset.file_manager.add_file(file_name=p.d_mean_map_png.format(dataset.tag),                     file_tag='d_mean_map_png',                   dir_tag='images')
            dataset.file_manager.add_file(file_name=p.z_map_naive_png.format(dataset.tag),                    file_tag='z_map_naive_png',                  dir_tag='images')
            dataset.file_manager.add_file(file_name=p.z_map_naive_norm_png.format(dataset.tag),               file_tag='z_map_naive_normalised_png',       dir_tag='images')
            dataset.file_manager.add_file(file_name=p.z_map_uncertainty_png.format(dataset.tag),              file_tag='z_map_uncertainty_png',            dir_tag='images')
            dataset.file_manager.add_file(file_name=p.z_map_uncertainty_norm_png.format(dataset.tag),         file_tag='z_map_uncertainty_normalised_png', dir_tag='images')
            dataset.file_manager.add_file(file_name=p.z_map_corrected_png.format(dataset.tag),                file_tag='z_map_corrected_png',              dir_tag='images')
            dataset.file_manager.add_file(file_name=p.z_map_corrected_norm_png.format(dataset.tag),           file_tag='z_map_corrected_normalised_png',   dir_tag='images')
            dataset.file_manager.add_file(file_name=p.z_map_qq_plot_png.format(dataset.tag),                  file_tag='z_map_qq_plot_png',                dir_tag='images')
            dataset.file_manager.add_file(file_name=p.bdc_est_png.format(dataset.tag, '{!s}'),                file_tag='bdc_est_png',                     dir_tag='images')
            dataset.file_manager.add_file(file_name=p.unc_qqplot_png.format(dataset.tag),                     file_tag='unc_qqplot_png',                   dir_tag='images')
            dataset.file_manager.add_file(file_name=p.obs_qqplot_sorted_png.format(dataset.tag),              file_tag='obs_qqplot_sorted_png',            dir_tag='images')
            dataset.file_manager.add_file(file_name=p.obs_qqplot_unsorted_png.format(dataset.tag),            file_tag='obs_qqplot_unsorted_png',          dir_tag='images')

            # Analysis files
            dataset.file_manager.add_file(file_name=f.z_peaks_csv.format(dataset.tag), file_tag='z_peaks_csv')

            # Scripts
            dataset.file_manager.add_dir(dir_name='scripts', dir_tag='scripts', top_dir_tag='root')
            dataset.file_manager.add_file(file_name=f.pymol_script,     file_tag='pymol_script',    dir_tag='scripts')
            dataset.file_manager.add_file(file_name=f.ccp4mg_script,    file_tag='ccp4mg_script',   dir_tag='scripts')

            # Output blobs
            dataset.file_manager.add_dir(dir_name='blobs', dir_tag='blobs', top_dir_tag='root')
            dataset.file_manager.add_file(file_name=f.ccp4mg_png,       file_tag='ccp4mg_png',      dir_tag='blobs')

            # Pickled objects
            dataset.file_manager.add_dir(dir_name='pickles', dir_tag='pickles', top_dir_tag='root')
            dataset.file_manager.add_file(file_name=f.dataset_pickle,   file_tag='dataset_pickle',  dir_tag='pickles')

            ##############################################################################################################

            # Links for the dataset input files
            link_pdb = dataset.file_manager.get_file('input_structure')
            link_mtz = dataset.file_manager.get_file('input_data')
            # Link the input files to the output folder
            if not os.path.exists(link_pdb): rel_symlink(orig=dataset.model.filename, link=link_pdb)
            if not os.path.exists(link_mtz): rel_symlink(orig=dataset.data.filename, link=link_mtz)

            # Search for ligand files and link them to the output ligands folder
            lig_glob  = os.path.join(os.path.dirname(dataset.model.filename), lig_style)
            lig_files = glob.glob(lig_glob)
            for lig_file in lig_files:
                # Find all files with the same basename but allowing for different extensions. Then link to output folder.
                lig_base = os.path.splitext(lig_file)[0] + '.*'
                lig_matches = glob.glob(lig_base)
                for lig in lig_matches:
                    out_path = os.path.join(dataset.file_manager.get_dir('ligand'), os.path.basename(lig))
                    if os.path.exists(lig) and (not os.path.exists(out_path)):
                        try: shutil.copy(lig, out_path)
                        except: pass

            # Lastly: Update the pointer to the new path (relative to the pandda directory)
            dataset.model.filename = os.path.relpath(link_pdb, start=self.out_dir)
            dataset.data.filename = os.path.relpath(link_mtz, start=self.out_dir)

        self.datasets.add(loaded_datasets)
        self.log('{!s} Datasets Loaded (New).          '.format(len(loaded_datasets), True))
        self.log('{!s} Datasets Loaded (Total).        '.format(self.datasets.size(), True))

    #########################################################################################################
    #                                                                                                       #
    #                                     Dataset loading / processing                                      #
    #                                                                                                       #
    #########################################################################################################

    def load_diffraction_data(self, datasets=None):
        """Extract amplitudes and phases for creating map"""

        if datasets is None:
            datasets=self.datasets.mask(mask_name='rejected - total', invert=True)

        all_cols = self.args.params.maps.structure_factors
        ref_cols = self.args.input.reference.structure_factors if self.args.input.reference.structure_factors else all_cols

        t1 = time.time()
        self.log('----------------------------------->>>', True)
        self.log('Loading structure factors: {}'.format(all_cols), True)
        self.log('----------------------------------->>>', True)

        self.datasets.reference().data.get_structure_factors(columns=ref_cols)

        for dataset in datasets:
            if all_cols in dataset.data.structure_factors.keys():
                self.log('\rAlready Loaded: Dataset {!s}                           '.format(dataset.tag))
            else:
                print('\rLoading Diffraction Data: Dataset {!s}                 '.format(dataset.tag), end=''); sys.stdout.flush()
                dataset.data.get_structure_factors(columns=all_cols)

        t2 = time.time()
        self.log('\r> Structure Factors Extracted > Time Taken: {!s} seconds'.format(int(t2-t1)), True)

        return all_cols,ref_cols

    def align_datasets(self, method):
        """Align each structure the reference structure"""

        assert method in ['local','global'], 'METHOD NOT DEFINED: {!s}'.format(method)

        # Select the datasets for alignment
        datasets_for_alignment = [d for d in self.datasets.mask(mask_name='rejected - total', invert=True) if not d.model.alignment]

        self.log('', True)
        self.log('----------------------------------->>>', True)

        if not datasets_for_alignment:
            self.log('All datasets are already aligned/No datasets to align')
            return
        else:
            print('Generating Alignments (using {!s} cores) for {} datasets'.format(self.settings.cpus, len(datasets_for_alignment)))

        # Create a shifted version of the reference model for alignment
        ref_model = copy.deepcopy(self.datasets.reference().model)
        ref_model.hierarchy.atoms().set_xyz(ref_model.alignment.nat2ref(ref_model.hierarchy.atoms().extract_xyz()))
        # Generate the alignments for each structure
        t1 = time.time()
        arg_list = [DatasetAligner(model=d.model, other=ref_model, method=method, id=d.tag) for d in datasets_for_alignment]
        dataset_alignments = easy_mp.pool_map(func=wrapper_run, args=arg_list, processes=self.settings.cpus)
        t2 = time.time()
        # Post-process the alignments and output the aligned structures
        self.log('----------------------------------->>>', True)
        self.log('Alignment Summaries:', True)
        for alignment in dataset_alignments:
            # Get the master copy of the dataset
            dataset = self.datasets.get(tag=alignment.id)
            dataset.model.alignment = alignment
            # Output an aligned copy of the structure
            aligned_struc = dataset.model.hierarchy.deep_copy()
            aligned_struc.atoms().set_xyz(dataset.model.alignment.nat2ref(coordinates=dataset.model.hierarchy.atoms().extract_xyz()))
            aligned_struc.write_pdb_file(file_name=os.path.join(self.file_manager.get_dir('aligned_structures'), '{!s}-aligned.pdb'.format(dataset.tag)))
            aligned_struc.write_pdb_file(file_name=dataset.file_manager.get_file('aligned_structure'))
            # Write alignment summary to log
            self.log('------------------------->>>')
            self.log(dataset.model.alignment.summary())

        t3 = time.time()
        self.log('----------------------------------->>>', True)
        self.log('> Generating Alignments > Time Taken: {!s} seconds'.format(int(t2-t1)), True)
        self.log('> Aligning Structures   > Time Taken: {!s} seconds'.format(int(t3-t2)), True)
        self.log('----------------------------------->>>', True)

    def collate_dataset_variables(self):
        """Go through all of the datasets and collect lots of different characteristics of the datasets for identifying odd datasets"""

        self.log('----------------------------------->>>', True)
        self.log('Collating Dataset Structure/Crystal Variables', True)

        for d in self.datasets.all():
            # Resolution info
            self.tables.dataset_info.set_value(d.tag, 'high_resolution', numpy.round(d.data.summary.high_res,3))
            self.tables.dataset_info.set_value(d.tag, 'low_resolution',  numpy.round(d.data.summary.low_res,3))
            # Unit cell info
            self.tables.dataset_info.set_value(d.tag, ['uc_a','uc_b','uc_c','uc_alpha','uc_beta','uc_gamma'],   numpy.round(d.data.summary.unit_cell.parameters(),3))
            self.tables.dataset_info.set_value(d.tag, 'uc_vol',                                                 numpy.round(d.data.summary.unit_cell.volume()),3)
            # Spacegroup info
            self.tables.dataset_info.set_value(d.tag, 'space_group', d.data.summary.space_group.info().type().lookup_symbol())
            # Quality info
            self.tables.dataset_info.set_value(d.tag, 'r_work', round_no_fail(d.model.input.get_r_rfree_sigma().r_work,3))
            self.tables.dataset_info.set_value(d.tag, 'r_free', round_no_fail(d.model.input.get_r_rfree_sigma().r_free,3))
            # Alignment info
            if d.model.alignment:
                self.tables.dataset_info.set_value(d.tag, 'rmsd_to_reference', numpy.round(d.model.alignment.alignment_rmsd(),3))

        self.log('----------------------------------->>>', True)

    #########################################################################################################
    #                                                                                                       #
    #                                 Dataset variation analysis functions                                  #
    #                                                                                                       #
    #########################################################################################################

    def analyse_unit_cell_variation(self):
        pass

    def analyse_alignment_variation(self):
        """Look at all of the rotation matrices for the local alignments and calculate the rms between neighbours"""

        # TODO TODO TODO
        return
        # TODO TODO TODO

        assert self.params.alignment.method in ['global', 'local']

        if self.params.alignment.method == 'global':
            self.log('GLOBAL ALIGNMENT SELECTED - NOT ANALYSING ROTATION MATRICES')
            return

        if self.settings.plot_graphs:
            import matplotlib
            matplotlib.interactive(False)
            from matplotlib import pyplot

        # Select datasets to analyse
        used_datasets = self.datasets.mask(mask_name='rejected - total', invert=True)

        ref_calpha_atoms = calphas(self.datasets.reference().model.hierarchy).atoms()
        ref_calpha_sites = ref_calpha_atoms.extract_xyz()
        ref_calpha_labels = [make_label(a) for a in ref_calpha_atoms]

        # Array to hold the output data
        num_datasets = len(used_datasets)
        num_pairs =  len(ref_c_alpha_labels)-1
        output_diffs = numpy.zeros((num_datasets, num_pairs, 2))

        # Iterate through the datasets and pull out the alignment matrices
        for d_num, dataset in enumerate(used_datasets):
            # Extract and sort dataset alignments
            alignments = dataset.local_alignment_transforms()
            alignment_keys = sorted(alignments.keys())
            assert alignment_keys == ref_c_alpha_labels

            # Iterate through adjacent pairs of matrices
            for i in range(0, num_pairs):
                # Label and lsq fit for the current calpha
                calpha_1 = alignment_keys[i]
                rt_1 = alignments[calpha_1]
                # And for the next calpha
                calpha_2 = alignment_keys[i+1]
                rt_2 = alignments[calpha_2]

                assert calpha_1 == ref_c_alpha_labels[i]
                assert calpha_2 == ref_c_alpha_labels[i+1]

                # Calculate the difference in the angles of the alignment matrices
                theta_1 = scitbx.math.math.acos((rt_1.r.trace()-1)/2.0)
                theta_2 = scitbx.math.math.acos((rt_2.r.trace()-1)/2.0)
                # XXX Should we calculate the absolute of the difference?
                theta_rad = theta_2-theta_1
                theta_deg = theta_rad * 180.0/scitbx.math.math.pi
                # Calculate the difference in the translation
                t_shift = (rt_2.t-rt_1.t).norm_sq()**0.5

#                # Calculate the angles from the multiplication of one by the inverse of the other
#                rt_1_2 = rt_1 * rt_2.inverse()
#                # Calculate the angle of the rotation matrix
#                theta_rad = scitbx.math.math.acos((rt_1_2.r.trace()-1)/2.0)
#                theta_deg = theta_rad * 180.0/scitbx.math.math.pi
#                # Calculate the length of the shift
#                t_shift =  rt_1_2.t.norm_sq()**0.5

                # Append to the array
                output_diffs[d_num, i, :] = theta_deg, t_shift

        # Directory to write the output to
        var_out_dir = self.file_manager.get_dir('analyses')
        # Write out to file
        numpy.savetxt(  fname = os.path.join(var_out_dir, 'calpha_rt_r_variation.csv'), X=output_diffs[:,:,0], delimiter=',', newline='\n' )
        numpy.savetxt(  fname = os.path.join(var_out_dir, 'calpha_rt_t_variation.csv'), X=output_diffs[:,:,1], delimiter=',', newline='\n' )

        # Write out graphs
        if self.settings.plot_graphs:

            # Create labels
            labels = ['']*num_pairs
            for i in range(0, num_pairs, 5)+[num_pairs-1]:
                labels[i] = i+1
            # Clear the last n before the last one
            n = 4
            labels[-1-n:-1] = ['']*n

            # BOX PLOT OF ROTATION AND TRANSLATION SHIFTS
            fig = pyplot.figure()
            pyplot.title('Rotation-translation alignment matrix variation between adjacent C-alpha')
            # ADJACENT ANGLE VARIATION
            pyplot.subplot(2, 1, 1)
            pyplot.boxplot(x=output_diffs[:,:,0], notch=True, sym='.', widths=0.5, whis=[5,95], whiskerprops={'ls':'-'}, flierprops={'ms':1}, labels=labels) # whis='range'
            pyplot.xlabel('C-alpha index')
            pyplot.ylabel('Angle Difference\n(degrees)')
            # ADJACENT SHIFT VARIATION
            pyplot.subplot(2, 1, 2)
            pyplot.boxplot(x=output_diffs[:,:,1], notch=True, sym='.', widths=0.5, whis=[5,95], whiskerprops={'ls':'-'}, flierprops={'ms':1}, labels=labels) # whis='range'
            pyplot.xlabel('C-alpha Index')
            pyplot.ylabel('Translation Difference\n(angstroms)')
            # Apply tight layout to prevent overlaps
            fig.set_tight_layout(True)
            # Save both
            pyplot.savefig(os.path.join(var_out_dir, 'calpha_rt_variation.png'), format='png')
            pyplot.close(fig)

    #########################################################################################################
    #                                                                                                       #
    #                                             Dataset filtering                                         #
    #                                                                                                       #
    #########################################################################################################

    def filter_datasets_1(self, filter_dataset=None):
        """Filter out the datasets which contain different protein models (i.e. protein length, sequence, etc)"""

        self.log('', True)
        self.log('----------------------------------->>>', True)
        self.log('Filtering Datasets (before alignment). Potential Classes:', True)
        for failure_class in PanddaMaskNames.reject_mask_names:
            self.log('\t{!s}'.format(failure_class), True)
        self.log('----------------------------------->>>', True)

        # If no filtering dataset given, filter against the reference dataset
        if not filter_dataset: filter_dataset = self.datasets.reference()

        # Check that the same protein structure is present in each dataset - THIS MASK SHOULD INCLUDE ALL DATASETS AT FIRST
        for dataset in self.datasets.mask(mask_name='rejected - total', invert=True):

            print('\rFiltering Dataset {!s}          '.format(dataset.tag), end=''); sys.stdout.flush()
            if self.params.filtering.flags.same_space_group_only and (dataset.model.space_group.info().symbol_and_number() != filter_dataset.model.space_group.info().symbol_and_number()):
                self.log('\rRejecting Dataset: {!s}          '.format(dataset.tag))
                self.log('Different Space Group')
                self.log('Reference: {!s}, {!s}: {!s}'.format(filter_dataset.model.space_group.info().symbol_and_number(),
                                                        dataset.tag, dataset.model.space_group.info().symbol_and_number()))
                self.log('---------->>>', True)
                self.tables.dataset_info.set_value(dataset.tag, 'rejection_reason', 'Different Space Group to Reference')
                self.datasets.all_masks().set_value(name='rejected - different space group', id=dataset.tag, value=True)
            elif dataset.model.input.get_r_rfree_sigma().r_free > self.params.filtering.max_rfree:
                self.log('\rRejecting Dataset: {!s}          '.format(dataset.tag))
                self.log('RFree is higher than cutoff: {!s}'.format(self.params.filtering.max_rfree))
                self.log('High RFree: {!s}'.format(dataset.model.input.get_r_rfree_sigma().r_free))
                self.log('---------->>>', True)
                self.tables.dataset_info.set_value(dataset.tag, 'rejection_reason', 'R-free is too high')
                self.datasets.all_masks().set_value(name='rejected - rfree', id=dataset.tag, value=True)
            elif self.params.filtering.flags.similar_models_only and (not dataset.model.hierarchy.is_similar_hierarchy(filter_dataset.model.hierarchy)):
                self.log('\rRejecting Dataset: {!s}          '.format(dataset.tag))
                self.log('Non-Identical Structure (Structures do not contain the same atoms)')
                self.log('---------->>>', True)
                self.tables.dataset_info.set_value(dataset.tag, 'rejection_reason', 'Atoms present in the dataset are different to atoms present in the reference structure')
                self.datasets.all_masks().set_value(name='rejected - non-identical structures', id=dataset.tag, value=True)
            else:
                pass

        # Update the combined masks
        combined_reject_mask = self.datasets.all_masks().combine_masks(names=PanddaMaskNames.reject_mask_names)
        self.datasets.all_masks().add_mask(name='rejected - total', values=combined_reject_mask, overwrite=True)

        self.log('\rDatasets Filtered.                          ', True)
        self.log('----------------------------------->>>')
        self.log('Rejected Datasets (Total):     {!s}'.format(sum(self.datasets.all_masks().get_mask(name='rejected - total'))), True)
        self.log('----------------------------------->>>')

        reject_reasons = self.tables.dataset_info['rejection_reason'].value_counts().sort_index()
        if reject_reasons.any():
            self.log('Reasons for Rejection:')
            for reason, count in reject_reasons.iteritems():
                self.log('{} Dataset(s) - {}'.format(count, reason))

#        # Link all rejected datasets into the rejected directory
#        for dataset in self.datasets.mask(mask_name='rejected - total'):
#            reject_dir = os.path.join(self.file_manager.get_dir('rejected_datasets'), dataset.tag)
#            if not os.path.exists(reject_dir):
#                rel_symlink(orig=dataset.file_manager.get_dir('root'), link=reject_dir)

    def filter_datasets_2(self):
        """Filter out the non-isomorphous datasets"""

        self.log('', True)
        self.log('----------------------------------->>>', True)
        self.log('Filtering Datasets (after alignment). Potential Classes:', True)
        for failure_class in PanddaMaskNames.reject_mask_names:
            self.log('\t{!s}'.format(failure_class), True)
        self.log('----------------------------------->>>', True)

        # Check that each dataset is similar enough to be compared
        for dataset in self.datasets.mask(mask_name='rejected - total', invert=True):

            print('\rFiltering Dataset {!s}          '.format(dataset.tag), end=''); sys.stdout.flush()
            # Check the deviation from the average sites
            if dataset.model.alignment.alignment_rmsd() > self.params.filtering.max_rmsd_to_reference:
                self.log('\rRejecting Dataset: {!s}          '.format(dataset.tag))
                self.log('Alignment RMSD is too large')
                self.log('Aligned (Calpha) RMSD: {!s}'.format(dataset.model.alignment.alignment_rmsd()))
                self.log('---------->>>', True)
                self.tables.dataset_info.set_value(dataset.tag, 'rejection_reason', 'High RMSD to aligned reference structure')
                self.datasets.all_masks().set_value(name='rejected - rmsd to reference', id=dataset.tag, value=True)
            else:
                pass

        # Combine all of the masks
        combined_reject_mask = self.datasets.all_masks().combine_masks(names=PanddaMaskNames.reject_mask_names)
        self.datasets.all_masks().add_mask(name='rejected - total', values=combined_reject_mask, overwrite=True)

        self.log('\rDatasets Filtered.               ', True)
        self.log('----------------------------------->>>')
        self.log('Rejected Datasets (Total):     {!s}'.format(sum(self.datasets.all_masks().get_mask(name='rejected - total'))), True)

        reject_reasons = self.tables.dataset_info['rejection_reason'].value_counts().sort_index()
        if reject_reasons.any():
            self.log('Reasons for Rejection:')
            for reason, count in reject_reasons.iteritems():
                self.log('{} Dataset(s) - {}'.format(count, reason))

#        # Link all rejected datasets into the rejected directory
#        for dataset in self.datasets.mask(mask_name='rejected - total'):
#            reject_dir = os.path.join(self.file_manager.get_dir('rejected_datasets'), dataset.tag)
#            if not os.path.exists(reject_dir):
#                rel_symlink(orig=dataset.file_manager.get_dir('root'), link=reject_dir)

    #########################################################################################################
    #                                                                                                       #
    #                              Dataset utility functions (e.g. reset/sync)                              #
    #                                                                                                       #
    #########################################################################################################

    def reset_loaded_datasets(self):
        """Check that pickled datasets are ready for reprocessing, etc, if required"""

        if self.args.method.reprocess_existing_datasets:   datasets_for_reprocessing = self.datasets.all()
        elif self.args.method.reprocess_selected_datasets: datasets_for_reprocessing = [self.datasets.get(tag=t) for t in self.args.method.reprocess_selected_datasets.split(',')]
        else:                                              datasets_for_reprocessing = []

        for dataset in datasets_for_reprocessing:
            # Reset the meta objects
            dataset.meta.analysed = False
            dataset.meta.dataset_map_info = None
            # Delete events from before
            dataset.events = []
            # Reset the map information for the dataset
            self.tables.dataset_map_info.loc[dataset.tag] = numpy.nan

    def check_loaded_datasets(self, datasets):
        """Check that the datasets are analysable (have the right mtz columns, etc)"""

        self.log('----------------------------------->>>', True)
        self.log('Performing checks on the loaded datasets', True)

        self.log('Checking for structure factors:', True)
        sf_cols = self.params.maps.structure_factors.split(',')
        self.log('> '+'\n> '.join(sf_cols), True)
        for dataset in datasets:
            for c in sf_cols:
                if not dataset.data.mtz_object().has_column(c):
                    raise Sorry('Structure factor column "{}" was not found in the reflection data for dataset {}. You may need to change the pandda.params.maps.structure_factors option.'.format(c, dataset.tag))

    def sync_datasets(self, datasets=None, overwrite_dataset_meta=False):
        """Sync the loaded datasets and the pandda dataset tables"""

        if not datasets: datasets = self.datasets.all()

        for dataset in datasets:
            # Copy data from pandda dataset tables to dataset
            if (dataset.meta.dataset_info is None) or overwrite_dataset_meta:
                dataset.meta.dataset_info = self.tables.dataset_info.loc[dataset.tag]
            # Copy data from dataset to pandda dataset tables
            else:
                for col,val in dataset.meta.dataset_info.iteritems():
                    self.tables.dataset_info.set_value(index=dataset.tag, col=col, value=val)

            # Copy data from pandda dataset tables to dataset
            if (dataset.meta.dataset_map_info is None) or overwrite_dataset_meta:
                dataset.meta.dataset_map_info = self.tables.dataset_map_info.loc[dataset.tag]
            # Copy data from dataset to pandda dataset tables
            else:
                for col,val in dataset.meta.dataset_map_info.iteritems():
                    self.tables.dataset_map_info.set_value(index=dataset.tag, col=col, value=val)

    #########################################################################################################
    #                                                                                                       #
    #                                Analysis functions (Dataset selection)                                 #
    #                                                                                                       #
    #########################################################################################################

    def select_for_building_distributions(self, high_res_cutoff):
        """Select all datasets with resolution better than high_res_cutoff"""

        building_mask_name = 'characterisation @ {!s}A'.format(high_res_cutoff)

        # Create empty mask
        self.datasets.all_masks().add_mask(name=building_mask_name, values=False)
        # Counter for the number of datasets to select
        total_build = 0
        # Select from the datasets that haven't been rejected
        for dataset in self.datasets.mask(mask_name='rejected - total', invert=True):
            # Check the resolution of the dataset
            if self.tables.dataset_info.get_value(index=dataset.tag, col='high_resolution') > high_res_cutoff:
                continue
            # Check to see if this has been excluded from building
            elif self.datasets.all_masks().get_value(name='no_build', id=dataset.tag) == True:
                self.log('Rejecting Dataset {!s}: Excluded from building'.format(dataset.tag))
                continue
            else:
                self.datasets.all_masks().set_value(name=building_mask_name, id=dataset.tag, value=True)
                # Check to see if the number of datasets to use in building has been reached
                total_build += 1
                if total_build >= self.params.analysis.max_build_datasets:
                    self.log('Maximum number of datasets for building reached: {!s}={!s}'.format(total_build, self.params.analysis.max_build_datasets))
                    break
        return building_mask_name, self.datasets.all_masks().get_mask(name=building_mask_name)

    def select_for_analysis(self, high_res_large_cutoff, high_res_small_cutoff):
        """Select all datasets with resolution between high and low limits"""

        assert high_res_large_cutoff > high_res_small_cutoff, '{!s} must be larger than {!s}'.format(high_res_large_cutoff, high_res_small_cutoff)

        analysis_mask_name = 'analysis @ {!s}A'.format(high_res_large_cutoff)

        if self.args.method.reprocess_selected_datasets: datasets_for_reprocessing = self.args.method.reprocess_selected_datasets.split(',')
        else:                                            datasets_for_reprocessing = []

        # Create empty mask
        self.datasets.all_masks().add_mask(name=analysis_mask_name, values=False)
        # Select from the datasets that haven't been rejected
        for dataset in self.datasets.mask(mask_name='rejected - total', invert=True):
            # Check the resolution of the dataset (is not too low)
            if self.tables.dataset_info.get_value(index=dataset.tag, col='high_resolution') > high_res_large_cutoff:
                continue
            # Check the resolution of the dataset (is not too high)
            elif self.tables.dataset_info.get_value(index=dataset.tag, col='high_resolution') <= high_res_small_cutoff:
                continue
            # Check to see if this has been excluded from building
            elif self.datasets.all_masks().get_value(name='no_analyse', id=dataset.tag) == True:
                self.log('Rejecting Dataset {!s}: Excluded from analysis'.format(dataset.tag))
                continue
            elif self.datasets.all_masks().get_value(name='old datasets', id=dataset.tag) and (not self.args.method.reprocess_existing_datasets) and (dataset.tag not in datasets_for_reprocessing):
                self.log('Rejecting Dataset {!s}: Already Processed (Old Dataset)'.format(dataset.tag))
                continue
            else:
                self.datasets.all_masks().set_value(name=analysis_mask_name, id=dataset.tag, value=True)
        return analysis_mask_name, self.datasets.all_masks().get_mask(name=analysis_mask_name)

    def truncate_scaled_data(self, datasets, res_truncate):
        """Truncate data at the same indices across all the datasets"""

        self.log('----------------------------------->>>', True)
        self.log('Truncating Reflection Data', True)

        # Column names
        all_cols = self.args.params.maps.structure_factors
        ref_cols = self.args.input.reference.structure_factors if self.args.input.reference.structure_factors else all_cols

        # Calculate which reflections are present in the reference dataset
        ref_size = self.datasets.reference().data.structure_factors[ref_cols].set().size()
        self.log('Number of Reflections in Reference Dataset: {!s}'.format(ref_size))

        # TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
        #
        # NEED TO REMOVE TRUNCATION STEP IF SPACEGROUPS ARE NOT IDENTICAL (OR UNIT CELLS ISOMORPHOUS)
        # 1) ONLY TRUNCATE IF IDENTICAL SPACEGROUP
        # 2) IF NOT IDENTICAL SPACEGROUP ASSERT THAT EACH DATASET HAS A FULL (COMPLETE) SET OF MILLER OF INDICES
        #
        # TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO

        # Truncate reflections to the common set (not including the reference dataset)
        common_set = datasets[0].data.structure_factors[all_cols].set()
        for dataset in datasets[1:]:
            common_set = common_set.common_set(dataset.data.structure_factors[all_cols], assert_is_similar_symmetry=False)

        self.log('----------------------------------->>>', True)
        self.log('Number of Common Reflections between Datasets: {!s} ({!s}% of reference)'.format(common_set.size(), int(100.0*common_set.size()/ref_size)))
        self.log('After Truncation - Reflections per dataset: {!s}'.format(common_set.size()))

        # Create maps for all of the datasets (including the reference dataset)
        self.datasets.reference().data.structure_factors['truncated'] = self.datasets.reference().data.structure_factors[ref_cols].common_set(common_set, assert_is_similar_symmetry=False)
        for dataset in datasets:
            dataset.data.structure_factors['truncated'] = dataset.data.structure_factors[all_cols].common_set(common_set, assert_is_similar_symmetry=False)

        # Get resolution range of the truncated data and plot histogram
        reslns = [d.data.structure_factors['truncated'].d_min() for d in datasets]
        min_res, max_res = min(reslns), max(reslns)
        self.log('After Truncation - Resolution Range: {!s}-{!s}'.format(min_res, max_res))
        if self.settings.plot_graphs:
            f_name = os.path.join(self.file_manager.get_dir('m_graphs'), '{!s}A-truncated_dataset_resolutions.png'.format(res_truncate))
            self.log('Writing histogram to {}'.format(f_name))
            simple_histogram(filename=f_name, data=reslns, title='Truncated dataset resolutions', x_lab='Resolution (A)', n_bins=15)

    def load_reference_map(self, map_resolution=0):
        """Load the reference map, and calculate some map statistics"""

        # Get the reference dataset
        ref_dataset = self.datasets.reference()

        # Take the scaled diffraction data for the dataset and create fft
        fft_map = ref_dataset.data.get_fft_map(structure_factors='truncated',
                    resolution_factor=self.params.maps.resolution_factor, d_min=map_resolution)

        # Scale the map
        if   self.params.maps.scaling == 'none':   pass
        elif self.params.maps.scaling == 'sigma':  fft_map.apply_sigma_scaling()
        elif self.params.maps.scaling == 'volume': fft_map.apply_volume_scaling()

        # Extract the points for the map (in the grid frame)
        masked_cart = self.grid.grid2cart(self.grid.global_mask().outer_mask(), origin=True)
        # Transform to the frame of the reference dataset (this should become unnecessary in the future)
        masked_cart = self.datasets.reference().grid2nat(masked_cart)
        # Create map handler in the native frame and extract the map values
        ref_map_true = ref_dataset.data.get_electron_density_map(fft_map_name='truncated').as_map()
        masked_vals = ref_map_true.get_cart_values(masked_cart)

        # Create a new electron density map object for the "grid map"
        ref_map = ElectronDensityMap(map_data=masked_vals, unit_cell=self.grid.unit_cell(),
                        map_indices=self.grid.global_mask().outer_mask_indices(),
                        map_size=self.grid.grid_size(), sparse=True)

        # Add some meta for debugging, etc
        ref_map.meta.type = 'reference-map'
        ref_map.meta.resolution = map_resolution
        ref_map.meta.map_mean = ref_map.data.min_max_mean().mean
        ref_map.meta.map_rms = ref_map.data.standard_deviation_of_the_sample()

        return ref_map

    def load_and_morph_maps(self, datasets, ref_map, map_resolution=0):
        """Create map from miller arrays. Transform map into the reference frame by sampling at the given points."""

        assert ref_map.is_sparse(), 'Reference map is not in sparse form'

        self.log('----------------------------------->>>', True)
        self.log('Loading electron density maps @ {!s}A'.format(map_resolution), True)

        # Create holder for the output map objects
        map_list = MapHolderList()

        self.log('----------------------------------->>>', True)
        self.log('Converting structure factors to electron density maps')
        start = time.time()
        for dataset in datasets:
            dataset.data.get_fft_map(structure_factors='truncated', resolution_factor=self.params.maps.resolution_factor, d_min=map_resolution)
        finish = time.time()
        self.log('> FFT-ing structure factors ({!s} Datasets) > Time Taken: {!s} seconds'.format(len(datasets), int(finish-start)), True)

        self.log('----------------------------------->>>', True)
        print('Loading maps (using {!s} cores)'.format(self.settings.cpus))
        start = time.time()
        arg_list = [MapLoader(dataset=d, grid=self.grid, reference_map=ref_map, args=self.args, verbose=self.settings.verbose) for d in datasets]
        # Print a sort of progress bar
        print('1'+''.join(['{:<5}'.format(i) for i in range(0,len(arg_list)+5,5)])[2:])
        print(' '*len(arg_list)+'|\r', end=''); sys.stdout.flush()
        dataset_maps = easy_mp.pool_map(func=wrapper_run, args=arg_list, processes=self.settings.cpus, chunksize=1)
        print('|')
        map_list.add(dataset_maps)

        # Assign parents and clear maps to save memory
        for map in map_list.all():
            map.parent = self.datasets.get(tag=map.meta.tag)
            map.parent.data.fft_map['truncated'] = None

        finish = time.time()
        self.log('> Loading maps ({!s} Datasets) > Time Taken: {!s} seconds'.format(map_list.size(), int(finish-start)), True)
        self.log('----------------------------------->>>')

        # Set the mask lengths and entry ids
        map_list.all_masks().set_index_ids(ids=map_list.all_tags())

        return map_list

    def collate_event_counts(self):
        """Collate events from all of the datasets"""

        self.log('----------------------------------->>>', False)
        self.log('Collating Clusters', False)

        # List of points to be returned
        all_dataset_events = dict([(d.tag, d.events) for d in self.datasets.all()])

        # Print Cluster Summaries
        event_num = [(k, len(all_dataset_events[k])) for k in sorted(all_dataset_events.keys()) if all_dataset_events[k]]
        event_total = sum([a[1] for a in event_num])

        return event_total, event_num, all_dataset_events

    def cluster_events_and_update(self, events=[], update_tables=True, update_output=True):
        """Cluster events to sites and add information to the pandda tables"""

        if not events:
            print('No Events Found')
            return None

        self.log('----------------------------------->>>', True)
        self.log('Clustering identified events: {} Event(s)'.format(len(events)))
        site_list = cluster_events(events=events, cutoff=15.0/self.grid.grid_spacing(), linkage='average')
        site_list.sort(key=lambda s: (s.info.num_events, max([e.cluster.max for e in s.children])), reverse=True).renumber()
        # Add meta to the site list TODO implement this function -- blank at the moment TODO
        [s.find_protein_context(hierarchy=self.datasets.reference().model.hierarchy) for s in site_list.children]
        # Update the pandda tables?
        if update_tables:
            self.update_site_table(site_list=site_list, clear_table=True)
            self.update_event_table_site_info(events=events)
        # Generate output images and graphs?
        if update_output:
            # Plot output graph of site list
            self.log('Deleting old images: ')
            delete_with_glob(glob_str=self.file_manager.get_file('analyse_site_graph_mult').format('*'))
            bar.multiple_bar_plot_over_several_images(
                                    f_template = self.file_manager.get_file('analyse_site_graph_mult'),
                                    plot_vals  = [sorted([e.cluster.max for e in s.children],reverse=True) for s in site_list.children]   )
            # Create pictures of the sites on the protein
            self.make_pymol_site_image_and_scripts(site_list=site_list, make_images=True)

        return site_list

#    def image_blob(self, script, image, dataset, point, point_no, towards=[10,10,10]):
#        """Take pictures of the maps with ccp4mg"""
#
#        from giant.graphics import calculate_view_quaternion, multiply_quaternions
#
#        # Get the template to be filled in
#        template = PANDDA_HTML_ENV.get_template('ccp4mg-pic.py')
#
#        orientation = calculate_view_quaternion(towards, point)
#        rotate_1 = multiply_quaternions(orientation, (0.0, 0.5**0.5, 0.0, 0.5**0.5))
#        rotate_2 = multiply_quaternions(orientation, (0.5**0.5, 0.0, 0.0, 0.5**0.5))
#
#        for view_no, view in enumerate([orientation, rotate_1, rotate_2]):
#
#            view_script = script.format(point_no, view_no)
#            view_image  = image.format(point_no, view_no)
#
#            ccp4mg_script = template.render({
#                                                'view'  :{
#                                                                'camera_centre' : [-1*c for c in point],
#                                                                'orientation'   : list(view)
#                                                            },
#                                                'mol'   :{
#                                                                'path'  : dataset.file_manager.get_file('aligned_structure'),
#                                                                'name'  : 'aligned_structure'
#                                                            },
#                                         #       'map'   :{
#                                         #                       'path'    : dataset.file_manager.get_file('sampled_map'),
#                                         #                       'name'    : 'sampled_map',
#                                         #                       'contour' : [1]
#                                         #                   },
#                                                'diff_map' :{
#                                                                'path'    : dataset.file_manager.get_file('z_map_corrected_normalised'),
#                                                                'name'    : 'diff_map',
#                                                            #    'neg-contour' : -3,
#                                                                'pos-contour' : [2,3,4,5]
#                                                            }
#                                            })
#
#            # Write out the ccp4mg script to the dataset's scripts folder
#            with open(view_script, 'w') as fh:
#                fh.write(ccp4mg_script)
#
#            # Make the images
#            c = CommandManager('ccp4mg')
#            c.SetArguments(['-norestore','-picture', view_script, '-R', view_image, '-RO', """'{"size":"1500x1500"}'""", '-quit'])
#            c.Run()
#
#            if not os.path.exists(view_image):
#                print('FAILED TO MAKE IMAGES')
#                print(c.err)

    def write_map_analyser_maps(self, map_analyser):
        """Write statistical maps for a map_analyser object"""

        map_res = map_analyser.meta.resolution

        self.log('----------------------------------->>>')
        self.log('=> Writing characterised statistical maps @ {!s}A'.format(map_res))

        write_array_as_map(grid     = self.grid,
                           array    = map_analyser.statistical_maps.mean_map.copy().as_dense().data,
                           f_name   = self.file_manager.get_file('mean_map').format(map_res))
        write_array_as_map(grid     = self.grid,
                           array    = map_analyser.statistical_maps.stds_map.copy().as_dense().data,
                           f_name   = self.file_manager.get_file('stds_map').format(map_res))
        write_array_as_map(grid     = self.grid,
                           array    = map_analyser.statistical_maps.sadj_map.copy().as_dense().data,
                           f_name   = self.file_manager.get_file('sadj_map').format(map_res))
        write_array_as_map(grid     = self.grid,
                           array    = map_analyser.statistical_maps.skew_map.copy().as_dense().data,
                           f_name   = self.file_manager.get_file('skew_map').format(map_res))
        write_array_as_map(grid     = self.grid,
                           array    = map_analyser.statistical_maps.kurt_map.copy().as_dense().data,
                           f_name   = self.file_manager.get_file('kurt_map').format(map_res))
        write_array_as_map(grid     = self.grid,
                           array    = map_analyser.statistical_maps.bimo_map.copy().as_dense().data,
                           f_name   = self.file_manager.get_file('bimo_map').format(map_res))

    def write_output_csvs(self):
        """Write CSV file of dataset variables"""

        self.log('----------------------------------->>>')
        self.log('Writing Dataset + Dataset Map Summary CSV')

        # Write the dataset information to csv file
        self.tables.dataset_info.to_csv(path_or_buf=self.file_manager.get_file('dataset_info'), index_label='dtag')
        self.tables.dataset_map_info.to_csv(path_or_buf=self.file_manager.get_file('dataset_map_info'), index_label='dtag')
        self.datasets.all_masks().table.to_csv(path_or_buf=self.file_manager.get_file('dataset_masks'), index_label='id')

        self.log('----------------------------------->>>')
        self.log('Writing COMBINED Dataset Summary CSV')

        # Join the tables on the index of the main table
        comb_tab = (
            self.tables.dataset_info
                .join(self.tables.dataset_map_info, how='outer')
                .join(self.datasets.all_masks().table[PanddaMaskNames.write_mask_names], how='outer')
        )
        comb_tab.to_csv(path_or_buf=self.file_manager.get_file('dataset_combined_info'), index_label='dtag')

        self.log('----------------------------------->>>')
        self.log('Writing Event+Site Summary CSVs')

        # Sort the event data by z-peak and write out
        sort_eve = self.tables.event_info.sort_values(by=['site_idx',self.args.results.events.order_by], ascending=[1,0])
        sort_eve = sort_eve.join(comb_tab, how='right')
        sort_eve.to_csv(path_or_buf=self.file_manager.get_file('event_info'))
        # Sort the sites by number of events and write out
        sort_sit = self.tables.site_info.sort_values(by=[self.args.results.sites.order_by],ascending=[0])
        sort_sit.to_csv( path_or_buf=self.file_manager.get_file('site_info'))

    def update_site_table(self, site_list, clear_table=True):
        """Add site entries to the site table"""

        # Clear an existing table
        if clear_table:
            self.tables.site_info = pandas.DataFrame(data    = None,
                                                     index   = self.tables.site_info.index.reindex([])[0],
                                                     columns = self.tables.site_info.columns)
        # Go through and update the site information
        for site in site_list.children:
            self.tables.site_info.loc[site.id,:] = None
            centroid_cart = tuple(flex.double(site.info.centroid)*self.grid.grid_spacing())
            self.tables.site_info.set_value(site.id, 'centroid', centroid_cart)
            self.tables.site_info.set_value(site.id, 'native_centroid', tuple(self.datasets.reference().model.alignment.ref2nat(coordinates=[centroid_cart])[0]))

    def make_pymol_site_image_and_scripts(self, site_list, make_images=True):
        """Generate pymol script to mark the location of identified sites"""

        pymol_str =  '# Mark the identified sites on the protein\n'
        pymol_str += 'from pymol import cmd\n'
        pymol_str += 'from pymol.cgo import *\n'
        pymol_str += 'cmd.load("{}", "reference")\n'.format(os.path.relpath(self.file_manager.get_file('reference_on_origin'), start=self.file_manager.get_dir('output_summaries')))
        pymol_str += 'cmd.show_as("cartoon", "reference")\n'
        pymol_str += 'cmd.color("cyan", "reference")\n'
        # Add sphere at each of the sites
        for site in site_list.children:
            # Only print the site if it has more than one event
            if len(site.children) > 1:
                lab = 'site_{}'.format(site.id)
                com = tuple(flex.double(site.info.centroid)*self.grid.grid_spacing())
                pymol_str += 'cmd.pseudoatom("{}", pos={}, vdw=2.5)\n'.format(lab, com)
                pymol_str += 'cmd.show("sphere", "{}")\n'.format(lab)
                pymol_str += 'cmd.label("{}", "{}")\n'.format(lab, site.id)
                pymol_str += 'cmd.color("deepteal", "{}")\n'.format(lab)
                pymol_str += 'cmd.set("label_color", "white", "{}")\n'.format(lab)
            # Label events as smaller spheres
            for event in site.children:
                lab = 'event'
                com = tuple(flex.double(event.cluster.centroid)*self.grid.grid_spacing())
                pymol_str += 'cmd.pseudoatom("{}", pos={}, vdw=0.5)\n'.format(lab, com)
                pymol_str += 'cmd.show("sphere", "{}")\n'.format(lab)
                pymol_str += 'cmd.color("blue", "{}")\n'.format(lab)
        # Set label things...
        pymol_str += 'cmd.set("label_size", 25)\n'
        pymol_str += 'cmd.set("label_position", (0,0,4))\n'
        pymol_str += 'cmd.bg_color(color="white")\n'
        # Write as python script
        with open(self.file_manager.get_file(file_tag='pymol_sites_py'), 'w') as fh:
            fh.write(pymol_str)

        # Run Pymol to generate images and output to pngs
        if make_images:
            pymol_str =  '# Load the protein representation and output images of sites\n'
            pymol_str += 'run {}\n'.format(os.path.relpath(self.file_manager.get_file(file_tag='pymol_sites_py'), start=self.file_manager.get_dir('output_summaries')))
            pymol_str += 'set ray_opaque_background, off\n'
            pymol_str += 'set specular, off\n'
            pymol_str += 'orient\n'
            pymol_str += 'png {}, width=1200, dpi=300, ray=1\n'.format(os.path.relpath(self.file_manager.get_file(file_tag='pymol_sites_png_1'), start=self.file_manager.get_dir('output_summaries')))
            pymol_str += 'rotate y, 180\n'
            pymol_str += 'png {}, width=1200, dpi=300, ray=1\n'.format(os.path.relpath(self.file_manager.get_file(file_tag='pymol_sites_png_2'), start=self.file_manager.get_dir('output_summaries')))
            pymol_str += 'quit'

            with open(self.file_manager.get_file(file_tag='pymol_sites_pml'), 'w') as fh:
                fh.write(pymol_str)

            # Change into directory as script runs off of relative paths
            os.chdir(self.file_manager.get_dir('output_summaries'))
            c = CommandManager('pymol')
            c.add_command_line_arguments(['-k', '-q', '-c', self.file_manager.get_file(file_tag='pymol_sites_pml')])
            try:    c.run()
            except: print("Failed to start pymol - maybe it's not available?")
            # Change back to top directory
            os.chdir(self.out_dir)

        #os.remove(self.file_manager.get_file(file_tag='pymol_sites_pml'))
        #os.remove(self.file_manager.get_file(file_tag='pymol_sites_py'))

    def add_event_to_event_table(self, dataset, event):
        """Add event entries to the event table"""

        # Check event has not been added previously
        assert event.id not in self.tables.event_info.index.values.tolist(), 'Event Already Added!: {!s}'.format(event.id)
        # Add values to a new row in the table
        self.tables.event_info.loc[event.id,:] = None
        # Default to site_idx of 0 if no site given
        if event.parent:    site_idx = event.parent.id
        else:               site_idx = 0
        self.tables.event_info.set_value(event.id, 'site_idx', site_idx)
        # Event and cluster information
        self.tables.event_info.set_value(event.id, '1-BDC',  round(1.0-event.info.estimated_bdc,2))
        self.tables.event_info.set_value(event.id, 'z_peak', round(event.cluster.max,2))
        self.tables.event_info.set_value(event.id, 'z_mean', round(event.cluster.mean,2))
        self.tables.event_info.set_value(event.id, 'cluster_size', event.cluster.size)
        self.tables.event_info.set_value(event.id, ['refx','refy','refz'], list(self.grid.grid2cart([event.cluster.peak],origin=False)[0]))
        self.tables.event_info.set_value(event.id, ['x','y','z'], list(dataset.model.alignment.ref2nat(coordinates=self.grid.grid2cart([event.cluster.peak],origin=False))[0]))

    def update_event_table_site_info(self, events):
        """Update the event table for pre-existing events"""
        for e in events:
            assert e.id, 'NO ID GIVEN: {!s}'.format(e.id)
            assert e.parent, 'EVENT HAS NO PARENT: {!s}'.format(e.parent)
            self.tables.event_info.set_value(e.id, 'site_idx', e.parent.id)


class PanddaMapAnalyser(object):


    def __init__(self, dataset_maps, meta=None, statistical_maps=None, parent=None, log=None):
        """Class to hold dataset maps, statistical maps and meta data for a set of related maps. Also holds functions for analysing the maps."""

        # Validate the dataset maps
        assert (dataset_maps is None) or isinstance(dataset_maps, MapHolderList), 'dataset_maps must be stored in a MapHolderList. Type given: {!s}'.format(type(dataset_maps))
        # Validate the meta
        if meta:
            assert isinstance(meta, Meta), 'meta must be of type Meta. Type given: {!s}'.format(type(meta))
        else:
            meta = Meta()
        # Validate the statistical maps
        if statistical_maps:
            assert isinstance(statistical_maps, PanddaStatMapList), 'statistical_maps must be of type MapList. Type given: {!s}'.format(type(statistical_maps))
        else:
            statistical_maps = PanddaStatMapList()
        # Validate the parent object (main pandda object)
        if parent:
            assert isinstance(parent, PanddaMultiDatasetAnalyser), 'parent must be of type PanddaMultiDatasetAnalyser. Type given: {!s}'.format(type(parent))
        # Create log as appropriate
        if log:
            pass
        elif parent:
            log = parent.log
        else:
            log = Log(verbose=False)

        self.dataset_maps = dataset_maps
        self.meta = meta
        self.statistical_maps = statistical_maps
        self.parent = parent
        self.log = log

        if self.dataset_maps is not None: self._validate()

    def _validate(self):
        """Check that all of the added maps are the same size etc..."""

        m_ref = None
        for m in self.dataset_maps.all():
            if m_ref is None: m_ref = m

            assert m.data.all()  == m_ref.data.all(), 'Not all maps are the same shape'
            assert m.data.size() == m_ref.data.size(), 'Not all maps are the same size'

            if not hasattr(m.meta, 'map_uncertainty'):
                m.meta.map_uncertainty = None

        self.meta.map_data_size  = m_ref.data.size()
        self.meta.map_data_shape = m_ref.data.all()

    def calculate_mean_map(self, mask_name=None):
        """Calculate the mean map from all of the different observations"""

        self.log('----------------------------------->>>', True)
        self.log('Calculating Mean Map', True)

        self.log('Mean Map Size: {}'.format(self.meta.map_data_size))

        # Chunk the points into groups - Compromise between cpu time and memory usage - ~200 dataset -> chunksize of 5000
        chunk_size = 500*iceil(1000.0/self.dataset_maps.size(mask_name=mask_name))
        chunk_idxs = [i for i in range(0, self.meta.map_data_size, chunk_size)]
        num_chunks = len(chunk_idxs)

        self.log('Iterating through {!s} points in {!s} chunks'.format(self.meta.map_data_size, num_chunks), True)
        t1 = time.time()

        mean_map_vals = numpy.zeros(self.meta.map_data_size)
        medn_map_vals = numpy.zeros(self.meta.map_data_size)

        for i_chunk, chunk_start in enumerate(chunk_idxs):
            status_bar_2(n=i_chunk, n_max=num_chunks)

            tmp_map_vals = numpy.array([m.data[chunk_start:chunk_start+chunk_size] for m in self.dataset_maps.mask(mask_name=mask_name)])

            # Check that the output values are the expected dimensions
            if i_chunk+1 < num_chunks:
                assert len(tmp_map_vals) == self.dataset_maps.size(mask_name=mask_name)
                assert len(tmp_map_vals.T) == chunk_size

            tmp_map_means = numpy.mean(tmp_map_vals, axis=0)
            mean_map_vals[chunk_start:chunk_start+chunk_size] = tmp_map_means
            tmp_map_medns = numpy.median(tmp_map_vals, axis=0)
            medn_map_vals[chunk_start:chunk_start+chunk_size] = tmp_map_medns

        status_bar_2(n=num_chunks, n_max=num_chunks)

        t2 = time.time()
        self.log('> Calculation of mean map > Time Taken: {!s} seconds'.format(int(t2-t1)))

        self.statistical_maps.mean_map = m.new_from_template(map_data=flex.double(mean_map_vals.flatten()), sparse=m.is_sparse())
        self.statistical_maps.medn_map = m.new_from_template(map_data=flex.double(medn_map_vals.flatten()), sparse=m.is_sparse())

        return self.statistical_maps.mean_map

    def calculate_map_uncertainties(self, masked_idxs=None, mask_name=None, q_cut=1.5, cpus=1):
        """Calculate the uncertainty in each of the different maps"""

        if masked_idxs is None:
            masked_idxs = flex.size_t(range(0, self.meta.map_data_size))
        else:
            assert max(masked_idxs) < self.meta.map_data_size, 'masked_idxs out of range of map'
            masked_idxs = flex.size_t(masked_idxs)

        # Extract masked map values from the mean map... and sort them
        mean_vals = self.statistical_maps.mean_map.data.select(masked_idxs)

        self.log('----------------------------------->>>')
        self.log('Selecting datasets for uncertainty calculation')

        arg_list = []

        for i_m, m in enumerate(self.dataset_maps.mask(mask_name=mask_name)):

            if m.meta.map_uncertainty is not None:
                print('SKIPPING Dataset {!s} ({!s}/{!s})'.format(m.meta.tag, i_m+1, self.dataset_maps.size(mask_name=mask_name)))
                arg_list.append(None)
                continue

            if m.parent:
                file_manager = m.parent.file_manager
            else:
                file_manager = None

            u = UncertaintyCalculator(query_values=m.data.select(masked_idxs), ref_values=mean_vals, file_manager=file_manager)
            arg_list.append(u)

        self.log('Calculating uncertainties of {!s} maps (using {!s} cores)'.format(len(arg_list), cpus))
        self.log('----------------------------------->>>')

        t1 = time.time()
        print('1'+''.join(['{:<5}'.format(i) for i in range(0,len(arg_list)+5,5)])[2:])
        print(' '*len(arg_list)+'|\r', end=''); sys.stdout.flush()
        map_uncertainties = easy_mp.pool_map(func=wrapper_run, args=arg_list, processes=cpus, chunksize=1)
        print('|')
        t2 = time.time()
        self.log('> Calculation of map uncertainties > Time Taken: {!s} seconds'.format(int(t2-t1)))
        self.log('----------------------------------->>>')

        for i_m, m in enumerate(self.dataset_maps.mask(mask_name=mask_name)):
            map_unc = map_uncertainties[i_m]
            if m.meta.map_uncertainty is not None:
                assert map_unc is None
            else:
                assert map_unc is not None
                m.meta.map_uncertainty = map_unc
                self.log('Dataset {}: Uncertainty {:.4f}'.format(m.meta.tag, m.meta.map_uncertainty))

        self.log('----------------------------------->>>')

        return [m.meta.map_uncertainty for m in self.dataset_maps.mask(mask_name=mask_name)]

    def calculate_statistical_maps(self, mask_name=None, ignore_warnings=True, cpus=1):
        """Take the sampled maps and calculate statistics for each grid point across the datasets"""

        # Create statistics objects for each grid point
        self.log('----------------------------------->>>')
        self.log('Calculating statistics of {} grid points (using {!s} cores)'.format(self.meta.map_data_size, cpus))

        if ignore_warnings:
            self.log('Suppressing Numpy Iteration Warnings... (Nothing to worry about)')
            warnings.simplefilter('ignore', category=RuntimeWarning)

        # Extract the map uncertainties
        uncertainties = [m.meta.map_uncertainty for m in self.dataset_maps.mask(mask_name=mask_name)]
        assert uncertainties.count(None) == 0, 'some maps have not got associated uncertainties'

        # Chunk the points into groups - Compromise between cpu time and memory usage - 1000 per cpu at 50 datasets
        chunk_size = iceil(1000.0*cpus*50.0/self.dataset_maps.size(mask_name=mask_name))
        chunk_idxs = [i for i in range(0, self.meta.map_data_size, chunk_size)]
        num_chunks = len(chunk_idxs)

        # Second level of iteration - split the first chunk level between the cpus
        chunk_size_2 = iceil(1.0*chunk_size/cpus)
        chunk_idxs_2 = [i for i in range(0, chunk_size, chunk_size_2)]
        num_chunks_2 = len(chunk_idxs_2)

        self.log('Iterating through {!s} points in blocks of {!s} ({!s} chunks)'.format(self.meta.map_data_size, chunk_size, num_chunks))

        t1 = time.time()

        # Output array of the 5 statistics for each map point
        point_statistics = numpy.zeros((self.meta.map_data_size, 5))

        tot = 0
        for i_chunk, chunk_start in enumerate(chunk_idxs):
            status_bar_2(n=i_chunk, n_max=num_chunks)

            # Argument list for multiprocessing
            arg_list = []

            # Loop through the secondary chunks and send for multi-core processing
            for i_chunk_2, chunk_start_2 in enumerate(chunk_idxs_2):

#                print('\n--------------------->')
#                print('Chunk {}-{}'.format(i_chunk, i_chunk_2))
#                print('Getting {} - {}'.format(chunk_start+chunk_start_2, chunk_start+chunk_start_2+chunk_size_2))

                # Lower limit - always the beginning of the chunk
                l1 = chunk_start+chunk_start_2
                # Upper limit - full chunk size, limited by the larger chunk size, or by map size
                l2 = min(chunk_start+chunk_start_2+chunk_size_2, chunk_start+chunk_size, self.meta.map_data_size)

                if l1 >= l2:
                    continue

                # Extract map values from the maps
                map_vals = [m.data[l1:l2] for m in self.dataset_maps.mask(mask_name=mask_name)]
                # Want to iterate over grid points not datasets
                map_vals = numpy.transpose(map_vals)
                assert map_vals.shape[1] == self.dataset_maps.size(mask_name=mask_name)

                # Create DensityStatistics object for analysis of the density variation
                arg_list.append(DensityStatistics(observations_array=map_vals, uncertainties=uncertainties))

            if not arg_list: continue

            # Calculate the statistis of the grid points
            tmp_point_statistics = easy_mp.pool_map(func=wrapper_run, args=arg_list, processes=cpus)

            # Put values into the output array
            offset = 0
            for point_vals in tmp_point_statistics:
                assert point_vals.shape[1] == 5
                l1 = chunk_start+offset
                l2 = l1+point_vals.shape[0]
                if not (point_statistics[l1:l2,:] == 0.0).all():
                    print('Overwriting data?!')
                    print(point_statistics[l1-10:l2+10,:])
                    assert point_statistics[l1:l2,:].shape == point_vals.shape, '{} != {}'.format(point_statistics[l1:l2,:].shape, point_vals.shape)
                point_statistics[l1:l2,:] = point_vals
                offset += point_vals.shape[0]
            tot += offset

        status_bar_2(n=num_chunks, n_max=num_chunks)

        # Check that we've calculated the right number of things
        assert tot == self.meta.map_data_size, 'tot {}, map size {}'.format(tot, self.meta.map_data_size)

        t2 = time.time()
        self.log('> Calculation of map variation statistics > Time Taken: {!s} seconds'.format(int(t2-t1)))

        # Use the mean map as the template for the other maps
        mean_map = self.statistical_maps.mean_map
        # Create the other statistical maps
        self.statistical_maps.stds_map = mean_map.new_from_template(map_data=flex.double(point_statistics[:,0].tolist()), sparse=mean_map.is_sparse())
        self.statistical_maps.sadj_map = mean_map.new_from_template(map_data=flex.double(point_statistics[:,1].tolist()), sparse=mean_map.is_sparse())
        self.statistical_maps.skew_map = mean_map.new_from_template(map_data=flex.double(point_statistics[:,2].tolist()), sparse=mean_map.is_sparse())
        self.statistical_maps.kurt_map = mean_map.new_from_template(map_data=flex.double(point_statistics[:,3].tolist()), sparse=mean_map.is_sparse())
        self.statistical_maps.bimo_map = mean_map.new_from_template(map_data=flex.double(point_statistics[:,4].tolist()), sparse=mean_map.is_sparse())

        return self.statistical_maps

    def calculate_z_map(self, method, map=None, tag=None, uncertainty=None):
        """Calculate the z-map relative to the mean and std map"""

        assert method in ['naive','adjusted','uncertainty','adjusted+uncertainty']
        assert [tag, map].count(None) == 1, 'Must provide tag OR map'

        if tag:
            map = self.dataset_maps.get(tag=tag)

        if uncertainty is not None:
            uncertainty = map.meta.map_uncertainty

        if 'uncertainty' in method:
            assert uncertainty is not None

        # Calculate Z-values
        if method == 'naive':
            z_map = (map - self.statistical_maps.mean_map.data)*(1.0/self.statistical_maps.stds_map.data)
        elif method == 'adjusted':
            z_map = (map - self.statistical_maps.mean_map.data)*(1.0/self.statistical_maps.sadj_map.data)
        elif method == 'uncertainty':
            z_map = (map - self.statistical_maps.mean_map.data)*(1.0/uncertainty)
        elif method == 'adjusted+uncertainty':
            z_map = (map - self.statistical_maps.mean_map.data)*(1.0/flex.sqrt(self.statistical_maps.sadj_map.data**2 + uncertainty**2))
        else:
            raise Exception('method not found: {!s}'.format(method))

        return z_map


import os

from libtbx import easy_pickle

from bamboo.common.logs import Log

class Program(object):
    """Class meant to provide basic functionality for programs and pipelines"""

    _NAME = None
    _TEXT = None
    _VERSION = None

    _allowed_statuses = ['running','done','errored']

    log = Log()

    def write_running_parameters_to_log(self):
        self.log('===================================>>>', True)
        self.log('PROCESSED ARGUMENTS', True)
        self.log('===================================>>>', True)
        self.log(self.master_phil.format(python_object=self._input_params).as_str(), True)

        self.log('===================================>>>', True)
        self.log('ARGUMENTS DIFFERENT TO DEFAULTS', True)
        self.log('===================================>>>', True)
        self.log(self.master_phil.fetch_diff(source=self.master_phil.format(python_object=self._input_params)).as_str(), True)

    def check_for_matplotlib(self):
        """Check to see whether we can load matplotlib"""
        try:
            import matplotlib
            matplotlib.use('Agg')
            # Setup so that we can write without a display connected
            matplotlib.interactive(0)
            # Find the backend
            default_backend, validate_function = matplotlib.defaultParams['backend']
            self.log('===================================>>>')
            self.log('MATPLOTLIB LOADED. Using Backend: {!s}'.format(default_backend))
            from matplotlib import pyplot
            # Test non-interactive-ness
            assert not pyplot.isinteractive(), 'PYPLOT IS STILL IN INTERACTIVE MODE'
            # Set the style to ggplot (the prettiest)
            pyplot.style.use('ggplot')
            self.log('PYPLOT loaded successfully')
            self.log('===================================>>>', True)
            self.log('', True)
            return True
        except:
            self.log('===================================>>>', True)
            self.log('>> COULD NOT IMPORT MATPLOTLIB. CANNOT GENERATE GRAPHS.', True)
            self.log('===================================>>>', True)
            self.log('', True)
            return False

    def update_status(self, status):
        """Set log files to indicate the status of the program"""

        assert status in self._allowed_statuses
        # Delete any that may exist
        existing_files = [self.output_handler.get_file('status').format(f) for f in self._allowed_statuses]
        [os.remove(f) for f in existing_files if os.path.exists(f)]
        # Create the new  status file
        with open(self.output_handler.get_file('status').format(status), 'w') as fh: fh.write('')

    def pickle(self, pickle_file, pickle_object, overwrite=True):
        """Takes an object and pickles it"""
        if os.path.exists(pickle_file) and not overwrite:
            self.log('NOT PICKLING: {!s}'.format(os.path.relpath(pickle_file, start=self.out_dir)))
        else:
            self.log('Pickling Object: {!s}'.format(os.path.relpath(pickle_file, start=self.out_dir)))
            easy_pickle.dump(pickle_file, pickle_object)

    def unpickle(self, pickle_file):
        """Takes an object and unpickles it"""
        self.log('Unpickling File: {!s}'.format(os.path.relpath(pickle_file, start=self.out_dir)))
        return easy_pickle.load(pickle_file)


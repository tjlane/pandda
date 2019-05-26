import os

from libtbx import adopt_init_args
from bamboo.common.path import easy_directory

class PandemicAdpFileSystem:


    def __init__(self,
            output_directory,
            ):
        output_directory    = self.make(output_directory)
        hierarchy_directory = self.make(output_directory, 'hierarchical_model')
        optimisation_directory = self.make(hierarchy_directory, 'optimisation')
        structure_directory = self.make(output_directory, 'output_structures')
        analysis_directory  = self.make(output_directory, 'output_analysis')
        adopt_init_args(self, locals())

    def make(self, *args):
        return easy_directory(os.path.join(*args))

class TidyOutputFolder:


    def __init__(self,
            delete_mtzs = False,
            compress_logs = False,
            compress_pdbs = False,
            log = None,
            ):
        adopt_init_args(self, locals())

    def __call__(self,
            output_directory,
            ):

        """Remove MTZ files, etc"""
        if not actions: return
        self.log.heading('Tidying/Cleaning output folder')
        if 'delete_mtzs' in actions:
            self.log.subheading('Deleting unnecessary MTZ files')
            for i_mod, m in enumerate(self.models):
                for f in [m.o_mtz, m.r_mtz]:
                    if (f is not None) and os.path.exists(f):
                        self.log('Removing file: {}'.format(f))
                        os.remove(f)
        if 'compress_logs' in actions:
            self.log.subheading('Compressing log files')
            from bamboo.common.file import compress_file
            main_log = os.path.abspath(self.log.log_file().path)
            for log in glob.glob(os.path.join(self.file_manager.get_dir('logs'),'*.log')):
                # Do not compress main log file
                if os.path.abspath(log) == main_log:
                    continue
                if self.params.settings.verbose:
                    self.log('Compressing: {}'.format(log))
                zip_file = compress_file(log, delete_original=True)

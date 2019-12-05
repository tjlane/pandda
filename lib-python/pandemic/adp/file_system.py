import os

from libtbx import adopt_init_args
from bamboo.common.path import easy_directory
from bamboo.common.file import compress_file


class PandemicAdpFileSystem:


    def __init__(self, output_directory):

        output_directory    = self.make(output_directory)
        hierarchy_directory = self.make(output_directory, 'hierarchical_model')
        optimisation_directory = self.make(output_directory, 'optimisation')
        structure_directory = self.make(output_directory, 'structures')
        analysis_directory  = self.make(output_directory, 'analysis')
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

    def __call__(self, output_directory):
        """Remove MTZ files, etc"""

        self.log.heading('Tidying/Cleaning output folder')

        if (self.delete_mtzs is True):

            self.log.subheading('Deleting unnecessary MTZ files')

            for r, d, f in os.walk(output_directory):
                for file in f:
                    if file.endswith('.mtz'):
                        filepath = os.path.join(r, file)
                        self.log('Removing {}'.format(filepath))
                        os.remove(filepath)

        if (self.compress_logs is True):

            self.log.subheading('Compressing log files')

            for r, d, f in os.walk(output_directory):
                for file in f:
                    if file.endswith('.log'):
                        filepath = os.path.join(r, file)
                        self.log('Compressing: {}'.format(filepath))
                        zip_file = compress_file(filepath, delete_original=True)

        if (self.compress_pdbs is True): 

            self.log.subheading('Compressing pdb files')

            for r, d, f in os.walk(output_directory):
                for file in f:
                    if file.endswith('.pdb'):
                        filepath = os.path.join(r, file)
                        self.log('Compressing: {}'.format(filepath))
                        zip_file = compress_file(filepath, delete_original=True)
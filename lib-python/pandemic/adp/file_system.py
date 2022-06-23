import giant.logs as lg
logger = lg.getLogger(__name__)

import os

from libtbx import adopt_init_args
from giant.paths import easy_directory, compress_file


class PandemicAdpFileSystem(object):


    def __init__(self, output_directory):

        output_directory    = self.make(output_directory)
        partition_directory = self.make(output_directory, 'level_partitions')
        hierarchy_directory = self.make(output_directory, 'hierarchical_model')
        optimisation_directory = self.make(output_directory, 'optimisation')
        structure_directory = self.make(output_directory, 'structures')
        analysis_directory  = self.make(output_directory, 'analysis')
        adopt_init_args(self, locals())

    def make(self, *args):
        return easy_directory(os.path.join(*args))


class TidyOutputFolder(object):


    def __init__(self,
            delete_mtzs = False,
            compress_logs = False,
            compress_pdbs = False,
            ):
        adopt_init_args(self, locals())

    def __call__(self, output_directory):
        """Remove MTZ files, etc"""

        logger.heading('Tidying/Cleaning output folder')

        if (self.delete_mtzs is True):

            logger.subheading('Deleting unnecessary MTZ files')

            for r, d, f in os.walk(output_directory):
                for file in f:
                    if file.endswith('.mtz'):
                        filepath = os.path.join(r, file)
                        logger('Removing {}'.format(filepath))
                        os.remove(filepath)

        if (self.compress_logs is True):

            logger.subheading('Compressing log files')

            for r, d, f in os.walk(output_directory):
                for file in f:
                    if file.endswith('.log'):
                        filepath = os.path.join(r, file)
                        logger('Compressing: {}'.format(filepath))
                        try:
                            zip_file = compress_file(filepath, delete_original=True)
                        except Exception as e:
                            logger('...failed to compress.')

        if (self.compress_pdbs is True):

            logger.subheading('Compressing pdb files')

            for r, d, f in os.walk(output_directory):
                for file in f:
                    if file.endswith('.pdb'):
                        filepath = os.path.join(r, file)
                        logger('Compressing: {}'.format(filepath))
                        try:
                            zip_file = compress_file(filepath, delete_original=True)
                        except Exception as e:
                            logger('...failed to compress.')

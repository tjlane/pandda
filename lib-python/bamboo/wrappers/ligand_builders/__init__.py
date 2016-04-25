import os

from bamboo.common.command import CommandManager
from bamboo.constants import DEFAULT_LIGAND_NAMES
from bamboo.wrappers.program_selector import allowed_builder_args

class BuilderObject(object):
    """Template Object for creating an Object to generate ligand restraints"""

    def __init__(self, time=True, verbose=True):

        self.allowedArgs = allowed_builder_args
        self.ligname = DEFAULT_LIGAND_NAMES[0]
        # Settings
        self.time = time
        self.verbose = verbose
        self.runtime = -1.0
        self.timeout = 1800
        # Custom Init
        self._custom_init()

    def generate_ligand(self, ligsmile, outdir, outfile, flags=[]):
        """Generate ligand PDB and CIF from smile string"""

        # Prepare standard settings
        self._prepare_to_generate_ligand(ligsmile, outdir, outfile, flags)
        # Check to see if already run, and if not, run the custom part of the object - different for each program
        if os.path.exists(self.outpdb) and os.path.exists(self.outcif) and os.path.exists(self.outlog):
            if self.verbose:
                print('\tLigand already generated - DOING NOTHING')
        else:
            cmd_line_args, std_inpt_args = self._create_program_arguments(ligsmile, outdir, outfile, flags)
            # Initialise CommandManager
            self.Builder = CommandManager(self.program)
            # Set Command-Line Args
            if cmd_line_args: self.Builder.SetArguments(cmd_line_args)
            # Set Standard Input
            if std_inpt_args: self.Builder.SetInput(std_inpt_args)
            # Set Parameters
            self.Builder.SetParameters(timeout=self.timeout)
            # RUN
            self.Builder.Run()
            # Get runtime
            self.runtime = self.Builder.runtime

            try:
                # Could add postprocess here
                pass

            finally:
                self.write_log_file()

        return self.outpdb, self.outcif, self.outlog

    def _prepare_to_generate_ligand(self, ligsmile, outdir, outfile, flags):
        """Set up the generic file names"""

        # Process outputfile
        if '/' in outfile:
            raise ValueError('outfile must be a file, not a path')
        # Record Template and Outdir
        self.outtemplate = os.path.join(outdir,outfile)
        self.outdir      = outdir
        # Record Filenames
        self.outpdb = self.outtemplate+'.pdb'
        self.outcif = self.outtemplate+'.cif'
        self.outlog = self.outtemplate+'.log'

    def write_log_file(self):
        """Write the log file"""

        with open(self.outlog,'w') as logfile:
            # Write out the input command
            logfile.write('\nCOMMAND\n\n')
            logfile.write('\n'.join(self.Builder.command))
            logfile.write('\nINPUT\n\n')
            logfile.write(self.Builder.inp)
            # Write out & err
            logfile.write('STDOUT\n\n')
            logfile.write(self.Builder.out)
            logfile.write('STDERR\n\n')
            logfile.write(self.Builder.err)


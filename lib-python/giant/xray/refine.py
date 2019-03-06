import os,sys,re,glob,copy,shutil,tempfile

from bamboo.common.command import CommandManager
from bamboo.common.path import easy_directory
from bamboo.common.logs import Log

from giant.xray.tls import phenix_find_tls_groups

class _refiner(object):

    program = None

    _kw_list = []

    def __init__(self, pdb_file, mtz_file=None, cif_files=None, out_prefix=None, log=None, **kw_args):

        # Set defaults if not given
        if mtz_file is None:
            mtz_file = pdb_file.replace('.pdb','.mtz')
        if out_prefix is None:
            out_prefix = os.path.splitext(pdb_file)[0] + '-refined'

        if log is None: log=Log()
        self.log = log

        # Main files
        self.pdb_file = pdb_file
        self.mtz_file = mtz_file
        self.cif_files = cif_files

        # Eventual output prefix
        self.out_prefix = out_prefix
        self.out_pdb_file = self.out_prefix+'.pdb'
        self.out_mtz_file = self.out_prefix+'.mtz'
        self.out_log_file = self.out_prefix+'.log'

        # Validate input
        assert os.path.exists(self.pdb_file)
        assert os.path.exists(self.mtz_file)
        assert not os.path.exists(self.out_pdb_file), 'Output file already exists: {}'.format(self.out_pdb_file)
        assert not os.path.exists(self.out_mtz_file), 'Output file already exists: {}'.format(self.out_mtz_file)

        # Create temporary folder for refinement
        self.tmp_dir = tempfile.mkdtemp(prefix='refine-model-')
        self.tmp_pre = os.path.join(self.tmp_dir, 'refine')

        # Command object for refinement
        self.cmd = CommandManager(self.program, log=self.log)
        self.kw_args = kw_args

        self.validate()

    def validate(self):
        """Validate kw_args"""
        if set(self.kw_args.keys()).difference(self._kw_list):
            raise Sorry('Invalid/unsupported arguments provided: {}'.format(set(self.kw_args.keys()).difference(self._kw_list)))

    def print_settings(self):
        self.log.bar()
        self.log("Refinement Command:")
        self.log.bar()
        self.cmd.print_settings()
        self.log.bar()
        self.log('Writing output to:\n\t{}\n\t{}'.format(self.out_pdb_file,self.out_mtz_file))
        self.log.bar()

    def run(self):
        """...run refinement amd export files"""
        ret = self.refine()
        self.export()
        return ret

    def setup(self):
        raise Exception('Dummy class -- not implemented')

    def refine(self):
        """...run the refinement"""
        return self.cmd.run()

    def export(self):
        """Copy files to output destination"""

        try:
            # Find pdb file in the output folder
            tmp_pdb = glob.glob(self.tmp_pre+'*.pdb')
            assert tmp_pdb, 'No refined files found: {}'.format(self.tmp_dir)
            tmp_pdb = tmp_pdb[0]
            tmp_mtz = tmp_pdb.replace('.pdb', '.mtz')
            assert os.path.exists(tmp_pdb)
            assert os.path.exists(tmp_mtz)
            # Copy to output folder
            shutil.copy(tmp_pdb, self.out_pdb_file)
            shutil.copy(tmp_mtz, self.out_mtz_file)
            assert os.path.exists(self.out_pdb_file)
            assert os.path.exists(self.out_mtz_file)
            # Write the log to the output log file
            self.cmd.write_output(self.out_log_file)
        except Exception as e:
            print '------------>'
            print e
            print '------------>'
            print self.cmd.output
            print '------------>'
            print self.cmd.error
            print '------------>'
            raise Exception('Failed during refinement')
        finally:
            # Delete temporary directory
            shutil.rmtree(self.tmp_dir)


class refine_phenix(_refiner):

    program = 'phenix.refine'

    _kw_list = [
            'cmd_line_args',
            'n_cycles',
            'strategy',
            ]

    def setup(self):
        """Prepare command object"""
        kw = self.kw_args

        self.validate()

        self.cmd.add_command_line_arguments(self.pdb_file, self.mtz_file)

        if self.cif_files is not None:
            self.cmd.add_command_line_arguments(self.cif_files)

        self.cmd.add_command_line_arguments('output.prefix={}'.format(self.tmp_pre))

        if 'cmd_line_args' in kw:
            self.cmd.add_command_line_arguments(kw['cmd_line_args'])

        if 'strategy' in kw:
            self.cmd.add_command_line_arguments('refine.strategy={}'.format(kw['strategy']))

        if 'n_cycles' in kw:
            self.cmd.add_command_line_arguments('main.number_of_macro_cycles={:d}'.format(kw['n_cycles']))

        return self

    def refine_coordinates_only(self):
        self.kw_args['strategy'] = 'individual_sites'
        return self

    def refine_b_factors_only(self):
        self.kw_args['strategy'] = 'individual_adp'
        return self


class refine_refmac(_refiner):

    program = 'refmac5'

    _kw_list = [
            'cmd_line_args',
            'std_in',
            'n_cycles',
            ]

    def setup(self):
        """Prepare command object"""
        kw = self.kw_args

        self.validate()

        # Input files
        self.cmd.add_command_line_arguments('XYZIN', self.pdb_file)
        self.cmd.add_command_line_arguments('HKLIN', self.mtz_file)

        # Output files
        self.cmd.add_command_line_arguments('XYZOUT', self.tmp_pre+'.pdb')
        self.cmd.add_command_line_arguments('HKLOUT', self.tmp_pre+'.mtz')

        if self.cif_files is not None:
            if len(self.cif_files) > 1:
                cif_file = merge_cif_libraries(self.cif_files, outcif=os.path.join(self.tmp_dir, 'combined.cif'))
            else:
                cif_file = self.cif_files[0]
            assert os.path.exists(cif_file)
            self.cmd.add_command_line_arguments('LIBIN', cif_file)

        # Process STDIN arguments
        if 'std_in' in kw:
            for l in kw['std_in']:
                self.cmd.add_standard_input(l)
        if 'n_cycles' in kw:
            self.cmd.add_standard_input('NCYC {}'.format(n_cycles))

        ###################################
        #          MUST BE AT END         #
        self.cmd.add_standard_input('END')
        #          MUST BE AT END         #
        ###################################

        return self

    def refine_b_factors_only(self):
        self.kw_args.setdefault('std_in',[]).append("REFI BONLY")
        return self

    def refine_coordinates_only(self):
        self.kw_args.setdefault('std_in',[]).append("BREF OVER")
        return self

    def zero_cycles(self):
        self.kw_args['n_cycles'] = 0
        return self

    def complete_model(self):
        self.kw_args.setdefault('std_in',[]).append("BUIL Y")
        return self


class BFactorRefinementFactory(object):

    _refine = refine_phenix

    def __init__(self, pdb_file, mtz_file, out_dir, cif_files=[], tag=None, tls_selections=None, prefix='refined'):

        self.pdb_file = pdb_file
        self.mtz_file = mtz_file
        self.cif_files = cif_files
        self.out_dir = easy_directory(out_dir)
        self.tag = tag
        self.tls_selections = []
        self.tls_matrices = None

        self.initial_pdb = os.path.join(self.out_dir, 'initial.pdb')
        self.out_template = os.path.join(self.out_dir, prefix)

        shutil.copy(self.pdb_file, self.initial_pdb)

        self.log = Log()

        if not tls_selections:
            tls_selections = self.determine_tls_groups(pdb_file=pdb_file)

        # Sanitise the tls selections
        for tls in tls_selections:
            if tls.startswith('"') and tls.endswith('"'):
                tls=tls[1:-1]
            assert '\"' not in tls, 'TLS selection cannot include \": {}'.format(tls)
            self.tls_selections.append(tls)

    def determine_tls_groups(self, pdb_file):
        """Use phenix.find_tls_groups to generate tls groups for the structure"""
        self.log.subheading('Determining TLS groups for: {}'.format(pdb_file))

        tls_selections = phenix_find_tls_groups(pdb_file)

        self.log.subheading('Identified TLS Selections:')
        for s in tls_selections:
            self.log(s)

        return tls_selections

    def refine_b_factors(self, mode='tls', suffix=None):
        """Refine the model with phenix.refine, including the TLS model"""

        assert mode in ['isotropic', 'tls', 'anisotropic']

        if suffix is None: suffix = mode

        strategy = "individual_sites+individual_adp+occupancies"

        if mode == 'isotropic':
            strategy += ''
            args = [r'convert_to_isotropic=True']
        elif mode == 'tls':
            strategy += '+tls'
            args = [r'refinement.refine.adp.tls="{}"'.format(t) for t in self.tls_selections]
        else:
            strategy += ''
            args = [r'refinement.refine.adp.individual.anisotropic="{}"'.format(' or '.join(['('+t+')' for t in self.tls_selections]))]

        self.log.subheading('Refining B-factor model with {}'.format(self._refine.program))
        obj = self._refine(pdb_file=self.pdb_file, mtz_file=self.mtz_file, cif_files=self.cif_files,
                           out_prefix=self.out_template+'-'+suffix,
                           strategy=strategy, n_cycles=3, cmd_line_args=args).run()

        return obj.out_pdb_file, obj.out_mtz_file

    @staticmethod
    def extract_tls_from_pdb(pdb_file):
        ih = iotbx.pdb.hierarchy.input(pdb_file)
        tls_params = ih.input.extract_tls_params(ih.hierarchy)
        return tls_params

    def show_tls_params(self, tls_params=None, pdb_file=None):
        if pdb_file: tls_params=self.extract_tls_from_pdb(pdb_file=pdb_file)
        T = tls_params.tls_params[0].t
        L = tls_params.tls_params[0].l
        S = tls_params.tls_params[0].s

        o = ""
        for tls in tls_params.tls_params:
            o += '\n'
            o += 'selection: {}\n'.format(tls.selection_string)
            o += 'origin: {}\n'.format(tls.origin)
            o += 'T: '+str(tls.t)+'\n'
            o += 'L: '+str(tls.l)+'\n'
            o += 'S: '+str(tls.s)+'\n'
        o += '\n'
        self.log(o)


import giant.logs as lg
logger = lg.getLogger(__name__)

import os

from giant.paths import easy_directory
from giant.exceptions import Sorry, Failure
from giant.dispatcher import Dispatcher
from giant.structure.tls import phenix_find_tls_groups


class RefinementError(Exception):
    pass


class _refiner(object):

    program = None

    def __init__(self, pdb_file, mtz_file=None, cif_files=None, out_prefix=None):

        # Set defaults if not given
        if mtz_file is None:
            mtz_file = pdb_file.replace('.pdb','.mtz')
        if out_prefix is None:
            out_prefix = os.path.splitext(pdb_file)[0] + '-refined'

        # Convert cif files to list
        if isinstance(cif_files, str):
            cif_files = [cif_files]

        # Main files
        self.input_files = dict(
            pdb = pdb_file,
            mtz = mtz_file,
            cif = cif_files,
        )

        # Eventual output prefix and output files
        self.out_prefix = out_prefix
        self.output_files = dict(
            pdb = self.out_prefix+'.pdb',
            mtz = self.out_prefix+'.mtz',
            log = self.out_prefix+'.log',
        )

        # Validate input
        if not os.path.exists(self.input_files['pdb']):
            raise Sorry('File does not exist: {}'.format(self.input_files['pdb']))
        if not os.path.exists(self.input_files['mtz']):
            raise Sorry('File does not exist: {}'.format(self.input_files['mtz']))
        if os.path.exists(self.output_files['pdb']):
            raise Sorry('Output file already exists: {}'.format(self.output_files['pdb']))
        if os.path.exists(self.output_files['mtz']):
            raise Sorry('Output file already exists: {}'.format(self.output_files['mtz']))

        # Create temporary folder for refinement
        import tempfile
        self.tmp_dir = tempfile.mkdtemp(prefix='refine-model-')
        self.tmp_pre = os.path.join(self.tmp_dir, 'refine')

        # Command object for refinement
        self.dispatcher = Dispatcher(self.program)

        self.setup()

    def __str__(self):
        out_str = (
            "* Refinement Parameters *\n\n" +
            self.dispatcher.as_string() + '\n' +
            'Writing output to:\n\t{}\n\t{}\n'.format(
                self.output_files['pdb'],
                self.output_files['mtz'],
            ) +
            "\n* Refinement Command *\n\n" +
            self.dispatcher.as_command()
        )
        return out_str

    def as_string(self):
        return str(self)

    def run(self):
        """...run refinement amd export files"""
        result = self.refine()
        try:
            self.export()
        except RefinementError:
            logger(result.stdout)
            logger(result.stderr)
            raise
        return result

    def add_parameter_file(self):
        raise NotImplementedError('Dummy class -- not implemented')

    def setup(self):
        raise NotImplementedError('Dummy class -- not implemented')

    def finalise(self):
        raise NotImplementedError('Dummy class -- not implemented')

    def refine(self):
        """...run the refinement"""
        self.finalise()
        self.dispatcher.run()
        return self.dispatcher.result

    def export(self):
        """Copy files to output destination"""

        import glob, shutil

        # Find pdb file in the output folder
        tmp_pdbs = glob.glob(self.tmp_pre+'*.pdb')
        if not tmp_pdbs:
            raise RefinementError('No refined files found: {}'.format(self.tmp_dir))
        tmp_pdb = tmp_pdbs[0]
        tmp_mtz = tmp_pdb.replace('.pdb', '.mtz')

        if not os.path.exists(tmp_pdb):
            raise RefinementError('Output PDB file does not exist: {}'.format(tmp_pdb))
        if not os.path.exists(tmp_mtz):
            raise RefinementError('Output MTZ file does not exist: {}'.format(tmp_mtz))

        # Copy to output folder
        shutil.copy(tmp_pdb, self.output_files['pdb'])
        shutil.copy(tmp_mtz, self.output_files['mtz'])

        # Write the log to the output log file
        self.dispatcher.write_output(self.output_files['log'])

        if not (
            os.path.exists(self.output_files['pdb']) and
            os.path.exists(self.output_files['mtz'])
            ):
            raise RefinementError(
                'Output PDB & MTZ files do not exist: {} / {}'.format(
                    self.output_files['pdb'],
                    self.output_files['mtz'],
                )
            )

        if not os.path.exists(self.output_files['log']):
            raise RefinementError(
                'Output log file does not exist: {}'.format(
                    self.output_files['log'],
                )
            )

        # Delete temporary directory
        shutil.rmtree(self.tmp_dir)


class refine_phenix(_refiner):

    program = 'phenix.refine'

    def setup(self):
        """Prepare command object"""

        self.dispatcher.extend_args([
            self.input_files['pdb'],
            self.input_files['mtz'],
        ])

        if self.input_files['cif'] is not None:
            self.dispatcher.extend_args(
                self.input_files['cif'],
            )

        self.dispatcher.append_arg(
            'output.prefix={}'.format(self.tmp_pre),
        )

        return self

    def finalise(self):
        return self

    def add_parameter_file(self, parameter_file):

        # Suffixes recognised by phenix
        if not (
            parameter_file.endswith('.eff') or
            parameter_file.endswith('.def') or
            parameter_file.endswith('.params')
        ):
            raise Sorry('Parameter files must end with .eff or .def or .params')

        # Files can just be passed straight to the program
        self.dispatcher.append_arg(parameter_file)

        return self

    def set_n_cycles(self, n):
        self.dispatcher.append_arg(
            'main.number_of_macro_cycles={:d}'.format(n)
        )
        return self

    def set_refine_coordinates_only(self):
        self.dispatcher.append_arg(
            'refine.strategy=individual_sites'
        )
        return self

    def set_refine_b_factors_only(self):
        self.dispatcher.append_arg(
            'refine.strategy=individual_adp'
        )
        return self


class refine_refmac(_refiner):

    program = 'refmac5'

    def setup(self):
        """Prepare command object"""

        # Input files
        self.dispatcher.extend_args([
            'XYZIN', self.input_files['pdb'],
            'HKLIN', self.input_files['mtz'],
            'XYZOUT', self.tmp_pre+'.pdb',
            'HKLOUT', self.tmp_pre+'.mtz',
        ])

        cif_files = self.input_files['cif']

        if (cif_files is not None) and (len(cif_files) > 0):

            # Refmac can only take one cif file -- merge
            if len(cif_files) > 1:
                from giant.refinement.restraints.merge_cifs import merge_cif_libraries
                cif_file = os.path.splitext(self.output_files['pdb'])[0] + '-combined-restraints.cif'
                merge_cif_libraries(cif_files, outcif=cif_file)
            else:
                cif_file = cif_files[0]
            assert os.path.exists(cif_file)

            self.dispatcher.extend_args(
                ['LIBIN', cif_file]
            )

        return self

    def finalise(self):

        ###################################
        #          MUST BE AT END         #
        self.dispatcher.append_stdin('END')
        #          MUST BE AT END         #
        ###################################

        return self

    def add_parameter_file(self, parameter_file):
        with open(parameter_file, 'r') as fh:
            for line in fh.readlines():
                if line.strip() == 'END':
                    logger('END line found in parameter file -- ignoring lines from this point on')
                    break
                self.dispatcher.append_stdin(line)
        return self

    def set_n_cycles(self, n):
        self.dispatcher.append_stdin(
            'NCYC {:d}'.format(n)
        )
        return self

    def set_refine_coordinates_only(self):
        self.dispatcher.append_stdin(
            "BREF OVER"
        )
        return self

    def set_refine_b_factors_only(self):
        self.dispatcher.append_stdin(
            "REFI BONLY"
        )
        return self

    ################# just refmac ##################

    def zero_cycles(self):
        self.set_n_cycles(0)
        return self

    def complete_model(self):
        self.dispatcher.append_stdin(
            "BUIL Y"
        )
        return self

    ################# just refmac ##################

refiners = {
    'refmac' : refine_refmac,
    'phenix' : refine_phenix,
}

def get_refiner(name):

    try:
        return refiners[name]
    except KeyError:
        raise Failure('Invalid refinement program selected: {}'.format(name))
    raise Failure('Invalid refinement program selected: {}'.format(name))


class BFactorRefinementFactory(object):

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

        import shutil
        shutil.copy(self.pdb_file, self.initial_pdb)

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
        logger.subheading('Determining TLS groups for: {}'.format(pdb_file))

        tls_selections = phenix_find_tls_groups(pdb_file)

        logger.subheading('Identified TLS Selections:')
        for s in tls_selections:
            logger(s)

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

        # Strategy command line arg
        args.append('refine.strategy={}'.format(strategy))

        logger.subheading('Refining B-factor model')
        obj = refine_phenix(
            pdb_file=self.pdb_file,
            mtz_file=self.mtz_file,
            cif_files=self.cif_files,
            out_prefix=self.out_template+'-'+suffix,
        )
        obj.set_n_cycles(3)
        obj.dispatcher.extend_args(args)

        obj.run()

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
        logger(o)


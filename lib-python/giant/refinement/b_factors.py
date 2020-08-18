import giant.logs as lg
logger = lg.getLogger(__name__)

import os

from giant.paths import easy_directory

from giant.exceptions import Sorry, Failure
from giant.structure.tls import phenix_find_tls_groups


class BFactorRefinementFactory(object):

    _n_cycles = 5

    def __init__(self, pdb_file, mtz_file, out_dir, cif_files=[], tls_selections=None, prefix='refined'):

        self.pdb_file = pdb_file
        self.mtz_file = mtz_file
        self.cif_files = cif_files
        self.out_dir = easy_directory(out_dir)
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
        from giant.refinement.wrappers import refine_phenix
        obj = refine_phenix(
            pdb_file=self.pdb_file,
            mtz_file=self.mtz_file,
            cif_files=self.cif_files,
            out_prefix=self.out_template+'-'+suffix,
        )
        obj.set_n_cycles(self._n_cycles)
        obj.dispatcher.extend_args(args)

        obj.run()

        return obj.output_files['pdb'], obj.output_files['mtz']

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


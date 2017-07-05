#!/usr/bin/env ccp4-python

import os, sys, re, glob, shutil, copy, tempfile, gc
import math, re

import pandas, numpy

import libtbx.phil
import libtbx.easy_mp
import iotbx.pdb
import mmtbx.tls.tools

from libtbx.utils import Sorry, Failure
from scitbx.array_family import flex

from bamboo.common.path import easy_directory
from bamboo.common.command import CommandManager
from bamboo.common.logs import Log
from bamboo.stats.cluster import generate_group_idxs

from giant.dataset import CrystallographicModel
from giant.structure.b_factors import occupancy_weighted_average_b_factor
from giant.structure.select import protein, backbone, sidechains

numpy.set_printoptions(threshold=numpy.nan)

EIGHT_PI_SQ = 8*math.pi*math.pi

############################################################################

PROGRAM = 'giant.datasets.b_factor_refine'

DESCRIPTION = """
    Refine the B-factors of a series of datasets using different protocols.
"""

############################################################################

blank_arg_prepend = {'.pdb':'pdb=', '.cif':'cif='}

master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .help = "input pdb file"
        .multiple = True
        .type = str
    cif = None
        .type = str
        .multiple = True
    tls_selections = tls_selections.params
        .help = "define the tls groups (used for all structures)"
        .type = path
        .multiple = False
    labelling = filename *foldername
        .type = choice
        .multiple = False
}
output {
    out_dir = refined-b-factors
        .help = "output directory"
        .type = str
    log_file = refined-b-factors.log
        .type = str
}
settings {
    cpus = 48
        .type = int
        .multiple = False
}
""")

############################################################################

def extract_tls_from_pdb(pdb_file):
    ih = iotbx.pdb.hierarchy.input(pdb_file)
    tls_params = ih.input.extract_tls_params(ih.hierarchy)
    return tls_params

############################################################################


class BFactorRefiner(object):

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

        self.log = Log(verbose=True)

        if not tls_selections:
            tls_selections = self.determine_tls_groups(pdb_file=pdb_file)

        # Sanitise the tls selections
        for tls in tls_selections:
            if tls.startswith('"') and tls.endswith('"'):
                tls=tls[1:-1]
            assert '\"' not in tls, 'TLS selection cannot include \": {}'.format(tls)
            self.tls_selections.append(tls)

    def determine_tls_groups(self, pdb_file):

        self.log.subheading('Determining TLS groups for: {}'.format(pdb_file))

        cmd = CommandManager('phenix.find_tls_groups')
        cmd.add_command_line_arguments(pdb_file)

        cmd.print_settings()
        ret_code = cmd.run()

        if ret_code != 0:
            self.log(cmd.output)
            self.log(cmd.error)
            raise Exception('Failed to determine TLS groups: {}'.format(' '.join(cmd.program)))

        regex = re.compile("refinement\.refine\.adp \{([\s\S]*?)\}")
        tls_command = regex.findall(cmd.output)[0]
        tls_selections = [s.strip() for s in tls_command.split('tls =') if s.strip()]

        self.log.subheading('Identified TLS Selections:')
        for s in tls_selections:
            self.log(s)

        return tls_selections

#    def initial_tls_parameters(self):
#        """Characterise TLS with phenix.tls - legacy function"""
#
#        self.log.subheading('Fitting TLS Matrices to selections')
#        self.log('writing to output file: {}'.format(self.tls_initial_pdb))
#
#        cmd = CommandManager('phenix.tls')
#        cmd.add_command_line_arguments(self.pdb_file)
#        cmd.add_command_line_arguments(self.cif_files)
#        cmd.add_command_line_arguments('extract_tls=True')
#        cmd.add_command_line_arguments([r'selection="{}"'.format(s) for s in self.tls_selections if s is not None])
#        cmd.add_command_line_arguments('output_file_name={}'.format(self.tls_initial_pdb))
#
#        cmd.print_settings()
#        ret_code = cmd.run()
#        cmd.write_output(self.tls_initial_pdb.replace('.pdb', '.log'))
#
#        if ret_code != 0:
#            self.log(cmd.output)
#            self.log(cmd.error)
#            raise Exception('Failed to determine TLS parameters: {}'.format(' '.join(cmd.program)))
#
#        return self.tls_initial_pdb, self.extract_tls_from_pdb(self.tls_initial_pdb)

    def refine_b_factors(self, mode='tls', suffix=None):
        """Refine the model with phenix.refine, including the TLS model"""

        assert mode in ['isotropic', 'tls', 'anisotropic']

        if suffix is None: suffix = mode

        tmp_dir = tempfile.mkdtemp(prefix='b-factor-refine-')
        tmp_pre = os.path.join(tmp_dir, 'refine')

        out_pre = self.out_template+'-'+suffix
        out_pdb = out_pre+'.pdb'
        out_mtz = out_pre+'.mtz'
        out_log = out_pre+'.log'

        if os.path.exists(out_pdb) and os.path.exist(out_mtz):
            self.log('refined structure already exists: {} - skipping'.format(out_pdb))
            return out_pdb, out_mtz

        if mode == 'isotropic':
            strategy = 'individual_adp'
            params = [r'convert_to_isotropic=True']
        elif mode == 'tls':
            strategy = 'tls+individual_adp'
            params = [r'refinement.refine.adp.tls="{}"'.format(t) for t in self.tls_selections]
        else:
            strategy = 'individual_adp'
            params = [r'refinement.refine.adp.individual.anisotropic="{}"'.format(t) for t in self.tls_selections]

        self.log.subheading('Refining B-factor model with phenix.refine')
        self.log('writing to temporary output directory: {}'.format(tmp_dir))

        cmd = CommandManager('phenix.refine')
        cmd.add_command_line_arguments(self.pdb_file, self.mtz_file)
        cmd.add_command_line_arguments(self.cif_files)
        cmd.add_command_line_arguments('refine.strategy='+strategy)
        cmd.add_command_line_arguments('main.number_of_macro_cycles=3')
        cmd.add_command_line_arguments(params)
        cmd.add_command_line_arguments('output.prefix={}'.format(tmp_pre))

        cmd.print_settings()
        ret_code = cmd.run()
        cmd.write_output(out_log)

        if ret_code != 0:
            self.log(cmd.output)
            self.log(cmd.error)
            raise Exception('Failed to determine refine model: {}'.format(' '.join(cmd.program)))

        self.log('copying output to output directory: {}'.format(self.out_dir))

        tmp_pdb = glob.glob(tmp_pre+'*.pdb')[0]
        tmp_mtz = tmp_pdb.replace('.pdb', '.mtz')

        self.log('copying output pdb: {} -> {}'.format(tmp_pdb, self.out_dir))
        shutil.copy(tmp_pdb, out_pdb)
        shutil.copy(tmp_mtz, out_mtz)

        shutil.rmtree(tmp_dir)

        return out_pdb, out_mtz

    @staticmethod
    def extract_tls_from_pdb(pdb_file):
        return extract_tls_from_pdb(pdb_file)

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

def wrapper_run(tls_fit):
    tls_fit.log.heading('Processing: {}'.format(tls_fit.tag))
    tls_fit.refine_b_factors(mode='isotropic')
    tls_fit.refine_b_factors(mode='tls')
    tls_fit.refine_b_factors(mode='anisotropic')
    return tls_fit

def run(params):

    log = Log(verbose=True)
    log.heading('Validating input parameters')

    out_dir = easy_directory(params.output.out_dir)

    if os.path.exists(params.input.tls_selections):
        tls_selections = open(params.input.tls_selections, 'r').read().strip().split('\n')
        log('Using existing TLS selections:')
        log('\t'+'\n\t'.join(tls_selections))
    else:
        tls_selections = None

    log.heading('Fitting TLS Parameters')

    all_fits = []

    for p in params.input.pdb:

        if params.input.labelling == 'foldername':
            tag = os.path.basename(os.path.dirname(p))
        elif params.input.labelling == 'filename':
            tag = os.path.basename(os.path.splitext(p)[0])

        fit = BFactorRefiner(pdb_file=p,
                             mtz_file=p.replace('.pdb', '.mtz'),
                             cif_files=params.input.cif,
                             out_dir=os.path.join(out_dir,tag),
                             tag=tag,
                             tls_selections=tls_selections)

        if tls_selections is None:
            tls_selections=fit.tls_selections

            assert not os.path.exists(params.input.tls_selections)
            with open(params.input.tls_selections, 'w') as fh: fh.write('\n'.join(tls_selections))

        all_fits.append(fit)

    all_fits = libtbx.easy_mp.pool_map(fixed_func=wrapper_run, args=all_fits, processes=params.settings.cpus, chunksize=1)

    log.heading('DONE')

############################################################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
        program             = PROGRAM,
        description         = DESCRIPTION)

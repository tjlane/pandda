#!/usr/bin/env ccp4-python

import os, sys, re, glob, shutil, copy, tempfile

import libtbx.phil
import libtbx.easy_mp
import iotbx.pdb

import mdp
import pandas, numpy

numpy.set_printoptions(threshold=numpy.nan)

from bamboo.common.path import easy_directory
from bamboo.common.command import CommandManager
from bamboo.common.logs import Log

from matplotlib import pyplot
pyplot.interactive(0)
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d

from matplotlib.patches import FancyArrowPatch

############################################################################

PROGRAM = 'giant.datasets.tls_analysis'

DESCRIPTION = """
    Analyse the TLS models of a set of related structures.
"""

############################################################################

blank_arg_prepend = {'.pdb' : 'pdb='}

master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .help = "input pdb file"
        .multiple = True
        .type = str
    tls_selection = None
        .help = "define the tls groups (used for all structures)"
        .type = str
        .multiple = True
}
output {
    out_dir = multi-dataset-tls-characterisation
        .help = "output directory"
        .type = str
    log_file = multi-dataset-tls-characterisation.log
        .type = str
}
settings {
    cpus = 48
        .type = int
        .multiple = False
}
""")

############################################################################

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

class TLSAnalyser(object):

    def __init__(self, pdb_file, mtz_file, out_dir, dataset_id=None, tls_selections=None):

        self.pdb_file = pdb_file
        self.mtz_file = mtz_file
        self.out_dir = out_dir
        self.dataset_id = dataset_id
        self.tls_selections = tls_selections
        self.tls_matrices = None

        self.tls_initial_pdb = os.path.join(self.out_dir, 'initial.pdb')
        self.tls_refined_pdb = os.path.join(self.out_dir, 'refined.pdb')

        self.log = Log(verbose=True)

        if not tls_selections:
            self.tls_selections = self.determine_tls_groups(pdb_file=pdb_file)

        if not os.path.exists(self.out_dir): os.mkdir(self.out_dir)

    def determine_tls_groups(self, pdb_file):

        self.log.heading('Determining TLS groups for: {}'.format(pdb_file))

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

    def initial_tls_parameters(self):
        """Characterise TLS with phenix.tls"""

        self.log.heading('Fitting TLS Matrices to selections')
        self.log('writing to output file: {}'.format(self.tls_initial_pdb))

        cmd = CommandManager('phenix.tls')
        cmd.add_command_line_arguments(self.pdb_file)
        cmd.add_command_line_arguments('extract_tls=True')
        cmd.add_command_line_arguments([r'selection={}'.format(s) for s in self.tls_selections if s is not None])
        cmd.add_command_line_arguments('output_file_name={}'.format(self.tls_initial_pdb))

        cmd.print_settings()
        ret_code = cmd.run()
        cmd.write_output(self.tls_initial_pdb.replace('.pdb', '.log'))

        if ret_code != 0:
            self.log(cmd.output)
            self.log(cmd.error)
            raise Exception('Failed to determine TLS parameters: {}'.format(' '.join(cmd.program)))

        return self.tls_initial_pdb, self.extract_tls_from_pdb(self.tls_initial_pdb)

    def refined_tls_parameters(self):
        """Refine the model with phenix.refine, including the TLS model"""

        out_dir = tempfile.mkdtemp(prefix='tls_refine_')
        out_pre = os.path.join(out_dir, 'tls-refine')

        self.log.heading('Refining TLS model with phenix.refine')
        self.log('writing to output directory: {}'.format(out_dir))

        cmd = CommandManager('phenix.refine')
        cmd.add_command_line_arguments(self.pdb_file, self.mtz_file)
        cmd.add_command_line_arguments('refine.strategy=individual_adp+tls')
        cmd.add_command_line_arguments('main.number_of_macro_cycles=1')
        cmd.add_command_line_arguments([r'refinement.refine.adp.tls={}'.format(t) for t in self.tls_selections])
        cmd.add_command_line_arguments('output.prefix={}'.format(out_pre))

        cmd.print_settings()
        ret_code = cmd.run()
        cmd.write_output(out_pre+'.log')

        if ret_code != 0:
            self.log(cmd.output)
            self.log(cmd.error)
            raise Exception('Failed to determine refine model: {}'.format(' '.join(cmd.program)))

        out_pdb = glob.glob(out_pre+'*.pdb')[0]

        self.log('copying output pdb: {} -> {}'.format(out_pdb, self.tls_refined_pdb))
        shutil.copy(out_pdb, self.tls_refined_pdb)

        #os.rmdir(out_dir)

        return self.tls_refined_pdb, self.extract_tls_from_pdb(self.tls_refined_pdb)

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


class MultiTLS(object):

    _t_name = ['T11','T22','T33','T12','T13','T23']
    _l_name = ['L11','L22','L33','L12','L13','L23']
    _s_name = ['S11','S12','S13','S21','S22','S23','S31','S32','S33']

    def __init__(self, csv_base='tls-params-'):

        self.csv_base = csv_base
        self.tables = {}

        self.log = Log(verbose=True)

    def add(self, tls):

        tls_params = tls.extract_tls_from_pdb(tls.tls_refined_pdb)

        for tls_fit in tls_params.tls_params:
            # Extract data table for this selection
            tls_table = self.tables.setdefault(tls_fit.selection_string, pandas.DataFrame(columns=self._t_name+self._l_name+self._s_name))

            tls_table.loc[tls.dataset_id]=None
            tls_table.loc[tls.dataset_id,self._t_name] = tls_fit.t
            tls_table.loc[tls.dataset_id,self._l_name] = tls_fit.l
            tls_table.loc[tls.dataset_id,self._s_name] = tls_fit.s

    def show(self):
        for selection in sorted(self.tables.keys()):
            self.log.subheading(selection)
            print self.tables[selection]

    def write(self):
        with open(self.csv_base+'selections.log', 'w') as csv_log:
            for i, selection in enumerate(sorted(self.tables.keys())):
                n = i+1
                csv_log.write('selection {:03d} : {}\n'.format(n, selection))
                self.tables[selection].to_csv(self.csv_base+'selection-{:03d}'.format(n)+'.csv')

    def run_pca(self):

        t = self.table.loc[:,self._t_name]
        l = self.table.loc[:,self._l_name]
        s = self.table.loc[:,self._s_name]

        #self.plot_3d(s.values[:,0:3], block=False)
        #self.plot_3d(s.values[:,3:6], block=False)
        #self.plot_3d(s.values[:,5:8], block=False)
        #self.plot_3d(s.values[:,[0,4,8]],block=True)

        data = self.table.values

        pca = mdp.nodes.PCANode()
        pca.train(data)
        pca.stop_training()

        print repr(pca)

        print 'average: ', pca.avg

        print 'projection matrix: '
        print pca.get_projmatrix()

        proj_vals = pca.execute(data)

        print proj_vals.shape

        print pca.get_explained_variance()
        print pca.d

        self.plot_bar(pca.d)
        self.plot_3d(proj_vals[:,:3])
        self.plot_2d(proj_vals[:,:2])

        help(pca)

    def plot_2d(self, data, block=True):
        fig = pyplot.figure(figsize=(8,8))
        ax = fig.add_subplot(111)
        pyplot.rcParams['legend.fontsize'] = 10
        ax.plot(data[:,0], data[:,1], 'o', markersize=8, color='blue', alpha=0.5)
        pyplot.show(block)

    def plot_3d(self, data, block=True):
        fig = pyplot.figure(figsize=(8,8))
        ax = fig.add_subplot(111, projection='3d')
        pyplot.rcParams['legend.fontsize'] = 10
        ax.plot(data[:,0], data[:,1], data[:,2], 'o', markersize=8, color='blue', alpha=0.5)
        pyplot.show(block)

    def plot_bar(self, vals, block=True):
        fig = pyplot.figure(figsize=(8,8))
        ax = fig.add_subplot(111)
        pyplot.rcParams['legend.fontsize'] = 10
        ax.plot(range(1,len(vals)+1), vals, 'o', markersize=8, color='blue', alpha=0.5)
        pyplot.show(block)

def wrapper_run(tls_fit):
    tls_fit.log.heading('Processing: {}'.format(tls_fit.dataset_id))
    if not os.path.exists(tls_fit.tls_initial_pdb):
        pdb1 = tls_fit.initial_tls_parameters()
    if not os.path.exists(tls_fit.tls_refined_pdb):
        pdb2 = tls_fit.refined_tls_parameters()
    return tls_fit

def run(params):

    log = Log(verbose=True)
    log.heading('Validating input parameters')

    out_dir = easy_directory(params.output.out_dir)

    if params.input.tls_selection == [None]:
        params.input.tls_selection = None

    all_fits = []

    for p in params.input.pdb:
        print 'fitting: {}'.format(p)
        fit = TLSAnalyser(pdb_file=p,
                          mtz_file=p.replace('.pdb', '.mtz'),
                          out_dir=os.path.join(out_dir,os.path.basename(os.path.splitext(p)[0])+'-tls'),
                          dataset_id=os.path.basename(os.path.splitext(p)[0]),
                          tls_selections=copy.copy(params.input.tls_selection))

        if not params.input.tls_selection:
            params.input.tls_selection=fit.tls_selections

        all_fits.append(fit)

    all_fits = libtbx.easy_mp.pool_map(fixed_func=wrapper_run, args=all_fits, processes=params.settings.cpus, chunksize=1)

    multi_tls = MultiTLS()

    for fit in all_fits:

        fit.log.heading('Results for {}'.format(fit.dataset_id))
        fit.log.subheading('Initial Fit')
        fit.show_tls_params(pdb_file=fit.tls_initial_pdb)
        fit.log.subheading('Refined Fit')
        fit.show_tls_params(pdb_file=fit.tls_refined_pdb)

        multi_tls.add(fit)

    multi_tls.show()
    multi_tls.write()
#    multi_tls.run_pca()

    from IPython import embed; embed()

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

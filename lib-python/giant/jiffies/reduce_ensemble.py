import os, sys, glob, re

#################################
import matplotlib
matplotlib.use('Agg')
matplotlib.interactive(0)
from matplotlib import pyplot
pyplot.style.use('ggplot')
#################################

import libtbx.phil

import numpy, pandas

import iotbx.pdb

from scitbx.array_family import flex

from giant.structure.ensembles import StructureCollection

#######################################

bar = '=======================++>'

#######################################

blank_arg_prepend = {'.pdb':'pdb1=', '.mtz':'mtz1='}

master_phil = libtbx.phil.parse("""""")

#######################################

def run(params):

    pdbs = sorted(glob.glob('*.pdb'))
    out_dir = 'reduce_ensembles'

    e = StructureCollection.from_files(pdbs)
    e.load_all()

    e.write_tables(out_dir=out_dir)
    e.write_graphs(out_dir=out_dir)

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run(None)
#    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)

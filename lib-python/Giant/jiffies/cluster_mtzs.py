import os, sys, copy, re

#import scipy.cluster
import libtbx.cluster
import libtbx.phil

import numpy

from scitbx.array_family import flex
from Giant.Xray.Data import crystalSummary

blank_arg_prepend = 'mtz='

master_phil = libtbx.phil.parse("""
input {
    mtz = None
        .type = path
        .multiple = True
    mtz_regex = '(x.*?)\/'
        .type = str
}
summary {
    res_cutoff = 4
        .type = float
    cluster_cutoff = 10
        .type = float
}
output {
    out_dir = None
        .type = path
}
""")

def comp_angles(m1, m2):
    uc1 = m1.unit_cell.parameters()[3:]
    uc2 = m2.unit_cell.parameters()[3:]
    return (flex.double(uc1) - flex.double(uc2)).norm()
def comp_axes(m1, m2):
    uc1 = m1.unit_cell.parameters()[:3]
    uc2 = m2.unit_cell.parameters()[:3]
    return (flex.double(uc1) - flex.double(uc2)).norm()

def run(params):

    mtz_summaries = []

    if params.input.mtz_regex: mtz_regex = params.input.mtz_regex
    else:                      mtz_regex = None

    if params.output.out_dir:
        if not os.path.exists(params.output.out_dir): os.mkdir(params.output.out_dir)

    #########################################################################################################
    print '=====================================>'
    print '=====================================>'
    print 'Reading Input Data'
    print '=====================================>'
    print '=====================================>\n'
    #########################################################################################################

    # Find resolutions of the datasets
    for i, m in enumerate(params.input.mtz):

        print 'READING:', m

        # Create a name for the dataset
        if mtz_regex: m_id = re.findall(mtz_regex, m)[0]
        else:         m_id = 'MTZ-{:06d}'.format(i)
        # Create a summary of the dataset
        summ = crystalSummary(m, id=m_id)

        # Append to list
        mtz_summaries.append(summ)

    #########################################################################################################

    print '=====================================>'
    print 'CLUSTERING UNIT CELLS'

    print '=====================================>'
    print 'CLUSTERING UNIT CELL ANGLES'
    cl = libtbx.cluster.HierarchicalClustering(mtz_summaries, lambda x,y: comp_angles(x,y))
    clusters = cl.getlevel(10)

    cl.data[0].display()

    print len(clusters), 'CLUSTERS'
    print 'SIZES', map(len, clusters)
    #assert len(clusters) == 1

    print 'CLUSTERING UNIT CELL AXES'
    cl = libtbx.cluster.HierarchicalClustering(mtz_summaries, lambda x,y: comp_axes(x,y))
    clusters = cl.getlevel(params.summary.cluster_cutoff)

    cl.data[0].display()

    print len(clusters), 'CLUSTERS'
    print 'SIZES', map(len, clusters)

    #########################################################################################################

    print '=====================================>'
    print 'CREATING OUTPUT CLUSTER DIRECTORIES'

    if params.output.out_dir:

        for i_clust, clust in enumerate(clusters):
            print '=====================================>'
            print 'OUTPUTTING CLUSTER {}'.format(i_clust)
            sub_dir = os.path.join(params.output.out_dir, 'CLUSTER_{}'.format(i_clust))
            if not os.path.exists(sub_dir): os.mkdir(sub_dir)
            for s in clust:
                sub_sub_dir = os.path.join(sub_dir, s.id)
                if not os.path.exists(sub_sub_dir): os.mkdir(sub_sub_dir)

                os.symlink(os.path.abspath(s.mtz_file), os.path.join(sub_sub_dir, s.id+'.mtz'))
                potential_pdb = s.mtz_file.replace('.mtz', '.pdb')
                if os.path.exists(potential_pdb):
                    os.symlink(os.path.abspath(potential_pdb), os.path.join(sub_sub_dir, s.id+'.pdb'))







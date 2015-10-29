import os, sys, copy, re

import libtbx.phil

import numpy
import pandas

from Giant.Xray.Data.cluster import sort_mtzs_by_spacegroup

blank_arg_prepend = 'mtz='

master_phil = libtbx.phil.parse("""
input {
    mtz = None
        .type = path
        .multiple = True
    mtz_regex = None
        .type = str
}
output {
    out_dir = None
        .type = path
}
""")

def run(params):

    if params.output.out_dir and (not os.path.exists(params.output.out_dir)):
        os.mkdir(params.output.out_dir)
    # Dump images into the current directory if no directory is given
    if not params.output.out_dir: img_dir = './'
    else: img_dir = params.output.out_dir

    #########################################################################################################

    if params.input.mtz_regex: labels = [re.findall(params.input.mtz_regex, m)[0] for m in params.input.mtz]
    else:                      labels = ['MTZ-{:06d}'.format(i) for i in range(len(params.input.mtz))]

    #########################################################################################################

    print '============================================>'
    print 'GROUPING BY SPACE GROUP'

    crystal_groups = sort_mtzs_by_spacegroup(mtz_files=params.input.mtz, labels=labels)

    #########################################################################################################

    for cg in crystal_groups:

        sg_name = 'SG-{}'.format(cg.space_groups[0].replace(' ','').replace('(','-').replace(')',''))

        print '==============================================================================>'
        print 'SPACE GROUP {}: {} DATASETS'.format(cg.space_groups[0], len(cg.crystals))

        #########################################################################################################

        res_sorted = sorted(cg.crystals, key=lambda c: c.high_res)
        high_res_crystal = res_sorted[0]
        low_res_crystal  = res_sorted[-1]

        print '===========================>'
        print 'HIGHEST RESOLUTION: {:.2f}'.format(high_res_crystal.high_res)
        print 'DATASET + PATH: {} - {}'.format(high_res_crystal.id, high_res_crystal.mtz_file)
        print 'LOWEST RESOLUTION: {:.2f}'.format(low_res_crystal.high_res)
        print 'DATASET + PATH: {} - {}'.format(low_res_crystal.id, low_res_crystal.mtz_file)

        #########################################################################################################

        print '===========================>'
        print 'UNIT CELL VARIATION ANALYSIS'
        print numpy.round(cg.uc_stats.as_pandas_table().T, 2)

        #########################################################################################################

        print '===========================>'
        print 'UNIT CELL DENDROGRAMS'
        if len(cg.crystals) > 1:
            cg.dendrogram(  fname = os.path.join(img_dir,'{}-dendrogram.png'.format(sg_name)),
                            xlab  = 'Crystal',
                            ylab  = 'Linear Cell Variation')

        #########################################################################################################

        for i_cg2, cg2 in enumerate(cg.by_unit_cell(cg.crystals,cutoff=0.2)):

            cluster_name = '{}-Cluster-{}'.format(sg_name, i_cg2+1)

            print '===========================>'
            print cluster_name

            #########################################################################################################

            if len(cg2.crystals) > 1:
                cg2.dendrogram( fname = os.path.join(img_dir, '{}-dendrogram.png'.format(cluster_name)),
                                xlab  = 'Crystal',
                                ylab  = 'Linear Cell Variation',
                                ylim  = (0,0.2),
                                annotate_y_min = 0.2)

            #########################################################################################################

            if params.output.out_dir:
                print '===========================>'
                print 'CREATING OUTPUT DIRECTORIES'

                # Go through and link the datasets for each of the spacegroups into a separate folder
                sub_dir = os.path.join(params.output.out_dir, cluster_name)
                if not os.path.exists(sub_dir): os.mkdir(sub_dir)

                for c in cg2.crystals:
                    print 'Linking {}'.format(c.mtz_file)

                    sub_sub_dir = os.path.join(sub_dir, c.id)
                    if not os.path.exists(sub_sub_dir): os.mkdir(sub_sub_dir)
                    os.symlink(os.path.abspath(c.mtz_file), os.path.join(sub_sub_dir, c.id+'.mtz'))

                    potential_pdb = c.mtz_file.replace('.mtz','.pdb')
                    if os.path.exists(potential_pdb):
                       os.symlink(os.path.abspath(potential_pdb), os.path.join(sub_sub_dir, c.id+'.pdb'))

def do_a_graph():
    from ascii_graph import Pyasciigraph
    g = Pyasciigraph()
    resolution_shells = numpy.arange(min(all_resolutions)-0.15, max(all_resolutions)+0.15, 0.1)
    graph_data = [('{!s}-{!s}'.format(shell, shell+0.1), sum([1 for r in all_resolutions if (r<shell+0.1 and r>shell)])) for shell in resolution_shells]
    for l in g.graph(label='Resolution Distribution', data=graph_data, sort=0):
        print(l)

    print '-------------------------------->'
    print 'Total Datasets:', len(all_resolutions)
    print '-------------------------------->'
    print len(high_res_datasets), 'Datasets above', min_res
    print '-------------------------------->'
    print ','.join(high_res_datasets)
    print '-------------------------------->'
    print 'Range:', min(all_resolutions), '->', max(all_resolutions)
    print '-------------------------------->'
    print 'Highest Resolution:', highest_res
    print 'Dataset:', highest_res_dataset
    print '-------------------------------->'

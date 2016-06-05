import os, sys, copy, re

import libtbx.phil

import numpy
import pandas

from giant.xray.data import CrystalSummary
from giant.xray.data.cluster import CrystalGroup

#######################################

blank_arg_prepend = {'.mtz':'mtz=', '.pdb':'pdb='}

master_phil = libtbx.phil.parse("""
input {
    mtz = None
        .type = path
        .multiple = True
    mtz_regex = None
        .type = str

    pdb = None
        .type = path
        .multiple = True
    pdb_regex = None
        .type = str
}
clustering{
    lcv_cutoff = 0.2
        .type = float
    label_nodes_above = 0.2
        .help = "Label nodes above this value"
        .type = float
}
output {
    out_dir = None
        .type = path
    split_pdbs_and_mtzs = True
        .type = bool
    html_out = None
        .type = path
}
""")

#######################################

bar = '============================================>'

#######################################

def run(params):

    if not (params.input.pdb or params.input.mtz): raise Exception('No pdb/mtz files have been provided.')

    # Make output directory if required
    if params.output.out_dir and (not os.path.exists(params.output.out_dir)):
        os.mkdir(params.output.out_dir)
    # Dump images into the current directory if no directory is given
    if not params.output.out_dir: img_dir = '.'
    else: img_dir = params.output.out_dir

    #########################################################################################################

    print bar
    print 'Making pdb/mtz dataset labels'

    try:
        if params.input.pdb_regex: p_labels = [re.findall(params.input.pdb_regex, f)[0] for f in params.input.pdb]
        else:                      p_labels = ['PDB-{:06d}'.format(i) for i in range(len(params.input.pdb))]
        if params.input.mtz_regex: m_labels = [re.findall(params.input.mtz_regex, f)[0] for f in params.input.mtz]
        else:                      m_labels = ['MTZ-{:06d}'.format(i) for i in range(len(params.input.mtz))]
    except:
        print 'Error parsing: {}'.format(f)
        raise

    #########################################################################################################

    print bar
    print 'Reading {} pdb(s) and {} mtz(s)'.format(len(params.input.pdb), len(params.input.mtz))

    if params.input.pdb: pdb_summaries = [CrystalSummary.from_pdb(pdb_file=f, id=lab) for f,lab in zip(params.input.pdb, p_labels)]
    else:                pdb_summaries = []
    if params.input.mtz: mtz_summaries = [CrystalSummary.from_mtz(mtz_file=f, id=lab) for f,lab in zip(params.input.mtz, m_labels)]
    else:                mtz_summaries = []

    print bar
    print 'Grouping crystals by space group...'

    # Group by SpaceGroup
    crystal_groups = CrystalGroup.by_space_group(crystals=pdb_summaries+mtz_summaries)

    print '...grouped crystals into {} groups'.format(len(crystal_groups))

    #########################################################################################################

    print bar
    print 'Analysing variation of unit cells for each space group'

    for cg in crystal_groups:

        sg_name = 'SG-{}'.format(cg.space_groups[0].replace(' ','').replace('(','-').replace(')',''))

        print bar
        print ''
        print bar
        print 'Space Group {}: {} dataset(s)'.format(cg.space_groups[0], len(cg.crystals))

        #########################################################################################################

        print '===========================>'
        print 'Unit Cell Variation:'
        print numpy.round(cg.uc_stats.as_pandas_table().T, 2)

        #########################################################################################################

        print '===========================>'
        print 'Making unit cell dendrogram for all crystals with this spacegroup'
        if len(cg.crystals) > 1:
            cg.dendrogram(  fname = os.path.join(img_dir,'{}-dendrogram.png'.format(sg_name)),
                            xlab  = 'Crystal',
                            ylab  = 'Linear Cell Variation')

        #########################################################################################################

        print '===========================>'
        print 'Clustering unit cells...'

        sg_crystal_groups = cg.by_unit_cell(cg.crystals, cutoff=params.clustering.lcv_cutoff)

        print '...clustered crystals into {} groups'.format(len(sg_crystal_groups))

        for i_cg2, cg2 in enumerate(sg_crystal_groups):

            cluster_name = '{}-Cluster-{}'.format(sg_name, i_cg2+1)

            print '===========================>'
            print 'Processing: {}'.format(cluster_name)

            #########################################################################################################

            print 'Making unit cell dendrogram for this cluster of crystals'

            if len(cg2.crystals) > 1:
                cg2.dendrogram( fname = os.path.join(img_dir, '{}-dendrogram.png'.format(cluster_name)),
                                xlab  = 'Crystal',
                                ylab  = 'Linear Cell Variation',
                                ylim  = (0, params.clustering.lcv_cutoff),
                                annotate_y_min = params.clustering.label_nodes_above)

            #########################################################################################################

            if params.output.out_dir:
                print 'Making and populating output directory'

                # Go through and link the datasets for each of the spacegroups into a separate folder
                sub_dir = os.path.join(params.output.out_dir, cluster_name)
                if not os.path.exists(sub_dir): os.mkdir(sub_dir)

                if params.output.split_pdbs_and_mtzs:
                    # Split the mtzs and pdbs into separate directories
                    mtz_dir = os.path.join(sub_dir, 'mtzs')
                    if not os.path.exists(mtz_dir): os.mkdir(mtz_dir)
                    pdb_dir = os.path.join(sub_dir, 'pdbs')
                    if not os.path.exists(pdb_dir): os.mkdir(pdb_dir)
                else:
                    mtz_dir = pdb_dir = sub_dir

                for c in cg2.crystals:
                    if c.mtz_file:
                        # Create subdirectory
                        sub_sub_dir = os.path.join(mtz_dir, c.id)
                        if not os.path.exists(sub_sub_dir): os.mkdir(sub_sub_dir)
                        os.symlink(os.path.abspath(c.mtz_file), os.path.join(sub_sub_dir, c.id+'.mtz'))
                        # Link PDB as well if filenames are the same
                        potential = c.mtz_file.replace('.mtz','.pdb')
                        if os.path.exists(potential):
                            os.symlink(os.path.abspath(potential), os.path.join(sub_sub_dir, c.id+'.pdb'))
                    if c.pdb_file:
                        # Create subdirectory
                        sub_sub_dir = os.path.join(pdb_dir, c.id)
                        if not os.path.exists(sub_sub_dir): os.mkdir(sub_sub_dir)
                        os.symlink(os.path.abspath(c.pdb_file), os.path.join(sub_sub_dir, c.id+'.pdb'))
                        # Link PDB as well if filenames are the same
                        potential = c.pdb_file.replace('.pdb','.mtz')
                        if os.path.exists(potential):
                            os.symlink(os.path.abspath(potential), os.path.join(sub_sub_dir, c.id+'.mtz'))

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)

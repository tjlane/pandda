import os, sys, copy, re

import libtbx.phil

import numpy

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
}
output {
    out_dir = None
        .type = path
}
""")

def run(params):

    all_mtz_summaries = []

    if params.input.mtz_regex: mtz_regex = params.input.mtz_regex
    else:                      mtz_regex = None

    if params.output.out_dir:
        if not os.path.exists(params.output.out_dir): os.mkdir(params.output.out_dir)

    #########################################################################################################

    print '============================================>'
    print 'Reading Input Data'
    print '============================================>'

    # Find resolutions of the datasets
    for i, m in enumerate(params.input.mtz):

        print 'READING:', m

        # Create a name for the dataset
        if mtz_regex: m_id = re.findall(mtz_regex, m)[0]
        else:         m_id = 'MTZ-{:06d}'.format(i)
        # Create a summary of the dataset
        summ = crystalSummary(m, id=m_id)

        # Append to list
        all_mtz_summaries.append(summ)

    #########################################################################################################

    print '============================================>'
    print 'GROUPING BY SPACE GROUP'

    # Pull out all of the space groups
    space_groups = [s.space_group.info().symbol_and_number() for s in all_mtz_summaries]
    print '============================================>'
    print '{!s} UNIQUE SPACE GROUPS (FROM {} DATASETS)'.format(len(set(space_groups)), len(all_mtz_summaries))
    print '============================================>'
    for sg in sorted(set(space_groups)):
        print '{:>20s} -> {:<3d}'.format(sg, space_groups.count(sg))
    # Go through and show summaries of each space group
    for sel_sg in sorted(set(space_groups)):
        # Get all of the mtzs with this space group
        sel_mtzs = [s for s in all_mtz_summaries if s.space_group.info().symbol_and_number()==sel_sg]
        print '==============================================================================>'
        print 'SPACE GROUP {}: {} DATASETS'.format(sel_sg, len(sel_mtzs))

        #########################################################################################################

        # Summarise the resolutions of mtzs for this sg
        high_res_limits = [s.high_res for s in sel_mtzs]
        # Find highest
        highest_res_dataset_idx = high_res_limits.index(min(high_res_limits))
        highest_res_dataset = sel_mtzs[highest_res_dataset_idx]
        # Find lowest
        lowest_res_dataset_idx = high_res_limits.index(max(high_res_limits))
        lowest_res_dataset = sel_mtzs[lowest_res_dataset_idx]

        print '===========================>'
        print 'HIGHEST RESOLUTION: {:.2f}'.format(highest_res_dataset.high_res)
        print 'DATASET + PATH: {} - {}'.format(highest_res_dataset.id, highest_res_dataset.mtz_file)
        print 'LOWEST RESOLUTION: {:.2f}'.format(lowest_res_dataset.high_res)
        print 'DATASET + PATH: {} - {}'.format(lowest_res_dataset.id, lowest_res_dataset.mtz_file)

        #########################################################################################################

        print '===========================>'
        print 'UNIT CELL VARIATION ANALYSIS'

        # Unit Cell Variation
        unit_cell_params = zip(*[s.unit_cell.parameters() for s in sel_mtzs])

        cols = ['VAR', 'MEAN', 'STD', '(%MEAN)', 'RANGE', 'MIN', '(dMEAN)', 'MAX', '(dMEAN)']
        bars = [8*'-']*len(cols)
        head = ['{:^8s}']*len(cols)

        out_table = []
        # Print top of table
        out_table.append(('--'+'---'.join(head)+'--').format(*bars))
        # Print headers
        out_table.append(('| '+' | '.join(head)+' |').format(*cols))
        # Print bottom of headers
        out_table.append(('| '+' | '.join(head)+' |').format(*bars))

        d_pr = '{:>8.2f}'
        p_pr = '{:>8.2%}'

        for param_vals, param_name in zip(unit_cell_params, ['a','b','c','alpha','beta','gamma']):
            param_mean = numpy.mean(param_vals)
            param_var  = numpy.std(param_vals)
            param_var_frac = param_var/param_mean
            param_min = min(param_vals)
            param_min_frac = param_min/param_mean
            param_min_diff = param_min - param_mean
            param_max = max(param_vals)
            param_max_frac = param_max/param_mean
            param_max_diff = param_max - param_mean
            param_range = param_max - param_min

            vars = [param_name, param_mean, param_var, param_var_frac, param_range, param_min, param_min_diff, param_max, param_max_diff]
            form = ['{:>8s}',   d_pr,       d_pr,      p_pr,           d_pr,        d_pr,      d_pr,           d_pr,      d_pr]
            out_table.append(('| '+' | '.join(form)+' |').format(*vars))
        # Print bottom of table
        out_table.append(('--'+'---'.join(head)+'--').format(*bars))

        print '\n'.join(out_table)
        print ''

        #########################################################################################################

        print '===========================>'
        print 'CREATING OUTPUT DIRECTORIES'

        if params.output.out_dir:
            # Go through and link the datasets for each of the spacegroups into a separate folder
            sub_dir = os.path.join(params.output.out_dir, sel_sg.replace(' ','').replace('(','-').replace(')',''))
            if not os.path.exists(sub_dir): os.mkdir(sub_dir)
            for s in sel_mtzs:
                sub_sub_dir = os.path.join(sub_dir, s.id)
                if not os.path.exists(sub_sub_dir): os.mkdir(sub_sub_dir)

                os.symlink(os.path.abspath(s.mtz_file), os.path.join(sub_sub_dir, s.id+'.mtz'))
                potential_pdb = s.mtz_file.replace('.mtz','.pdb')
                if os.path.exists(potential_pdb):
                   os.symlink(os.path.abspath(potential_pdb), os.path.join(sub_sub_dir, s.id+'.pdb'))

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

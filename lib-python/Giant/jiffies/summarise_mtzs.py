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
    out_dir = ./
        .type = path
}
""")

def run(params):

    mtz_summaries = []

    if params.input.mtz_regex: mtz_regex = params.input.mtz_regex
    else:                      mtz_regex = None

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
    print '\n=====================================>'
    print '=====================================>'
    print 'Processing Input Data'
    print '=====================================>'
    print '=====================================>\n'
    #########################################################################################################

    print '\n=====================================>'
    print 'RESOLUTION ANALYSIS'
    print '=====================================>\n'

    # Resolution Variation
    high_res_limits = [s.high_res for s in mtz_summaries]
    # Find highest
    highest_res_dataset_idx = high_res_limits.index(min(high_res_limits))
    highest_res_dataset = mtz_summaries[highest_res_dataset_idx]
    # Find lowest
    lowest_res_dataset_idx = high_res_limits.index(max(high_res_limits))
    lowest_res_dataset = mtz_summaries[lowest_res_dataset_idx]

    print 'HIGHEST RESOLUTION: {:.2f}'.format(highest_res_dataset.high_res)
    print 'HIGH RES DATASET:', highest_res_dataset.id
    print 'HIGH RES PATH:', highest_res_dataset.mtz_file
    print ''
    print 'LOWEST RESOLUTION: {:.2f}'.format(lowest_res_dataset.high_res)
    print 'LOW RES DATASET:', lowest_res_dataset.id
    print 'LOW RES PATH:', lowest_res_dataset.mtz_file

    #########################################################################################################
    # Calculate the variations in the datasets
    #########################################################################################################

    print '\n=====================================>'
    print 'UNIT CELL VARIATION ANALYSIS'
    print '=====================================>\n'

    # Unit Cell Variation
    unit_cell_params = zip(*[s.unit_cell.parameters() for s in mtz_summaries])

    cols = ['VAR', 'MEAN', 'STD', '(%MEAN)', 'MIN', '(dMEAN)', 'MAX', '(dMEAN)']
    bars = [8*'-']*len(cols)
    head = ['{:^8s}']*len(cols)

    print ('--'+'---'.join(head)+'--').format(*bars)
    print ('| '+' | '.join(head)+' |').format(*cols)
    print ('| '+' | '.join(head)+' |').format(*bars)

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

        vars = [param_name, param_mean, param_var, param_var_frac, param_min, param_min_diff, param_max, param_max_diff]
        form = ['{:>8s}',   d_pr,       d_pr,      p_pr,           d_pr,      d_pr,           d_pr,      d_pr]
        print ('| '+' | '.join(form)+' |').format(*vars)

    print ('--'+'---'.join(head)+'--').format(*bars)

    #########################################################################################################

    print '\n=====================================>'
    print 'SPACE GROUP VARIATION ANALYSIS'
    print '=====================================>\n'

    # Space Group Variation
    space_groups = [s.space_group.info().symbol_and_number() for s in mtz_summaries]

    print '=====================================>'
    print '{!s} UNIQUE SPACE GROUPS'.format(len(set(space_groups)))
    print '=====================================>'
    for sg in sorted(set(space_groups)):
        print '{:>20s} -> {:<3d}'.format(sg, space_groups.count(sg))


    print ''
    print ''
    print ''
    print ''
    print ''


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

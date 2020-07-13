#!/usr/bin/env cctbx.python

import giant.logs as lg
logger = lg.getLogger(__name__)

from giant.exceptions import Sorry, Failure

import os, sys, json, collections
import math

import numpy as np

from pandemic.adp.plots import PlotHelper
helper = PlotHelper()

#######################################

blank_arg_prepend = {
    '.pdb' : 'input.pdb=',
    '.json' : 'input.json=',
}

import libtbx.phil
master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .type = path
        .multiple = True
    json = None
        .type = path
        .multiple = False
}
""")

#######################################

def splice_chain_id(pdb_str, chain_id):
    chain_id = str(chain_id)
    assert len(chain_id) == 1
    return ( pdb_str[:9] + chain_id + pdb_str[10:] )

def get_chains_from_labels(labels):
    return [l[9] for l in labels]

def get_residue_number_from_labels(labels):
    return [int(l[10:14]) for l in labels]

def sort_and_average(labels, values):

    uniq_x, y_arr = np.unique(labels, return_inverse=True)
    uniq_y = np.zeros_like(uniq_x, dtype=float)
    for i in range(len(uniq_y)):
        vals = values[y_arr == i]
        uniq_y[i] = np.mean(vals)
    return uniq_x, uniq_y

def make_plot(
    data_frame,
    output_file,
    ):

    import numpy as np
    from matplotlib import pyplot

    logger('> {}'.format(output_file))

    fig, axes = pyplot.subplots(nrows=1, ncols=2)

    ############

    ax1 = axes[0]

    x = data_frame['b_values']
    y = data_frame['ref_values']

    ax1.scatter(
        x = x,
        y = y,
        s = 1,
    )

    non_zero_bool = np.logical_not((x==0) | (y==0))
    corr = round(np.corrcoef(x[non_zero_bool],y[non_zero_bool])[0][1], 3)
    min_x, max_x = ax1.get_xlim()
    min_y, max_y = ax1.get_ylim()

    ax1.text(
        x = (0.05*min_x + 0.95*max_x),
        y = (0.05*min_y + 0.95*max_y),
        s = corr,
        horizontalalignment = 'right',
        verticalalignment = 'bottom',
    )

    ############

    ax2_1 = axes[1]
    ax2_2 = ax2_1.twinx()

    new_x = get_residue_number_from_labels(data_frame['labels'])

    ax2_1.scatter(
        x = new_x,
        y = data_frame['ref_values'],
        s = 1,
        color = 'b',
    )

    uniq_x_1, uniq_y_1 = sort_and_average(
        labels = new_x,
        values = data_frame['ref_values'].values,
    )
    ax2_1.plot(uniq_x_1, uniq_y_1, 'b-')

    ax2_2.scatter(
        x = new_x,
        y = data_frame['b_values'],
        s = 1,
        color = 'r',
    )
    uniq_x_2, uniq_y_2 = sort_and_average(
        labels = new_x,
        values = data_frame['b_values'].values,
    )
    ax2_2.plot(uniq_x_2, uniq_y_2, 'r-')

    non_zero_bool = np.logical_not((uniq_y_1==0) | (uniq_y_2==0))
    corr = round(np.corrcoef(uniq_y_1[non_zero_bool],uniq_y_2[non_zero_bool])[0][1], 3)
    min_x, max_x = ax2_1.get_xlim()
    min_y, max_y = ax2_1.get_ylim()

    ax2_1.text(
        x = (0.95*min_x + 0.05*max_x),
        y = (0.05*min_y + 0.95*max_y),
        s = corr,
        horizontalalignment = 'left',
        verticalalignment = 'bottom',
    )


    ############

    helper.write_and_close_fig(
        fig = fig,
        filename = output_file,
    )

def run(params):

    lg.setup_logging(
        name = __name__,
    )

    with open(params.input.json, 'r') as fh:
        json_dict = json.loads(fh.read())

    json_chains = sorted(set(get_chains_from_labels(json_dict.keys())))

    for pdb_file in params.input.pdb:

        logger.heading(pdb_file)

        pdb_root = os.path.splitext(pdb_file)[0]

        import iotbx.pdb
        h = iotbx.pdb.hierarchy.input(pdb_file).hierarchy

        for json_chain_id in json_chains:

            for c in h.chains():

                if not (c.is_protein() or c.is_na()):
                    continue

                plot_values = []

                for a in c.atoms():

                    a_lab = a.pdb_label_columns()
                    a_b = a.b

                    json_lab = splice_chain_id(a_lab, chain_id=json_chain_id)

                    a_val = json_dict.get(json_lab)

                    if (a_val is None):
                        continue

                    plot_values.append((a_lab, a_b, a_val))

                import pandas
                data_frame = pandas.DataFrame(
                    data = plot_values,
                    columns = ['labels', 'ref_values', 'b_values'])

                logger(data_frame)

                make_plot(
                    data_frame = data_frame,
                    output_file = (pdb_root + '_pdb_chain_{}_json_chain_{}'.format(c.id, json_chain_id) + '.png'),
                )

if __name__ == '__main__':
    from giant.jiffies import run_default
    run_default(
        run = run,
        master_phil = master_phil,
        args = sys.argv[1:],
        blank_arg_prepend = blank_arg_prepend,
    )

import giant.logs as lg
logger = lg.getLogger(__name__)

import os, glob

from pandemic.adp.plots import PlotHelper
helper = PlotHelper()

def get_values_from_structure(
    hierarchy,
    ):

    b_values = []

    for c in hierarchy.chains():
        if not c.is_protein():
            continue

        b_values.extend(list(c.atoms().extract_b()))

    return numpy.array(b_values)

def get_level_bmeans(pandemic_dir):

    import pandas

    csv = os.path.join(pandemic_dir, 'optimisation', 'tracking_levels.csv')

    csv_data = pandas.read_csv(csv)

    last_cycle = max(csv_data.cycle)

    last_cycle_data = csv_data[(csv_data.cycle == last_cycle) & (csv_data.step == 'end')]

    b_vals_tab = last_cycle_data[['level', 'b_avg']].set_index('level').T
    b_vals_dict = collections.OrderedDict([(l, tab) for l, tab in list(b_vals_tab.iteritems())])

    return b_vals_dict

def plot_stuff(data_dict):

    import pandas

    d_frame = pandas.DataFrame(columns=('Protein', 'Level', 'Delta-B', 'Delta-B%'))

    for prot_name, prot_data in data_dict.iteritems():

        pdb_ids = [d[0] for d in prot_data]
        prot_values = [d[1] for d in prot_data]

        levels = prot_values[0].keys()

        # Totals
        b1 = sum([v['b_avg'] for v in prot_values[0].values()])
        b2 = sum([v['b_avg'] for v in prot_values[1].values()])

        d_frame.loc[len(d_frame)] = {
            'Protein' : prot_name,
            'Level' : 'Total',
            'Delta-B' : b2 - b1,
            'Delta-B%' : (b2 - b1) / b1,
        }

        # Levels
        for l in levels:
            b1 = prot_values[0][l]['b_avg']
            b2 = prot_values[1][l]['b_avg']

            d_frame.loc[len(d_frame)] = {
                'Protein' : prot_name,
                'Level' : l,
                'Delta-B' : b2 - b1,
                'Delta-B%' : (b2 - b1) / b1,
            }

    import seaborn as sns
    ax = sns.barplot(
        y = "Protein",
        x = "Delta-B",
        hue = "Level",
        data = d_frame,
    )

    ax.hlines(
        y = [n + 0.5 for n in range(len(data_dict.keys())-1)],
        xmin = ax.get_xlim()[0],
        xmax = ax.get_xlim()[1],
        linewidth = 0.25,
    )

    ax.legend(
        bbox_to_anchor = (1.0, 1.0),
        loc = 'lower right',
        fontsize = 6,
        ncol = 2,
    )

    ax.set_xticks(range(-5, 1+int(ax.get_xlim()[1]), 1), minor=True)
    ax.set_xticks(range(-5, 1+int(ax.get_xlim()[1]), 5))

    ax.yaxis.grid(True, which='minor')

    ax.tick_params(axis='x', labelsize=6)
    ax.tick_params(axis='y', labelsize=6)

    ax.set_aspect(2/1.)

    from matplotlib import pyplot as plt
    #plt.setp(ax.get_xticklabels(), rotation=90)
    helper.write_and_close_fig(
        fig = plt.gcf(),
        filename = 'deltab.png',
    )

######################

cryo_datasets = [
    '1wmd',
    '2wlt',
    '1b2x',
    '2ftl',
    '1a6g',
    '1jxt',
    '1gcs',
    '3k0p',
    '3k0m',
    '1a6n',
    '2fdn',
    '3djj',
    '3lzt',
    '1bx7',
    '1a6k',
    '1x6q',
    '1pnc',
    '1lni',
    '1i0v',
    '3kyv',
    '1jcv',
    '1dqi',
    '1ly0',
    '1tgt',
    '1gdq',
    '1i1w',
    '2dfb',
    '1lw9',
    '3dmv',
    '1ctq',
]

logger = lg.setup_logging(
    name = __name__,
)

in_dir = './temperature-comparison'
out_dir = './temperature-comparison-analysis'

import collections
data_dict = collections.OrderedDict()

for d in sorted(
        glob.glob(os.path.join(in_dir, '*'))
    ):

    prot_label = os.path.basename(d)

    logger.subheading(prot_label)

    pandemic_dirs = sorted(
        glob.glob(
            os.path.join(d, '*', 'pdb_redo', 'pandemic-adp')
        )
    )

    if len(pandemic_dirs) != 2:
        continue

    prot_data = []

    try:
        for run_dir in pandemic_dirs:
            pdb_id = os.path.basename(run_dir.replace('/pdb_redo/pandemic-adp', ''))
            logger(pdb_id)
            mean_b_values = get_level_bmeans(run_dir)
            if pdb_id in cryo_datasets:
                prot_data.insert(0, (pdb_id, mean_b_values))
            else:
                prot_data.append((pdb_id, mean_b_values))
    except:
        continue

    data_dict[prot_label] = prot_data

plot_stuff(data_dict)


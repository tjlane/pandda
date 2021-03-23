import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys, copy, itertools

import numpy as np

from giant.exceptions import Sorry, Failure

from giant.plot import setup_plotting
setup_plotting()

#################################

PROGRAM = 'giant.score_model'
DESCRIPTION = """
    A tool to quickly score residues (e.g. ligands) against crystallographic electron density.

    Single Structures:

    1) Simple usage (for ligand called LIG or UNL):
        > giant.score_model <filename>.pdb <filename>.mtz

    2) Score particular residues in the file (residues XX1, XX2 and XX3)
        > giant.score_model ... resnames=XX1,XX2,XX3

    3) Define a "reference" model and data to compare the input model to. Providing the MTZ is optional
        > giant.score_model ... ref_pdb=<filename>.pdb ref_mtz=<filename>.mtz
        > giant.score_model ... ref_pdb=<filename>.pdb

    Multiple Structures:

    1) Simple usage:
        Mtz files can be provided explicitly:
        > giant.score_model filename1.pdb filename1.mtz filename2.pdb filename2.mtz filename3.pdb filename3.mtz
        or can be omitted when the same same as the pdb files (e.g. filename1.mtz, filename2.mtz, etc)
        > giant.score_model filename1.pdb filename2.pdb filename3.pdb

    See giant.score_model for more information.
"""

#######################################

blank_arg_prepend = {
    '.pdb' : 'pdb=',
    '.mtz' : 'mtz=',
    None : 'directory=',
}

residue_plot_phil =  """
plot {
    mode = *separate combined
        .help = "Combine residues into one plot or plot as separate images?"
        .type = choice(multi=True)
    limits = manual *automatic
        .help = "Rescale plot limits to the minimum/maximum values (automatic) or use preset values."
        .type = choice(multi=False)
    remove_blank_entries = False
        .type = bool
    axis_labels = none *simplified all
        .type = choice(multi=False)
    parameters {
        rscc {
            title = '\n\nModel\nQuality\n(RSCC)'
                .type = str
            axis_min = 0.60
                .type = float
            axis_max = 0.85
                .type = float
            axis_invert = True
                .type = bool
        }
        rszd {
            title = 'Model\nAccuracy\n(RSZD)'
                .type = str
            axis_min = 1.50
                .type = float
            axis_max = 4.00
                .type = float
            axis_invert = False
                .type = bool
        }
        rszo {
            title = 'Density\nPrecision\n(RSZO/OCC)'
                .type = str
            axis_min = 0.00
                .type = float
            axis_max = 2.00
                .type = float
            axis_invert = True
                .type = bool
        }
        b_factor_ratio {
            title = 'B-Factor\nStability\n(B-factor Ratio)'
                .type = str
            axis_min = 1.00
                .type = float
            axis_max = 3.00
                .type = float
            axis_invert = False
                .type = bool
        }
        rmsd {
            title = 'Coordinate\nStability\n(RMSD)'
                .type = str
            axis_min = 0.00
                .type = float
            axis_max = 1.50
                .type = float
            axis_invert = False
                .type = bool
        }
    }
}
"""

import libtbx.phil
master_phil = libtbx.phil.parse("""
input {

    pdb = None
        .type = str
        .multiple = True
    mtz = None
        .type = str
        .multiple = True

    ref_pdb = None
        .type = str
        .multiple = True
    ref_mtz = None
        .type = str
        .multiple = True

    label = None
        .help = "what is the label for the input structure -- only valid for scoring of a single structure"
        .type = str
        .multiple = False
    label_func = *filename foldername
        .help = "labelling function for each input structure -- only valid for scoring multiple structures"
        .type = choice(multi=False)

    f_label = None
        .help = 'Column label for experimental amplitudes in input mtz files. If left blank will try to guess the labels.'
        .type = str
        .multiple = False

    directory = None
        .type = str
        .multiple = True
    pdb_style = None
        .type = str
    mtz_style = None
        .type = str
    ref_pdb_style = None
        .type = str
    ref_mtz_style = None
        .type = str

}
options {
    resnames = LIG,UNL,DRG
        .type = str
        .help = "Comma-separated list of residue names to score"
}
output {
    out_dir = ./
        .type = path
}
include scope giant.phil.settings_phil
"""+residue_plot_phil, process_includes=True)

#######################################

table_columns = [
    'Residue',
    'Model RMSD',
    'RSCC',
    'RSZD',
    'RSZO',
    'RSR',
    'RSZO/OCC',
    'Occupancy',
    'dRSCC',
    'dRSZD',
    'dRSZO',
    'dRSR',
    'dRSZO/OCC',
    'dOccupancy',
    'Surroundings B-factor Ratio',
    'Surroundings Residue Labels',
    'Average B-factor (Residue)',
    'Average B-factor (Surroundings)',
    'PDB',
    'MTZ',
    'PDB-2',
    'MTZ-2',
]

def prepare_output_directory(params):
    from giant.paths import easy_directory
    out_dir = easy_directory(params.output.out_dir)
    img_dir = easy_directory(os.path.join(out_dir, 'residue_plots'))
    return out_dir, img_dir

def score_model_wrapper(kw_args):
    try:
        return score_model(**kw_args)
    except Exception as e:
        return e

def score_model(params, pdb1, mtz1, pdb2=None, mtz2=None, label_prefix=''):
    """
    Score residues against density, and generate other model quality indicators.
    Identified residues in pdb1 are scored against mtz1 (and mtz2, if provided) using edstats.
    Identified residues in pdb1 are compared to the equivalent residues in pdb2, if provided.
    B-factors ratios of identified residues to surrounding sidechains are calculated.
    """

    from giant.io.pdb import strip_pdb_to_input
    from giant.xray.edstats import EdstatsFactory
    from giant.structure import sanitise_hierarchy, calculate_paired_conformer_rmsds
    from giant.structure.b_factors import calculate_residue_group_bfactor_ratio
    from giant.structure.occupancy import calculate_residue_group_occupancy
    from giant.structure.formatting import ShortLabeller
    from giant.structure.select import non_h, protein, backbone, sidechains

    if (label_prefix) and (not label_prefix.endswith('-')):
        label_prefix = label_prefix + '-'

    # Extract variables from params
    res_names = params.options.resnames.split(',')
    f_label = params.input.f_label

    logger('Reading main input structure: {}'.format(pdb1))

    # Extract Structure and strip hydrogens etc
    h1_all = non_h(
        strip_pdb_to_input(
            pdb1,
            remove_ter = True,
            remove_end = True,
        ).hierarchy
    )

    # Normalise hierarchy (standardise atomic naming, etc...)
    sanitise_hierarchy(h1_all)

    # Extract subselections
    h1_pro = protein(h1_all)
    h1_bck = backbone(h1_all)
    h1_sch = sidechains(h1_all)

    # Pull out residues to analyse
    if res_names:
        rg_for_analysis = [
            rg for rg in h1_all.residue_groups()
            if [n for n in rg.unique_resnames() if n in res_names]
        ]
        logger(
            'Selecting residues named {}: {} residue(s)'.format(
                ' or '.join(res_names),
                len(rg_for_analysis),
            )
        )
    else:
        rg_for_analysis = h1_all.residue_groups()
        logger(
            'Analysing all residues ({} residues)'.format(
                len(rg_for_analysis),
            )
        )

    # Check residues to analyse or skip
    if (not rg_for_analysis):
        raise Sorry('There are no residues called {} in {}'.format(' or '.join(res_names), pdb1))

    # Extract PDB2
    if (pdb2 is not None):
        logger('Reading secondary (comparison) input structure: {}'.format(pdb2))
        h2_all = non_h(
            strip_pdb_to_input(
                pdb2,
                remove_ter = True,
                remove_end = True,
            ).hierarchy
        )
        sanitise_hierarchy(h2_all)

    # Score structures against MTZ(s)
    logger('Scoring model(s) against mtz file(s) with EDSTATs')
    run_edstats = EdstatsFactory()
    mtz1_edstats_scores = mtz2_edstats_scores = None

    # Score MTZ1
    if (mtz1 is not None):
        pdb = pdb1
        logger('Scoring {} >>> {}'.format(pdb, mtz1))
        mtz1_edstats_scores = run_edstats(
            pdb_file = pdb,
            mtz_file = mtz1,
            f_label = f_label,
        )

    # Score MTZ2
    if (mtz2 is not None):
        pdb = (pdb2 if (pdb2 is not None) else pdb1)
        logger('Scoring {} >>> {}'.format(pdb, mtz2))
        mtz2_edstats_scores = run_edstats(
            pdb_file = pdb,
            mtz_file = mtz2,
            f_label = f_label,
        )

    # Prepare output table
    import pandas
    data_table = pandas.DataFrame(columns=table_columns)

    for rg_sel in rg_for_analysis:

        # label for the current residue
        rg_label = ShortLabeller.format(rg_sel)

        # Create label for the output table
        full_label = (label_prefix + rg_label).replace(' ','')

        if len(rg_sel.unique_resnames()) != 1:
            raise Sorry(rg_label+': More than one residue name associated with residue group in the input structure -- cannot process')

        # Residue occupancy
        rg_occ = calculate_residue_group_occupancy(
            residue_group = rg_sel,
        )

        # B-factor ratios to surroundings
        b_factor_scores = calculate_residue_group_bfactor_ratio(
            residue_group = rg_sel,
            hierarchy = h1_sch,
        )

        # Extract structural comparative values
        if (pdb2 is None):
            # Compare to pdb1 instead
            rg_occ_2 = rg_occ
            model_rmsd = None
        else:

            # Extract the equivalent residue in pdb2
            rg_sel_2 = [
                rg for rg in h2_all.residue_groups()
                if (
                    ShortLabeller.format(rg) == rg_label
                )
            ]

            if (not rg_sel_2):
                raise ValueError('Residue is not present in pdb file: {} not in {}'.format(rg_label, pdb2))
            if len(rg_sel_2) >= 1:
                raise ValueError('More than one residue has been selected for {} in {}'.format(rg_label, pdb2))
            rg_sel_2 = rg_sel_2[0]

            # Extract occupancy
            rg_occ_2 = calculate_residue_group_occupancy(
                residue_group = rg_sel_2,
            )

            # Calculate the RMSD between the models
            confs1, confs2, rmsds = zip(
                *calculate_paired_conformer_rmsds(
                    conformers_1 = rg_sel.conformers(),
                    conformers_2 = rg_sel_2.conformers(),
                )
            )
            model_rmsd = min(rmsds)

        # Extract Density Scores - MTZ 1
        ed_dict = {} # initialise optional values to None
        if (mtz1_edstats_scores is not None):
            ed_dict = mtz1_edstats_scores.extract_residue_group_scores(
                residue_group  = rg_sel,
            )
            # Normalise the RSZO by the Occupancy of the ligand
            ed_dict['rszo_norm'] = ed_dict['rszo'] / rg_occ

        # Extract Density Scores - MTZ 2
        ed_dict_2 = {}
        d_ed_dict = {}
        if (mtz2_edstats_scores is not None):
            ed_dict_2 = mtz2_edstats_scores.extract_residue_group_scores(
                residue_group  = rg_sel,
            )
            # Normalise the RSZO by the Occupancy of the ligand
            ed_dict_2['rszo_norm'] = ed_dict_2['rszo'] / rg_occ_2
            # Delta scores
            d_ed_dict.update({
                'dRSCC' : ed_dict.get('rscc') - ed_dict_2.get('rscc'),
                'dRSZD' : ed_dict.get('rszd') - ed_dict_2.get('rszd'),
                'dRSZO' : ed_dict.get('rszo') - ed_dict_2.get('rszo'),
                'dRSR'  : ed_dict.get('rsr') - ed_dict_2.get('rsr'),
                'dRSZO/OCC' : ed_dict.get('rszo_norm') - ed_dict_2.get('rszo_norm'),
            })

        rg_dict = {
            'Residue' : rg_label,
            'Model RMSD' : model_rmsd,
            'Occupancy' : rg_occ,
            'dOccupancy' : rg_occ - rg_occ_2,
            'RSCC' : ed_dict.get('rscc'),
            'RSZD' : ed_dict.get('rszd'),
            'RSZO' : ed_dict.get('rszo'),
            'RSR'  : ed_dict.get('rsr'),
            'RSZO/OCC' : ed_dict.get('rszo_norm'),
            'dRSCC' : d_ed_dict.get('dRSCC'),
            'dRSZD' : d_ed_dict.get('dRSZD'),
            'dRSZO' : d_ed_dict.get('dRSZO'),
            'dRSR'  : d_ed_dict.get('dRSR'),
            'dRSZO/OCC' : d_ed_dict.get('dRSZO/OCC'),
            'Surroundings B-factor Ratio' : b_factor_scores.b_factor_ratio,
            'Surroundings Residue Labels' : b_factor_scores.surroundings_names,
            'Average B-factor (Residue)' : b_factor_scores.selection_av_bfactor,
            'Average B-factor (Surroundings)' : b_factor_scores.surroundings_av_bfactor,
            'PDB' : pdb1,
            'MTZ' : mtz1,
            'PDB-2' : pdb2,
            'MTZ-2' : mtz2,
        }

        # Append to table
        data_table.loc[full_label] = rg_dict

    return data_table

def label_and_value(v):
    try:
        l = round(v, 2)
    except Exception as e:
        l = 'n/a'
        v = None
    return l, v

def make_residue_radar_plot(path, data, columns_dict, linetype=None, remove_blank_entries=False, axis_labels='simplified'):
    "Plot radar graph. data is a list of (label, scores) tuples. label is a string. scores is a pandas Series."

#    column_limit = [(0.6,0.85),(1.5,4),(0,2),(1,3),(0,1.5)]
#    column_names = ['RSCC','RSZD','RSZO/OCC','Surroundings B-factor Ratio','Model RMSD']
#    column_title = ['Model\nQuality\n(RSCC)', 'Model\nAccuracy\n(RSZD)', 'Model\nPrecision\n(RSZO/OCC)', 'B-Factor\nRatio', 'Model\nRMSD']
#    column_invse = [1,0,1,0,0]

    from giant.plot import Radar

    assert isinstance(columns_dict, dict)

    column_title = columns_dict.get('titles') # titles of the axes
    column_names = columns_dict.get('names')  # keys for data frame
    column_limit = columns_dict.get('limits') # axis limits
    column_invse = columns_dict.get('invert') # which axes are inverted

    assert (column_title is not None)
    assert (column_names is not None)
    assert len(column_names) == len(column_title)
    if (column_limit is not None):
        assert len(column_limit) == len(column_title)
    if (column_invse is not None):
        assert len(column_invse) == len(column_title)

    # Extract the plot data from the data_frame
    plot_data = data[column_names]
    n_lines = len(data.index)

    if (remove_blank_entries is True):

        # Filter the entries based on whether there is at least one values for each column
        data_mask_i = [i for i, c in enumerate(column_names) if data[c].any()]
        filter_func = lambda values: [values[i] for i in data_mask_i]

        # Filter against the mask
        column_title = filter_func(column_title)
        column_names = filter_func(column_names)
        if column_invse:
            column_invse = filter_func(column_invse)
        if column_limit:
            column_limit = filter_func(column_limit)

        # Reselect the plot_data
        plot_data = data[column_names]

    # Format tick values for plotting
    col_tick_vals = [[] for c in column_names] # blank
    col_tick_labs = [[] for c in column_names] # blank

    if (axis_labels == 'all'):
        for i_col, col in enumerate(column_names):
            for v in plot_data[col]:
                l, v = label_and_value(v)
                col_tick_vals[i_col].append(v)
                col_tick_labs[i_col].append(l)
    elif (axis_labels == 'simplified'):
        for i_col, col in enumerate(column_names):
            vals = plot_data[col]
            if len(vals) == 1:
                l, v = label_and_value(vals[0])
                col_tick_vals[i_col].append(v)
                col_tick_labs[i_col].append(l)
    else:
        pass # leave blank

    # Colors of lines
    colors = ['r', 'g', 'b', 'y']

    # Prepare markers
    if (plot_data.index.size > 1):
        markers = ['o','^','s','D','*','+'][0:len(colors)]
        markersize = 5
    else:
        markers = [None]
        markersize = None

    # Prepare linetypes
    if (linetype is None):
        linetype = np.repeat(['-','--'], len(markers)).tolist()
    elif isinstance(linetype, str):
        linetype = [linetype]
    else:
        pass

    # Make cycles for simplicity
    colors_iter = itertools.cycle(colors)
    linetype_iter = itertools.cycle(linetype)
    markers_iter = itertools.cycle(markers)

    # Add each row of the data frame as a separate line
    r = Radar(titles=column_title)
    for label, row in plot_data.iterrows():
        r.add(
            row[column_names].tolist(),
            color = colors_iter.next(),
            linestyle = linetype_iter.next(),
            marker = markers_iter.next(),
            markersize = markersize,
            markeredgecolor = 'k',
            label = label,
            lw = 2,
        )
    # Set axes to invert
    if column_invse:
        r.set_inversion(column_invse)
    # Set manual limits
    if column_limit:
        r.set_limits(column_limit)
    # Axis values to plot
    r.set_ticks(values=col_tick_vals, labels=col_tick_labs)

    # Plot, modify and save
    r.plot()
    r.legend()
    r.axis_limits()
    r.savefig(path)
    r.close()

    return

def format_parameters_for_plot(params):
    """Convert plot scope parameters to parameter dict"""

    p = params
    columns = {}
    columns['titles'] = [
        p.rscc.title,
        p.rszd.title,
        p.rszo.title,
        p.b_factor_ratio.title,
        p.rmsd.title,
    ]
    columns['names']  = [
        'RSCC',
        'RSZD',
        'RSZO/OCC',
        'Surroundings B-factor Ratio',
        'Model RMSD',
    ]
    columns['invert'] = [
        p.rscc.axis_invert,
        p.rszd.axis_invert,
        p.rszo.axis_invert,
        p.b_factor_ratio.axis_invert,
        p.rmsd.axis_invert,
    ]
    columns['limits'] = [
        (p.rscc.axis_min,           p.rscc.axis_max),
        (p.rszd.axis_min,           p.rszd.axis_max),
        (p.rszo.axis_min,           p.rszo.axis_max),
        (p.b_factor_ratio.axis_min, p.b_factor_ratio.axis_max),
        (p.rmsd.axis_min,           p.rmsd.axis_max),
    ]

    return columns

#######################################

def replace_extension(filename, current_extension, new_extension):
    assert current_extension.startswith('.')
    assert new_extension.startswith('.')
    assert filename.endswith(current_extension)
    base, ext = os.path.splitext(filename)
    assert ext == current_extension
    return (base + new_extension)

def process_and_validate_input_files(params):

    if (params.input.pdb and params.input.directory):
        raise Sorry('Must provide input structures OR input directories')

    if (not params.input.pdb) and (not params.input.directory):
        raise Sorry('No input structures provided')

    if (params.input.pdb):

        pdb_files = params.input.pdb
        mtz_files = params.input.mtz

        if (not mtz_files):

            mtz_files = [replace_extension(f, '.pdb', '.mtz') for f in pdb_files]

            logger(
                'No MTZ files provided. Looking for files with the same naming as pdb files: \n\t{}'.format(
                    '\n\t'.join(['{} -> {}'.format(p,m) for p,m in zip(pdb_files, mtz_files)]),
                )
            )

    elif (params.input.directory):

        from giant.paths import resolve_glob

        pdb_style = params.input.pdb_style
        mtz_style = params.input.mtz_style

        if (pdb_style is None):
            raise Sorry('Must provide pdb_style when providing directories as input')
        if (mtz_style is None):
            raise Sorry('Must provide mtz_style when providing directories as input')

        pdb_files = []
        mtz_files = []

        logger(
            'Looking for pdb & mtz files with names "{}" and "{}" in directories: \n\t{}'.format(
                pdb_style,
                mtz_style,
                '\n\t'.join(params.input.directory),
            )
        )

        for in_d in params.input.directory:

            try:
                pdb = resolve_glob(os.path.join(in_d, pdb_style), n=1)
                mtz = resolve_glob(os.path.join(in_d, mtz_style), n=1)
            except:
                logger('Failed to find the expected files in directory "{}"'.format(in_d))
                raise

            pdb_files.append(pdb)
            mtz_files.append(mtz)

    # Check number of files
    n_pdb = len(pdb_files)
    n_mtz = len(mtz_files)

    if (n_mtz != n_pdb):
        raise Sorry(
            'Different numbers of PDB and MTZ files have been provided ({} pdbs != {} mtzs)'.format(
                n_pdb, n_mtz,
            )
        )

    for f in (pdb_files + mtz_files):

        if not os.path.exists(f):
            raise IOError(
                'Input file does not exist ({}). Please check provided command line arguments'.format(
                    f,
                )
            )

    # Labels
    if (n_pdb == 1) and (params.input.label is not None):
        labels = [params.input.label]
    else:
        from giant.paths import filename, foldername
        if params.input.label_func == 'filename':
            label_func = filename
        else:
            label_func = foldername
        labels = [label_func(f) for f in pdb_files]

        if len(labels) != len(set(labels)):
            raise ValueError(
                'Labels generated using "{}" lead to duplicate labels.\n\t{}.\nChange labelling function or rename input files.'.format(
                    params.input.label_func,
                    '\n\t'.join(sorted(labels)),
                )
            )

    return pdb_files, mtz_files, labels

def process_reference_files(params, n):

    if (params.input.pdb):

        ref_pdbs = params.input.ref_pdb
        ref_mtzs = params.input.ref_mtz

        if (not ref_pdbs):
            ref_pdbs = [None] * n
        elif len(ref_pdbs) == 1:
            ref_pdbs = list(ref_pdbs) * n

        if (not ref_mtzs):
            ref_mtzs = [None] * n
        elif len(ref_mtzs) == 1:
            ref_mtzs = list(ref_mtzs) * n

    elif (params.input.directory):

        from giant.paths import resolve_glob

        ref_pdb_style = params.input.ref_pdb_style
        ref_mtz_style = params.input.ref_mtz_style

        if (ref_mtz_style is not None) and (not ref_pdb_style is None):
            raise Sorry('Cannot provide ref_mtz_style without providing ref_pdb_style')

        ref_pdbs = []
        ref_mtzs = []

        logger(
            'Looking for reference pdb & mtz files with names "{}" and "{}" in directories: \n\t{}'.format(
                ref_pdb_style,
                ref_mtz_style,
                '\n\t'.join(params.input.directory),
            )
        )

        for in_d in params.input.directory:

            try:
                if (ref_pdb_style is not None):
                    pdb = resolve_glob(os.path.join(in_d, ref_pdb_style), n=1)
                else:
                    pdb = None
                if (ref_mtz_style is not None):
                    mtz = resolve_glob(os.path.join(in_d, ref_mtz_style), n=1)
                else:
                    mtz = None
            except:
                logger('Failed to find the expected files in directory "{}"'.format(in_d))
                raise

            ref_pdbs.append(pdb)
            ref_mtzs.append(mtz)

    if len(ref_pdbs) != n:
        raise IOError(
            '{} reference pdbs have been provided -- expected {}'.format(
                len(ref_pdbs),
                n,
            )
        )

    if len(ref_mtzs) != n:
        raise IOError(
            '{} reference mtzs have been provided -- expected {}'.format(
                len(ref_mtzs),
                n,
            )
        )

    for f in (ref_pdbs + ref_mtzs):

        if (f is None):
            continue

        if not os.path.exists(f):
            raise IOError(
                'Input file does not exist ({}). Please check provided command line arguments'.format(
                    f,
                )
            )

    return ref_pdbs, ref_mtzs

def run(params):

    # REMOVE THIS WHEN IMPLEMENTED
    if not params.options.resnames:
        raise ValueError('Must provide resnames')

    logger.heading('Processing input files')
    pdb_files, mtz_files, labels = process_and_validate_input_files(params)
    ref_pdbs, ref_mtzs = process_reference_files(params, n=len(pdb_files))

    output_dir, images_dir = prepare_output_directory(params)
    scores_file = os.path.join(output_dir, 'residue_scores.csv')
    summary_file = os.path.join(output_dir, 'residue_scores.html')

    # Iterate through input files and generate sets of arguments
    input_args = []
    for pdb, mtz, lab, ref_pdb, ref_mtz in zip(
        pdb_files, mtz_files, labels, ref_pdbs, ref_mtzs,
        ):

        input_args.append(
            dict(
                params = params,
                pdb1 = pdb,
                mtz1 = mtz,
                pdb2 = ref_pdb,
                mtz2 = ref_mtz,
                label_prefix = lab,
            )
        )

        #try:

        #    # Generate scores as table
        #    logger.heading('Scoring model')
        #    data_table = score_model(
        #        params = params,
        #        pdb1 = pdb,
        #        mtz1 = mtz,
        #        pdb2 = ref_pdb,
        #        mtz2 = ref_mtz,
        #        label_prefix = lab,
        #    )

        #    # Append to output list
        #    output_tables.append(data_table)
        #    logger.subheading('done!')

        #except Sorry as e:
        #    logger.heading('Program exited with an error')
        #    logger(str(e))
        #    sys.exit(1)

        #except Failure as e:
        #    raise

    logger.heading('Scoring models')

    if params.settings.cpus == 1:
        output_tables = []
        for kw_args in input_args:
            output_tables.append(score_model(**kw_args))
    else:
        import multiprocessing
        workers = multiprocessing.Pool(processes=params.settings.cpus)
        output_tables = workers.map(score_model_wrapper, input_args)
        workers.close()

    import pandas
    combined_table = pandas.concat(output_tables)

    logger.heading('Writing output csv')

    combined_table.dropna(axis=1, how='all').to_csv(scores_file)
    logger('Output written to {}'.format(scores_file))

    logger.heading('Generating output plots')

    # Create dictionary of plot parameters
    columns_dict_master = format_parameters_for_plot(params=params.plot.parameters)

    # Create ouptut lists to be optionally populated
    import collections
    separate_images = collections.OrderedDict()
    combined_images = collections.OrderedDict()

    if 'separate' in params.plot.mode:

        logger.subheading('Generating separate plots')

        # Create copy of dict to manipulate
        columns_dict = copy.deepcopy(columns_dict_master)

        # Input data
        graph_table_data = combined_table

        # Iterate through combined table and generate one plot for each row
        for label, row in graph_table_data.iterrows():

            image_path = os.path.join(
                images_dir,
                '{}.png'.format(
                    label.replace(' ','')
                )
            )

            logger('{} >> {}'.format(label, image_path))

            make_residue_radar_plot(
                path = image_path,
                data = row.to_frame().T,
                columns_dict = columns_dict,
                remove_blank_entries = params.plot.remove_blank_entries,
                axis_labels = params.plot.axis_labels,
            )

            separate_images[label] = image_path

    if ('combined' in params.plot.mode) and (len(pdb_files) > 1):

        logger.subheading('Generating comparative plots')

        # Create copy of dict to manipulate
        columns_dict = copy.deepcopy(columns_dict_master)

        # Remove limits information if requested
        if params.plot.limits == 'automatic':
            columns_dict.pop('limits', None)
        elif params.plot.limits == 'manual':
            pass

        # Input data
        graph_table_data = combined_table

        # Sort into groups based on residue labels
        rg_labels = graph_table_data['Residue']

        from giant.stats.cluster import generate_group_idxs
        for common_label, idxs in generate_group_idxs(rg_labels):

            # Create output label
            out_label = 'residue-' + common_label.replace(' ','')

            image_path = os.path.join(
                images_dir,
                'compare-{}.png'.format(
                    out_label.replace(' ','')
                ),
            )

            logger('{} >> {}'.format(out_label, image_path))

            make_residue_radar_plot(
                path = image_path,
                data = graph_table_data.iloc[idxs],
                columns_dict = columns_dict,
                remove_blank_entries = params.plot.remove_blank_entries,
                axis_labels = params.plot.axis_labels,
            )

    if True:

        logger.subheading('Collating all images and outputting to output HTML...')

        # Get template to be filled in
        from giant.html import HTML_ENV
        template = HTML_ENV.get_template('summary_page.html')

        # Output directory (for relative symlinks)
        out_dir = os.path.abspath(os.path.dirname(summary_file))

        # ===========================================================>
        # Construct the data object to populate the template
        output_data = {}
        output_data['header'] = 'Residue Score Summaries'
        output_data['title'] = 'Residue Score Summaries'
        output_data['introduction'] = 'Model Quality and Validation checks.'
        # ===========================================================>
        # Header Images
        output_data['small_images'] = []
        for img in (combined_images.values() + separate_images.values()):
            output_data['small_images'].append({
                'path': './'+os.path.relpath(path=img, start=out_dir),
                'title': 'Scores for {}'.format(os.path.splitext(os.path.basename(img))[0]),
            })
        # ===========================================================>
        # Write Output
        with open(summary_file, 'w') as out_html:
            out_html.write(template.render(output_data))

    logger.heading('finished normally!')

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
        program             = PROGRAM,
        description         = DESCRIPTION,
    )

import os, sys, copy, re, glob

#################################
import matplotlib
matplotlib.use('Agg')
matplotlib.interactive(0)
from matplotlib import pyplot
pyplot.style.use('ggplot')
#################################

import libtbx.phil

import numpy, pandas

from bamboo.plot import Radar
from bamboo.edstats import Edstats

from giant.utils.pdb import strip_pdb_to_input
from giant.structure import calculate_residue_group_occupancy, calculate_residue_group_rmsd, calculate_residue_group_bfactor_ratio
from giant.xray.edstats import extract_residue_group_density_scores

#######################################

bar = '=======================++>'

#######################################

blank_arg_prepend = {'.pdb':'pdb1=', '.mtz':'mtz1='}

master_phil = libtbx.phil.parse("""
input {
    pdb1 = None
        .type = str
        .multiple = False
    mtz1 = None
        .type = str
        .multiple = False

    pdb2 = None
        .type = str
        .multiple = False
    mtz2 = None
        .type = str
        .multiple = False

    label = ''
        .type = str
        .multiple = False
}
selection {
    res_names = LIG,UNL,DRG
        .type = str
        .help = "Comma-separated list of residue names to score -- if None then scores all residues"
}
output {
    out_dir = ./
        .type = path
}
""")

#######################################

def sanitise_hierarchy(hierarchy):
    hierarchy.atoms().set_chemical_element_simple_if_necessary()

def prepare_table():

    columns_fix = ['Model RMSD']
    columns_num = ['RSCC','RSZD','RSZO','OCC','RSR',
                    'Average B-factor (Residue)',
                    'Average B-factor (Surroundings)',
                    'Surroundings B-factor Ratio']
    columns_end = ['PDB','MTZ']

    columns = []
    columns.extend(columns_fix)
    columns.extend(columns_num)
    columns.extend([c+'-2' for c in columns_num])
    columns.extend(columns_end)
    columns.extend([c+'-2' for c in columns_end])

    return pandas.DataFrame(index=[], columns=columns)

def prepare_output_directory(params):
    if not os.path.exists(params.output.out_dir): os.mkdir(params.output.out_dir)
    images_dir =   os.path.join(params.output.out_dir, 'residue_plots')
    if not os.path.exists(images_dir): os.mkdir(images_dir)
    return params.output.out_dir, images_dir

def score_model(params, pdb1, mtz1, pdb2=None, mtz2=None, label_prefix=''):
    """
    Score residues against density, and generate other model quality indicators.
    Identified residues in pdb1 are scored against mtz1 (and mtz2, if provided) using edstats.
    Identified residues in pdb1 are compared to the equivalent residues in pdb2, if provided.
    B-factors ratios of identified residues to surrounding sidechains are calculated.
    """

    if label_prefix: label_prefix = label_prefix + '-'

    # Extract the residues to look for
    res_names = params.selection.res_names_list

    # Extract Structure
    h1_all = strip_pdb_to_input(pdb1, remove_ter=True, remove_end=True).hierarchy
    sanitise_hierarchy(h1_all)
    h1_all = h1_all.select(h1_all.atom_selection_cache().selection('not element H'), copy_atoms=True)
    h1_pro = h1_all.select(h1_all.atom_selection_cache().selection('pepnames'), copy_atoms=True)
    h1_bck = h1_pro.select(h1_pro.atom_selection_cache().selection('(name CA or name C or name O or name N)'), copy_atoms=True)
    h1_sch = h1_pro.select(h1_pro.atom_selection_cache().selection('not (name CA or name C or name O or name N)'), copy_atoms=True)

    # Pull out residues to analyse
    if res_names: rg_for_analysis = [rg for rg in h1_all.residue_groups() if [n for n in rg.unique_resnames() if n in res_names]]
    else:         rg_for_analysis = h1_all.residue_groups()

    # Check residues to analyse or skip
    if not rg_for_analysis:
        raise Exception('No residues to analyse: {}'.format(pdb1))

    # Extract PDB2
    if pdb2 is not None:
        h2_all = strip_pdb_to_input(pdb2, remove_ter=True, remove_end=True).hierarchy
        sanitise_hierarchy(h2_all)
        h2_all = h2_all.select(h2_all.atom_selection_cache().selection('not element H'), copy_atoms=True)

    # Score MTZ1
    if mtz1 is not None: mtz1_edstats_scores = Edstats(mtz_file=mtz1, pdb_file=pdb1)
    else:                mtz1_edstats_scores = None
    # Score MTZ2
    if mtz2 is not None: mtz2_edstats_scores = Edstats(mtz_file=mtz2, pdb_file=pdb1)
    else:                mtz2_edstats_scores = None

    # Prepare output table
    data_table = prepare_table()

    for rg_sel in rg_for_analysis:

        # Create label for the output table
        rg_label = (label_prefix+rg_sel.unique_resnames()[0]+'-'+rg_sel.parent().id+'-'+rg_sel.resseq+rg_sel.icode).replace(' ','')

        if len(rg_sel.unique_resnames()) != 1:
            raise Exception(rg_label+': More than one residue name associated with residue group -- cannot process')

        # Append empty row to output table
        data_table = data_table.append(pandas.DataFrame(index=[rg_label]), verify_integrity=True)

        data_table.set_value(   index = rg_label, col='PDB', value=pdb1 )
        data_table.set_value(   index = rg_label,
                                col   = 'Occupancy',
                                value = calculate_residue_group_occupancy(residue_group=rg_sel) )

        data_table = calculate_residue_group_bfactor_ratio( residue_group   = rg_sel,
                                                            hierarchy       = h1_sch,
                                                            data_table      = data_table,
                                                            rg_label        = rg_label
                                                        )

        if pdb2 is not None:
            data_table.set_value(   index = rg_label, col='PDB-2', value=pdb2 )

            # Extract the equivalent residue in pdb2
            rg_sel_2 = [rg for rg in h2_all.residue_groups() if (rg.unique_resnames() == rg_sel.unique_resnames())
                                                            and (rg.parent().id       == rg_sel.parent().id)
                                                            and (rg.resseq            == rg_sel.resseq)
                                                            and (rg.icode             == rg_sel.icode)              ]
            assert rg_sel_2, 'Residue is not present in pdb file: {} not in {}'.format(rg_label, pdb2)
            assert len(rg_sel_2) == 1, 'More than one residue has been selected for {} in {}'.format(rg_label, pdb2)

            # Calculate the RMSD between the models
            data_table.set_value(   index = rg_label,
                                    col   = 'Model RMSD',
                                    value = calculate_residue_group_rmsd(residue_group_1=rg_sel, residue_group_2=rg_sel_2[0])   )

        # Extract Density Scores - MTZ 1
        if mtz1 is not None:
            data_table.set_value(   index = rg_label, col='MTZ', value=mtz1)
        if mtz1_edstats_scores is not None:
            data_table = extract_residue_group_density_scores(  residue_group  = rg_sel,
                                                                edstats_scores = mtz1_edstats_scores,
                                                                data_table     = data_table,
                                                                rg_label       = rg_label )
            # Normalise the RSZO by the Occupancy of the ligand
            data_table['RSZO/OCC'] = data_table['RSZO']/data_table['Occupancy']

        # Extract Density Scores - MTZ 2
        if mtz2 is not None:
            data_table.set_value(   index = rg_label, col='MTZ-2', value=mtz2)
        if mtz2_edstats_scores is not None:
            data_table = extract_residue_group_density_scores(  residue_group  = rg_sel,
                                                                edstats_scores = mtz2_edstats_scores,
                                                                data_table     = data_table,
                                                                rg_label       = rg_label,
                                                                column_suffix  = '-2' )
            # Normalise the RSZO by the Occupancy of the ligand
            data_table['RSZO/OCC-2'] = data_table['RSZO-2']/data_table['Occupancy']

    return data_table

def make_residue_radar_plot(path, data, columns=None):
    "Plot radar graph. data is a list of (label, scores) tuples. label is a string. scores is a pandas Series."

    column_limit = [(0.6,0.85),(1.5,4),(0,2),(1,3),(0,1.5)]
    column_names = ['RSCC','RSZD','RSZO/OCC','Surroundings B-factor Ratio','Model RMSD']
    column_title = ['Model Quality (RSCC)', 'Model Accuracy (RSZD)', 'Model Precision (RSZO/OCC)', 'B-Factor Ratio', 'Model RMSD']
    column_invse = [1,0,1,0,0]

    if columns:
        assert isinstance(columns, dict)
        if 'limits' in columns: column_limit = columns['limits']
        if 'names'  in columns: column_names = columns['names']
        if 'titles' in columns: column_title = columns['titles']
        if 'invert' in columns: column_invse = columns['invert']

    # Round column values
    col_tick_vals = []
    col_tick_labs = []
    for col in column_names:
        tvs=[]; tls=[]
        for i, v in enumerate(data[col]):
            try: l = round(v, 2)
            except:
                l = 'Error'; v = None;
            # Add point marker
            tvs.append(v); tls.append(l)
        # Add to column ticks
        col_tick_vals.append(tvs); col_tick_labs.append(tls)

    r = Radar(titles=column_title)
    # Add each row of the data frame as a separate line
    for label, row in data.iterrows():
        r.add(row[column_names].tolist(), "-", label=label, lw=2)
    # Set axis meta manually
    if column_invse: r.set_inversion(column_invse)
    if column_limit: r.set_limits(column_limit)
    r.set_ticks(values=col_tick_vals, labels=col_tick_labs)
    # Plot, modify and save
    r.plot()
    r.ax.legend(loc='upper center', fancybox=True, bbox_to_anchor=(0.5, 0.0))
    r.savefig(path)
    r.close()

    return

#######################################

def run(params):

    assert params.input.pdb1, 'No pdb1 provided'
    assert params.input.mtz1, 'No mtz1 provided'
    # REMOVE THIS WHEN IMPLEMENTED
    assert params.selection.res_names
    # REMOVE THIS WHEN IMPLEMENTED
    if params.selection.res_names:      params.selection.__inject__("res_names_list", params.selection.res_names.split(','))
    else:                               params.selection.__inject__("res_names_list", None)

    output_dir, images_dir = prepare_output_directory(params)
    scores_file = os.path.join(output_dir, 'residue_scores.csv')

    print bar
    print 'Scoring model...'
    data_table = score_model(   params = params,
                                pdb1   = params.input.pdb1,
                                mtz1   = params.input.mtz1,
                                pdb2   = params.input.pdb2,
                                mtz2   = params.input.mtz2,
                                label_prefix = params.input.label
                            )
    print '...Done'
    print bar

    data_table.to_csv(scores_file)
    print 'Output written to {}'.format(scores_file)
    print bar

    # Output Images
    all_images = []
    print 'Generating Output Images...'
    for label, row in data_table.iterrows():
        print 'Making: {}...'.format(label)
        image_path = os.path.join(images_dir,'{}.png'.format(label))
        make_residue_radar_plot(path=image_path, data=row.to_frame().T)
        all_images.append(image_path)
    print '...Done.'
    print bar

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)

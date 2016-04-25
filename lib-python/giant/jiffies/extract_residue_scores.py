import os, sys, copy, re, glob

import libtbx.phil
import libtbx.easy_mp

import iotbx.pdb
import numpy, pandas

from scitbx.array_family import flex

from bamboo.plot import Radar
from bamboo.html import BAMBOO_HTML_ENV
from bamboo.edstats import edstats

from giant.utils.pdb import strip_pdb_to_input
from giant.maths.geometry import is_within

#################################
import matplotlib
matplotlib.interactive(0)
from matplotlib import pyplot
pyplot.style.use('ggplot')
#################################

blank_arg_prepend = {'.pdb':'pdb=', '.mtz':'mtz=', None:'dir='}

master_phil = libtbx.phil.parse("""
input {
    dir = None
        .type = path
        .multiple = True

    pdb_style = 'refine.pdb'
        .type = str
    mtz_style = 'refine.mtz'
        .type = str

    pdb_2_style = '*-ensemble-model.pdb'
        .type = str
    mtz_2_style = '*-input.mtz'
        .type = str

    labels = basename *folder_name
        .type = choice
        .multiple = False
}
selection {
    res_names = LIG,UNL,DRG
        .type = str
        .help = "Comma-separated list of residue names to score -- if None then scores all"
}
output {
    out_dir = residue_scores
        .type = path
    generate_summary = True
        .type = bool
}
settings{
    cpus = 1
        .type = int
}
""")

def sanitise_hierarchy(hierarchy):
    hierarchy.atoms().set_chemical_element_simple_if_necessary()

def prepare_table():

    columns_fix = ['RMSD']
    columns_num = ['RSCC','RSZD','RSZO','OCC','RSR','B_RATIO','B_AV_QUERY','B_AV_CONTACTS']
    columns_end = ['PDB','MTZ']

    columns = []
    columns.extend(columns_fix)
    columns.extend(columns_num)
    columns.extend([c+'-2' for c in columns_num])
    columns.extend(columns_end)
    columns.extend([c+'-2' for c in columns_end])

    return pandas.DataFrame(index=[], columns=columns)

def process_folder_mp(args):
    try:
        return process_folder(*args)
    except Exception as e:
        raise
        return e

def process_folder(dir, params, data_table=None):
    """Score the files in a folder"""

    if data_table is None: data_table = pandas.DataFrame(index=[], column=[])

    try:    pdb  = glob.glob(os.path.join(dir, params.input.pdb_style))[0]
    except: raise Exception('FAILED TO FIND FILE: {} in {}'.format(params.input.pdb_style, dir))
    try:
        if params.input.pdb_2_style:    pdb2 = glob.glob(os.path.join(dir, params.input.pdb_2_style))[0]
        else:                           pdb2 = None
    except: raise Exception('FAILED TO FIND FILE: {} in {}'.format(params.input.pdb_2_style, dir))
    try:
        if params.input.mtz_style:      mtz1 = glob.glob(os.path.join(dir, params.input.mtz_style))[0]
        else:                           mtz1 = None
    except: raise Exception('FAILED TO FIND FILE: {} in {}'.format(params.input.mtz_style, dir))
    try:
        if params.input.mtz_2_style:    mtz2 = glob.glob(os.path.join(dir, params.input.mtz_2_style))[0]
        else:                           mtz2 = None
    except: raise Exception('FAILED TO FIND FILE: {} in {}'.format(params.input.mtz_2_style, dir))

    if   params.input.labels == 'basename':     lab = os.path.split_ext(os.path.basename(pdb))[0]
    elif params.input.labels == 'folder_name':  lab = os.path.basename(dir)

    # Extract the residues to look for
    res_names = params.selection.res_names_list

    print 'Processing: {}'.format(lab)

    # Extract Structure
    h_all = strip_pdb_to_input(pdb, remove_ter=True, remove_end=True).hierarchy
    sanitise_hierarchy(h_all)
    h_all = h_all.select(h_all.atom_selection_cache().selection('not element H'), copy_atoms=True)
    h_pro = h_all.select(h_all.atom_selection_cache().selection('pepnames'), copy_atoms=True)
    h_bck = h_pro.select(h_pro.atom_selection_cache().selection('(name CA or name C or name O or name N)'), copy_atoms=True)
    h_sch = h_pro.select(h_pro.atom_selection_cache().selection('not (name CA or name C or name O or name N)'), copy_atoms=True)

    # Pull out residues to analyse
    if res_names: rg_for_analysis = [rg for rg in h_all.residue_groups() if [n for n in rg.unique_resnames() if n in res_names]]
    else:         rg_for_analysis = h_all.residue_groups()

    # Check residues to analyse or skip
    if not rg_for_analysis:
        raise Exception('No residues to analyse, skipping: {}'.format(pdb))

    # Extract PDB2
    if pdb2 is not None:
        pdb2_h_all = strip_pdb_to_input(pdb2, remove_ter=True, remove_end=True).hierarchy
        sanitise_hierarchy(pdb2_h_all)
        pdb2_h_all = pdb2_h_all.select(pdb2_h_all.atom_selection_cache().selection('not element H'), copy_atoms=True)

    # Score MTZ1
    if mtz1 is not None: mtz1_edstats_scores = edstats(mtz_file=mtz1, pdb_file=pdb)
    else:                mtz1_edstats_scores = None
    # Score MTZ2
    if mtz2 is not None: mtz2_edstats_scores = edstats(mtz_file=mtz2, pdb_file=pdb)
    else:                mtz2_edstats_scores = None

    for rg_sel in rg_for_analysis:

        # Create label for the output table
        row_label = (lab+'-'+rg_sel.unique_resnames()[0]+'-'+rg_sel.parent().id+'-'+rg_sel.resseq+rg_sel.icode).replace(' ','')

        if len(rg_sel.unique_resnames()) != 1:
            raise Exception(row_label+': More than one residue name associated with residue group -- cannot process')

        # Append empty row to output table
        data_table = data_table.append(pandas.DataFrame(index=[row_label]), verify_integrity=True)
        data_table.set_value(index=row_label, col='PDB', value=pdb)

        # Extract for residue
        rg_atoms   = rg_sel.atoms()
        rg_coords  = rg_atoms.extract_xyz()
        rg_confs   = [c.altloc for c in rg.conformers()]

        # Calculate Occupancies of atoms
        if rg_confs == ['']: res_occ = max(rg.atoms().extract_occ())
        else:                res_occ = sum([max(c.atoms().extract_occ()) for c in rg.conformers() if c.altloc])

        if res_occ > 1: res_occ = None
        data_table.set_value(index=row_label, col='OCC',    value=res_occ )

        # Select nearby atom_groups
        near_ags = [ag for ag in h_sch.atom_groups() if (is_within(4, rg_coords, ag.atoms().extract_xyz()) and (ag not in rg_sel.atom_groups()))]
        if near_ags:
            # Extract atoms from nearby groups
            near_ats = iotbx.pdb.hierarchy.af_shared_atom()
            [near_ats.extend(ag.detached_copy().atoms()) for ag in near_ags]
            # Calculate B-factors of atoms
            res_mean_b = flex.mean_weighted(rg_atoms.extract_b(), rg_atoms.extract_occ())
            data_table.set_value(index=row_label, col='B_AV_QUERY',    value=res_mean_b )
            sch_mean_b = flex.mean_weighted(near_ats.extract_b(), near_ats.extract_occ())
            data_table.set_value(index=row_label, col='B_AV_CONTACTS', value=sch_mean_b )
            data_table.set_value(index=row_label, col='B_RATIO', value=res_mean_b/sch_mean_b )

        # Extract Changes between models
        if pdb2 is not None:
            data_table.set_value(index=row_label, col='PDB-2', value=pdb2)
            # Extract the equivalent residue in pdb2
            rg_sel2 = [rg for rg in pdb2_h_all.residue_groups() if  (rg.unique_resnames() == rg_sel.unique_resnames())
                                                                and (rg.parent().id       == rg_sel.parent().id)
                                                                and (rg.resseq            == rg_sel.resseq)
                                                                and (rg.icode             == rg_sel.icode)              ]
            try:
                assert rg_sel2
                assert len(rg_sel2) == 1
                rg_sel2 = rg_sel2[0]

                ats_1 = rg_sel.atoms()
                ats_2 = rg_sel2.atoms()

                assert ats_1.size()                       == ats_2.size()
                assert list(ats_1.extract_name()   )      == list(ats_2.extract_name()   )
                assert list(ats_1.extract_element())      == list(ats_2.extract_element())
                assert [a.parent().altloc for a in ats_1] == [a.parent().altloc for a in ats_2]

                rmsd = (ats_1.extract_xyz()-ats_2.extract_xyz()).rms_length()
                data_table.set_value(index=row_label, col='RMSD', value=rmsd)

            except Exception as e:
                print row_label+': Failed to calculate RMSD between models'
                data_table.set_value(index=row_label, col='RMSD', value=None)

        # Extract Density Scores - MTZ 1
        if mtz1 is not None:
            data_table.set_value(index=row_label, col='MTZ', value=mtz1)
        if mtz1_edstats_scores is not None:
            ed_scores = mtz1_edstats_scores.scores[(rg_sel.unique_resnames()[0], rg_sel.parent().id, rg_sel.resseq_as_int(), rg_sel.icode)]
            data_table.set_value(index=row_label, col='RSCC', value=ed_scores['CCSa'])
            data_table.set_value(index=row_label, col='RSR',  value=ed_scores['Ra']  )
            data_table.set_value(index=row_label, col='B_AV', value=ed_scores['BAa'] )
            data_table.set_value(index=row_label, col='RSZO', value=ed_scores['ZOa'] )
            data_table.set_value(index=row_label, col='RSZD', value=ed_scores['ZDa'] )

        # Extract Density Scores - MTZ 2
        if mtz2 is not None:
            data_table.set_value(index=row_label, col='MTZ-2', value=mtz2)
        if mtz2_edstats_scores is not None:
            ed_scores = mtz2_edstats_scores.scores[(rg_sel.unique_resnames()[0], rg_sel.parent().id, rg_sel.resseq_as_int(), rg_sel.icode)]
            data_table.set_value(index=row_label, col='RSCC-2', value=ed_scores['CCSa'])
            data_table.set_value(index=row_label, col='RSR-2',  value=ed_scores['Ra']  )
            data_table.set_value(index=row_label, col='B_AV-2', value=ed_scores['BAa'] )
            data_table.set_value(index=row_label, col='RSZO-2', value=ed_scores['ZOa'] )
            data_table.set_value(index=row_label, col='RSZD-2', value=ed_scores['ZDa'] )

    return data_table

def make_image(lab, row, path):
    "Plot radar graph of data"

    cols = row[['RSCC','RSZD','RSZO','B_RATIO','RMSD']].tolist()
    col_names = []
    col_names.append('Quality (RSCC)')
    col_names.append('Accuracy (RSZD)')
    if row['OCC'] is None:
        cols = cols[0:2]+[cols[2]]+cols[3:5]
        col_names.append('Precision (RSZO)')
    else:
        cols = cols[0:2]+[cols[2]/row['OCC']]+cols[3:5]
        col_names.append('Precision (RSZO/OCC)')
    col_names.append('B-Factor Ratio')
    col_names.append('RMSD')

    limits = [(0.6,0.85),(1.5,4),(0,2),(1,3),(0,1.5)]

    # Round column values
    col_tick_vals = []
    col_tick_labs = []
    for i,v in enumerate(cols):
        tvs=[]; tls=[]
        try:
            l = round(v, 2)
        except:
            v = None
            l = 'Error'
        # Get axis limits
        lmin, lmax = limits[i]
        # More than lower limit include min
#        if (v is not None) and (v > 1.1*lmin): tvs.append(lmin); tls.append(lmin)
        # Add point marker
        tvs.append(v); tls.append(l)
        # Less than upper limit include max
#        if (v is not None) and (v < 0.9*lmax): tvs.append(lmax); tls.append(lmax)
        # Add to column ticks
        col_tick_vals.append(tvs); col_tick_labs.append(tls)

    r = Radar(titles=col_names)
    r.add(cols, "-", label=lab, lw=2)
    r.set_inversion([1,0,1,0,0])
    r.set_limits(limits)
    r.set_ticks(values=col_tick_vals, labels=col_tick_labs)
    r.plot()
    r.ax.legend(loc='upper center', fancybox=True, bbox_to_anchor=(0.5, 0.0))
    r.savefig(path)
    r.close()

    return

def run(params):

    assert params.input.dir, 'No directories provided'
    assert params.input.pdb_style
    # REMOVE THIS WHEN IMPLEMENTED
    assert params.selection.res_names
    if params.selection.res_names:      params.selection.__inject__("res_names_list", params.selection.res_names.split(','))
    else:                               params.selection.__inject__("res_names_list", None)
    if not os.path.exists(params.output.out_dir): os.mkdir(params.output.out_dir)
    scores_file =  os.path.join(params.output.out_dir, 'residue_scores.csv')
    images_dir =   os.path.join(params.output.out_dir, 'images')
    if not os.path.exists(images_dir): os.mkdir(images_dir)
    summary_file = os.path.join(params.output.out_dir, 'residue_scores.html')

    proc_dirs = []
    # Filter directories
    for dir in params.input.dir:
        if not os.path.isdir(dir): continue
        proc_dirs.append(dir)
    # Build arglist for mapping
    arg_list = []
    for dir in proc_dirs:
        new_table = prepare_table()
        arg_list.append((dir, params, new_table))
    # Score files in parallel
    returned_objs = libtbx.easy_mp.pool_map(fixed_func=process_folder_mp, args=arg_list, processes=params.settings.cpus)
    # Extract results
    data_table = None
    for obj in returned_objs:
        if isinstance(obj, pandas.DataFrame):
            if data_table is None:  data_table = obj
            else:                   data_table = data_table.append(obj, verify_integrity=True)
        else:
            print '{}: {}'.format(obj.__class__, obj.message)

    print '#=======================++>'
    data_table.to_csv(scores_file)
    print 'Output written to {}'.format(scores_file)

    if params.output.generate_summary:
        # Output Images
        all_images = []
        print 'Generating Output Images'
        for lab, row in data_table.iterrows():
            print 'Making: {}...'.format(lab)
            image_path = os.path.join(images_dir,'{}.png'.format(lab))
            make_image(lab=lab, row=row, path=image_path)
            all_images.append(image_path)

        print 'Generating output HTML'
        # Get template to be filled in
        template = BAMBOO_HTML_ENV.get_template('summary_page.html')
        # Output directory (for relative symlinks)
        out_dir  = os.path.abspath(os.path.dirname(summary_file))

        # ===========================================================>
        # Construct the data object to populate the template
        output_data = {}
        output_data['header'] = 'Residue Score Summaries'
        output_data['title'] = 'Residue Score Summaries'
        output_data['introduction'] = 'Model Quality and Validation checks.'
        # ===========================================================>
        # Header Images
        output_data['small_images'] = []
        for img in all_images:
            output_data['small_images'].append({ 'path': './'+os.path.relpath(path=img, start=out_dir),
                                                 'title': 'Scores for {}'.format(os.path.splitext(os.path.basename(img))[0]) })
        # ===========================================================>
        # Write Output
        with open(summary_file, 'w') as out_html:
            out_html.write(template.render(output_data))





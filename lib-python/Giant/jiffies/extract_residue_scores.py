import os, sys, copy, re, glob

import libtbx.phil

import iotbx.pdb
import pandas

from scitbx.array_family import flex
from Bamboo.Density.Edstats import edstats
from Giant.Utils.pdb import strip_pdb_to_input
from IPython import embed

#################################
import matplotlib
try:    matplotlib.use('Agg')
except: print 'Failed to load Agg as Backend for Matplotlib'
matplotlib.interactive(0)
from matplotlib import pyplot
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
    prefix = residue_scores.csv
        .type = path
    plot_graphs = True
        .type = bool
}
""")

def is_within(dist, coords_1, coords_2):
    """Checks if any of coords_1 is within dist of coords_2"""

    dist_sq = dist**2
    if len(coords_2) < len(coords_1):
        i1 = coords_2; i2 = coords_1
    else:
        i1 = coords_1; i2 = coords_2

    for c1 in i1:
        diffs = i2 - c1
        min_d_sq = flex.min(diffs.dot())
        if min_d_sq < dist_sq:
            return True
    return False

def sanitise_hierarchy(hierarchy):
    print 'Sanitising:'
    print hierarchy.atoms().set_chemical_element_simple_if_necessary()

def prepare_table(mtz_1=True, mtz_2=True, pdb_2=True):

    columns_1 = ['PDB','OCC','B-AV-QUERY','B-AV-CONTACTS','B-RATIO']
    columns_2 = ['MTZ','RSCC','RSR','B-AV','RSZO','RSZD']

    columns = []

    # PDB/MODEL Columns
    columns.extend(columns_1)
    if pdb_2:
        columns.extend([c+'-2' for c in columns_1])

    # MTZ/DENSITY Columns
    if mtz_1:
        columns.extend(columns_2)
    if mtz_2:
        columns.extend([c+'-2' for c in columns_2])

    empty_table = pandas.DataFrame(index=[], columns=columns)

    return empty_table

def run(params):

    assert params.input.dir, 'No directories provided'
    assert params.input.pdb_style
#    assert params.input.mtz_style

    # REMOVE THIS WHEN IMPLEMENTED
    assert params.selection.res_names

    if params.selection.res_names:      res_names = params.selection.res_names.split(',')
    else:                               res_names = None

    data_table = prepare_table( mtz_1=(True if params.input.mtz_style   else False),
                                mtz_2=(True if params.input.mtz_2_style else False),
                                pdb_2=(True if params.input.pdb_2_style else False))

    for dir in params.input.dir:

        if not os.path.isdir(dir): continue

        try:
            pdb  = glob.glob(os.path.join(dir, params.input.pdb_style))[0]
        except: print 'FAILED TO FIND FILE: {} in {}'.format(params.input.pdb_style, dir); continue
        try:
            if params.input.pdb_2_style:    pdb2 = glob.glob(os.path.join(dir, params.input.pdb_2_style))[0]
            else:                           pdb2 = None
        except: print 'FAILED TO FIND FILE: {} in {}'.format(params.input.pdb_2_style, dir); continue
        try:
            if params.input.mtz_style:      mtz1 = glob.glob(os.path.join(dir, params.input.mtz_style))[0]
            else:                           mtz1 = None
        except: print 'FAILED TO FIND FILE: {} in {}'.format(params.input.mtz_style, dir); continue
        try:
            if params.input.mtz_2_style:    mtz2 = glob.glob(os.path.join(dir, params.input.mtz_2_style))[0]
            else:                           mtz2 = None
        except: print 'FAILED TO FIND FILE: {} in {}'.format(params.input.mtz_2_style, dir); continue

        if   params.input.labels == 'basename':     lab = os.path.split_ext(os.path.basename(pdb))[0]
        elif params.input.labels == 'folder_name':  lab = os.path.split(pdb)[-2]

        print 'Processing: {}'.format(lab)
        print 'PDB:     {}'.format(pdb)

        # Extract Structure
        h_all = strip_pdb_to_input(pdb, remove_ter=True).hierarchy
        sanitise_hierarchy(h_all)
        h_all = h_all.select(h_all.atom_selection_cache().selection('not element H'))
        h_pro = h_all.select(h_all.atom_selection_cache().selection('pepnames'))
        h_bck = h_pro.select(h_pro.atom_selection_cache().selection('(name CA or name C or name O or name N)'))
        h_sch = h_pro.select(h_pro.atom_selection_cache().selection('not (name CA or name C or name O or name N)'))

        # Pull out residues to analyse
        if res_names:
            rg_for_analysis = [rg for rg in h_all.residue_groups() if [n for n in rg.unique_resnames() if n in res_names]]
        else:
            rg_for_analysis = h_all.residue_groups()

        # Check residues to analyse or skip
        if not rg_for_analysis:
            print 'No residues to analyse, skipping: {}'.format(pdb)
            continue

        # Extract PDB2
        if pdb2 is not None:
            print 'READING PDB2: {}'.format(pdb2)
            pdb2_h_all = iotbx.pdb.hierarchy.input(pdb2).hierarchy
            sanitise_hierarchy(pdb2_h_all)
            pdb2_h_all = pdb2_h_all.select(pdb2_h_all.atom_selection_cache().selection('not element H'))

        # Score MTZ1
        if mtz1 is not None:
            print 'SCORING MTZ1: {}'.format(mtz1)
            mtz1_edstats_scores = edstats(mtz_file=mtz1, pdb_file=pdb)
        else:
            mtz1_edstats_scores = None
        # Score MTZ2
        if mtz2 is not None:
            print 'SCORING MTZ2: {}'.format(mtz2)
            mtz2_edstats_scores = edstats(mtz_file=mtz2, pdb_file=pdb)
        else:
            mtz2_edstats_scores = None

        for rg_sel in rg_for_analysis:

            if len(rg_sel.unique_resnames()) != 1:
                print 'More than one residue name associated with residue group -- cannot process'

            # Create label for the output table
            row_label = (lab, rg_sel.unique_resnames()[0], rg_sel.parent().id, rg_sel.resseq_as_int(), rg_sel.icode)
            # Append empty row to output table
            data_table = data_table.append(pandas.DataFrame(index=[row_label]), verify_integrity=True)
            data_table.set_value(index=row_label, col='PDB', value=pdb)

            # Extract Model Scores - MTZ 1
            rg_atoms   = rg_sel.atoms()
            rg_coords  = rg_atoms.extract_xyz()
            rg_confs   = [c.altloc for c in rg.conformers()]
            # Select nearby atom_groups
            near_ags = [ag for ag in h_sch.atom_groups() if (   is_within(4, rg_coords, ag.atoms().extract_xyz())
                                                                and (ag not in rg_sel.atom_groups())                )]
            print '{} Atom Groups Selected {}'.format(len(near_ags), row_label)
#            for ag in near_ags:
#                print ag.id_str(), ag.altloc, [at.name for at in ag.atoms()]
            if near_ags:
                # Extract atoms from nearby groups
                near_ats = iotbx.pdb.hierarchy.af_shared_atom()
                [near_ats.extend(ag.atoms()) for ag in near_ags]
                # Calculate Occupancies of atoms
                res_occ = sum(set(rg_atoms.extract_occ()))
                if res_occ > 1: res_occ=None
                data_table.set_value(index=row_label, col='OCC',    value=res_occ )
                # Calculate B-factors of atoms
                res_mean_b = flex.mean_weighted(rg_atoms.extract_b(), rg_atoms.extract_occ())
                data_table.set_value(index=row_label, col='B-AV-QUERY',    value=res_mean_b )
                sch_mean_b = flex.mean_weighted(near_ats.extract_b(), near_ats.extract_occ())
                data_table.set_value(index=row_label, col='B-AV-CONTACTS', value=sch_mean_b )
                data_table.set_value(index=row_label, col='B-AV-RATIO', value=res_mean_b/sch_mean_b )

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

                    if not rg_sel2.detached_copy().is_similar_hierarchy(rg_sel.detached_copy()):
                        print [c.altloc for c in rg_sel.conformers()], [c.altloc for c in rg_sel2.conformers()]
                        assert rg_sel2.detached_copy().is_similar_hierarchy(rg_sel.detached_copy())

                    assert rg_sel.atoms().size()                  == rg_sel2.atoms().size()
                    assert list(rg_sel.atoms().extract_name()   ) == list(rg_sel2.atoms().extract_name()   )
                    assert list(rg_sel.atoms().extract_element()) == list(rg_sel2.atoms().extract_element())

                    rmsd = (rg_sel.atoms().extract_xyz() - rg_sel2.atoms().extract_xyz()).rms_length()
                    data_table.set_value(index=row_label, col='RES-RMSD', value=rmsd)
                except:
                    print 'Failed to calculate RMSD between models'

            # Extract Density Scores - MTZ 1
            if mtz1 is not None:
                data_table.set_value(index=row_label, col='MTZ', value=mtz1)
            if mtz1_edstats_scores is not None:
                ed_scores = mtz1_edstats_scores.scores[(rg_sel.unique_resnames()[0], rg_sel.parent().id, rg_sel.resseq_as_int(), rg_sel.icode)]

                data_table.set_value(index=row_label, col='RSCC', value=ed_scores['CCSa'])
                data_table.set_value(index=row_label, col='RSR',  value=ed_scores['Ra']  )
                data_table.set_value(index=row_label, col='B-AV', value=ed_scores['BAa'] )
                data_table.set_value(index=row_label, col='RSZO', value=ed_scores['ZOa'] )
                data_table.set_value(index=row_label, col='RSZD', value=ed_scores['ZDa'] )

            # Extract Density Scores - MTZ 2
            if mtz2 is not None:
                data_table.set_value(index=row_label, col='MTZ-2', value=mtz2)
            if mtz2_edstats_scores is not None:
                ed_scores = mtz2_edstats_scores.scores[(rg_sel.unique_resnames()[0], rg_sel.parent().id, rg_sel.resseq_as_int(), rg_sel.icode)]

                data_table.set_value(index=row_label, col='RSCC-2', value=ed_scores['CCSa'])
                data_table.set_value(index=row_label, col='RSR-2',  value=ed_scores['Ra']  )
                data_table.set_value(index=row_label, col='B-AV-2', value=ed_scores['BAa'] )
                data_table.set_value(index=row_label, col='RSZO-2', value=ed_scores['ZOa'] )
                data_table.set_value(index=row_label, col='RSZD-2', value=ed_scores['ZDa'] )

            print data_table.loc[[row_label]].T
            print '#=======================++>'
            data_table.to_csv(params.output.prefix)

    print '#=======================++>'
    data_table.to_csv(params.output.prefix)
    print 'Output written to {}'.format(params.output.prefix)

    if params.output.plot_graphs:

        fig = pyplot.figure()
        pyplot.title('RSCC HISTOGRAM')
        pyplot.hist(x=data_table['RSCC'], bins=30, normed=True)
        pyplot.xlabel('RSCC')
        pyplot.ylabel('DENSITY')
        pyplot.tight_layout()
        pyplot.savefig(os.path.splitext(params.output.prefix)[0]+'-rscc-hist.png')
        pyplot.close(fig)

        fig = pyplot.figure()
        pyplot.title('RSZD HISTOGRAM')
        pyplot.hist(x=data_table['RSZD'], bins=30, normed=True)
        pyplot.xlabel('RSZD')
        pyplot.ylabel('DENSITY')
        pyplot.tight_layout()
        pyplot.savefig(os.path.splitext(params.output.prefix)[0]+'-rszd-hist.png')
        pyplot.close(fig)

        fig = pyplot.figure()
        pyplot.title('RSZD HISTOGRAM')
        pyplot.plot(data_table['OCC'], data_table['RSCC'], 'go')
        pyplot.xlabel('OCCUPANCY')
        pyplot.ylabel('RSCC')
        pyplot.tight_layout()
        pyplot.savefig(os.path.splitext(params.output.prefix)[0]+'-rscc-v-occ.png')
        pyplot.close(fig)

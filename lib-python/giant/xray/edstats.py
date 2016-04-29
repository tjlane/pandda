import pandas

def extract_residue_group_density_scores(residue_group, edstats_scores, data_table=None, rg_label=None, column_suffix=''):
    """Extract density quality metrics for a residue group from precalculated edstats scores"""

    rg_sel = residue_group
    # Set defaults
    if rg_label is None:    rg_label = (rg_sel.unique_resnames()[0]+'-'+rg_sel.parent().id+'-'+rg_sel.resseq+rg_sel.icode).replace(' ','')
    if data_table is None:  data_table = pandas.DataFrame(index=[rg_label], column=[])
    # Check validity
    if len(rg_sel.unique_resnames()) != 1: raise Exception(rg_label+': More than one residue name associated with residue group -- cannot process')

    # Extract residue scores
    ed_scores = edstats_scores.scores[(rg_sel.unique_resnames()[0], rg_sel.parent().id, rg_sel.resseq_as_int(), rg_sel.icode)]
    # Append scores to data_table
    data_table.set_value(index=rg_label, col='RSCC'+column_suffix, value=ed_scores['CCSa'])
    data_table.set_value(index=rg_label, col='RSR' +column_suffix, value=ed_scores['Ra']  )
    data_table.set_value(index=rg_label, col='B_AV'+column_suffix, value=ed_scores['BAa'] )
    data_table.set_value(index=rg_label, col='RSZO'+column_suffix, value=ed_scores['ZOa'] )
    data_table.set_value(index=rg_label, col='RSZD'+column_suffix, value=ed_scores['ZDa'] )

    return data_table






from bamboo.macro.utils import get_residue_labels

def extract_scores(scores, label):
    """Pull out the score specified by 'label'"""

    # Validate label (NAME-CHAIN-NUMBER-INSCODE)
    if label[3] == '':
        label = list(label)
        label[3] = ' '
        label = tuple(label)

    # Extract Scores
    try:
        labscores = scores[label]
    except KeyError:
        if label[1] == ' ':
            reflabel=list(label)
            reflabel[1] = '_'
            if tuple(reflabel) in scores:
                labscores = scores[tuple(reflabel)]
            else:
                reflabel[1] = 'X'
                labscores = scores[tuple(reflabel)]
        else:
            print 'Scores not found in EDSTATS output!'
            raise

    return labscores

def score_against_density_by_subset_file(mtz_file, pdb_file, subpdb):
    """Scores all the residues with edstats and returns the scores for the residues present in `subpdb`"""

    # Get the labels
    labels = get_residue_labels(subpdb)
    # Score by labels
    labels, summary, scores = score_against_density_by_res_labels(mtz_file, pdb_file, labels)

    return labels, summary, scores

def score_against_density_by_res_labels(mtz_file, pdb_file, labels):
    """Scores all the residues with edstats and returns the scores for the residues identified by `labels`"""

    # Get the scores
    edstats, summary = score_file_with_edstats(mtz_file, pdb_file)
    # Score by label (list of tuples)
    scores = dict([(label,extract_scores(edstats.scores, label)) for label in labels])
    # Add summary percentiles for the labelled residue
    summary.percentiles = dict([(label,get_labelled_residue_percentiles(label, edstats.scores, summary)) for label in labels])

    return labels, summary, scores

def get_labelled_residue_percentiles(label, edstats_scores):
    """Takes the scores and summary object from edstats and compares the residue given by 'label' to the rest of the protein"""

    # Get the res scores given by the label
    res_score = extract_scores(edstats_scores, label)

    # Get the total number of residues in the file
    num_res = edstats_scores.columns.size

    # Get the list of ZD scores
    all_ZD =       edstats_scores.loc['ZDa']
    # Get the list of ZD+ scores
    all_ZD_plus =  edstats_scores.loc['ZD+a']
    # Get the list of ZD- scores
    all_ZD_minus = edstats_scores.loc['ZD-a']
    # Get the list of ZO scores
    all_ZO =       edstats_scores.loc['ZOa']
    # Get the list of RSR Scores
    all_RSR =      edstats_scores.loc['Ra']
    # Get the list of RSCC Scores
    all_RSCC =     edstats_scores.loc['CCPa']
    # Get the list of B-factor Scores
    all_B =        edstats_scores.loc['BAa']

    # Get the percentiles of the ligand scores
    percentiles = {}
    percentiles['%ile_ZD']  =  (len([i for i in all_ZD       if i < res_score['ZDa']])      + len([i for i in all_ZD       if i <= res_score['ZDa']])) * 50 / float(num_res)
    percentiles['%ile_ZD+']  = (len([i for i in all_ZD_plus  if i < res_score['ZD+a']])     + len([i for i in all_ZD_plus  if i <= res_score['ZD+a']])) * 50 / float(num_res)
    percentiles['%ile_ZD-']  = (len([i for i in all_ZD_minus if i > res_score['ZD-a']])     + len([i for i in all_ZD_minus if i >= res_score['ZD-a']])) * 50 / float(num_res)
    percentiles['%ile_ZO']   = (len([i for i in all_ZO       if i > res_score['ZOa']])      + len([i for i in all_ZO       if i >= res_score['ZOa']])) * 50 / float(num_res)
    percentiles['%ile_RSR']  = (len([i for i in all_RSR      if i < res_score['Ra']])       + len([i for i in all_RSR      if i <= res_score['Ra']])) * 50 / float(num_res)
    percentiles['%ile_RSCC'] = (len([i for i in all_RSCC     if i > res_score['CCPa']])     + len([i for i in all_RSCC     if i >= res_score['CCPa']])) * 50 / float(num_res)
    percentiles['%ile_B']    = (len([i for i in all_B        if i < res_score['BAa']])      + len([i for i in all_B        if i <= res_score['BAa']])) * 50 / float(num_res)

    percentiles['%ile_ZOs']   = (len([i for i in all_ZO       if i > res_score['ZOs']])      + len([i for i in all_ZO       if i >= res_score['ZOs']])) * 50 / float(num_res)
    percentiles['%ile_RSCCs'] = (len([i for i in all_RSCC     if i > res_score['CCPs']])     + len([i for i in all_RSCC     if i >= res_score['CCPs']])) * 50 / float(num_res)
    percentiles['%ile_Bs']    = (len([i for i in all_B        if i < res_score['BAs']])      + len([i for i in all_B        if i <= res_score['BAs']])) * 50 / float(num_res)

    return percentiles


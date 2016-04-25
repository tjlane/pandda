from bamboo.macro.utils import get_residue_labels
from bamboo.edstats import score_with_edstats_to_dict

class EdstatsLogSummary:
    """Class to process, store and present the output from edstats.pl"""

    def __init__(self, logtext):
        self.logtext = logtext
        self.mean_biso = None
        self.rms_scores = {'atm_all_ZD+':None,'atm_all_ZD-':None,'atm_all_ZO':None,'atm_main_ZD+':None,'atm_main_ZD-':None,'atm_main_ZO':None,'atm_side_ZD+':None,'atm_side_ZD-':None,'atm_side_ZO':None,'het_all_ZD+':None,'het_all_ZD-':None,'het_all_ZO':None,'het_main_ZD+':None,'het_main_ZD-':None,'het_main_ZO':None,'het_side_ZD+':None,'het_side_ZD-':None,'het_side_ZO':None}

        # Populate fields
        self._parselogtext()

    def _parselogtext(self):
        """Read and sort through log input"""

        output = [line for line in self.logtext.split('\n') if line]
        for i, line in enumerate(output):
            if line.startswith('Overall mean Biso:'):
                self.mean_biso = float(line.replace('Overall mean Biso:',''))
            elif line.startswith('RMS Z-scores for protein residues:'):
                assert ['Main','Side','All'] == output[i+1].strip().split(), 'COLUMN HEADINGS ARE NOT WHAT IS EXPECTED! {!s}'.format(line)
                for incr in [2,3,4]:
                    score, main_val, side_val, all_val = output[i+incr].strip().split()
                    self.rms_scores['_'.join(['atm','main',score])] = float(main_val)
                    self.rms_scores['_'.join(['atm','side',score])] = float(side_val)
                    self.rms_scores['_'.join(['atm','all',score])] = float(all_val)
            elif line.startswith('RMS Z-scores for hetero residues:'):
                assert ['Main','Side','All'] == output[i+1].strip().split(), 'COLUMN HEADINGS ARE NOT WHAT IS EXPECTED! {!s}'.format(line)
                for incr in [2,3,4]:
                    score, main_val, side_val, all_val = output[i+incr].strip().split()
                    self.rms_scores['_'.join(['het','main',score])] = float(main_val)
                    self.rms_scores['_'.join(['het','side',score])] = float(side_val)
                    self.rms_scores['_'.join(['het','all',score])] = float(all_val)

def extract_scores(scores, label):
    """Pull out the score specified by 'label'"""

    # Validate label
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

def score_against_density_by_subset_file(mtzpath, pdbpath, subpdb):
    """Scores all the residues with edstats and returns the scores for the residues present in `subpdb`"""

    # Get the labels
    labels = get_residue_labels(subpdb)
    # Score by labels
    labels, summary, scores = score_against_density_by_res_labels(mtzpath, pdbpath, labels)

    return labels, summary, scores

def score_against_density_by_res_labels(mtzpath, pdbpath, labels):
    """Scores all the residues with edstats and returns the scores for the residues identified by `labels`"""

    # Get the scores
    summary, rawscores = score_file_with_edstats(mtzpath, pdbpath)
    # Score by label (list of tuples)
    scores = dict([(label,extract_scores(rawscores, label)) for label in labels])
    # Add summary percentiles for the labelled residue
    summary.percentiles = dict([(label,get_labelled_residue_percentiles(label, rawscores, summary)) for label in labels])

    return labels, summary, scores

def score_file_with_edstats(mtzpath, pdbpath):
    """Score pdb file against electron density"""

    # Score the complex with Edstats
    scores, command = score_with_edstats_to_dict(mtzpath, pdbpath)
    if (not scores) and command.err:
        raise ScoringError('EDSTATS ERROR : {!s}'.format(command.err))

    # Process the std out of the program
    summary = edstatsLogSummary(command.out)

    return summary, scores

def get_labelled_residue_percentiles(label, edstats_scores, edstats_summary):
    """Takes the scores and summary object from edstats and compares the residue given by 'label' to the rest of the protein"""

    # Get the res scores given by the label
    res_score = extract_scores(edstats_scores, label)

    # Get the total number of residues in the file
    num_res = len(edstats_scores)

    # Get the list of ZD scores
    all_ZD =       [score['ZDa'] for score in edstats_scores.values()]
    # Get the list of ZD+ scores
    all_ZD_plus =  [score['ZD+a'] for score in edstats_scores.values()]
    # Get the list of ZD- scores
    all_ZD_minus = [score['ZD-a'] for score in edstats_scores.values()]
    # Get the list of ZO scores
    all_ZO =       [score['ZOa'] for score in edstats_scores.values()]
    # Get the list of RSR Scores
    all_RSR =      [score['Ra'] for score in edstats_scores.values()]
    # Get the list of RSCC Scores
    all_RSCC =     [score['CCPa'] for score in edstats_scores.values()]
    # Get the list of B-factor Scores
    all_B =        [score['BAa'] for score in edstats_scores.values()]

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
    percentiles['$ile_Bs']    = (len([i for i in all_B        if i < res_score['BAs']])      + len([i for i in all_B        if i <= res_score['BAs']])) * 50 / float(num_res)

    return percentiles


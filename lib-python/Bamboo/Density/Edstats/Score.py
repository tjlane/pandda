#! /usr/local/python/python2.7.3-64bit/bin/python

from Bamboo.Macro.Utils import get_residue_labels
from Bamboo.Density.Edstats.Utils import score_with_edstats_to_dict, extract_scores, edstatsSummary, get_labelled_residue_percentiles

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
    summary = edstatsSummary(command.out)

    return summary, scores


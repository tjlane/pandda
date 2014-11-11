#! /usr/local/python/python2.7.3-64bit/bin/python

# Adapted from Code by Sebastian Kelm

import os, sys, tempfile

from Bamboo.Common.Command import CommandManager
from Bamboo.Utils.MTZ import mtzFile
from Bamboo.Macro.Utils import get_residue_labels

class edstatsSummary:
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

def score_with_edstats_to_dict(mtzpath, pdbpath):
    """Scores residues against density, then converts to dict"""

    # Generate the list of the EDSTATS scores for each residue
    scores, header, command = score_with_edstats_to_list(mtzpath, pdbpath)
    # Create dict for the mapping between residue_id and scores
    mapping = {}
    for label, values in scores:
        hash = {}
        for k,v in zip(header,values):
            hash[k] = v
        mapping[label] = hash

    return mapping, command

def score_with_edstats_to_list(mtzpath, pdbpath):
    """Scores residues against density, then returns list"""

    assert os.path.exists(mtzpath), 'MTZ FILE FOR EDSTATS DOES NOT EXIST! {!s}'.format(mtzpath)
    assert os.path.exists(pdbpath), 'PDB FILE FOR EDSTATS DOES NOT EXIST! {!s}'.format(mtzpath)

    # Create a file handle and path for the output
    temp_handle, temp_path = tempfile.mkstemp(suffix='.table', prefix='edstats_')

    # Find the labels in the mtzfile
    file_obj = mtzFile(mtzpath)
    f_label = file_obj.label.f
    if not f_label:
        raise ReflectionException('MTZ Summary ERROR: No F Label Found in MTZ File: {!s}'.format(mtzpath))

    # Run EDSTATS on the files
    try:
        # Initialise Command Manager to run edstats
        command = CommandManager('edstats.pl')
        args = ['-hklin',mtzpath,'-xyzin',pdbpath,'-output',temp_path,'-noerror','-flabel',f_label]
        command.SetArguments(args)
        command.SetParameters(timeout=600)
        command.Run()
        # Read the output
        with os.fdopen(temp_handle) as f:
            output = f.read().strip().replace('\r\n','\n').replace('\r','\n').splitlines()
    finally:
        os.remove(temp_path)

    # Process the output header
    if output:
        header = output.pop(0).split()
        assert header[:3] == ['RT', 'CI', 'RN'], 'EDSTATS OUTPUT HEADERS ARE NOT AS EXPECTED! {!s}'.format(output)
        num_fields = len(header)
        header = header[3:]
    else:
        header = []

    # List to be returned
    outputdata = []

    # Process the rest of the data
    for line in output:
        line = line.strip()
        if not line:
            continue

        fields = line.split()
        if len(fields) != num_fields:
            raise ValueError("Error Parsing EDSTATS output: Header & Data rows have different numbers of fields")

        # Get and process the residue information
        residue, chain, resnum = fields[:3]
        try:
            resnum = int(resnum)
            inscode = ' '
        except ValueError:
            inscode = resnum[-1:]
            resnum = int(resnum[:-1])

        # Remove the processed columns
        fields = fields[3:]

        # Process the other columns (changing n/a to None and value to int)
        for i, x in enumerate(fields):
            if x == 'n/a':
                fields[i] = None
            else:
                try:
                    fields[i] = int(x)
                except ValueError:
                    fields[i] = float(x)

        outputdata.append([(residue, chain, resnum, inscode),fields])

    return outputdata, header, command

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

def get_z_values_for_scores(edstats_scores, edstats_summary):
    pass

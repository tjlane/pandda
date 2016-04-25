# Adapted from Code by Sebastian Kelm

import os, sys, tempfile
import pandas

from bamboo.common.command import CommandManager
from bamboo.utils.mtz import MtzFile

class edstats(object):
    def __init__(self, mtz_file, pdb_file):
        scores, command = score_with_edstats_to_dict(mtzpath=mtz_file, pdbpath=pdb_file)
        self.scores = pandas.DataFrame.from_dict(scores)
        self._command = command

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
        command.add_command_line_arguments(['-hklin',mtzpath,'-xyzin',pdbpath,'-output',temp_path,'-noerror','-flabel',f_label])
        command.set_timeout(timeout=600)
        command.run()
        # Read the output
        with os.fdopen(temp_handle) as f:
            output = f.read().strip().replace('\r\n','\n').replace('\r','\n').splitlines()
        command.file_output = output
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



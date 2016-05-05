import os,sys,re

from bamboo.common.file import FileObj
from bamboo.common.command import CommandManager

################################################################################################################
### CLASSES                                                                                                    #
################################################################################################################

class MtzMeta(object):
    """Object to hold the meta data for an MTZ file"""

    def __init__(self):
        self.reslow = None
        self.reshigh = None
        self.spacegroup = None
        self.spacegroupno = None
        self.cell = None
        self.cryst = None
    def __str__(self):
        return str(self.__dict__)

class ColumnLabels(object):
    """Object to contain the column labels of an MTZ file"""

    def __init__(self):
        self.f = None
        self.sigf = None
        self.phas = None
        self.free = None
        self.comp_f = None
        self.comp_p = None
        self.diff_f = None
        self.diff_p = None
    def __str__(self):
        return str(self.__dict__)

class MtzSummary(object):
    """Class for summarising an MTZ file"""

    def __init__(self, mtz_file):
        self.path = os.path.abspath(mtz_file)
        self.file = FileObj(mtz_file)
        # Read the mtz_file and create summary
        self.summary = get_mtz_summary_dict(mtz_file)
        # Record information from summary
        self.label = ColumnLabels()
        self.label.f    = self.summary['f_labels'][0]    if self.summary['f_labels']    else None
        self.label.sigf = self.summary['sigf_labels'][0] if self.summary['sigf_labels'] else None
        self.label.i    = self.summary['i_labels'][0]    if self.summary['i_labels']    else None
        self.label.sigi = self.summary['sigi_labels'][0] if self.summary['sigi_labels'] else None
        self.label.phas = self.summary['p_labels'][0]    if self.summary['p_labels']    else None
        self.label.free = self.summary['r_flags'][0]     if self.summary['r_flags']     else None
        self.label.comp_f = self.summary['wtmap_f_comp'][0] if self.summary['wtmap_f_comp'] else None
        self.label.comp_p = self.summary['wtmap_p_comp'][0] if self.summary['wtmap_p_comp'] else None
        self.label.diff_f = self.summary['wtmap_f_diff'][0] if self.summary['wtmap_f_diff'] else None
        self.label.diff_p = self.summary['wtmap_p_diff'][0] if self.summary['wtmap_p_diff'] else None
        # Get some Meta
        self.data = MtzMeta()
        self.data.reslow        = self.summary['reslow']        if self.summary['reslow']       else None
        self.data.reshigh       = self.summary['reshigh']       if self.summary['reshigh']      else None
        self.data.spacegroup    = self.summary['spacegroup']    if self.summary['spacegroup']   else None
        self.data.spacegroupno  = self.summary['spacegroupno']  if self.summary['spacegroupno'] else None
        self.data.cell          = self.summary['cell']
        self.data.cryst         = dict(zip(['A','B','C','ALPHA','BETA','GAMMA','SG'],self.data.cell+[self.data.spacegroup]))
        self.data.a,     self.data.b,    self.data.c     = self.summary['cell'][0:3]
        self.data.alpha, self.data.beta, self.data.gamma = self.summary['cell'][3:6]

################################################################################################################
### FUNCTIONS                                                                                                  #
################################################################################################################

def get_mtz_summary_dict(mtz_file):
    """Get an MTZ Summary"""

    # Extract Contents of MTZ
    MTZDMP = CommandManager('mtzdmp')
    MTZDMP.add_command_line_arguments(mtz_file)
    MTZDMP.run()
    # Check for errors
    if MTZDMP.process.returncode != 0:
        raise RuntimeError('mtzdmp failed to read file {!s}:\nReturn: {!s}\nOut: {!s}'.format(mtz_file,MTZDMP.process.returncode,MTZDMP.output))

    # Create empty dict to contain the summary
    summary = {}

    # Get the resolution range
    regex = re.compile('\*  Resolution Range :.*\n.*\n.*\((.*)A \)\n')
    matches = regex.findall(MTZDMP.output)
    assert matches, 'No Resolution Range found in MTZFile {!s}'.format(mtz_file)
    assert len(matches)==1, 'Too many matching lines found for Resolution Range in MTZFile {!s}\n\t{!s}'.format(mtz_file,matches)
    summary['reslow'],summary['reshigh'] = map(float,matches[0].replace(' ','').split('-'))

    # Get the Number of Columns
    regex = re.compile('\* Number of Columns =(.*)\n')
    matches = regex.findall(MTZDMP.output)
    assert matches, 'No Number of Columns found for {!s}'.format(mtz_file)
    assert len(matches)==1, 'Too many matching lines found for Number of Columns in MTZFile {!s}\n\t{!s}'.format(mtz_file,matches)
    summary['numcols'] = int(matches[0].strip())

    # Get the Number of Reflections
    regex = re.compile('\* Number of Reflections =(.*)\n')
    matches = regex.findall(MTZDMP.output)
    assert matches, 'No Number of Reflections found for {!s}'.format(mtz_file)
    assert len(matches)==1, 'Too many matching lines found for Number of Reflections in MTZFile {!s}\n\t{!s}'.format(mtz_file,matches)
    summary['numreflections'] = int(matches[0].strip())

    # Get the Column Labels
    regex = re.compile('\* Column Labels :.*\n.*\n(.*)\n')
    matches = regex.findall(MTZDMP.output)
    assert matches, 'No Column Labels found for {!s}'.format(mtz_file)
    assert len(matches)==1, 'Too many matching lines found for Column Headings in MTZFile {!s}\n\t{!s}'.format(mtz_file,matches)
    summary['colheadings'] = matches[0].strip().split()

    # Get the Column Types
    regex = re.compile('\* Column Types :.*\n.*\n(.*)\n')
    matches = regex.findall(MTZDMP.output)
    assert matches, 'No Column Types found for {!s}'.format(mtz_file)
    assert len(matches)==1, 'Too many matching lines found for Column Types in MTZFile {!s}\n\t{!s}'.format(mtz_file,matches)
    summary['coltypes'] = matches[0].strip().split()

    # Get the different datasets
    regex = re.compile('\* Associated datasets :.*\n.*\n(.*)\n')
    matches = regex.findall(MTZDMP.output)
    assert matches, 'No Dataset Numbers found for {!s}'.format(mtz_file)
    assert len(matches)==1, 'Too many matching lines found for Dataset Numbers in MTZFile {!s}\n\t{!s}'.format(mtz_file,matches)
    summary['coldatasets'] = matches[0].strip().split()

    # Get the Spacegroup
    regex = re.compile('\* Space group =.*\'(.*)\'.\(number(.*)\)')
    matches = regex.findall(MTZDMP.output)
    assert matches, 'No Space Group found for {!s}'.format(mtz_file)
    assert len(matches)==1, 'Too many matching lines found for Spacegroup in MTZFile {!s}\n\t{!s}'.format(mtz_file,matches)
    summary['spacegroup'] = matches[0][0].strip()
    summary['spacegroupno'] = int(matches[0][1].strip())

    # Get the Cell Dimensions
    regex = re.compile('\* Cell Dimensions :.*\n.*\n(.*)\n')
    matches = regex.findall(MTZDMP.output)
    assert matches, 'No Cell Dimensions found for {!s}'.format(mtz_file)
    assert len(matches)==1, 'Too many matching lines found for Cell Dimensions in MTZFile {!s}\n\t{!s}'.format(mtz_file,matches)
    summary['cell'] = map(float,matches[0].split())

    # Get the Cell Dimensions

    # Sort the flags
    # F: Amplitudes, H: Miller Indices, I: Integers, P: Phases, Q: Standard Deviations, J: Intensities
    sfac = [summary['colheadings'][i] for i,type in enumerate(summary['coltypes']) if type=='F']
    inty = [summary['colheadings'][i] for i,type in enumerate(summary['coltypes']) if type=='J']
    mill = [summary['colheadings'][i] for i,type in enumerate(summary['coltypes']) if type=='H']
    ints = [summary['colheadings'][i] for i,type in enumerate(summary['coltypes']) if type=='I']
    phas = [summary['colheadings'][i] for i,type in enumerate(summary['coltypes']) if type=='P']
    sdev = [summary['colheadings'][i] for i,type in enumerate(summary['coltypes']) if type=='Q']

    # 2FOFC cols
    wt_f_map_opts = ['2FOFCWT','FWT']
    wt_p_map_opts = ['PH2FOFCWT','PHWT','PHFWT']
    # FOFC cols
    wt_f_dif_opts = ['FOFCWT','DELFWT']
    wt_p_dif_opts = ['PHFOFCWT','DELPHWT','PHDELWT']
    # Record 2FOFC pair (composite map)
    summary['wtmap_f_comp'] = [s for s in sfac if s in wt_f_map_opts]
    summary['wtmap_p_comp'] = [p for p in phas if p in wt_p_map_opts]
    # Record FOFC pair (different map)
    summary['wtmap_f_diff'] = [s for s in sfac if s in wt_f_dif_opts]
    summary['wtmap_p_diff'] = [p for p in phas if p in wt_p_dif_opts]

    # Find the main F col
    summary['f_labels'] = [s for s in sfac if (s in ['F','FP','FCTR','FOSC','FOBS'] \
                            or s.startswith('F_') or s.startswith('FP_') \
                            or s.upper().startswith('F-OBS') or s.upper().startswith('FOUT_'))]
    if not summary['f_labels']:
        summary['f_labels'] = [s for s in sfac if s.startswith('FP')]
    # Find the main SIGF col
    summary['sigf_labels'] = [s for s in sdev if (s in ['SIGF','SIGFP','SIGFOBS'] \
                            or s.startswith('SIGF_') or s.startswith('SIGFP_') \
                            or s.upper().startswith('SIGF-OBS') or s.upper().startswith('SIGFOUT_'))]
    if not summary['sigf_labels']:
        summary['sigf_labels'] = [s for s in sdev if s.startswith('SIGFP')]
    # Find the I cols
    summary['i_labels'] = [s for s in inty if s.startswith('I')]
    # Find the SIGI cols
    summary['sigi_labels'] = [s for s in sdev if s.startswith('SIGI')]
    # Find the F_calc
    summary['f_calcs'] = [s for s in sfac if (s=='FC' or s=='FMODEL' or s.startswith('FC_') or s.upper().startswith('F-MODEL'))]
    # Find the PHI_calcs
    summary['p_calcs'] = [p for p in phas if (p=='PHIC' or p=='PHIFMODEL' or p.startswith('PHIC_') or p.upper().startswith('PHIF-MODEL'))]
    # Find the main phase col
    summary['p_labels'] = summary['p_calcs']
    # Find the RFree Flag
    summary['r_flags'] = [r for r in ints if ('RFREE' in r.upper() or 'FREER' in r.upper() or r=='FREE' or r.upper().startswith('R-FREE'))]

    # XXX Will probably delete this later
    summary['unknown'] = [c for c in summary['colheadings'] if c not in summary['f_labels']+summary['sigf_labels']+summary['f_calcs']+summary['p_calcs']+summary['p_labels']+summary['r_flags']+summary['wtmap_f_comp']+summary['wtmap_p_comp']+summary['wtmap_f_diff']+summary['wtmap_p_diff']]

    return summary

def get_mtz_resolution(mtz_file):
    """Gets the max resolution from the file"""

    # Extract Contents of MTZ
    MTZDMP = CommandManager('mtzdmp')
    MTZDMP.add_command_line_arguments(mtz_file)
    MTZDMP.run()
    # Check for errors
    if MTZDMP.process.returncode != 0:
        raise RuntimeError('mtzdmp failed to read file {!s}:\nReturn: {!s}\nOut: {!s}'.format(mtz_file,MTZDMP.process.returncode,MTZDMP.output))
    # Search for the Column Headings
    regex = re.compile('\*  Resolution Range :.*\n.*\n.*\((.*)A \)\n')
    matches = regex.findall(MTZDMP.output)
    # Check for validity of matches
    assert matches, 'No Resolution Range found in MTZFile {!s}'.format(mtz_file)
    assert len(matches)==1, 'Too many matching lines found for Column Headings in MTZFile {!s}\n\t{!s}'.format(mtz_file,matches)
    # Return
    return map(float,matches[0].replace(' ','').split('-'))


import os,sys,re

from bamboo.common.file import FileObj
from bamboo.macro.molecule import MacroMol

################################################################################################################
### CLASSES                                                                                                    #
################################################################################################################

class PdbMeta(object):
    """Object to hold the meta data for a PDB file"""

    def __init__(self):
        self.reslow = None
        self.reshigh = None
        self.spacegroup = None
        self.spacegroupno = None
        self.cell = None
        self.cryst = None
        self.rfree = None
        self.rwork = None
    def __str__(self):
        return str(self.__dict__)

class BFactorSummary(object):
    """Class for summarising the b-factors of a residue"""
    def __init__(self, chain, resid, main_bav, main_brms, side_bav, side_brms, all_bav, all_brms):
        self.chain = chain
        self.resid = resid
        self.main_bav  = main_bav
        self.main_brms = main_brms
        self.side_bav  = side_bav
        self.side_brms = side_brms
        self.all_bav   = all_bav
        self.all_brms  = all_brms

class PdbSummary(object):
    """Class for summarising an PDBfile"""
    def __init__(self, pdb_file, parse=False):
        self.path = os.path.abspath(pdb_file)
        self.file = FileObj(pdb_file)
        # Read the pdb_file and create a summary
        self.summary = get_pdb_summary_dict(pdb_file)
        # Get some Meta
        self.data = PdbMeta()
        self.data.reslow        = self.summary['reslow']        if self.summary['reslow']       else None
        self.data.reshigh       = self.summary['reshigh']       if self.summary['reshigh']      else None
        self.data.spacegroup    = self.summary['spacegroup']    if self.summary['spacegroup']   else None
        self.data.spacegroupno  = self.summary['spacegroupno']  if self.summary['spacegroupno'] else None
        self.data.cell          = self.summary['cell']
        self.data.cryst         = dict(zip(['A','B','C','ALPHA','BETA','GAMMA','SG'],self.data.cell+[self.data.spacegroup]))
        self.data.a,     self.data.b,    self.data.c     = self.summary['cell'][0:3]
        self.data.alpha, self.data.beta, self.data.gamma = self.summary['cell'][3:6]
        # Create MacroMol object
        if parse:
            self.macromol = MacroMol(pdbpath)
        else:
            self.macromol = None

################################################################################################################
### FUNCTIONS                                                                                                  #
################################################################################################################

def get_pdb_summary_dict(pdb_file):
    """Get a pdb summary from the header information"""

    # Extract contents of PDB
    pdbcontents = open(pdb_file, 'r').read()

    # Create empty dict to contain the summary
    summary = {}

    # Get the Refinement High Resolution Limit
    regex = re.compile('\nREMARK.*RESOLUTION RANGE HIGH \(ANGSTROMS\) :(.*)\n')
    matches = regex.findall(pdbcontents)
#    assert matches, 'No High Resolution Limit found for {!s}'.format(pdb_file)
    assert len(matches)<2, 'Too many matching lines found for High Res Limit in PDBFile {!s}\n\t{!s}'.format(pdb_file,matches)
    summary['reshigh'] = float(matches[0]) if matches else None

    # Get the Refinement Low Resolution Limit
    regex = re.compile('\nREMARK.*RESOLUTION RANGE LOW  \(ANGSTROMS\) :(.*)\n')
    matches = regex.findall(pdbcontents)
#    assert matches, 'No Low Resolution Limit found for {!s}'.format(pdb_file)
    assert len(matches)<2, 'Too many matching lines found for Low Res Limit in PDBFile {!s}\n\t{!s}'.format(pdb_file,matches)
    summary['reslow'] = float(matches[0]) if matches else None

    # Get the R Free
    regex = re.compile('\nREMARK.*\ \ FREE R VALUE\ *:(.*)\n')
    matches = regex.findall(pdbcontents)
#    assert matches, 'No R FREE found for {!s}'.format(pdb_file)
    assert len(matches)<2, 'Too many matching lines found for R FREE in PDBFile {!s}\n\t{!s}'.format(pdb_file,matches)
    summary['rfree'] = float(matches[0]) if matches else None

    # Get the R Work
    regex = re.compile('\nREMARK.*\ \ R VALUE.*\(WORKING SET\) :(.*)\n')
    matches = regex.findall(pdbcontents)
#    assert matches, 'No R WORK found for {!s}'.format(pdb_file)
    assert len(matches)<2, 'Too many matching lines found for R WORK in PDBFile {!s}\n\t{!s}'.format(pdb_file,matches)
    summary['rwork'] = float(matches[0]) if matches else None

    # Get the Cryst Line
    regex = re.compile('\n(CRYST1.*)\n')
    matches = regex.findall(pdbcontents)
    assert matches, 'No CRYST LINE found for {!s}'.format(pdb_file)
    assert len(matches)<2, 'Too many matching lines found for CRYST LINE in PDBFile {!s}\n\t{!s}'.format(pdb_file,matches)
    assert len(matches[0]) > 54, 'Cryst Line is not long enough in PDBFile {!s}\n\t{!s}'.format(pdb_file, matches)
    summary['cell'] = map(float,[matches[0][6:15], matches[0][15:24], matches[0][24:33], matches[0][33:40], matches[0][40:47], matches[0][47:54]])
    summary['spacegroup'] = matches[0][55:66].strip()
    summary['spacegroupno'] = None

    return summary

def get_pdb_cryst_line(pdb_file):
    """Get the cryst Line from the PDB File"""

    cryst = [l for l in open(pdb_file,'r').readlines() if l.startswith('CRYST')]
    # There should only be one CRYST line
    if len(cryst)!=1:
        return ('More than one cryst line found',cryst)
    else:
        cryst = cryst[0]
    # Check for validity of CRYST line
    if not cryst[0:6]=='CRYST1':
        return ('CRYST line does not begin with \'CRYST1\'',cryst)
    elif len(cryst)<56:
        return ('CRYST line is too short for PDB format',cryst)
    # Unit Cell Dimensions
    a = float(cryst[6:15])
    b = float(cryst[15:24])
    c = float(cryst[24:33])
    # Unit Cell Angles
    alpha = float(cryst[33:40])
    beta  = float(cryst[40:47])
    gamma = float(cryst[47:54])
    # Space Group and Z Value
    spacegroup = cryst[55:66].strip()
    if len(cryst) > 66:
        z_value = cryst[66:70]
    else:
        z_value = ''

    # Format
    cols = ['A', 'B', 'C', 'ALPHA', 'BETA', 'GAMMA', 'SG', 'Z']
    cryst = [a,b,c,alpha,beta,gamma,spacegroup,z_value]
    # Turn into dict with Column Headings as keys
    returncryst = dict(zip(cols,cryst))

    return (None,returncryst)


#! /usr/local/python/python2.7.3-64bit/bin/python

import os, urllib, datetime

from Bamboo.Utils.Constants import LONGESTPEPTIDE, LONGESTOLIGO
from Bamboo.Utils.Constants import aaModsDict, aaCodesDict, ntCodesDict, waterNames, ionNames, solventNames
from Bamboo.Utils.Constants import pdb_headers_to_keep, pdb_footers_to_keep
from Bamboo.Rdkit.Smile import get_smile_from_block, match_smile_to_list

# RDKit Error Suppression
from rdkit import RDLogger
lg = RDLogger.logger()

class MacroMol:
    """Class to handle a MacroMolecule"""

    def __init__(self, pdbin, keepLines=-1, headerOnly=-1):
        """Initialise pdbfile <pdbin> object"""

        self.pdbin = pdbin
        self.pdbid = ''

        if os.path.exists(pdbin) and ('/' in pdbin or '\\' in pdbin or '.pdb' in pdbin):
            self.intype = 'file'
            self.pdbsource = pdbin
            self._setState(1,0)
        elif len(pdbin)==4:
            self.intype = 'web'
            self.pdbsource = 'http://www.rcsb.org/pdb/files/{!s}.pdb'.format(pdbin)
            self._setState(1,0)
        elif 'ATOM' in pdbin or 'HETATM' in pdbin:
            self.intype = 'block'
            self.pdbsource = 'block'
            self._setState(1,0)
        else:
            raise TypeError("Invalid PDBIN given. Must be a '.pdb' file or a valid PDB code: \n\n{!s}".format(pdbin))
        # Update state
        self._setState(keepLines,headerOnly)

        # Read & Parse File
        self._readPdb()
        self._parsePdb()
        self._wipeLines()

    def __iter__(self):
        return iter(self.getChains())

    def __getitem__(self, id):
        if isinstance(id,int):
            return self.getChains()[id]
        if isinstance(id,str):
            return self.getChain(id)

    def _setState(self, keepLines, headerOnly):
        if keepLines>=0: self.keepLines=keepLines
        if headerOnly>=0: self.headerOnly=headerOnly

    def _readPdb(self):
        """Get and Read PDB information"""
        if self.intype=='file':
            self.pdbLines = open(self.pdbsource,'r').readlines()
        elif self.intype=='web':
            self.pdbLines = urllib.urlopen(self.pdbsource).readlines()
        elif self.intype=='block':
            self.pdbLines = self.pdbin.replace('\n','\nZ\rX\rY\r').strip('Z\rX\rY\r').split('Z\rX\rY\r')
            self.pdbin = ''

        if len(self.pdbLines) > 100000:
            raise MacroMolError('LARGE PDB FILE: {!s}'.format(self.pdbsource))

    def _parsePdb(self):
        """Process PDB Information"""
        self.chainList = []
        self.residueList = []
        self.authorLines = []
        self.sourceLines = []
        self.titleLines = []
        self.depositionDate = datetime.datetime.today()
        self.organism = ''
        self.experiment = ''
        # Record Lines that can be easily echoed to output
        self.easy_header = []
        self.easy_footer = []
        # Initialise current residue values
        currChainID = None
        currChain = None
        currResnum = None
        currRes = None
        # Build the lists by going through lines
        for lineNum,line in enumerate(self.pdbLines):
            fields = line.split()
            # Record easy-header and easy-footer
            if [flag for flag in pdb_headers_to_keep if line.startswith(flag)]:
                self.easy_header.append(line)
            if [flag for flag in pdb_footers_to_keep if line.startswith(flag)]:
                self.easy_footer.append(line)
            # Parse meta-data
            if line[:4]!='ATOM' and line[:6]!='HETATM':
                # Source
                if line[:6]=='SOURCE':
                    self.sourceLines.append(line)
                    if 'ORGANISM_COMMON' in line:
                        self.organism = line.split(':')[1].strip()
                # Authors
                if line[:6]=='AUTHOR':
                    self.authorLines.append(line)
                # Deposition Date - on first line of PDB entry
                if line[:6]=='HEADER':
                    try:
                        self.depositionDate = datetime.datetime.strptime(line[50:59],'%d-%b-%y')
                    except ValueError:
                        self.depositionDate = None
                # Release date - on REVDAT record, the one flagged "0" on column 31
                if line[:6]=='REVDAT':
                    if len(line)>31 and line[31]=='0':
                        self.releaseDate = datetime.datetime.strptime(line[13:22],'%d-%b-%y')
                    else:
                        self.releaseDate = None
                # Title
                if line[:5]=='TITLE':
                    self.titleLines.append(line)
                # Experiment type
                if line[12:27]=='EXPERIMENT TYPE':
                    self.experiment = line[45:50].lower().replace('-','')
                # (No other metadata just yet)
                continue
            # Parse atoms
            else:
                # Bail if headerOnly
                if self.headerOnly:
                    break
                # Parse res stats
                resname = line[17:20].strip()
                chainID = line[21].strip()
                try:
                    resnum = int(line[22:26])
                except ValueError:
                    resnum = 0
                inscode = line[26].strip()
                # New chain as needed
                if not chainID:
                    chainID = ' '
                if chainID != currChainID:
                    if isinstance(currChain,Chain) and isinstance(currRes,Residue):
                        # New Chain! - add the last residue to the old chain before moving on...
                        currChain.addRes(currRes)
                        currRes = None
                        currResnum = None
                    # Create a new chain
                    currChainID=chainID
                    if not self.chainList:
                        currChain = Chain(chainID)
                        self.chainList.append(currChain)
                    elif chainID not in self.getChainIDs():
                        currChain = Chain(chainID)
                        self.chainList.append(currChain)
                    else:
                        currChain = self.getChain(chainID)
                # New residue if needed (incl. if new chain)
                if (resnum != currResnum) or (inscode != currInscode):
                    if isinstance(currRes,Residue):
                        # New Residue - add the old one to the chain
                        currChain.addRes(currRes)
                    # Change resnum and create a new residue
                    currResnum = resnum
                    currInscode = inscode
                    currRes = Residue(resname,chainID,resnum,inscode)
                    self.residueList.append(currRes)
                # Add line to residue
                currRes.addLine(line,lineNum)
        # At the end of the file, add the LAST Residue to the chain
        if isinstance(currChain,Chain) and isinstance(currRes,Residue):
            currChain.addRes(currRes)
        # Finish off
        for chn in self.chainList:
            chn.updateTypes()

    def _wipeLines(self):
        """Deletes pdb text if required"""
        if not self.keepLines:
            print('Wiping...')
            self.pdbLines = []

    def _forceFullRead(self):
        """Force full parsing"""
        # Don't proceed if already in this loop
        if '_forcing' in self.__dict__ and self._forcing: return
        # Get the file lines if necessary, but then restore the state
        self._forcing=1
        if not self.chainList:
            state=(self.keepLines,self.headerOnly)
            self._setState(-1,0)
            self._readPdb()
            self._parsePdb()
            self._wipeLines()
            self._setState(*state)
        self._forcing=0

    def getChains(self):
        if not self.chainList:
            self._forceFullRead()
        return self.chainList

    def getChain(self, chainID):
        """Returns <Chain> Object for <chainID>"""
        subset = [c for c in self.getChains() if c.chainID==chainID]
        if len(subset)==1:
            return subset[0]
        else:
            raise Exception('More than one chain found with chainid: {!s}'.format(chainID))

    def getChainIDs(self):
        return [c.chainID for c in self.getChains()]

    def getResidues(self):
        return [res for chn in self.getChains() for res in chn.residues]

    def getResidue(self, resid):
        residues = [res for res in self.getResidues() if res.get_res_id()==resid]
        if len(residues)>1: raise MacroMolError('Too Many Residues Found!')
        if len(residues)<1: raise KeyError('No Residue Found!')
        return residues[0]

    def getAminos(self):
        return [res for chn in self.getChains() for res in chn.aminos]

    def getWaters(self):
        return [res for chn in self.getChains() for res in chn.waters]

    def getIons(self):
        return [res for chn in self.getChains() for res in chn.ions]

    def getSolvent(self):
        return [res for chn in self.getChains() for res in chn.solvent]

    def getUnknowns(self):
        return [res for chn in self.getChains() for res in chn.unknowns]

    def getAtoms(self):
        return [atm for res in self.getResidues() for atm in res.atoms]

    def printChains(self):
        """Prints Chain IDs"""
        for ch in self.getChains():
            print('{!s}: {!s}'.format(ch.chainID, ch.types))

    def getType(self):
        """Returns types of molecule present"""
        typesList = list(set([cTyp for c in self.getChains() for cTyp in ['AMINO','NUCLEO'] if cTyp in c.types]))
        return ''.join(typesList)

    def isProtein(self):
        for chn in self.getChains():
            if chn.isProtein():
                return True
        return False

    def isProteinNucleoComplex(self):
        if ('AMINO' in self.getType()) and ('NUCLEO' in self.getType()):
            return True
        else:
            return False

    def getTitle(self):
        return ''.join([l[10:].strip() for l in self.titleLines])

    def getAuthors(self):
        return ''.join([l[10:].strip() for l in self.authorLines])

class Chain:
    """Class to handle PDB chains"""

    def __init__(self, chainID=''):
        self.chainID = chainID.strip()
        self.residues = []
        self.aminos = []
        self.waters = []
        self.ions = []
        self.solvent = []
        self.unknowns = []
        self.types = []

    def __iter__(self):
        return iter(self.residues)

    def __getitem__(self, id):
        if isinstance(id,int):
            return self.residues[id]

    def __str__(self):
        str = 'Chain {!s}:\n\tNRes = {:d} ({:d} to {:d}), type {!s}'.format(self.chainID, len(self.residues), self.firstResnum(), self.lastResnum(), '-'.join([t for t in self.types]))
        str += '\n\tSequence: {!s}'.format(self.getSequence())
        if len(self.residues)<5:
            str += ' ==> '+''.join([n for n in self.getNames()])
        str += '\n\tAminos: {!s}'.format(len(self.aminos))
        str += '\n\tWaters: {!s}'.format(len(self.waters))
        str += '\n\tIons: {!s}'.format(len(self.ions))
        return str

    def __len__(self):
        return len(self.residues)

    def printInfo(self):
        print(str(self))

    def addRes(self, residue):
        # Check the residue is up to date
        residue.assignTypeAndCode()
        self.residues.append(residue)
        # Now sort it
        if residue.code == 'i':
            self.ions.append(residue)
        elif residue.code == 'o':
            self.waters.append(residue)
        elif 'SOLVENT' in residue.types:
            self.solvent.append(residue)
        elif 'UNKNOWN' in residue.types:
            self.unknowns.append(residue)
        else:
            self.aminos.append(residue)
        [self.types.append(t) for t in residue.types if t not in self.types]

    def getRes(self, resnum):
        subset = [r for r in self.residues if r.resnum == resnum]
        if len(subset)==1:
            return subset[0]
        else:
            raise Exception('More than one residue found with residue number {!s}'.format(resnum))

    def firstResidue(self):
        if self.residues:
            return self.residues[0]

    def lastResidue(self):
        if self.residues:
            return self.residues[-1]

    def firstResnum(self):
        return min(self.getNumbers())

    def lastResnum(self):
        return max(self.getNumbers())

    def getNumbers(self):
        return [r.resnum for r in self.residues]

    def getNames(self):
        return [r.resname for r in self.residues]

    def getSequence(self):
        sequence = '-'
        for res in self.aminos:
            if len(res.code)==1:
                sequence += res.code
            else:
                if (sequence[-1]!='-'):
                    sequence += '-'
                sequence += '('+res.code+')-'
        return sequence.strip('-')

    def getLongSequence(self):
        return '-'.join([r.resname for r in self.aminos])

    def updateTypes(self):
        [self.types.append(t) for r in self.residues for t in r.types if t not in self.types]

    def isProtein(self):
        if self.getAACount() > LONGESTPEPTIDE:
            return True
        else:
            return False

    def isPeptide(self):
        aaCount = self.getAACount()
        if aaCount > 0 and aaCount <= LONGESTPEPTIDE:
            return True
        else:
            return False

    def getAACount(self):
        return len([r for r in self.residues if (('AMINO' in r.types) or ('AMINO-MOD' in r.types))])

class Residue:
    """Class to handle PDB residues"""

    def __init__(self, resname, chain, resnum, inscode):
        self.resname = resname.strip()
        self.chain = chain.strip()
        self.resnum = resnum
        self.inscode = inscode
        self.firstLine = -1
        self.lastLine = -1
        self.atoms = []
        self.smile = ''
        self.assignTypeAndCode()

    def __iter__(self):
        return iter(self.atoms)

    def __len__(self):
        return len(self.atoms)

    def getPdbString(self):
        return ''.join([atm.line for atm in self.atoms])

    def printLines(self):
        print(self.getPdbString())

    def addLine(self, pdbLine, lineNum):
        pdbLine = pdbLine.strip('\n')
        # Make the right length for PDB format
        if len(pdbLine)<80:
            pdbLine = pdbLine + ' '*(80-len(pdbLine))
        # Add new atom
        self.atoms.append(Atom(pdbLine))
        # Check for line numbers
        if self.firstLine<0:
            self.firstLine=lineNum
        self.lastLine=lineNum

    def assignTypeAndCode(self):
        types = []
        altCodes = []
        # Amino Acids
        if self.resname in aaCodesDict:
            types.append('AMINO')
            code = aaCodesDict[self.resname]
        # Nucleotides
        elif self.resname in ntCodesDict:
            types.append('NUCLEO')
            code = ntCodesDict[self.resname]
        # Waters
        elif self.isWater():
            types.append('WATER')
            code = 'o'
        # Ions
        elif self.isIon():
            types.append('ION')
            code = 'i'
        # Solvent/Buffer
        elif self.isSolventOrBuffer():
            types.append('SOLVENT')
            code = self.resname
            self.getSmiles()
        # Modified amino acids or ligands (or unknown)
        else:
            code = self.resname
            self.getSmiles()
            if self.smile:
                matches = match_smile_to_list(self.smile, aaModsDict.keys(), assignbonds=True)
                if matches:
                    types.append('AMINO-MOD')
                    [altCodes.extend(aaModsDict[s]) for s in matches]
            if not (types and code):
                types.append('UNKNOWN')

        self.types = types
        self.code = code
        if altCodes:
            self.altCodes = altCodes
        else:
            self.altCodes = [code]

    def isWater(self):
        if self.resname.upper() in waterNames:
            return True
        else:
            return False

    def isIon(self):
        """Test if the residue is a non-water ion"""

        # Count number of atoms (accounting for differing alternate locations)
        if len([at for at in self.atoms if at.altloc in ['', 'A']])!=1:
            return False
        elif self.isWater():
            return False
        else:
            return True

    def isSolventOrBuffer(self):
        """Test is the residue is a (common) solvent or buffer molecule"""

        # Count number of atoms (must be greater than one)
        if len([at for at in self.atoms if at.altloc in ['', 'A']])==1:
            return False
        elif self.resname.upper() in solventNames:
            return True
        else:
            return False

    def getSmiles(self):
        if self.atoms:
# THIS IS A COMMON PLACE FOR RDKIT ERROR MESSAGES
#            print('Making Noise.')
            lg.setLevel(RDLogger.CRITICAL)
            self.smile = get_smile_from_block(self.getPdbString())
            lg.setLevel(RDLogger.ERROR)
#            print("I'll be quiet now...")

    def get_res_id(self):
        """Get a tuple that identifies the residue (used as a key in other functions to refer to the residue)"""
        return (self.resname, self.chain, self.resnum, self.inscode)

    def get_non_h_atoms(self):
        """Get non-hydrogen atoms of the residue"""
        if not self.atoms:
            raise ValueError('No atoms in residue! Cannot calculate average occupancy!')
        return [atm for atm in self.atoms if atm.element != 'H']

    def get_centroid(self):
        """Return the centre of the residue, giving equal weight to each of the atoms (no H)"""
        non_h = self.get_non_h_atoms()
        x_cent = sum([atm.x for atm in non_h])/len(non_h)
        y_cent = sum([atm.y for atm in non_h])/len(non_h)
        z_cent = sum([atm.z for atm in non_h])/len(non_h)
        return (x_cent, y_cent, z_cent)

    def get_mean_occupancy(self):
        """Get the average occupancy for atoms in the residue (excluding hydrogens)"""
        if not self.atoms:
            raise ValueError('No atoms in residue! Cannot calculate average occupancy!')
        non_h = self.get_non_h_atoms()
        mean_occupancy = sum([atm.occupancy for atm in non_h])/len(non_h)
        return mean_occupancy

    def get_mean_bfactor(self):
        """Get the average bfactor for atoms in the residue (excluding hydrogens)"""
        if not self.atoms:
            raise ValueError('No atoms in residue! Cannot calculate average bfactor!')
        non_h = self.get_non_h_atoms()
        mean_bfactor = sum([atm.bfactor for atm in non_h])/len(non_h)
        return mean_bfactor

    def update_inscode(self, inscode):
        """Change the inscode of atom and update atom line"""
        if not self.atoms:
            raise ValueError('No atoms in residue! Cannot change insertion code!')
        new_ins = str(inscode.upper())
        assert len(new_ins)<2, 'New inscode too long! {!s}'.format(new_ins)
        for atm in self.atoms:
            # Set the inscode
            atm.inscode = new_ins
            # Change the pdbline
            newline = list(atm.line)
            newline[26] = new_ins
            atm.line = ''.join(newline)
            # Sanity check line length
            assert len(atm.line)==81, 'ATOM LINE IS NOT THE RIGHT LENGTH FOR PDB FORMAT! {!s}'.format(atm.line)

    def update_atom_numbers(self,atom_numbering):
        """Update the atom numbers to those in atom_numbering - Requires a dict with atom names as the keys for the new numbers"""
        if not self.atoms:
            raise ValueError('No atoms in residue! Cannot change atom numbers!')
        for atm in self.atoms:
            # Set the new serial num
            atm.serialnum = atom_numbering[atm.atomname]
            # Change the pdbline
            newline = list(atm.line)
            newline[6:11] = '{:>5}'.format(str(atom_numbering[atm.atomname]))
            atm.line = ''.join(newline)
            # Sanity check line length
            assert len(atm.line)==81, 'ATOM LINE IS NOT THE RIGHT LENGTH FOR PDB FORMAT! {!s}'.format(atm.line)

class Atom:
    """Atom Class"""

    def __init__(self, line):

        # Process input line and populate record
        line = line.strip('\n').upper()
        # Make the right length
        if len(line) < 80:
            line = line+' '*(80-len(line))
        # Save the line
        self.line = line+'\n'
        # Check the ATOM or HETATM records
        try:
            self.recordtype = line[0:6].strip()
            self.serialnum  = line[6:11].strip()
            self.atomname   = line[12:16].strip()
            self.altloc     = line[16].strip()
            self.resname    = line[17:20].strip()
            self.chnname    = line[21].strip()
            self.resnum     = line[22:26].strip()
            self.inscode    = line[26].strip()
            self.x          = float(line[30:38])
            self.y          = float(line[38:46])
            self.z          = float(line[46:54])
            self.occupancy  = float(line[54:60])
            self.bfactor    = float(line[60:66])
            self.segment    = line[72:76].strip()
            self.element    = line[76:78].strip()
            self.charge     = line[78:80].strip()
        except ValueError as Err:
            raise MacroMolError('Error in PDB File Atom Coordinates. Maybe Wrong File Format?\nOriginal Error is: {!s}'.format(Err))


import os, urllib, datetime

from bamboo.constants import LONGEST_PEPTIDE, LONGEST_OLIGO, \
                                    AA_MODIFICATIONS_DICT, AA_3_TO_1_CODES_DICT, NUCLEOTIDE_DICT, \
                                    WATER_NAMES, ION_NAMES, SOLVENT_NAMES, \
                                    PDB_HEADERS_TO_KEEP, PDB_FOOTERS_TO_KEEP
from bamboo.rdkit_utils.smile import get_smile_from_block, match_smile_to_list

# RDKit Error Suppression
from rdkit import RDLogger
lg = RDLogger.logger()

class MacroMolError(Exception):
    pass

class MacroMol:
    """Class to handle a MacroMolecule"""

    def __init__(self, pdbin, keep_lines=-1, header_only=-1):
        """Initialise pdbfile <pdbin> object"""

        self.pdbin = pdbin
        self.pdbid = ''

        if os.path.exists(pdbin) and ('/' in pdbin or '\\' in pdbin or '.pdb' in pdbin):
            self.intype = 'file'
            self.pdbsource = pdbin
            self._set_state(1,0)
        elif len(pdbin)==4:
            self.intype = 'web'
            self.pdbsource = 'http://www.rcsb.org/pdb/files/{!s}.pdb'.format(pdbin)
            self._set_state(1,0)
        elif 'ATOM' in pdbin or 'HETATM' in pdbin:
            self.intype = 'block'
            self.pdbsource = 'block'
            self._set_state(1,0)
        else:
            raise TypeError("Invalid PDBIN given. Must be a '.pdb' file or a valid PDB code: \n\n{!s}".format(pdbin))
        # Update state
        self._set_state(keep_lines,header_only)

        # Read & Parse File
        self._read_pdb()
        self._parse_pdb()
        self._wipe_lines()

    def __iter__(self):
        return iter(self.get_chains())

    def __getitem__(self, id):
        if isinstance(id,int):
            return self.get_chains()[id]
        if isinstance(id,str):
            return self.get_chain(id)

    def _set_state(self, keep_lines, header_only):
        if keep_lines>=0: self.keep_lines=keep_lines
        if header_only>=0: self.header_only=header_only

    def _read_pdb(self):
        """Get and Read PDB information"""
        if self.intype=='file':
            self.pdb_lines = open(self.pdbsource,'r').readlines()
        elif self.intype=='web':
            self.pdb_lines = urllib.urlopen(self.pdbsource).readlines()
        elif self.intype=='block':
            self.pdb_lines = self.pdbin.replace('\n','\nZ\rX\rY\r').strip('Z\rX\rY\r').split('Z\rX\rY\r')
            self.pdbin = ''

        if len(self.pdb_lines) > 100000:
            raise MacroMolError('LARGE PDB FILE: {!s}'.format(self.pdbsource))

    def _parse_pdb(self):
        """Process PDB Information"""
        self.chain_list = []
        self.residue_list = []
        self.author_lines = []
        self.source_lines = []
        self.title_lines = []
        self.deposition_date = datetime.datetime.today()
        self.organism = ''
        self.experiment = ''
        # Record Lines that can be easily echoed to output
        self.easy_header = []
        self.easy_footer = []
        # Initialise current residue values
        curr_chain_id = None
        curr_chain = None
        curr_resnum = None
        curr_res = None
        # Build the lists by going through lines
        for line_num,line in enumerate(self.pdb_lines):
            fields = line.split()
            # Record easy-header and easy-footer
            if [flag for flag in PDB_HEADERS_TO_KEEP if line.startswith(flag)]:
                self.easy_header.append(line)
            if [flag for flag in PDB_FOOTERS_TO_KEEP if line.startswith(flag)]:
                self.easy_footer.append(line)
            # Parse meta-data
            if line[:4]!='ATOM' and line[:6]!='HETATM':
                # Source
                if line[:6]=='SOURCE':
                    self.source_lines.append(line)
                    if 'ORGANISM_COMMON' in line:
                        self.organism = line.split(':')[1].strip()
                # Authors
                if line[:6]=='AUTHOR':
                    self.author_lines.append(line)
                # Deposition Date - on first line of PDB entry
                if line[:6]=='HEADER':
                    try:
                        self.deposition_date = datetime.datetime.strptime(line[50:59],'%d-%b-%y')
                    except ValueError:
                        self.deposition_date = None
                # Release date - on REVDAT record, the one flagged "0" on column 31
                if line[:6]=='REVDAT':
                    if len(line)>31 and line[31]=='0':
                        self.release_date = datetime.datetime.strptime(line[13:22],'%d-%b-%y')
                    else:
                        self.release_date = None
                # Title
                if line[:5]=='TITLE':
                    self.title_lines.append(line)
                # Experiment type
                if line[12:27]=='EXPERIMENT TYPE':
                    self.experiment = line[45:50].lower().replace('-','')
                # (No other metadata just yet)
                continue
            # Parse atoms
            else:
                # Bail if headerOnly
                if self.header_only:
                    break
                # Parse res stats
                resname = line[17:20].strip()
                chain_id = line[21].strip()
                try:
                    resnum = int(line[22:26])
                except ValueError:
                    resnum = 0
                inscode = line[26].strip()
                # New chain as needed
                if not chain_id:
                    chain_id = ' '
                if chain_id != curr_chain_id:
                    if isinstance(curr_chain,Chain) and isinstance(curr_res,Residue):
                        # New Chain! - add the last residue to the old chain before moving on...
                        curr_chain.add_res(curr_res)
                        curr_res = None
                        curr_resnum = None
                    # Create a new chain
                    curr_chain_id=chain_id
                    if not self.chain_list:
                        curr_chain = Chain(chain_id)
                        self.chain_list.append(curr_chain)
                    elif chain_id not in self.get_chain_ids():
                        curr_chain = Chain(chain_id)
                        self.chain_list.append(curr_chain)
                    else:
                        curr_chain = self.get_chain(chain_id)
                # New residue if needed (incl. if new chain)
                if (resnum != curr_resnum) or (inscode != curr_inscode):
                    if isinstance(curr_res,Residue):
                        # New Residue - add the old one to the chain
                        curr_chain.add_res(curr_res)
                    # Change resnum and create a new residue
                    curr_resnum = resnum
                    curr_inscode = inscode
                    curr_res = Residue(resname,chain_id,resnum,inscode)
                    self.residue_list.append(curr_res)
                # Add line to residue
                curr_res.add_line(line,line_num)
        # At the end of the file, add the LAST Residue to the chain
        if isinstance(curr_chain,Chain) and isinstance(curr_res,Residue):
            curr_chain.add_res(curr_res)
        # Finish off
        for chn in self.chain_list:
            chn.update_types()

    def _wipe_lines(self):
        """Deletes pdb text if required"""
        if not self.keep_lines:
            print('Wiping...')
            self.pdb_lines = []

    def _force_full_read(self):
        """Force full parsing"""
        # Don't proceed if already in this loop
        if '_forcing' in self.__dict__ and self._forcing: return
        # Get the file lines if necessary, but then restore the state
        self._forcing=1
        if not self.chain_list:
            state=(self.keep_lines,self.header_only)
            self._set_state(-1,0)
            self._read_pdb()
            self._parse_pdb()
            self._wipe_lines()
            self._set_state(*state)
        self._forcing=0

    def get_chains(self):
        if not self.chain_list:
            self._force_full_read()
        return self.chain_list

    def get_chain(self, chain_id):
        """Returns <Chain> Object for <chainID>"""
        subset = [c for c in self.get_chains() if c.chain_id==chain_id]
        if len(subset)==1:
            return subset[0]
        else:
            raise Exception('More than one chain found with chainid: {!s}'.format(chain_id))

    def get_chain_ids(self):
        return [c.chain_id for c in self.get_chains()]

    def get_residues(self):
        return [res for chn in self.get_chains() for res in chn.residues]

    def get_residue(self, resid):
        residues = [res for res in self.get_residues() if res.get_res_id()==resid]
        if len(residues)>1: raise MacroMolError('Too Many Residues Found!')
        if len(residues)<1: raise KeyError('No Residue Found!')
        return residues[0]

    def get_aminos(self):
        return [res for chn in self.get_chains() for res in chn.aminos]

    def get_waters(self):
        return [res for chn in self.get_chains() for res in chn.waters]

    def get_ions(self):
        return [res for chn in self.get_chains() for res in chn.ions]

    def get_solvent(self):
        return [res for chn in self.get_chains() for res in chn.solvent]

    def get_unknowns(self):
        return [res for chn in self.get_chains() for res in chn.unknowns]

    def get_atoms(self):
        return [atm for res in self.get_residues() for atm in res.atoms]

    def print_chains(self):
        """Prints Chain IDs"""
        for ch in self.get_chains():
            print('{!s}: {!s}'.format(ch.chain_id, ch.types))

    def get_type(self):
        """Returns types of molecule present"""
        types_list = list(set([c_typ for c in self.get_chains() for c_typ in ['AMINO','NUCLEO'] if c_typ in c.types]))
        return ''.join(types_list)

    def is_protein(self):
        for chn in self.get_chains():
            if chn.is_protein():
                return True
        return False

    def is_protein_nucleo_complex(self):
        if ('AMINO' in self.get_type()) and ('NUCLEO' in self.get_type()):
            return True
        else:
            return False

    def get_title(self):
        return ''.join([l[10:].strip() for l in self.title_lines])

    def get_authors(self):
        return ''.join([l[10:].strip() for l in self.author_lines])

class Chain:
    """Class to handle PDB chains"""

    def __init__(self, chain_id=''):
        self.chain_id = chain_id.strip()
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
        str = 'Chain {!s}:\n\tNRes = {:d} ({:d} to {:d}), type {!s}'.format(self.chain_id, len(self.residues), self.first_resnum(), self.last_resnum(), '-'.join([t for t in self.types]))
        str += '\n\tSequence: {!s}'.format(self.get_sequence())
        if len(self.residues)<5:
            str += ' ==> '+''.join([n for n in self.get_names()])
        str += '\n\tAminos: {!s}'.format(len(self.aminos))
        str += '\n\tWaters: {!s}'.format(len(self.waters))
        str += '\n\tIons: {!s}'.format(len(self.ions))
        return str

    def __len__(self):
        return len(self.residues)

    def print_info(self):
        print(str(self))

    def add_res(self, residue):
        # Check the residue is up to date
        residue.assign_type_and_code()
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

    def get_res(self, resnum):
        subset = [r for r in self.residues if r.resnum == resnum]
        if len(subset)==1:
            return subset[0]
        else:
            raise Exception('More than one residue found with residue number {!s}'.format(resnum))

    def first_residue(self):
        if self.residues:
            return self.residues[0]

    def last_residue(self):
        if self.residues:
            return self.residues[-1]

    def first_resnum(self):
        return min(self.get_numbers())

    def last_resnum(self):
        return max(self.get_numbers())

    def get_numbers(self):
        return [r.resnum for r in self.residues]

    def get_names(self):
        return [r.resname for r in self.residues]

    def get_sequence(self):
        sequence = '-'
        for res in self.aminos:
            if len(res.code)==1:
                sequence += res.code
            else:
                if (sequence[-1]!='-'):
                    sequence += '-'
                sequence += '('+res.code+')-'
        return sequence.strip('-')

    def get_long_sequence(self):
        return '-'.join([r.resname for r in self.aminos])

    def update_types(self):
        [self.types.append(t) for r in self.residues for t in r.types if t not in self.types]

    def is_protein(self):
        if self.get_aa_count() > LONGEST_PEPTIDE:
            return True
        else:
            return False

    def is_peptide(self):
        aa_count = self.get_aa_count()
        if aa_count > 0 and aa_count <= LONGEST_PEPTIDE:
            return True
        else:
            return False

    def get_aa_count(self):
        return len([r for r in self.residues if (('AMINO' in r.types) or ('AMINO-MOD' in r.types))])

class Residue:
    """Class to handle PDB residues"""

    def __init__(self, resname, chain, resnum, inscode):
        self.resname = resname.strip()
        self.chain = chain.strip()
        self.resnum = resnum
        self.inscode = inscode
        self.first_line = -1
        self.last_line = -1
        self.atoms = []
        self.smile = ''
        self.assign_type_and_code()

    def __iter__(self):
        return iter(self.atoms)

    def __len__(self):
        return len(self.atoms)

    def get_pdb_string(self):
        return ''.join([atm.line for atm in self.atoms])

    def print_lines(self):
        print(self.get_pdb_string())

    def add_line(self, pdb_line, line_num):
        pdb_line = pdb_line.strip('\n')
        # Make the right length for PDB format
        if len(pdb_line)<80:
            pdb_line = pdb_line + ' '*(80-len(pdb_line))
        # Add new atom
        self.atoms.append(Atom(pdb_line))
        # Check for line numbers
        if self.first_line<0:
            self.first_line=line_num
        self.last_line=line_num

    def assign_type_and_code(self):
        types = []
        alt_codes = []
        # Amino Acids
        if self.resname in AA_3_TO_1_CODES_DICT:
            types.append('AMINO')
            code = AA_3_TO_1_CODES_DICT[self.resname]
        # Nucleotides
        elif self.resname in NUCLEOTIDE_DICT:
            types.append('NUCLEO')
            code = NUCLEOTIDE_DICT[self.resname]
        # Waters
        elif self.is_water():
            types.append('WATER')
            code = 'o'
        # Ions
        elif self.is_ion():
            types.append('ION')
            code = 'i'
        # Solvent/Buffer
        elif self.is_solvent_or_buffer():
            types.append('SOLVENT')
            code = self.resname
            self.get_smiles()
        # Modified amino acids or ligands (or unknown)
        else:
            code = self.resname
            self.get_smiles()
            if self.smile:
                matches = match_smile_to_list(self.smile, AA_MODIFICATIONS_DICT.keys(), assignbonds=True)
                if matches:
                    types.append('AMINO-MOD')
                    [alt_codes.extend(AA_MODIFICATIONS_DICT[s]) for s in matches]
            if not (types and code):
                types.append('UNKNOWN')

        self.types = types
        self.code = code
        if alt_codes:
            self.alt_codes = alt_codes
        else:
            self.alt_codes = [code]

    def is_water(self):
        if self.resname.upper() in WATER_NAMES:
            return True
        else:
            return False

    def is_ion(self):
        """Test if the residue is a non-water ion"""

        # Count number of atoms (accounting for differing alternate locations)
        if len([at for at in self.atoms if at.altloc in ['', 'A']])!=1:
            return False
        elif self.is_water():
            return False
        else:
            return True

    def is_solvent_or_buffer(self):
        """Test is the residue is a (common) solvent or buffer molecule"""

        # Count number of atoms (must be greater than one)
        if len([at for at in self.atoms if at.altloc in ['', 'A']])==1:
            return False
        elif self.resname.upper() in SOLVENT_NAMES:
            return True
        else:
            return False

    def get_smiles(self):
        if self.atoms:
            lg.setLevel(RDLogger.CRITICAL)

            altlocs = sorted(set([at.altloc for at in self.atoms if at.altloc]))
            for alt in altlocs:
                atms = [atm for atm in self.atoms if ((not atm.altloc) or (atm.altloc == alt))]
                pdb_string = ''.join([atm.line[0:16]+' '+atm.line[17:] for atm in atms])
                smile = get_smile_from_block(pdb_string)
                if smile:
                    self.smile = smile
                    break
#            self.smile = get_smile_from_block(self.get_pdb_string())
            lg.setLevel(RDLogger.ERROR)

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


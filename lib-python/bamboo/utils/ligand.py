from rdkit import Chem

from bamboo.common.file import FileObj
from bamboo.macro.molecule import MacroMol
from bamboo.rdkit_utils.mol import check_pdb_readable
from bamboo.constants import ATOM_NUMBER_DICT

################################################################################################################
### CLASSES                                                                                                    #
################################################################################################################

class LigFile(object):
    """Class for summarising a ligand PDBFile"""
    def __init__(self, ligpath, parse=False):
        self.file = FileObj(ligpath)
        # Convert to mol, checking that's it's valid chemically
        self.mol = check_pdb_readable(ligpath)
        # Create MacroMol object
        if parse:
            self.macromol = MacroMol(ligpath)
        else:
            self.macromol = None

    def get_total_atom_count(self):
        return self.mol.GetNumAtoms()

    def get_heavy_atom_count(self):
        return self.mol.GetNumHeavyAtoms()

    def get_light_atom_count(self):
        return self.mol.GetNumAtoms() - self.mol.GetNumHeavyAtoms()

    def get_bond_count(self):
        return self.mol.GetNumBonds()

    def get_conformer_count(self):
        return self.mol.GetNumConformers()

    def get_atom_counts(self):
        """Gets a summary of the counts of each atom in the ligand"""

        atomcounts = {}
        ats = self.mol.GetAtoms()
        for a in ats:
            atomic_symbol = ATOM_NUMBER_DICT[a.GetAtomicNum()]
            if atomic_symbol not in atomcounts:
                atomcounts[atomic_symbol] = 1
            else:
                atomcounts[atomic_symbol] += 1

        return atomcounts

    def get_atom_counts_by_fragment(self):
        """Gets atom counts by residue"""

        fragments = {}

        inscodes = sorted(list(set([at.GetMonomerInfo().GetInsertionCode() for at in self.mol.GetAtoms()])))

        for ins in inscodes:
            ats = [at for at in self.mol.GetAtoms() if at.GetMonomerInfo().GetInsertionCode()==ins]

            atomcounts = {}
            for a in ats:
                atomic_symbol = ATOM_NUMBER_DICT[a.GetAtomicNum()]
                if atomic_symbol not in atomcounts:
                    atomcounts[atomic_symbol] = 1
                else:
                    atomcounts[atomic_symbol] += 1

            fragments[ins] = atomcounts

        return fragments

    def get_total_atom_count_by_fragment(self):

        fragments = {}

        inscodes = sorted(list(set([at.GetMonomerInfo().GetInsertionCode() for at in self.mol.GetAtoms()])))

        for ins in inscodes:
            ats = [at for at in self.mol.GetAtoms() if at.GetMonomerInfo().GetInsertionCode()==ins]
            fragments[ins] = len(ats)

        return fragments


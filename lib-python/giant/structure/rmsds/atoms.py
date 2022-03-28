import giant.logs as lg
logger = lg.getLogger(__name__)

import giant.common.geometry as gm


class CalculatePairedAtomSetsRMSD(object):
    """
    Calculate the RMSD between two sets of equivalent atoms (i.e. atoms with the same names)
    - sort                   - sort atoms prior to rmsd calculation (so that equivalent atoms are compared)
    - truncate_to_common_set - only calculate rmsd over the common set of atoms, return None if different atoms are provided.
    Raises error if atom names are present more than once in each list (e.g. alternate conformers of atoms)
    """

    def __init__(self,
        truncate_to_common_set = True,
        remove_hydrogens = True,
        ):

        self.truncate_to_common_set = bool(
            truncate_to_common_set
            )

        self.remove_hydrogens = bool(
            remove_hydrogens
            )

    def __call__(self,
        atoms_1,
        atoms_2,
        ):

        atoms_1 = self.filter_atoms(atoms_1)
        atoms_2 = self.filter_atoms(atoms_2)

        atom_dict_1 = atoms_1.build_dict()
        atom_dict_2 = atoms_2.build_dict()

        atom_keys_1 = set(atom_dict_1.keys())
        atom_keys_2 = set(atom_dict_2.keys())

        if self.truncate_to_common_set is False:

            if atom_keys_1.symmetric_difference(atom_keys_2):
                return None

        atom_keys = list(atom_keys_1.intersection(atom_keys_2)) # ordering does not matter

        rmsd = gm.rmsd_coordinates(
            points_1 = [atom_dict_1[k].xyz for k in atom_keys],
            points_2 = [atom_dict_2[k].xyz for k in atom_keys],
            )

        return rmsd

    def filter_atoms(self, atoms):

        if self.remove_hydrogens is True:

            atoms = atoms.select(atoms.extract_element() != ' H')

        return atoms


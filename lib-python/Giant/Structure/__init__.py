import iotbx.pdb


class residueGroupSummary(object):
    def __init__(self, residue_group):
        self.residue_group = residue_group

class atomGroupSummary(object):

    def __init__(self, atom_group):

        self.atom_group = atom_group

        # Meta about atom group
        self.residue_class = iotbx.pdb.common_residue_names_get_class(atom_group.resname)

        # Occupancy summary
        occ_stats = atom_group.atoms().extract_occ().min_max_mean()
        self.occ_min = occ_stats.min
        self.occ_max = occ_stats.max
        self.occ_mean = occ_stats.mean

        # B-factor summary
        b_stats = atom_group.atoms().extract_b().min_max_mean()
        self.b_min = b_stats.min
        self.b_max = b_stats.max
        self.b_mean = b_stats.mean

        # Coordinate summary
        xyz_coords = atom_group.atoms().extract_xyz()
        self.xyz_min = xyz_coords.min()
        self.xyz_max = xyz_coords.max()
        self.xyz_mean = xyz_coords.mean()


import numpy
import iotbx.mtz, iotbx.pdb

class CrystalInfo(object):

    def __init__(self,
        crystal_symmetry = None,
        space_group = None,
        unit_cell = None,
        resolution_high = None,
        resolution_low = None,
        r_free = None,
        r_work = None,
        column_labels = None,
        name = None,
        ):

        if (crystal_symmetry is not None):
            assert space_group is None
            assert unit_cell is None
            space_group = crystal_symmetry.space_group()
            unit_cell = crystal_symmetry.unit_cell()

        self.crystal_symmetry = crystal_symmetry
        self.space_group = space_group
        self.unit_cell = unit_cell

        self.resolution_high = resolution_high
        self.resolution_low = resolution_low

        self.r_free = r_free
        self.r_work = r_work

        self.column_labels = column_labels

        self.name = name

    def __str__(self):
        return self.summary()

    def summary(self):

        lines = []

        lines.append(
            'Name: {}'.format(
                str(self.name)
                )
            )

        lines.append(
            'Space group: {}'.format(
                self.space_group.type().lookup_symbol()) 
                if hasattr(self.space_group, "type") 
                else str(self.space_group)
            )

        lines.append(
            'Unit Cell: {}'.format(
                tuple([round(p, 3) for p in self.unit_cell.parameters()])
                if hasattr(self.unit_cell, "parameters")
                else str(self.unit_cell)
                )
            )

        lines.append(
            'Resolution (high): {}'.format(
                round(self.resolution_high, 3)
                if (self.resolution_high is not None)
                else str(resolution_high)
                )
            )
        
        lines.append(
            'Resolution (low): {}'.format(
                round(self.resolution_low, 3)
                if (self.resolution_low is not None)
                else str(self.resolution_low)
                )
            )

        lines.append(
            'R-free: {}'.format(
                round(self.r_free, 3)
                if (self.r_free is not None)
                else str(self.r_free)
                )
            )
        
        lines.append(
            'R-work: {}'.format(
                round(self.r_work,3)
                if (self.r_work is not None)
                else str(self.r_work)
                )
            )
    
        lines.append(
            'Columns: {}'.format(
                ', '.join(self.column_labels)
                if (self.column_labels is not None)
                else str(self.column_labels)
                )
            )

        return '\n'.join(lines)

    @classmethod
    def from_mtz(cls, mtz_file=None, mtz_object=None, name=None):

        assert [mtz_file, mtz_object].count(None)==1,'Provide mtz_file OR mtz_object'

        if (mtz_object is None):
            mtz_object = iotbx.mtz.object(mtz_file)

        # Use first crystal
        crystal = mtz_object.crystals()[0]
        crystal_symmetry = crystal.crystal_symmetry()
        resolution_low, resolution_high = mtz_object.max_min_resolution()
        column_labels = mtz_object.column_labels()

        return cls(
            crystal_symmetry = crystal_symmetry,
            resolution_high = resolution_high,
            resolution_low = resolution_low,
            column_labels = column_labels,
            name = name,
        )

    @classmethod
    def from_pdb(cls, pdb_file=None, pdb_input=None, name=None):

        assert [pdb_file, pdb_input].count(None)==1,'Provide pdb_file OR pdb_input'

        if (pdb_input is None):
            pdb_input = iotbx.pdb.input(pdb_file)

        try:
            info = pdb_input.get_r_rfree_sigma()
        except TypeError:
            info = pdb_input.get_r_rfree_sigma(pdb_file)

        return cls(
            crystal_symmetry = pdb_input.crystal_symmetry(),
            resolution_high = info.high,
            resolution_low = info.low,
            r_free = info.r_free,
            r_work = info.r_work,
            name = name,
        )

    def as_cryst(self):

        cell = self.unit_cell.parameters()
        sg = self.space_group

        return 'CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f}{:>11}{:4}'.format(
            cell[0],
            cell[1],
            cell[2],
            cell[3],
            cell[4],
            cell[5],
            sg.type().lookup_symbol(),
            ' ',
        )


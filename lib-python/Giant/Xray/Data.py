import iotbx.mtz

class crystalSummary(object):

    def __init__(self, mtz_file=None, id=None):
        # Store filename
        self.mtz_file = mtz_file
        self.id = id
        # Create mtz object of reflection data
        refl_data = self.mtz_object()

        # Extract the resolution limits
        self.low_res, self.high_res = refl_data.max_min_resolution()
        self.space_group = refl_data.space_group()

        # Extract unit cell from the crystal
        crystal = refl_data.crystals()[0]
        self.unit_cell = crystal.unit_cell()
        self.symmetry = crystal.crystal_symmetry()

        # Column information
        self.column_labels = refl_data.column_labels()

    def mtz_object(self):
        return iotbx.mtz.object(self.mtz_file)

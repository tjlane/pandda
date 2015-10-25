import numpy
import iotbx.mtz
from Bamboo.Common import Info

class crystalSummary(object):
    def __init__(self, mtz_file=None, mtz_object=None, pdb_file=None, id=None):

        assert [mtz_file, mtz_object].count(None)==1,'Provide mtz_file OR mtz_object'

        # Store filename
        self.mtz_file = mtz_file
        self._mtz_object = mtz_object
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
        if self._mtz_object: return self._mtz_object
        else: return iotbx.mtz.object(self.mtz_file)

class unitCellVariation(Info):
    def __init__(self, unit_cells):
        unit_cell_params = numpy.array([uc.parameters() for uc in unit_cells])
        self.mean  = tuple(unit_cell_params.mean(axis=0))
        self.std   = tuple(unit_cell_params.std(axis=0))
        self.max   = tuple(unit_cell_params.max(axis=0))
        self.min   = tuple(unit_cell_params.min(axis=0))
        self.range = tuple(numpy.array(self.max)-numpy.array(self.min))
        self.std_perc = tuple(100.0*numpy.array(self.std)/numpy.array(self.mean))
        self.max_perc = tuple(100.0*numpy.array(self.min)/numpy.array(self.mean))
        self.min_perc = tuple(100.0*numpy.array(self.max)/numpy.array(self.mean))
        self._initialized = True
    def as_pandas_table(self):
        """Return the variation as a pandas table"""
        import pandas
        pd = pandas.DataFrame(index=['a','b','c','alpha','beta','gamma'])
        pd['min']   = self.min
        pd['mean']  = self.mean
        pd['max']   = self.max
        pd['std']   = self.std
        pd['range'] = self.range
        pd['%std']  = self.std_perc
        pd['%min']  = self.min_perc
        pd['%max']  = self.max_perc
        return pd


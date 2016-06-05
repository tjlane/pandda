import numpy
import iotbx.mtz, iotbx.pdb
from bamboo.common import Info

class CrystalSummary(object):

    def __init__(self, id=None):
        self.id = id
        # File-specific
        self.mtz_file = None
        self._mtz_object = None
        self.pdb_file = None
        self._pdb_input = None
        # Crystal-specific
        self.space_group = None
        self.unit_cell = None
        self.symmetry = None
        # MTZ-Specific
        self.column_labels = None
        self.high_res = None
        self.low_res = None
        # PDB-Specific

    def mtz_object(self):
        if self._mtz_object: return self._mtz_object
        else: return iotbx.mtz.object(self.mtz_file)

    def pdb_input(self):
        if self._pdb_input: return self._pdb_input
        else: return iotbx.pdb.input(self.pdb_file)

    ######################################################
    @classmethod
    def from_mtz(cls, mtz_file=None, mtz_object=None, id=None):
        assert [mtz_file, mtz_object].count(None)==1,'Provide mtz_file OR mtz_object'
        # Initialise class (populate manually)
        new_cls = cls(id=id)
        # Store filename
        new_cls.mtz_file = mtz_file
        new_cls._mtz_object = mtz_object
        # Create mtz object of reflection data
        refl_data = new_cls.mtz_object()
        # Extract the resolution limits
        new_cls.low_res, new_cls.high_res = refl_data.max_min_resolution()
        new_cls.space_group = refl_data.space_group()
        # Crystal Properties
        crystal = refl_data.crystals()[0]
        new_cls.unit_cell = crystal.unit_cell()
        new_cls.symmetry = crystal.crystal_symmetry()
        # Column information
        new_cls.column_labels = refl_data.column_labels()
        return new_cls

    @classmethod
    def from_pdb(cls, pdb_file=None, pdb_input=None, id=None):
        assert [pdb_file, pdb_input].count(None)==1,'Provide pdb_file OR pdb_input'
        # Initialise class (populate manually)
        new_cls = cls(id=id)
        # Store filename
        new_cls.pdb_file = pdb_file
        new_cls._pdb_input = pdb_input
        # Create pdb input of structure
        pdb_input = new_cls.pdb_input()
        # Crystal Properties
        new_cls.symmetry = pdb_input.crystal_symmetry()
        new_cls.space_group = new_cls.symmetry.space_group()
        new_cls.unit_cell = new_cls.symmetry.unit_cell()
        return new_cls

class UnitCellVariation(Info):
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


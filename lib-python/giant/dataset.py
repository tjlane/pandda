import os, copy

import numpy

import iotbx.pdb, iotbx.mtz
from libtbx.utils import Sorry, Failure

from bamboo.common import Meta, Info
from bamboo.common.file import FileManager
from bamboo.common.path import easy_directory

from giant.io.pdb import strip_pdb_to_input
from giant.xray.data import extract_structure_factors
from giant.xray.crystal import CrystalSummary
from giant.xray.symmetry import get_crystal_contact_operators, apply_symmetry_operators, combine_hierarchies
from giant.structure.align import align_structures_rigid, align_structures_flexible


class _DatasetObj(object):


    def label(self, num=-1, tag=None):
        self.num = num
        if tag: self.tag = str(tag)
        return self


class ModelAndData(_DatasetObj):


    def __init__(self, model, data):
        """Convenience object to hold model-data pairs"""
        self.model = model
        self.data  = data
        self.num = self.tag = None
        self.file_manager = None
        self.meta = Meta()
        self.parent = None
        self.child = None
        self.children = []

    @classmethod
    def from_file(cls, model_filename=None, data_filename=None):
        model = data = None
        if model_filename:
            assert os.path.exists(model_filename), 'Model file does not exist!'
            model = CrystallographicModel.from_file(filename=model_filename)
        if data_filename:
            assert os.path.exists(data_filename),  'Experimental Data file does not exist!'
            data = ExperimentalData.from_file(filename=data_filename)
        return cls(model=model, data=data)

    def initialise_output_directory(self, dir):
        """Initialise a dataset output directory"""
        # Create a file and directory organiser
        self.file_manager = FileManager(rootdir=easy_directory(dir))

    def get_pickle_copy(self):
        """Get copy of self that can be pickled - some cctbx objects cannot be pickled..."""

        return self


class AtomicModel(_DatasetObj):


    def __init__(self, input, hierarchy):
        self.input = input
        self.hierarchy = hierarchy
        self.filename = None
        self.alignment = None

    @classmethod
    def from_file(cls, filename):
#        ih = iotbx.pdb.hierarchy.input(filename)
        ih = strip_pdb_to_input(filename, remove_ter=True)
        c = cls(input=ih.input, hierarchy=ih.hierarchy)
        c.filename = filename
        return c

    @classmethod
    def from_other(cls, other):
        return cls(input=other.input, hierarchy=other.hierarchy)

    def align_to(self, other_hierarchy, method='local', **kwargs):
        """Align this model to another"""
        assert not self.alignment, 'Already aligned!'
        assert isinstance(other_hierarchy, iotbx.pdb.hierarchy.root)
        assert method in ['global','local'], 'alignment method not supported'
        if method == 'global':
            self.alignment = align_structures_rigid(mov_hierarchy=self.hierarchy, ref_hierarchy=other_hierarchy, **kwargs)
        else:
            self.alignment = align_structures_flexible(mov_hierarchy=self.hierarchy, ref_hierarchy=other_hierarchy, **kwargs)
        return self.alignment


class CrystallographicModel(AtomicModel):


    def __init__(self, input, hierarchy):
        super(CrystallographicModel, self).__init__(input=input, hierarchy=hierarchy)
        self.crystal_symmetry = self.input.crystal_symmetry()
        if self.crystal_symmetry is None:
            raise Sorry('There is no crystal symmetry for this structure')
        self.unit_cell = self.crystal_symmetry.unit_cell()
        if self.unit_cell is None:
            raise Sorry('There is no unit cell information for this structure')
        self.space_group = self.crystal_symmetry.space_group()
        if self.space_group is None:
            raise Sorry('There is no space group information for this structure')

        self._crystal_contacts_operators = None

    def crystal_contact_operators(self, distance_cutoff=10):
        """Return the symmetry operations to obtain the crystallographic copies within buffer distance of the atomic model"""
        return get_crystal_contact_operators(
                            hierarchy=self.hierarchy,
                            crystal_symmetry=self.crystal_symmetry,
                            distance_cutoff=distance_cutoff)

    def crystal_contacts(self, distance_cutoff=10, combine_copies=False):
        """Return the crystallographic symmetry copies within buffer distance of the atomic model"""
        ops = self.crystal_contact_operators(distance_cutoff=distance_cutoff)
        sym_hierarchies, chain_mappings = apply_symmetry_operators(
                                                    hierarchy=self.hierarchy,
                                                    crystal_symmetry=self.crystal_symmetry,
                                                    sym_ops_mat=ops)

        if combine_copies: sym_hierarchies = combine_hierarchies(sym_hierarchies)
        return sym_hierarchies


class ExperimentalData(_DatasetObj):


    @classmethod
    def from_file(cls, filename):
        if filename.endswith('.mtz'):
            return XrayData.from_file(filename=filename)


class XrayData(ExperimentalData):


    def __init__(self, mtz_object):
        self._mtz_object = mtz_object
        self.filename = None
        self.summary = CrystalSummary.from_mtz(mtz_object=mtz_object)
        self.miller_arrays = {}
        self.fft_maps = {}

    def _remove_mtz_objects(self):
        self._mtz_object = None
        self.summary._mtz_object = None

    @classmethod
    def from_file(cls, filename):
        assert filename.endswith('.mtz'), 'Given filename is not an mtz file'
        c = cls(mtz_object=iotbx.mtz.object(filename))
        c.filename = filename
        # Have the filename so no need to hold the object in memory -- allows pickling of object
        c._remove_mtz_objects()
        return c

    def mtz_object(self):
        if self._mtz_object:
            return self._mtz_object
        elif self.filename:
            return iotbx.mtz.object(self.filename)
        else:
            raise Exception('No filename to load data from')

    def get_structure_factors(self, columns):
        """Extract a let of structure factors from the mtz_object"""
        assert columns.count(',') == 1
        return extract_structure_factors(self.mtz_object(), ampl_label=columns.split(',')[0], phas_label=columns.split(',')[1])


class ModelAndMap(_DatasetObj):


    def __init__(self, model=None, map=None):
        """Object defining a subsection of a Dataset (such as a chain, or subunit)"""
        self.model = model
        self.map = map



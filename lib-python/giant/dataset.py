import os, copy

import iotbx.pdb, iotbx.mtz

from giant.exceptions import Sorry


class _DatasetObj(object):

    def label(self, num=None, tag=None):
        if (num is not None):
            self.num = num
        if (tag is not None):
            self.tag = str(tag)
        return self


class AtomicModel(_DatasetObj):

    def __init__(self, input, hierarchy):
        self.input = input
        self.hierarchy = hierarchy
        self.filename = None

    @classmethod
    def from_file(cls, filename):
        from giant.io.pdb import strip_pdb_to_input
        ih = strip_pdb_to_input(
            filename,
            remove_ter = True,
        )
        c = cls(
            input = ih.input,
            hierarchy = ih.hierarchy,
        )
        c.filename = filename
        return c

    def align_to(self,
        other_hierarchy,
        method = 'local',
        **kwargs
        ):
        """Align this model to another"""

        from giant.structure.align import align_structures_rigid, align_structures_flexible

        assert isinstance(other_hierarchy, iotbx.pdb.hierarchy.root)
        if method not in ['global','local']:
            raise NotImplementedError('alignment method not supported: {}'.format(method))

        if (method == 'global'):
            alignment = align_structures_rigid(
                mov_hierarchy = self.hierarchy,
                ref_hierarchy = other_hierarchy,
                **kwargs
            )
        elif (method == 'local'):
            alignment = align_structures_flexible(
                mov_hierarchy = self.hierarchy,
                ref_hierarchy = other_hierarchy,
                **kwargs
            )

        return alignment


class CrystallographicModel(AtomicModel):

    def __init__(self, input, hierarchy):

        super(CrystallographicModel, self).__init__(
            input = input,
            hierarchy = hierarchy,
        )

        from giant.xray.crystal import CrystalInfo
        self.crystal = CrystalInfo.from_pdb(
            pdb_input = input,
        )

    def crystal_contact_operators(self, distance_cutoff=10):
        """Return the symmetry operations to obtain the crystallographic copies within buffer distance of the atomic model"""

        if self.crystal.crystal_symmetry is None:
            raise Sorry('Cannot generate crystal contact operators: there is no crystal symmetry for this structure')

        from giant.xray.symmetry import get_crystal_contact_operators

        return get_crystal_contact_operators(
            hierarchy = self.hierarchy,
            crystal_symmetry = self.crystal.crystal_symmetry,
            distance_cutoff = distance_cutoff,
        )

    def crystal_contacts(self, distance_cutoff=10, combine_copies=False):
        """Return the crystallographic symmetry copies within buffer distance of the atomic model"""

        ops = self.crystal_contact_operators(
            distance_cutoff = distance_cutoff,
        )

        from giant.xray.symmetry import apply_symmetry_operators, combine_hierarchies

        sym_hierarchies, chain_mappings = apply_symmetry_operators(
            hierarchy = self.hierarchy,
            crystal_symmetry = self.crystal.crystal_symmetry,
            sym_ops_mat = ops,
        )

        if (combine_copies is not False):
            sym_hierarchies = combine_hierarchies(sym_hierarchies)

        return sym_hierarchies


class ExperimentalData(_DatasetObj):

    def __init__(self):

        self.filename = None

    @classmethod
    def from_file(cls, filename):

        if filename.endswith('.mtz'):
            return XrayData.from_file(filename=filename)
        else:
            raise NotImplementedError('No data object for this filetype')


class XrayData(ExperimentalData):

    def __init__(self, mtz_object):

        super(XrayData, self).__init__()

        self.mtz_object = mtz_object

        from giant.xray.crystal import CrystalInfo
        self.crystal = CrystalInfo.from_mtz(
            mtz_object = mtz_object,
        )

    @classmethod
    def from_file(cls, filename):
        assert filename.endswith('.mtz'), 'Given filename is not an mtz file'
        c = cls(
            mtz_object = iotbx.mtz.object(filename),
        )
        c.filename = filename
        return c

    def get_structure_factors(self, columns):
        """Extract a let of structure factors from the mtz_object"""
        from giant.xray.data import extract_structure_factors
        assert columns.count(',') == 1
        amplitude, phase = columns.split(',')
        return extract_structure_factors(
            self.mtz_object,
            ampl_label = amplitude,
            phas_label = phase,
        )


class ModelAndData(_DatasetObj):

    ModelClass = CrystallographicModel
    DataClass = ExperimentalData

    def __init__(self, model, data):
        """
        Convenience object to hold model-data pairs
        Model: an atomic / crystallographic model object
            containing model with coordinates and other related data
        Data: a map / xray / em data object
            containing reflections / maps and other related data
        """
        self.model = model
        self.data  = data
        self.num = None
        self.tag = None
        self.parent = None
        self.children = []

    @classmethod
    def from_file(cls,
        model_filename = None,
        data_filename = None,
        ):

        model = None
        data = None

        if (model_filename is not None):

            if not os.path.exists(model_filename):
                raise IOError('Model file does not exist: {}'.format(model_filename))

            model = cls.ModelClass.from_file(
                filename = model_filename,
            )

        if (data_filename is not None):

            if not os.path.exists(data_filename):
                raise IOError('Experimental Data file does not exist: {}'.format(data_filename))

            data = cls.DataClass.from_file(
                filename = data_filename,
            )

        return cls(
            model = model,
            data = data,
        )


#class EmData:


#class MapData:

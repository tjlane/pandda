import os, copy

import iotbx.pdb, iotbx.mtz

from giant.exceptions import Sorry


class _DatasetObj(object):
    
    name = 'DummyDatasetClass'

    def __init__(self):

        self.tag = None
        self.num = None
        self.filename = None

    def __repr__(self):

        return '<{name} {tag} from file {filename}>'.format(
            name = self.name,
            tag = self.tag,
            filename = self.filename,
            )

    def __str__(self):

        s_ = (
            'Object Type: {name}\n'
            '| Label: {tag}\n'
            '| From file: {filename}\n'
            '`---->'
            ).format(
            name = self.name,
            tag = self.tag,
            filename = self.filename,
            )

        return s_.strip()

    def label(self, num=None, tag=None):
        if (num is not None):
            self.num = num
        if (tag is not None):
            self.tag = str(tag)
        return self


class AtomicModel(_DatasetObj):

    name = "AtomicModel"

    def __init__(self, input, hierarchy, filename=None):    

        super(AtomicModel, self).__init__()

        self.input = input
        self.hierarchy = hierarchy
        self.filename = (
            str(filename)
            if filename is not None
            else None
            )

    @classmethod
    def from_file(cls, filename):
        from giant.io.pdb import strip_pdb_to_input
        ih = strip_pdb_to_input(
            str(filename),
            remove_ter = True,
        )
        c = cls(
            input = ih.input,
            hierarchy = ih.hierarchy,
            filename = filename,
        )
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

    name = "CrystallographicModel"

    def __init__(self, input, hierarchy, filename=None):

        super(CrystallographicModel, self).__init__(
            input = input,
            hierarchy = hierarchy,
            filename = filename,
        )

        from giant.xray.crystal import CrystalInfo
        self.crystal = CrystalInfo.from_pdb(
            pdb_input = input,
        )

    def __str__(self):

        s_ = super(CrystallographicModel, self).__str__()

        s_ += (
            '\n'
            '| Crystal information: \n'
            '|\t{crystal}\n'
            '`---->'
            ).format(
            crystal = str(self.crystal).strip().replace('\n','\n|\t'),
            )

        return s_.strip()

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
    
    name = "ExperimentalData"

    def __init__(self):

        super(ExperimentalData, self).__init__()

    @classmethod
    def from_file(cls, filename):

        if filename.endswith('.mtz'):
            return CrystallographicData.from_file(
                mtz_filename = filename,
                )
        else:
            raise NotImplementedError('No data object for this filetype')


class CrystallographicData(ExperimentalData):

    name = "CrystallographicData"

    def __init__(self,
        mtz_object = None,
        mtz_filename = None, # standardise! 
        ):

        super(CrystallographicData, self).__init__()

        self._mtz_object = mtz_object
        self.filename = (
            str(mtz_filename)
            if mtz_filename is not None
            else None
            )

        from giant.xray.crystal import CrystalInfo
        self.crystal = CrystalInfo.from_mtz(
            mtz_object = self.mtz_object(),
        )

    def __str__(self):

        s_ = super(CrystallographicData, self).__str__()

        s_ += (
            '\n'
            '| Crystal information: \n'
            '|\t{crystal}\n'
            '`---->'
            ).format(
            crystal = str(self.crystal).strip().replace('\n','\n|\t'),
            )

        return s_.strip()

    @classmethod
    def from_file(cls, filename):
        assert str(filename).endswith('.mtz'), 'Given filename is not an mtz file'
        c = cls(
            mtz_filename = filename,
            #mtz_object = iotbx.mtz.object(filename),
        )

        return c

    def get_structure_factors(self, columns):
        """Extract a let of structure factors from the mtz_object"""
        if isinstance(columns, str):
            assert columns.count(',') == 1
            columns = columns.split(',')
        assert len(columns) == 2
        from giant.xray.data import extract_structure_factors
        return extract_structure_factors(
            self.mtz_object(),
            ampl_label = columns[0],
            phas_label = columns[1],
        )

    def mtz_object(self):

        if self._mtz_object is not None:
            return self._mtz_object

        return iotbx.mtz.object(
            self.filename
            )


class ModelAndData(_DatasetObj):

    name = "ModelAndData"

    ModelClass = None
    DataClass = None

    def __init__(self, model, data):
        """
        Convenience object to hold model-data pairs
        Model: an atomic / crystallographic model object
            containing model with coordinates and other related data
        Data: a map / xray / em data object
            containing reflections / maps and other related data
        """

        super(ModelAndData, self).__init__()

        self.model = model
        self.data  = data

    def __str__(self):

        s_ = (
            'Object Type: {name}\n'
            '| Label: {tag}\n'
            '| From file: {filename}\n'
            '| Model:\n'
            '|\t{model}\n'
            '| Data:\n'
            '|\t{data}\n'
            '`---->'
            ).format(
            name = self.name,
            tag = self.tag,
            filename = self.filename,
            model = str(
                self.model
                ).strip('\n').replace('\n','\n|\t'),
            data = str(
                self.data
                ).strip('\n').replace('\n','\n|\t'),
            )

        return s_.strip()

    @classmethod
    def from_file(cls,
        model_filename = None,
        data_filename = None,
        ):

        model = None
        data = None

        if (model_filename is not None):

            if not os.path.exists(str(model_filename)):
                raise IOError('Model file does not exist: {!s}'.format(model_filename))

            model = cls.ModelClass.from_file(
                filename = model_filename,
            )

        if (data_filename is not None):

            if not os.path.exists(str(data_filename)):
                raise IOError('Experimental Data file does not exist: {}'.format(data_filename))

            data = cls.DataClass.from_file(
                filename = data_filename,
            )

        return cls(
            model = model,
            data = data,
        )

class CrystallographicDataset(ModelAndData):

    name = "CrystallographicDataset"

    ModelClass = CrystallographicModel
    DataClass = CrystallographicData


#class EmData:


#class MapData:

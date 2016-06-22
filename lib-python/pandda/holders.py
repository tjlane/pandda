import cctbx.sgtbx
import cctbx.uctbx

from scitbx.array_family import flex

from bamboo.common.masks import mask_collection
from bamboo.common import Meta

from pandda.handlers import DatasetHandler

########################################################################################################
#
#   HOLDER CLASSES
#
########################################################################################################

class MapHolder(object):
    """Class to hold map values and meta data"""
    def __init__(self, num, tag, map, unit_cell, space_group, meta, parent=None, child=None):
        assert isinstance(num, int), 'Num must be int. Type given: {!s}'.format(type(num))
        self.num = num
        assert isinstance(tag, str), 'Tag must be str. Type given: {!s}'.format(type(tag))
        self.tag = tag
        assert isinstance(map, flex.double), 'Map data must be flex.double. Type given: {!s}'.format(type(map))
        self.map = map
        assert isinstance(unit_cell, cctbx.uctbx.ext.unit_cell), 'Unit cell must be of type unit_cell. Type given: {!s}'.format(type(unit_cell))
        self.unit_cell = unit_cell
        assert isinstance(space_group, cctbx.sgtbx.ext.space_group), 'Space group must be of type space_group. Type given: {!s}'.format(type(space_group))
        self.space_group = space_group
        assert isinstance(meta, Meta), 'Meta must be dict. Type given: {!s}'.format(type(meta))
        self.meta = meta
        if parent:
            assert isinstance(parent, DatasetHandler) or isinstance(parent, MapHolder), 'parent must be of type DatasetHandler or MapHolder. Type given: {!s}'.format(type(parent))
            self.parent = parent
        if child:
            self.child = child

########################################################################################################
#
#   HOLDER LIST CLASSES
#
########################################################################################################

class HolderList(object):
    """Class for grouping many holders together"""
    _holder_class = None

    def __init__(self, *args):
        self._holder_list = []
        self._masks = mask_collection()
        # Allow for subclassing
        self.__class_init__(args)

    def __class_init__(self, args):
        pass

    def __getitem__(self, idx):
        if   isinstance(idx, int): return self.get(num=idx)
        elif isinstance(idx, str): return self.get(tag=idx)
        else: raise Exception('CANNOT INDEX EXCEPT BY int OR str. TYPE GIVEN: {!s}'.format(type(idx)))

    def __call__(self):
        """Return all holders"""
        return self.all()

    def all(self):
        """Return all holders"""
        return self._holder_list
    def all_nums(self):
        """Return the list of holder ids"""
        return [h.num for h in self.all()]
    def all_tags(self):
        """Return the list of holder tags"""
        return [h.tag for h in self.all()]

    def all_masks(self):
        """Return the mask object"""
        return self._masks
    def mask(self, mask_name, invert=False):
        """Retrieve a masked list of datasets"""
        if not mask_name:
            return self.all()
        else:
            return self._masks.mask(mask_name=mask_name, input_list=self.all(), invert=invert)

    def size(self, mask_name=None, invert=False):
        """Return the number of holders in the list (with optional mask applied)"""
        if mask_name:   return len(self.mask(mask_name=mask_name, invert=invert))
        else:           return len(self.all())

    def add(self, new_holders):
        """Add new datasets"""

        for new_h in new_holders:
            # Check all added datasets are the right class
            if self._holder_class:
                assert isinstance(new_h, self._holder_class), 'OBJECTS MUST BE OF TYPE: {!s}\n(ADDED OF TYPE: {!s})'.format(self._holder_class, type(new_h))
            # Check all added dataset id tags are strs
            assert isinstance(new_h.tag, str), 'TAG MUST BE str. Type given: {!s}'.format(type(new_h.tag))
            assert new_h.tag not in self.all_tags(), 'HOLDER WITH TAG ALREADY EXISTS: {!s}'.format(new_h.tag)
            # Check all added dataset id nums are ints
            assert isinstance(new_h.num, int), 'NUM MUST BE int. Type given: {!s}'.format(type(new_h.num))
            assert new_h.num not in self.all_nums(), 'HOLDER WITH NUM ALREADY EXISTS: {!s}'.format(new_h.num)
            # No problems, add to list
            self._holder_list.append(new_h)

    def get(self, tag=None, num=None):
        """Get a dataset by tag or num"""

        assert [num, tag].count(None) == 1, 'Must give EITHER num OR tag'
        if num: matching = [m for m in self.all() if m.num == num]
        else:   matching = [m for m in self.all() if m.tag == tag]
        if len(matching) == 0: raise Exception('NO MATCHING HOLDER FOUND - NUM: {!s}, TAG: {!s}'.format(num, tag))
        if len(matching) != 1: raise Exception('MORE THAN ONE MATCHING HOLDER FOUND - NUM: {!s}, TAG: {!s}'.format(num, tag))
        return matching[0]

########################################################################################################
#
#   HOLDER LIST SUBCLASSES
#
########################################################################################################

class DatasetHandlerList(HolderList):
    """Class for grouping many dataset handlers together"""
    _holder_class = DatasetHandler

    def __class_init__(self, *args):
        pass

class MapHolderList(HolderList):
    """Class for grouping many MapHolder objects together"""
    _holder_class = MapHolder

    def __class_init__(self, *args):
        pass


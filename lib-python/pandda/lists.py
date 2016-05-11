
from bamboo.common import Meta, Info

########################################################################################################
#
#   MAP LIST CLASSES
#
########################################################################################################

class MapList(Info):
    _map_names = []
    def __init__(self, map_names=None):
        if map_names is None: map_names=[]
        assert self._map_names+map_names, 'No Map Names defined'
        for m in self._map_names+map_names:
            self.__dict__[m] = None
        # Intialise meta dict
        self.meta = Meta()
        # Set initialised so attributes can't be added
        self._initialized = True

class PanddaStatMapList(MapList):
    _map_names = ['mean_map','medn_map','stds_map','sadj_map','skew_map','kurt_map','bimo_map']

class PanddaMultipleStatMapList(object):
    def __init__(self):
        """Store lots of statistical maps"""
        self.map_lists = {}
    def get_resolutions(self):
        return sorted(self.map_lists.keys())
    def add(self, stat_map_list, resolution, overwrite=False):
        assert isinstance(resolution, float), 'Resolution of map must be of type float. Type given: {!s}'.format(type(resolution))
        if not overwrite: assert resolution not in self.map_lists.keys(), 'MAPS OF THIS RESOLUTION ALREADY ADDED'
        assert isinstance(stat_map_list, PanddaStatMapList), 'stat_map_list must be of type PanddaMultipleStatMapList. Type given: {!s}'.format(type(stat_map_list))
        self.map_lists[resolution] = stat_map_list
    def get(self, resolution):
        return self.map_lists[resolution]



from scipy import stats as scipy_stats
import numpy

class data_collection(object):
    def __init__(self):
        """Stores observations across a set of datasets"""
        self.data = {}
        self.entry_ids = None
        self._data_length = None

    def set_data_length(self, data_length):
        """Set the mask length"""
        if self._data_length: assert self._data_length == data_length, "CAN'T CHANGE THE DATASET LENGTH"
        else:                 self._data_length = data_length

    def set_entry_ids(self, entry_ids):
        """List of ids that label the different observations"""
        assert not self.entry_ids, 'DATASET IDS HAVE ALREADY BEEN SET'
        if self._data_length: assert len(entry_ids) == self._data_length
        else:                 self._data_length = len(entry_ids)
        self.entry_ids = entry_ids

    def add_empty_data(self, data_name):
        """Add a new empty set of data"""
        if data_name in self.data.keys():
            print('PROGRAM WARNING: OVERWRITING DATA IN DATA SUMMARY: {!s}'.format(data_name))
        self.data[data_name] = [numpy.nan]*self._data_length

    def add_data(self, data_name, data_values):
        """Add a new set of data"""
        if data_name in self.data.keys():
            print('PROGRAM WARNING: OVERWRITING DATA IN DATA SUMMARY: {!s}'.format(data_name))
        assert len(data_values) == self._data_length
        self.data[data_name] = data_values

    def get_data(self, data_name):
        """Get a single set of data"""
        return self.data[data_name]

    def has_data(self, data_name):
        """Check if it has a single mask"""
        return data_name in self.data.keys()

    def set_data_value(self, data_name, entry_id, value):
        """Set a particular entry in the data list `data_name`, corresponding to `entry_id` to `value`"""
        assert data_name in self.data.keys()
        assert entry_id in self.entry_ids
        self.data[data_name][self.entry_ids.index(entry_id)] = value

    def get_data_value(self, data_name, entry_id):
        """Retrieve an observation for a particular entry_id"""
        assert data_name in self.data.keys()
        assert entry_id in self.entry_ids
        return self.data[data_name][self.entry_ids.index(entry_id)]

    def get_data_as_z_values(self, data_name):
        """Convert a set of data to z-values"""
        data_values = self.get_data(data_name=data_name)
        z_values = scipy_stats.zscore(a=data_values).tolist()
        return zip(self.entry_ids, z_values)

    def get_data_as_dict(self, data_name):
        """Get a single set of data - returned as a dictionary with the entry_ids as keys"""
        return dict(zip(self.entry_ids, self.get_data(data_name=data_name)))

class multiple_data_collection(object):
    def __init__(self):
        """Stores multiple related data_collection objects - Allows the same observation to be pulled from all of the data_collections"""

        # Contains keys for the dictionaries
        self.collection_ids = []
        # Contains the collection objects
        self.collections = {}
        # Contains information about each collection object
        self.collection_info = {}
        # Number of collection objects
        self._collection_length = None

    def set_collection_ids(self, collection_ids):
        """Set the ids of the different collections"""
        if self._collection_length: assert len(collection_ids) == self._collection_length
        else:                       self._collection_length = len(collection_ids)
        assert not self.collection_ids, 'COLLECTION IDS HAVE ALREADY BEEN SET'
        self.collection_ids = collection_ids

    def add_collection(self, collection_id, collection):
        """Add a data collection"""
        assert collection_id in self.collection_ids, 'INVALID COLLECTION ID'
        self.collections[collection_id] = collection

    def get_collection(self, collection_id):
        """Get a data collection"""
        return self.collections[collection_id]

    def add_collection_info(self, info_name, info_values):
        """Store info about the different data collections"""
        assert self._collection_length == len(info_values)
        self.collection_info[info_name] = dict(zip(self.collection_ids, info_values))

    def get_collection_info(self, info_name):
        """Get info about the different data collections"""
        return self.collection_info[info_name]

    def get_data_from_collections(self, data_name):
        """Retrieve a zipped list of collection ids, and the data contained within them referenced by data_name"""
        return [(col_id, self.get_collection(col_id).get_data(data_name)) for col_id in self.collection_ids]

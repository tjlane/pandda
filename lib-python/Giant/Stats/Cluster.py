import scipy.spatial
import scipy.cluster
import numpy

class cluster_data(object):
    def __init__(self, input_dict):
        """
        Object for processing the clusters of points on a grid.
        Takes a dictionary as input of the form:
            {cluster_number/name: [(3d_point_tuple, value), ...], ...}

        """

        # Names/Numbers of Clusters (sorted for repeatability)
        self._cluster_keys = sorted(input_dict.keys())
        self._cluster_data = [input_dict[key] for key in self._cluster_keys]
        # Translate between cluster key and index in lists
        self._key_dict = dict([(key, i) for i, key in enumerate(self._cluster_keys)])

        # Number of clusters
        self._num_clusters = len(input_dict.keys())

        # Extract sets of points and values that make up the clusters
        self._cluster_sizes = [len(input_dict[key]) for key in self._cluster_keys]
        self._cluster_points = [[d[0] for d in input_dict[key]] for key in self._cluster_keys]
        self._cluster_values = [[d[1] for d in input_dict[key]] for key in self._cluster_keys]

        # Store raw data
        self.raw_data = input_dict

        # Cluster summary statistics
        self._cluster_means = [numpy.mean(g) for g in self._cluster_values]
        self._cluster_maxs  = [max(g) for g in self._cluster_values]

        # Representative points
        self._cluster_centroids = [tuple(numpy.mean(g, axis=0).tolist()) for g in self._cluster_points]
        self._cluster_peaks     = [g[v.index(v_max)] for v_max, v, g in zip(self._cluster_maxs, self._cluster_values, self._cluster_points)]

    def get_indices_from_keys(self, keys):
        """Returns the list of indices for a list of keys"""
        return [self._key_dict[key] for key in keys]

    def get_keys(self, indices=None):
        """Gets the list of sorted cluster keys"""
        if not indices: return self._cluster_keys
        else:           return [self._cluster_keys[i] for i in indices]

    def get_data(self, indices=None):
        """Gets the list of sorted (raw) cluster data"""
        if not indices: return self._cluster_data
        else:           return [self._cluster_data[i] for i in indices]

    def get_sizes(self, indices=None):
        """Get the sizes of the clusters"""
        if not indices: return self._cluster_sizes
        else:           return [self._cluster_sizes[i] for i in indices]

    def get_points(self, indices=None):
        """Get groups of points in each cluster"""
        if not indices: return self._cluster_points
        else:           return [self._cluster_points[i] for i in indices]

    def get_values(self, indices=None):
        """Get groups of values in each cluster"""
        if not indices: return self._cluster_values
        else:           return [self._cluster_values[i] for i in indices]

    def get_means(self, indices=None):
        """Get mean values of each cluster"""
        if not indices: return self._cluster_means
        else:           return [self._cluster_means[i] for i in indices]

    def get_maxima(self, indices=None):
        """Get max values of each cluster"""
        if not indices: return self._cluster_maxs
        else:           return [self._cluster_maxs[i] for i in indices]

    def get_centroids(self, indices=None):
        """Get centroids of points in each cluster"""
        if not indices: return self._cluster_centroids
        else:           return [self._cluster_centroids[i] for i in indices]

    def get_peaks(self, indices=None):
        """Get locations of cluster peaks"""
        if not indices: return self._cluster_peaks
        else:           return [self._cluster_peaks[i] for i in indices]

    def sort(self, sorting_data='values', sorting_function=None, decreasing=True):
        """Return sorted indices for clusters ordered by sorting_function applied to the cluster values or sizes - i.e. mean(values)"""

        if sorting_data=='values':
            sorting_vals = [sorting_function(v) for v in self._cluster_values]
        elif sorting_data=='sizes':
            sorting_vals = self._cluster_sizes
        else:
            raise Exception()

        # Indices of the sorted clusters
        sorted_indices = sorted(range(len(sorting_vals)), key=lambda k: sorting_vals[k], reverse=decreasing)

        return sorted_indices

    def create_new_cluster_from_mask(self, mask):
        """Returns selection of the cluster data based on a binary mask - can be Used to create new clusters from subsets of the data"""

        return cluster_data(dict([(self._cluster_keys[i], self._cluster_data[i]) for i, m in enumerate(mask) if m==1]))

def combine_clusters(clusters, ids=None):
    """Combine multiple clusters into one object. Keys will be given the prefix in `names` if given, or integers - clusters can be None (these will be ignored)"""

    if ids: assert len(ids)==len(clusters)
    else: ids = range(0,len(clusters))

    new_keys = []
    new_data = []

    for i_clust, old_cluster in enumerate(clusters):
        if old_cluster is None:
            continue
        new_keys.extend([(ids[i_clust], old_key) for old_key in old_cluster.get_keys()])
        new_data.extend(old_cluster.get_data())

    return cluster_data(dict(zip(new_keys, new_data)))

def find_connected_groups(connection_matrix):
    """Take a matrix of connected elements and output groups"""

    assert not set(connection_matrix.flatten()).difference(set([0,1])), 'CONNECTION MATRIX MUST CONSIST OF 1s AND 0s ONLY'
    assert connection_matrix.ndim == 2, 'CONNECTION MATRIX MUST BE OF DIMENSION 2'
    assert connection_matrix.shape[0] > 1, 'MATRIX MUST BE LARGER THAN 1x1'
    assert connection_matrix.shape[0] == connection_matrix.shape[1], 'CONNECTION MATRIX MUST BE SQUARE'

    # Make it symmetrical - if 1 is connected to 2, 2 is connected to 1
    sym_connection_matrix = ((connection_matrix + connection_matrix.T) != 0).astype(int)
    # Convert from a connection array to a distance array (0 for connected, 1 for unconnected)
    dist_mx = 1 - sym_connection_matrix
    # Convert to condensed distance matrix
    cd_dist_mx = scipy.spatial.distance.squareform(dist_mx)
    # Calculate the linkage matrix
    l = scipy.cluster.hierarchy.linkage(cd_dist_mx, method='single', metric='euclidean')
    # Cluster with very low cutoff
    connected_groups = scipy.cluster.hierarchy.fcluster(Z=l, t=0.01).tolist()

    return connected_groups


import numpy as np


class Cluster(object):

    def __init__(self, points, values):

        self.points = np.array(points)
        self.values = np.array(values)

        self.validate()

    def validate(self):

        assert self.values.ndim == 1
        assert self.points.ndim == 2
        assert self.points.shape == self.values.shape + (3,)

    def size(self):
        return len(self.points)

    def min(self):
        return np.min(self.values)

    def mean(self):
        return np.mean(self.values)

    def max(self):
        return np.max(self.values)

    def peak(self):
        idx = np.where(
            self.values == self.values.max()
            )[0][0]
        return self.points[idx]

    def centroid(self):
        return self.points.mean(axis=0)


class BasicClusterFinder(object):

    def __init__(self, 
        grid_clustering_cutoff, 
        negative_values, 
        cluster_method, 
        contour_level, 
        n_error = 10000,
        ):

        assert cluster_method in [
            "agglomerative_hierarchical",
            ], "still need to implement alternatives!"

        self.negative_values = negative_values
        self.contour_level = contour_level

        self.grid_clustering_cutoff = grid_clustering_cutoff
        self.cluster_method = cluster_method

        self.n_error = n_error

    def __call__(self, 
        map_data,
        ):

        map_data = self.truncate_map(
            map_data = map_data,
            )

        clusters = self.find_clusters(
            map_data = map_data,
            )

        return clusters

    def truncate_map(self,
        map_data,
        mask_value = 0.0,
        ):

        map_copy = map_data.copy()

        if (self.negative_values is True):
            map_copy = np.abs(map_copy)

        # Points to be masked
        zero_mask = (map_copy <= self.contour_level)

        # Mask map
        map_copy[zero_mask] = mask_value

        return map_copy

    def find_clusters(self,
        map_data,
        ):

        map_points, map_values = self.extract_points_and_values(
            map_data = map_data,
            )

        # Raise Error
        if (len(map_points) >= self.n_error):
            raise ValueError("Number of points exceeds allowed maximum")

        # Special case! 
        if len(map_points) == 0:
            return []

        # Special case! 
        if len(map_points) == 1: 
            cluster_labels = [0]
        else:
            cluster_labels = self.cluster_points(
                grid_points = map_points,
                )

        clusters = self.make_clusters(
            cluster_labels = cluster_labels,
            map_points = map_points,
            map_values = map_values,
            )

        return clusters

    def extract_points_and_values(self,
        map_data,
        ):

        assert (map_data.ndim == 3)

        all_points = np.transpose(
            np.indices(
                map_data.shape, 
                dtype = int,
                ),
            axes = (1,2,3,0),
            )

        sel_mask = (map_data != 0.0)

        return (
            all_points[sel_mask], 
            map_data[sel_mask],
            )

    def cluster_points(self,
        grid_points,
        ):

        if self.cluster_method == "agglomerative_hierarchical":

            return self.cluster_hierarcical(
                grid_points = grid_points,
                )

        raise NotImplemented()

    def cluster_hierarcical(self, 
        grid_points,
        ):

        import scipy.cluster

        cluster_ids = scipy.cluster.hierarchy.fclusterdata(
            X = grid_points,
            t = self.grid_clustering_cutoff,
            criterion = 'distance',
            metric = 'euclidean',
            method = 'single',
            )

        cluster_ids = np.array(cluster_ids)

        return cluster_ids

    def make_clusters(self, 
        cluster_labels,
        map_points,
        map_values,
        ):

        from giant.common.clustering import generate_group_idxs

        clusters = []

        for c_lab, c_idxs in generate_group_idxs(cluster_labels):

            # c_lab is not used!

            clusters.append(
                Cluster(
                    points = map_points[c_idxs],
                    values = map_values[c_idxs],
                    )
                )

        return clusters

# def cluster_hdbscan(above_gps):
#     sample_by_feature = above_gps.to_numpy()
#     print(sample_by_feature.shape)
#     clusterer = HDBSCAN()
#     clusterer.fit(sample_by_feature)
#     cluster_ids = list(clusterer.labels_)
#     return cluster_ids



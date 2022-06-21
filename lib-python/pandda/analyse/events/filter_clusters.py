import numpy as np 

from giant.common.clustering import (
    find_connected_groups, 
    generate_group_idxs,
    )

from .find_clusters import Cluster

# TODO: replace with giant.common.geometry.pairwise_distances
def pairwise_distances_between_sites(sites1, sites2):

    n1 = len(sites1)
    n2 = len(sites2)

    dm = (
        sites1.reshape((n1, 1, 3)).repeat(n2, axis=1) - sites2
        )

    return np.sqrt(
        (dm * dm).sum(axis=-1)
        )

def minimum_distance_between_sites(sites1, sites2):

    dm = pairwise_distances_between_sites(
        sites1 = sites1, 
        sites2 = sites2,
        )

    min_dist = dm.min()

    return min_dist

def minimum_distance_between_clusters(cluster1, cluster2):

    return minimum_distance_between_sites(
        sites1 = cluster1.points,
        sites2 = cluster2.points,
        )


class ClusterFilterList(object):

    def __init__(self, filters):
        
        self.filters = filters

    def __call__(self, 
        clusters,
        **kw_args
        ): 

        for f in self.filters:

            if len(clusters) == 0:
                break

            clusters = f(
                clusters = clusters,
                **kw_args
                )

        return clusters


class PeakAndSizeClusterFilter(object):

    def __init__(self,
        min_blob_z_peak,
        min_blob_volume,
        ):

        self.min_blob_z_peak = min_blob_z_peak
        self.min_blob_volume = min_blob_volume

    def __call__(self,
        clusters,
        map_grid,
        **kw_args
        ):

        clusters = self.filter_size_and_peak(
            clusters = clusters,
            min_peak = (
                self.min_blob_z_peak
                ),
            min_size = (
                self.min_blob_volume / (map_grid.grid_spacing()**3)
                ),
            )

        return clusters

    def filter_size_and_peak(self, 
        clusters,
        min_peak,
        min_size,
        ):

        out_clusters = []

        for c in clusters: 

            if c.size() < min_size:
                continue

            if c.max() < min_peak: 
                continue

            out_clusters.append(c)

        return out_clusters


class GroupNearbyClustersFilter(object):

    def __init__(self, 
        distance_cutoff,
        ):

        self.distance_cutoff = distance_cutoff

    def __call__(self,
        clusters,
        map_grid,
        **kw_args
        ):

        if len(clusters) == 1: 
            return clusters

        clusters = self.group_clusters(
            clusters = clusters,
            grid_distance_cutoff = (
                self.distance_cutoff / map_grid.grid_spacing()
                ),
            )

        return clusters

    def group_clusters(self,
        clusters,
        grid_distance_cutoff,
        ):

        # Record which clusters are to be joined
        connect_array = np.eye(
            len(clusters), 
            dtype = int,
            )

        for i1, c1 in enumerate(clusters):

            for i2, c2 in enumerate(clusters):

                if i1 == i2: 
                    break

                min_dist = minimum_distance_between_clusters(
                    cluster1 = c1,
                    cluster2 = c2,
                    )

                if (min_dist < grid_distance_cutoff): 
                    connect_array[i1,i2] = 1
                    connect_array[i2,i1] = 1

        # Cluster the connection array
        cluster_groupings = find_connected_groups(
            connection_matrix = connect_array,
            )

        # Concatenate smaller clusters into larger clusters
        out_clusters = []

        for g_lab, c_idxs in generate_group_idxs(cluster_groupings):

            new_cluster = Cluster(
                points = np.concatenate([clusters[i].points for i in c_idxs]),
                values = np.concatenate([clusters[i].values for i in c_idxs]),
                )

            assert sum([clusters[i].size() for i in c_idxs]) == new_cluster.size()

            out_clusters.append(new_cluster)

        return out_clusters


class ContactsClusterFilter(object):

    def __init__(self,
        distance_cutoff,
        selection_string = "not hetero", 
        ):

        self.distance_cutoff = distance_cutoff
        self.selection_string = selection_string

    def __call__(self,
        clusters,
        map_grid,
        dataset,
        ):

        clusters = self.find_clusters_with_contacts(
            clusters = clusters, 
            map_grid = map_grid,
            dataset = dataset,
            grid_distance_cutoff = (
                self.distance_cutoff / map_grid.grid_spacing()
                ),
            )

        return clusters

    def find_clusters_with_contacts(self, 
        clusters,
        map_grid,
        dataset,
        grid_distance_cutoff,
        ):

        # Extract the sites to the reference frame

        ref_sites_cart = dataset.model.alignment.nat2ref(
            self.get_filter_sites(dataset)
            )

        ref_sites_grid = map_grid.cart2grid(
            ref_sites_cart
            )

        # Remove clusters not within distance of sites

        out_clusters = []

        for c in clusters:

            min_dist = minimum_distance_between_sites(
                sites1 = ref_sites_grid,
                sites2 = c.points,
                )

            if (min_dist < grid_distance_cutoff): 

                out_clusters.append(c)

        return out_clusters

    def get_filter_sites(self, dataset):

        cache = dataset.model.hierarchy.atom_selection_cache()

        selection = cache.selection(
            self.selection_string,
            )

        atoms = dataset.model.hierarchy.atoms()

        xyz = atoms.select(selection).extract_xyz()

        return xyz


class SymmetryClusterFilter(object):

    def __init__(self,
        cluster_distance_cutoff,
        contact_distance_cutoff,
        selection_string = "not hetero",
        ):

        self.cluster_distance_cutoff = cluster_distance_cutoff
        self.contact_distance_cutoff = contact_distance_cutoff
        self.selection_string = selection_string

    def __call__(self,
        clusters,
        map_grid,
        dataset,
        ):

        if len(clusters) == 1: 
            return clusters

        clusters = self.remove_symmetry_related_clusters(
            clusters = clusters,
            map_grid = map_grid,
            dataset = dataset,
            )

        return clusters

    def remove_symmetry_related_clusters(self,
        clusters,
        map_grid,
        dataset,
        ):

        from giant.xray.symmetry import (
            find_symmetry_equivalent_groups,
            )

        d_unit_cell = dataset.model.crystal.unit_cell

        # Cartesianise and fractionalise the points in each of the clusters (in the crystallographic frame)

        cluster_points_cart = [
            dataset.model.alignment.ref2nat(
                map_grid.grid2cart(c.points)
                )
            for c in clusters
            ]

        cluster_points_frac = [
            d_unit_cell.fractionalize(points)
            for points in cluster_points_cart
            ]

        # Find the sets of clusters that are symmetry related

        sym_equiv_groups = find_symmetry_equivalent_groups(
            points_frac = cluster_points_frac,
            sym_ops = dataset.model.crystal_contact_operators(),
            unit_cell = d_unit_cell,
            cutoff_cart = self.cluster_distance_cutoff,
            )

        # Choose which cluster to use by the numebr of contacts with the protein

        d_sites_cart = self.get_filter_sites(dataset)

        # Iterate through and chose one from each group to keep

        out_clusters = []

        for c_lab, c_idxs in generate_group_idxs(sym_equiv_groups):

            if len(c_idxs) == 1: 
                out_clusters.append(
                    clusters[c_idxs[0]]
                    )

                continue

            # Count the number of contact for each cluster in the group
            c_contacts = [
                self.find_number_of_contacts(
                    sites1 = d_sites_cart,
                    sites2 = cluster_points_cart[i_c],
                    ) 
                for i_c in c_idxs
                ]

            # Find the cluster with the most contacts
            max_contacts = max(c_contacts)

            if (max_contacts == 0):

                # TODO just choose the first one for now -- revisit

                out_clusters.append(
                    clusters[c_idxs[0]]
                    )

            else:

                i_best = c_idxs[c_contacts.index(max_contacts)]

                out_clusters.append(
                    clusters[i_best]
                    )

        return out_clusters

    def get_filter_sites(self, dataset):

        cache = dataset.model.hierarchy.atom_selection_cache()

        selection = cache.selection(
            self.selection_string,
            )

        atoms = dataset.model.hierarchy.atoms()

        xyz = atoms.select(selection).extract_xyz()

        return np.array(xyz)

    def find_number_of_contacts(self, sites1, sites2):

        pairwise_distances = pairwise_distances_between_sites(
            sites1 = sites1,
            sites2 = sites2,
            )

        n_links = (pairwise_distances < self.contact_distance_cutoff).sum()

        return n_links

import numpy
import scipy.cluster

from scitbx.array_family import flex

from giant.stats.cluster import find_connected_groups, generate_group_idxs
from giant.xray.symmetry import find_symmetry_equivalent_groups


class PanddaZMapAnalyser(object):


    def __init__(self, params, grid_spacing, log):
        self.log = log
        self.params = params
        self.grid_spacing = grid_spacing
        self.grid_clustering_cutoff = 1.1 * numpy.math.sqrt(3)
        self.real_clustering_cufoff = self.grid_clustering_cutoff * grid_spacing
        self.grid_minimum_volume = int(self.params.min_blob_volume/(grid_spacing**3))

    def print_settings(self):
        self.log('----------------------------------->>>', True)
        self.log('Z-Scores Clustering', True)
        self.log('----------------------------------->>>', True)
        self.log('Clustering Points with Z-Scores > {!s}'.format(self.params.contour_level), True)
        self.log('----------------------------------->>>', True)
        self.log('Clustering Cutoff (A):       {!s}'.format(self.real_clustering_cufoff), True)
        self.log('----------------------------------->>>', True)
        self.log('Minimum Cluster Z-Peak:      {!s}'.format(self.params.min_blob_z_peak), True)
        self.log('Minimum Cluster Volume (A):  {!s}'.format(self.params.min_blob_volume), True)
        self.log('Minimum Cluster Size:        {!s}'.format(self.grid_minimum_volume), True)

    def cluster_high_z_values(self, z_map_data, point_mask):
        """Finds all the points in the z-map above `z_cutoff`, points will then be clustered into groups of cutoff `clustering_cutoff` angstroms"""

        point_mask_gps = flex.vec3_double(point_mask)
        # Convert the gps to idxs
        point_mask_idx = flex.size_t(map(z_map_data.accessor(), point_mask))
        # Select these values from the map
        point_mask_val = z_map_data.select(point_mask_idx)
        # Find values above cutoff
        above_bool = (point_mask_val >= self.params.contour_level)
        above_val = point_mask_val.select(above_bool)
        above_gps = point_mask_gps.select(above_bool)
        above_len = len(above_val)

        # No Cluster points found
        if   above_len == 0:
            return 0, []
        # One Cluster point found
        elif above_len == 1:
            return 1, [(above_gps, above_val)]
        # Can't cluster if there are too many points
        elif above_len > 10000:
            return -1, [(above_gps, above_val)]
        # Cluster points if we have found them
        else:
            self.log('> Clustering {!s} Points.'.format(above_len))
            # Cluster the extracted points
            t1 = time.time()
            cluster_ids = scipy.cluster.hierarchy.fclusterdata( X = above_gps,
                                                                t = self.grid_clustering_cutoff,
                                                                criterion = 'distance',
                                                                metric    = 'euclidean',
                                                                method    = 'single' )
            cluster_ids = list(cluster_ids)
            t2 = time.time()
            self.log('> Clustering > Time Taken: {!s} seconds'.format(int(t2-t1)))

            # Get the number of clusters
            num_clusters = max(cluster_ids)
            # Group the values by cluster id
            z_clusters = []
            for c_id, c_idxs in generate_group_idxs(cluster_ids):
                c_idxs = flex.size_t(c_idxs)
                c_gps = above_gps.select(c_idxs)
                c_val = above_val.select(c_idxs)
                z_clusters.append((c_gps, c_val))
            assert num_clusters == len(z_clusters)
            return num_clusters, z_clusters

    def validate_clusters(self, z_clusters):
        for i, (gps, vals) in enumerate(z_clusters):
            #print('Cluster {}'.format(i))
            #print('Points: {}'.format(len(gps)))
            #print('Values: {}'.format(len(vals)))
            assert len(gps) == len(vals)

    def filter_z_clusters_1(self, z_clusters):
        """Filter the z-clusters on a variety of criteria (size, peak value)"""

        self.log('----------------------------------->>>')
        self.log('Filtering by blob size and peak value')

        filt_z_clusters = z_clusters
        # Filter out small clusters - get numbers of clusters satisfying the minimum cluster size
        large_clusters = (flex.int([x[1].size() for x in filt_z_clusters]) >= self.grid_minimum_volume).iselection()
        if large_clusters.size() == 0:  return 0, []
        filt_z_clusters = [filt_z_clusters[i] for i in large_clusters]
        # Filter out weak clusters - get numbers of clusters satisfying the minimum z_peak value
        strong_clusters = (flex.double([x[1].min_max_mean().max for x in filt_z_clusters]) >= self.params.min_blob_z_peak).iselection()
        if strong_clusters.size() == 0: return 0, []
        filt_z_clusters = [filt_z_clusters[i] for i in strong_clusters]

        self.log('Filtered {!s} Clusters to {!s} Clusters'.format(len(z_clusters), len(filt_z_clusters)))
        return len(filt_z_clusters), filt_z_clusters

    def filter_z_clusters_2(self, z_clusters, grid_origin_cart, ref_structure, min_contact_dist=6):
        """Find and remove clusters more than a minimum distance from the protein"""

        # min_contact_dist - blobs are rejected if they are more than this distance from the protein

        self.log('----------------------------------->>>')
        self.log('Filtering by minimum distance from protein')

        # Extract the protein sites
        ref_sites_cart = ref_structure.select(ref_structure.atom_selection_cache().selection('pepnames')).atoms().extract_xyz()
        # Save time - calculate the square of the contact distance
        min_contact_dist_sq = min_contact_dist**2

        # Remove any clusters that are more than min_contact_dist from the protein
        filtered_c_idxs = []
        for c_idx, (c_gps, c_val) in enumerate(z_clusters):
            # Extract points in cluster
            cluster_points_cart = (c_gps * self.grid_spacing) - grid_origin_cart
            # Calculate minimum distance to protein
            for r_site_cart in ref_sites_cart:
                diff_vecs_cart = cluster_points_cart - r_site_cart
                # Keep cluster if minimum distance is less than min_contact_dist
                if min(diff_vecs_cart.dot()) < min_contact_dist_sq:
                    filtered_c_idxs.append(c_idx)
                    break
            # Report
#            if self.log.verbose:
#                if filtered_c_idxs and (filtered_c_idxs[-1] == c_idx):
#                    print('KEEPING CLUSTER:', c_idx)
#                else:
#                    print('REJECTING CLUSTER:', c_idx)
        # Select filtered clusters
        filt_z_clusters = [z_clusters[i] for i in filtered_c_idxs]

        self.log('Filtered {!s} Clusters to {!s} Clusters'.format(len(z_clusters), len(filt_z_clusters)))
        return len(filt_z_clusters), filt_z_clusters

    def filter_z_clusters_3(self, z_clusters, dataset, max_contact_dist=8):
        """Find and remove symmetry equivalent clusters"""

        if len(z_clusters) == 1:
            return 1, z_clusters
        else:
            self.log('----------------------------------->>>')
            self.log('Filtering symmetry equivalent clusters')

        # Extract the protein sites
        ref_sites_cart = ref_structure.select(ref_structure.atom_selection_cache().selection('pepnames')).atoms().extract_xyz()

        # Cartesianise and fractionalise the points in each of the clusters (in the crystallographic frame)
        points_cart = [None]*len(z_clusters)
        points_frac = [None]*len(z_clusters)
        for c_idx, (c_gps, c_val) in enumerate(z_clusters):
            # Extract points in cluster
            points_cart[c_idx] = (c_gps * self.grid_spacing) - grid_origin_cart
            # Fractionalise them to the unit cell of the reference structure
            points_frac[c_idx] = ref_unit_cell.fractionalize(points_cart[c_idx])
        # Find the sets of clusters that are symmetry related
        sym_equiv_groups = find_symmetry_equivalent_groups( points_frac = points_frac,
                                                            sym_ops     = ref_sym_ops,
                                                            unit_cell   = ref_unit_cell,
                                                            cutoff_cart = 1.05*1.7321*self.grid_spacing )
        # max_contact_dist - a point contacts an atom if the atoms is within this distance of it
        # Save time - calculate the square of the contact distance
        max_contact_dist_sq = max_contact_dist**2
        # Iterate through and chose one from each group to keep
        filt_z_clusters = []
        for g_id, g_idxs in generate_group_idxs(sym_equiv_groups):
            # Count the number of contact for each cluster in the group
            c_contacts = []
            # Iterate through cluster in the group
            for c_idx in g_idxs:
                # Initialise contact counter
                contacts = 0
                # Get the cartesian points for the cluster
                c_points_cart = points_cart[c_idx]
                # Again, use the brute force all-v-all method
                for rp in ref_sites_cart:
                    diffs_cart = c_points_cart - rp
                    # Check to see if site closer to cluster than minimum
                    if min(diffs_cart.dot()) < max_contact_dist_sq:
                        contacts += 1
                # Record the number of contacts (over size of cluster)
                c_contacts.append(1.0*contacts/len(c_points_cart))
#                if self.log.verbose:
#                    print('CLUSTER:', c_idx, ', CONTACTS PER POINT:', round(c_contacts[-1],3))

            # Find the cluster with the most contacts
            max_contacts = max(c_contacts)
            if max_contacts == 0:
                raise Exception('MAX CONTACTS IS 0!')
            else:
                cluster_to_keep = g_idxs[c_contacts.index(max_contacts)]
                filt_z_clusters.append(z_clusters[cluster_to_keep])
#                if self.log.verbose:
#                    print('KEEPING CLUSTER', cluster_to_keep)
        assert len(filt_z_clusters) == max(sym_equiv_groups), 'NUMBER OF UNIQUE GROUPS AND GROUPS TO BE RETURNED NOT THE SAME'

        self.log('Filtered {!s} Clusters to {!s} Clusters'.format(len(z_clusters), len(filt_z_clusters)))
        return len(filt_z_clusters), filt_z_clusters

    def group_clusters(self, z_clusters, separation_cutoff=5):
        """Join clusters that are separated by less than max_separation"""

        if len(z_clusters) == 1:
            return 1, z_clusters
        else:
            self.log('----------------------------------->>>')
            self.log('Grouping Nearby Clusters')

        # Minimum distance between grid points to be joined (squared)
        grid_cutoff_sq = (separation_cutoff/self.grid_spacing)**2

        # Record which clusters are to be joined
        connect_array = numpy.zeros((len(z_clusters),len(z_clusters)), dtype=int)
        for i_clust_1, (c_gps_1, c_val_1) in enumerate(z_clusters):
            for i_clust_2, (c_gps_2, c_val_2) in enumerate(z_clusters):
                # Skip if this is the same blob
                if i_clust_1 == i_clust_2:
                    connect_array[(i_clust_1, i_clust_2)] = 1
                    continue
                # Extract the minimum separation of the grid points
                min_dist_sq = min([min((c_gps_2 - gp).dot()) for gp in c_gps_1])
                # Check to see if they should be joined
                if min_dist_sq < grid_cutoff_sq:
                    connect_array[(i_clust_1, i_clust_2)] = 1
        # Cluster the connection array
        cluster_groupings = find_connected_groups(connection_matrix=connect_array)
        # Concatenate smaller clusters into larger clusters
        grouped_clusters = []
        for g_id, g_idxs in generate_group_idxs(cluster_groupings):
            g_gps = []; [g_gps.extend(z_clusters[i][0]) for i in g_idxs]
            g_gps = flex.vec3_double(g_gps)
            g_val = []; [g_val.extend(z_clusters[i][1]) for i in g_idxs]
            g_val = flex.double(g_val)
            grouped_clusters.append((g_gps, g_val))

        assert len(grouped_clusters) == max(cluster_groupings)

        self.log('Grouped {!s} Clusters together to form {!s} Clusters'.format(len(z_clusters), len(grouped_clusters)))
        return len(grouped_clusters), grouped_clusters


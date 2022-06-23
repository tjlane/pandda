import collections
import numpy
from scitbx.matrix import sym
from libtbx import adopt_init_args, group_args
from pandemic.adp.echt.analysis.cluster_tls_groups.utils import \
    make_selection_strings, tls_group_comparison_matrix, tls_group_adjacency_matrix
from pandemic.adp.echt.analysis.cluster_tls_groups.pymol import base_pymol_script

def uij_helliger_distance(a,b):
    """Calculate helliger distance between two variance tensors"""
    sa = sym(sym_mat3=a)
    sb = sym(sym_mat3=b)
    da = max(0.0, sa.determinant())
    db = max(0.0, sb.determinant())
    dab = max(0.0, (0.5*(sa+sb)).determinant())
    fa = (da*db)**0.25
    fb = (dab)**0.5
    if (fa == 0.0) and (fb == 0.0): return 0.0
    d = max(0.0, ( 1.0 -  ( fa / fb  )  ))
    return d ** 0.5

def uij_helliger_distance_function(uijs_a, uijs_b):
    """Calculate average helliger distance between two arrays of variance tensors"""
    u_dists = [uij_helliger_distance(a,b) for a,b in zip(uijs_a,uijs_b)]
    u_dist = numpy.mean(u_dists)
    return u_dist

def make_pymol_script(
    tls_groups,
    connectivity,
    distance_matrix,
    indices_sets,
    model_structure,
    output_filename,
    ):
    """Make pymol script to show results of helliger distance clustering"""

    from giant.pymol_utils import shapes

    minimum = 0.
    maximum = 1.
    midpoint = 0.5

    def map_value_to_colour(v, midpoint=0.5):
        midpoint = max(midpoint, 1e-6)
        assert 0. <  midpoint <= 1.0
        assert 0. <= v <= 1.0
        r = min(
            (0. + v)/(0. + midpoint),
            1.0,
            )
        g = min(
            (1. - v)/(1. - midpoint),
            1.0,
            )
        b = min(
            (1. - v)/(1. - midpoint),
            (0. + v)/(0. + midpoint),
            )
        return (r,g,b)

    def map_value_to_radius(v, max_radius=0.5):
        r = max_radius*(1.0-v)
        return r

    sc = base_pymol_script(
        tls_groups = tls_groups,
        connectivity_matrix = connectivity,
        model_structure = model_structure,
        output_filename = output_filename,
        )

    sc.add_comparison_matrix(
        comparison_matrix = distance_matrix,
        colour_function = map_value_to_colour,
        radius_function = map_value_to_radius,
        obj_name = 'tls-group-similarities',
        min_max_values = (1.0, 0.0),
        )

    s = sc.script

    ##################################################

    # Iterate through different clustering thresholds
    for i_iter, cluster_sets in enumerate(cluster_sets):
        # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        break # XXXXXXXXXXXXXXXXXXXXXXXXXX
        # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        group_shapes = []
        # Iterate through the clusters
        for c_indices in cluster_sets:
            # Select pairs of clusters
            for i in c_indices:
                for j in c_indices:
                    if i == j:
                        continue
                    if not connectivity[i,j]:
                        continue

                    cyl = shapes.Cylinder(
                        start=sc.centroids[i], end=sc.centroids[j], colors=[sc.bg_colour,sc.bg_colour],
                        radius=map_value_to_radius(0),
                    )
                    group_shapes.append(cyl)

        s.add_shapes(group_shapes, obj='clusters', state=i_iter+2) # states 2+

    sc.write()

    return


class ClusterTLSGroups_HelligerDistance(object):


    _comparison_functions = {
        'helliger_distance' : uij_helliger_distance_function,
        }

    _linkage_functions = {
        'single' : numpy.min,
        'average': numpy.mean,
        'complete': numpy.max,
    }

    def __init__(self,
        xyz_cutoff_distance = 4.0,
        comparison = 'helliger_distance',
        linkage = 'average',
        ):
        adopt_init_args(self, locals())

        if linkage not in self._linkage_functions:
            raise Sorry('Invalid linkage type ({}). Must be one of: {}'.format(linkage, ', '.join(self._linkage_functions.keys())))

        if self.comparison not in self._comparison_functions:
            raise Sorry('Invalid comparison function ({}). Must be one of: {}'.format(self.comparison, ', '.join(self._comparison_functions.keys())))

    def run(self,
        tls_groups,
        ):

        uij_distance_function = self._comparison_functions[self.comparison]

        tls_distance = tls_group_comparison_matrix(
            tls_groups = tls_groups,
            comparison_function = uij_distance_function,
            make_symmetric = 'min',
            )
        xyz_distance = tls_group_adjacency_matrix(
            tls_groups = tls_groups,
            )

        connectivity = (xyz_distance < self.xyz_cutoff_distance)

        clustering_info = self.do_hierarchical_clustering(
            connectivity = connectivity,
            distance_matrix = tls_distance,
            )

        processed_clusters = self.resolve_clustering(
            clustering_info = clustering_info,
            )

        self.result = group_args(
            distance_matrix = tls_distance,
            connectivity = connectivity,
            thresholded_module_lists = None,
            thresholded_hierarchy = None,
            output_hierarchy = None,
            )

        return self.result

    def do_hierarchical_clustering(self,
        distance_matrix,
        connectivity,
        ):

        n_points = len(distance_matrix)

        linkage_func = self._linkage_functions[self.linkage]

        # Dictionary of contents of each cluster (ALL)
        cluster_children = {i:[i] for i in range(n_points)}
        # Currently unmerged clusters
        current_clusters = list(range(n_points))
        # Calculated linkages between clusters
        cluster_linkages = {}
        # list of merges
        merging_steps = []
        # Neighbours of clusters
        cluster_neighbours = {}
        for i in range(n_points):
            j_indices = list(numpy.where(connectivity[i])[0])
            assert i in j_indices, 'Connectivity matrix must connect points to themselves (matrix must have unitary diagonal)'
            j_indices.remove(i)
            cluster_neighbours[i] = set(j_indices)

        # Calculate hierarchical clustering (max steps == n_points-1)
        for k_iter in range(n_points-1):

            # New cluster number
            k_clust = n_points + k_iter

            min_pair = None
            min_linkage = 1.0

            for j_clust in sorted(current_clusters):
                # Skip if not connected
                if j_clust not in cluster_neighbours:
                    continue
                for i_clust in sorted(current_clusters):
                    # Skip if not connected to anything
                    if i_clust not in cluster_neighbours:
                        continue
                    # Make sure clusters are compared only once
                    if i_clust >= j_clust:
                        break
                    # Clusters connected?
                    if not (cluster_neighbours[i_clust].intersection([j_clust]) or cluster_neighbours[j_clust].intersection([i_clust])):
                        continue

                    # Use precalculated linkage where possible
                    linkage = cluster_linkages.get((i_clust, j_clust), None)

                    if linkage is None:
                        i_children = cluster_children[i_clust]
                        j_children = cluster_children[j_clust]
                        i_j_distances = distance_matrix[i_children][:,j_children]
                        linkage = linkage_func(i_j_distances)
                        cluster_linkages[(i_clust, j_clust)] = linkage

                    #print linkage
                    if linkage < min_linkage:
                        min_linkage = linkage
                        min_pair = (i_clust, j_clust)

            # No pairs to merge
            if min_pair is None:
                break

            logger.debug('Merge {} and {} into {} with linkage {:.3f}'.format(min_pair[0], min_pair[1], k_clust, min_linkage))

            merging_steps.append((min_linkage, min_pair))

            i_merge, j_merge = min_pair
            i_children = cluster_children[i_merge]
            j_children = cluster_children[j_merge]
            k_children = sorted(i_children+j_children)

            # Update existing clusters
            current_clusters.remove(i_merge)
            current_clusters.remove(j_merge)

            # Create new cluster
            current_clusters.append(k_clust)
            cluster_children[k_clust] = k_children
            cluster_neighbours[k_clust] = set()

            # Done?
            if len(current_clusters) == 1:
                break

            # Calculate new neighbours
            for l_clust in current_clusters:
                if l_clust == k_clust:
                    continue
                l_children = cluster_children[l_clust]
                if connectivity[l_children][:,k_children].any():
                    # Append to l - SHOULD ALREADY EXIST
                    cluster_neighbours.get(l_clust).add(k_clust)
                    # Append to k
                    cluster_neighbours.get(k_clust).add(l_clust)

        return group_args(
            n_points = n_points,
            linkage = linkage,
            children = cluster_children,
            #neighbours = cluster_neighbours,
            #linkages = cluster_linkages,
            merging_steps = merging_steps,
            )

    def resolve_clustering(self,
        clustering_info,
        selection_strings,
        max_thresholds = None,
        ):

        # Output list of clusterings for different linkages
        output_linkages = []
        output_clusters = []
        output_selections = []

        # Extract from input
        cluster_children = clustering_info.children
        merging_steps = clustering_info.merging_steps

        n_merges = len(merging_steps)
        n_points = clustering_info.n_points

        if max_thresholds is None:
            max_thresholds = n_merges
        n_steps = min(max_thresholds, n_merges)

        n_print = max(1, min(10, n_merges))
        n_delta = int(float(n_steps) / float(n_print))

        for i_step in range(1, n_steps+1):

            stop_step = int(float(n_merges) * float(i_step) / float(n_steps))

            active_steps = merging_steps[:stop_step+1] # inclusive

            final_clusters = list(range(n_points))
            for i_merge, (merge_linkage, merge_pair) in enumerate(active_steps):
                i_new = n_points + i_merge
                final_clusters.remove(merge_pair[0])
                final_clusters.remove(merge_pair[1])
                final_clusters.append(i_new)

            max_linkage = max(list(zip(*active_steps))[0])
            clusters = [cluster_children[f] for f in final_clusters]
            strings = make_selection_strings(
                selections = selection_strings,
                clusters = clusters,
                )

            logger.debug('Clusters at step {} (linkage {:.3f})'.format(i_step, max_linkage))
            for s in strings:
                logger.debug(s)

            output_linkages.append(max_linkage)
            output_clusters.append(clusters)
            output_selections.append(strings)

        return group_args(
            cluster_sets = output_clusters,
            maximimum_linkages = output_linkages,
            selection_strings = output_selections,
            )

    def write_output(self,
        tls_groups,
        model_structure,
        output_prefix,
        output_directory = './',
        ):

        from pandemic.adp.echt.analysis.cluster_tls_groups.plots import make_linkage_graphs, Dendrogram

        cluster_info = self.result

        output_files = collections.OrderedDict()

        filename = str(output_directory / (output_prefix+'-clustering-process.png'))
        ret = make_linkage_graphs(
            cluster_info = cluster_info.clustering_info,
            output_filename = filename,
            )
        output_files['linkage_graph'] = filename

        filename = str(output_directory / (output_prefix+'-dendro.png'))
        dendro = Dendrogram(
            n_points = cluster_info.clustering_info.n_points,
            children_hash = cluster_info.clustering_info.children,
            merging_steps = cluster_info.clustering_info.merging_steps,
            )
        ret = dendro(output_filename = filename)
        output_files['dendrogram'] = filename

        filename = str(output_directory / (output_prefix+'-pymol.py'))
        ret = make_pymol_script(
            tls_groups = tls_groups,
            connectivity = cluster_info.connectivity,
            distance_matrix = cluster_info.distance_matrix,
            indices_sets = cluster_info.processed_clusters.cluster_sets,
            output_filename = filename,
            model_structure = model_structure,
            )
        output_files['pymol_script'] = filename

        return output_files


import giant.logs as lg
logger = lg.getLogger(__name__)

from libtbx import adopt_init_args, group_args

import copy, collections
import numpy

from sklearn.cluster import dbscan


class IdentifyModules(object):

    debug = False

    def __init__(self,
        threshold_start = 1.0,
        threshold_delta = 1.0,
        comparison_metric = 'similarity',
        ):

        """
        Identify groups of indices that cluster together by a certain "comparison" metric.

        comparions_metric: 'similarity' or 'distance'
        similarity -> modules are identified from low similarity (low thresholds) -> high similarity (high thresholds)
        distance -> modules are identified from low similarity (high thresholds) -> high similarity (low thresholds)
        The order of output modules is therefore the same for both in terms of similiarity,
        even though in one case the thresholds increase and in the other they decrease.
        """

        assert comparison_metric in ['similarity','distance'], 'Invalid comparison_metric ({}) provided'.format(comparison_metric)
        if comparison_metric == 'similarity':
            # Similarity : Increasing threshold
            assert threshold_start >= 0.0, 'start_threshold must be positive: similarity thresholds must be positive'
            assert threshold_delta > 0.0,  'delta_threshold must be positive: search is from small thresholds to large thresholds'
            threshold_function = self._threshold_above
        else:
            # Distance : Decreasing threshold
            assert threshold_start >= 0.0, 'start_threshold must be positive: distance thresholds must be positive'
            assert threshold_delta < 0.0,  'delta_threshold must be negative: search is from large thresholds to small thresholds'
            threshold_function = self._threshold_below

        adopt_init_args(self, locals())

    def __call__(self,
        connectivity,
        comparison_matrix,
        filter_complete_sets = False,
        ):

        n = len(comparison_matrix)

        threshold = self.threshold_start
        delta = self.threshold_delta

        logger.debug('Threshold jumps: {} (start), {} (delta)'.format(threshold, delta))

        # "raw" modules identified at each threshold
        threshold_modules = collections.OrderedDict()

        # Break condition for decreasing thresholds (won't affect increasing thresholds -- different break condition)
        while threshold > 0.0:

            modules = self.find_modules(
                connectivity = connectivity,
                comparison_matrix = comparison_matrix,
                threshold = threshold,
                )
            modules = self.filter_subsets(modules)

            # Remove trivial sets of all
            if (filter_complete_sets is True):
                modules = self.filter_complete_sets(modules, n_total=n)
                # Skip to next if still finding everything to avoid breaking below
                if not modules: continue

            # Stop once we stop finding modules
            if not modules: break

            # Store output modules
            threshold_modules[threshold] = modules

            # Update threshold
            threshold += delta

        ##############################################

        # Modules added at each level
        threshold_unique_modules = collections.OrderedDict()
        threshold_unique_hierarchies = collections.OrderedDict()
        # Cumulative hierarchy for each threshold
        threshold_cumulative_hierarchies = collections.OrderedDict()
        threshold_reduced_hierarchies = collections.OrderedDict()

        # Cumulative lists for all levels
        cumulative_modules = []
        cumulative_hierarchy = []

        # Need to iterate through in reverse order (remember -- using sorted dict!)
        threshold_keys = list(threshold_modules.keys())
        threshold_keys.reverse()
        # Create output levels for each threshold -- REMEMBER NOW USING REVERSED ORDER
        for threshold in threshold_keys:

            # Extract modules for this threshold
            modules = threshold_modules[threshold]

            # Filter modules that have been seen at a higher level
            unique_modules = self.filter_modules_by_reference(filter_modules=modules, reference_modules=cumulative_modules)
            # Skip if no new modules
            if not unique_modules: continue

            # Sort the new modules
            self.sort_modules_in_place(unique_modules)

            # Create hierarchy from the unique modules
            unique_hierarchy = self.create_hierarchy_from_modules(unique_modules, resort=False) # already sorted

            # Prepend modules to total list
            cumulative_modules = unique_modules + cumulative_modules
            # Prepend hierarchy to total list
            cumulative_hierarchy = unique_hierarchy + cumulative_hierarchy

            # Reduce hierarchy (sort modules)
            reduced_hierarchy = self.create_hierarchy_from_modules(cumulative_modules, resort=True) # already sorted

            # Store output
            threshold_unique_modules[threshold] = [tuple(sorted(m)) for m in unique_modules]
            threshold_unique_hierarchies[threshold] = self.sets_to_tuples(unique_hierarchy)
            threshold_cumulative_hierarchies[threshold] = self.sets_to_tuples(cumulative_hierarchy)
            threshold_reduced_hierarchies[threshold] = self.sets_to_tuples(reduced_hierarchy)

        # Output hierarchy is the last module hierarchy
        reduced_hierarchy = self.create_hierarchy_from_modules(cumulative_modules, resort=True)

        if (self.debug is True):
            self.show('threshold_unique_modules', threshold_unique_modules.items())
            self.show('threshold_cumulative_hierarchies', threshold_cumulative_hierarchies.items())
            for threshold, r_hierarchy in threshold_reduced_hierarchies.items():
                self.show('Reduced Hierarchy (@ Threshold {}):'.format(threshold), r_hierarchy)
            self.show('Total Hierarchy (all thresholds):', reduced_hierarchy)

        return group_args(
            threshold_unique_modules     = threshold_unique_modules,
            threshold_unique_hierarchies = threshold_unique_hierarchies,
            threshold_cumulative_hierarchies = threshold_cumulative_hierarchies,
            threshold_reduced_hierarchies = threshold_reduced_hierarchies,
            cumulative_hierarchy = self.sets_to_tuples(cumulative_hierarchy),
            reduced_hierarchy    = self.sets_to_tuples(reduced_hierarchy),
            )

    def _threshold_above(self, matrix, threshold_value):
        return (matrix >= threshold_value)

    def _threshold_below(self, matrix, threshold_value):
        return (matrix <= threshold_value)

    def show(self, message, item_list):
        logger.subheading(message)
        for i in item_list:
            logger(i)

    def filter_modules_by_reference(self, filter_modules, reference_modules):

        output_modules = []

        for m in filter_modules:
            found = False
            for m2 in reference_modules:
                if not m2.symmetric_difference(m):
                    found = True
                    break
            if found is False:
                output_modules.append(m)

        return output_modules

    def find_modules(self,
        connectivity,
        comparison_matrix,
        threshold,
        min_neighbours = 3,
        ):

        # How many nodes is each node connected to?
        max_node_connectivity = connectivity.sum(axis=0)

        modules = []

        for i, node_connectivity in enumerate(connectivity):

            # Threshold the comparison matrix values
            sim_mask = self.threshold_function(comparison_matrix[i], threshold).astype(int)
            # Ensure core point is always selected
            sim_mask[i] = 1

            # Copy connectivity matrix
            connectivity_copy = numpy.zeros_like(connectivity).astype(int)
            connectivity_copy[i] = (sim_mask * node_connectivity)

            # Create set of connected nodes
            module_nodes = set([i])
            prev_module_nodes = set([i])
            # Iterate over connectivity to find full network
            while True:
                # A connected node is one with at least one connection
                connection_counts = connectivity_copy.sum(axis=0)
                # Extract connected groups
                connections = set(numpy.where(connection_counts)[0])
                # Find only new connections
                new_connections = list(connections.difference(module_nodes))
                # Break if no new nodes
                if (not new_connections):
                    break
                # Add connections from connected nodes
                module_nodes.update(new_connections)
                for i_new in new_connections:
                    connectivity_copy[i_new] = (sim_mask * connectivity[i_new])

                # Prune nodes with low connectivity
                #
                # len(module_nodes)-1   -> overrides min_neighbours for small modules
                # min_neighbours+1      -> +1 accounts for self connectivity
                #
                curr_min_neighbours = min(len(module_nodes)-1, min_neighbours+1)
                # Count connections for each point
                connection_counts = connectivity_copy.sum(axis=0)
                # Which nodes have fewer required connections?
                prune_bool = (connection_counts < curr_min_neighbours)
                # Unselect any nodes for pruning that have fewer than curr_min_neighbours connections (as impossible to ever fulfil)
                prune_bool[(max_node_connectivity < curr_min_neighbours)] = False
                # Extract nodes for pruning
                prune_nodes = list(numpy.where(prune_bool)[0])

                if prune_nodes:
                    # Remove connections from pruned nodes
                    connectivity_copy[prune_nodes] = 0
                    # Remove from module nodes
                    module_nodes = module_nodes.difference(prune_nodes)
                    #module_nodes.add(i)

                if not module_nodes:
                    break
                # Check convergence
                if not prev_module_nodes.symmetric_difference(module_nodes):
                    break

                prev_module_nodes = copy.deepcopy(module_nodes)

            # Skip singletons
            if (not module_nodes) or (len(module_nodes) <= 1):
                continue

            modules.append(module_nodes)

        return modules

    def create_hierarchy_from_modules(self, modules, resort=False):

        # Flatten modules into single list (maintains order)
        modules = self.filter_duplicate_modules(modules)

        if resort is True:
            self.sort_modules_in_place(modules)

        output_levels = []
        for m1 in modules:
            inserted = False
            for level_modules in output_levels:
                intersects = False
                for m2 in level_modules:
                    if m2.intersection(m1):
                        intersects = True
                if intersects is False:
                    level_modules.append(m1)
                    inserted = True
                    break
            if inserted is False:
                output_levels.append([m1])

        return output_levels

    def filter_subsets(self, module_list):
        modules_out = []
        while module_list:
            m1 = module_list.pop()
            found = False
            for m2 in module_list:
                if m1.issubset(m2) or m2.issubset(m1):
                    m2.update(m1)
                    found = True
                    break
            if found is False:
                modules_out.append(m1)
        self.sort_modules_in_place(modules_out)
        return modules_out

    def filter_complete_sets(self, modules, n_total):
        modules_out = []
        for m in modules:
            if len(m) == n_total:
                continue
            modules_out.append(m)
        return modules_out

    def filter_duplicate_modules(self, modules):
        modules_out = []

        for m1 in modules:
            duplicate = False
            for m2 in modules_out:
                # If has no different elements
                if not m1.symmetric_difference(m2):
                    duplicate = True
                    break
            if duplicate is False:
                modules_out.append(m1)
        return modules_out

    def sets_to_tuples(self, modules_list):

        output_list = []

        for l in modules_list:
            o_l = [tuple(sorted(m)) for m in l]
            output_list.append(o_l)

        return output_list

    def sort_modules_in_place(self, module_list, presort_by_index=False):
        if presort_by_index is True:
            module_list.sort(key=lambda x: min(x))
        module_list.sort(key=lambda x: len(x), reverse=True)

import giant.logs as lg
logger = lg.getLogger(__name__)

from libtbx import adopt_init_args, group_args


class LevelGroupTree(object):


    def __init__(self,
        links,
        cache_sets = True,
        ):

        # Cache for saving results
        set_cache = {}
        # Max number of possible links
        max_links = len(list(links.keys()))

        adopt_init_args(self, locals())

        # Find the roots of the tree
        self.roots = self._find_roots()
        # Reverse the directional links to find parents of each node
        self.node_parents = self._find_node_parents()

    def _find_roots(self):
        """Find roots of the tree (nodes which are not children)"""
        roots = set(self.iter_nodes())
        for l, g in self.iter_nodes():
            l_g_pairs = self.get_children(l,g)
            # Remove anything that is a child
            roots.difference_update(l_g_pairs)
        return sorted(roots)

    def _find_node_parents(self):
        """Find the parents of each node and return as dictionary"""
        # Dictionary of parents for each node
        node_parents = {}
        for l, g in self.iter_nodes():
            for lg_child in self.get_children(l,g):
                node_parents.setdefault(lg_child, []).append((l,g))
        return node_parents

    def iter_nodes(self):
        for l_label, l_links in sorted(self.links.items()):
            for g_label in sorted(l_links.keys()):
                yield (l_label, g_label)

    def iter_nodes_by_level(self):
        for l_label, l_links in sorted(self.links.items()):
            yield (l_label, list(l_links.keys()))

    def is_convergent_node(self, l, g):
        """Check if a node has more than one parent"""
        return (len(self.node_parents.get((l,g),[])) > 1)

    def get_children(self, l, g):
        l_g_links = self.links.get(l,{}).get(g,{})
        children = []
        for l_child, g_children in sorted(l_g_links.items()):
            children.extend([(l_child, g_child) for g_child in g_children])
        return children

    def get_parents(self, l, g):
        return self.node_parents.get((l,g),[])

    def get_children_recursive(self,
        l, g, max_recursions,
        get_parents_for_convergent_nodes = True,
        ):
        """Get the children of a node recursively, going through `max_recursions` number of links"""
        children = []
        # If one recursion, just return the current children
        if max_recursions > 0:
            # Get the direct children of this node
            children += self.get_children(l,g)
            # If more than one recursion, recurse to the level below (get chidren of children)
            if max_recursions > 1:
                for cl, cg in children:
                    # Get children of this node
                    children += self.get_children_recursive(
                        l=cl, g=cg,                         # One "level" lower than original (l,g)
                        max_recursions=max_recursions-1,    # !!! Decrement this as we have moved down a level
                        get_parents_for_convergent_nodes=get_parents_for_convergent_nodes,
                        )
                    # Get parents of this node?
                    if (get_parents_for_convergent_nodes is True) and self.is_convergent_node(cl,cg):
                        parents = self.get_parents(cl,cg)
                        # Should be a convergent node with more than one parent
                        assert len(parents) > 1
                        # Remove the seed parent
                        parents.remove((l,g))
                        # Append to output list
                        children += parents
                        # Recursively add all children
                        for pl, pg in parents:
                            children += self.get_children_recursive(
                                l=pl, g=pg,                     # Same "level" as original (l,g)
                                max_recursions=max_recursions,  # !!! Do not decrement this as we have not descended a level
                                get_parents_for_convergent_nodes=get_parents_for_convergent_nodes,
                                )
        # Make non-redundant and sort
        return sorted(set(children))

    def get_recursive_sets(self, max_recursions):

        # Look for previous solutions if present
        result = self.set_cache.get(max_recursions, None)
        if result is not None:
            return result

        # List of recursions starting from ALL levels
        output_sets = []
        # Iter through groups in each level
        for l, groups in self.iter_nodes_by_level():
            # List of recursions starting from THIS level
            level_sets = []
            # Find recursive sets for each group
            for g in groups:
                # Group starts with itself
                level_group_set = [(l,g)] + self.get_children_recursive(
                        l = l,
                        g = g,
                        max_recursions = max_recursions,
                        get_parents_for_convergent_nodes = True,
                        )
                level_sets.append(level_group_set)
            # Merge overlapping groups (paths may converge)
            output_sets += self.merge_overlapping_groups(level_sets)

        # Cache result?
        if self.cache_sets is True:
            self.set_cache[max_recursions] = output_sets

        return output_sets

    def get_complete_sets(self):

        cache_label = 'complete'

        result = self.set_cache.get(cache_label, None)
        if result is not None:
            return result

        max_recursions = self.max_links

        # List of recursions starting from each root
        root_sets = []
        # Iter through each root
        for l, g in self.roots:
            root_set = [(l,g)] + self.get_children_recursive(
                    l = l,
                    g = g,
                    max_recursions = max_recursions,
                    get_parents_for_convergent_nodes = False, # Do not need to find parents as will be merged later
                    )
            root_sets.append(root_set)

        output_sets = self.merge_overlapping_groups(root_sets)

        # Cache result?
        if self.cache_sets is True:
            self.set_cache[cache_label] = output_sets

        return output_sets

    def merge_overlapping_groups(self, list_of_l_g_pairs):

        input_sets = [set(l_g_pairs) for l_g_pairs in list_of_l_g_pairs]
        output_sets = []

        # Any overlaps found between the sets?
        global_overlaps_found = False

        for i_s in input_sets:
            # local flag for this set
            this_overlaps_found = False
            # Iterate through sets and merge into those with overlaps
            for o_s in output_sets:
                if i_s.intersection(o_s):
                    o_s.update(i_s)
                    this_overlaps_found = True
            # If no overlaps, append to output list
            if this_overlaps_found is False:
                output_sets.append(i_s)
            else:
                global_overlaps_found = True

        # Run recursively until no overlaps are returned
        if global_overlaps_found is True:
            return self.merge_overlapping_groups(output_sets)

        return list(map(sorted, output_sets))


class BuildLevelArrayAsTreeTask(object):


    def __init__(self,
            ):
        """Create the selection array for the supplied levels"""
        adopt_init_args(self, locals())

    def run(self,
            level_group_array,
            ):
        """Identify the tree of groups that form each atom"""

        logger.subheading('Converting selections into hierarchy')

        import numpy

        # Dictionary pointing to "owned" groups in the next level(s)
        # Goes from lower levels to higher levels: graph[chain_level] -> groups on ss_level
        graph = {}

        # Search downwards
        for i_l_above in range(0,level_group_array.shape[0]):
            # Get the groups in the current level
            l_above = level_group_array[i_l_above]
            # Search group by group
            for i_g_above in sorted(numpy.unique(l_above)):
                # Skip in not in a group
                if (i_g_above == -1): continue
                # populate empty to ensure level + groups are definitely known in tree (but connect to nothing)
                graph.setdefault(i_l_above, {}).setdefault(i_g_above, {})
                # Special case -- first level -- populate empty
                if i_l_above == 0:
                    continue
                # Select atoms in this group
                g_sel = (l_above == i_g_above)

                # Find the nearest group "below" (may not be on the adjacent level)
                i_l_below = i_l_above-1
                while i_l_below >= 0:
                    l_below = level_group_array[i_l_below]
                    if (l_below[g_sel] != -1).sum() > 0:
                        break
                    i_l_below -= 1
                if i_l_below < 0:
                    # No group found
                    continue

                # Create dictionary for the lower level, pointing to higher levels
                below_dict = graph.setdefault(i_l_below, {})
                # Iterate through the groups in the lower level
                for i_g_below in sorted(numpy.unique(l_below[g_sel])):
                    if (i_g_below == -1): continue
                    # Initialise a dictionary for this group
                    below_group_dict = below_dict.setdefault(i_g_below, {})
                    # Create list of groups in this level that correspond to higher levels
                    above_list = below_group_dict.setdefault(i_l_above, [])
                    # Add the group
                    above_list.append(i_g_above)

        self.result = group_args(
            tree = LevelGroupTree(links=graph),
            )
        self.show_summary()

        return self.result

    def show_summary(self):

        graph = self.result.tree.links

        logger.subheading('Tree summary of hierarchy:')
        for i in sorted(graph.keys()):
            logger.bar()
            logger('Ownership of level {} (Related groups that will be co-optimised)'.format(i))
            logger.bar()
            for g, vals in graph[i].items():
                groups = ['(Level {}, Groups {})'.format(l, ', '.join(map(str,vs))) for l, vs in vals.items()]
                if not groups: groups = ['none']
                logger('Group {} -> {}'.format(g, ', '.join(groups)))



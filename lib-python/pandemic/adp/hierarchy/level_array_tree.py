from libtbx import adopt_init_args, group_args
from bamboo.common.logs import Log

from pandemic.adp.echt.optimise.inter_level import LevelGroupTree


class BuildLevelArrayAsTreeTask:


    def __init__(self,
            log = None,
            ):
        """Create the selection array for the supplied levels"""
        if log is None: log = Log()
        adopt_init_args(self, locals())

    def run(self,
            level_group_array,
            ):
        """Identify the tree of groups that form each atom"""

        self.log.subheading('Converting selections into hierarchy')

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

        log = self.log
        log.subheading('Tree summary of hierarchy:')
        for i in sorted(graph.keys()):
            log.bar()
            log('Ownership of level {} (Related groups that will be co-optimised)'.format(i))
            log.bar()
            for g, vals in graph[i].items():
                groups = ['(Level {}, Groups {})'.format(l, ', '.join(map(str,vs))) for l, vs in vals.items()]
                if not groups: groups = ['none']
                log('Group {} -> {}'.format(g, ', '.join(groups)))



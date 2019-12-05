from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry
from bamboo.common.logs import Log

import numpy


class BuildLevelArrayTask:


    def __init__(self,
        overall_selection = None,
        remove_duplicate_groups = None,
        warnings = None,
        log = None,
        ):
        """Create the selection array for the supplied levels"""
        if log is None: log = Log()
        adopt_init_args(self, locals())

    def run(self,
        hierarchy,
        level_group_selection_strings,
        ):

        level_group_array = self.generate_array(
            hierarchy = hierarchy,
            selection_strings = level_group_selection_strings,
            )

        self.validate_array_and_selections(
            level_group_array = level_group_array,
            selection_strings = level_group_selection_strings,
            assert_complete = False,
            )

        if (self.remove_duplicate_groups is not None):
            self.remove_duplicate_groups_from_array(
                level_group_array = level_group_array,
                keep_which_group = self.remove_duplicate_groups,
                )

        # Check for empty levels
        invalid_levels = [i+1 for i,l in enumerate(level_group_array) if (l == -1).all()]
        if len(invalid_levels) > 0:
            raise Sorry('Levels have been created that contain no atoms!\nInvalid Levels: {}'.format(invalid_levels))

        # Filter the selection strings to those actually used (e.g. not masked by overall_selection)
        filtered_selections = self.filter_selection_strings_by_array(
            level_group_array = level_group_array,
            selection_strings = level_group_selection_strings,
            )

        # Renumber so each level starts from 1
        level_group_array = self.renumber_array(level_group_array)

        # Find the actual global selection bool
        overall_atom_mask = (level_group_array!=-1).any(axis=0)

        # Filter and store array by overall mask
        level_group_array = level_group_array[:,overall_atom_mask]

        self.validate_array_and_selections(
            level_group_array = level_group_array,
            selection_strings = filtered_selections,
            overall_mask = overall_atom_mask,
            assert_complete = True,
            )

        self.result = group_args(
            level_group_array = level_group_array,
            level_group_selection_strings = filtered_selections,
            overall_atom_mask = overall_atom_mask,
            )
        self.show_summary()

        return self.result

    def validate_array_and_selections(self,
        level_group_array,
        selection_strings,
        overall_mask = None,
        assert_complete = False,
        ):

        assert len(level_group_array) == len(selection_strings)

        if overall_mask is not None:
            n_atoms = sum(overall_mask)

        for level_values, level_strings in zip(level_group_array, selection_strings):

            if overall_mask is not None:
                assert len(level_values) == n_atoms

            assert min(level_values) >= -1
            assert max(level_values) < len(level_strings)

            if assert_complete is True:
                assert range(len(level_strings)) == sorted(set(level_values).difference({-1}))

    def generate_array(self,
        hierarchy,
        selection_strings,
        ):

        log = self.log
        log.subheading('Converting selections to array')

        # Extract the atoms for each tls group
        atom_cache = hierarchy.atom_selection_cache()

        # Array of which atoms are in which group at which level
        idx_array = -1 * numpy.ones((len(selection_strings), hierarchy.atoms_size()), dtype=int)

        # Apply the overall filter if provided
        if self.overall_selection is not None:
            global_selection = numpy.array(atom_cache.selection(self.overall_selection), dtype=bool)
        else:
            global_selection = numpy.ones(hierarchy.atoms().size(), dtype=bool)

        # List of any selections that result in no atoms
        warnings = []
        errors = []

        # Convert selection strings to boolean selections
        for i_level, group_selections in enumerate(selection_strings):
            # Groups with no atoms
            no_atoms = 0
            # Iterate through selections and create array
            for i_group, group_sel_str in enumerate(group_selections):
                # Create boolean selection
                group_sel_bool = numpy.array(atom_cache.selection(group_sel_str), dtype=bool)
                # Print warning if there are no atoms selected
                if sum(group_sel_bool) == 0:
                    warnings.append('Level {}, Group "{}": selects no atoms ({} atoms)'.format(i_level+1, group_sel_str, sum(group_sel_bool)))
                    no_atoms += 1
                    continue
                elif sum(group_sel_bool[global_selection]) == 0:
                    warnings.append('Level {}, Group "{}" is omitted by the global overall mask'.format(i_level+1, group_sel_str))
                    no_atoms += 1
                    continue
                # Check that atoms are not already allocated to another group
                if (idx_array[i_level, group_sel_bool] != -1).any():
                    for i_other in numpy.unique(idx_array[i_level, group_sel_bool]):
                        if i_other == -1: continue
                        errors.append('Selection {} and {} overlap in Level {} (atoms can only be in one group per level)'.format(group_selections[i_other], group_sel_str))
                # Set level to this group
                idx_array[i_level, group_sel_bool] = i_group
            # Sanity check that the expected number of groups have been produced
            n_groups = len(set(idx_array[i_level]).difference({-1}))
            assert n_groups+no_atoms == len(group_selections)

        # Apply the global mask (as technically not restrictive above)
        if self.overall_selection is not None:
            ivt_selection = numpy.logical_not(global_selection)
            idx_array[:,ivt_selection] = -1

        # Raise errors
        if errors:
            for e in errors: log(e)
            raise Sorry('Errors raised during level generation. See above messages.')

        # List warnings
        if warnings:
            msg = 'WARNING: One or more group selections do not select any atoms: \n\t{}'.format('\n\t'.join(warnings))
            if self.warnings is not None:
                self.warnings.append(msg)
                self.warnings.flush()
            else:
                log(msg)

        return idx_array

    def remove_duplicate_groups_from_array(self,
        level_group_array,
        keep_which_group,
        ):

        assert keep_which_group in ['keep_highest_group', 'keep_lowest_group']

        report_strings = []

        for i_level, i_level_values in enumerate(level_group_array):

            for i_group in sorted(set(i_level_values)):

                if i_group == -1:
                    continue

                # Get the selection for this group in this level
                i_group_sel = (i_level_values == i_group)

                for j_level, j_level_values in enumerate(level_group_array[i_level+1:]):

                    # Need to apply offset because of slicing of array
                    j_level = j_level + (i_level + 1)

                    # Get the j_values for i_group in j_level
                    j_group_values = j_level_values[i_group_sel]

                    # If more than one index in the group on the j_level then cannot be the same
                    if len(set(j_group_values)) > 1:
                        continue

                    # Get the j_index for the j_group in the j_level
                    j_group = j_group_values[0]

                    # Get the full selection for the j_group on the j_level
                    j_group_sel = (j_level_values == j_group)

                    # If the selections are the same then these are the same group
                    duplicate = numpy.array_equal(i_group_sel, j_group_sel)

                    # Duplicate groups!
                    if duplicate is True:
                        if keep_which_group == 'keep_lowest_group':
                            # Remove the group at the higher level
                            level_group_array[j_level, j_group_sel] = -1
                            # Report
                            report_strings.append('Removing Group {} from Level {} (also exists as Group {} on Level {})'.format(j_group+1, j_level+1, i_group+1, i_level+1))
                        else:
                            # Remove the group at the lower level
                            level_group_array[i_level, i_group_sel] = -1
                            # Report
                            report_strings.append('Removing Group {} from Level {} (also exists as Group {} on Level {})'.format(i_group+1, i_level+1, j_group+1, j_level+1))
                        # Don't need to check any further for this group
                        break

        if report_strings:
            msg = 'Removed {} duplicated groups from hierarchy: \n\t{}'.format(len(report_strings), '\n\t'.join(report_strings))
            if self.warnings is not None:
                self.warnings.append(msg)
                self.warnings.flush()
            else:
                self.log(msg)

        return level_group_array

    def filter_selection_strings_by_array(self,
        level_group_array,
        selection_strings,
        ):

        new_selection_strings = [list() for i in selection_strings]

        for i_level, level_values in enumerate(level_group_array):

            # Iterate through groups in the level
            for i_group in sorted(set(level_values)):

                # Skip null group
                if (i_group == -1): continue

                # Transfer the appropriate string
                new_selection_strings[i_level].append(selection_strings[i_level][i_group])

        return new_selection_strings

    def renumber_array(self,
        level_group_array,
        ):

        # Iterate through levels
        for i_level, level_values in enumerate(level_group_array):
            # Counter for the new groups
            i_group_new = 0
            # Iterate through the old group indices in order
            for i_group_old in sorted(set(level_values)):
                # Don't care about the null group
                if i_group_old == -1: continue
                # Get the selection for the group
                i_group_sel = (level_values == i_group_old)
                # Assign the new group number - in place!
                level_group_array[i_level, i_group_sel] = i_group_new
                # Increment
                i_group_new += 1

        return level_group_array

    def show_summary(self):

        log = self.log
        n_levels, n_atoms = self.result.level_group_array.shape
        log('> Constructed {} levels containing {} atoms'.format(n_levels, n_atoms))
        for i_level, level_values in enumerate(self.result.level_group_array):
            group_nums, group_counts = numpy.unique(level_values, return_counts=True)
            log.bar(True, True)
            log('> Level {} - covers {}/{} atoms -- partitioned into {} groups\n'.format(i_level+1, sum(level_values!=-1), n_atoms, len(group_nums)))
            for i_grp, count in zip(group_nums, group_counts):
                if i_grp == -1: continue
                label = self.result.level_group_selection_strings[i_level][i_grp]
                log('> Group {:>5d}: {:>5d} atoms ({})'.format(i_grp+1, count, label))
        log.bar(True, False)



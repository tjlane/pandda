import giant.logs as lg
logger = lg.getLogger(__name__)

from libtbx import adopt_init_args, group_args
from libtbx import easy_mp
from libtbx.utils import Sorry

import numpy


class SelectionStringConverter(object):


    def __init__(self,
        atom_cache,
        global_selection,
        ):
        adopt_init_args(self, locals())

    def __call__(self, selection_string):
        # Create boolean selection
        sel_bool = numpy.array(self.atom_cache.selection(selection_string), dtype=bool)
        # Print warning if there are no atoms selected
        if sum(sel_bool) == 0:
            return '"{}" selects no atoms ({} atoms)'.format(selection_string, sum(sel_bool))
        elif sum(sel_bool[self.global_selection]) == 0:
            return '"{}" is omitted by the global overall mask'.format(selection_string)

        return sel_bool


class BuildLevelArrayTask(object):


    def __init__(self,
        overall_selection = None,
        remove_duplicate_groups = None,
        n_cpus = 1,
        ):
        """Create the selection array for the supplied levels"""
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
                assert list(range(len(level_strings))) == sorted(set(level_values).difference({-1}))

    def generate_array(self,
        hierarchy,
        selection_strings,
        ):

        logger.subheading('Converting selections to array')

        # Array of which atoms are in which group at which level
        idx_array = -1 * numpy.ones((len(selection_strings), hierarchy.atoms_size()), dtype=int)

        # Extract the atoms for each tls group
        atom_cache = hierarchy.atom_selection_cache()

        # Apply the overall filter if provided
        if self.overall_selection is not None:
            global_selection = numpy.array(atom_cache.selection(self.overall_selection), dtype=bool)
        else:
            global_selection = numpy.ones(hierarchy.atoms().size(), dtype=bool)

        # Function for converting selection strings to boolean selections
        extract_selection = SelectionStringConverter(
            atom_cache = atom_cache,
            global_selection = global_selection,
            )

        # Create multi-process wrapper
        from pandemic.adp.parallel import RunParallelWithProgressBarUnordered
        extract_selections_parallel = RunParallelWithProgressBarUnordered(
            function = extract_selection,
            n_cpus = self.n_cpus,
            max_chunksize = 100,
        )

        # Extract selection strings
        arg_dicts = []
        for i_level, group_selections in enumerate(selection_strings):
            # Iterate through selections and create array
            for i_group, group_sel_str in enumerate(group_selections):
                arg_dicts.append(dict(selection_string=group_sel_str))

        # Convert selection strings to boolean selections
        logger('Converting selection strings to boolean selections')
        results = extract_selections_parallel(arg_dicts=arg_dicts)

        # List of any selections that result in no atoms
        errors = []; warnings = []

        # Extract results and place in overall array
        for i_level, group_selections in enumerate(selection_strings):
            # Groups with no atoms/other error
            no_sel = 0
            # Iterate through selections and create array
            for i_group, group_sel_str in enumerate(group_selections):
                # Extract from mp -- results should be in order
                group_sel_bool = results.pop(0)
                # Check if error extracting selection
                if isinstance(group_sel_bool, str):
                    # warning message is the result
                    warning_msg = 'Level {level}, Group {group}: '.format(
                        level = i_level+1,
                        group = i_group+1,
                        ) + str(group_sel_bool)
                    warnings.append(warning_msg)
                    # mark this group as not selecting anything
                    no_sel += 1
                    continue
                # Check that atoms are not already allocated to another group
                if (idx_array[i_level, group_sel_bool] != -1).any():
                    # Check if these atoms already assigned to another group
                    for i_other in numpy.unique(idx_array[i_level, group_sel_bool]):
                        # If == -1, not overlapping
                        if i_other == -1:
                            continue
                        error_msg = 'Selection "{selection1}" and "{selection2}" overlap in Level {level} (atoms can only be in one group per level)'.format(
                            selection1 = group_selections[i_other],
                            selection2 = group_sel_str,
                            level = i_level+1,
                            )
                        errors.append(error_msg)
                # Set level to this group
                idx_array[i_level, group_sel_bool] = i_group
            # Sanity check that the expected number of groups have been produced
            n_groups = len(set(idx_array[i_level]).difference({-1}))
            assert n_groups+no_sel == len(group_selections)

        # Apply the global mask (as technically not restrictive above)
        if self.overall_selection is not None:
            ivt_selection = numpy.logical_not(global_selection)
            idx_array[:,ivt_selection] = -1

        # List warnings
        if len(warnings) > 0:
            msg = 'One or more group selections do not select any atoms: \n\t{}'.format('\n\t'.join(warnings))
            logger.warning(msg)

        # Raise errors
        if len(errors) > 0:
            msg = 'Errors raised during level generation: \n\t{}'.format('\n\t'.join(errors))
            raise Sorry(msg)

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

        if len(report_strings) > 0:
            msg = 'Removed {} duplicated groups from hierarchy: \n\t{}'.format(len(report_strings), '\n\t'.join(report_strings))
            logger.warning(msg)

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

        n_levels, n_atoms = self.result.level_group_array.shape
        logger('> Constructed {} levels containing {} atoms'.format(n_levels, n_atoms))
        for i_level, level_values in enumerate(self.result.level_group_array):
            group_nums, group_counts = numpy.unique(level_values, return_counts=True)
            logger.bar(True, True)
            logger('> Level {} - covers {}/{} atoms -- partitioned into {} groups\n'.format(i_level+1, sum(level_values!=-1), n_atoms, len(group_nums)))
            for i_grp, count in zip(group_nums, group_counts):
                if i_grp == -1: continue
                label = self.result.level_group_selection_strings[i_level][i_grp]
                logger('> Group {:>5d}: {:>5d} atoms ({})'.format(i_grp+1, count, label))
        logger.bar(True, False)



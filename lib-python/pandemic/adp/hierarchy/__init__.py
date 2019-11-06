import os, glob, collections
import numpy

from bamboo.common.logs import Log
from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure

from giant.structure.pymol import auto_residue_images, auto_chain_images, selection_images

def translate_phenix_selections_to_pymol_selections_simple(selections, verbose=False):
    """Convert simplex phenix selections into simplex pymol selections"""

    easy_regexes = [
            ("chain '?[A-Za-z]{1,2}'?",
                lambda x: x
                ),
            ("chain '?[A-Za-z]{1,2}'? and resid '?[0-9]*?[A-Z]?'? through '?[0-9]*[A-Z]?'?",
                lambda x: x.replace("resid", "resi").replace(" through ",":")
                ),
                    ]

    from bamboo.common.logs import ScreenLogger
    log = ScreenLogger(stdout=verbose)

    import re
    output_selections = []
    for s in selections:
        # Start
        o = None
        # Remove double spaces and padding
        while '  ' in s:
            s = s.replace('  ',' ')
        s = s.strip(' ')
        log('Phenix Selection String: {}'.format(s))
        for rgx_s, trans_func in easy_regexes:
            rgx = re.compile(rgx_s)
            mat = rgx.findall(s)
            # No matches or too many matches
            if len(mat) != 1:
                log('> No matches or multiple matches to {}'.format(rgx_s))
                continue
            # If the match is the same as input string, then process
            if mat[0] == s:
                log('Exact match for {}'.format(rgx_s))
                o = trans_func(s)
                log('Translating to pymol: \n\t{} -> {}'.format(s, o))
                break
        # Append processed selection or none
        output_selections.append(o)

    assert len(output_selections) == len(selections)
    return output_selections


class GenerateLevelSelectionsTask:


    def __init__(self,
            auto_levels = [],
            custom_levels = None,
            overall_selection = None,
            cbeta_in_backbone = True,
            verbose = False,
            log = None,
            ):
        if log is None: log = Log()
        adopt_init_args(self, locals())

    def run(self,
        hierarchy,
        ):
        """Build the levels for the hierarchical fitting"""

        from giant.structure.formatting import PhenixSelection
        from giant.structure.select import \
            protein, backbone, sidechains, \
            default_secondary_structure_selections_filled

        log = self.log

        log.subheading('Building selections for levels')
        levels = []; labels=[];

        filter_h = hierarchy
        cache = filter_h.atom_selection_cache()

        log('Input hierarchy contains {} atoms'.format(filter_h.atoms().size()))

        # Filter the hierarchy by the overall selection
        if self.overall_selection:
            choice = cache.selection(self.overall_selection)
            filter_h = filter_h.select(choice, copy_atoms=True)
            log('Overall selection ({}) selects {} atoms'.format(self.overall_selection, filter_h.atoms().size()))

        # Whether cbeta is in backbone
        cbeta_flag = self.cbeta_in_backbone
        backbone_atoms = ['C','CA','N','O']+(['CB']*cbeta_flag)
        backbone_atoms_sel = '({})'.format(' or '.join(['name {}'.format(a) for a in backbone_atoms]))
        back_sel = ' and '+backbone_atoms_sel
        side_sel = ' and not '+backbone_atoms_sel

        log.bar(True, False)
        log('Creating automatic levels:')

        if 'chain' in self.auto_levels:
            log('Level {}: Creating level with groups for each chain'.format(len(levels)+1))
            groups = [PhenixSelection.format(c) for c in filter_h.chains()]
            levels.append(sorted(set(groups))) # Chains can be present multiple times
            labels.append('chain')
        if 'auto_group' in self.auto_levels:
            from giant.structure.tls import phenix_find_tls_groups
            log('Level {}: Creating level with groups determined by phenix.find_tls_groups'.format(len(levels)+1))
            groups = [s.strip('"') for s in phenix_find_tls_groups(hierarchy=hierarchy)]
            levels.append([g for g in groups if not cache.selection(g).all_eq(False)])
            labels.append('groups')
        if ('secondary_structure' in self.auto_levels) or ('ss' in self.auto_levels):
            log('Level {}: Creating level with groups based on secondary structure'.format(len(levels)+1))
            groups = [s.strip('"') for s in default_secondary_structure_selections_filled(hierarchy=filter_h, verbose=self.verbose)]
            levels.append([g for g in groups if not cache.selection(g).all_eq(False)])
            labels.append('sec. struct.')
        if 'residue' in self.auto_levels:
            log('Level {}: Creating level with groups for each residue'.format(len(levels)+1))
            levels.append([PhenixSelection.format(r) for r in filter_h.residue_groups()])
            labels.append('residue')
        if 'backbone_sidechain' in self.auto_levels:
            log('Level {}: Creating level with groups for each residue backbone/sidechain'.format(len(levels)+1))
            b_gps = backbone(filter_h, cbeta=cbeta_flag).atom_groups()
            b_sels = [PhenixSelection.format(r)+back_sel for r in b_gps if (r.resname not in ['ALA','GLY','PRO'])]
            s_gps = sidechains(filter_h, cbeta=(not cbeta_flag)).atom_groups()
            s_sels = [PhenixSelection.format(r)+side_sel for r in s_gps if (r.resname not in ['ALA','GLY','PRO'])]
            levels.append(sorted(b_sels+s_sels))
            labels.append('backbone/sidechain')
        if 'atom' in self.auto_levels:
            log('Level {}: Creating level with groups for each atom'.format(len(levels)+1))
            levels.append([PhenixSelection.format(a) for a in filter_h.atoms()])
            labels.append('atom')
        log.bar()

        # Insert custom levels
        if self.custom_levels:
            # Print auto levels
            log('> {} automatic levels created:'.format(len(levels)))
            for i_l, level in enumerate(levels):
                log('\tLevel {} ({}) - {} groups'.format(i_l+1, labels[i_l], len(level)))
            log.bar()
            # Insert custom levels
            log.bar(True, False)
            log('Inserting custom levels:')
            for l_params in self.custom_levels:

                # Skip blank levels that might be inserted
                if (l_params.selection == []):
                    continue

                # Only ONE can be given
                if [l_params.depth, l_params.insert_before, l_params.insert_after].count(None) < 2:
                    msg = ""
                    msg += "For each custom level, you must define EITHER depth OR insert_before OR insert_after (OR none of them). "
                    msg += "\nLevel {} is currently defined with:".format(l_params.label)
                    msg += "\n\tdepth:         {}".format(l_params.depth)
                    msg += "\n\tinsert_before: {}".format(l_params.insert_before)
                    msg += "\n\tinsert_after:  {}".format(l_params.insert_after)
                    raise Sorry(msg)

                # label MUST be given (?)
                if l_params.label is None:
                    raise Sorry("Must provide label for each custom level")

                # List index to insert level
                if l_params.depth is not None:
                    if l_params.depth < 1:
                        msg = 'Custom level depths cannot be less that 1! (input depth is {} for level with label {})'.format(l_params.depth, l_params.label)
                        raise Sorry(msg)
                    if l_params.depth > (len(levels) + 1):
                        msg = 'Trying to add group "{label}" at position {depth}, but only {length} levels currently exist!'.format(
                            label=l_params.label,
                            depth=l_params.depth,
                            length=len(levels),
                            )
                        msg += '\nThe maximum possible depth that can be added at this point is {} given the current levels (shown below):'.format(len(levels)+1)
                        msg += '\nYou might try changing the depth of this level or adding groups in a different order?'
                        msg += '\nLevels added until this point: \n\t{}'.format('\n\t'.join(['Depth {}: {}'.format(i+1, l) for i,l in enumerate(labels)]))
                        raise Sorry(msg)
                    # insertion index is just the depth - 1
                    idx = l_params.depth - 1
                elif l_params.insert_before is not None:
                    if l_params.insert_before not in labels:
                        msg = 'Trying to insert level "{new}" before level "{ref}" but "{ref}" does not exist yet!'.format(new=l_params.label, ref=l_params.insert_before)
                        msg += '\nYou might try adding groups in a different order? (the reference group must be added before this group).'
                        raise Sorry(msg)
                    # insertion index is where the reference group is (will be inserted at this point, shifting reference down)
                    idx = labels.index(l_params.insert_before)
                elif l_params.insert_after is not None:
                    if l_params.insert_after not in labels:
                        msg = 'Trying to insert level "{new}" after level "{ref}" but "{ref}" does not exist yet!'.format(new=l_params.label, ref=l_params.insert_after)
                        msg += '\nYou might try adding groups in a different order? (the reference group must be added before this group).'
                        raise Sorry(msg)
                    # insertion index is one after the reference group
                    idx = labels.index(l_params.insert_after) + 1
                else:
                    # Just append to the hierarchy
                    idx = len(levels)
                    #raise Failure('Should not reach here! Please contact developer (send log file!)')

                log('Inserting level: \n\tLabel: {label}\n\tPosition: {pos}\n\tNumber of groups: {n_groups}'.format(
                    label = l_params.label,
                    pos = idx+1,
                    n_groups = len(l_params.selection),
                    ))
                assert len(l_params.selection) > 0, 'No selections provided for this group!'
                for g in l_params.selection:
                    if cache.selection(g).all_eq(False):
                        raise Sorry('Selection "{}" does not select any atoms'.format(g))
                levels.insert(idx, l_params.selection)
                labels.insert(idx, l_params.label)
            log.bar()

        # Report
        log.subheading('Hierarchy summary: {} levels created'.format(len(levels)))
        for i_l, level in enumerate(levels):
            log.bar()
            log('Level {} ({}) - {} groups'.format(i_l+1, labels[i_l], len(level)))
            log.bar()
            for i, l in enumerate(level): log('\t{:>5d} : {}'.format(i+1,l))

        self.result = group_args(
            level_group_selection_strings = levels,
            level_labels = labels,
            )
        return self.result


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
                    #raise Exception('Broken! Couldn\'t find a group in a lower level')

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
            graph = graph,
            )
        self.show_summary()

        return self.result

    def show_summary(self):

        graph = self.result.graph

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


class CreateHierarchicalModelTask:


    def __init__(self,
        auto_levels,
        custom_levels,
        overall_selection,
        cbeta_in_backbone,
        remove_duplicate_groups = False,
        warnings = None,
        verbose = False,
        log = None,
        ):
        if log is None: log = Log()
        adopt_init_args(self, locals())

        # Selection strings for each group for each level
        self.selections_constructor = GenerateLevelSelectionsTask(
            auto_levels = self.auto_levels,
            custom_levels = self.custom_levels,
            overall_selection = self.overall_selection,
            cbeta_in_backbone = self.cbeta_in_backbone,
            log = self.log,
            )

        # Construct 2-d integer array of group indices for each atom on each level
        self.array_constructor = BuildLevelArrayTask(
            overall_selection = self.overall_selection,
            remove_duplicate_groups = self.remove_duplicate_groups,
            warnings = self.warnings,
            log = self.log,
            )

        self.array_as_tree = BuildLevelArrayAsTreeTask(
            log = self.log,
            )

    def run(self,
        hierarchy,
        ):

        self.log.heading('Generating hierarchical model')

        selections = self.selections_constructor.run(
            hierarchy = hierarchy,
            )

        arrays = self.array_constructor.run(
            hierarchy = hierarchy,
            level_group_selection_strings = selections.level_group_selection_strings,
            )

        graphs = self.array_as_tree.run(
            level_group_array = arrays.level_group_array,
            )

        self.result = group_args(
            level_labels                    = selections.level_labels,
            level_group_array               = arrays.level_group_array,
            level_group_selection_strings   = arrays.level_group_selection_strings,
            level_group_tree                = graphs.graph,
            overall_atom_mask               = arrays.overall_atom_mask,
            )

        return self.result


class WriteHierarchySummaryTask:

    debug = False

    level_atoms_pdb = 'level_{:04d}-partitions.pdb'
    level_atoms_pymol_png = 'level_{:04d}-partitions-chain_{}.png'
    level_partitions_png = 'all-level-partitions-chain_{}.png'

    pymol_script_py = 'pymol_script.py'

    def __init__(self,
        output_directory,
        pymol_images = None,
        warnings = None,
        verbose = False,
        log = None,
        ):
        if log is None: log = Log()
        adopt_init_args(self, locals())

    def show(self, file_dict, indent=0):
        log = self.log
        s = '  '
        for k, v in file_dict.iteritems():
            if isinstance(v, dict):
                log(s*indent + '> {}'.format(k))
                self.show(v, indent+1)
            elif isinstance(v, str):
                log(s*indent + '> {}: {}'.format(k, v))
            else:
                log(s*indent + '> {}'.format(k))
                for vv in v:
                    log(s*(indent+1)+vv)

    def filepath(self, filename):
        return os.path.join(self.output_directory, filename)

    def run(self,
        reference_hierarchy,
        level_group_array,
        level_group_selection_strings,
        level_labels,
        overall_atom_mask,
        plotting_object,
        ):
        """Write out the composition of the hierarchical model"""

        log = self.log

        log.subheading('Writing summary of hierarchical model')

        # Object to generate structures with different b-factors etc.
        from pandemic.adp.utils import StructureFactory
        s_fac = StructureFactory(master_h=reference_hierarchy)

        output_files = collections.OrderedDict()

        of, level_hierarchies = self.write_level_atoms_structures(
            level_group_array = level_group_array,
            level_labels = level_labels,
            overall_atom_mask = overall_atom_mask,
            structure_factory = s_fac,
            )
        self.show(of)
        output_files.update(of)

        of = self.make_level_atoms_plots(
            level_hierarchies = level_hierarchies,
            level_labels = level_labels,
            structure_factory = s_fac,
            plotting_object = plotting_object,
            )
        self.show(of)
        output_files.update(of)

        # Generate images of each chain of each level coloured by group
        if self.pymol_images is not None:
            of = self.make_pymol_images_of_level_atoms(
                level_atoms_pdb = output_files['level_atoms_pdb'],
                level_group_selection_strings = level_group_selection_strings,
                level_labels = level_labels,
                structure_factory = s_fac,
                )
            self.show(of)
            output_files.update(of)

        of = {'pymol_script' : self.make_pymol_script(file_dict=output_files)}
        self.show(of)
        output_files.update(of)

        self.result = group_args(
            output_files = output_files,
            )

        return self.result

    def write_level_atoms_structures(self,
        level_group_array,
        level_labels,
        overall_atom_mask,
        structure_factory,
        ):

        file_dict = collections.OrderedDict()

        # Generate hierarchy for each level with groups as b-factors
        level_hierarchies = []
        for i_level, g_vals in enumerate(level_group_array):
            level_lab = level_labels[i_level]
            all_values = -1 * numpy.ones_like(overall_atom_mask)
            all_values[overall_atom_mask] = g_vals # initially indexed from zero!
            m_h = structure_factory.custom_copy(uij=None, iso=all_values, mask=None, blank_copy=True)
            m_f = self.filepath(self.level_atoms_pdb.format(i_level+1))
            m_h.write_pdb_file(m_f)
            file_dict.setdefault('level_atoms_pdb',collections.OrderedDict())[level_lab] = m_f
            # Append to list for plotting
            level_hierarchies.append(m_h)

        return file_dict, level_hierarchies

    def make_level_atoms_plots(self,
        level_hierarchies,
        level_labels,
        structure_factory,
        plotting_object,
        ):

        file_dict = collections.OrderedDict()

        # Write hierarchy plot for each chain
        b_h = structure_factory.blank_copy()
        b_c = b_h.atom_selection_cache()
        for c in b_h.chains():
            chain_sel = b_c.selection('chain {}'.format(c.id))
            hierarchies = [h.select(chain_sel, copy_atoms=True) for h in level_hierarchies]
            # Skip if no partitions in this chain
            #if (numpy.array([h.atoms().extract_b() for h in hierarchies]) == -1).all():
            #    continue
            filename = self.filepath(self.level_partitions_png.format(c.id))
            plotting_object.level_plots(
                filename=filename,
                hierarchies=hierarchies,
                labels = ['Level {}\n({})'.format(i_l+1, l) for i_l, l in enumerate(level_labels)],
                title='chain {}'.format(c.id),
                )
            file_dict.setdefault('level_partitions_png',collections.OrderedDict())[c.id] = filename

        return file_dict

    def make_pymol_images_of_level_atoms(self,
        level_atoms_pdb,
        level_group_selection_strings,
        level_labels,
        structure_factory,
        ):

        file_dict = collections.OrderedDict()

        from giant.structure.formatting import PymolSelection
        for level_lab, structure_filename in level_atoms_pdb.iteritems():
            i_level = level_labels.index(level_lab)
            # Images for each chain (of the partitions) - coloured by group
            filenames_glob   = self.filepath(self.level_atoms_pymol_png.format(i_level+1, '*'))
            filenames_prefix = filenames_glob.replace('-chain_*.png','')
            # Choose the style based on whether interested in atoms or larger regions
            styles = 'cartoon' if level_lab in ['chain','groups','sec. struct.'] else 'lines+spheres'
            # Create selections for each group in each level
            pymol_selections = translate_phenix_selections_to_pymol_selections_simple(level_group_selection_strings[i_level])
            # If returned selection is none, create atom-by-atom pymol selection
            blank_h = structure_factory.blank_copy()
            cache_h = blank_h.atom_selection_cache()
            selections = [s1 if (s1 is not None) else PymolSelection.join_or([PymolSelection.format(a) for a in blank_h.atoms().select(cache_h.selection(s2))]) for s1,s2 in zip(pymol_selections, level_group_selection_strings[i_level])]
            # Create image of different selections
            auto_chain_images(structure_filename = structure_filename,
                              output_prefix = filenames_prefix,
                              style = styles,
                              het_style = 'lines+spheres',
                              colours = ['red','green','blue'],
                              colour_selections = selections,
                              settings = [('cartoon_oval_length', 0.5),
                                          ('cartoon_discrete_colors', 'on'),
                                          ('sphere_scale', 0.25)],
                              width=1000, height=750,
                              delete_script = (not self.debug))
            output_images = glob.glob(filenames_glob)
            if not output_images:
                self.warnings.append('no plots have been generated! ({})'.format(filenames_glob))
            else:
                # Identify the chain for each image
                output_hash = collections.OrderedDict([(s.split('_')[-1].replace('.png',''), s) for s in output_images])
                # Store in output dictionary
                file_dict.setdefault('level_atoms_pymol_png', collections.OrderedDict())[level_lab] = output_hash

        return file_dict

    def make_pymol_script(self, file_dict):

        from bamboo.pymol_utils import PythonScript

        s = PythonScript(pretty_models=False, pretty_maps=False)

        s.change_into_directory(path=os.path.abspath(self.output_directory))

        for f in file_dict.get('level_atoms_pdb',{}).values():
            if f is None: continue
            obj = os.path.basename(f)
            s.load_pdb(
                f_name = os.path.relpath(os.path.abspath(f), start=self.output_directory), 
                obj = os.path.basename(f),
                )
            s.colour(obj=obj, colour='grey')
            s.custom('spectrum', expression='b%10', palette="blue_white_green", selection="{} and (b>-1)".format(obj))
        
        s.show_as(obj='all', style='spheres')
        s.show(obj='all', style='sticks')
        s.set('sphere_scale', 0.25)
        s.set('grid_mode', 1)
        s.orient(obj='all')

        filename = self.filepath(self.pymol_script_py)
        s.write_script(filename)

        return filename



import giant.logs as lg
logger = lg.getLogger(__name__)

from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure


class GenerateLevelSelectionsTask(object):


    def __init__(self,
            auto_levels = [],
            custom_levels = None,
            overall_selection = None,
            cbeta_in_backbone = True,
            assign_het_residues_to_nearest_ss_groups = True,
            assign_het_residues_to_nearest_custom_groups = True,
            assignment_distance_cutoff = 10,
            ):
        adopt_init_args(self, locals())

    def run(self,
        hierarchy,
        ):
        """Build the levels for the hierarchical fitting"""

        from giant.structure.formatting import PhenixSelection
        from giant.structure.select import \
            protein, backbone, sidechains, \
            default_secondary_structure_selections_filled

        logger.subheading('Building selections for levels')
        levels = []; labels=[];

        filter_h = hierarchy
        cache = filter_h.atom_selection_cache()

        logger('Input hierarchy contains {} atoms'.format(filter_h.atoms().size()))

        # Filter the hierarchy by the overall selection
        if self.overall_selection:
            choice = cache.selection(self.overall_selection)
            filter_h = filter_h.select(choice, copy_atoms=True)
            logger('Overall selection ({}) selects {} atoms'.format(self.overall_selection, filter_h.atoms().size()))

        # Whether cbeta is in backbone
        cbeta_flag = self.cbeta_in_backbone
        backbone_atoms = ['C','CA','N','O']+(['CB']*cbeta_flag)
        backbone_atoms_sel = '({})'.format(' or '.join(['name {}'.format(a) for a in backbone_atoms]))
        back_sel = ' and '+backbone_atoms_sel
        side_sel = ' and not '+backbone_atoms_sel

        logger.bar(True, False)
        logger('Creating automatic levels:')

        if 'chain' in self.auto_levels:
            logger('Level {}: Creating level with groups for each chain'.format(len(levels)+1))
            groups = [PhenixSelection.format(c) for c in filter_h.chains()]
            levels.append(sorted(set(groups))) # Chains can be present multiple times
            labels.append('chain')

        if 'phenix_find_tls_groups' in self.auto_levels:
            from giant.structure.tls import phenix_find_tls_groups
            logger('Level {}: Creating level with groups determined by phenix.find_tls_groups'.format(len(levels)+1))
            unfiltered_groups = [s.strip('"') for s in phenix_find_tls_groups(hierarchy=filter_h)]
            groups = [g for g in unfiltered_groups if not cache.selection(g).all_eq(False)]
            levels.append(groups)
            labels.append('groups')

        if ('secondary_structure' in self.auto_levels) or ('ss' in self.auto_levels):
            logger('Level {}: Creating level with groups based on secondary structure'.format(len(levels)+1))
            try:
                unfiltered_groups = [s.strip('"') for s in default_secondary_structure_selections_filled(hierarchy=filter_h)]
            except Exception as e:
                import traceback
                logger.debug(traceback.format_exc())
                logger.warning('\nError during secondary structure identification: {}\n'.format(str(e)))
                raise Sorry('DSSP algorithm failed to identify secondary structure elements -- you will have to provide secondary structure selections manually as a custom level')
            groups = [g for g in unfiltered_groups if not cache.selection(g).all_eq(False)]
            # Assign het molecules to each ss group?
            if (self.assign_het_residues_to_nearest_ss_groups is True):
                groups = self.assign_unselected_het_molecules_to_nearest_group(
                    hierarchy = filter_h,
                    selections = groups,
                )
            levels.append(groups)
            labels.append('sec. struct.')

        if 'residue' in self.auto_levels:
            logger('Level {}: Creating level with groups for each residue'.format(len(levels)+1))
            groups = [PhenixSelection.format(r) for r in filter_h.residue_groups()]
            levels.append(groups)
            labels.append('residue')

        if 'backbone_sidechain' in self.auto_levels:
            logger('Level {}: Creating level with groups for each residue backbone/sidechain'.format(len(levels)+1))
            b_gps = backbone(filter_h, cbeta=cbeta_flag).atom_groups()
            b_sels = [PhenixSelection.format(r)+back_sel for r in b_gps if (r.resname not in ['ALA','GLY','PRO'])]
            s_gps = sidechains(filter_h, cbeta=(not cbeta_flag)).atom_groups()
            s_sels = [PhenixSelection.format(r)+side_sel for r in s_gps if (r.resname not in ['ALA','GLY','PRO'])]
            groups = sorted(b_sels+s_sels)
            levels.append(groups)
            labels.append('backbone/sidechain')

        if 'atom' in self.auto_levels:
            logger('Level {}: Creating level with groups for each atom'.format(len(levels)+1))
            groups = [PhenixSelection.format(a) for a in filter_h.atoms()]
            levels.append(groups)
            labels.append('atom')
        logger.bar()

        # Insert custom levels
        if self.custom_levels:

            # Print auto levels
            logger('> {} automatic levels created:'.format(len(levels)))
            for i_l, level in enumerate(levels):
                logger('\tLevel {} ({}) - {} groups'.format(i_l+1, labels[i_l], len(level)))
            logger.bar()
            # Insert custom levels
            logger.bar(True, False)
            logger('Processing custom levels:')

            for l_params in self.custom_levels:

                # Extract selections
                label = l_params.label
                groups = l_params.selection

                # Skip blank levels that might be inserted
                if (groups == []):
                    continue

                # Only ONE can be given
                if [l_params.depth, l_params.insert_before, l_params.insert_after].count(None) < 2:
                    msg = ""
                    msg += "For each custom level, you must define EITHER depth OR insert_before OR insert_after (OR none of them). "
                    msg += "\nLevel {} is currently defined with:".format(label)
                    msg += "\n\tdepth:         {}".format(l_params.depth)
                    msg += "\n\tinsert_before: {}".format(l_params.insert_before)
                    msg += "\n\tinsert_after:  {}".format(l_params.insert_after)
                    raise Sorry(msg)

                # label MUST be given (?)
                if label is None:
                    raise Sorry("Must provide label for each custom level")

                # List index to insert level
                if l_params.depth is not None:
                    if l_params.depth < 1:
                        msg = 'Custom level depths cannot be less that 1! (input depth is {} for level with label {})'.format(l_params.depth, label)
                        raise Sorry(msg)
                    if l_params.depth > (len(levels) + 1):
                        msg = 'Trying to add group "{label}" at position {depth}, but only {length} levels currently exist!'.format(
                            label=label,
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
                        msg = 'Trying to insert level "{new}" before level "{ref}" but "{ref}" does not exist yet!'.format(new=label, ref=l_params.insert_before)
                        msg += '\nYou might try adding groups in a different order? (the reference group must be added before this group).'
                        raise Sorry(msg)
                    # insertion index is where the reference group is (will be inserted at this point, shifting reference down)
                    idx = labels.index(l_params.insert_before)
                elif l_params.insert_after is not None:
                    if l_params.insert_after not in labels:
                        msg = 'Trying to insert level "{new}" after level "{ref}" but "{ref}" does not exist yet!'.format(new=label, ref=l_params.insert_after)
                        msg += '\nYou might try adding groups in a different order? (the reference group must be added before this group).'
                        raise Sorry(msg)
                    # insertion index is one after the reference group
                    idx = labels.index(l_params.insert_after) + 1
                else:
                    # Just append to the hierarchy
                    idx = len(levels)

                logger('Inserting level: \n\tLabel: {label}\n\tPosition: {pos}\n\tNumber of groups: {n_groups}'.format(
                    label = label,
                    pos = idx+1,
                    n_groups = len(groups),
                    ))

                assert len(groups) > 0, 'No selections provided for this group!'

                for g in groups:
                    if cache.selection(g).all_eq(False):
                        raise Sorry('Selection "{}" does not select any atoms'.format(g))

                if (self.assign_het_residues_to_nearest_custom_groups is True):
                    groups = self.assign_unselected_het_molecules_to_nearest_group(
                        hierarchy = filter_h,
                        selections = groups,
                    )

                levels.insert(idx, groups)
                labels.insert(idx, label)

            logger.bar()

        # Report
        logger.subheading('Hierarchy summary: {} levels created'.format(len(levels)))
        for i_l, level in enumerate(levels):
            logger.bar()
            logger('Level {} ({}) - {} groups'.format(i_l+1, labels[i_l], len(level)))
            logger.bar()
            for i, l in enumerate(level): logger('\t{:>5d} : {}'.format(i+1,l))

        self.result = group_args(
            level_group_selection_strings = levels,
            level_labels = labels,
            )
        return self.result

    def assign_unselected_het_molecules_to_nearest_group(self,
        hierarchy,
        selections,
        ):

        h = hierarchy
        a = h.atoms()
        cache = h.atom_selection_cache()

        # Select HET atoms from the hierarchy
        het_sel = cache.selection('hetero')
        # Return selections if there are no het atoms
        if sum(het_sel) == 0:
            return selections

        # Extract bool selections
        sel_bool = [cache.selection(s) for s in selections]
        # Extract atoms and coordinates for each selection
        group_atoms = [a.select(s_b) for s_b in sel_bool]
        group_xyz = [a.extract_xyz() for a in group_atoms]

        # Calculate total selection from the groups
        cum_sel = sel_bool[0]
        for s_b in sel_bool[1:]:
            cum_sel = (cum_sel | s_b)

        # Invert the total selection and combine with het_sel
        rest_sel = het_sel & (cum_sel == False)
        if sum(rest_sel) == 0:
            return selections

        # Extract atoms
        rest_h = h.select(rest_sel)

        from giant.structure.formatting import PhenixSelection
        formatter = PhenixSelection()

        import numpy

        # Hash to map groups to the selections to be added
        assigned_groups_hash = {}
        # Iterate through unclaimed groups
        for ag in rest_h.atom_groups():
            # Create selection for atom group and extract coordinates
            ag_sel = formatter.format(ag)
            xyz = ag.atoms().extract_xyz()
            # Find the minimum distance to another group
            dists = [xyz.min_distance_between_any_pair(xyz2) for xyz2 in group_xyz]
            min_dist = numpy.min(dists)
            # If over minimum distance, do not assign
            if min_dist > self.assignment_distance_cutoff:
                continue
            idx = numpy.where(numpy.array(dists)==min_dist)[0][0]
            # Add to hash
            assigned_groups_hash.setdefault(idx, []).append(ag_sel)

        # Start with copy of input
        import copy
        out_selections = copy.deepcopy(selections)

        for i, ag_sels in assigned_groups_hash.items():
            orig_sel = selections[i]
            new_sel = '(' + ') or ('.join([orig_sel]+ag_sels) + ')'
            out_selections[i] = new_sel

        return out_selections


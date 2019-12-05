from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure
from bamboo.common.logs import Log


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



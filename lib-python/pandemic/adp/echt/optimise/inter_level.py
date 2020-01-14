from libtbx import adopt_init_args, group_args
from scitbx.array_family import flex

import numpy

def uij_modulus(uij_array):
    sh = uij_array.shape
    assert sh[-1] == 6
    sh_ = sh[:-1]
    uij_mods = uij_array[...,0:3].mean(axis=len(sh_))
    #uijs_flex = flex.sym_mat3_double(uij_array.reshape((numpy.product(sh_),6)))
    #uijs_eigs = numpy.array(uij_eigenvalues(uijs_flex)).reshape(sh_+(3,))
    #uijs_mods = numpy.mean(uijs_eigs, axis=len(sh_))
    assert uij_mods.shape == sh_
    return uij_mods


class LevelGroupTree:


    def __init__(self,
        links,
        cache_sets = True,
        ):

        # Cache for saving results
        set_cache = {}
        # Max number of possible links
        max_links = len(links.keys())

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
            yield (l_label, l_links.keys())

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

        return map(sorted, output_sets)


class OptimisationSetGenerator:


    def __init__(self, level_group_tree):
        adopt_init_args(self, locals())

    def __call__(self,
        i_datasets,
        max_recursions = None,
        optimise_atomic_adp_amplitudes = True):

        optimisation_sets = self.generate_optimisation_sets(
            max_recursions = max_recursions,
            )

        optimisation_jobs = self.generate_optimisation_jobs_from_sets(
            set_list = optimisation_sets,
            dataset_indices = i_datasets,
            optimise_atomic_adp_amplitudes = optimise_atomic_adp_amplitudes,
            )

        return optimisation_jobs

    def generate_optimisation_sets(self, max_recursions):

        if max_recursions is None:
            return self.level_group_tree.get_complete_sets()
        else:
            return self.level_group_tree.get_recursive_sets(max_recursions)

    def generate_optimisation_jobs_from_sets(self,
        set_list,
        dataset_indices,
        optimise_atomic_adp_amplitudes = True,
        ):
        """Generate control list for optimisation -- if datasets are independent, can be optimised separately, etc.."""

        output = []
        for s in set_list:
            if optimise_atomic_adp_amplitudes is True:
                output.append((dataset_indices, s))                 # One optimisation job
            else:
                output.extend([([i], s) for i in dataset_indices])  # One optimisation job for each dataset
        return output


class OptimiseInterLevelAmplitudes:


    def __init__(self,
        convergence_tolerance,
        optimise_atomic_adp_amplitudes = True,
        optimisation_weights = None,
        verbose = False,
        log = None,
        ):
        if optimisation_weights is not None:
            from mmtbx.tls.optimise_amplitudes import OptimisationWeights
            optimisation_weights = OptimisationWeights.from_other(optimisation_weights)

        adopt_init_args(self, locals())

    def __call__(self,
        uij_target,
        uij_target_weights,
        uij_isotropic_mask,
        model_object,
        level_group_tree,
        adp_optimisation_dataset_mask = None, # TODO
        max_recursions = None,
        ):

        # Extract all groups and amplitudes from the input model object
        group_amplitudes, group_uij_values = self.extract_group_values_for_optimisation(
            model_object = model_object,
            uij_isotropic_mask = uij_isotropic_mask,
            )

        if self.optimise_atomic_adp_amplitudes is True:
            adp_amplitudes, adp_uij_values = self.extract_adp_values_for_optimisation(
                adp_values = model_object.adp_values,
                uij_isotropic_mask = uij_isotropic_mask,
                )
        else:
            # Subtract the adp values from the target as they will not be optimised
            uij_target = uij_target - numpy.array(model_object.adp_values)
            adp_amplitudes, adp_uij_values = None, None

        # Initialise this always, even if not used
        if adp_optimisation_dataset_mask is None:
            adp_optimisation_dataset_mask = flex.bool(model_object.n_datasets, True)
        assert adp_optimisation_dataset_mask.size() == model_object.n_datasets

        get_optimisation_indices_sets = OptimisationSetGenerator(level_group_tree=level_group_tree)

        optimisation_sets = get_optimisation_indices_sets(
            i_datasets = range(model_object.n_datasets),
            max_recursions = max_recursions,
            optimise_atomic_adp_amplitudes = self.optimise_atomic_adp_amplitudes,
            )

        ##############
        # Report
        if self.verbose:
            from itertools import groupby
            log = self.log
            log.bar()
            log('Inter-level optimisations:')
            for i_j, (i_ds, l_g_pairs) in enumerate(optimisation_sets):
                log.bar()
                log('Optimising amplitudes for datasets:\n\t{}'.format(', '.join(map(str,[i+1 for i in i_ds]))))
                log('Optimising amplitudes for groups:')
                for l, g in groupby(sorted(l_g_pairs), lambda x: x[0]):
                    start_str = "Level {}, Groups ".format(l+1)
                    padding = '\n\t'+(' '*len(start_str))
                    group_str = ', '.join([padding*(int(i)%20==19)+str(i_g+1) for i, i_g in enumerate(zip(*g)[1])])
                    log('\t'+start_str+group_str)
            log.bar()
        # Report
        ##############

        print('')
        import tqdm
        pbar = tqdm.tqdm(total=len(optimisation_sets), ncols=100)

        for dataset_indices, level_group_indices in optimisation_sets:

            opt_parameters = self.get_optimisation_variables(
                dataset_indices = dataset_indices,
                level_group_indices = level_group_indices,
                uij_target = uij_target,
                uij_weights = uij_target_weights,
                group_amplitudes = group_amplitudes,
                group_uij_values = group_uij_values,
                group_selections = model_object.tls_selections,
                adp_uij_values = adp_uij_values,
                adp_amplitudes = adp_amplitudes,
                )

            group_amplitudes_sel, adp_amplitudes_sel = self.run_optimisation(
                uij_target = opt_parameters.uij_target,
                uij_weights = opt_parameters.uij_weights,
                base_uij_values = opt_parameters.base_uij_values,
                base_amplitudes = opt_parameters.base_amplitudes,
                base_atom_indices = opt_parameters.base_atom_indices,
                base_dataset_hash = opt_parameters.base_dataset_hash,
                adp_uij_values = opt_parameters.adp_uij_values,
                adp_amplitudes = opt_parameters.adp_amplitudes,
                adp_optimisation_dataset_mask = adp_optimisation_dataset_mask.select(flex.size_t(dataset_indices)),
                )

            self.insert_values(
                dataset_indices = dataset_indices,
                level_group_indices = level_group_indices,
                atom_selection = opt_parameters.atom_selection,
                group_amplitudes_all = group_amplitudes,
                group_amplitudes_new = group_amplitudes_sel,
                adp_amplitudes_all = adp_amplitudes,
                adp_amplitudes_new = adp_amplitudes_sel,
                )

            pbar.update(1)

        pbar.close()
        print('')

        self.apply_group_multipliers(
            model_object = model_object,
            group_multipliers = group_amplitudes,
            )

        if self.optimise_atomic_adp_amplitudes:
            self.apply_adp_multipliers(
                model_object = model_object,
                adp_multipliers = adp_amplitudes,
                adp_uij_values = adp_uij_values,
                )

        return None
        #return (model_object, adp_uij_values)

    @staticmethod
    def extract_group_values_for_optimisation(
        model_object,
        uij_isotropic_mask = None,
        ):

        # One array of uij values
        group_uij_values = numpy.zeros((
            model_object.n_tls_levels,
            model_object.n_modes,
            model_object.n_datasets,
            model_object.n_atoms,
            6,
            ), dtype=float)

        # One array of amplitudes/multipliers
        group_multipliers = [numpy.zeros((
            len(groups),
            model_object.n_modes,
            model_object.n_datasets,
            )) for groups in model_object.tls_objects]

        # Extract amplitudes and uijs group-by-group
        for i_level in xrange(model_object.n_tls_levels):
            for i_group, group in enumerate(model_object.tls_objects[i_level]):
                selection = model_object.tls_selections[i_level][i_group]
                for i_mode, (amps, uijs) in enumerate(group.uijs_unmultiplied()):
                    # Extract the uijs for each dataset
                    for i_dataset in xrange(group.n_datasets):
                        group_uij_values[i_level, i_mode, i_dataset, selection] = uijs[i_dataset]
                    # Extract the amplitudes for this mode/group
                    group_multipliers[i_level][i_group, i_mode, :] = amps

        # Apply isotropic mask if provided
        if uij_isotropic_mask is not None:
            for i_level in xrange(model_object.n_tls_levels):
                for i_mode in xrange(model_object.n_modes):
                    group_uij_values[i_level, i_mode] = uij_isotropic_mask(group_uij_values[i_level, i_mode])

        return (group_multipliers, group_uij_values)

    @staticmethod
    def extract_adp_values_for_optimisation(
        adp_values,
        uij_isotropic_mask = None,
        ):
        # Apply isotropic mask if provided
        if uij_isotropic_mask is not None:
            adp_values = uij_isotropic_mask(adp_values)

        # Convert to numpy array and calculate modulus (size) of each uij
        adp_values = numpy.array(adp_values).reshape((adp_values.all()+(6,)))
        adp_multipliers = uij_modulus(adp_values)

        # Prevent dividing by zero
        mask = (adp_multipliers != 0.0)
        # Calculate normalised adps
        adp_uij_values = numpy.zeros_like(adp_values)
        adp_uij_values[mask] = adp_values[mask] / adp_multipliers.reshape(adp_multipliers.shape+(1,))[mask]

        return (adp_multipliers, adp_uij_values)

    def get_optimisation_variables(self,
        dataset_indices,
        level_group_indices,
        uij_target,
        uij_weights,
        group_amplitudes,
        group_uij_values,
        group_selections,
        adp_uij_values,
        adp_amplitudes,
        ):

        opt_parameters = self.get_fitted_values(
            dataset_indices = dataset_indices,
            level_group_indices = level_group_indices,
            group_amplitudes = group_amplitudes,
            group_uij_values = group_uij_values,
            group_selections = group_selections,
            )

        sel_uij_target, sel_uij_weights = self.get_target_values(
            dataset_indices = dataset_indices,
            level_group_dict_mask = opt_parameters.level_group_dict,
            atom_selection = opt_parameters.atom_selection,
            uij_target = uij_target,
            uij_weights = uij_weights,
            group_amplitudes = group_amplitudes,
            group_uij_values = group_uij_values,
            group_selections = group_selections,
            )

        if adp_uij_values is not None:
            assert adp_amplitudes is not None
            sel_adp_uij_values = flex.sym_mat3_double(adp_uij_values[opt_parameters.atom_selection])
            sel_adp_amplitudes = flex.double(adp_amplitudes[opt_parameters.atom_selection])
        else:
            assert adp_amplitudes is None
            sel_adp_uij_values = flex.sym_mat3_double(opt_parameters.atom_selection.sum())
            sel_adp_amplitudes = flex.double(sel_adp_uij_values.size(), 0.0)

        result = group_args(
            uij_target = sel_uij_target,
            uij_weights = sel_uij_weights,
            atom_selection = opt_parameters.atom_selection,
            base_amplitudes = opt_parameters.base_amplitudes,
            base_uij_values = opt_parameters.base_uij_values,
            base_atom_indices = opt_parameters.base_atom_indices,
            base_dataset_hash = opt_parameters.base_dataset_hash,
            adp_uij_values = sel_adp_uij_values,
            adp_amplitudes = sel_adp_amplitudes,
            )

        return result

    @staticmethod
    def get_fitted_values(
        dataset_indices,
        level_group_indices,
        group_amplitudes,
        group_uij_values,
        group_selections,
        ):

        n_groups = len(level_group_indices)
        n_modes = group_uij_values.shape[1]
        n_datasets = len(dataset_indices)
        n_base = (n_groups * n_modes * n_datasets)

        # Extracted values for selected groups (1d arrays)
        base_amplitudes = flex.double(n_base, 0.0)
        base_uij_values = []
        base_atom_indices = []
        i_base_amplitudes = 0

        # Objects to be returned along with extracted values
        level_group_dict = {}
        base_dataset_hash = flex.size_t(range(n_datasets) * n_groups * n_modes)
        atom_selection = numpy.zeros((group_uij_values.shape[-2],), dtype=bool)

        # Generate necessary reindexing operations and masks
        for (i_l, i_g) in level_group_indices:
            atom_selection[group_selections[i_l][i_g]] = True
            # Add group to dict
            level_group_dict.setdefault(i_l, {}).setdefault(i_g, True)

        # Extract the contributions for optimised selections
        for (i_l, i_g) in level_group_indices:
            # Atom selection for this group
            g_sel = group_selections[i_l][i_g]
            # Reindex this on the subselection
            this_group_atom_indices = flex.bool(g_sel[atom_selection]).iselection()
            # Extract uijs for each mode for each dataset
            for i_m in xrange(group_uij_values.shape[1]):
                for i_d in dataset_indices:
                    # Extract uij values for mode of group in dataset
                    this_group_uij_values = flex.sym_mat3_double(group_uij_values[i_l, i_m, i_d, g_sel])
                    # Append to output lists
                    base_amplitudes[i_base_amplitudes] = group_amplitudes[i_l][i_g,i_m,i_d]
                    base_uij_values.append(this_group_uij_values)
                    base_atom_indices.append(this_group_atom_indices)
                    i_base_amplitudes += 1
        assert i_base_amplitudes == n_base

        return group_args(
            level_group_dict = level_group_dict,
            atom_selection = atom_selection,
            base_amplitudes = base_amplitudes,
            base_uij_values = base_uij_values,
            base_atom_indices = base_atom_indices,
            base_dataset_hash = base_dataset_hash,
            )

    @staticmethod
    def get_target_values(
        dataset_indices,
        level_group_dict_mask,
        atom_selection,
        uij_target,
        uij_weights,
        group_amplitudes,
        group_uij_values,
        group_selections,
        ):

        n_datasets = len(dataset_indices)
        n_atoms = atom_selection.sum()

        # Extract the contributions for non-optimised selections
        l_mults = numpy.zeros_like(group_uij_values)
        for i_l in xrange(len(group_selections)):
            for i_g in xrange(len(group_selections[i_l])):
                # Skip the groups being optimised
                if level_group_dict_mask.get(i_l, {}).get(i_g):
                    continue
                g_sel = group_selections[i_l][i_g]
                # Skip groups outside the atom selection
                if not (g_sel * atom_selection).any():
                    continue
                # Extract multipliers
                for i_t in xrange(group_uij_values.shape[1]):
                    for i_d in dataset_indices: # Only need to iterate over selected datasets
                        l_mults[i_l, i_t, i_d, g_sel] = group_amplitudes[i_l][i_g, i_t, i_d]
        # Multiply the fitted uijs and subtract them from target
        sel_fitted_uij_values = (l_mults * group_uij_values).sum(axis=(0,1))

        # Extract the ones required
        sel_uij_target = (uij_target - sel_fitted_uij_values)
        sel_uij_target = sel_uij_target[dataset_indices]   # select for datasets
        sel_uij_target = sel_uij_target[:, atom_selection]  # select for atoms
        assert sel_uij_target.shape == (n_datasets, n_atoms, 6)
        # Reshape for output
        sel_uij_target = flex.sym_mat3_double(sel_uij_target.reshape((n_datasets*n_atoms,6)))
        sel_uij_target.reshape(flex.grid((n_datasets,n_atoms)))

        # Extract optimisation weights
        sel_uij_weights = uij_weights[dataset_indices]
        sel_uij_weights = uij_weights[:,atom_selection]
        # Renormalise weights
        sel_uij_weights = sel_uij_weights / sel_uij_weights.mean()
        # Reshape for output
        sel_uij_weights = flex.double(sel_uij_weights.reshape((n_datasets*n_atoms,)))
        sel_uij_weights.reshape(flex.grid((n_datasets,n_atoms)))

        return sel_uij_target, sel_uij_weights

    def run_optimisation(self,
        uij_target,
        uij_weights,
        base_uij_values,
        base_amplitudes,
        base_atom_indices,
        base_dataset_hash,
        adp_uij_values,
        adp_amplitudes,
        adp_optimisation_dataset_mask = None,
        ):

        if adp_optimisation_dataset_mask is None:
            adp_optimisation_dataset_mask = flex.bool(uij_target.all()[0], True)

        assert uij_target.all() == uij_weights.all()
        assert base_amplitudes.all() == (len(base_uij_values),)
        assert base_amplitudes.all() == (len(base_atom_indices),)
        assert base_amplitudes.all() == base_dataset_hash.all()
        assert adp_uij_values.all() == (uij_target.all()[1],)
        assert adp_amplitudes.all() == (uij_target.all()[1],)
        assert adp_optimisation_dataset_mask.all() == (uij_target.all()[0],)

        from mmtbx.tls.optimise_amplitudes import OptimiseAmplitudes
        optimiser = OptimiseAmplitudes(
            target_uijs         = uij_target,
            target_weights      = uij_weights,
            base_amplitudes     = base_amplitudes,
            base_uijs           = base_uij_values,
            base_atom_indices   = base_atom_indices,
            base_dataset_hash   = base_dataset_hash,
            atomic_uijs         = adp_uij_values,
            atomic_amplitudes   = adp_amplitudes,
            atomic_optimisation_mask = adp_optimisation_dataset_mask,
            optimisation_weights  = self.optimisation_weights,
            convergence_tolerance = self.convergence_tolerance,
            )
        optimiser.run()

        amplitudes_opt = numpy.array(optimiser.result)
        assert len(amplitudes_opt) == len(base_amplitudes) + len(adp_amplitudes)

        base_amplitudes_new = flex.double(amplitudes_opt[0:len(base_amplitudes)])
        adp_amplitudes_new = flex.double(amplitudes_opt[len(base_amplitudes):])

        assert base_amplitudes_new.all() == base_amplitudes.all()
        assert adp_amplitudes_new.all() == adp_amplitudes.all()

        return base_amplitudes_new, adp_amplitudes_new

    @staticmethod
    def insert_values(
        dataset_indices,
        level_group_indices,
        atom_selection,
        group_amplitudes_all,
        group_amplitudes_new,
        adp_amplitudes_all,
        adp_amplitudes_new,
        ):

        n_modes = group_amplitudes_all[0].shape[1]

        i = 0
        for (i_l, i_g) in level_group_indices:
            for i_m in xrange(n_modes):
                for i_d in dataset_indices:
                    group_amplitudes_all[i_l][i_g,i_m,i_d] = group_amplitudes_new[i]
                    i += 1
        assert i == len(group_amplitudes_new)

        if adp_amplitudes_all is not None:
            adp_amplitudes_all[atom_selection] = adp_amplitudes_new

    @staticmethod
    def apply_group_multipliers(
        model_object,
        group_multipliers,
        ):
        """Apply the optimised multipliers to the input levels"""
        for i_lvl, level_groups in enumerate(model_object.tls_objects):
            for i_g, g in enumerate(level_groups):
                assert i_g == g.index
                mults = group_multipliers[i_lvl][i_g]
                for i_tls in range(model_object.n_modes):
                    mode = g.tls_parameters.get(index=i_tls)
                    mode.amplitudes.set(mults[i_tls])

    @staticmethod
    def apply_adp_multipliers(
        model_object,
        adp_multipliers,
        adp_uij_values,
        ):
        new_adp_values = adp_multipliers.reshape(adp_multipliers.shape+(1,)) * adp_uij_values
        new_adp_values = flex.sym_mat3_double(new_adp_values)
        model_object.adp_values = new_adp_values

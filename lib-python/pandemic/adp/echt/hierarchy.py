import giant.logs as lg
logger = lg.getLogger(__name__)

from libtbx import adopt_init_args, group_args
import numpy
from scitbx.array_family import flex


class CreateEchtModelTask(object):


    adp_level_name = 'atomic'

    def __init__(self,
        n_tls_modes,
        matrix_decimals,
        amplitude_decimals,
        matrix_tolerance,
        amplitude_tolerance,
        ):
        adopt_init_args(self, locals())

        # Construct tasks for later
        from pandemic.adp.echt.hierarchy import CreateMultiDatasetTLSGroupHierarchyTask
        self.tls_objects_constructor = CreateMultiDatasetTLSGroupHierarchyTask(
            n_tls_modes = self.n_tls_modes,
            matrix_decimals = self.matrix_decimals,
            amplitude_decimals = self.amplitude_decimals,
            matrix_tolerance = self.matrix_tolerance,
            amplitude_tolerance = self.amplitude_tolerance,
            )

        self.model_object = None

    def run(self,
            models,
            level_group_array,
            level_group_labels,
            level_labels,
            overall_atom_mask,
            ):

        self.tls_objects_constructor.run(
                models = models,
                level_group_array = level_group_array,
                level_group_labels = level_group_labels,
                level_labels = level_labels,
                overall_atom_mask = overall_atom_mask,
                )

        from pandemic.adp.echt.tls import EchtBFactorModel
        self.model_object = EchtBFactorModel(
                tls_objects     = self.tls_objects_constructor.result.tls_objects,
                tls_selections  = self.tls_objects_constructor.result.tls_selections,
                adp_values      = flex.sym_mat3_double(
                    int(overall_atom_mask.sum()),
                    (0.0,0.0,0.0,0.0,0.0,0.0),
                    ),
                tls_level_names = level_labels,
                adp_level_name  = self.adp_level_name,
                tls_selection_strings = level_group_labels,
                dataset_labels = [m.tag for m in models],
                )

        return self.model_object


class CreateMultiDatasetTLSGroupHierarchyTask(object):


    def __init__(self,
            n_tls_modes,
            matrix_decimals,
            amplitude_decimals,
            matrix_tolerance,
            amplitude_tolerance,
            ):
        adopt_init_args(self, locals())
        # Set precisions and tolerances
        from mmtbx.tls.utils import TLSMatrices, TLSAmplitudes
        TLSMatrices.set_precision(matrix_decimals)
        TLSAmplitudes.set_precision(amplitude_decimals)
        TLSMatrices.set_tolerance(matrix_tolerance)
        TLSAmplitudes.set_tolerance(amplitude_tolerance)

    def run(self,
            models,
            level_group_array,
            level_group_labels = None,
            level_labels = None,
            overall_atom_mask = None,
            ):

        # Create level labels if they don't exist for reporting purposes
        level_labels = (level_labels if level_labels else ['Level {}'.format(i+1) for i in range(len(level_group_array))])

        # Extract coordinates
        coordinates = self.extract_coordinates_from_models(
            models = models,
            atom_mask = overall_atom_mask,
            )

        # Check all compatible
        self.validate_input(
                models = models,
                coordinates = coordinates,
                level_group_array = level_group_array,
                level_group_labels = level_group_labels,
                level_labels = level_labels,
                )

        # Generate origins for the TLS groups that are valid across all levels
        tls_origins_hash, tls_origins = self.get_non_segmenting_tls_origins_and_partitions(
                level_group_array = level_group_array,
                level_labels = level_labels,
                atom_xyzs = coordinates,
                )

        level_group_objects, level_group_selections = self.get_tls_objects(
            models = models,
            n_tls_modes = self.n_tls_modes,
            level_group_array = level_group_array,
            level_group_labels = level_group_labels,
            level_labels = level_labels,
            coordinates = coordinates,
            tls_origins = tls_origins,
            tls_origins_hash = tls_origins_hash,
            )

        self.result = group_args(
                tls_objects = level_group_objects,
                tls_selections = level_group_selections,
                )

        return self.result

    @staticmethod
    def validate_input(
            models,
            coordinates,
            level_group_array,
            level_group_labels,
            level_labels,
            ):
        assert coordinates.ndim == 3
        assert coordinates.shape[0] == len(models)
        assert coordinates.shape[1] == level_group_array.shape[1]
        assert coordinates.shape[2] == 3
        assert len(level_labels) == len(level_group_array)
        assert len(level_labels) == len(level_group_labels)
        for i, l in enumerate(level_group_array):
            assert len(set(l).difference({-1})) == len(level_group_labels[i]), 'number of groups does not match labels'

    @staticmethod
    def get_non_segmenting_tls_origins_and_partitions(
        level_group_array,
        level_labels,
        atom_xyzs,
        ):

        logger.subheading('Identifying partitioning of atoms that does not split any TLS groups')

        array = level_group_array
        labels = level_labels

        atom_hash = None
        atom_coms = []

        # Find the levels with the largest coverage (number of atoms in groups in this level)
        logger('Calculating the number of atoms that are in groups for each level:')
        n_atoms = array.shape[1]
        n_atoms_per_level = [(l>=0).sum() for l in array]
        logger('Atoms covered per level: \n\t{}'.format('\n\t'.join(['Level {} -> {} atoms'.format(i+1, n) for i,n in enumerate(n_atoms_per_level)])))
        levels_with_max_coverage = [i_l for i_l, n in enumerate(n_atoms_per_level) if n==n_atoms]

        if levels_with_max_coverage:
            logger('Only considering partitionings that cover all atoms ({} atoms).'.format(n_atoms))
            logger('Considering partitionings from levels: \n\t{}'.format('\n\t'.join(['{} ({})'.format(l+1, labels[l]) for l in levels_with_max_coverage])))
            # Find how many groups per level, and find level indices with fewest groups
            logger('Looking for which of these levels has the smallest number of groups')
            n_groups_per_level = [len(numpy.unique(l[l>=0])) for l in array[levels_with_max_coverage]]
            logger('Groups per level: \n\t{}'.format('\n\t'.join(['Level {} -> {}'.format(levels_with_max_coverage[i]+1, n) for i,n in enumerate(n_groups_per_level)])))
            levels_with_min_ngroups = [levels_with_max_coverage[i] for i, n in enumerate(n_groups_per_level) if n==min(n_groups_per_level)]
            logger('Only considering partitionings with the smallest number of groups ({} groups).'.format(min(n_groups_per_level)))
            logger('Considering partitionings from levels: \n\t{}'.format('\n\t'.join(['{} ({})'.format(i_l+1, labels[i_l]) for i_l in levels_with_min_ngroups])))
            # Set to False just because (even though it is immediately overwritten)
            found_valid_level = False
            # Test each level to see if it splits tls groups in other levels
            for i_l in levels_with_min_ngroups:
                logger.bar()
                logger('Checking to see if the partitioning on level {} ({}) cuts any TLS groups on other levels'.format(i_l+1, labels[i_l]))
                # Assume level is valid until proven otherwise
                found_valid_level = True
                # Extract groups for this level
                level_groups = array[i_l]
                # Check to see if groups on other levels occupy the non_sel positions (the "null space" of this level)
                level_sel = (level_groups >= 0)
                have_no_group_this_level = array[:,numpy.bitwise_not(level_sel)]
                have_group_other_levels = (have_no_group_this_level >= 0)
                # group members outside selection
                if have_group_other_levels.sum():
                    found_valid_level = False
                    logger('> There are atoms in groups on other levels that do not have a group on this level.')
                    logger('> This partitioning is not appropriate.')
                    continue
                else:
                    logger('> All atoms in other levels are contained by groups on this level.')
                # Get the group idxs (labels/names) for this level
                g_idxs = numpy.unique(level_groups)
                # Iterate through groups and get selection for each
                logger('> Checking that every TLS group on other levels is contained within only one group on this level')
                for g in g_idxs:
                    if g < 0: continue
                    logger('  ...checking this condition for atoms in group {} of this level'.format(g))
                    # Get the selection for this group
                    g_sel = (level_groups==g)
                    # Get the chunk of values for other levels -- for this selection, and NOT this selection
                    g_this = array[:,g_sel]
                    g_other = array[:,numpy.bitwise_not(g_sel)]
                    # Check condition for all other levels
                    for i_l_test in range(len(array)):
                        # Skip checking against itself
                        if i_l==i_l_test: continue
                        logger('     ...checking this condition for the partitioning on level {}.'.format(i_l_test+1))
                        # Check if groups in this selection in this level are present in the level outside this selection
                        overlap = set(g_this[i_l_test]).intersection(set(g_other[i_l_test])).difference({-1})
                        if overlap:
                            logger('      ...there is a group on level {} whose atoms are split between this group and another'.format(i_l_test+1))
                            logger('      ...(groups: {})'.format(', '.join(map(str,sorted(overlap)))))
                            logger('> This partitioning is not appropriate.')
                            found_valid_level = False
                            break
                    # If invalid it's not necessary to check the rest
                    if found_valid_level is False:
                        break
                # All checks passed - return this
                if found_valid_level is True:
                    logger.bar()
                    logger('This level -- level {} ({}) -- has a partitioning that is compatible with all other levels.'.format(i_l+1,labels[i_l]))
                    logger('Using this partitioning.')
                    import copy
                    atom_hash = copy.copy(array[i_l])
                    break
        else:
            logger('No levels that cover all atoms')
            found_valid_level = False

        # No valid levels -- create over-arching partition
        if found_valid_level is False:
            logger('No suitable level found -- using one partition containing all atoms.')
            atom_hash = numpy.ones(array.shape[1], dtype=int)
        # Renumber the array to start from 0
        for i,v in enumerate(sorted(numpy.unique(atom_hash))):
            if v < 0: continue
            atom_hash[atom_hash==v] = i

        # Iterate through and calculate origins for each dataset for each group
        logger.subheading('Calculating Centres-of-mass for each partition')
        for i in sorted(numpy.unique(atom_hash)):
            if i<0: continue
            # Get the selection for this group
            i_sel = (atom_hash == i)
            # Get the atoms in this group
            xyzs = atom_xyzs[:,i_sel,:]
            # Calculate centre of mass for each dataset
            com = numpy.mean(xyzs, axis=1)
            atom_coms.append(com)
        # Convert to numpy array and transpose
        atom_coms = numpy.array(atom_coms)

        # Report
        for i_p, coms in enumerate(atom_coms):
            logger.bar()
            logger('Centres-of-mass for atoms in Partition {}'.format(i_p+1))
            for j, c in enumerate(coms):
                logger('\tDataset {:>3d}: {}'.format(j+1, c))
        logger.bar()

        assert atom_hash.shape == (atom_xyzs.shape[1],)            # n_atoms
        assert atom_coms.shape[0] == len(numpy.unique(atom_hash))  # n_partitions
        assert atom_coms.shape[1] == atom_xyzs.shape[0]            # n_datasets
        assert atom_coms.shape[2] == 3                             # ... n_dim

        return atom_hash, atom_coms

    @staticmethod
    def extract_coordinates_from_models(
            models,
            atom_mask=None,
            ):
        coordinates = numpy.array([m.hierarchy.atoms().extract_xyz() for m in models])
        assert coordinates.ndim == 3
        assert coordinates.shape[0] == len(models)
        assert coordinates.shape[2] == 3
        if atom_mask is not None:
            coordinates = coordinates[:,atom_mask,:]
        return coordinates

    @staticmethod
    def get_tls_objects(
        models,
        n_tls_modes,
        level_group_array,
        level_group_labels,
        level_labels,
        coordinates,
        tls_origins,
        tls_origins_hash,
        ):
        from mmtbx.tls.utils import TLSMatricesAndAmplitudesList
        from pandemic.adp.echt.tls import MultiDatasetTLSGroup
        # Generate tls hierarchy classes
        all_level_group_selections = []
        all_level_group_objects = []
        for i_level, group_num_array in enumerate(level_group_array):
            # Names for the groups
            group_labels = level_group_labels[i_level]
            # Number of the group
            group_nums = sorted(set(group_num_array).difference({-1}))
            this_level_group_selections = []
            this_level_group_objects = []
            for i_group, group_number in enumerate(group_nums):
                if group_number==-1: continue
                # Atomic selection
                group_selection = (group_num_array == group_number)
                # Which origin do we use?
                i_origin = numpy.unique(tls_origins_hash[group_selection])
                assert len(i_origin) == 1
                i_origin = i_origin[0]
                # Create tls parameters object (c++ class)
                tls_parameters = TLSMatricesAndAmplitudesList(
                        length = n_tls_modes,
                        n_amplitudes = len(models),
                        )
                # Initalise to isotropic
                for mode in tls_parameters:
                    mode.amplitudes.set([0.,]*mode.amplitudes.size())
                    mode.matrices.set(values=(1.,1.,1.,0.,0.,0.), component_string='T')
                # Create convenience class for generating uijs from parameters
                group = MultiDatasetTLSGroup(
                    index = group_number,
                    label = group_labels[i_group],
                    tls_parameters = tls_parameters,
                    coordinates = coordinates[:,group_selection],
                    origins = tls_origins[i_origin],
                    )
                # Store the new group object, and the original atoms it refers to
                this_level_group_selections.append(group_selection)
                this_level_group_objects.append(group)
                #level_groups.append((group_selection, group))
            assert [g.index for g in this_level_group_objects] == list(range(len(this_level_group_objects)))
            all_level_group_selections.append(this_level_group_selections)
            all_level_group_objects.append(this_level_group_objects)

        return all_level_group_objects, all_level_group_selections


import giant.logs as lg
logger = lg.getLogger(__name__)

from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure

from pandemic.adp.hierarchy.summary import WriteHierarchicalModelSummaryTask
from pandemic.adp.hierarchy.utils import StructureFactory, MaskedStructureFactory


class CreateHierarchicalModelTask(object):


    def __init__(self,
        auto_levels,
        custom_levels,
        overall_selection,
        cbeta_in_backbone,
        assign_het_residues_to_nearest_ss_groups,
        assign_het_residues_to_nearest_custom_groups,
        remove_duplicate_groups = None,
        n_cpus = 1,
        ):
        adopt_init_args(self, locals())

        # Selection strings for each group for each level
        from pandemic.adp.hierarchy.level_selections import GenerateLevelSelectionsTask
        self.selections_constructor = GenerateLevelSelectionsTask(
            auto_levels = self.auto_levels,
            custom_levels = self.custom_levels,
            overall_selection = self.overall_selection,
            cbeta_in_backbone = self.cbeta_in_backbone,
            assign_het_residues_to_nearest_ss_groups = self.assign_het_residues_to_nearest_ss_groups,
            assign_het_residues_to_nearest_custom_groups = self.assign_het_residues_to_nearest_custom_groups,
            )

        # Construct 2-d integer array of group indices for each atom on each level
        from pandemic.adp.hierarchy.level_array import BuildLevelArrayTask
        self.array_constructor = BuildLevelArrayTask(
            overall_selection = self.overall_selection,
            remove_duplicate_groups = self.remove_duplicate_groups,
            n_cpus = n_cpus,
            )

        from pandemic.adp.hierarchy.level_array_tree import BuildLevelArrayAsTreeTask
        self.array_as_tree = BuildLevelArrayAsTreeTask(
            )

    def run(self,
        hierarchy,
        ):

        logger.heading('Generating hierarchical model')

        selections = self.selections_constructor.run(
            hierarchy = hierarchy,
            )

        arrays = self.array_constructor.run(
            hierarchy = hierarchy,
            level_group_selection_strings = selections.level_group_selection_strings,
            )

        array_tree = self.array_as_tree.run(
            level_group_array = arrays.level_group_array,
            )

        self.result = group_args(
            level_labels                    = selections.level_labels,
            level_group_array               = arrays.level_group_array,
            level_group_selection_strings   = arrays.level_group_selection_strings,
            level_group_tree                = array_tree.tree,
            overall_atom_mask               = arrays.overall_atom_mask,
            # Convenience outputs
            structure_factory               = StructureFactory(master_h=hierarchy),
            masked_structure_factory        = MaskedStructureFactory(master_h=hierarchy, mask=arrays.overall_atom_mask),
            )

        return self.result

import giant.logs as lg
logger = lg.getLogger(__name__)

import math, copy, collections
import numpy, pandas

from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure

from pandemic.adp import constants
from pandemic.adp.utils import show_file_dict

import pathlib as pl


class AssessHierarchyGroupsTask(object):


    show_file_dict = show_file_dict

    def __init__(self,
        output_directory,
        plotting_object,
        min_b_factor = 1.0,
        ):

        adopt_init_args(self, locals())

    def run(self,
        model_object,
        model_hierarchy_info,
        model_hierarchy_files,
        reference_hierarchy,
        overall_atom_selection,
        write_levels_function,
        ):

        if not self.output_directory.exists():
            self.output_directory.mkdir(parents=True)

        # Extract from hierarchy info (generic input object)
        level_labels = model_hierarchy_info.level_labels
        level_group_array = model_hierarchy_info.level_group_array

        output_files = collections.OrderedDict()

        level_group_array_filt = self.remove_null_groups(
            model_uijs = model_object.uijs(),
            level_group_array = level_group_array,
            b_factor_threshold = self.min_b_factor,
            )

        # Get/plot images of the input hierarchy
        input_hierarchy_images = model_hierarchy_files.get('level_partitions_png')
        if (input_hierarchy_images is None):
            input_hierarchy_images = self.plot_level_group_array(
                level_group_array = level_group_array,
                level_labels = level_labels,
                reference_hierarchy = reference_hierarchy,
                overall_atom_selection = overall_atom_selection,
                template_filename = 'input-level-partitions-chain-{}.png',
                )
        output_files['input_partitions_png'] = input_hierarchy_images

        # Plot images of output hierarchy
        output_hierarchy_images = self.plot_level_group_array(
            level_group_array = level_group_array_filt,
            level_labels = level_labels,
            reference_hierarchy = reference_hierarchy,
            overall_atom_selection = overall_atom_selection,
            template_filename = 'output-level-partitions-chain-{}.png',
            )
        output_files['output_partitions_png'] = output_hierarchy_images

        filename =  str(self.output_directory / 'input_hierarchy.eff')
        write_levels_function(
            indices_hierarchy = self.array_to_indices(
                level_group_array = level_group_array,
                ),
            input_level_names = level_labels,
            output_level_names = level_labels,
            output_filename = filename,
            )
        output_files['input_eff_file'] = filename

        filename =  str(self.output_directory / 'output_hierarchy.eff')
        write_levels_function(
            indices_hierarchy = self.array_to_indices(
                level_group_array = level_group_array_filt,
                ),
            input_level_names = level_labels,
            output_level_names = level_labels,
            output_filename = filename,
            )
        output_files['output_eff_file'] = filename

        self.show_file_dict(output_files)

        self.result = group_args(
            output_files = output_files,
            )

        return self.result

    def remove_null_groups(self,
        model_uijs,
        level_group_array,
        b_factor_threshold,
        ):

        logger('Filtering groups with a B-factor less than {}A^2.'.format(b_factor_threshold))

        level_group_array = copy.deepcopy(level_group_array)
        u_threshold = float(b_factor_threshold) / constants.EIGHTPISQ

        # Calculate the maximum iso uij for each atom (potentially across datasets)
        non_atomic_levels_uijs = model_uijs[:-1]
        # mean over diagonal of u to get B, then max over datasets
        uijs_iso_max = non_atomic_levels_uijs[..., 0:3].mean(axis=-1).max(axis=1)
        assert (uijs_iso_max.shape == level_group_array.shape)

        for i_level, level_values in enumerate(level_group_array):

            for i_group in sorted(set(level_values)):
                if i_group == -1:
                    continue
                # Find the atoms in this group
                group_sel = (level_values == i_group)
                # Extract the values for the group
                u_group = uijs_iso_max[i_level, group_sel]

                # Zero group if no atoms above threshold
                if not (u_group > u_threshold).any():
                    level_group_array[i_level, group_sel] = -1

        return level_group_array

    def plot_level_group_array(self,
        level_group_array,
        level_labels,
        reference_hierarchy,
        overall_atom_selection,
        template_filename = 'level-partitions-chain-{}.png',
        ):

        # Create copy of input hierarchy and set all b-factors to -1
        from scitbx.array_family import flex
        reference_hierarchy = reference_hierarchy.deep_copy()
        reference_hierarchy.atoms().set_b(flex.double(reference_hierarchy.atoms_size(), -1))
        # Create structure factory
        from pandemic.adp.hierarchy.utils import StructureFactory
        structure_factory = StructureFactory(master_h=reference_hierarchy)

        # Create hierarchies of the groups for each level
        level_hierarchies = [structure_factory.custom_copy(
            uij = None,
            iso = v,
            mask = overall_atom_selection,
            blank_copy = False) for v in level_group_array]

        output_files = collections.OrderedDict()

        # Write hierarchy plot for each chain
        b_h = structure_factory.blank_copy()
        b_c = b_h.atom_selection_cache()
        for c_id in sorted([c.id for c in b_h.chains()]):
            chain_sel = b_c.selection('chain {}'.format(c_id))
            hierarchies = [h.select(chain_sel, copy_atoms=True) for h in level_hierarchies]
            # Skip if no partitions in this chain
            #if (numpy.array([h.atoms().extract_b() for h in hierarchies]) == -1).all():
            #    continue
            filename = str(self.output_directory / template_filename.format(c_id))
            self.plotting_object.level_plots(
                filename=filename,
                hierarchies=hierarchies,
                labels = ['Level {}\n({})'.format(i_l+1, l) for i_l, l in enumerate(level_labels)],
                title='chain {}'.format(c_id),
                )
            output_files[c_id] = filename

        return output_files

    def array_to_indices(self,
        level_group_array,
        ):

        level_group_indices = []

        for level_values in level_group_array:
            indices = sorted(set(level_values).difference({-1}))
            level_group_indices.append([[i] for i in indices])

        return level_group_indices

    def as_html_summary(self):
        from pandemic.adp.analysis.hierarchy.html import AssessHierarchyGroupsHtmlSummary
        return AssessHierarchyGroupsHtmlSummary(self)

import giant.logs as lg
logger = lg.getLogger(__name__)

import collections

import numpy as np
import pandas as pd
import pathlib as pl

from scitbx.array_family import flex

from giant.utils import merge_dicts

from .pymol import (
    WriteEchtPymolImages,
    WriteEchtPymolScript,
    )

from .graphs import (
    WriteEchtModelGraphs,
    )

from .tables import (
    WriteEchtModelTables,
    )

from .statistics import (
    CalculateBFactorStatistics,
    )

from .structures import (
    WriteEchtAverageStructures,
    )


class AverageEchtModelValues(object):

    def __init__(self, model_object, isotropic_mask):
        """Extract average values"""

        from giant.structure.uij import uij_to_b

        # AVERAGE over *datasets*
        # (n_tls_level, n_mode, *n_dataset*, n_atom, 6)
        self.tls_levels_modes_uij = np.array(
            isotropic_mask( # only need to apply this once
                model_object.uijs_tls().mean(axis=2)
                )
            )
        # (n_tls_level, n_mode, n_atom)
        self.tls_levels_modes_b = np.array(
            [list(map(uij_to_b, u)) for u in self.tls_levels_modes_uij]
            )

        # Now SUM over *modes*
        # (n_tls_level, *n_mode*, n_atom, 6)
        self.tls_levels_uij = self.tls_levels_modes_uij.sum(axis=1)
        # (n_tls_level, *n_mode*, n_atom)
        self.tls_levels_b = self.tls_levels_modes_b.sum(axis=1)

        # (n_atoms, 6)
        self.atomic_level_uij = np.array(
            isotropic_mask(model_object.adp_values)
            )
        # (n_atoms, )
        self.atomic_level_b = np.array(
            uij_to_b(self.atomic_level_uij)
            )

        # (n_levels, n_atom, 6)
        self.all_levels_uij = np.concatenate([
            self.tls_levels_uij,
            self.atomic_level_uij.reshape((1,)+self.atomic_level_uij.shape)
            ])
        # (n_levels, n_atom)
        self.all_levels_b = np.concatenate([
            self.tls_levels_b,
            self.atomic_level_b.reshape((1,)+self.atomic_level_b.shape)
            ])

        # (n_atom, 6)
        self.total_uij = (
            self.all_levels_uij.sum(axis=0)
            )
        # (n_atom, )
        self.total_b = (
            self.all_levels_b.sum(axis=0)
            )

        self.n_levels = model_object.n_levels
        self.n_tls_levels = model_object.n_tls_levels
        self.n_tls_modes = model_object.n_modes
        self.n_atoms = model_object.n_atoms

        self.all_level_names = model_object.all_level_names
        self.tls_level_names = model_object.tls_level_names
        self.adp_level_name = model_object.adp_level_name


class AverageTargetValues(object): 

    def __init__(self, uij_target):

        from giant.structure.uij import uij_to_b
        
        self.total_uij = uij_target.mean(axis=0)
        self.total_b = np.array(uij_to_b(self.total_uij))


class EchtModelSummaryOutput(object):

    def __init__(self,
        output_files,
        level_b_factor_statistics_table,
        ):

        self.output_files = output_files
        self.level_b_factor_statistics_table = level_b_factor_statistics_table


class WriteEchtModelSummary(object):

    def __init__(self,
        output_directory,
        pymol_images = None,
        distribution_images = False,
        ):
        
        self.pymol_images = pymol_images

        self.setup(
            output_directory = output_directory,
            )

    def setup(self, output_directory):

        output_directory = pl.Path(output_directory)

        self.output_directory = output_directory

        self.write_structures = WriteEchtAverageStructures(
            output_directory = output_directory / 'structures',
            )

        self.write_model_graphs = WriteEchtModelGraphs(
            output_directory = output_directory,
            )

        self.write_model_tables = WriteEchtModelTables(
            output_directory = output_directory / 'tables',
            )

        self.write_pymol_images = WriteEchtPymolImages(
            output_directory = output_directory / 'images',
            )

        self.write_pymol_script = WriteEchtPymolScript(
            output_directory = output_directory,
            )

        self.calculate_statistics = CalculateBFactorStatistics()

    def __call__(self,
        overall_atom_mask,
        reference_model,
        model_object,
        isotropic_mask,
        uij_target,
        level_group_array,
        plotting_object,
        ):

        from pandemic.adp.hierarchy.utils import (
            PartitionBordersFactory,
            )

        if not self.output_directory.exists():
            self.output_directory.mkdir(parents=True)

        ###

        # Get average uijs, etc
        model_values = AverageEchtModelValues(
            model_object = model_object,
            isotropic_mask = isotropic_mask,
            )
        
        target_values = AverageTargetValues(
            uij_target = uij_target,
            )

        # Structure-writer
        structure_factory = PartitionBordersFactory(
            master_h = reference_model.hierarchy,
            )

        # Ensure flex
        overall_atom_mask = flex.bool(overall_atom_mask)

        ###

        output_files = collections.OrderedDict()

        #

        of = self.write_structures(
            model_values = model_values,
            target_values = target_values,
            structure_factory = structure_factory,
            atom_mask = overall_atom_mask,
            )

        merge_dicts(master_dict=output_files, merge_dict=of)

        # 

        of = self.write_model_graphs(
            target_values = target_values,
            level_group_array = level_group_array,
            model_values = model_values,
            structure_factory = structure_factory,
            atom_selection = overall_atom_mask,
            plotting_object = plotting_object,
            )

        merge_dicts(master_dict=output_files, merge_dict=of)

        #

        of = self.write_model_tables(
            model_object = model_object,
            )

        merge_dicts(master_dict=output_files, merge_dict=of)

        #

        if self.pymol_images: 

            of = self.write_pymol_images(
                level_pdbs_dict = output_files['level_uijs_pdb'],
                model_object = model_object,
                reference_hierarchy = reference_model.hierarchy,
                overall_atom_mask = overall_atom_mask,
                )

            merge_dicts(master_dict=output_files, merge_dict=of)

        # 

        of = self.write_pymol_script(
            file_dict = output_files,
            )

        merge_dicts(master_dict=output_files, merge_dict=of)    

        # 

        level_b_factor_statistics_table = self.calculate_statistics(
            level_b_values = model_values.all_levels_b,
            level_names = model_values.all_level_names,
            structure_factory = structure_factory,
            overall_atom_mask = overall_atom_mask,
        )

        ###

        self.result = EchtModelSummaryOutput(
            output_files = output_files,
            level_b_factor_statistics_table = level_b_factor_statistics_table,
        )

        return self.result


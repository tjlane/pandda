import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy as np
import pathlib as pl

from giant.structure.formatting import (
    short_labeller,
    )

from giant.structure.common import (
    GetInterestingResnames,
    )

from giant.refinement import (
    get_refiner
    )


class MakeAltlocOmitMap(object):

    def __init__(self,
        refinement_program = 'refmac',
        output_prefix = "omit_altlocs_",
        ):

        self.output_prefix = output_prefix

        self.refine = get_refiner(refinement_program)

    def __call__(self, dataset, altlocs, dir_path=pl.Path()):

        hierarchy = self.get_new_hierarchy(dataset, altlocs)

        out_base = (
            self.output_prefix + ''.join(sorted(altlocs))
            )

        out_prefix = (
            dir_path / out_base
            )

        input_pdb_path = out_prefix.with_suffix('.input.pdb')

        self.write_pdb(
            dataset = dataset,
            hierarchy = hierarchy,
            pdb_path = input_pdb_path,
            )

        result = self.make_omit_map(
            input_pdb_path = input_pdb_path,
            input_mtz_path = dataset.data.filename,
            output_prefix = out_prefix,
            )

        return result

    def get_new_hierarchy(self, dataset, altlocs):

        h = dataset.model.hierarchy

        selection_string = "not ({selection})".format(
            selection = " or ".join([
                "altloc '{altloc}'".format(
                    altloc = a,
                    )
                for a in altlocs
                ]),
            )

        logger(selection_string)

        ac = h.atom_selection_cache()

        selection_bool = ac.selection(selection_string)

        hierarchy = dataset.model.hierarchy.select(selection_bool)

        return hierarchy

    def write_pdb(self, dataset, hierarchy, pdb_path):

        hierarchy.write_pdb_file(
            file_name = str(pdb_path),
            crystal_symmetry = dataset.model.crystal.crystal_symmetry,
            )

    def make_omit_map(self,
        input_pdb_path,
        input_mtz_path,
        output_prefix,
        ):

        refine = self.refine(
            pdb_file = str(input_pdb_path),
            mtz_file = str(input_mtz_path),
            out_prefix = str(output_prefix),
            )

        refine.set_n_cycles(0)

        result = refine.run()

        return (
            refine.output_files['pdb'], refine.output_files['mtz']
            )


class MakeLeaveOneOutOmitMaps(object):

    def __init__(self,
        get_interesting_resnames = None,
        make_altloc_omit_map = None,
        ):

        self.get_interesting_resnames = (
            get_interesting_resnames
            if get_interesting_resnames is not None
            else GetInterestingResnames()
            )

        self.make_altloc_omit_map = (
            make_altloc_omit_map
            if make_altloc_omit_map is not None
            else MakeAltlocOmitMap()
            )

    def __call__(self, dataset, dir_path=pl.Path()):

        resname_altloc_hash = self.get_resname_altloc_sets(dataset)

        unique_altloc_combinations = list(map(
            tuple, sorted(set(resname_altloc_hash.values()))
            ))

        altloc_result_dict = {}

        for altlocs in unique_altloc_combinations:

            logger('Making OMIT map for altlocs {}'.format(altlocs))

            results = self.make_altloc_omit_map(
                dataset = dataset,
                altlocs = altlocs,
                dir_path = dir_path,
                )

            altloc_result_dict[altlocs] = results

        logger(altloc_result_dict)

        residue_result_dict = {
            rg_label : altloc_result_dict[altlocs]
            for rg_label, altlocs
            in resname_altloc_hash.items()
        }

        logger(residue_result_dict)

        return residue_result_dict

    def get_focus_residue_groups(self, hierarchy):

        resnames = self.get_interesting_resnames(hierarchy)

        logger.debug('Selected resnames: {s}'.format(s=str(resnames)))

        if len(resnames) == 0:

            logger('No interesting residue names identified')

            return None

        resnames_set = set(resnames)

        residue_groups = [
            rg
            for rg in hierarchy.residue_groups()
            if resnames_set.intersection(rg.unique_resnames())
            ]

        return residue_groups

    def get_resname_altloc_sets(self, dataset):

        residue_groups = self.get_focus_residue_groups(dataset.model.hierarchy)

        # Get conformers for each residue

        conformer_sets = {
            short_labeller(rg) : tuple(sorted(c.altloc for c in rg.conformers()))
            for rg in residue_groups
        }

        # Get the conformers not used by any residue

        all_conformers = set([
            a for a in dataset.model.hierarchy.altloc_indices() if a.strip(' ')
            ])

        all_conformers_used = set(
            np.concatenate(list(conformer_sets.values())).tolist()
            )

        unused_conformers = (
            all_conformers.difference(all_conformers_used)
            )

        if len(unused_conformers) > 0:

            conformer_sets[None] = tuple(sorted(unused_conformers))

        return conformer_sets


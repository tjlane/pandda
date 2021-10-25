import giant.logs as lg
logger = lg.getLogger(__name__)

from giant.structure.formatting import (
    PhenixSelection,
    )

#
# Phenix
#

class PhenixFormatter(object):

    selection = PhenixSelection


class PhenixDistanceRestraintFormatter(PhenixFormatter):

    _string = (
        "bond {{\n"
        "    action = *add\n"
        "    atom_selection_1 = {string1}\n"
        "    atom_selection_2 = {string2}\n"
        "    distance_ideal = {distance}\n"
        "    sigma = {sigma}\n"
        "    slack = None\n"
        "}}"
        )

    def __init__(self):

        pass

    def __call__(self, distance_restraint, add=True):

        r = distance_restraint

        return self._string.format(
            string1 = self.selection.format(r.atom1),
            string2 = self.selection.format(r.atom2),
            distance = r.length,
            sigma = r.sigma,
            )


class PhenixOccupancyRestraintFormatter(PhenixFormatter):

    _group_string_1 = 'selection = {selection}'
    _group_string_2 = 'constrained_group {{\n    {selection_block}\n}}'

    def __init__(self):

        pass

    def __call__(self, occupancy_restraint):

        group_strings = []

        for group in occupancy_restraint.occupancy_groups:

            formatted_string = self.format_group(
                objects = group.objects,
                )

            group_strings.append(formatted_string)

        return self.format_block(
            group_strings = group_strings,
            )

    def format_group(self, objects):

        group_string = self._group_string_1.format(
            selection = self.selection.join_or(
                [
                    self.selection.format(obj)
                    for obj in objects
                    ],
                extra_join = '\\\n    ',
                ),
            )

        return group_string

    def format_block(self, group_strings):

        return self._group_string_2.format(
            selection_block = (
                '\n'.join(group_strings).replace('\n','\n    ')
                ),
            )


class PhenixRestraintCollectionFormatter(object):

    _distance_string = "refinement.geometry_restraints.edits {{\n    {restraints_string}\n}}"
    _occupancy_string = "refinement.refine.occupancies {{\n    {restraints_string}\n}}"

    def __init__(self):

        self.format_distance_restraint = (
            PhenixDistanceRestraintFormatter()
            )

        self.format_occupancy_restraint = (
            PhenixOccupancyRestraintFormatter()
            )

    def __call__(self,
        restraints_collection,
        filepath = None,
        ):

        rc = restraints_collection
        fp = filepath

        blocks = []

        #

        distance_restraints_block = self.format_distance_restraints(
            distance_restraints = rc.distance_restraints,
            )

        if distance_restraints_block is not None:
            blocks.append(distance_restraints_block)

        #

        occupancy_restraints_block = self.format_occupancy_restraints(
            occupancy_restraints = rc.occupancy_restraints,
            )

        if occupancy_restraints_block is not None:
            blocks.append(occupancy_restraints_block)

        #

        output_block = self.format_blocks(
            blocks = blocks,
            )

        if (fp is not None) and (output_block is not None):

            with open(str(fp), 'a') as fh:
                fh.write(
                    output_block
                    )

        return output_block

    def format_distance_restraints(self, distance_restraints):

        if (distance_restraints is None) or (len(distance_restraints) == 0):
            return None

        restraints_strings = [
            self.format_distance_restraint(r)
            for r in distance_restraints
            ]

        return self._distance_string.format(
            restraints_string = (
                '\n'.join(restraints_strings).replace('\n','\n    ')
                ),
            )

    def format_occupancy_restraints(self, occupancy_restraints):

        if (occupancy_restraints is None) or (len(occupancy_restraints) == 0):
            return None

        restraints_strings = [
            self.format_occupancy_restraint(r)
            for r in occupancy_restraints
            ]

        return self._occupancy_string.format(
            restraints_string = (
                '\n'.join(restraints_strings).replace('\n','\n    ')
                ),
            )

    def format_blocks(self, blocks):

        if (blocks is None) or (len(blocks) == 0):
            return None

        return '\n'.join(blocks)

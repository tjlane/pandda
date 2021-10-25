import giant.logs as lg
logger = lg.getLogger(__name__)

from giant.structure.formatting import (
    RefmacSelection,
    )

#
# Refmac
#


class RefmacFormatter(object):

    selection = RefmacSelection


class RefmacDistanceRestraintFormatter(RefmacFormatter):

    _string = (
        "exte dist "
        "first {string1} "
        "second {string2} "
        "value {distance} "
        "sigma {sigma} "
        "type {type}"
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
            type = (
                1 if add else 0
                ),
            )


class RefmacOccupancyRestraintFormatter(RefmacFormatter):

    _group_string_1 = 'occupancy group id {idx} {selection}'
    _group_string_2 = 'occupancy group alts {complete} {group_idxs}'

    def __init__(self):

        self.idx = 0

    def __call__(self, occupancy_restraint):

        group_idxs = []
        group_strings = []

        for group in occupancy_restraint.occupancy_groups:

            group_idx, formatted_strings = self.format_group(
                objects = group.objects,
                )

            group_idxs.append(group_idx)
            group_strings.extend(formatted_strings)

        return self.format_block(
            group_idxs = group_idxs,
            group_strings = group_strings,
            complete = occupancy_restraint.complete,
            )

    def format_group(self, objects):

        self.idx += 1

        group_strings = [
            self._group_string_1.format(
                idx = str(self.idx),
                selection = str(self.selection.format(obj)),
                )
            for obj in objects
            ]

        return (self.idx, group_strings)

    def format_block(self, group_idxs, group_strings, complete):

        return '\n'.join([
            '\n'.join(group_strings),
            self._group_string_2.format(
                complete = (
                    'complete' if (complete is True) else 'incomplete'
                    ), 
                group_idxs = ' '.join(map(str,group_idxs)),
                ),
            ])


class RefmacRestraintCollectionFormatter(object):

    _occupancy_string = "{restraints_string} \noccupancy refine "

    def __init__(self):

        self.format_distance_restraint = (
            RefmacDistanceRestraintFormatter()
            )

        self.format_occupancy_restraint = (
            RefmacOccupancyRestraintFormatter()
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

        return '\n'.join([
            self.format_distance_restraint(r)
            for r in distance_restraints
            ])

    def format_occupancy_restraints(self, occupancy_restraints):

        if (occupancy_restraints is None) or (len(occupancy_restraints) == 0):
            return None

        restraints_strings = [
            self.format_occupancy_restraint(r)
            for r in occupancy_restraints
            ]

        return self._occupancy_string.format(
            restraints_string = '\n'.join(restraints_strings),
            )

    def format_blocks(self, blocks):

        if (blocks is None) or (len(blocks) == 0):
            return None

        return '\n'.join(blocks)

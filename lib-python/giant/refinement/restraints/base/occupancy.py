import copy

from giant.structure.formatting import (
    Labeller,
    GenericSelection,
    )

#
# Occupancy Restraints
#


class OccupancyGroup(object):

    def __init__(self,
        objects,
        label = None,
        ):

        self.objects = objects
        self.label = label

    def __str__(self):

        formatted_string = '\n'.join([
            str(Labeller.format(o)).strip()
            for o in self.objects
            ])

        s_ = (
            'Occupancy Group:\n'
            '| Label: {name}\n'
            '| Elements:\n'
            '|\t{object_string}\n'
            '`---->'
            ).format(
                name = str(self.label),
                object_string = str(
                    formatted_string
                    ).strip('\n').replace('\n','\n|\t'),
                )

        return s_.strip()


class OccupancyRestraint(object):

    def __init__(self,
        occupancy_groups,
        complete = True,
        ):

        self.occupancy_groups = occupancy_groups
        self.complete = complete

    def __str__(self):

        formatted_groups = '\n'.join([
            str(g).strip()
            for g in self.occupancy_groups
            ])

        s_ = (
            'Occupancy Restraint:\n'
            '| Complete: {complete}\n'
            '| Groups:\n'
            '|\t{groups}\n'
            '`---->'
            ).format(
                complete = str(self.complete),
                groups = str(
                    formatted_groups
                    ).strip('\n').replace('\n','\n|\t'),
                )

        return s_.strip()

    def add_group(self, objects, label=None):

        self.occupancy_groups.append(
            OccupancyGroup(
                objects = objects,
                label = label,
                )
            )

    @classmethod
    def from_list_of_lists(cls,
        objects_list,
        *args, **kwargs
        ):

        return cls(
            occupancy_groups = [
                OccupancyGroup(
                    objects = objects,
                    )
                for objects in objects_list
                ],
            *args, **kwargs
            )

    def set_complete(self, complete=True):
        self.complete = bool(complete)


class OccupancyRestraintList(object):

    name = "OccupancyRestraintList"

    def __init__(self,
        occupancy_restraints = None,
        ):

        self.occupancy_restraints = []

        if occupancy_restraints is not None:
            self.add(
                occupancy_restraints = occupancy_restraints,
                )

    def __str__(self):

        formatted_restraints = '\n'.join([
            str(r).strip('\n')
            for r in self.occupancy_restraints
            ])

        s_ = (
            'Object Type: {name}\n'
            '| Contents:\n'
            '|\t{restraints}\n'
            '`---->'
            ).format(
            name = str(self.name),
            restraints = str(
                formatted_restraints
                ).strip('\n').replace('\n','\n|\t'),
            )

        return s_.strip()

    def __iter__(self):

        return iter(self.occupancy_restraints)

    def __len__(self):

        return len(self.occupancy_restraints)

    def add(self, occupancy_restraints):

        if hasattr(occupancy_restraints, 'occupancy_restraints'):
            occupancy_restraints = copy.deepcopy(
                occupancy_restraints.occupancy_restraints
                )
            #occupancy_restraints = occupancy_restraints.occupancy_restraints

        self.occupancy_restraints.extend(
            occupancy_restraints
            )

        return self

    def remove(self, restraint):

        self.occupancy_restraints.remove(restraint)

    def set_complete(self, complete=True):

        for r in self.occupancy_restraints:
            r.set_complete(complete)


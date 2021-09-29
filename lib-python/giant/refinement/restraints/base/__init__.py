import copy

from .distance import (
    DistanceRestraint,
    DistanceRestraintList,
    )

from .occupancy import (
    OccupancyGroup,
    OccupancyRestraint,
    OccupancyRestraintList,
    )

#
# Collection Objects
#


class RestraintsCollection(object):

    """Collections of distance restraints, occupancy groups, etc."""

    name = "RestraintsCollection"

    def __init__(self,
        distance_restraints = None,
        occupancy_restraints = None,
        label = None,
        ):

        self.label = (
            str(label)
            if label is not None
            else None
            )

        self.distance_restraints = DistanceRestraintList(
            distance_restraints = distance_restraints,
            )

        self.occupancy_restraints = OccupancyRestraintList(
            occupancy_restraints = occupancy_restraints,
            )

    def __str__(self):

        s_ = (
            'Object Type: {name}\n'
            '| Label: {label}\n'
            '| Distance Restraints:\n'
            '|\t{distance_restraints}\n'
            '| Occupancy Restraints:\n'
            '|\t{occupancy_restraints}\n'
            '`---->'
            ).format(
                name = str(self.name),
                label = str(self.label),
                distance_restraints = str(
                    self.distance_restraints
                    ).strip('\n').replace('\n','\n|\t'),
                occupancy_restraints = str(
                    self.occupancy_restraints
                    ).strip('\n').replace('\n','\n|\t'),
            )

        return s_.strip()

    def add(self, other):

        self.distance_restraints.add(
            distance_restraints = other.distance_restraints,
            )

        self.occupancy_restraints.add(
            occupancy_restraints = other.occupancy_restraints,
            )

    def format(self, program='refmac', filepath=None):

        formatter = (
            RefmacFormatter()
            if (program=='refmac') else
            PhenixFormatter()
            if (program=='phenix') else
            NotImplemented()
            )

        output = formatter(
            restraints_collection = self,
            filepath = str(filepath),
            )

        return output


#
# Dummy
#


class DummyRestraintMaker(object):

    name = "DummyRestraintMaker"

    def __init__(self, name=None, *args, **kwargs):

        if name is not None:
            self.name = str(name)

    def __str__(self):

        s_ = (
            'Task: {name}\n'
            '| This task does nothing.\n'
            '`---->'
            ).format(
                name = self.name,
                )

        return s_.strip()

    def __call__(self, hierarchy):
        return RestraintsCollection()

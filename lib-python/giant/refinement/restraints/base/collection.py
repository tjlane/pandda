from .distance import (
    DistanceRestraintList,
    )

from .occupancy import (
    OccupancyRestraintList,
    )


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

        from ..format import (
            get_restraint_formatter,
            )

        formatter = get_restraint_formatter(
            program = program,
            )

        output = formatter(
            restraints_collection = self,
            filepath = str(filepath),
            )

        return output

    def prune_conflicting(self, method='simple'):

        # can't do imports at top because of circular referencing

        if method == 'simple':

            from .prune import (
                PruneConflictingRestraintsSimple,
                )

            prune = PruneConflictingRestraintsSimple(
                in_place = True,
                )

        prune(self)

        return self


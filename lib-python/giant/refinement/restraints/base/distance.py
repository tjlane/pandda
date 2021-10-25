import copy

from giant.structure.formatting import (
    Labeller,
    GenericSelection,
    )

#
# Distance Restraints
#


class DistanceRestraint(object):

    def __init__(self,
        atom1,
        atom2,
        length,
        sigma,
        ):

        self.atom1 = GenericSelection.to_dict(atom1)
        self.atom2 = GenericSelection.to_dict(atom2)

        self.length = length
        self.sigma = sigma

    def __str__(self):

        s_ = (
            'Distance Restraint:\n'
            '| {atom1}\n'
            '| {atom2}\n'
            '| length {length} (sigma {sigma})\n'
            '`---->'
            ).format(
                atom1 = Labeller.format(self.atom1),
                atom2 = Labeller.format(self.atom2),
                length = self.length,
                sigma = self.sigma,
            )

        return s_.strip()

    def set_sigma(self, sigma):
        self.sigma = float(sigma)


class DistanceRestraintList(object):

    name = "DistanceRestraintList"

    def __init__(self,
        distance_restraints = None,
        ):

        self.distance_restraints = []

        if distance_restraints is not None:
            self.add(
                distance_restraints = distance_restraints,
                )

    def __str__(self):

        formatted_restraints = '\n'.join([
            str(r).strip()
            for r in self.distance_restraints
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

        return iter(self.distance_restraints)

    def __len__(self):

        return len(self.distance_restraints)

    def add(self, distance_restraints):

        if hasattr(distance_restraints, 'distance_restraints'):
            distance_restraints = copy.deepcopy(
                distance_restraints.distance_restraints
                )

        self.distance_restraints.extend(
            distance_restraints
            )

        return self

    def remove(self, restraint):

        self.distance_restraints.remove(restraint)

    def set_sigma(self,
        sigma,
        ):

        for r in self.distance_restraints:
            r.set_sigma(sigma)



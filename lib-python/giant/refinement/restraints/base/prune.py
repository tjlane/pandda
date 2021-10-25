import giant.logs as lg
logger = lg.getLogger(__name__)

from . import (
    RestraintsCollection
    )

from giant.structure.formatting import (
    short_labeller,
    )


class PruneConflictingRestraintsSimple(object):

    def __init__(self, in_place=True):

        self.in_place = bool(in_place)

    def __call__(self, restraints_collection):

        logger.debug('** Pruning redundant restraints **')

        pruned = RestraintsCollection()

        if self.in_place is False:

            restraints_collection = copy.deepcopy(restraints_collection)

        pruned_distance_restraints = self.prune_distance_restraints(
            restraints_collection = restraints_collection,
            )

        if pruned_distance_restraints is not None:
            pruned.add(pruned_distance_restraints)

        pruned_occupancy_restraints = self.prune_occupancy_restraints(
            restraints_collection = restraints_collection,
            )

        if pruned_occupancy_restraints is not None:
            pruned.add(pruned_occupancy_restraints)

        logger.debug("** Pruned restraints **\n"+str(pruned))

        return restraints_collection

    def prune_distance_restraints(self, restraints_collection):

        prune_list = []

        r_hash = {}

        for r in restraints_collection.distance_restraints:

            label = (
                short_labeller(r.atom1),
                short_labeller(r.atom2),
                )

            # Get any previous group with this label
            previous = r_hash.get(label, None)

            if previous is not None:
                prune_list.append(previous)

            # Replace with the current
            r_hash[label] = r

        for r in prune_list:

            restraints_collection.distance_restraints.remove(r)

        return RestraintsCollection(
            distance_restraints = prune_list,
            )

    def prune_occupancy_restraints(self, restraints_collection):

        prune_list = []

        r_hash = set()

        for r in restraints_collection.occupancy_restraints:

            prune = False

            r_labels = set()

            for g in r.occupancy_groups:

                for obj in g.objects:

                    label = short_labeller(obj)

                    if label in r_hash:
                        prune = True
                        break

                    r_labels.add(label)

                if prune is True:
                    break

            if prune is True:
                prune_list.append(r)
            else:
                r_hash.update(r_labels)

        for r in prune_list:

            restraints_collection.occupancy_restraints.remove(r)

        return RestraintsCollection(
            occupancy_restraints = prune_list,
            )

import giant.logs as lg
logger = lg.getLogger(__name__)

from .base import (
    RestraintsCollection,
    )

from .occupancy import (
    MakeSimpleOccupancyRestraints,
    MakeMultiStateOccupancyRestraints,
    )

from .conformers import (
    MakeIntraConformerRestraints,
    MakeDuplicateConformerRestraints,
    )


class MakeMultiStateRestraints(object):

    def __init__(self,
        #prune_conflicting_restraints = True,
        make_intra_conformer_restraints = None,
        make_duplicate_conformer_restraints = None,
        make_simple_occupancy_restraints = None,
        make_multi_state_occupancy_restraints = None,
        ):

        self.make_intra_conformer_restraints = (
            make_intra_conformer_restraints
            if (make_intra_conformer_restraints is not None)
            else MakeIntraConformerRestraints()
            )

        self.make_duplicate_conformer_restraints = (
            make_duplicate_conformer_restraints
            if (make_duplicate_conformer_restraints is not None)
            else MakeDuplicateConformerRestraints()
            )

        self.make_simple_occupancy_restraints = (
            make_simple_occupancy_restraints
            if (make_simple_occupancy_restraints is not None)
            else MakeSimpleOccupancyRestraints()
            )

        self.make_multi_state_occupancy_restraints = (
            make_multi_state_occupancy_restraints
            if (make_multi_state_occupancy_restraints is not None)
            else MakeMultiStateOccupancyRestraints()
            )

    def __call__(self, hierarchy):

        rc = RestraintsCollection()

        logger.heading('Making multi-state restraints')

        for restraint_maker in [
            self.make_intra_conformer_restraints,
            self.make_duplicate_conformer_restraints,
            self.make_simple_occupancy_restraints,
            self.make_multi_state_occupancy_restraints,
            ]:

            if restraint_maker is None:
                continue

            logger.subheading(
                'Making restraints for {}'.format(
                    restraint_maker.name
                    )
                )
            logger(str(restraint_maker))

            restraints = restraint_maker(
                hierarchy = hierarchy,
                )

            logger('\nOutput Restraints:\n')
            logger(str(restraints))

            if restraints is not None:
                rc.add(restraints)

        # TODO!
        #rc.prune_conflicting()

        return rc



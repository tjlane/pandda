import giant.logs as lg
logger = lg.getLogger(__name__)

from iotbx.pdb.hierarchy import (
    common_residue_names_get_class,
    )

common_molecules = [
        'ACE',
        'EDO',
        'DMS',
        'GOL',
        'PO4',
        'SO4',
    ]

class GetInterestingResnames(object):

    name = "GetInterestingResnames"

    common_molecules = list(common_molecules)

    def __init__(self,
        ignore_common_molecules = False,
        include_resnames_list = None,
        ignore_resnames_list = None,
        ):

        self.ignore_common_molecules = (
            bool(ignore_common_molecules)
            )

        self.include_resnames_list = (
            list(include_resnames_list)
            if (include_resnames_list is not None)
            else None
            )

        self.ignore_resnames_list = (
            list(ignore_resnames_list)
            if (ignore_resnames_list is not None)
            else None
            )

    def __str__(self):

        s_ = (
            'Task: {name}\n'
            '| ignore_common_molecules: {ignore_common_molecules}\n'
            '| include_resnames_list: {include_resnames_list}\n'
            '| ignore_resnames_list: {ignore_resnames_list}\n'
            '`---->'
            ).format(
                name = self.name,
                ignore_common_molecules = str(self.ignore_common_molecules),
                include_resnames_list = str(self.include_resnames_list),
                ignore_resnames_list = str(self.ignore_resnames_list),
                )

        return s_.strip()

    def __call__(self, hierarchy):

        unq_resnames = list(hierarchy.overall_counts().resnames.keys())

        logger.debug(
            'Residues in structure (and classes): \n\t{}'.format(
                '\n\t'.join(
                    map(str,sorted(
                        [
                            (n,common_residue_names_get_class(n))
                            for n in unq_resnames
                            ],
                        key = lambda t: t[1],
                        ))
                    )
                )
            )

        int_resnames = set([
            r for r in unq_resnames
            if 'common_' not in common_residue_names_get_class(r)
            ])

        logger.debug(
            "Interesting residue names -> {}".format(
                str(int_resnames),
                )
            )

        # Add back in small molecules
        if self.ignore_common_molecules is False:
            #
            int_resnames.update([
                r for r in unq_resnames
                if common_residue_names_get_class(r) not in [
                    'common_amino_acid',
                    'common_element',
                    'common_water',
                    ]
                ])
            #
            logger.debug(
                "After adding small/common molecules -> {}".format(
                    str(int_resnames),
                    )
                )
        else:
            # Make sure other (non-phenix-defined) common molecules are removed
            int_resnames.difference_update(self.common_molecules)
            #
            logger.debug(
                "After removing small molecules -> {}".format(
                    str(int_resnames),
                    )
                )

        # Override -- must be last
        if self.include_resnames_list is not None:
            #
            int_resnames.update(self.include_resnames_list)
            #
            logger.debug(
                "After adding selected molecules ({}) -> {}".format(
                    str(self.include_resnames_list),
                    str(int_resnames),
                    )
                )

        # Override -- must be last
        if self.ignore_resnames_list is not None:
            int_resnames.difference_update(self.ignore_resnames_list)
            logger.debug(
                "After removing selected molecules ({}) -> {}".format(
                    str(self.ignore_resnames_list),
                    str(int_resnames),
                    )
                )

        return sorted(int_resnames)


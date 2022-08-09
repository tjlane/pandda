import giant.logs as lg
logger = lg.getLogger(__name__)

from .cif_tools import CifMerger
from .links import UpdateLinksFromCif


class ConsolidationResult:

    def __init__(self,
        model,
        cif_manager,
        link_records,
        ):

        self.model = model
        self.cif_manager = cif_manager
        self.link_records = link_records

    def write_pdb(self, path):

        logger('Writing {}'.format(str(path)))

        with open(str(path), 'w') as fh:
            fh.write("HEADER    ----                                                XXXX\n")
            fh.write("TITLE     ---\n")
            fh.write("COMPND    ---\n")

        self.model.hierarchy.write_pdb_file(
            str(path),
            open_append = True,
            crystal_symmetry = self.model.crystal.crystal_symmetry,
            link_records = '\n'.join(map(str, self.link_records)),
            )

    def write_cif(self, path):

        logger('Writing {}'.format(str(path)))

        self.cif_manager.write_cif(
            str(path),
            )


class ConsolidateRestraintsAndUpdateModel:

    def __init__(self):

        self.merge_cifs = CifMerger()

        self.update_link_records = UpdateLinksFromCif()

    def __call__(self,
        model,
        cif_managers,
        ):

        merged_cif_obj = (
            self.merge_cifs(
                cif_managers = cif_managers,
                )
            if len(cif_managers) > 1
            else cif_managers[0]
            )

        logger.subheading(
            'Consolidated Cif Object'
            )

        logger(
            str(merged_cif_obj)
            )

        if model is None:

            link_records = None

        else:

            logger.subheading(
                'Updating link records'
                )

            link_records = self.update_link_records(
                model = model,
                cif_manager = merged_cif_obj,
                )

            logger.subheading(
                'Updated link records from model'
                )

            logger(
                '\n'.join(map(str,link_records))
                )

        return ConsolidationResult(
            model = model,
            cif_manager = merged_cif_obj,
            link_records = link_records,
            )


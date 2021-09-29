import giant.logs as lg
logger = lg.getLogger(__name__)

from .phenix import (
    PhenixRestraintCollectionFormatter,
    )

from .refmac import (
    RefmacRestraintCollectionFormatter,
    )

#
# General / collators
#

def get_restraint_formatter(program):

    if str(program) == 'refmac':
        return RefmacRestraintCollectionFormatter()

    if str(program) == 'phenix':
        return PhenixRestraintCollectionFormatter()

    raise NotImplementedError()


class WriteRestraints(object):

    def __init__(self,
        formats,
        output_path_dict = None,
        ):

        self.formats = formats

        self.output_path_dict = (
            output_path_dict
            if output_path_dict is not None
            else dict()
            )

        self.formatters = {
            k : get_restraint_formatter(k)
            for k in formats
        }

    def __call__(self, restraint_collection):

        formatted = {}

        logger.heading('Writing Restraints')

        for k, formatter in sorted(self.formatters.items()):

            logger.subheading('Writing restraints for {}'.format(k))

            formatted[k] = formatter(
                restraint_collection = restraint_collection,
                filepath = self.output_path_dict.get(k),
                )

            logger(str(formatted[k]))

        return formatted

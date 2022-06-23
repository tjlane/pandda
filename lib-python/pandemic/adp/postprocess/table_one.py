import giant.logs as lg
logger = lg.getLogger(__name__)

import os, copy
from libtbx import adopt_init_args
from libtbx.utils import Sorry, Failure

from giant.jiffies import multi_table_ones


class CalculateTableOnes(object):

    debug = False

    def __init__(self,
        output_directory,
        table_one_options,
        n_cpus = 1,
        ):

        # Populate table one phil
        base_phil = multi_table_ones.master_phil.extract()
        base_phil.input.dir        = []
        base_phil.input.labelling  = 'foldername'
        base_phil.options          = table_one_options
        base_phil.settings.cpus    = n_cpus
        base_phil.settings.verbose = self.debug

        adopt_init_args(self, locals())

    def __call__(self,
        output_prefix,
        structures,
        ):

        output_prefix = os.path.join(self.output_directory, output_prefix)
        output_csv = output_prefix + '.csv'
        output_eff = output_prefix + '.eff'

        if os.path.exists(output_csv):
            raise Failure('Output table ones already exist.')

        phil = copy.deepcopy(self.base_phil)
        phil.input.pdb = structures
        phil.output.parameter_file = output_eff
        phil.output.output_basename = output_prefix

        logger.subheading("Making Table One parameter file")

        # Generate parameter file!
        multi_table_ones.run(params=phil)

        from giant.dispatcher import Dispatcher
        prog = Dispatcher('phenix.table_one')
        prog.append_arg(output_eff)

        logger.subheading("Making Table One")
        logger(prog.as_string()+'\n')

        prog.run()
        prog.write_output(output_eff.replace('.eff','.log'))

        if not os.path.exists(output_csv):
            raise Failure('Failed to make table one: {}'.format(output_csv))

        return output_csv


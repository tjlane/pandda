import os, copy
from libtbx import adopt_init_args
from libtbx.utils import Sorry, Failure
from bamboo.common.logs import Log

from giant.jiffies import multi_table_ones


class CalculateTableOnes:


    def __init__(self,
        output_directory,
        table_one_options,
        n_cpus = 1,
        verbose = False,
        log = None,
        ):

        if log is None: log = Log()

        # Populate table one phil
        base_phil = multi_table_ones.master_phil.extract()
        base_phil.input.dir        = []
        base_phil.input.labelling  = 'foldername'
        base_phil.options          = table_one_options
        base_phil.settings.cpus    = n_cpus
        base_phil.settings.verbose = verbose

        adopt_init_args(self, locals())

    def __call__(self,
        output_prefix,
        structures,
        ):

        from bamboo.common.command import CommandManager

        output_prefix = os.path.join(self.output_directory, output_prefix)
        output_csv = output_prefix + '.csv'
        output_eff = output_prefix + '.eff'

        if os.path.exists(output_csv):
            raise Failure('Output table ones already exist.')

        phil = copy.deepcopy(self.base_phil)
        phil.input.pdb = structures
        phil.output.parameter_file = output_eff
        phil.output.output_basename = output_prefix

        # Generate parameter file!
        multi_table_ones.run(params=phil)

        cmd = CommandManager('phenix.table_one')
        cmd.add_command_line_arguments([output_eff])
        self.log.bar()
        self.log("Table One Command:")
        self.log.bar()
        cmd.print_settings()
        cmd.run()
        cmd.write_output(output_eff.replace('.eff','.log'))

        if not os.path.exists(output_csv):
            raise Failure('Failed to make table one: {}'.format(output_csv))

        return output_csv



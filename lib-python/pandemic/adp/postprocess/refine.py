import os, collections
from libtbx import adopt_init_args
from libtbx.utils import Sorry, Failure
from bamboo.common.logs import Log
from bamboo.common.path import easy_directory

from pandemic.adp.parallel import RunParallelWithProgressBarUnordered

import traceback
def wrapper_run(obj):
    try:
        ret = obj.run()
    except:
        tr = traceback.format_exc()
        return tr
    return obj


class RefineStructures:


    def __init__(self,
            output_directory,
            refinement_program = 'phenix',
            #refinement_params = None,
            n_cpus = 1,
            verbose = False,
            log = None,
            ):

        if log is None: log = Log()

        from giant.xray.refine import refine_phenix, refine_refmac
        if  refinement_program == 'refmac':
            refine = refine_refmac
        elif refinement_program == 'phenix':
            refine = refine_phenix
        else:
            raise Failure('Invalid refinement program selected: {}'.format(refinement_program))

        run_parallel = RunParallelWithProgressBarUnordered(
            function = wrapper_run,
            n_cpus = n_cpus,
        )

        adopt_init_args(self, locals())

    def __call__(self,
        input_structures,
        input_reflection_data,
        labels,
        cif_files = None,
        output_suffix = '-refined',
        ):

        n = len(labels)

        assert n == len(input_structures)
        assert n == len(input_reflection_data)

        output_directories = [easy_directory(os.path.join(self.output_directory, l)) for l in labels]
        output_prefixes = [os.path.join(d, l+output_suffix) for d,l in zip(output_directories, labels)]

        arg_list = []
        for i in xrange(n):

            i_pdb = input_structures[i]
            i_mtz = input_reflection_data[i]

            if not os.path.exists(i_mtz):
                raise Failure('Something has gone wrong -- Trying to refine output structures but mtz has not been provided/has been lost.')

            # Create refinement object
            obj = self.refine(
                    pdb_file = i_pdb,
                    mtz_file = i_mtz,
                    cif_files = (cif_files if cif_files else None),
                    out_prefix = output_prefixes[i],
                    log = self.log,
                    ) \
            .set_refine_coordinates_only() \
            .set_n_cycles(5)

            obj.tag = labels[i]

            if os.path.exists(obj.out_pdb_file):
                raise Failure('Refined PDB already exists! (model {})'.format(obj.tag))

            arg_list.append(obj)

        for o in arg_list:
            o.log('')
            o.print_settings()

        # Refine all of the models
        self.log.subheading('Running {} refinements'.format(len(arg_list)))
        results = self.run_parallel(arg_dicts=arg_list)

        output_structures = collections.OrderedDict()
        for r in results:
            if isinstance(r, str):
                self.log.subheading('Failed to refine structure -- continuing anyway')
                self.log(r)
                continue

            output_structures[r.tag] = (r.out_pdb_file, r.out_mtz_file)

        return output_structures



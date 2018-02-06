#!/usr/bin/env ccp4-python

import os, sys, copy, re, shutil

import libtbx.phil
import libtbx.easy_mp

from bamboo.common.path import splice_ext
from bamboo.common.command import CommandManager
from giant.xray.crystal import CrystalSummary

############################################################################

PROGRAM = 'giant.datasets.prepare'

DESCRIPTION = """
    Take input MTZs, reindex, transfer R-free flags from a reference, fill missing reflections, and optionally run refinement.
"""

############################################################################

blank_arg_prepend = {'.mtz':'mtz=', '.pdb':'reference_pdb='}

master_phil = libtbx.phil.parse("""
input {
    mtz = None
        .type = path
        .multiple = True
    reference_mtz = None
        .type = path
    reference_pdb = None
        .type = path
}
options {
    free_r_flag = FreeR_flag
        .type = str
}
settings {
    cpus = 1
        .help = "number of cpus to use for processing"
        .type = int
}
""")

############################################################################

def reindex_mtz_to_reference(AAA):

        reind_output_mtz = os.path.splitext(mtz_file)[0]+'.reindexed.mtz'
        cad_output_mtz   = os.path.splitext(mtz_file)[0]+'.cad.mtz'
        free_output_mtz  = os.path.splitext(mtz_file)[0]+'.free.mtz'
        free_output_log  = os.path.splitext(free_output_mtz)[0]+'.log'

        cmds_list = []
        commands.append(cmds_list)

        ######################################
        # REINDEX
        ######################################
        cm = CommandManager('pointless')
        cm.add_command_line_arguments([ 'hklin', mtz_file,
                                        'hklref', params.input.reference_mtz,
                                        'hklout', reind_output_mtz   ])
        cm.add_standard_input([ 'tolerance 5' ])
        cm.print_settings()

        cmds_list.append(cm)

        ######################################
        # Transfer Free-R Flags
        ######################################
        cm = CommandManager('cad')
        cm.add_command_line_arguments([ 'hklin1', reind_output_mtz,
                                        'hklin2', params.input.reference_mtz,
                                        'hklout', cad_output_mtz   ])
        cm.add_standard_input([ 'labin file_number 1 E1=IMEAN E2=SIGIMEAN',
                                'labin file_number 2 E1={}'.format(params.options.free_r_flag),
                                'labout file_number 2 E1=FreeR_flag',
                                'END'     ])
        cm.print_settings()

        cmds_list.append(cm)

        ######################################
        # Fill in missing reflections (nans)
        ######################################
        cm = CommandManager('uniqueify')
        cm.add_command_line_arguments([ '-f', 'FreeR_flag',
                                        cad_output_mtz,
                                        free_output_mtz   ])
        cm.print_settings()

        cmds_list.append(cm)

def wrapper_run(cmds):
    for c in cmds:
        ret_code = c.run()
        if ret_code != 0:
            print '============================>'
            print 'Process returned with an error'
            print '============================>'
            c.print_settings()
            print '============================>'
            print c.output
            print '============================>'
            print c.error
            print '============================>'
            break
    return ret_code

def wrapper_reprocess(arg_dict):

    params = arg_dict['params']
    mtz    = arg_dict['mtz']

    out_files = []

    ###########################################
    #               Reindexing                #
    ###########################################

    try:
        new_file = splice_ext(path=in_mtz, new='reindexed')
        ret = reindex_mtz_to_reference(in_mtz  = mtz,
                                       out_mtz = new_file,
                                       reference_mtz=params.input.reference_mtz)
        assert os.path.exists(new_file), 'reindexing has failed -- {} does not exist'.format(new_file)
        out_files.append(new_file)
    except Exception as e:
        return Exception('Failed during reindexing: {} ({})'.format(mtz, str(e)))


def run(params):

    assert params.input.reference_mtz is not None, 'no reference mtz supplied (provided with reference_mtz=XXX)'
    assert len(params.input.mtz) > 0, 'no mtz files supplied'

    arg_list = [{'params':params, 'mtz':m} for m in params.input.mtz]

    print '============================>'
    print 'Running re-processing for {} mtz files'.format(len(arg_list))
    returned = libtbx.easy_mp.pool_map(fixed_func   = wrapper_reprocess,
                                       args         = arg_list,
                                       processes    = params.settings.cpus,
                                       chunksize    = 1)

    print '============================>'
    print 'Return code total: {}'.format(sum(returned))


    # Reindex MTZs and transfer R-free flags
    reindex_mtzs(reference_mtz  = params.input.reference_mtz,
                 target_mtzs    = params.input.mtzs)



    commands = []

    for mtz_file in params.input.mtz:

    # Clean up
    for mtz_file in params.input.mtz:

        reind_output_mtz = os.path.splitext(mtz_file)[0]+'.reindexed.mtz'
        cad_output_mtz   = os.path.splitext(mtz_file)[0]+'.cad.mtz'
        free_output_mtz  = os.path.splitext(mtz_file)[0]+'.free.mtz'
        free_output_log  = os.path.splitext(free_output_mtz)[0]+'.log'

        if os.path.exists(os.path.basename(free_output_log)) and (os.path.basename(free_output_log) != free_output_log):
            shutil.move(os.path.basename(free_output_log), free_output_log)

        if os.path.exists(free_output_mtz):
            os.remove(reind_output_mtz)
            os.remove(cad_output_mtz)
            os.remove(free_output_log)

############################################################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
        program             = PROGRAM,
        description         = DESCRIPTION)

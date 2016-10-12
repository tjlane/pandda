import os, sys, glob

import libtbx.phil

from bamboo.common.command import CommandManager
from bamboo.common.path import rel_symlink

#######################################

blank_arg_prepend = {   '.mtz':'mtz='   }

master_phil = libtbx.phil.parse("""
input {
    mtz = 'refine.mtz'
        .multiple = True
        .type = str
}
settings {
    f_obs = F,SIGF
        .type = str
    f_calc = F-model,PHIF-model
        .type = str
}

""")

#######################################

def run(params):

    assert params.input.mtz, 'No MTZs given for comparison'

    assert not os.path.exists('merged.mtz'), 'A file called merged.mtz already exists. Please delete it before re-running.'

    merge = CommandManager('giant.mtz.merge')
    merge.add_command_line_arguments( 'label_suffix=incremental' )
    merge.add_command_line_arguments( params.input.mtz )
    merge.add_command_line_arguments( params.settings.f_obs.split(',') )
    merge.add_command_line_arguments( params.settings.f_calc.split(',') )
    merge.run()
    merge.write_output('merged.log')

    assert os.path.exists('merged.mtz')

    fo1, fo2 = params.settings.f_obs.split(',')
    fc1, fc2 = params.settings.f_calc.split(',')

    for i_2 in range(1, len(params.input.mtz)+1):
        for i_1 in range(1, len(params.input.mtz)+1):
            if i_1 == i_2:
                break
            match = CommandManager('cphasematch')
            match.add_command_line_arguments( ['-mtzin', 'merged.mtz'] )
            match.add_command_line_arguments( ['-mtzout', 'phasematch.mtz'] )
            match.add_command_line_arguments( ['-colin-fo',   '{}-{},{}-{}'.format(fo1,i_1,fo1,i_1)] )
            match.add_command_line_arguments( ['-colin-fc-1', '{}-{},{}-{}'.format(fc1,i_1,fc1,i_1)] )
            match.add_command_line_arguments( ['-colin-fc-2', '{}-{},{}-{}'.format(fc2,i_2,fc2,i_2)] )
            match.run()
            march.write_output('phasematch-{}-{}.log'.format(i_1.i_2))


#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)

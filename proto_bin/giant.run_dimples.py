#!/usr/bin/env pandda.python

import os, sys, copy, re

import libtbx.phil
from Bamboo.Common.command import commandManager
from Giant.jiffies import parse_phil_args
from Giant.Xray.Data import crystalSummary

blank_arg_prepend = 'mtz='

master_phil = libtbx.phil.parse("""
input {
    reference_pdb = None
        .type = path
        .multiple = True
    mtz = None
        .type = path
        .multiple = True
}
output {
    dir_suffix = '_dimple'
        .type = str
}
""")

def run(params):

    for ref_pdb_file in params.input.reference_pdb:
        ref_mtz_file = ref_pdb_file.replace('.pdb','.mtz')
        assert os.path.exists(ref_mtz_file), 'Reference MTZ does not exist: {}'.format(ref_mtz_file)

        output_dirname = os.path.basename(ref_pdb_file).replace('.pdb', params.output.dir_suffix)
        print 'Placing all dimple runs for {} into directories called {}'.format(ref_pdb_file, output_dirname)

        for mtz_file in params.input.mtz:

            output_dir = os.path.join(os.path.dirname(mtz_file), output_dirname)

            cm = commandManager('dimple')
            cm.add_command_line_arguments([ mtz_file,
                                            ref_pdb_file,
                                            output_dir,
                                            '--free-r-flags', ref_mtz_file
                                        ])
            print '============================>'
            print 'Running Dimple:'
            cm.print_settings()
            ret_code = cm.run()
            if ret_code != 0:
                print '============================>'
                print 'Dimple returned with an error'
                print '============================>'
                print cm.output
                print '============================>'
                print cm.error
                print '============================>'

if __name__ == '__main__':

    # Show Defaults (just values)
    if '--show-defaults' in sys.argv:
        master_phil.show(attributes_level=0)
    # Show Defaults (including information)
    elif '--help' in sys.argv:
        master_phil.show(attributes_level=2)
    # ... or just run ...
    elif '--expert' in sys.argv:
        master_phil.show(attributes_level=4)
    # ... or just run ...
    else:
        working_phil = parse_phil_args(master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)
        out = run(params=working_phil.extract())
    # Exit (unnecessary, but eh)
    sys.exit()


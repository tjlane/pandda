import os, sys, copy, glob
import shutil

import libtbx.phil

import numpy

############################################################################

master_phil = libtbx.phil.parse("""
input {
    pandda_dir = './pandda'
        .help = 'Path to the pandda directory to use as a template'
        .type = str
}
output {
    out_dir = './pandda-migrate'
        .help = 'output folder name'
        .type = str
}

verbose = False
    .type = bool

""")

############################################################################

def copy_directory_contents(in_dir, out_dir, regex='*'):

    if isinstance(regex, str):
        regex = [regex]
    else:
        assert isinstance(regex, list), 'regex must either be str or list of str'
        for r in regex: assert isinstance(r, str), 'regex must either be str or list of str'

    if not os.path.exists(out_dir): os.mkdir(out_dir)

    for file_regex in regex:
        for file_in in sorted(glob.glob(os.path.join(in_dir, file_regex))):
            # Create output filename
            file_out = os.path.join(out_dir, os.path.basename(file_in))
            # Copy file
            print 'Copying {} -> {}'.format(file_in, file_out)
            shutil.copy(file_in, file_out)

def migrate_pandda(params):

    # Create output directory
    os.mkdir(params.output.out_dir)

    # Copy across reference structures (copying the file if it's a softlink)
    copy_directory_contents(
        in_dir  = os.path.join(params.input.pandda_dir, 'reference'),
        out_dir = os.path.join(params.output.out_dir,   'reference'),
        regex   = ['*.pdb', '*.mtz', '*.ccp4']
    )

    # Copy across pickles (copying the file if it's a softlink)
    copy_directory_contents(
        in_dir  = os.path.join(params.input.pandda_dir, 'pickled_panddas'),
        out_dir = os.path.join(params.output.out_dir,   'pickled_panddas'),
        regex   = ['reference*.pickle', 'statistical*.pickle']
    )

def run(params):

    assert params.input.pandda_dir, 'Must specify pandda directory'
    assert params.output.out_dir,   'Must specify output directory'
    assert not os.path.exists(params.output.out_dir), 'Output directory already exists'

    print 'Migrating Pandda from {} to {}'.format(params.input.pandda_dir, params.output.out_dir)
    migrate_pandda(params)

#######################################

if __name__ == '__main__':
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)

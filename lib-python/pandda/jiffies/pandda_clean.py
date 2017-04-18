import os, sys, copy, glob
import shutil

import libtbx.phil

import numpy

from bamboo.common.path import delete_with_glob

############################################################################

PROGRAM = 'pandda.clean'

DESCRIPTION = """
Cleaning Levels:
BASIC       - Remove unnecessary files
NORMAL      - Remove easily regenerated dataset files - re-run pandda to remake: quick
RUTHLESS    - Remove almost all dataset files         - re-run pandda to remake: slow
DEVASTATING - Remove almost everything                - re-run pandda to remake: very slow

Running with DEVASTATING essentially requires the pandda to be run again from scratch.

Modelled Stuctures (in the modelled_structures folders) are NEVER deleted by pandda.clean.
Files in the "analyses" folders are NEVER deleted by pandda.clean. Directories are also
NEVER deleted (only files are deleted).
"""

############################################################################

master_phil = libtbx.phil.parse("""
input {
    pandda_dir = None
        .help = 'Path to the pandda directory'
        .type = str
}
options {
    level = *basic normal ruthless devastating
        .type = choice
        .multiple = False
    skip_modelled = True
        .help = 'Clean interesting datasets as well?'
        .type = bool
}
settings {
    verbose = False
        .type = bool
    YES_I_AM_SURE = False
        .type = bool
}
""")

blank_arg_prepend = {None:'pandda_dir='}

############################################################################

def clean_misc(top_dir, level=0):
    print '============================================================>>>'
    print 'Cleaning Miscellaneous:'

    # Level 2 cleaning
    if level >= 2:
        delete_with_glob(os.path.join(top_dir, 'aligned_structures', '*.pdb'))

def clean_datasets(top_dir, level=0, skip_modelled=True):
    print '============================================================>>>'
    print 'Cleaning Dataset Folders:'

    for d_dir in sorted(glob.glob(os.path.join(top_dir, 'processed_datasets', '*'))):

        d_name = os.path.basename(d_dir)
        print '=======================================>>>'

        if skip_modelled and glob.glob(os.path.join(d_dir, 'modelled_structures', '*.pdb')):
            print 'Not cleaning modelled dataset (modelled_structures contains a pdb file): {}'.format(d_name)
            continue
        else:
            print 'Cleaning folder: {}'.format(d_dir)

        # Level 1 cleaning - clean maps
        if level >= 1:
            delete_with_glob(os.path.join(d_dir, '*.ccp4'))

        # Level 2 cleaning - clean structures and pickles as well
        if level >= 2:
            delete_with_glob(os.path.join(d_dir, '*.csv'))
            delete_with_glob(os.path.join(d_dir, '*.pdb'))
            delete_with_glob(os.path.join(d_dir, '*.mtz'))
            delete_with_glob(os.path.join(d_dir, 'pickles', '*.pickle'))
            delete_with_glob(os.path.join(d_dir, 'images', '*.png'))
            delete_with_glob(os.path.join(d_dir, 'graphs', '*.png'))
            delete_with_glob(os.path.join(d_dir, 'scripts', '*'))

    return

def clean_statistical_maps(top_dir, level=0):
    print '============================================================>>>'
    print 'Cleaning Statistical Maps:'

    # Level 0 cleaning - clean skew and kurt maps
    if level >= 0:
        delete_with_glob(os.path.join(top_dir, 'reference', 'statistical_maps', '*-skew_map.ccp4'))
        delete_with_glob(os.path.join(top_dir, 'reference', 'statistical_maps', '*-kurt_map.ccp4'))
        delete_with_glob(os.path.join(top_dir, 'reference', 'statistical_maps', '*-bimo_map.ccp4'))
        delete_with_glob(os.path.join(top_dir, 'reference', 'statistical_maps', '*-stds_map.ccp4'))

    # Level 1 cleaning - clean mean and sadj maps
    if level >= 1:
        delete_with_glob(os.path.join(top_dir, 'reference', 'statistical_maps', '*-mean_map.ccp4'))
        delete_with_glob(os.path.join(top_dir, 'reference', 'statistical_maps', '*-sadj_map.ccp4'))

    return

def clean_pickles(top_dir, level=0):
    print '============================================================>>>'
    print 'Cleaning Pickles:'

    # Level 2 cleaning - clean meta
    if level >= 2:
        delete_with_glob(os.path.join(top_dir, 'pickled_data', 'dataset_meta.pickle'))
        delete_with_glob(os.path.join(top_dir, 'pickled_data', 'reference_dataset.pickle'))

    # Level 3 cleaning - delete everthing
    if level >= 3:
        delete_with_glob(os.path.join(top_dir, 'pickled_data', '*.pickle'))

def run(params):

    assert params.input.pandda_dir, 'Must specify pandda directory'

    if params.options.level   == 'basic':
        ferocity = 0
    elif params.options.level == 'normal':
        ferocity = 1
    elif params.options.level == 'ruthless':
        ferocity = 2
    elif params.options.level == 'devastating':
        ferocity = 3

    if ferocity > 1 and (not params.settings.YES_I_AM_SURE):
        raise SystemExit('level set to {} -- are you sure? if yes, set YES_I_AM_SURE=True and re-run'.format(params.options.level))

    print 'Cleaning Pandda in folder: {}'.format(params.input.pandda_dir)
    print 'Running at Ferocity Level of {} ({})'.format(ferocity, params.options.level)

    clean_misc(
        top_dir=params.input.pandda_dir,
        level=ferocity
    )

    clean_datasets(
        top_dir=params.input.pandda_dir,
        level=ferocity,
        skip_modelled=params.options.skip_modelled
    )

    clean_statistical_maps(
        top_dir=params.input.pandda_dir,
        level=ferocity
    )

    clean_pickles(
        top_dir=params.input.pandda_dir,
        level=ferocity
    )

#######################################

if __name__ == '__main__':
    from pandda.jiffies import run_default
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
        program             = PROGRAM,
        description         = DESCRIPTION)

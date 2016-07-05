import os, sys, copy, glob
import shutil

import libtbx.phil

import numpy

from bamboo.common.path import delete_with_glob

############################################################################

master_phil = libtbx.phil.parse("""
input {
    pandda_dir = None
        .help = 'Path to the pandda directory to export files from'
        .type = str
}

level = *basic normal ruthless devastating
    .type = choice
    .multiple = False

skip_interesting = True
    .help = 'Clean interesting datasets as well?'
    .type = bool

verbose = False
    .type = bool

YES_I_AM_SURE = False
    .type = bool

""")

blank_arg_prepend = {None:'pandda_dir='}

help = """
Cleaning Levels:
BASIC       - Remove unnecessary files
NORMAL      - Remove easily regenerated dataset files - re-run pandda to remake: quick
RUTHLESS    - Remove almost all dataset files         - re-run pandda to remake: slow
DEVASTATING - Remove almost everything                - re-run pandda to remake: very slow

Running with DEVASTATING essentially requires the pandda to be run again from scratch.

Modelled Stuctures (in the modelled_structures folders) are NEVER deleted by pandda.clean.
Files in the "analyses" and "results_summaries" folders are NEVER deleted by pandda.clean.
Directories are also NEVER deleted (only files are deleted).
"""

############################################################################

def clean_misc(top_dir, level=0):
    print '============================================================>>>'
    print 'Cleaning Miscellaneous:'

    # Level 2 cleaning
    if level >= 2:
        delete_with_glob(os.path.join(top_dir, 'aligned_structures', '*.pdb'))

def clean_datasets(top_dir, level=0, skip_interesting=True):
    print '============================================================>>>'
    print 'Cleaning Dataset Folders:'

    for d_dir in sorted(glob.glob(os.path.join(top_dir, 'processed_datasets', '*'))):

        d_name = os.path.basename(d_dir)
        print '=======================================>>>'

        if skip_interesting and os.path.exists(os.path.join(top_dir, 'interesting_datasets', d_name)):
            print 'Not Cleaning Interesting Dataset: {}'.format(d_name)
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
            delete_with_glob(os.path.join(d_dir, '*_refine.params'))
            delete_with_glob(os.path.join(d_dir, 'pickles', '*.pickle'))
            delete_with_glob(os.path.join(d_dir, 'output_images', '*.png'))
            delete_with_glob(os.path.join(d_dir, 'ligand_files',  '*'))

    return

def clean_links(top_dir, level=0, skip_interesting=True):
    print '============================================================>>>'
    print 'Cleaning Soft Links:'

    # Level 0 cleaning
    if level >= 0:
        delete_with_glob(os.path.join(top_dir, 'resolutions', '*'))

    # Level 1 cleaning
    if level >= 1:
        pass

    # Level 2 cleaning
    if level >= 2:
        delete_with_glob(os.path.join(top_dir, 'empty_directories', '*'))
#        if not skip_interesting:
#            delete_with_glob(os.path.join(top_dir, 'interesting_datasets', '*'))

    return

def clean_statistical_maps(top_dir, level=0):
    print '============================================================>>>'
    print 'Cleaning Statistical Maps:'

    # Level 0 cleaning - clean skew and kurt maps
    if level >= 0:
        delete_with_glob(os.path.join(top_dir, 'statistical_maps', '*-skew_map.ccp4'))
        delete_with_glob(os.path.join(top_dir, 'statistical_maps', '*-kurt_map.ccp4'))
        delete_with_glob(os.path.join(top_dir, 'statistical_maps', '*-bimo_map.ccp4'))
        delete_with_glob(os.path.join(top_dir, 'statistical_maps', '*-stds_map.ccp4'))

    # Level 1 cleaning - clean mean and sadj maps
    if level >= 1:
        delete_with_glob(os.path.join(top_dir, 'statistical_maps', '*-mean_map.ccp4'))
        delete_with_glob(os.path.join(top_dir, 'statistical_maps', '*-sadj_map.ccp4'))

    return

def clean_pickles(top_dir, level=0):
    print '============================================================>>>'
    print 'Cleaning Pickles:'

    # Level 2 cleaning - clean meta
    if level >= 2:
        delete_with_glob(os.path.join(top_dir, 'pickled_panddas', 'dataset_meta.pickle'))
        delete_with_glob(os.path.join(top_dir, 'pickled_panddas', 'reference_dataset.pickle'))

    # Level 3 cleaning - delete everthing
    if level >= 3:
        delete_with_glob(os.path.join(top_dir, 'pickled_panddas', '*.pickle'))

def run(params):

    assert params.input.pandda_dir, 'Must specify pandda directory'

    if params.level   == 'basic':
        ferocity = 0
    elif params.level == 'normal':
        ferocity = 1
    elif params.level == 'ruthless':
        ferocity = 2
    elif params.level == 'devastating':
        ferocity = 3

    if ferocity > 1 and (not params.YES_I_AM_SURE):
        raise SystemExit('level set to {} -- are you sure? if yes, set YES_I_AM_SURE=True and re-run'.format(params.level))

    print 'Cleaning Pandda in folder: {}'.format(params.input.pandda_dir)
    print 'Running at Ferocity Level of {} ({})'.format(ferocity, params.level)

    clean_misc(
        top_dir=params.input.pandda_dir,
        level=ferocity
    )

    clean_datasets(
        top_dir=params.input.pandda_dir,
        level=ferocity,
        skip_interesting=params.skip_interesting
    )

    clean_links(
        top_dir=params.input.pandda_dir,
        level=ferocity,
        skip_interesting=params.skip_interesting
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
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)

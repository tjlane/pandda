import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys, copy, re, shutil

from giant.paths import easy_directory
from giant.exceptions import Sorry, Failure

from giant.plot import setup_plotting
setup_plotting()

############################################################################

PROGRAM = 'giant.datasets.cluster'

DESCRIPTION = """
    A tool to cluster sets of pdb and/or mtz files.

    1) Simple usage:
        > giant.datasets.cluster *.pdb labelling=filename
    or
        > giant.datasets.cluster */model.pdb labelling=foldername

"""

############################################################################

blank_arg_prepend = {
    '.pdb':'pdb=',
    '.mtz':'mtz=',
}

import libtbx.phil
master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .type = path
        .multiple = True
    mtz = None
        .type = path
        .multiple = True
    labelling = *filename foldername
        .type = choice(multi=False)
    labelling_regex = None
        .type = str
}
clustering {
    lcv_cutoff = 0.2
        .type = float
    label_nodes_above = 0.2
        .help = "Label nodes above this value"
        .type = float
}
output {
    out_dir = clustered_datasets
        .type = path
    file_mode = copy *symlink
        .type = choice
}
""")

def regex_label(pattern, string):
    import re
    matches = re.findall(pattern, string)
    if len(matches) == 0:
        raise Sorry('Regex did not match any part of string: \n\tPattern {}\n\tString {}'.format(pattern, string))
    if len(matches) > 1:
        raise Sorry('Regex matched multiple parts of string: \n\tPattern {}\n\tString {}\n\tMatches {}'.format(pattern, string, matches))
    match = matches[0]
    if not match:
        raise Sorry('Matched pattern is blank: \n\tPattern {}\n\tString {}\n\tMatch {}'.format(pattern, string, match))
    return match

def validate_labels(labels, filenames):
    if len(filenames) != len(labels):
        raise Failure('Different number of filenames and labels')
    logger(
        'Input filenames/labels: \n\t{}'.format(
            '\n\t'.join(
                ['{} : {}'.format(l,f) for l,f in zip(labels, filenames)]
            ),
        )
    )
    for l, f in zip(labels, filenames):
        if (not l):
            raise Sorry(
                'Invalid label generated: {} (from {}). \n\tConsider changing the labelling function.'.format(
                    l, f,
                )
            )
    if len(labels) != len(set(labels)):
        dups = [(l, f) for l, f in zip(labels, filenames) if labels.count(l) > 1]
        raise Sorry(
            'One or more labels are the same for the input files: \n\t{}\n\tConsider using a different labelling function'.format(
                '\n\t'.join(['{} : {}'.format(l, f) for l, f in sorted(dups)]),
            )
        )

def validate_params(params):
    # Validate input files
    if not (params.input.pdb or params.input.mtz):
        raise Sorry('No pdb/mtz files have been provided: specify with input.pdb or input.mtz')
    # Check and create output directory
    if os.path.exists(params.output.out_dir):
        raise Sorry('Output directory already exists!')
    if not params.output.out_dir:
        raise Sorry('No output directory has been specified: specify with output.out_dir')
    easy_directory(params.output.out_dir)

############################################################################

def run(params):

    validate_params(params)

    # Define and create image directory
    img_dir = easy_directory(os.path.join(params.output.out_dir, 'dendrograms'))

    # Setup logging
    logger = lg.setup_logging(
        name = __name__,
        log_file = os.path.join(params.output.out_dir, 'clustering.log'),
    )

    # Define output_file_function to copy or symlink files as needed
    if params.output.file_mode == 'symlink':
        out_file_func = os.symlink
    elif params.output.file_mode == 'copy':
        out_file_func = shutil.copy

    from giant.paths import filename, foldername
    if params.input.labelling_regex:
        import re
        label_func = lambda s: regex_label(pattern=params.input.labelling_regex, string=s)
    elif params.input.labelling == 'filename':
        label_func = filename
    elif params.input.labelling == 'foldername':
        label_func = foldername

    logger.heading('Processing input pdb/mtz files')

    pdb_files = params.input.pdb
    mtz_files = params.input.mtz

    n_pdb = len(pdb_files)
    n_mtz = len(mtz_files)
    n = n_pdb + n_mtz

    if (n_pdb > 0) and (n_mtz > 0):
        raise Sorry('Must provide either PDB files or MTZ files but not both')

    input_files = sorted(pdb_files + mtz_files)

    logger('Input: {} pdb(s) / {} mtz(s)\n'.format(n_pdb, n_mtz))

    labels = [label_func(f) for f in input_files]
    validate_labels(labels, input_files)

    label_hash = dict(zip(labels, input_files))

    # Load crystal summaries
    logger('\nReading input files')

    from giant.xray.crystal import CrystalInfo
    from giant.xray.crystal.cluster import CrystalGroup

    if (n_pdb > 0):
        crystal_summaries = [CrystalInfo.from_pdb(pdb_file=f, name=lab) for f,lab in zip(input_files, labels)]
    else:
        crystal_summaries = [CrystalInfo.from_mtz(mtz_file=f, name=lab) for f,lab in zip(input_files, labels)]

    logger.subheading('Crystal Summaries')

    for c in crystal_summaries:
        logger(str(c)+'\n')

    # Group by SpaceGroup
    logger.heading('Grouping crystals by space group')
    crystal_groups = CrystalGroup.by_space_group(
        crystals = crystal_summaries,
    )
    logger('Grouped crystals into {} space group(s)'.format(len(crystal_groups)))

    logger.heading('Analysing variation of unit cells for each space group')

    import numpy

    for cg in crystal_groups:

        sg_name = 'sg-{}'.format(cg.space_groups[0].split(' (')[0].replace(' ','_'))

        logger.subheading('Space Group {}: {} dataset(s)'.format(cg.space_groups[0], len(cg.crystals)))

        logger('Unit Cell Variation:')
        logger(numpy.round(cg.uc_stats.as_pandas_table().T, 2))

        logger('\nMaking unit cell dendrogram for all crystals with this spacegroup')
        if len(cg.crystals) > 1:
            filename = os.path.join(img_dir,'{}-all.png'.format(sg_name))
            logger('> {}'.format(filename))
            cg.dendrogram(
                fname = filename,
                xlab  = 'Crystal',
                ylab  = 'Linear Cell Variation',
                annotate_y_min = params.clustering.label_nodes_above,
            )

        logger('\nClustering {} unit cells using LCV cutoff of {}'.format(len(cg.crystals), params.clustering.lcv_cutoff))
        sg_crystal_groups = cg.by_unit_cell(cg.crystals, cutoff=params.clustering.lcv_cutoff)
        logger('Clustered crystals into {} group(s)'.format(len(sg_crystal_groups)))

        for i_cg2, cg2 in enumerate(sg_crystal_groups):

            cluster_name = '{}-cluster-{}'.format(sg_name, i_cg2+1)

            logger.subheading('Spacegroup {}, Cluster {}'.format(sg_name, i_cg2+1))

            logger('Unit Cell Variation:')
            logger(numpy.round(cg2.uc_stats.as_pandas_table().T, 2))

            logger('\nMaking unit cell dendrogram for this cluster of crystals')
            if len(cg2.crystals) > 1:
                filename = os.path.join(img_dir, '{}.png'.format(cluster_name))
                logger('> {}'.format(filename))
                cg2.dendrogram(
                    fname = filename,
                    xlab  = 'Crystal',
                    ylab  = 'Linear Cell Variation',
                    ylim  = (0, params.clustering.lcv_cutoff),
                    annotate_y_min = params.clustering.label_nodes_above,
                )

            logger('Copying files to output directory')

            # Go through and link the datasets for each of the spacegroups into a separate folder
            sub_dir = easy_directory(os.path.join(params.output.out_dir, cluster_name))

            for c in cg2.crystals:
                # Output directory for crystal
                c_dir = easy_directory(os.path.join(sub_dir, c.name))
                # Root path of input file without extension
                in_base = os.path.splitext(os.path.abspath(label_hash[c.name]))[0]
                # Output file base template
                out_base = os.path.join(c_dir, c.name)
                # Export file(s)
                for ext in ['.pdb', '.mtz']:
                    in_f = (in_base + ext)
                    out_f = (out_base + ext)
                    if os.path.exists(in_f):
                        out_file_func(in_f, out_f)

    logger.heading('finished')

############################################################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
        program             = PROGRAM,
        description         = DESCRIPTION,
    )

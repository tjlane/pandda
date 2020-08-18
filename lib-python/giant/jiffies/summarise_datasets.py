import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys, copy, re

from giant.exceptions import Sorry, Failure

############################################################################

PROGRAM = 'giant.summarise_datasets'

DESCRIPTION = """
    A tool to summarise the variation in sets of pdb and/or mtz files.
"""

############################################################################

blank_arg_prepend = {
    '.pdb':'pdb=',
    '.mtz':'mtz=',
}

import libtbx.phil
master_phil = libtbx.phil.parse("""
input {
    mtz = None
        .type = path
        .multiple = True
    pdb = None
        .type = path
        .multiple = True
    labelling = *filename foldername
        .type = choice
}
options {
    column_label = None
        .type = str
        .multiple = True
}
output {
    log_file = dataset_summary.log
        .type = path
}
""")

############################################################################

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

def format_column(col, col_width, position, fmt):
    COL_STR = '{:'+position+str(col_width)+fmt+'}'
    s = COL_STR.format(col)
    return s

def format_row(cols, row='', col_width=20, header=False):
    if header:
        fmt = ''
    else:
        fmt = '.3f'
    s = ('| '+format_column(row,col_width,'',''))*bool(row)+'|'+'|'.join([format_column(c,col_width,'^',fmt) for c in cols])+'|'
    if header:
        s = '|-'+('-'*col_width+'|')*(len(cols)+bool(row))+'\n'+s+'\n|-'+('-'*col_width+'|')*(len(cols)+bool(row))
    return s

def crystal_statistics(title, crystals, value_func, header=False, footer=False):
    from scitbx.array_family import flex
    from scitbx.python_utils import robust_statistics
    sorted_crystals = sorted(crystals, key=value_func)
    sorted_values = flex.double([value_func(c) for c in crystals])
    min_max_mean = sorted_values.min_max_mean()
    stddev = sorted_values.sample_standard_deviation()
    hinges = robust_statistics.hinges(sorted_values)
    median = robust_statistics.median(sorted_values)
    s_out = []
    if header:
        s_out.append(
            format_row(cols=['Min', 'Lower Quartile', 'Median', 'Upper Quartile', 'Max', 'Mean', 'Standard Dev'], row=' ', header=True)
        )
    s_out.append(
        format_row(cols=[min_max_mean.min, hinges[0], median, hinges[1], min_max_mean.max, min_max_mean.mean, stddev], row=title, header=False)
    )
    if footer:
        s_out.append(
            format_row(cols=['Min', 'Lower Quartile', 'Median', 'Upper Quartile', 'Max', 'Mean', 'Standard Dev'], row=' ', header=True)
        )
    return '\n'.join(s_out)

def crystal_min_max(title, crystals, label_hash, value_func):
    sorted_crystals = sorted(crystals, key=value_func)
    s_out = []
    smallest = sorted_crystals[0]
    largest = sorted_crystals[-1]
    s_out.append(
        'Smallest {}: {} - {} / {}'.format(title, value_func(smallest), smallest.name, label_hash[smallest.name])
    )
    s_out.append(
        'Largest {}:  {} - {} / {}'.format(title, value_func(largest), largest.name, label_hash[largest.name])
    )
    return '\n'.join(s_out)

############################################################################

def run(params):

    import iotbx.mtz
    from giant.xray.crystal import CrystalInfo
    from giant.xray.crystal.cluster import CrystalGroup

    # Create logging object
    logger = lg.setup_logging(
        name = __name__,
        log_file = params.output.log_file,
    )

    from giant.paths import filename, foldername
    if params.input.labelling =='filename':
        label_func = filename
    elif params.input.labelling =='foldername':
        label_func = foldername
    else:
        raise Failure(
            'labelling function not supported: {}'.format(
                params.input.labelling,
            )
        )

    # Process PDBs
    if params.input.pdb:

        filenames = params.input.pdb
        n = len(filenames)

        logger.heading('Processing {} PDB Files'.format(n))

        labels = [label_func(f) for f in filenames]
        validate_labels(labels=labels, filenames=filenames)
        label_hash = dict(zip(labels, filenames))

        crystals = [
            CrystalInfo.from_pdb(
                pdb_file = f,
                name = l,
            ) for f,l in zip(filenames, labels)
        ]

        logger.subheading('Crystal information')
        for c in crystals:
            logger(str(c)+'\n')

        logger.heading('Grouping {} pdb files by space group'.format(n))

        crystal_groups = CrystalGroup.by_space_group(crystals=crystals)

        logger('> Clustered into {} space group(s)'.format(len(crystal_groups)))

        for cg in crystal_groups:

            logger.subheading(
                'Space group: {} - {} datasets'.format(
                    ','.join(cg.space_groups),
                    len(cg.crystals),
                )
            )

            logger(crystal_statistics('Resolution', cg.crystals, value_func=lambda c: c.resolution_high, header=True))

            logger(crystal_statistics('R-work', cg.crystals, value_func=lambda c: c.r_work, header=True))
            logger(crystal_statistics('R-free', cg.crystals, value_func=lambda c: c.r_free, footer=True))

            logger('\nSmallest + Largest Values\n')

            logger(crystal_min_max('Resolution', cg.crystals, label_hash=label_hash, value_func=lambda c: c.resolution_high))
            logger(crystal_min_max('R-work', cg.crystals, label_hash=label_hash, value_func=lambda c: c.r_work))
            logger(crystal_min_max('R-free', cg.crystals, label_hash=label_hash, value_func=lambda c: c.r_free))

    # Process MTZs
    if params.input.mtz:

        filenames = params.input.mtz
        n = len(filenames)

        logger.heading('Processing {} MTZ Files'.format(n))

        labels = [label_func(f) for f in filenames]
        validate_labels(labels=labels, filenames=filenames)
        label_hash = dict(zip(labels, filenames))

        crystals = [
            CrystalInfo.from_mtz(
                mtz_file = f,
                name = l,
            ) for f,l in zip(filenames, labels)
        ]

        logger.subheading('Crystal information')
        for c in crystals:
            logger(str(c)+'\n')

        logger.heading('Grouping {} mtz files by space group'.format(n))

        crystal_groups = CrystalGroup.by_space_group(crystals=crystals)

        logger('> Clustered into {} space group(s)'.format(len(crystal_groups)))

        for cg in crystal_groups:

            logger.subheading(
                'Space group {} - {} datasets'.format(
                    ','.join(cg.space_groups),
                    len(cg.crystals),
                )
            )

            error = False
            for c in cg.crystals:
                for label in params.options.column_label:
                    if label is None:
                        continue
                    if label not in c.column_labels:
                        logger('Checking: column "{}" not in diffraction data of {}. columns present are {}'.format(label, label_hash[c.name], c.column_labels))
                        error = True
            if (error is True):
                raise Sorry('There are datasets that do not contain the right columns.')

            logger('\nCrystal Statistics\n')

            logger(crystal_statistics('Resolution (high)',  cg.crystals, value_func=lambda c: c.resolution_high, header=True))
            logger(crystal_statistics('Resolution (low)',   cg.crystals, value_func=lambda c: c.resolution_low))
            logger(crystal_statistics('Unit cell - vol',    cg.crystals, value_func=lambda c: c.unit_cell.volume(), header=True))
            logger(crystal_statistics('Unit cell - a',      cg.crystals, value_func=lambda c: c.unit_cell.parameters()[0]))
            logger(crystal_statistics('Unit cell - b',      cg.crystals, value_func=lambda c: c.unit_cell.parameters()[1]))
            logger(crystal_statistics('Unit cell - c',      cg.crystals, value_func=lambda c: c.unit_cell.parameters()[2]))
            logger(crystal_statistics('Unit cell - alpha',  cg.crystals, value_func=lambda c: c.unit_cell.parameters()[3]))
            logger(crystal_statistics('Unit cell - beta',   cg.crystals, value_func=lambda c: c.unit_cell.parameters()[4]))
            logger(crystal_statistics('Unit cell - gamma',  cg.crystals, value_func=lambda c: c.unit_cell.parameters()[5]))

            wavelength_func = lambda c: iotbx.mtz.object(label_hash[c.name]).crystals()[1].datasets()[0].wavelength()
            logger(crystal_statistics('Wavelength',         cg.crystals, value_func=wavelength_func, header=True, footer=True))

            if params.options.column_label:
                column_labels = params.options.column_label
                n_lab = len(column_labels)
                for i, label in enumerate(column_labels):
                    if label is None:
                        continue
                    if (i == 0):
                        logger('\nNumber of reflections in each column:\n')
                    column_func = lambda c: iotbx.mtz.object(label_hash[c.name]).get_column(label).n_valid_values()
                    logger(crystal_statistics(
                        label,
                        cg.crystals,
                        value_func = column_func,
                        header = (i==0),
                        footer = (i+1==n_lab),
                    ))

            logger('\nSmallest + Largest Values\n')

            logger(crystal_min_max('Resolution', cg.crystals, label_hash=label_hash, value_func=lambda c: c.resolution_high))

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
        description         = DESCRIPTION)

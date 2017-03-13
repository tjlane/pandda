import os, sys, copy, re

import libtbx.phil

import numpy
import pandas

from giant.xray.crystal import CrystalSummary
from giant.xray.crystal.cluster import CrystalGroup

from scitbx.array_family import flex
from scitbx.python_utils import robust_statistics

#######################################

blank_arg_prepend = {'.mtz':'mtz=', '.pdb':'pdb='}

master_phil = libtbx.phil.parse("""
input {
    mtz = None
        .type = path
        .multiple = True
    pdb = None
        .type = path
        .multiple = True
    file_label = *filename folder
        .type = choice
}
check_for{
    label = None
        .type = str
        .multiple = True
}
summary {
    observations = F,SIGF
        .type = str
}
output {
    out_dir = None
        .type = path
}
""")

#######################################

def print_table_row(cols, row='', col_width=20, header=False):
    def format_column(col, col_width, position, format):
        COL_STR = '{:'+position+str(col_width)+format+'}'
        return COL_STR.format(col)
    if header: format = ''
    else:      format = '.3f'
    if header: print '|-'+('-'*col_width+'|')*(len(cols)+bool(row))
    print ('| '+format_column(row,col_width,'',''))*bool(row)+'|'+'|'.join([format_column(c,col_width,'^',format) for c in cols])+'|'
    if header: print '|-'+('-'*col_width+'|')*(len(cols)+bool(row))

def print_crystal_statistics(label, crystals, value_func, header=False, footer=False):
    sorted_crystals = sorted(crystals, key=value_func)
    sorted_values   = flex.double([value_func(c) for c in crystals])
    min_max_mean = sorted_values.min_max_mean()
    stddev       = sorted_values.sample_standard_deviation()
    hinges       = robust_statistics.hinges(sorted_values)
    median       = robust_statistics.median(sorted_values)
    if header: print_table_row(cols=['Min', 'Lower Quartile', 'Median', 'Upper Quartile', 'Max', 'Mean', 'Standard Dev'], row=' ',   header=True)
    print_table_row(cols=[min_max_mean.min, hinges[0], median, hinges[1], min_max_mean.max, min_max_mean.mean, stddev],   row=label, header=False)
    if footer: print_table_row(cols=['Min', 'Lower Quartile', 'Median', 'Upper Quartile', 'Max', 'Mean', 'Standard Dev'], row=' ',   header=True)

def print_crystal_min_max(label, crystals, value_func):
    sorted_crystals = sorted(crystals, key=value_func)
    print label+' Smallest:', sorted_crystals[0].id,  sorted_crystals[0].pdb_file,  sorted_crystals[0].mtz_file
    print label+' Largest: ', sorted_crystals[-1].id, sorted_crystals[-1].pdb_file, sorted_crystals[-1].mtz_file
    print '----------------------------------->'

#######################################

def run(params):

    if params.output.out_dir and (not os.path.exists(params.output.out_dir)):
        os.mkdir(params.output.out_dir)
    # Dump images into the current directory if no directory is given
    if not params.output.out_dir: img_dir = './'
    else: img_dir = params.output.out_dir

    #########################################################################################################

    if params.input.mtz:

        print ''
        print '====================================================================================================>'
        print '====================================================================================================>'
        print '====================================================================================================>'
        print 'Processing {} MTZ Files'.format(len(params.input.mtz))
        print '====================================================================================================>'
        print '====================================================================================================>'
        print '====================================================================================================>'
        print ''

        if   params.input.file_label=='filename': labels = [os.path.basename(os.path.splitext(f)[0]) for f in params.input.mtz]
        elif params.input.file_label=='folder':   labels = [os.path.basename(os.path.dirname(f)) for f in params.input.mtz]
        else: raise Exception('MTZ labelling function not supported: {}'.format(params.input.file_label))

        print '============================================>'
        print 'GROUPING BY SPACE GROUP'
        crystal_groups = CrystalGroup.by_space_group(crystals=[CrystalSummary.from_mtz(mtz_file=f, id=lab) for f,lab in zip(params.input.mtz, labels)])
        print '> {} SPACE GROUP(S)'.format(len(crystal_groups))
        print '============================================>'

        for cg in crystal_groups:

            print ''
            print '===============================================>'
            print '======================================================================>'
            print '> SPACE GROUPS: {} - {} DATASETS'.format(','.join(cg.space_groups), len(cg.crystals))
            print '======================================================================>'
            print '===============================================>'
            print ''

            if params.check_for.label:
                for c in cg.crystals:
                    for lab in params.check_for.label:
                        if lab not in c.column_labels:
                            print 'COLUMN NOT IN FILE ({}): {}'.format(lab, c.mtz_file)

            print_crystal_statistics('WAVELENGTH',         cg.crystals, value_func=lambda c: c.mtz_object().crystals()[1].datasets()[0].wavelength(), header=True)
            print_crystal_statistics('RESOLUTION (HIGH)',  cg.crystals, value_func=lambda c: c.high_res,                                              header=False)
            print_crystal_statistics('RESOLUTION (LOW)',   cg.crystals, value_func=lambda c: c.low_res,                                               header=False)
            print_crystal_statistics('UNIT CELL - VOL',    cg.crystals, value_func=lambda c: c.unit_cell.volume(),                                    header=False)
            print_crystal_statistics('UNIT CELL - A',      cg.crystals, value_func=lambda c: c.unit_cell.parameters()[0],                             header=False)
            print_crystal_statistics('UNIT CELL - B',      cg.crystals, value_func=lambda c: c.unit_cell.parameters()[1],                             header=False)
            print_crystal_statistics('UNIT CELL - C',      cg.crystals, value_func=lambda c: c.unit_cell.parameters()[2],                             header=False)
            print_crystal_statistics('UNIT CELL - Alpha',  cg.crystals, value_func=lambda c: c.unit_cell.parameters()[3],                             header=False)
            print_crystal_statistics('UNIT CELL - Beta',   cg.crystals, value_func=lambda c: c.unit_cell.parameters()[4],                             header=False)
            print_crystal_statistics('UNIT CELL - Gamma',  cg.crystals, value_func=lambda c: c.unit_cell.parameters()[5],                             header=False, footer=True)

            if params.summary.observations:
                obs,sig = params.summary.observations.split(',')
                print_crystal_statistics('UNIQUE REFLECTIONS', cg.crystals, value_func=lambda c: c.mtz_object().extract_observations(obs,sig).data.size(),     header=False, footer=True)

            print ''
            print '================================================>'
            print 'Smallest + Largest Values'
            print '================================================>'

            print_crystal_min_max('RESOLUTION', cg.crystals, value_func=lambda c: c.high_res)

    #########################################################################################################
    #                                                                                                       #
    #########################################################################################################
    #                                                                                                       #
    #########################################################################################################

    if params.input.pdb:

        print ''
        print '====================================================================================================>'
        print '====================================================================================================>'
        print '====================================================================================================>'
        print 'Processing {} PDB Files'.format(len(params.input.pdb))
        print '====================================================================================================>'
        print '====================================================================================================>'
        print '====================================================================================================>'
        print ''

        if   params.input.file_label=='filename': labels = [os.path.basename(os.path.splitext(f)[0]) for f in params.input.pdb]
        elif params.input.file_label=='folder':   labels = [os.path.basename(os.path.dirname(f)) for f in params.input.pdb]
        else: raise Exception('PDB labelling function not supported: {}'.format(params.input.file_label))

        print '============================================>'
        print 'GROUPING BY SPACE GROUP'
        crystal_groups = CrystalGroup.by_space_group(crystals=[CrystalSummary.from_pdb(pdb_file=f, id=lab) for f,lab in zip(params.input.pdb, labels)])
        print '> {} SPACE GROUP(S)'.format(len(crystal_groups))
        print '============================================>'

        for cg in crystal_groups:

            print ''
            print '===============================================>'
            print '======================================================================>'
            print '> SPACE GROUPS: {} - {} DATASETS'.format(','.join(cg.space_groups), len(cg.crystals))
            print '======================================================================>'
            print '===============================================>'
            print ''

#            from IPython import embed; embed()

            print_crystal_statistics('R-WORK',         cg.crystals, value_func=lambda c: c.pdb_input().get_r_rfree_sigma().r_work, header=True)
            print_crystal_statistics('R-FREE',         cg.crystals, value_func=lambda c: c.pdb_input().get_r_rfree_sigma().r_free, header=False, footer=True)

            print ''
            print '================================================>'
            print 'Smallest + Largest Values'
            print '================================================>'

#            print_crystal_min_max('RESOLUTION', cg.crystals, value_func=lambda c: c.high_res)

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)

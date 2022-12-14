import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys

############################################################################

PROGRAM = 'multi.table_ones'

DESCRIPTION = """
    Make a parameter eff file for phenix.table_one
"""

############################################################################

blank_arg_prepend = {
    '.pdb' : 'pdb=',
    None : 'dir=',
}

options_phil = """
    column_labels = 'F-obs,SIGF-obs'
        .type = str
        .multiple = False
    r_free_label = 'R-free-flags'
        .type = str
        .multiple = False
"""

import libtbx.phil
master_phil = libtbx.phil.parse("""
input  {
    pdb = None
        .type = path
        .multiple = True
    dir = None
        .type = path
        .multiple = True
    pdb_style = None
        .type = str
        .multiple = False
    mtz_style = None
        .type = str
        .multiple = False
    labelling = *foldername filename
        .type = choice
        .multiple = False
}
options {{options_phil}
}
output {
    parameter_file = table_one.eff
        .type = str
        .multiple = False
    output_basename = table_one_multiple
        .type = str
        .multiple = False
}
settings {
    cpus = 1
        .type = int
    verbose = True
        .type = bool
}
""".replace(
    '{options_phil}',
    options_phil,
))

############################################################################

structure_block = """
  structure {{
    name = "{label}"
    pdb_file = "{pdb}"
    mtz_file = "{mtz}"
    data_labels = "{columns}"
    r_free_flags_label = "{rfree}"
    wavelength = 1.0000
    cif_file = None
    cif_directory = "{folder}"
    data_type = *xray neutron
    unmerged_data = None
    unmerged_labels = None
    use_internal_variance = False
    count_anomalous_pairs_separately = False
  }}"""

output_eff = """table_one {{
  {structure_effs}
  processing {{
    re_compute_r_factors = True
    n_bins = 10
    ligand_selection = None
  }}
  multiprocessing {{
    nproc = {cpus}
    technology = *multiprocessing sge lsf pbs condor slurm
    qsub_command = None
  }}
  output {{
    directory = "./"
    job_title = "multi.table_ones"
    show_missing_fields = True
    format = *txt *csv *rtf
    base_name = "{out_base}"
    verbose = "True"
    text_field_separation = 2
  }}
}}
"""

############################################################################

def run(params):

    if (not params.input.pdb) and (not params.input.dir):
        raise IOError('Need to supply input.pdb OR input.dir')

    if (params.input.pdb and params.input.dir):
        raise IOError('Cannot supply input.pdb AND input.dir')

    if (params.input.pdb_style) and (not params.input.mtz_style):
        params.input.mtz_style = params.input.pdb_style.replace('.pdb','')+'.mtz'

    out_file = params.output.parameter_file
    if os.path.exists(out_file):
        raise IOError('Output file already exists: {}'.format(out_file))

    from giant.paths import foldername, filename
    if params.input.labelling == 'foldername':
        lab_func = foldername
    else:
        lab_func = filename

    # Build file list
    file_list = []

    # From PDB files directly
    for pdb in params.input.pdb:
        logger('Processing {}'.format(pdb))
        mtz = pdb.replace('.pdb','.mtz')
        if not os.path.exists(pdb):
            raise IOError('PDB does not exist: {}'.format(pdb))
        if not os.path.exists(mtz):
            raise IOError('MTZ does not exist: {}\nMTZs must be named the same as the pdb files (except for .mtz extension) when using input.pdb option'.format(mtz))
        file_list.append((pdb,mtz))

    # From directories
    from giant.paths import resolve_glob
    for d in sorted(params.input.dir):
        logger('Processing directory: {}'.format(d))
        try:
            pdb = resolve_glob(os.path.join(d, params.input.pdb_style), n=1)
            mtz = resolve_glob(os.path.join(d, params.input.mtz_style), n=1)
        except ValueError:
            logger('PDB style/MTZ style did not resolve to the correct number of files')
            raise
        file_list.append((pdb,mtz))

    out_str = ''
    for pdb,mtz in file_list:
        out_str += structure_block.format(
            label = lab_func(pdb),
            pdb = os.path.abspath(pdb),
            mtz = os.path.abspath(mtz),
            columns = params.options.column_labels,
            rfree = params.options.r_free_label,
            folder = os.path.dirname(os.path.abspath(pdb)),
        )

    output = output_eff.format(
        structure_effs = out_str,
        out_base = params.output.output_basename,
        cpus = params.settings.cpus,
    ).replace(
        '{{','{'
    ).replace(
        '}}','}'
    )

    logger.subheading('Table One parameters')
    logger(output)

    with open(out_file, 'w') as fh:
        fh.write(output)

#######################################

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

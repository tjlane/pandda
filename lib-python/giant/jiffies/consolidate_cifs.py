import giant.logs as lg
logger = lg.getLogger(__name__)

import sys
import pathlib as pl

from giant.phil import (
    log_running_parameters,
    )

from giant.mulch.dataset import (
    CrystallographicModel,
    )

from giant.refinement.restraints.cif_tools import (
    CifManager,
    )

from giant.refinement.restraints.consolidate import (
    ConsolidateRestraintsAndUpdateModel,
    )

############################################################################

PROGRAM = 'giant.consolidate_cifs'

DESCRIPTION = """
    A script to make using multiple CIF files and LINKs a little bit less painful.
"""

############################################################################

blank_arg_prepend = {
    '.pdb' : 'input.pdb=',
    '.cif' : 'input.cif=',
}

import libtbx.phil
master_phil = libtbx.phil.parse("""
input {
    cif = None
        .type = str
        .multiple = True
    pdb = None
        .type = str
}
output {
    cif = 'consolidated.cif'
        .help = 'output cif file'
        .type = str
    pdb = 'consolidated.pdb'
        .help = 'output pdb file'
        .type = str
    log = None
        .help = 'output log file'
        .type = str
}
settings {
    overwrite = False
        .type = bool
    verbose = False
        .type = bool
}
""")

############################################################################

def set_get_log(params):

    if params.output.log is not None:
        pass
    elif params.output.cif is not None:
        params.output.log = str(
            pl.Path(params.output.cif).with_suffix('.log')
            )
    else:
        # probably going to error, so set some default
        params.output.log = "consolidated.log"

    return params.output.log

def validate_params(params):

    for f in [
        params.output.pdb,
        params.output.cif,
        ]:
        if pl.Path(f).exists():
            raise Exception('Output file already exists: {}'.format(str(f)))

def run(params):

    logger = lg.setup_logging(
        name = __name__,
        log_file = set_get_log(params),
        debug = params.settings.verbose,
        )

    # Report
    logger.subheading(
        'Validating input parameters and input files'
        )

    validate_params(params)

    log_running_parameters(
        params = params,
        master_phil = master_phil,
        logger = logger,
    )

    logger.subheading(
        'Setting up...'
        )

    consolidate_restraints = ConsolidateRestraintsAndUpdateModel()

    #####

    logger.heading(
        'Reading input files'
        )

    cif_managers = [
        CifManager.from_file(p)
        for p in params.input.cif
        ]

    for i_cif, cif in enumerate(cif_managers):

        logger.subheading(
            'Cif {n}: {s}'.format(
                n = i_cif + 1,
                s = cif.cif_obj.source,
                )
            )

        logger(str(cif))

    model = (
        CrystallographicModel.from_file(
            params.input.pdb
            )
        if params.input.pdb is not None
        else None
        )

    if model is not None:

        logger.subheading(
            'Input Atomic Model'
            )

        try:
            logger(
                str(model)
                )
        except:
            pass

    #####

    logger.subheading(
        'Consolidating and updating restraints'
        )

    result = consolidate_restraints(
        cif_managers = cif_managers,
        model = model,
        )

    # Write output file
    logger.subheading(
        'Writing output structures'
        )

    result.write_pdb(params.output.pdb)
    result.write_cif(params.output.cif)

    logger.heading('cif consolidator finished normally')

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

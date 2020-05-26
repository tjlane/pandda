import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys, copy

from giant.exceptions import Sorry, Failure

############################################################################

PROGRAM = 'giant.standardise_model'

DESCRIPTION = """
    A tool to standardise a multi-conformer model so that no backbone discontinuities are created, etc.

    1) Simple usage:
        > giant.standardise_model my.pdb
"""

############################################################################

blank_arg_prepend = {'.pdb': 'input.pdb='}

import libtbx.phil
master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .help = 'Input PDB files - multiple can be provided'
        .type = str
        .multiple = True
}
options {
    pruning_rmsd = 0.1
        .help = 'rmsd at which to prune all alternate conformations'
        .type = float
}
output {
    suffix = '.standardised'
        .help = 'output pdb file suffix'
        .type = str
    log = 'standardise-model.log'
        .help = 'output log file'
        .type = str
}
settings {
    overwrite = True
        .type = bool
    verbose = False
        .type = bool
}
""")

############################################################################

from giant.structure.altlocs import expand_alternate_conformations, prune_redundant_alternate_conformations

def standardise_multiconformer_model(hierarchy, pruning_rmsd=0.1, in_place=False, verbose=False):
    """Standardise hierarchies by expanding alternate model conformations, and then trimming alternate conformations where possible"""

    # Alter the original files?
    if (in_place is False):
        # Copy the hierarchies
        hierarchy = hierarchy.deep_copy()

    # Sort the atoms
    hierarchy.sort_atoms_in_place()

    logger.subheading('Explicitly expanding model to all conformations of the crystal')
    expand_alternate_conformations(
        hierarchy = hierarchy,
        in_place = True,
        verbose = verbose,
    )

    logger.subheading('Pruning unneccessary multi-conformer residues in the expanded structure')
    prune_redundant_alternate_conformations(
        hierarchy = hierarchy,
        required_altlocs = hierarchy.altloc_indices(),
        rmsd_cutoff = pruning_rmsd,
        in_place = True,
        verbose = verbose,
    )

    return hierarchy

############################################################################

def run(params):

    logger = lg.setup_logging(
        name = __name__,
        log_file = params.output.log,
    )

    logger('Input Parameters:')
    logger(master_phil.format(params).as_str())

    # Report
    logger.heading('Validating input parameters and input files')

    # Check one or other have been provided
    if not params.input.pdb:
        raise IOError('No pdb files have been provided')
    for pdb in params.input.pdb:
        if not os.path.exists(pdb):
            raise IOError('pdb does not exist: {}'.format(pdb))

    for pdb in params.input.pdb:

        logger.heading('Standardising {}'.format(pdb))

        from giant.io.pdb import strip_pdb_to_input
        obj = strip_pdb_to_input(pdb, remove_ter=True)
        if obj.hierarchy.models_size() > 1:
            raise Sorry('Input structures may only have one model (current structure has {})'.format(obj.hierarchy.models_size()))

        # Merge the hierarchies
        final =  standardise_multiconformer_model(
            hierarchy = obj.hierarchy,
            pruning_rmsd = params.options.pruning_rmsd,
            in_place = True,
            verbose = params.settings.verbose,
        )

        # Update the atoms numbering
        final.sort_atoms_in_place()

        # Write output file
        filename = os.path.splitext(pdb)[0]+params.output.suffix+'.pdb'
        logger('Writing output structure to {}'.format(filename))
        final.write_pdb_file(file_name=filename, crystal_symmetry=obj.crystal_symmetry())

    logger.heading('finished')

    return

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

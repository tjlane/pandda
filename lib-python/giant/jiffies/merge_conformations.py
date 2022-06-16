import giant.logs as lg
logger = lg.getLogger(__name__)

import sys
import pathlib as pl

from giant.phil import (
    log_running_parameters,
    )

from giant.mulch.dataset import (
    AtomicModel,
    )

from giant.structure.conformers import (
    MakeMultiStateModel,
    )

from giant.structure.occupancy import (
    ScaleOccupancies,
    SanitiseOccupancies,
    ResetOccupancies,
    )

from giant.jiffies import (
    make_restraints,
    )

############################################################################

PROGRAM = 'giant.merge_conformations'

DESCRIPTION = """
    A tool to merge multiple models of the same crystal into one structure.
        The alternate conformers of the two models are generated and merged.
        Any common sections are then converted back into a single conformer
        model. This tool is therefore particularly useful when one model is
        an edited version of the other as many of the residues will be
        identical and the output structure will not contain too many
        alternate conformations.

    1) Simple usage:
        > giant.merge_conformations structure_1.pdb structure_2.pdb
"""

############################################################################

blank_arg_prepend = {
    '.pdb' : 'input.pdb=',
}

import libtbx.phil
master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .type = str
        .multiple = True
}
output {
    pdb = 'merge_conformations.pdb'
        .help = 'output pdb file'
        .type = str
    log = None
        .help = 'output log file'
        .type = str
}
options {
    prune_duplicates_rmsd = 0.05
        .help = 'RMSD below which to remove duplicate residues'
        .type = float
    occupancy = None
        .help = 'comma-separated list of occupancies for the different states. Used to scale the input occupancies of each structure.'
        .type = str
        .multiple = False
    sanitise_occupancies = True
        .type = bool
        .help = '(Re)scale occupancies so that group occupancies in the output structure are between 0 and 1.'
    reset_all_occupancies = False
        .help = 'Set all conformer occupancies to the same value (1.0/<number of conformers>). overrides all other occupancy arguments.'
        .type = bool
}
restraints {
    make_restraints = True
        .help = 'generate refinement restraints for the output model'
        .type = bool
    output {
    include scope giant.jiffies.make_restraints.output_phil
    }
    include scope giant.jiffies.make_restraints.options_phil
}
settings {
    overwrite = False
        .type = bool
    verbose = False
        .type = bool
}
""", process_includes=True)

############################################################################

def set_get_log(params):

    if params.output.log is not None:
        pass
    elif params.output.pdb is not None:
        params.output.log = str(
            pl.Path(params.output.pdb).with_suffix('.log')
            )
    else:
        # probably going to error, so set some default
        params.output.log = "merge_conformations.log"

    return params.output.log

def validate_params(params):

    params.input.pdb = [
        p for p in params.input.pdb
        if p is not None
        ]

    if len(params.input.pdb) == 0:
        raise IOError('No input PDBs provided.')

    if len(params.input.pdb) == 1:
        raise IOError('Only one input PDB provided.')

    for p in params.input.pdb:

        if not pl.Path(p).exists():

            raise IOError(
                'File does not exist: {path}'.format(
                    str(p)
                    )
                )

    if (params.output.pdb is None) or (not params.output.pdb.strip()):

        raise IOError('No output PDB file provided.')

    # Check existence of output pdb and delete as necessary
    if pl.Path(params.output.pdb).exists() and (not params.settings.overwrite):

        raise IOError(
            'Output file already exists: {path}.\nRun with overwrite=True to remove this file'.format(
                path = params.output.pdb,
                )
            )

    if (params.options.occupancy is not None) and len(params.options.occupancy) > 0:

        occupancies = list(map(float, params.options.occupancy.split(',')))

        if len(occupancies) != len(params.input.pdb):

            raise IOError(
                "Must provide the same number of occupancies as input PDB files: \n{pdbs}\n{occupancies}".format(
                    pdbs = len(params.input.pdb),
                    occupancies = len(params.options.occupancy),
                    )
                )

        for o in occupancies:

            if (o < 0.0) or (o > 1.0):

                raise IOError(
                    "Invalid occupancy (must be between 0.0 and 1.0): {o}".format(
                        o = o,
                        )
                    )

        occ_sum = sum(occupancies)

        if occ_sum > 1.0:

            raise IOError(
                "Occupancy sum must be less than 1.0: {occs}={occ_sum}".format(
                    occs = '+'.join(map(str, occupancies)),
                    occ_sum = occ_sum,
                    )
                )


class MergeConformations(object):

    def __init__(self,
        prune_rmsd_cutoff,
        sanitise_occupancies = False,
        reset_all_occupancies = False,
        ):

        self.scale_occupancies = ScaleOccupancies(
            in_place = True,
            )

        self.merge_hierarchies = MakeMultiStateModel(
            prune_rmsd_cutoff = prune_rmsd_cutoff,
            in_place = True,
            )

        self.sanitise_occupancies = (
            SanitiseOccupancies()
            if sanitise_occupancies
            else None
            )

        self.reset_occupancies = (
            ResetOccupancies()
            if reset_all_occupancies
            else None
            )

    def __call__(self,
        hierarchies,
        occupancies = None,
        ):

        # Copy as will be modifying the list
        hierarchies = list(hierarchies)

        if (occupancies is not None):

            logger.subheading(
                'Applying input occupancies prior to merging: \n\t{occs}'.format(
                    occs = str(occupancies),
                    )
                )

            assert len(occupancies) == len(hierarchies)

            for h, o in zip(hierarchies, occupancies):

                self.scale_occupancies(
                    atoms = h.atoms(),
                    multiplier = o,
                    )

        #####

        logger.heading(
            'Merging models'
            )

        hierarchy = self.merge_hierarchies(
            hierarchies = hierarchies,
            )

        #####

        logger.subheading(
            'Post-processing occupancies'
            )

        # Reset occupancies
        if self.reset_occupancies is not None:

            hierarchy = self.reset_occupancies(hierarchy)

        # Always do this last
        if self.sanitise_occupancies is not None:

            hierarchy = self.sanitise_occupancies(hierarchy)

        return hierarchy


class MakeOutputRestraints(object):

    def __init__(self,
        params,
        ):

        self.make_restraints_master_phil = (
            make_restraints.master_phil
            )

        self.scope = self.get_scope(
            params = params,
            )

    def __call__(self, pdb_path):

        params = self.scope.copy().extract()

        params.input.pdb = str(pdb_path)

        # Report
        logger.heading(
            'Making restraints for output structure'
            )
        #

        output = make_restraints.run(
            params
            )

        return output

    def get_scope(self, params):

        # Transfer the other phil objects from the master phil
        #
        master_phil = self.make_restraints_master_phil
        #
        merged_params = (
            master_phil.fetch(
                master_phil.format(
                    params.restraints
                    )
                )
            ).extract()

        # Transfer from the main params
        merged_params.settings.overwrite = params.settings.overwrite
        merged_params.settings.verbose = params.settings.verbose

        # Set some different defaults
        merged_params.output.log = None

        # Reformat to a scope (can be copied)
        merged_scope = master_phil.format(
            merged_params
            )

        return merged_scope


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

    merge_conformations = MergeConformations(
        prune_rmsd_cutoff = params.options.prune_duplicates_rmsd,
        sanitise_occupancies = params.options.sanitise_occupancies,
        reset_all_occupancies = params.options.reset_all_occupancies,
        )

    make_output_restraints = (
        MakeOutputRestraints(
            params = params,
            )
        if params.restraints.make_restraints
        else None
        )

    #####

    logger.subheading(
        'Reading input files'
        )

    models = [
        AtomicModel.from_file(p)
        for p in params.input.pdb
        ]

    logger(
        '\n\n'.join(map(str, models))
        )

    #####

    hierarchy = merge_conformations(
        hierarchies = [m.hierarchy for m in models],
        occupancies = (
            list(map(float, params.options.occupancy.split(',')))
            if params.options.occupancy
            else None
            ),
        )

    # Update the atom numbering
    hierarchy.sort_atoms_in_place()
    hierarchy.atoms_reset_serial()

    # Write output file
    logger.subheading(
        'Writing output structure'
        )

    logger(
        'Writing output structure to {}'.format(params.output.pdb)
        )

    model = models[0]
    hierarchy.write_pdb_file(
        file_name = params.output.pdb,
        crystal_symmetry = (
            model.input.crystal_symmetry()
            ),
    )

    if make_output_restraints is not None:

        restraints = make_output_restraints(
            pdb_path = params.output.pdb,
            )

    logger.heading('merge_conformations finished normally')

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

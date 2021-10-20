import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys
import pathlib as pl

from giant.phil import (
    log_running_parameters,
    )

from giant.mulch.dataset import (
    AtomicModel
    )

from giant.refinement.restraints import (
    DummyRestraintMaker,
    )

from giant.refinement.restraints.multi_state import (
    MakeMultiStateRestraints,
    MakeSimpleOccupancyRestraints,
    MakeMultiStateOccupancyRestraints,
    MakeIntraConformerRestraints,
    MakeDuplicateConformerRestraints,
    )

from giant.refinement.restraints.format import (
    WriteRestraints,
    )


############################################################################

PROGRAM = 'giant.make_restraints'

DESCRIPTION = """
    A tool to simplify the generation of restraints for use during refinement with REFMAC or PHENIX.

    The output ".params" files may be passed to giant.quick_refine for use in refinement.

    1) Simple usage:
        > giant.make_restraints input.pdb

    2) Just make simple occupancy restraints
        > giant.make_restraints input.pdb modes=simple_occupancy_groups
"""

############################################################################

blank_arg_prepend = {
    '.pdb' : 'pdb=',
}

input_phil = """
    pdb = None
        .help = 'Protein model'
        .type = str
"""
output_phil = """
    output_prefix = 'restraints'
        .help = 'output file root'
        .type = str
    output_formats = *refmac *phenix
        .type = choice(multi=True)
    log = 'restraints.log'
        .help = 'log file name'
        .type = path
"""
options_phil = """
modes = *local_altloc_restraints *duplicated_atom_restraints *simple_occupancy_groups *multi_state_occupancy_groups
    .type = choice(multi=True)
    .help = "Which restraints should be generated"
atom_selection = None
    .help = "Select a subset of atoms to generate restraints for"
    .type = str
exclude_hydrogens = True
    .help = "Exclude all hydrogen restraints"
    .type = bool
local_altloc_restraints {
    altloc_selection = None
        .help = "Select a subset of altlocs to generate restraints for e.g. ABC"
        .type = str
    max_distance = 4.2
        .help = "Maximum distance to create local restraints between atoms"
        .type = float
    min_distance = 0.1
        .help = "Minimum distance to create local restraints between atoms"
        .type = float
    sigma_xyz = 0.1
        .help = "Sigma of the distance restraint controlling how strongly the restraint is enforced"
        .type = float
    atom_selection = None
        .help = "Select a subset of atoms to generate restraints for"
        .type = str
}
duplicated_atom_restraints {
    rmsd_cutoff = 0.1
        .help = "Cutoff at which two conformers are considered to be duplicated"
        .type = float
    sigma_xyz = 0.02
        .help = "Coordinate restraint term controlling how strongly the restraint is enforced"
        .type = float
    atom_selection = None
        .help = "Select a subset of atoms to generate restraints for"
        .type = str
}
simple_occupancy_groups {
    single_atom_groups = True
        .type = bool
    single_conformer_groups = True
        .type = bool
    atom_selection = None
        .help = "Select a subset of atoms to generate restraints for"
        .type = str
}
multi_state_occupancy_groups {
    group_dist = 6.0
        .type = float
        .help = 'Distance to use when clustering groups of atoms that should have the SAME occupancy'
    overlap_dist = 3.0
        .type = float
        .help = 'Distance to use when clustering groups of atoms that should have occupancies that SUM TO (LESS THAN) ONE'
    ignore_common_solvent_molecules = True
        .type = bool
    include_resname = None
        .type = str
        .multiple = True
    ignore_resname = None
        .type = str
        .multiple = True
    set_group_completeness_to = None
        .help = 'Generate a set of fully constrained groups (that sum to unitary occupancy) when True. Generate a set of weaker constraints for overlapping atoms when False.'
        .type = bool
    atom_selection = None
        .help = "Select a subset of atoms to generate restraints for"
        .type = str
}
"""

import libtbx.phil
master_phil = libtbx.phil.parse("""
input {{input_phil}}
output {{output_phil}}
options {{options_phil}}
settings {
    overwrite = True
        .type = bool
    verbose = True
        .type = bool
}
""".replace(
    '{input_phil}',
    input_phil,
).replace(
    '{output_phil}',
    output_phil,
).replace(
    '{options_phil}',
    options_phil,
))

def validate_params(params):

    if not params.input.pdb:
        raise IOError('No PDB File Provided')

    if not pl.Path(params.input.pdb).exists():
        raise IOError(
            'File does not exist: {p}'.format(
                p = p,
                )
            )

def build_restraints_maker(options):

    lar = options.local_altloc_restraints
    dar = options.duplicated_atom_restraints
    sog = options.simple_occupancy_groups
    msog = options.multi_state_occupancy_groups

    modes = options.modes
    exclude_hydrogens = options.exclude_hydrogens

    maker = MakeMultiStateRestraints(
        make_intra_conformer_restraints = (
            MakeIntraConformerRestraints(
                min_distance_cutoff = lar.min_distance,
                max_distance_cutoff = lar.max_distance,
                distance_restraint_sigma = lar.sigma_xyz,
                select_altlocs = lar.altloc_selection,
                atom_selection = lar.atom_selection,
                exclude_hydrogens = exclude_hydrogens,
                )
            if ('local_altloc_restraints' in modes)
            else DummyRestraintMaker(name='< step skipped >')
            ),
        make_duplicate_conformer_restraints = (
            MakeDuplicateConformerRestraints(
                rmsd_cutoff = dar.rmsd_cutoff,
                distance_restraint_sigma = dar.sigma_xyz,
                atom_selection = dar.atom_selection,
                exclude_hydrogens = exclude_hydrogens,
                )
            if ('duplicated_atom_restraints' in modes)
            else DummyRestraintMaker(name='< step skipped >')
            ),
        make_simple_occupancy_restraints = (
            MakeSimpleOccupancyRestraints(
                single_atom_groups = sog.single_atom_groups,
                single_conformer_groups = sog.single_conformer_groups,
                atom_selection = sog.atom_selection,
                )
            if ('simple_occupancy_groups' in modes)
            else DummyRestraintMaker(name='< step skipped >')
            ),
        make_multi_state_occupancy_restraints = (
            MakeMultiStateOccupancyRestraints(
                group_distance_cutoff = msog.group_dist,
                overlap_distance_cutoff = msog.overlap_dist,
                ignore_common_solvent_molecules = msog.ignore_common_solvent_molecules,
                include_resnames_list = msog.include_resname,
                ignore_resnames_list = msog.ignore_resname,
                set_group_completeness_to = msog.set_group_completeness_to,
                atom_selection = msog.atom_selection,
                )
            if ('multi_state_occupancy_groups' in modes)
            else DummyRestraintMaker(name='< step skipped >')
            ),
        )

    return maker


def run(params):

    logger = lg.setup_logging(
        name = __name__,
        log_file = params.output.log,
        debug = params.settings.verbose,
        )

    validate_params(params)

    log_running_parameters(
        params = params,
        master_phil = master_phil,
        logger = logger,
    )

    ######################################################################
    # Prepare output and input
    ######################################################################

    overwrite = bool(params.settings.overwrite)
    input_pdb = params.input.pdb
    output_prefix = params.output.output_prefix
    output_formats = params.output.output_formats

    out_root = (
        pl.Path(
            str(output_prefix) + '.params'
            )
        if (output_prefix is not None)
        else pl.Path(
            input_pdb
            ).with_suffix('.params')
        )

    output_paths = {
        k : out_root.with_suffix('.'+str(k)+'.params')
        for k in output_formats
    }

    for k, p in output_paths.items():
        if p.exists():
            if (overwrite is True):
                os.remove(str(p))
            else:
                raise IOError(
                    'File already exists for {k}: {p}'.format(
                        k = str(k),
                        p = str(p),
                        )
                    )

    ######################################################################
    # Generate objects
    ######################################################################

    make_restraints = build_restraints_maker(params.options)

    write_restraints = WriteRestraints(
        formats = output_paths.keys(),
        output_path_dict = output_paths,
        )

    ######################################################################
    # Generate restraints
    ######################################################################

    model = AtomicModel.from_file(input_pdb)
    model.hierarchy.sort_atoms_in_place()

    restraint_collection = make_restraints(
        hierarchy = model.hierarchy,
        )

    logger.heading('Output Restraints')
    logger(str(restraint_collection))

    formatted_restraints = write_restraints(
        restraint_collection = restraint_collection,
        )

    logger.heading('done')

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

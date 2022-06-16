import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys
import pathlib as pl

from giant.phil import (
    log_running_parameters,
    )

from giant.mulch.dataset import (
    AtomicModel,
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

from giant.refinement.restraints.to_pymol import (
    WriteRestraintsPymolScript,
    )

from giant.refinement.restraints import (
    DummyRestraintMaker,
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
    output_root = None
        .help = 'manual output filepath root -- by default use the input pdb path as root'
        .type = str
    output_formats = *refmac *phenix
        .type = choice(multi=True)
    write_pymol_script = False
        .help = "Write restraints as pymol script for visualisation"
        .type = bool
    log = None
        .help = 'log file name'
        .type = path
"""
options_phil = """
modes = *multi_state_occupancy_groups *local_altloc_restraints *duplicated_atom_restraints *simple_occupancy_groups
    .type = choice(multi=True)
    .help = "Which restraints should be generated"
atom_selection = None
    .help = "Select a subset of atoms to generate restraints for"
    .type = str
exclude_hydrogens = True
    .help = "Exclude all hydrogen restraints"
    .type = bool
local_altloc_restraints {
    selection = all *multi_state_occupancy_groups
        .type = choice(multi=False)
        .help = "what atoms to generate restraints for (default just those that are selected by multi_state_occupancy_groups)"
    altloc_selection = None
        .help = "Select a subset of altlocs to generate restraints for e.g. ABC"
        .type = str
    max_distance = 4.0
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
    ignore_common_molecules = False
        .type = bool
    include_resname = None
        .type = str
        .multiple = True
    ignore_resname = None
        .type = str
        .multiple = True
    set_group_completeness_to = None
        .help = 'Generate a set of fully constrained groups (that sum to unitary occupancy) when True. Generate a set of weaker constraints for overlapping atoms when False. None will decide automatically.'
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
    overwrite = False
        .type = bool
    verbose = False
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

def set_get_log(params):

    if params.output.log is not None:
        pass
    elif params.output.output_root is not None:
        params.output.log = str(
            pl.Path(params.output.output_root).with_suffix('.log')
            )
    else:
        # set some sensible default
        params.output.log = "restraints.log"

    return params.output.log

def validate_params(params):

    if not params.input.pdb:
        raise IOError('No PDB File Provided')

    if not pl.Path(params.input.pdb).exists():
        raise IOError(
            'File does not exist: {p}'.format(
                p = params.input.pdb,
                )
            )

def build_restraints_maker(options):

    lar = options.local_altloc_restraints
    dar = options.duplicated_atom_restraints
    sog = options.simple_occupancy_groups
    msog = options.multi_state_occupancy_groups

    modes = options.modes
    exclude_hydrogens = options.exclude_hydrogens

    make_intra_conformer_restraints_all = DummyRestraintMaker(name='< step skipped >')
    make_intra_conformer_restraints_occ = None

    if 'local_altloc_restraints' in modes:

        make_intra_conformer_restraints = MakeIntraConformerRestraints(
            min_distance_cutoff = lar.min_distance,
            max_distance_cutoff = lar.max_distance,
            distance_restraint_sigma = lar.sigma_xyz,
            select_altlocs = lar.altloc_selection,
            atom_selection = lar.atom_selection,
            exclude_hydrogens = exclude_hydrogens,
            )

        if lar.selection == 'all':
            make_intra_conformer_restraints_all = make_intra_conformer_restraints
        elif lar.selection == 'multi_state_occupancy_groups':
            make_intra_conformer_restraints_occ = make_intra_conformer_restraints
        else:
            raise NotImplementedError()

    maker = MakeMultiStateRestraints(
        make_multi_state_occupancy_restraints = (
            MakeMultiStateOccupancyRestraints(
                group_distance_cutoff = msog.group_dist,
                overlap_distance_cutoff = msog.overlap_dist,
                ignore_common_molecules = msog.ignore_common_molecules,
                include_resnames_list = msog.include_resname,
                ignore_resnames_list = msog.ignore_resname,
                set_group_completeness_to = msog.set_group_completeness_to,
                atom_selection = msog.atom_selection,
                make_intra_conformer_distance_restraints = make_intra_conformer_restraints_occ,
                )
            if ('multi_state_occupancy_groups' in modes)
            else DummyRestraintMaker(name='< step skipped >')
            ),
        make_intra_conformer_restraints = (
            make_intra_conformer_restraints_all
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
        )

    return maker

def run(params):

    logger = lg.setup_logging(
        name = __name__,
        log_file = set_get_log(params),
        debug = params.settings.verbose,
        )

    log_running_parameters(
        params = params,
        master_phil = master_phil,
        logger = logger,
    )

    validate_params(params)

    ######################################################################
    # Prepare output and input
    ######################################################################

    overwrite = bool(
        params.settings.overwrite
        )

    input_pdb = pl.Path(
        params.input.pdb
        )

    output_root = (
        params.output.output_root
        )

    output_formats = (
        params.output.output_formats
        )

    ###

    output_root = (
        pl.Path(
            str(output_root)
            )
        if (
            output_root is not None
            )
        else (
            input_pdb
            ).with_name(
            input_pdb.stem+'-restraints'
            )
        )

    output_paths = {
        k : output_root.with_suffix(
            '.{k}.params'.format(k=k)
            )
        for k in output_formats
    }

    output_pymol_filepath = (
        output_root.with_suffix('.pymol.py')
        )
    if output_pymol_filepath.exists():
        if (overwrite is True):
            os.remove(str(output_pymol_filepath))

    for k, p in list(output_paths.items()):
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
        formats = list(output_paths.keys()),
        output_path_dict = output_paths,
        )

    if params.output.write_pymol_script is True:

        write_restraints_pymol_script = WriteRestraintsPymolScript(
            filepath = output_pymol_filepath,
            )

    else:

        write_restraints_pymol_script = None

    ######################################################################
    # Generate restraints
    ######################################################################

    model = AtomicModel.from_file(str(input_pdb))
    model.hierarchy.sort_atoms_in_place()

    restraints_collection = make_restraints(
        hierarchy = model.hierarchy,
        )

    logger.heading('Output Restraints')
    logger(str(restraints_collection))

    formatted_restraints = write_restraints(
        restraints_collection = restraints_collection,
        )

    if write_restraints_pymol_script is not None:

        write_restraints_pymol_script(
            hierarchy = model.hierarchy,
            restraints_collection = restraints_collection,
            )

    logger.heading('make_restraints finished normally')

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

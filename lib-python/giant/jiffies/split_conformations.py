import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys, copy
import pathlib as pl

from giant.phil import (
    log_running_parameters,
    )

from giant.mulch.dataset import (
    AtomicModel,
    )

from giant.structure.conformers import (
    SplitHierarchyByConformer,
    SplitHierarchyByConformerGroup,
    SplitHierarchyByResidueNames,
    SplitMultiStateModel,
    )

############################################################################

PROGRAM = 'giant.split_conformations'

DESCRIPTION = """
    A tool to split a multi-conformer model into its component states.
        Conformation groups to keep can be kept can be declared explicity (mode=by_conformer_group [ + options ])
        or implicitly, through keeping the conformations associated with certain residue names (mode=by_residue_name [ + options ]).

    1) Simple usage:
        > giant.split_conformations input1.pdb input2.pdb [options]
"""

############################################################################

blank_arg_prepend = {
    '.pdb':'input.pdb=',
}

import libtbx.phil
master_phil = libtbx.phil.parse("""
input  {
    pdb = None
        .help = 'A model containing multiple states/conformations (for example unbound and bound states)'
        .type = str
        .multiple = True
}
output {
    suffix_prefix = '-split'
        .help = "prefix for the suffix to be appended to the output structure -- appended to the basename of the input structure before any conformer-specific labels"
        .type = str
    log = None
        .help = "output log file"
        .type = path
}
options {
    mode = *by_residue_name by_conformer by_conformer_group
        .help = 'How to split the model:\n\tby_conformer : output a model for each conformer\n\tby_conformer_group : output a model for a selection of conformers\n\tby_residue_name : output a model containing all states containing a specifice residue name'
        .type = choice
    by_conformer_group {
        conformers = None
            .help = "Output these conformers in one model. Multiple can be provided. e.g. [conformers=ABC conformers=BCD]"
            .multiple = True
            .type = str
    }
    by_residue_name {
        ignore_common_molecules = True
            .type = bool
        combine_bound_states = False
            .type = bool
            .help = "Put all bound states into hierarchy or split by molecule"
        include_resname = None
            .help = "Split conformers based on residues with these names."
            .type = str
            .multiple = True
        ignore_resname = None
            .help = "Ignore these residue names when making decisions."
            .type = str
            .multiple = True
        selected_name = 'bound-state'
            .help = "Output filename component containing selected residues"
            .type = str
        unselected_name = 'ground-state'
            .help = "Output filename component containing residues
            not selected through 'rename'"
            .type = str
        atom_selection = None
            .help = "Select a set of atoms to use to decide (used in combination with other flags)."
            .type = str
    }
    pruning {
        prune_duplicates = True
            .help = 'Remove duplicated conformers in the output structure (convert to blank altloc)'
            .type = bool
        prune_duplicates_rmsd = 0.05
            .type = float
    }
    reset_altlocs = True
        .help = 'Relabel conformers of kept residues to begin with "A" (i.e. C,D,E -> A,B,C)'
        .type = bool
    rescale_occupancies = False
        .help = 'Normalise the occupancies so that the maximum occupancy of the output structure is 1.0 (all relative occupancies are maintained)'
        .type = bool
}
settings {
    overwrite = False
        .type = bool
    verbose = False
        .type = bool
}

""")

def set_get_log(params):

    if params.output.log is not None:
        pass
    elif len(params.input.pdb) == 1:
        params.output.log = str(
            pl.Path(params.input.pdb[0]).with_name(
                pl.Path(params.input.pdb[0]).stem+'-split_conformations.log'
                )
            )
    else:
        params.output.log = "split_conformations.log"

    return params.output.log

def validate_params(params):

    if not params.input.pdb:
        raise IOError('No PDB files given')

    for p in params.input.pdb:
        if not pl.Path(p).exists():
            raise IOError('PDB file {} does not exist'.format(p))

def build_splitter(options):

    if options.mode == 'by_residue_name':

        opt = options.by_residue_name

        split_hierarchy = SplitHierarchyByResidueNames(
            ignore_common_molecules = opt.ignore_common_molecules,
            include_resnames_list = opt.include_resname,
            ignore_resnames_list = opt.ignore_resname,
            selected_name = opt.selected_name,
            unselected_name = opt.unselected_name,
            atom_selection = opt.atom_selection,
            combine_split_states = opt.combine_bound_states,
            )

    elif options.mode == 'by_conformer':

        split_hierarchy = SplitHierarchyByConformer()

    elif options.mode == 'by_conformer_group':

        opt = options.by_conformer_group

        split_hierarchy = SplitHierarchyByConformerGroup(
            conformer_id_sets = opt.conformers,
            )

    else:

        raise NotImplementedError()

    split_states = SplitMultiStateModel(
        split_hierarchy = split_hierarchy,
        prune_duplicates_rmsd = (
            options.pruning.prune_duplicates_rmsd
            if options.pruning.prune_duplicates is True
            else None
            ),
        reset_altlocs = options.reset_altlocs,
        reset_occupancies = options.rescale_occupancies,
        )

    return split_states


class WriteStructures(object):

    def __init__(self, overwrite=False, filepath_append=None):

        self.overwrite = bool(overwrite)

        self.filepath_append = (
            str(filepath_append)
            if filepath_append is not None
            else None
            )

    def __call__(self, model, hierarchy_dict, output_path_root):

        stem = (
            output_path_root.stem
            )

        if self.filepath_append is not None:
            stem = (
                stem + str(self.filepath_append)
                )

        filepaths = {
            k : output_path_root.with_name(
                stem + '-' + k
                ).with_suffix(
                '.pdb'
                )
            for k in list(hierarchy_dict.keys())
            }

        if (self.overwrite is False):

            for k, p in list(filepaths.items()):

                if p.exists():
                    raise IOError(
                        'Output file already exists: {}'.format(
                            str(p)
                            )
                        )

        for k, h in sorted(hierarchy_dict.items()):

            p = filepaths[k]

            logger(
                'Writing {k} to {p}'.format(
                    k = k,
                    p = str(p),
                    )
                )

            h.write_pdb_file(
                file_name = str(p),
                crystal_symmetry = (
                    model.input.crystal_symmetry()
                    ),
                )

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

    #####

    split_states = build_splitter(
        options = params.options,
        )

    write_structures = WriteStructures(
        filepath_append = str(params.output.suffix_prefix),
        overwrite = bool(params.settings.overwrite),
        )

    #####

    input_pdbs = list(params.input.pdb)

    logger.heading('Processing structures')

    for pdb_path in input_pdbs:

        logger.subheading('Splitting {}'.format(str(pdb_path)))

        pdb_path = pl.Path(pdb_path)

        model = AtomicModel.from_file(str(pdb_path))
        model.hierarchy.sort_atoms_in_place()

        split_dict = split_states(
            hierarchy = model.hierarchy,
            )

        write_structures(
            model = model,
            hierarchy_dict = split_dict,
            output_path_root = (
                pdb_path.with_name(pdb_path.stem)
                ),
            )

    logger.heading('split_conformations finished normally')

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

import os, sys, copy

import iotbx.pdb
import libtbx.phil

from bamboo.common.logs import Log

from giant.io.pdb import strip_pdb_to_input, get_pdb_header
from giant.structure.formatting import Labeller
from giant.structure.altlocs import prune_redundant_alternate_conformations

############################################################################

PROGRAM = 'giant.strip_conformations'

DESCRIPTION = """
    A tool to remove unwanted conformations of a multi-conformer model.
        Conformations to keep can be kept can be declared explicity (using conf=...) or implicitly,
        through keeping the conformations associated with certain residue names (res=...).

    1) Simple usage:
        > giant.strip_conformations input.pdb
"""

############################################################################

blank_arg_prepend = {'.pdb':'pdb='}

master_phil = libtbx.phil.parse("""
input  {
    pdb = None
        .help = 'A model containing multiple states/conformations (for example unbound and bound states)'
        .type = str
        .multiple = True
}
output {
    suffix_prefix = 'split'
        .help = "prefix for the suffix to be appended to the output structure -- appended to the basename of the input structure before any conformer-specific labels"
        .type = str
    log = None
        .help = ""
        .type = path
    reset_altlocs = True
        .help = 'Relabel conformers of kept residues to begin with "A" (i.e. C,D,E -> A,B,C)'
        .type = bool
    prune_duplicates = True
        .help = 'Remove duplicated conformers in the output structure (convert to blank altloc)'
        .type = bool
}
options {
    mode = by_conformer by_conformer_group *by_residue_name
        .help = 'How to split the model:\n\tby_conformer : output a model for each conformer\n\tby_conformer_group : output a model for a selection of conformers\n\tby_residue_name : output a model containing all states containing a specifice residue name'
        .type = choice
    by_conformer_group {
        conformers = C,D,E,F,G
            .help = "Output these conformers in one model. Multiple can be provided."
            .multiple = True
            .type = str
    }
    by_residue_name {
        resname = DRG,FRG,LIG,UNK,UNL
            .help = "Group conformers containing any of these residue names"
            .type = str
        selected_name = 'bound-state'
            .help = "Output filename component containing selected residues"
            .type = str
        unselected_name = 'ground-state'
            .help = "Output filename component containing residues not selected through 'rename'"
            .type = str
    }
}
settings {
    overwrite = False
        .type = bool
    verbose = False
        .type = bool
}

""")

############################################################################

def split_conformations(filename, params, log=None):

    if log is None: log = Log(verbose=True)

    # Read the pdb header - for writing later...
    header_contents = get_pdb_header(filename)

    # Read in and validate the input file
    ens_obj = strip_pdb_to_input(filename, remove_ter=True)
    ens_obj.hierarchy.only_model()

    # Create a new copy of the structures
    new_ens = ens_obj.hierarchy.deep_copy()

    # Extract conformers from the structure as set
    all_confs = set(ens_obj.hierarchy.altloc_indices())
    all_confs.discard('')

    if params.options.mode == 'by_residue_name':
        sel_resnames = params.options.by_residue_name.resname.split(',')
        sel_confs = [ag.altloc for ag in new_ens.atom_groups() if (ag.resname in sel_resnames)]
        # List of conformers to output for each structure, and suffixes
        out_confs = map(sorted, [all_confs.intersection(sel_confs), all_confs.difference(sel_confs)])
        out_suffs = [params.options.by_residue_name.selected_name, params.options.by_residue_name.unselected_name]
    elif params.options.mode == 'by_conformer':
        sel_resnames = None
        sel_confs = None
        # One structure for each conformer
        out_confs = [[c] for c in sorted(all_confs)]
        out_suffs = [''.join(c) for c in out_confs]
    elif params.options.mode == 'by_conformer_group':
        sel_resnames = None
        sel_confs = None
        # One structure for each set of supplied conformer sets
        out_confs = [s.split(',') for s in params.options.by_conformer_group.conformers]
        out_suffs = [''.join(c) for c in out_confs]
    else:
        raise Exception('Invalid selection for options.mode: {}'.format(params.options.mode))

    print out_confs
    print out_suffs

    # Create paths from the suffixes
    out_paths = ['.'.join([os.path.splitext(filename)[0],params.output.suffix_prefix,suff,'pdb']) for suff in out_suffs]

    log.subheading('Pruning conformers in {}'.format(filename[-70:]))

    for this_confs, this_path in zip(out_confs, out_paths):

        if not this_confs: continue

        # Select atoms to keep - no altloc, or altloc in selection
        sel_string = ' or '.join(['altid " "']+['altid "{}"'.format(alt) for alt in this_confs])
        # Extract selection from the hierarchy
        sel_hiery = new_ens.select(new_ens.atom_selection_cache().selection(sel_string), copy_atoms=True)

        log.bar(True, False)
        log('Outputting conformer(s) {} to {}'.format(''.join(this_confs), this_path))
        log.bar()
        log('Keeping ANY atom with conformer id: {}'.format(' or '.join(['" "']+this_confs)))
        log('Selection: \n\t'+sel_string)
        log.bar()

        if params.output.prune_duplicates:
            # Remove an alternate conformers than are duplicated after selection
            prune_redundant_alternate_conformations(
                hierarchy           = sel_hiery,
                required_altlocs    = sel_hiery.altloc_indices(),
                rmsd_cutoff         = 0.1,
                in_place            = True,
                verbose             = params.settings.verbose)

        if params.output.reset_altlocs:
            # Change the altlocs so that they start from "A"
            if len(this_confs) == 1:
                conf_hash = {this_confs[0]: ' '}
            else:
                conf_hash = dict(zip(this_confs, iotbx.pdb.systematic_chain_ids()))
            log('Resetting structure altlocs:')
            for k in sorted(conf_hash.keys()):
                log('\t{} -> "{}"'.format(k, conf_hash[k]))
            log.bar()
            for ag in sel_hiery.atom_groups():
                if ag.altloc in this_confs:
                    if params.settings.verbose:
                        log('{} -> alt {}'.format(Labeller.format(ag), conf_hash[ag.altloc]))
                    ag.altloc = conf_hash[ag.altloc]
            if params.settings.verbose: log.bar()

        log('Writing structure: {}'.format(this_path))
        log.bar(False, True)

        # Write header contents
        with open(this_path, 'w') as fh: fh.write(header_contents)
        # Write output file
        sel_hiery.write_pdb_file(this_path, open_append=True)

    return out_paths

############################################################################

def run(params):

    # Create log file
    log = Log(log_file=params.output.log, verbose=True)

    log.heading('Validating input parameters')

    assert params.input.pdb, 'No PDB files given'

    print params.options.by_conformer_group.conformers

    log.heading('Splitting multi-state structures')

    # Iterate through the input structures and extract the conformation
    for pdb in params.input.pdb:
        split_conformations(filename=pdb, params=params, log=log)

    log.heading('FINISHED')

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
        program             = PROGRAM,
        description         = DESCRIPTION)

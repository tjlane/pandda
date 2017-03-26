import os, sys, copy

import iotbx.pdb
import libtbx.phil

from giant.io.pdb import strip_pdb_to_input
from giant.structure.utils import resolve_residue_id_clashes, transfer_residue_groups_from_other
from giant.structure.altlocs import expand_alternate_conformations, find_next_conformer_idx, increment_altlocs, prune_redundant_alternate_conformations
from giant.structure.occupancy import set_conformer_occupancy

############################################################################

PROGRAM = 'giant.merge_conformations'
DESCRIPTION = """
    A tool to merge two models of the same crystal into one structure.
        The alternate conformers of the two models are generated and merged. Any common sections are then
        converted back into a single conformer model. This tool is therefore particularly useful when one
        model is an edited version of the other as many of the residues will be identical and the output
        structure will not contain too many alternate conformations.

    1) Simple usage:
        > giant.merge_conformations major=structure_1.pdb minor=structure_2.pdb
"""

############################################################################

blank_arg_prepend = None

master_phil = libtbx.phil.parse("""
input {
    major = None
        .help = 'The major conformation of the protein (the conformation you want to have the lowest altlocs)'
        .type = str
    minor = None
        .help = 'The minor conformation of the protein (the conformation that will have the higher altlocs)'
        .type = str
}
output {
    pdb = 'merged.pdb'
        .help = 'output pdb file'
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

def merge_complementary_hierarchies(hierarchy_1, hierarchy_2, verbose=False):
    """Merge hierarchies that are alternate models of the same crystal by expanding alternate model conformations, merging, and then trimming alternate conformations where possible"""

    # Copy the hierarchies
    hierarchy_1 = hierarchy_1.deep_copy()
    hierarchy_2 = hierarchy_2.deep_copy()
    # Sort the atoms
    hierarchy_1.sort_atoms_in_place()
    hierarchy_2.sort_atoms_in_place()

    print '============================================ *** ============================================'
    print '                                 Preparing to merge structures'
    print '============================================ *** ============================================'
    print ''
    print '-------------------------------------------- *** --------------------------------------------'
    print 'Explicitly expanding models to all conformations of the crystal'
    print '-------------------------------------------- *** --------------------------------------------'
    print ''
    print 'Expanding alternate conformations in structure 1'
    expand_alternate_conformations(
        hierarchy   = hierarchy_1,
        in_place    = True,
        verbose     = verbose)
    print 'Expanding alternate conformations in structure 2'
    expand_alternate_conformations(
        hierarchy   = hierarchy_2,
        in_place    = True,
        verbose     = verbose)
    print ''
    print '-------------------------------------------- *** --------------------------------------------'
    print 'Applying conformer shift to the second structure before merging'
    print '-------------------------------------------- *** --------------------------------------------'
    print ''
    print 'Identifying the altloc shift required from the number of alternate conformers in structure 1'
    conf_offset = find_next_conformer_idx(
        hierarchy           = hierarchy_1,
        all_ids             = iotbx.pdb.systematic_chain_ids())
    print 'Incrementing all altlocs in structure 2 by', conf_offset
    increment_altlocs(
        hierarchy           = hierarchy_2,
        offset              = conf_offset,
        in_place            = True,
        verbose             = verbose)
    print ''
    print '-------------------------------------------- *** --------------------------------------------'
    print 'Renaming residues that do not align between structures'
    print '-------------------------------------------- *** --------------------------------------------'
    resolve_residue_id_clashes(
        fixed_hierarchy     = hierarchy_1,
        moving_hierarchy    = hierarchy_2,
        in_place            = True,
        verbose             = verbose)
    print ''
    print '============================================ *** ============================================'
    print '                                      Merging structures'
    print '============================================ *** ============================================'
    print ''
    print 'Transferring residues from Structure 2 to Structure 1'
    transfer_residue_groups_from_other(
        acceptor_hierarchy  = hierarchy_1,
        donor_hierarchy     = hierarchy_2,
        in_place            = True,
        verbose             = verbose)
    print ''
    print '============================================ *** ============================================'
    print '                                 Post-processing structure'
    print '============================================ *** ============================================'
    print ''
    print 'Pruning unneccessary multi-conformer residues in the merged structure'
    prune_redundant_alternate_conformations(
        hierarchy           = hierarchy_1,
        required_altlocs    = hierarchy_1.altloc_indices(),
        rmsd_cutoff         = 0.1,
        in_place            = True,
        verbose             = verbose)
    print ''
    altlocs = [a for a in hierarchy_1.altloc_indices() if a]
    new_occ = 1.0/len(altlocs)
    print 'Setting all conformer ({}) occupancies to {}'.format(','.join(altlocs), new_occ)
    hierarchy_1 = set_conformer_occupancy(
        hierarchy   = hierarchy_1,
        altlocs     = altlocs,
        occupancy   = new_occ,
        in_place    = True,
        verbose     = verbose)

    return hierarchy_1

############################################################################

def run(params):

    print ''
    print '============================================ *** ============================================'
    print '                         Validating input parameters and input files'
    print '============================================ *** ============================================'
    print ''

    assert params.input.major
    assert params.input.minor
    assert params.output.pdb
    if os.path.exists(params.output.pdb):
        if params.settings.overwrite: os.remove(params.output.pdb)
        else: raise Exception('File already exists: {}'.format(params.output.pdb))

    # Read in the ligand file and set each residue to the requested conformer
    maj_obj = strip_pdb_to_input(params.major, remove_ter=True)
    min_obj = strip_pdb_to_input(params.minor, remove_ter=True)

    # Check that ... something
    try:
        maj_obj.hierarchy.only_model()
        min_obj.hierarchy.only_model()
    except:
        raise Sorry('Structures may only have one model')

    final_struct = merge_complementary_hierarchies(
        hierarchy_1 = maj_obj.hierarchy,
        hierarchy_2 = min_obj.hierarchy,
        verbose     = params.settings.verbose)

    print 'Writing output structure'

    # Update the atoms numbering
    final_struct.sort_atoms_in_place()
    final_struct.atoms_reset_serial()
    # Write output file
    final_struct.write_pdb_file(file_name=params.output.pdb, crystal_symmetry=maj_obj.crystal_symmetry())

    print ''
    print '...finished!'
    print ''

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
        description         = DESCRIPTION)

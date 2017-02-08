import os, sys, copy

import numpy

import iotbx.pdb

import libtbx.phil

from giant.structure.restraints.peptides import generate_set_of_alternate_conformer_peptide_links, format_link_record
from giant.structure.restraints.conformers import find_duplicated_conformers_and_generate_atom_pairs
from giant.structure.restraints.external import find_atoms_around_alternate_conformers
from giant.structure.restraints.occupancy import overlapping_occupancy_groups, simple_occupancy_groups
from giant.structure.formatting import RefmacFormatter, PhenixFormatter

#######################################

PROGRAM = 'giant.make_restraints'
DESCRIPTION = """
    A tool to simplify the generation of restraints for use during refinement with REFMAC or PHENIX.

    The output ".params" files may be passed to giant.quick_refine for use in refinement.

    1) Simple usage:
        > giant.make_restraints input.pdb

    2) With all modes active
        > giant.make_restraints input.pdb all=True
"""

blank_arg_prepend = {'.pdb' : 'pdb='}

master_phil = libtbx.phil.parse("""
input {
    pdb = None
       .help = 'Protein model'
        .type = str
}

modes {
    all = False
        .help = "Turn on all funcationality"
        .multiple = False
        .type = bool
    peptide_bond_links = True
        .help = "Make continuity records for alternate conformations of the peptide backbone to ensure no discontinuities"
        .type = bool
    duplicated_atom_restraints = True
        .help = "Generate parameter restraints for duplicated residues (reduce free parameters in refinement)"
        .type = bool
    local_structure_restraints = False
        .help = "Generate a set of local distance restraints around certain regions to restrain them in refinement"
        .type = bool
    occupancy_groups = True
        .help = "Make occupancy constraints (and restraints) for groups in refinement"
        .type = bool
    b_factor_restraints = False
        .help = "NOT CURRENTLY IMPLEMENTED - Make B-factor restraints for ligands to surrounding residues"
        .expert_level = 1
        .type = bool
}

output {
    phenix = 'restraints_phenix.params'
        .help = 'Output restraints file for phenix'
        .type = path
    refmac = 'restraints_refmac.params'
        .help = 'Output restraints file for refmac'
        .type = path
}

duplicates {
    make_for = *all protein het
        .help = "Generate only restraints for protein residues?"
        .type = choice
        .multiple = False
    rmsd_cutoff = 0.1
        .help = "Cutoff at which two conformers are considered to be duplicated"
        .type = float
    sigma_xyz = 0.02
        .help = "Coordinate restraint term controlling how strongly the restraint is enforced"
        .type = float
}

peptide_bonds {
    suffix = '.link.pdb'
        .type = str
}

local_restraints {
    altlocs = C,D,E,F,G
        .help = "Which altlocs to generate local restraints for"
        .type = str
    max_distance = 4.2
        .help = "Maximum distance to create local restraints between atoms"
        .type = float
    min_distance = 1.6
        .help = "Minimum distance to create local restraints between atoms"
        .type = float
    sigma_xyz = 0.1
        .help = "Sigma of the distance restraint controlling how strongly the restraint is enforced"
        .type = float
}

occupancy {
    resname = DRG,FRG,LIG,UNK,UNL
        .help = 'Residues to generate constraint groups around for occupancy refinement (comma separated list of residue identifiers, i.e. resname=LIG or resname=LIG,UNL)'
        .type = str
    group_dist = 5
        .type = float
        .help = 'Distance to use when clustering atoms that should have the SAME occupancy'
    overlap_dist = 2
        .type = float
        .help = 'Distance to use when clustering atoms that should have occupancies that SUM TO LESS THAN ONE'
    complete_groups = True
        .help = 'Generate a set of fully constrained groups (that sum to unitary occupancy) when True. Generate a set of weaker constraints for overlapping atoms when False.'
        .type = bool
    simple_groups = False
        .help = 'Generate the set of default occupancy groups (conformers of the same residue, connected sets of residues)'
        .type = bool
    exclude_altlocs = None
        .help = 'Exclude certain altlocs from occupancy groups (e.g. A or A,B)'
        .type = str
}

b_factors
    .expert_level = 1
{
    resname = None
        .help = 'Residues to create B-factor restraints for. If None, defaults to occupancy.resname'
        .type = str
    contact_dist = 4
        .help = 'Cutoff distance for identifying contacting atoms'
        .type = float
    sigma_b = 1.0
        .help = "Sigma of the b-factor restraint controlling how strongly the restraint is enforced"
        .type = float
}

overwrite = True
    .type = bool
verbose = True
    .type = bool

""")

#######################################

def make_link_records(params, input_hierarchy, link_file):
    """Create link records to make a continuous peptide chain"""

    print '\n\n\n'
    print '============================================ *** ============================================'
    print '{:^93}'.format('Checking the continuity of the protein backbone')
    print '============================================ *** ============================================'
    print ''

    links, warnings = generate_set_of_alternate_conformer_peptide_links(hierarchy=input_hierarchy.hierarchy)

    if warnings:
        print '============================================>'
        print 'WARNINGS:'
        print '============================================>'
        print ''
        for w in warnings: print w
        print ''

    if (not links) and (not warnings):
        print 'No breaks in the backbone - hooray! (nothing needs to be done here)'
        return
    elif (not links):
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "!!! >>> There are breaks in the backbone but I'm not able to do anything to fix them    <<< !!!"
        print "!!! >>> You'll need to check them manually to see if these are going to be a problem... <<< !!!"
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        return

    link_block = '\n'.join([format_link_record(atom_1=a1,atom_2=a2,chain_id_1=c1,chain_id_2=c2,link_type=lt) for a1,a2,c1,c2,lt in links])

    print 'Need to apply {} links to make the backbone continuous:'.format(len(links))
    print ''
    print link_block
    print ''

    print 'Writing hierarchy with new link records to {}'.format(link_file)
    print '(This file can only be used for refinement with REFMAC)'
    print ''
    print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print '!!! ALTHOUGH THE FILE WITH BACKBONE LINKS HAS BEEN OUTPUT, IT SHOULD BE USED WITH CAUTION !!!'
    print '!!!   THE CONNECTION OF ALTERNATE CONFORMATIONS OF THE BACKBONE IS GENERALLY "INCORRECT"  !!!'
    print '!!!          THERE SHOULD BE A VERY GOOD REASON FOR THESE RESTRAINTS TO BE USED           !!!'
    print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

    input_hierarchy.hierarchy.write_pdb_file(file_name          = link_file,
                                             crystal_symmetry   = input_hierarchy.crystal_symmetry(),
                                             link_records       = link_block)

def make_duplication_restraints(params, input_hierarchy):
    """Create coordinate and b-factor restraints for duplicated conformers"""

    print '\n\n\n'
    print '============================================ *** ============================================'
    print '{:^93}'.format('Generating restraints for duplicated conformers')
    print '============================================ *** ============================================'
    print ''

    dup_groups = []

    for chn in input_hierarchy.hierarchy.chains():

        if (params.duplicates.make_for == 'protein') and not chn.is_protein():
            continue
        elif (params.duplicates.make_for == 'het') and chn.is_protein():
            continue

        for rg in chn.residue_groups():
            dup_groups += find_duplicated_conformers_and_generate_atom_pairs(residue_group=rg, rmsd_cutoff=params.duplicates.rmsd_cutoff)

    if not dup_groups:
        print 'No duplicated conformers (no restraints created)'
        return

    # Concatenate atoms into one list
    atom_pairs = []; [atom_pairs.extend(l) for l in dup_groups]

    print 'Found {} duplicated conformers consisting of {} atoms'.format(len(dup_groups), len(atom_pairs))
    print ''

    if params.output.refmac:
        restraint_list = [RefmacFormatter.make_distance_restraint(atm_1=a1, atm_2=a2, value=0.0, sigma=params.duplicates.sigma_xyz) for a1,a2 in atom_pairs]
        rest_block = RefmacFormatter.format_distance_restraints(restraint_list=restraint_list)
        with open(params.output.refmac, 'a') as fh: fh.write(rest_block+'\n')
        if params.verbose:
            print '-------------------------------------------- *** --------------------------------------------'
            print '{:^93}'.format('refmac duplicate conformer restraints (output to {})'.format(params.output.refmac))
            print '-------------------------------------------- *** --------------------------------------------'
            print ''
            print rest_block[:1000]+'...'*(len(rest_block)>1000)
            print ''

    if params.output.phenix:
        restraint_list = [PhenixFormatter.make_distance_restraint(atm_1=a1, atm_2=a2, value=0.0, sigma=params.duplicates.sigma_xyz) for a1,a2 in atom_pairs]
        rest_block = PhenixFormatter.format_distance_restraints(restraint_list=restraint_list)
        with open(params.output.phenix, 'a') as fh: fh.write(rest_block+'\n')
        if params.verbose:
            print '-------------------------------------------- *** --------------------------------------------'
            print '{:^93}'.format('phenix duplicate conformer restraints (output to {})'.format(params.output.phenix))
            print '-------------------------------------------- *** --------------------------------------------'
            print ''
            print rest_block[:1000]+'...'*(len(rest_block)>1000)
            print ''

def make_local_restraints(params, input_hierarchy):
    """Create local restraints for a hierarchy"""

    print '\n\n\n'
    print '============================================ *** ============================================'
    print '{:^93}'.format('Generating local structure restraints')
    print '============================================ *** ============================================'
    print ''

    atom_d_pairs = find_atoms_around_alternate_conformers(hierarchy     = input_hierarchy.hierarchy,
                                                          altlocs       = params.local_restraints.altlocs.split(',') if params.local_restraints.altlocs else None,
                                                          dist_cutoff   = params.local_restraints.max_distance)
    # Filter the 0-distance restraints
    atom_d_pairs = [(a1,a2,d) for a1,a2,d in atom_d_pairs if d>params.local_restraints.min_distance]

    print 'Created {} local restraints for {} conformers with distance cutoff of {}-{}A'.format(len(atom_d_pairs),
                                                                                               params.local_restraints.altlocs if params.local_restraints.altlocs else 'all',
                                                                                               params.local_restraints.min_distance,
                                                                                               params.local_restraints.max_distance)
    print ''

    if params.output.refmac:
        restraint_list = [RefmacFormatter.make_distance_restraint(atm_1=a1, atm_2=a2, value=d, sigma=params.local_restraints.sigma_xyz) for a1,a2,d in atom_d_pairs]
        rest_block = RefmacFormatter.format_distance_restraints(restraint_list=restraint_list)
        with open(params.output.refmac, 'a') as fh: fh.write(rest_block+'\n')
        if params.verbose:
            print '-------------------------------------------- *** --------------------------------------------'
            print '{:^93}'.format('refmac local structural restraints (output to {})'.format(params.output.refmac))
            print '-------------------------------------------- *** --------------------------------------------'
            print ''
            print rest_block[:1000]+'...'*(len(rest_block)>1000)
            print ''

    if params.output.phenix:
        restraint_list = [PhenixFormatter.make_distance_restraint(atm_1=a1, atm_2=a2, value=d, sigma=params.local_restraints.sigma_xyz) for a1,a2,d in atom_d_pairs]
        rest_block = PhenixFormatter.format_distance_restraints(restraint_list=restraint_list)
        with open(params.output.phenix, 'a') as fh: fh.write(rest_block+'\n')
        if params.verbose:
            print '-------------------------------------------- *** --------------------------------------------'
            print '{:^93}'.format('phenix duplicate conformer restraints (output to {})'.format(params.output.phenix))
            print '-------------------------------------------- *** --------------------------------------------'
            print ''
            print rest_block[:1000]+'...'*(len(rest_block)>1000)
            print ''

def make_occupancy_constraints(params, input_hierarchy):
    """Create occupancy groups for a hierarchy"""

    print '\n\n\n'
    print '============================================ *** ============================================'
    print '{:^93}'.format('Generating occupancy-constrained groups')
    print '============================================ *** ============================================'
    print ''

    # Ligand resname identifiers
    resnames = params.occupancy.resname.split(',')
    if params.verbose:
        print 'Looking for ligands with resname {!s}'.format(' or '.join(resnames))
        print ''

    # Make occupancy groups
    occupancy_groups = overlapping_occupancy_groups(hierarchy       = input_hierarchy.hierarchy,
                                                    resnames        = resnames,
                                                    group_dist      = params.occupancy.group_dist,
                                                    overlap_dist    = params.occupancy.overlap_dist,
                                                    complete_groups = params.occupancy.complete_groups,
                                                    exclude_altlocs = params.occupancy.exclude_altlocs.split(',') if params.occupancy.exclude_altlocs else [],
                                                    verbose         = params.verbose)
    # Record whether the occupancy groups are complete (occs sum to 1)
    if params.occupancy.complete_groups:
        occupancy_complete = [True]*len(occupancy_groups)
    else:
        occupancy_complete = [False]*len(occupancy_groups)

    if not occupancy_groups:
        print 'No matching residues were found (no occupancy constraints created)'
        return

    print '-------------------------------------->'
    print ''
    print 'Created {} occupancy groups for overlapping conformers'.format(len(occupancy_groups))
    print ''

    # Ref-make the default occupancy groups?
    if params.occupancy.simple_groups:
        print 'simple_groups=={}: Remaking default occupancy restraints for residues'.format(params.occupancy.simple_groups)
        if params.verbose: print ''
        simple_groups = simple_occupancy_groups(hierarchy   = input_hierarchy.hierarchy,
                                                verbose     = params.verbose)
        num_alts = len([a for a in input_hierarchy.hierarchy.altloc_indices() if a!=''])
        occupancy_complete += [True if len(g)==num_alts else False for g in simple_groups]
        occupancy_groups += simple_groups
        if params.verbose: print ''
        print 'Increased number of occupancy groups to {}'.format(len(occupancy_groups))
        print ''

    if params.output.refmac:
        restraint_list = RefmacFormatter.make_occupancy_restraints(list_of_lists_of_groups  = occupancy_groups,
                                                                   group_completeness       = occupancy_complete)
        rest_block = RefmacFormatter.format_occupancy_restraints(restraint_list=restraint_list)
        with open(params.output.refmac, 'a') as fh: fh.write(rest_block+'\n')
        if params.verbose:
            print '-------------------------------------------- *** --------------------------------------------'
            print '{:^93}'.format('refmac occupancy restraints (output to {})'.format(params.output.refmac))
            print '-------------------------------------------- *** --------------------------------------------'
            print ''
            print rest_block[:1000]+'...'*(len(rest_block)>1000)
            print ''

    if params.output.phenix:
        restraint_list = PhenixFormatter.make_occupancy_restraints(list_of_lists_of_groups  = occupancy_groups,
                                                                   group_completeness       = occupancy_complete)
        rest_block = PhenixFormatter.format_occupancy_restraints(restraint_list=restraint_list)
        with open(params.output.phenix, 'a') as fh: fh.write(rest_block+'\n')
        if params.verbose:
            print '-------------------------------------------- *** --------------------------------------------'
            print '{:^93}'.format('phenix occupancy restraints (output to {})'.format(params.output.phenix))
            print '-------------------------------------------- *** --------------------------------------------'
            print ''
            print rest_block[:1000]+'...'*(len(rest_block)>1000)
            print ''

def make_b_factor_restraints(params, input_hierarchy):
    pass

#######################################

def run(params):

    ######################################################################
    # Validate input
    ######################################################################

    assert params.input.pdb, 'No PDB File Provided'

    if params.modes.all:
        params.modes.peptide_bond_links         = True
        params.modes.duplicated_atom_restraints = True
        params.modes.local_structure_restraints = True
        params.modes.occupancy_groups           = True
        params.modes.b_factor_restraints        = True

    if params.modes.peptide_bond_links:
        link_file = os.path.splitext(params.input.pdb)[0]+params.peptide_bonds.suffix
    if params.modes.duplicated_atom_restraints:
        pass
    if params.modes.local_structure_restraints:
        pass
    if params.modes.occupancy_groups:
        pass
    if params.modes.b_factor_restraints:
        pass

    ######################################################################
    # Prepare output and input
    ######################################################################
    if params.output.phenix and os.path.exists(params.output.phenix):
        if params.overwrite: os.remove(params.output.phenix)
        else: raise Exception('File already exists: {}'.format(params.output.phenix))
    if params.output.refmac and os.path.exists(params.output.refmac):
        if params.overwrite: os.remove(params.output.refmac)
        else: raise Exception('File already exists: {}'.format(params.output.refmac))

    pdb_obj = iotbx.pdb.hierarchy.input(params.input.pdb)
    pdb_obj.hierarchy.sort_atoms_in_place()

    ######################################################################
    # Generate restraints
    ######################################################################

    if params.modes.peptide_bond_links:
        make_link_records(params=params, input_hierarchy=pdb_obj, link_file=link_file)

    if params.modes.duplicated_atom_restraints:
        make_duplication_restraints(params=params, input_hierarchy=pdb_obj)

    if params.modes.local_structure_restraints:
        make_local_restraints(params=params, input_hierarchy=pdb_obj)

    if params.modes.occupancy_groups:
        make_occupancy_constraints(params=params, input_hierarchy=pdb_obj)

    if params.modes.b_factor_restraints:
        make_b_factor_restraints(params=params, input_hierarchy=pdb_obj)

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

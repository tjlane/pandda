import os, sys, copy

import libtbx.phil

import numpy

import iotbx.pdb

from scitbx.array_family import flex

from giant.utils.pdb import strip_pdb_to_input, get_pdb_header

#######################################

blank_arg_prepend = {'.pdb':'pdb='}

master_phil = libtbx.phil.parse("""
pdb = None
    .help = 'The ensemble of protein conformations (normally the unbound and bound structures)'
    .type = str
    .multiple = True

res = LIG,UNL
    .help = 'Residues to be define selected conformations (comma separated list of residue names, i.e. res=LIG or res=LIG,UNL)'
    .type = str
conf = None
    .help = 'Define selected conformations explicitly (comma separated list of conformer IDs, i.e. conf=A or conf=A,B)'
    .type = str

suffix = '.stripped.pdb'
    .help = 'output suffix for input files'
    .type = str

reset_altlocs = True
    .help = 'Relabel conformers of kept residues'
    .type = bool

overwrite = False
    .type = bool
verbose = False
    .type = bool

log = None
    .type = path
    .multiple = False
""")

#######################################

def proc(ensemble_file, params, sel_resnames=None, sel_confs=None):

    # Create output file path
    output_file = os.path.splitext(ensemble_file)[0] + params.suffix

    # Read the pdb header
    header_contents = get_pdb_header(ensemble_file)

    # Read in the ligand file and set each residue to the requested conformer
    ens_obj = strip_pdb_to_input(ensemble_file, remove_ter=True)

    # Check that ... something
    ens_obj.hierarchy.only_model()

    # Create a new copy of the structures
    new_ens = ens_obj.hierarchy.deep_copy()

    confs_to_select = set()
    if params.res:
        confs_to_select = confs_to_select.union([ag.altloc for ag in new_ens.atom_groups() if (ag.resname in sel_resnames) and ag.altloc])
    if params.conf:
        confs_to_select = confs_to_select.union(sel_confs)
    confs_to_select = sorted(confs_to_select)

    print 'Keeping Conformer IDs: {}'.format(', '.join(confs_to_select))

    if not confs_to_select: raise Exception('No Conformers Selected')

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'PRUNING OUTPUT STRUCTURE'
    ######################################################################

    print '===========================================>>>'
    print 'Keeping ANY atom with conformer id: {}'.format(' or '.join(['" "']+confs_to_select))
    print 'Keeping ALL amino acid conformers with no conformer in {}'.format(' or '.join(confs_to_select))
    print 'Removing amino acid conformers where conformer not in {} but there is a conformer in this'.format(' or '.join(confs_to_select))
    print 'Removing any other atom where conformer id not in {}'.format(' or '.join(confs_to_select))
    print '===========================================>>>'

    # Select atoms to keep - no altloc, or altloc in selection
    sel_string_1 = ' or '.join(['altid " "']+['altid {}'.format(alt) for alt in confs_to_select])
    print 'Selection 1: \n\t'+sel_string_1

    # Select amino acids with no conformer overlapping with confs_to_select
    sel_ags = [ag for ag in new_ens.atom_groups() if (iotbx.pdb.common_residue_names_get_class(ag.resname)=='common_amino_acid') and not set(confs_to_select).intersection([a.altloc for a in ag.parent().atom_groups()])]
    sel_string_2 = ' or '.join(['(chain {!s} and resid {!s} and altid "{!s}")'.format(ag.parent().parent().id, ag.parent().resid(), ag.altloc) for ag in sel_ags])
    print 'Selection 2: \n\t'+sel_string_2.replace('or ','or\n\t')[:500]+'...'

    sel_cache = new_ens.atom_selection_cache()
    sel_string = ' or '.join(['({})'.format(s) for s in [sel_string_1, sel_string_2]])

    conf_sel = new_ens.select(sel_cache.selection(sel_string))

    if params.log:
        fh = open(params.log, 'w')
        fh.write('>>> Residues where the conformer will be reset\n')
    else:
        fh = None

    if params.reset_altlocs:
        ######################################################################
        if params.verbose: print '===========================================>>>'
        if params.verbose: print 'Resetting OUTPUT STRUCTURE'
        ######################################################################

        if len(confs_to_select)==1: hash = {confs_to_select[0]: ' '}
        else:                       hash = dict(zip(confs_to_select, iotbx.pdb.systematic_chain_ids()))

        if params.verbose:
            for it_pr in hash.items(): print 'Changing ALTLOCs:', ' -> '.join(it_pr)

        for ag in conf_sel.atom_groups():
            if ag.altloc in confs_to_select:
                # Reset the altloc
                ag.altloc = hash[ag.altloc]
                # Log this atom group
                if fh: fh.write('> Chain {}, Residue {}, (New) Conformer {}\n'.format(ag.parent().parent().id, ag.parent().resid(), ag.altloc))

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'WRITING OUTPUT STRUCTURE'
    ######################################################################

    # Update the atoms numbering
#        conf_sel.atoms_reset_serial()
    with open(output_file, 'w') as fh:
        fh.write(header_contents)
    # Write output file
    conf_sel.write_pdb_file(output_file, open_append=True)

#######################################

def run(params):

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'VALIDATING INPUT PARAMETERS'
    ######################################################################

    assert params.pdb, 'No PDB files given'
    assert params.res or params.conf, 'Must give either res or conf'

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'SPLITTING ENSEMBLES'
    ######################################################################

    if params.res:
        sel_resnames = params.res.split(',')
    else:
        sel_resnames = None

    if params.conf:
        sel_confs = params.conf.split(',')
    else:
        sel_confs = None

    # Iterate through the input structures and extract the conformation
    for ensemble_file in params.pdb:
        try:
            proc(ensemble_file, params, sel_resnames=sel_resnames, sel_confs=sel_confs)
        except:
            print "Failed:", ensemble_file
            raise

    return

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)

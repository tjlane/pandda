import os, sys, glob

import libtbx.phil

from bamboo.common.command import CommandManager
from bamboo.common.path import rel_symlink

#######################################

blank_arg_prepend = {   '.mtz':'mtz=',
                        '.pdb':'pdb=',
                        '.cif':'cif=',
                        '.params':'params=' }

master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .type = str
    mtz = 'refine.mtz'
        .type = str
    cif = None
        .type = str
        .multiple = True
    params = None
        .type = str
}
settings {
    program = *phenix refmac
        .type = choice
    create_occupancy_params = False
        .type = bool
    args = None
        .help = "Pass any additional arguments to the program? (Command Line arguments). Pass as quoted strings."
        .type = str
        .multiple = True
}
phenix {
    ordered_solvent = False
        .type = bool
    anisotropic = False
        .type = bool
    no_coords = False
        .type = bool
}
output {
    dir_prefix = 'Refine_'
        .type = str
    out_prefix = 'refine'
        .type = str
    link_prefix = 'refine'
        .type = str
}

""")

#######################################

def run(params):

    assert params.input.pdb is not None, 'No PDB given for refinement'
    assert params.input.mtz is not None, 'No MTZ given for refinement'
    if params.settings.create_occupancy_params: assert params.input.params, 'settings.create_occupancy_params is True, but you must also define input.params.'

    ########################
    current_dirs = sorted(glob.glob(params.output.dir_prefix+'*'))
    if not current_dirs:
        out_dir = params.output.dir_prefix+'1'
    else:
        current_dirs = [s.replace(params.output.dir_prefix, '') for s in current_dirs]
        out_dir = params.output.dir_prefix + str(sorted(map(int, current_dirs))[-1]+1)
    print 'Outputting to {}'.format(out_dir)
    os.mkdir(out_dir)

    ########################
    # Link input
    rel_symlink(params.input.pdb, os.path.abspath(os.path.join(out_dir, 'input.pdb')))
    rel_symlink(params.input.mtz, os.path.abspath(os.path.join(out_dir, 'input.mtz')))

    ########################
    output_prefix = os.path.join(out_dir, params.output.out_prefix)
    print 'Output Real Prefix = {}'.format(output_prefix)
    print 'Output Link Prefix = {}'.format(params.output.link_prefix)

    ########################
    # Create occupancy parameters
    if params.settings.create_occupancy_params and params.input.params:
        print 'Running giant.create_occupancy_params'
        print 'Generating occupancy model for the input pdb: {}'.format(params.input.pdb)
        if params.settings.program == 'phenix':
            os.system('giant.create_occupancy_params pdb={} phenix_occ_out={} refmac_occ_out=None verbose=False'.format(params.input.pdb, params.input.params))
        elif params.settings.program == 'refmac':
            os.system('giant.create_occupancy_params pdb={} phenix_occ_out=None refmac_occ_out={} verbose=False'.format(params.input.pdb, params.input.params))

    ########################
    # Phenix
    if params.settings.program == 'phenix':
        ########################
        # Build command string
        cm = CommandManager('phenix.refine')
        cm.add_command_line_arguments([ params.input.pdb, params.input.mtz ])
        cm.add_command_line_arguments([ 'output.prefix={}'.format(output_prefix) ])
        if params.input.cif:
            cm.add_command_line_arguments( params.input.cif )
        if params.input.params and os.path.exists(params.input.params):
            cm.add_command_line_arguments([ params.input.params ])
        if params.phenix.ordered_solvent:
            cm.add_command_line_arguments([ 'ordered_solvent=True' ])
        if params.phenix.anisotropic:
            cm.add_command_line_arguments([ 'adp.individual.anisotropic="not element H"' ])
        if params.phenix.no_coords:
            cm.add_command_line_arguments([ 'strategy=occupancies+individual_adp' ])

    elif params.settings.program == 'refmac':
        ########################
        # Build command string
        cm = CommandManager('refmac5')
        cm.add_command_line_arguments( ['xyzin', params.input.pdb, 'hklin', params.input.mtz] )
        cm.add_command_line_arguments( ['xyzout', output_prefix+'.pdb', 'hklout', output_prefix+'.mtz'] )
        if params.input.cif:
            for cif in params.input.cif:
              cm.add_command_line_arguments( ['libin', cif] )

        # Read in occupancy params
        if params.input.params:
            cm.add_standard_input( open(params.input.params).read().split('\n') )

        cm.add_standard_input( ['END'] )

    # Pass additional command line arguments?
    if params.settings.args:
        cm.add_command_line_arguments( params.input.args )

    ########################
    # Print + Run
    cm.print_settings()
    out = cm.run()
    print cm.output

    if out != 0:
        print cm.error

    ########################
    # Write log
    with open(output_prefix+'.quick.log', 'w') as fh:
        fh.write('\n\nSTDOUT\n\n')
        fh.write(cm.output)
        fh.write('\n\nSTDERR\n\n')
        fh.write(cm.error)

    ########################
    # Link output
    real_pdb = glob.glob(output_prefix+'*.pdb')
    link_pdb = params.output.link_prefix+'.pdb'
    if real_pdb:
        real_pdb = real_pdb[0]
        if os.path.exists(link_pdb) and os.path.islink(link_pdb):
            os.unlink(link_pdb)
        if (not os.path.exists(link_pdb)):
            os.symlink(real_pdb, link_pdb)

    real_mtz = glob.glob(output_prefix+'*.mtz')
    link_mtz = params.output.link_prefix+'.mtz'
    if real_mtz:
        real_mtz = real_mtz[0]
        if os.path.exists(link_mtz) and os.path.islink(link_mtz):
            os.unlink(link_mtz)
        if (not os.path.exists(link_mtz)):
            os.symlink(real_mtz, link_mtz)

    ########################
    # Create occupancy parameters
    if os.path.exists(link_pdb) and params.settings.create_occupancy_params and params.input.params:
        print 'Running giant.create_occupancy_params'
        print 'Generating occupancy model for the input pdb: {}'.format(link_pdb)
        os.system('giant.create_occupancy_params pdb={} verbose=False'.format(link_pdb))

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)

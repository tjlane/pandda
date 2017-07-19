#!/usr/bin/env ccp4-python

import os, sys, copy, glob

import numpy, pandas, json

import libtbx.phil, libtbx.easy_mp
import iotbx.pdb

from libtbx.utils import Sorry, Failure

from bamboo.common.logs import Log
from bamboo.common.path import easy_directory, rel_symlink
from bamboo.common.command import CommandManager

from giant.structure.select import protein, backbone, sidechains
from giant.structure.formatting import PhenixSelection

import matplotlib
matplotlib.interactive(False)
from matplotlib import pyplot
pyplot.switch_backend('agg')
pyplot.interactive(0)

numpy.set_printoptions(threshold=numpy.nan)

from IPython import embed

############################################################################

PROGRAM = 'pandemic.adp_sweep'

DESCRIPTION = """
    Fit a series of B-factor models to a series of B-factor refinements
"""

############################################################################

blank_arg_prepend = {None:'refinement_dir=', '.pdb':'pdb='}

master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .help = "input pdb files for refinement"
        .multiple = True
        .type = str
    cif = None
        .type = str
    refinement_dir = None
        .help = "input directory, containing different B-factor refinements"
        .optional = False
        .type = str
    labelling = filename *foldername
        .type = choice
}
output {
    out_dir = b-factor-fitting-sweep
        .help = "output directory"
        .type = str
}
options {
    tls_selections = tls_selections.params
        .type = path
        .multiple = False
    refinements = *isotropic *tls *anisotropic
        .type = choice(multi=True)
    fitting_groups = *tls *chain *residue *backbone_sidechain
        .type = choice(multi=True)
}
parameterisation {
    tbd = True
}
settings {
    cpus = 60
        .type = int
        .multiple = False
    cpus_per_job = 5
        .type = int
        .multiple = False

}
""", process_includes=True)

############################################################################

def wrapper_run(arg):
    return arg.run()

def validate_parameters(params):

    assert params.settings.cpus > params.settings.cpus_per_job
    assert params.settings.cpus % params.settings.cpus_per_job == 0, 'cpus must be integer multiple of cpus_per_job'

def run_refinements(params, log=None):

    if log is None: log = Log(verbose=True)
    log.heading('Creating B-factor refinement job')

    cmd = CommandManager('multi.adp_refine')
    cmd.add_command_line_arguments(params.input.pdb)
    cmd.add_command_line_arguments(params.input.cif)
    cmd.add_command_line_arguments(r'tls_selections={}'.format(params.options.tls_selections))
    cmd.add_command_line_arguments(r'labelling={}'.format(params.input.labelling))
    cmd.add_command_line_arguments(r'b_factor_models={}'.format('+'.join(params.options.refinements)))
    cmd.add_command_line_arguments(r'out_dir={}'.format(params.input.refinement_dir))
    cmd.add_command_line_arguments(r'log_file={}.log'.format(params.input.dir))
    cmd.add_command_line_arguments(r'cpus={}'.format(params.settings.cpus))
    cmd.print_settings()
    ret = cmd.run()
    if ret != 0:
        raise Exception('Failed in refinement')

def create_sweep_tls_group_choice(input_pdbs, out_dir, params, log=None):

    if log is None: log = Log(verbose=True)

    out_dir = easy_directory(out_dir)

    jobs = []

    base_manager = CommandManager("pandemic.adp")
    base_manager.add_command_line_arguments(input_pdbs)
    base_manager.add_command_line_arguments(r'labelling=foldername')
    base_manager.add_command_line_arguments(r'cpus={}'.format(params.settings.cpus_per_job))

    # Load first structure for creating TLS groups
    h = protein(iotbx.pdb.hierarchy.input(input_pdbs[0]).hierarchy)

    for mode in params.options.fitting_groups:
        # Compile command object
        cmd = copy.deepcopy(base_manager)
        jobs.append(cmd)
        # Create output directory
        sweep_dir = os.path.join(out_dir, 'group_{}'.format(mode))
        cmd.add_command_line_arguments(r'out_dir={}'.format(sweep_dir))
        # Create TLS groups for fitting
        tls_groups = []
        if   mode == "chain":
            for chain in h.chains():
                #tls_groups.append("chain '{}'".format(chain.id))
                tls_groups.append(PhenixSelection.format(chain))
        elif mode == "residue":
            for rg in h.residue_groups():
                tls_groups.append(PhenixSelection.format(rg))
        elif mode == "tls":
            for gp in open(params.options.tls_selections, 'r').read().strip().split('\n'):
                tls_groups.append(gp)
        elif mode == "backbone_sidechain":
            for ag in backbone(h).atom_groups():
                if ag.resname in ['ALA','GLY']:
                    tls_groups.append(PhenixSelection.format(ag))
                else:
                    tls_groups.append(PhenixSelection.format(ag)+" and (name C or name CA or name N or name O)")
            for ag in sidechains(h).atom_groups():
                if ag.resname in ['ALA','GLY']:
                    continue
                else:
                    tls_groups.append(PhenixSelection.format(ag)+" and not (name C or name CA or name N or name O)")

        assert tls_groups, 'No TLS groups found: mode={}'.format(mode)

        # Add groups to command object
        for g in sorted(tls_groups):
            cmd.add_command_line_arguments(r'tls_group="{}"'.format(g))

    assert len(params.options.fitting_groups) == len(jobs)

    for mode, cmd in zip(params.options.fitting_groups, jobs):
        log.subheading(mode)
        cmd.print_settings()

    return jobs

def create_sweep_refinement_type(params, log=None):

    if log is None: log = Log(verbose=True)

    jobs = []

    for ref_type in params.options.refinements:

        log.heading('Creating jobs to analyse {}-refined structures'.format(ref_type))

        pdbs = sorted(glob.glob(os.path.join(params.input.refinement_dir, '*', '*-{}.pdb'.format(ref_type))))
        log('Found {} files.'.format(len(pdbs)))
        out_dir = easy_directory(os.path.join(params.output.out_dir, ref_type))
        log('Outputting to directory: {}'.format(out_dir))

        jobs.extend(create_sweep_tls_group_choice(input_pdbs=pdbs, out_dir=out_dir, params=params))

    return jobs

def run_sweeps(jobs, params, log=None):

    if log is None: log = Log(verbose=True)

    n_jobs = params.settings.cpus // params.settings.cpus_per_job

    #if n_jobs > len(jobs):
    #    log('More process available than cores being used. Increasing cpus_per_job')

    log.heading('Running {} jobs ({} at a time with {} cores each)'.format(len(jobs), n_jobs, params.settings.cpus_per_job))

    ret = libtbx.easy_mp.pool_map(processes=n_jobs, func=wrapper_run, args=jobs, chunksize=1)

#    for i, j in enumerate(ret):
#        j.write_output('job-{}.log'.format(i+1))

############################################################################

def run(params):

    easy_directory(params.output.out_dir)

    validate_parameters(params)

    if params.input.pdb: run_refinements(params)

    jobs = create_sweep_refinement_type(params)

    run_sweeps(jobs, params)

    #embed()

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



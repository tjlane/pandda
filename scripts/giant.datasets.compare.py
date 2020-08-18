#!/usr/bin/env cctbx.python

import os, sys
import math
import numpy, pandas

import libtbx.phil
import iotbx.pdb

from libtbx.utils import Sorry, Failure

from scitbx.matrix import sym
from scitbx.array_family import flex

from bamboo.common.logs import Log
from bamboo.common.path import filename, foldername, easy_directory

from giant.dataset import AtomicModel
from giant.structure.align import align_chains_rigid, GlobalAlignment

############################################################################

PROGRAM = os.path.basename(__file__)

DESCRIPTION = """
    A tool to compare structure of proteins from similar crystals.
"""
############################################################################

blank_arg_prepend = {'.pdb':'pdb=', '.csv':'csv='}

master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .type = path
        .multiple = True
    labelling = *filename foldername
        .type = choice
}
output {
    out_dir = dataset_comparison
        .type = str
}
""")

############################################################################

def rms(vals, axis=None):
    return numpy.sqrt(numpy.mean(numpy.power(vals,2), axis=axis))

def transform_uij(r, sym_mat3_double):
    # Calculate inverse
    rt = r.inverse()
    # Convert to matrices for multiplication
    sym_mat3 = [sym(sym_mat3=v) for v in sym_mat3_double]
    # Apply matrix transformation
    trans = [r*v*rt for v in sym_mat3]
    # Convert back to sym_mat3_double
    out = flex.sym_mat3_double([v.as_sym_mat3() for v in trans])
    return out

def b_to_u(b):
    iso_u = b / (8*math.pi*math.pi)
    iso_u_arr = numpy.zeros((iso_u.size(),6))
    for i in range(3):
        iso_u_arr[:,i] = iso_u
    iso_u_arr = flex.sym_mat3_double(iso_u_arr.tolist())
    return iso_u_arr

def run(params):

    log = Log()

    #if os.path.exists(params.output.out_dir):
    #    raise Sorry('Output directory already exists, please delete or change output directory: {}'.format(params.output.out_dir))
    #
    #out_dir = easy_directory(params.output.out_dir)

    if params.input.labelling == 'filename':
        lab_func = filename
    else:
        lab_func = foldername

    reference_h = None

    for p in params.input.pdb:

        lab = lab_func(p)

        if (lab is None) or (lab == ''):
            raise Sorry('The selected labelling function ({}) does not return a valid label ({}) for input file {}'.format(params.input.labelling, lab, p))

        m = AtomicModel.from_file(p).label(tag=lab_func(p))

        log('Processing {}'.format(m.tag))

        if reference_h is None:
            log('Selected this as reference: {}'.format(m.tag))
            reference_h = m.hierarchy
            atom_hash = [a.id_str() for a in reference_h.atoms()]
            continue

        # Align each dataset to the reference
        log('Aligning to Reference')
        lsq_rt, alignment_sites, reference_sites = align_chains_rigid(
                mov_chain=m.hierarchy.models()[0].chains()[0],
                ref_chain=reference_h.models()[0].chains()[0]
                )
        m.alignment = GlobalAlignment(alignment_mx=lsq_rt, alignment_sites=alignment_sites, reference_sites=reference_sites, id=None)

        # Rotate xyz & Uij to new orientation
        rot_hierarchy = m.hierarchy.deep_copy()
        # Extract & transform xyz
        raw_xyz = rot_hierarchy.atoms().extract_xyz()
        rot_xyz = m.alignment.nat2ref(raw_xyz)
        # Extract & transform uij
        raw_uij = rot_hierarchy.atoms().extract_uij()
        # Atoms wth uijs or B-iso
        sel_u = flex.bool((numpy.array(raw_uij).reshape((raw_uij.size(), 6)) != -1.0).all(axis=1).tolist())
        sel_b = flex.bool(raw_uij.size(), True).set_selected(sel_u, False)
        log('{} atoms have U'.format(sum(sel_u)))
        log('{} atoms have B only'.format(sum(sel_b)))
        # Transform Uij
        u_trans = transform_uij(r=lsq_rt.r, sym_mat3_double=raw_uij.select(sel_u))
        rot_uij = raw_uij.deep_copy().set_selected(sel_u, u_trans)
        # Apply
        rot_hierarchy.atoms().set_xyz(rot_xyz)
        rot_hierarchy.atoms().set_uij(rot_uij)

        out_file = os.path.splitext(p)[0]+'-aligned.pdb'
        log('Writing output: {}'.format(out_file))
        rot_hierarchy.write_pdb_file(out_file)

        # Calculate differences to reference
        dif_hierarchy = rot_hierarchy.deep_copy()
        # Reset B+U values
        dif_atoms = dif_hierarchy.atoms()
        dif_atoms.set_b(flex.double(dif_atoms.size(), 0))
        dif_atoms.set_uij(flex.sym_mat3_double(dif_atoms.size(), (-1.0,)*6))

        pairs = []
        for i, a in enumerate(dif_atoms):
            try:
                idx = atom_hash.index(a.id_str())
                # Compare XYZ Coordinates?
                pairs.append((i, idx))
            except: pass
        # Common selections between structures -- THESE ARE ISELECTIONS
        sel_this, sel_ref = map(flex.size_t, zip(*pairs))

        # Atoms not to be compared
        not_sel_this = flex.bool(raw_uij.size(), True).set_selected(sel_this, False)
        # Atoms to be compared in B or U
        comp_b = sel_b.deep_copy().set_selected(not_sel_this, False)
        comp_u = sel_u.deep_copy().set_selected(not_sel_this, False)
        assert sum(comp_b) + sum(comp_u) == len(sel_this)

        # Masks of U/B on the selection of sel_this
        b_mask = comp_b.select(sel_this)
        u_mask = comp_u.select(sel_this)
        assert b_mask.size() == sel_this.size()

        dif_u = (
                reference_h.atoms().select(sel_ref).extract_uij() -
                rot_hierarchy.atoms().select(sel_this).extract_uij()
                ).select(u_mask)
        rms_uij = flex.double(numpy.abs(numpy.array(dif_u).reshape((dif_u.size(), 6))[:,0:3]).mean(axis=1))
        #rms_uij = flex.double(map(rms, dif_u))

        dif_b = flex.abs(
                reference_h.atoms().select(sel_ref).extract_b() -
                rot_hierarchy.atoms().select(sel_this).extract_b()
                ).select(b_mask)

        # Set these in the difference hierarchy
        dif_hierarchy.atoms().select(comp_u).set_b(8*math.pi*math.pi*rms_uij)
        dif_hierarchy.atoms().select(comp_b).set_b(dif_b)

        nrm_hierarchy = dif_hierarchy.deep_copy()
        nrm_b = nrm_hierarchy.atoms().extract_b() / rot_hierarchy.atoms().extract_b()
        nrm_b.set_selected((rot_hierarchy.atoms().extract_b() == 0.0), 0.0)
        nrm_hierarchy.atoms().set_b(nrm_b)

        iso_u = b_to_u(dif_hierarchy.atoms().extract_b())
        dif_hierarchy.atoms().set_uij(iso_u)

        out_file = os.path.splitext(p)[0]+'-difference.pdb'
        log('Writing output: {}'.format(out_file))
        dif_hierarchy.write_pdb_file(out_file)

        iso_u = b_to_u(nrm_hierarchy.atoms().extract_b())
        nrm_hierarchy.atoms().set_uij(iso_u)

        out_file = os.path.splitext(p)[0]+'-difference-normalised.pdb'
        log('Writing output: {}'.format(out_file))
        nrm_hierarchy.write_pdb_file(out_file)

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

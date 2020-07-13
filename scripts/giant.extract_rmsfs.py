#!/usr/bin/env cctbx.python

import giant.logs as lg
logger = lg.getLogger(__name__)

from giant.exceptions import Sorry, Failure

import os, sys, json, collections

#######################################

blank_arg_prepend = {
    '.pdb' : 'input.pdb=',
}

import libtbx.phil
master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .type = path
        .multiple = True
}
""")

#######################################

def run(params):

    lg.setup_logging(
        name = __name__,
    )

    for pdb_file in params.input.pdb:

        logger.heading(pdb_file)

        pdb_root = os.path.splitext(pdb_file)[0]

        import iotbx.pdb
        h = iotbx.pdb.hierarchy.input(pdb_file).hierarchy

        if h.models_size() == 1:
            logger('Only 1 model!')
            continue

        atom_coords = collections.OrderedDict()

        for model in h.models():

            logger('Reading model {}'.format(model.id))

            for a in model.atoms():

                lab = a.pdb_label_columns()
                xyz = a.xyz

                atom_coords.setdefault(lab, []).append(xyz)

        logger.subheading('Calculating RMSFs')
        from scitbx.array_family import flex
        rmsf_dict = collections.OrderedDict()
        for lab, coords in atom_coords.iteritems():

            coords = flex.vec3_double(coords)
            rmsf = (coords - coords.mean()).rms_length()

            rmsf_dict[lab] = rmsf

            logger('{}\t{:.1f}'.format(lab, rmsf))

        rmsf_json_file = (pdb_root + '-rmsf.json')
        with open(rmsf_json_file, 'w') as fh:
            fh.write(
                json.dumps(rmsf_dict, indent=2)
            )

        h_copy = h.deep_copy()
        h_copy.atoms().set_b(flex.double(h_copy.atoms_size(), 0.0))
        for a in h_copy.atoms():
            a.b = rmsf_dict.get(a.pdb_label_columns(), 0)

        rmsf_pdb_file = (pdb_root + '-rmsf.pdb')
        h_copy.write_pdb_file(rmsf_pdb_file)

if __name__ == '__main__':
    from giant.jiffies import run_default
    run_default(
        run = run,
        master_phil = master_phil,
        args = sys.argv[1:],
        blank_arg_prepend = blank_arg_prepend,
    )

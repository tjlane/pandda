import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys

from giant.exceptions import Sorry, Failure

#######################################

blank_arg_prepend = {
    '.pdb' : 'align_pdb=',
}

import libtbx.phil
master_phil = libtbx.phil.parse("""
input {
    reference_pdb = None
        .type = str
        .help = "Structure to align all other structures to."
        .optional = False
    align_pdb = None
        .type = str
        .multiple = True
    labelling = *filename foldername
        .type = choice(multi=False)
        .help = "how should a label be extracted from the input file(s)? Used for naming the output files."
}
output {
    out_dir = .
        .type = path
        .help = "output directory"
    prefix = ''
        .type = str
    suffix = '.aligned'
        .type = str
}
alignment {
    reference_selection = "chain A"
        .type = str
        .help = "Selection of the reference structure to align the other things to"
    align_selections = None
        .multiple = True
        .optional = True
    minimum_alignment = 30
        .type = int
        .help = 'Minimum number of residues to be used for each alignment'
    minimum_alignment_fraction = 0.50
        .type = float
        .help = 'Minimum fraaction of the sequences that must align for two things to be aligned'
    minimum_identity = 0.00
        .type = float
        .help = 'Minimum sequence identity for two things to be aligned'
}
settings {
    verbose = False
        .type = bool
}
""")


class Aligner(object):
    """Takes two chains and aligns them - return rt_mx"""

    def __init__(self,
        reference_hierarchy,
        reference_selection = None,
        minimum_alignment = 0,
        minimum_alignment_fraction = 0.00,
        minimum_identity = 0.00,
        ):

        if (reference_selection is not None):
            asc = reference_hierarchy.atom_selection_cache()
            selection = asc.selection(reference_selection)
            reference_hierarchy = reference_hierarchy.select(selection)

        r_seq, r_sites, r_flags = self.extract_sites_for_alignment(
            hierarchy = reference_hierarchy,
        )

        self.ref_hierarchy = reference_hierarchy
        self.ref_seq   = r_seq
        self.ref_sites = r_sites
        self.ref_flags = r_flags

        self.minimum_alignment = minimum_alignment
        self.minimum_alignment_fraction = minimum_alignment_fraction
        self.minimum_identity = minimum_identity

    def __call__(self, align_hierarchy, align_selection=None):

        if align_selection is not None:
            asc = align_hierarchy.atom_selection_cache()
            selection = asc.selection(align_selection)
            sel_align_hierarchy = align_hierarchy.select(selection)
        else:
            sel_align_hierarchy = align_hierarchy

        rt = self.get_alignment(
            hierarchy = sel_align_hierarchy,
            )

        if rt is None:
            raise Sorry('No alignment generated')

        rot_h = self.apply_alignment(
            hierarchy = align_hierarchy,
            rt = rt,
        )

        return rot_h

    def apply_alignment(self, hierarchy, rt):

        import numpy
        import scitbx.matrix
        from scitbx.array_family import flex

        r1 = rt.r
        r2 = rt.r.transpose()

        out_h = hierarchy.deep_copy()
        out_a = out_h.atoms()

        xyz_in = out_a.extract_xyz()
        uij_in = out_a.extract_uij()

        xyz_out = rt * xyz_in

        uij_d = uij_in.as_double()
        uij_d.reshape(flex.grid((uij_in.size(), 6)))
        uij_n = uij_d.as_numpy_array()
        uij_sel = flex.bool(numpy.logical_not((uij_n == -1).all(axis=1)))
        uij_rt = flex.sym_mat3_double([
            (r1 * scitbx.matrix.sym(sym_mat3=uu) * r2).as_sym_mat3()
            for uu in uij_in.select(uij_sel)
        ])

        uij_out = uij_in.deep_copy()
        uij_out.set_selected(uij_sel, uij_rt)

        out_a.set_xyz(xyz_out)
        out_a.set_uij(uij_out)

        return out_h

    def get_alignment(self, hierarchy):

        ref_seq, ref_sites, ref_flags = (self.ref_seq, self.ref_sites, self.ref_flags)
        mov_seq, mov_sites, mov_flags = self.extract_sites_for_alignment(
            hierarchy = hierarchy,
        )

        import mmtbx.alignment
        align_obj = mmtbx.alignment.align(
            seq_a = ref_seq,
            seq_b = mov_seq,
            gap_opening_penalty = 20,
            gap_extension_penalty = 2,
            similarity_function = 'blosum50',
            style = 'local',
        )

        # Extract the alignment
        alignment = align_obj.extract_alignment()
        # List of matches - '|' for exact match, '*' for good match
        matches = alignment.matches()
        equal = matches.count("|")
        similar = matches.count("*")
        alignment.pretty_print(
            matches = matches,
            block_size = 50,
            n_block = 1,
            top_name = "fixed",
            bottom_name = "moving",
            )

        identity = alignment.calculate_sequence_identity()
        if (identity < self.minimum_identity):
            logger(
                'Alignment sequence identity below threshold: {} < {}'.format(
                    identity,
                    self.minimum_identity,
                )
            )
            return None

        num_matches = equal + similar
        if (num_matches < self.minimum_alignment):
            logger(
                'Alignment length below threshold: {} < {}'.format(
                    num_matches,
                    self.minimum_alignment,
                )
            )
            return None

        fraction = float(equal+similar) / float(alignment.match_codes.count('m'))
        if (fraction < self.minimum_alignment_fraction):
            logger(
                'Alignment fraction below threshold: {} < {}'.format(
                    fraction,
                    self.minimum_alignment_fraction,
                )
            )
            return None

        # Create list of selected sites
        from scitbx.array_family import flex
        ref_sites_sel = flex.vec3_double()
        mov_sites_sel = flex.vec3_double()
        for ia,ib,m in zip(alignment.i_seqs_a, alignment.i_seqs_b, matches):
            if (m not in ["|", "*"]):
                continue
            # Check that the sites are flagged to be used
            if (ref_flags[ia] and mov_flags[ib]):
                # Append sites to list to align
                ref_sites_sel.append(ref_sites[ia])
                mov_sites_sel.append(mov_sites[ib])

        if (ref_sites_sel.size() == 0):
            raise Sorry("No matching C-alpha atoms.")

        from scitbx.math import superpose
        lsq_rt = superpose.least_squares_fit(
            reference_sites = ref_sites_sel,
            other_sites = mov_sites_sel,
        ).rt()

        return lsq_rt

    def extract_sites_for_alignment(self, hierarchy):
        """Extract sequence and sites of c-alphas - adapted from mmtbx.command_line.super"""

        import iotbx.pdb
        import iotbx.pdb.amino_acid_codes
        from scitbx.array_family import flex

        seq = []
        sites = flex.vec3_double()
        use_sites = flex.bool()
        model = hierarchy.models()[0]
        for chain in model.chains():
            for resi in chain.conformers()[0].residues():
                if (iotbx.pdb.common_residue_names_get_class(name=resi.resname) != "common_amino_acid"):
                    continue
                resn = resi.resname
                single = iotbx.pdb.amino_acid_codes.one_letter_given_three_letter[resn]
                seq.append(single)
                use = False
                xyz = (0,0,0)
                for atom in resi.atoms():
                    if (atom.name == " CA "):
                        xyz = atom.xyz
                        use = True
                        break
                sites.append(xyz)
                use_sites.append(use)

        return "".join(seq), sites, use_sites


#######################################

def validate_params(params):

    if (params.input.reference_pdb is None):
        raise Sorry('No reference structure provided (input.reference_pdb=XXX)')

    if (not params.input.align_pdb):
        raise Sorry('No structures provided for alignment (input.align_pdb=XXX)')

    for p in params.input.align_pdb:
        if not os.path.exists(p):
            raise IOError('Input file "{}" does not exists'.format(p))

def get_label(label_func, filename):

    label = str(label_func(filename))

    if (label is None) or (len(label) == 0):
        raise Sorry('No label generated for file: {}'.format(filename))

    return label

def run(params):

    from giant.paths import filename, foldername, easy_directory
    from giant.colours import pretty_string as ps

    logger = lg.setup_logging(
        name = __name__,
        log_file = 'giant_align_structures.log',
    )

    validate_params(params)

    if params.input.labelling == 'filename':
        label_func = filename
    elif params.input.labelling == 'foldername':
        label_func = foldername
    else:
        raise Failure('invalid input.labelling: {}'.format(params.input.labelling))

    out_dir = easy_directory(params.output.out_dir)

    logger(ps('giant.align_structures_with_adps\n').bold().blue())

    logger(ps('Aligning {} structure(s) to {}'.format(
        len(params.input.align_pdb),
        params.input.reference_pdb,
    )).bold())
    logger('\t{}'.format(
        '\n\t'.join(params.input.align_pdb),
    ))

    # =================================================>
    # Load structures and create aligned
    # =================================================>

    import iotbx.pdb

    h_ref = iotbx.pdb.hierarchy.input(params.input.reference_pdb).hierarchy

    aligner = Aligner(
        reference_hierarchy = h_ref,
        reference_selection = params.alignment.reference_selection,
        minimum_alignment = params.alignment.minimum_alignment,
        minimum_alignment_fraction = params.alignment.minimum_alignment_fraction,
        minimum_identity = params.alignment.minimum_identity,
        )

    for p_mov in params.input.align_pdb:

        logger.subheading('Aligning {}'.format(p_mov))

        h_mov = iotbx.pdb.hierarchy.input(p_mov).hierarchy

        mov_label = get_label(label_func, p_mov)

        if params.alignment.align_selections:
            align_selections = params.alignment.align_selections
            align_labels = ['sel{}'.format(i+1) for i in range(len(align_selections))]
        else:
            chain_ids = sorted(set([c.id for c in h_mov.chains()]))
            align_selections = ["chain {}".format(c_id) for c_id in chain_ids]
            align_labels = ['chain{}'.format(c_id) for c_id in chain_ids]

        for ali_label, ali_sel in zip(align_labels, align_selections):

            logger(ps('\nAligning by {}'.format(ali_label)).green())
            out_basename = params.output.prefix + mov_label + params.output.suffix + '.' + ali_label + '.pdb'
            out_filepath = os.path.join(params.output.out_dir, out_basename)

            if os.path.exists(out_filepath):
                raise IOError('file already exists! ({})'.format(out_filepath))

            try:
                aligned_hierarchy = aligner(
                    align_hierarchy = h_mov,
                    align_selection = ali_sel,
                )
            except Sorry as e:
                logger.warning(ps('Could not align selection {}: {}'.format(ali_sel, str(e))).red())
                continue

            if (aligned_hierarchy is not None):
                aligned_hierarchy.write_pdb_file(
                    file_name = out_filepath,
                )

    logger.subheading(ps('Finished normally!').green())

#######################################

if __name__ == '__main__':
    from giant.jiffies import run_default
    run_default(
        run = run,
        master_phil = master_phil,
        args = sys.argv[1:],
        blank_arg_prepend = blank_arg_prepend,
    )

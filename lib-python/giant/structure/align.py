import giant.logs as lg
logger = lg.getLogger(__name__)

from scitbx.array_family import flex

from giant.exceptions import Sorry, Failure
from giant.structure.select import protein, backbone, complete_backbone, common_residues, extract_atom


class AlignmentError(Exception):
    pass


class Alignment(object):

    def __init__(self,
        mov_label = None,
        ref_label = None,
        ):

        # IDs of moving and reference (e.g.) chains (if defined)
        self.mov_label = mov_label
        self.ref_label = ref_label

    def nat2ref(self, *args, **arg_dict):
        raise NotImplementedError('Function not implemented for class of this type: {}'.format(self.__class__))

    def ref2nat(self, *args, **arg_dict):
        raise NotImplementedError('Function not implemented for class of this type: {}'.format(self.__class__))

    def summary(self, *args, **arg_dict):
        raise NotImplementedError('Function not implemented for class of this type: {}'.format(self.__class__))

    def alignment_rmsd(self, *args, **arg_dict):
        raise NotImplementedError('Function not implemented for class of this type: {}'.format(self.__class__))


class CoordinateAlignment(Alignment):

    def __init__(self,
        alignment_sites = None,
        reference_sites = None,
        mov_label = None,
        ref_label = None,
        ):
        """An alignment derived from two sets of coordinates"""

        super(CoordinateAlignment, self).__init__(
            mov_label = mov_label,
            ref_label = ref_label,
        )

        self.n_sites = 0
        self.reference_sites = None
        self.alignment_sites_nat = None
        self.alignment_sites_ref = None

        if (alignment_sites is not None):
            # Sites of the moving atoms before and after alignment
            self.alignment_sites_nat = flex.vec3_double(alignment_sites)
            self.n_sites = len(self.alignment_sites_nat)

        # Sites used for alignment (defines reference frame)
        if (reference_sites is not None):
            self.reference_sites = flex.vec3_double(reference_sites)
            self.n_sites = len(self.reference_sites)

        if (reference_sites is not None) and (alignment_sites is not None):
            assert self.n_sites == len(self.reference_sites)
            assert self.n_sites == len(self.alignment_sites_nat)

    def alignment_rmsd(self):

        return self.reference_sites.rms_difference(self.alignment_sites_ref)


class GlobalAlignment(CoordinateAlignment):

    def __init__(self,
        alignment_mx,
        alignment_sites = None,
        reference_sites = None,
        mov_label = None,
        ref_label = None,
        ):
        """
        Object to hold the output of a global (rigid) alignment - performs coordinate transformations using a single alignment matrix

        alignment_mx    : alignment matrix.
        alignment_sites : sites used for the alignment (in the "native frame")
        reference_sites : corresponding reference sites used for the alignment (in the "reference" frame)
        """

        super(GlobalAlignment, self).__init__(
            alignment_sites = alignment_sites,
            reference_sites = reference_sites,
            mov_label = mov_label,
            ref_label = ref_label,
        )

        # Store alignment matrix
        self.alignment_mx = alignment_mx

        ### Should now be able to use inherited methods

        if (self.alignment_sites_nat is not None):
            self.alignment_sites_ref = flex.vec3_double(
                self.nat2ref(coordinates=self.alignment_sites_nat)
            )

    def nat2ref(self, coordinates, **kw_args):

        if not isinstance(coordinates, flex.vec3_double):
            coordinates = flex.vec3_double(coordinates)

        return self.alignment_mx * coordinates

    def ref2nat(self, coordinates, **kw_args):

        if not isinstance(coordinates, flex.vec3_double):
            coordinates = flex.vec3_double(coordinates)

        return self.alignment_mx.inverse() * coordinates

    def summary(self):

        s = (
            'Global Alignment:\n' +
            '  Alignment labels: {} (Moving) - {} (Reference)\n'.format(self.mov_label, self.ref_label) +
            '  Sites used for alignment: {}\n'.format(self.n_sites) +
            '  Post-alignment RMSD: {}\n'.format(round(self.alignment_rmsd(),2))
        )

        return s


class LocalAlignment(CoordinateAlignment):

    def __init__(self,
        alignment_mxs,
        alignment_sites,
        reference_sites = None,
        mov_label = None,
        ref_label = None,
        ):
        """
        Object to hold the output of a local alignment - performs coordinate transformations using the alignment matrices

        alignment       : list of alignment matrices.
        alignment_sites : list of sites associated with the alignment matrices
        """

        assert (alignment_sites is not None)

        super(LocalAlignment, self).__init__(
            alignment_sites = alignment_sites,
            reference_sites = reference_sites,
            mov_label = mov_label,
            ref_label = ref_label,
        )

        # List of alignment matrices attached to each site
        self.alignment_mxs = list(alignment_mxs)

        # Therefore must be same length
        assert len(alignment_mxs) == self.n_sites

        self._populate()

    def _populate(self):

        assert (self.alignment_sites_nat is not None)

        # Transform alignment sites to reference frame
        self.alignment_sites_ref = flex.vec3_double([
            m * s for m,s in zip(self.alignment_mxs, self.alignment_sites_nat)
        ])

        # Sites for searching (centres of voronoi cells)
        self._tree_sites = {
            'nat': self.alignment_sites_nat,
            'ref': self.alignment_sites_ref,
        }

    def _tree(self, coordinate_frame):

        import scipy.spatial

        assert coordinate_frame in list(self._tree_sites.keys())

        return scipy.spatial.KDTree(
            data = self._tree_sites[coordinate_frame],
        )

    def _get_alignment_indices_for_points(self, coordinates, coordinate_frame):
        """Get the index of the closes alignment_site for each point"""

        tree = self._tree(
            coordinate_frame = coordinate_frame,
        )

        return tree.query(coordinates)[1]

    def _validate_mappings(self, mappings):

        assert not (mappings<0).nonzero()[0].any(), 'Invalid alignment mapping: Cannot use mappings that are less than zero'
        assert not (mappings>self.n_sites).nonzero()[0].any(), 'Invalid alignment mapping: Trying to use more mappings than there are alignments'

    def nat2ref(self,
        coordinates,
        mappings = None,
        **kw_args
        ):

        if mappings is None:
            mappings = self._get_alignment_indices_for_points(
                coordinates = coordinates,
                coordinate_frame = 'nat',
            )
        else:
            self._validate_mappings(mappings)

        return transform_coordinates_with_multiple_alignments(
            coordinates = coordinates,
            alignments = self.alignment_mxs,
            mappings = mappings,
            inverse = False,
        )

    def ref2nat(self,
        coordinates,
        mappings = None,
        **kw_args
        ):

        if mappings is None:
            mappings = self._get_alignment_indices_for_points(
                coordinates = coordinates,
                coordinate_frame = 'ref',
            )
        else:
            self._validate_mappings(mappings)

        return transform_coordinates_with_multiple_alignments(
            coordinates = coordinates,
            alignments = self.alignment_mxs,
            mappings = mappings,
            inverse = True,
        )

    def summary(self):

        s = (
            'Local alignment:\n' +
            '  Alignment labels: {} (Moving) - {} (Reference)\n'.format(self.mov_label, self.ref_label) +
            '  Sites used for alignment: {}\n'.format(self.n_sites) +
            '  Post-alignment RMSD: {}\n'.format(round(self.alignment_rmsd(),2))
        )

        return s


class MultipleLocalAlignment(LocalAlignment):

    def __init__(self, local_alignments=[], id=None):
        """
        Object to hold the output of multiple local alignments (such as multiple chains per structure).
        Performs coordinate transformations using the alignment matrices supplied as if it were one LocalAlignment

        local_alignments       : list of LocalAlignment objects.
        """

        super(MultipleLocalAlignment, self).__init__(
            alignment_mxs = [],
            alignment_sites = flex.vec3_double(), # need to init to empty list
            reference_sites = None,
            mov_label = None,
            ref_label = None,
        )

        # List of child alignments
        self.local_alignments = []

        # Add alignments to the class
        for l in local_alignments:
            self.add(l, update=False)

        # Update internal objects
        self._update()

    def _update(self):
        """Refresh all of the attributes to make sure they are up-to-date"""

        # Reset these to their empty values
        self.alignment_mxs = []
        self.alignment_sites_nat = flex.vec3_double()
        self.alignment_sites_ref = flex.vec3_double()
        self.reference_sites = flex.vec3_double()


        # Add data from the alignments
        for l in self.local_alignments:
            self.alignment_mxs.extend(l.alignment_mxs)
            self.alignment_sites_nat.extend(l.alignment_sites_nat)
            self.alignment_sites_ref.extend(l.alignment_sites_ref)
            self.reference_sites.extend(l.reference_sites)

        # Class variables
        self.n_alignments = len(self.local_alignments)            
        self.n_sites = len(self.reference_sites)

        # Create new query trees
        self._tree_sites = {
            'nat': self.alignment_sites_nat,
            'ref': self.alignment_sites_ref,
        }

    def add(self, local_alignment, update=True):

        assert isinstance(local_alignment, LocalAlignment), 'Supplied object must be a LocalAlignment'

        assert (local_alignment not in self.local_alignments), 'Local alignment has already been added to object'

        self.local_alignments.append(local_alignment)

        if (update is True):
            self._update()

    def summary(self):

        s = (
            'Multiple local alignment:\n' +
            '  Number of alignments: {}\n'.format(self.n_alignments) +
            '  Alignment Summaries:\n'
        )
        for la in self.local_alignments:
            s += (
                '  ' +
                la.summary().replace('\n','\n  ') +
                '\n'
            )
        return s.strip('\n')


def nearby_coords_bool(query, coords, cutoff):
    """Find all points in coords within cutoff of query. Return boolean selection"""

    assert isinstance(coords, flex.vec3_double)
    return (coords-query).dot() < (cutoff**2)

def align_structures_flexible(
    mov_hierarchy,
    ref_hierarchy,
    altlocs = ['','A'],
    cutoff_radius = 15,
    sequence_identity_threshold = 0.95,
    one_to_one_mapping = True,
    require_hierarchies_identical = True,
    ):
    """
    Perform a flexible alignment on two hierarchies. Alignments are performed on a chain-by-chain basis.
    Each chain of mov_hierarchy is aligned
    """

    from giant.structure.sequence import \
            pairwise_chain_sequence_identity, \
            align_sequences_default

    # List of the alignments for each chain
    local_alignments = []

    # Trim to protein only
    mov_hierarchy = backbone(mov_hierarchy, copy=True)
    ref_hierarchy = backbone(ref_hierarchy, copy=True)

    # Check the structures only have one model
    try:
        mov_hierarchy.only_model()
        ref_hierarchy.only_model()
    except:
        raise Exception('Structures for alignment can only have one model!')

    # Check the structures are identical
    if require_hierarchies_identical:
        assert mov_hierarchy.is_similar_hierarchy(ref_hierarchy), 'Structures for alignment must have the same atoms (although atomic parameters can vary)'

    # Extract the chains from the structures
    c_mov = list(mov_hierarchy.chains())
    c_ref = list(ref_hierarchy.chains())

    # Match chains in the two structures (c_mov is first so the array is first indexed by the chains in mov)
    chn_ident = pairwise_chain_sequence_identity(
        c_mov,
        c_ref,
        seq_identity_threshold = None,
    )

    # Make the array boolean at the threshold value
    chn_sim = ( chn_ident > sequence_identity_threshold ).astype(int)

    # Create strings for use in case of errors/verbose printing
    logger.debug(
        '\n'.join([
            'Chain and sequences for aligment:',
            '{} chains in mov_hierarchy:'.format(len(c_mov)),
        ] + [
            '\t{}: {}'.format(
                c.id,
                ''.join(c.as_sequence()),
            ) for c in c_mov
        ] + [
            '{} chains in ref_hierarchy:'.format(len(c_ref)),
        ] + [
            '\t{}: {}'.format(
                c.id,
                ''.join(c.as_sequence()),
            ) for c in c_ref
        ] + [
            'Pairwise chain-by-chain sequence identities:',
            '     REFERENCE CHAINS --->',
            'MOV  {}'.format(
                ' '.join(['{:4}'.format(c.id) for c in c_ref]),
            ),
        ] + [
            '{:3}  {}'.format(
                i_c.id,
                ' '.join(['{:4}'.format(v) for v in chn_ident[i]]),
            ) for i,i_c in enumerate(c_mov)
        ] + [
            '',
            'Pairwise chain-by-chain sequence identities (thresholded at {}%):'.format(
                100 * sequence_identity_threshold,
            ),
            '     REFERENCE CHAINS --->',
            'MOV  {}'.format(
                ' '.join(['{:4}'.format(c.id) for c in c_ref]),
            ),
        ] + [
            '{:3}  {}'.format(
                i_c.id,
                ' '.join(['{:4}'.format(v) for v in chn_sim[i]]),
            ) for i,i_c in enumerate(c_mov)
        ])
    )

    # Iterate through and align the chains
    for i, chn_mov in enumerate(c_mov):

        # Skip if not protein
        if not chn_mov.is_protein():
            continue

        # Find the first chain in the reference structure that's "alignable"
        try:

            idx_ref = list(chn_sim[i]).index(1)
            chn_ref = c_ref[idx_ref]

            logger.debug(
                'Aligning chain {} of mov_hierarchy to chain {} in ref_hierarchy'.format(
                    chn_mov.id,
                    chn_ref.id,
                )
            )

            if (one_to_one_mapping is True):

                chn_sim[:,idx_ref] = 0

                logger.debug(
                    'Removing chain {} of ref_hierarchy from the pool of alignment chains (one_to_one_mapping is turned on)'.format(
                        chn_ref.id,
                    )
                )

        except ValueError:

            raise Failure(
                'Error raised during alignment.\n'
                'Unable to align chain {} from mov_hierarchy: there is no suitable chain in ref_hierarchy.\n'
                'This might be fixed by setting one_to_one_mapping to False or decreasing sequence_identity_threshold.\n'.format(chn_mov.id)
            )

        # Align the selected chains
        l_ali = align_chains_flexible(
            chn_mov = chn_mov,
            chn_ref = chn_ref,
            altlocs = altlocs,
            cutoff_radius = cutoff_radius,
        )

        # Add aligned chains as the ID of the LocalAlignment object
        l_ali.id = 'chain {} to chain {}'.format(
            chn_mov.id,
            chn_ref.id,
        )
        l_ali.mov_label = chn_mov.id
        l_ali.ref_label = chn_ref.id
        l_ali.sequence_alignment = align_sequences_default(
            seq_a = chn_ref.as_sequence(),
            seq_b = chn_mov.as_sequence(),
        )

        # Append to the alignments
        local_alignments.append(l_ali)

    # Print which chains were aligned to which
    logger.debug(
        '\n'.join(
            [
                'Alignment finished:'
            ] + [
                '\t(mov) chain {} aligned to (ref) chain {}'.format(
                    l_ali.mov_label,
                    l_ali.ref_label,
                ) for l_ali in local_alignments
            ]
        )
    )

    # Combine all of the local alignments
    return MultipleLocalAlignment(
        local_alignments = local_alignments,
    )

def align_chains_flexible(
    chn_mov,
    chn_ref,
    altlocs = ['','A'],
    cutoff_radius = 15,
    ):
    """
    Take two chains and perform flexible alignment on them.
    Only alternate conformations supplied in (e.g. altlocs=['','A']) will be used for alignment (maximum one conformer).
    Residues are removed that do not contain a full set of backbone atoms (N,CA,C,O) for the conformers selected (e.g. altlocs=['','A'])
    Chains will be truncated so that the chains contain an "aligned" set of residues (currently sequence-identical)

    returns LocalAlignment
    """

    import iotbx.pdb
    from scitbx.math import superpose

    # Trim both chains to residues with complete backbones
    chn_mov_cb = complete_backbone(chn_mov, altlocs=altlocs)
    chn_ref_cb = complete_backbone(chn_ref, altlocs=altlocs)

    # Trim both chains to the same set of residues
    chn_ref_cr, chn_mov_cr = common_residues(chn_ref_cb, chn_mov_cb)

    # Create new hierarchies to perform most processing
    h_mov = iotbx.pdb.hierarchy.new_hierarchy_from_chain(chn_mov_cr)
    h_ref = iotbx.pdb.hierarchy.new_hierarchy_from_chain(chn_ref_cr)

    h_mov.sort_atoms_in_place();
    h_ref.sort_atoms_in_place();

    # Extract new processed chain objects
    c_mov = h_mov.only_chain()
    c_ref = h_ref.only_chain()

    # Check that the chains contain the same atoms
    assert c_mov.atoms().extract_element() == c_ref.atoms().extract_element(), 'chn_mov and chn_ref must contain the same atoms'
    assert c_mov.atoms().extract_name()    == c_ref.atoms().extract_name(),    'chn_mov and chn_ref must contain the same atoms'

    # List of output alignments and alignment sites
    o_rts = []; o_xyz = []; r_xyz = []

    # Extract xyz coords
    xyz_mov = c_mov.atoms().extract_xyz()
    xyz_ref = c_ref.atoms().extract_xyz()

    # Iterate through and create an alignment for each C-alpha
    for rg_mov in c_mov.residue_groups():

        # Find the atoms near the C-alpha
        ca_atm = extract_atom(
            residue = rg_mov.conformers()[0].only_residue(),
            atom='CA',
        )
        nr_sel = nearby_coords_bool(
            query = ca_atm.xyz,
            coords = xyz_mov,
            cutoff = cutoff_radius,
        )

        # Select the sites from both chains
        xyz_mov_sel = xyz_mov.select(nr_sel)
        xyz_ref_sel = xyz_ref.select(nr_sel)

        # Calculate the alignment for this residue
        rt_atm = superpose.least_squares_fit(
            reference_sites = xyz_ref_sel,
            other_sites = xyz_mov_sel
        ).rt()

        # Save the rotation matrix and the coordinates of the c-alpha
        o_xyz.append(ca_atm.xyz)
        o_rts.append(rt_atm)
        r_xyz.append(xyz_ref_sel.select(((xyz_mov_sel-ca_atm.xyz).dot() == 0.0))[0])

    assert len(o_rts), o_rts
    assert len(o_xyz), o_xyz
    assert len(r_xyz), r_xyz

    # Return LocalAlignment object
    return LocalAlignment(
        alignment_mxs = o_rts,
        alignment_sites = o_xyz,
        reference_sites = r_xyz,
    )

def transform_coordinates_with_multiple_alignments(coordinates, alignments, mappings, inverse=False):
    """
    Apply the alignments to coordinates. Which alignment is used is controlled by mappings

    coordinates     : list of n coordinates
    alignments      : list of m alignments (or dict with m key/alignment pairs)
    mappings        : list of n indices (or keys), mapping each coordinate to an alignment matrix
    """

    import itertools
    import numpy

    assert len(coordinates) == len(mappings)

    if not isinstance(coordinates, flex.vec3_double):
        coordinates = flex.vec3_double(coordinates)

    # Sort the indices by the mapping values (allows coordinates to be grouped by the alignment used)
    num_coords      = len(coordinates)
    sorted_idxs     = flex.size_t(sorted(range(num_coords), key=lambda i: mappings[i]))
    sorted_coords   = coordinates.select(sorted_idxs)
    sorted_mappings = [mappings[i] for i in sorted_idxs]

    # Initialise output array
    out_coords = numpy.zeros(num_coords, dtype=[('x',float),('y',float),('z',float)])

    # Iterate through the coords in groups and transform
    for key_or_idx, sel_idxs in itertools.groupby(list(range(num_coords)), key=lambda i: sorted_mappings[i]):

        sel_idxs = flex.size_t(sel_idxs)
        sel_coords   = sorted_coords.select(sel_idxs)
        sel_mappings = [sorted_mappings[i] for i in sel_idxs]
        orig_idxs    = sorted_idxs.select(sel_idxs)

        assert max(sel_idxs)-min(sel_idxs) == len(sel_idxs)-1
        assert len(set(sel_mappings)) == 1

        rt = alignments[key_or_idx]

        if inverse:
            rt = rt.inverse()

        rt_coords = rt * sel_coords

        out_coords.put(orig_idxs, rt_coords)

    return flex.vec3_double(out_coords)

def extract_sites_for_alignment(chain_obj):
    """Extract sequence and sites of c-alphas - adapted from mmtbx.command_line.super"""

    import iotbx.pdb.amino_acid_codes

    seq = []
    sites = flex.vec3_double()
    use_sites = flex.bool()

    for resi in chain_obj.conformers()[0].residues():

        if (iotbx.pdb.common_residue_names_get_class(name=resi.resname) != "common_amino_acid"):
            continue

        single = iotbx.pdb.amino_acid_codes.one_letter_given_three_letter[resi.resname]
        xyz = (0,0,0)
        use = False

        for atom in resi.atoms():

            if (atom.name.strip() == "CA") and (atom.element.strip() == "C"):
              xyz = atom.xyz
              use = True
              break

        seq.append(single)
        sites.append(xyz)
        use_sites.append(use)

    return ("".join(seq), sites, use_sites)

def align_structures_rigid(
    mov_hierarchy, 
    ref_hierarchy,
    ):
    """Extract c-alpha sites from the structures and align"""

    lsq_rt, alignment_sites, reference_sites = align_chains_rigid(
        mov_chain = protein(mov_hierarchy, copy=True).models()[0].only_chain(),
        ref_chain = protein(ref_hierarchy, copy=True).models()[0].only_chain(),
    )

    return GlobalAlignment(
        alignment_mx = lsq_rt,
        alignment_sites = alignment_sites,
        reference_sites = reference_sites,
    )

def align_chains_rigid(mov_chain, ref_chain):
    """Takes two chains and aligns them - return rt_mx"""

    import mmtbx.alignment
    from scitbx.math import superpose

    mov_seq, mov_sites, mov_flags = extract_sites_for_alignment(mov_chain)
    ref_seq, ref_sites, ref_flags = extract_sites_for_alignment(ref_chain)

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
    total = len(alignment.a) - alignment.a.count("-")

    # Print
    out_catch = lg.ListStream()
    alignment.pretty_print(
        matches = matches,
        out = out_catch,
        block_size = 50,
        n_block = 1,
        top_name = "fixed",
        bottom_name = "moving",
    )

    logger.info(out_catch)

    # Create list of selected sites
    ref_sites_sel = flex.vec3_double()
    mov_sites_sel = flex.vec3_double()
    for ia,ib,m in zip(alignment.i_seqs_a, alignment.i_seqs_b, matches):
        # Ignore non-matches
        if (m not in ["|", "*"]):
            continue
        # Check that the sites are flagged to be used
        if (ref_flags[ia] and mov_flags[ib]):
            # Append sites to list to align
            ref_sites_sel.append(ref_sites[ia])
            mov_sites_sel.append(mov_sites[ib])

    if (ref_sites_sel.size() == 0):
        raise AlignmentError("No matching C-alpha atoms.")

    lsq_rt = superpose.least_squares_fit(
        reference_sites = ref_sites_sel,
        other_sites = mov_sites_sel
    ).rt()

    return lsq_rt, mov_sites_sel, ref_sites_sel


import os

from scitbx.array_family import flex

import iotbx

def generate_adjacent_symmetry_copies(ref_hierarchy, crystal_symmetry, buffer_thickness=0):
    """Find symmetry copies of the protein in contact with the asu and generate these copies"""

    sym_ops_mat, contacts_out = get_symmetry_operations_to_generate_crystal_contacts(ref_hierarchy=ref_hierarchy,
                                                                                     crystal_symmetry=crystal_symmetry,
                                                                                     buffer_thickness=buffer_thickness)

    sym_hierarchies, chain_mappings = generate_crystal_copies_from_operations(ref_hierarchy=ref_hierarchy,
                                                                              crystal_symmetry=crystal_symmetry,
                                                                              sym_ops_mat=sym_ops_mat)

    return sym_ops_mat, contacts_out, sym_hierarchies, chain_mappings

def get_symmetry_operations_to_generate_crystal_contacts(ref_hierarchy, crystal_symmetry, buffer_thickness):
    """Extract symmetry operations to generate symmetry copies of the protein that form crystal contacts"""

    # Extract the xray structure from the reference hierarchy
    ref_struc = ref_hierarchy.extract_xray_structure(crystal_symmetry=crystal_symmetry)

    # Extract the mappings that will tell us the adjacent symmetry copies
    asu_mappings = ref_struc.asu_mappings(buffer_thickness=buffer_thickness)
    uc = asu_mappings.unit_cell()
    # Symmetry operations for each atom
    mappings = asu_mappings.mappings()

    # There should be one mappings list per atom
    assert len(ref_struc.scatterers()) == len(mappings)

    # Filter out the non-identity transformations
    filtered_mappings = []
    for at, m in zip(ref_struc.scatterers(), mappings):
        for site_map in m:
            rt_mx = asu_mappings.get_rt_mx(site_map)
            if rt_mx.as_xyz() == 'x,y,z':
                assert site_map.mapped_site() == uc.orthogonalize(at.site), 'SYMMETRY OPERATION IS NOT THE IDENTITY'
            else:
                filtered_mappings.append(  ( at,
                                             rt_mx  ) )

    # Different symmetry operations that map to symmetry neighbours
    uniq_sym_op_xyz = sorted(list(set([m[1].as_xyz() for m in filtered_mappings])))

    # Returned list of symmatry contacts, and operations
    sym_ops_out = []
    contacts_out = []

    for sym_op_xyz in uniq_sym_op_xyz:

        # Select atom mappings with this sym op
        contacts, sym_op_maps = zip(*[x for x in filtered_mappings if x[1].as_xyz()==sym_op_xyz])

        uniq_sym_ops = list(set([x.as_xyz() for x in sym_op_maps]))

        assert len(uniq_sym_ops) == 1, 'MORE THAN ONE UNIQUE SYMMETRY OPERATION PRESENT FOR SAME OPERATION?!'
#        print 'NEW OPERATION'
#        print '\n'.join(uniq_sym_ops)
#        print 'SIZE OF CONTACT AREA:', len(sym_op_maps)

        sym_ops_out.append(sym_op_maps[0])
        contacts_out.append(contacts)

    return sym_ops_out, contacts_out

def generate_crystal_copies_from_operations(ref_hierarchy, crystal_symmetry, sym_ops_mat):
    """Take a list of symmetry operations and apply them to reference hierarchy"""

    # Extract the xray structure from the reference hierarchy
    ref_struc = ref_hierarchy.extract_xray_structure(crystal_symmetry=crystal_symmetry)

    # Create a list of chain ids for the new symmetry copies
    # Already existing chain ids
    ref_chains = [c.id for c in ref_hierarchy.chains()]
    # Number of new chains needed - one for each reference chain for each symmetry operation
    num_new_chains = len(ref_chains)*len(sym_ops_mat)
    # Select the new chains, ignoring the ones already in the reference
    new_chain_ids = [c for c in iotbx.pdb.systematic_chain_ids()[0:(num_new_chains+len(ref_chains))] if c not in ref_chains][0:num_new_chains]

    assert not [c for c in new_chain_ids if c in ref_chains], 'GENERATED CHAIN IDS ARE NOT UNIQUE'

    # Create combinations of the different symmetry operations and the reference chains and map them to new chain ids
    sym_op_ref_chain_combinations = []
    for sym_op in sym_ops_mat:
        for r_chain in ref_chains:
            sym_op_ref_chain_combinations.append((sym_op.as_xyz(), r_chain))
    # Dictionary to map symmetry operations and reference chain ids to new chain ids
    new_chain_ids_hash = dict(zip(sym_op_ref_chain_combinations, new_chain_ids))

    # Hierarchies to be returned
    sym_hierarchies = []
    chain_mappings = dict([(c, []) for c in ref_chains])

    for sym_op in sym_ops_mat:

        sym_op_rt = sym_op.as_rational().as_float()

        # Transform all of the coordinates in the reference structure
        transformed_coords = sym_op_rt * ref_struc.sites_frac()

        # Create copy of the xray structure to play with, and set the transformed coordinates
        new_struc = ref_struc.customized_copy()
        new_struc.set_sites_frac(transformed_coords)

        # Transfer the sites to a new hierarchy
        new_hierarchy = ref_hierarchy.deep_copy()
        new_hierarchy.adopt_xray_structure(new_struc)

        # Update the chain ids in the new structure
        for chain in new_hierarchy.chains():
            old_id = chain.id
            chain.id = new_chain_ids_hash[(sym_op.as_xyz(), chain.id)]
#            print 'CHANGE CHAIN ID: {!s} -> {!s}'.format(old_id, chain.id)
            chain_mappings[old_id].append(chain.id)

        # Add the hierarchy to the output dict, referenced by the symmetry operation
        sym_hierarchies.append(new_hierarchy)

    return sym_hierarchies, chain_mappings

def combine_hierarchies(list_of_hierarchies):
    """Combine a list of hierarchies into one hierarchy -- Requires all of the chain identifiers to be unique"""

    top_h = list_of_hierarchies[0].deep_copy()

    for next_h in list_of_hierarchies[1:]:
        top_h.transfer_chains_from_other(next_h.deep_copy())

    return top_h

#if __name__=='__main__':
#
#    from PANDDAs.Main import dataset_handler
#
#    d=dataset_handler(-1, './reference.pdb', './reference.mtz')
#
#    i = d.get_input()
#    h = d.get_hierarchy()
#
#    sym_ops, contact_mappings, sym_hierarchies, chain_mappings = generate_adjacent_symmetry_copies(    ref_hierarchy=h,
#                                                                                                       crystal_symmetry=i.crystal_symmetry(),
#                                                                                                       buffer_thickness=5)
#
#    print 'CHAIN MAPPINGS:'
#    for ch in chain_mappings.keys():
#        print '\tCHAIN {!s} maps to {!s}'.format(ch, chain_mappings[ch])
#    print 'SYMMETRY HEIRARCHIES:'
#    for x in sym_hierarchies:
#        print '\t',x
#
#    combined_sym_hierarchy = combine_hierarchies(sym_hierarchies)
#
#    output_file = './reference-symmetry-contacts.pdb'
#
#    assert not os.path.exists(output_file)
#    combined_sym_hierarchy.write_pdb_file(output_file)
#

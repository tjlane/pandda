import os

import scipy

from scitbx.array_family import flex

import iotbx.pdb
from cctbx import crystal

def generate_adjacent_symmetry_copies(ref_hierarchy, crystal_symmetry, buffer_thickness=0, method=2):
    """Find symmetry copies of the protein in contact with the asu and generate these copies"""

    sym_ops_mat = get_symmetry_operations_to_generate_crystal_contacts( ref_hierarchy=ref_hierarchy,
                                                                        crystal_symmetry=crystal_symmetry,
                                                                        buffer_thickness=buffer_thickness   )

    sym_hierarchies, chain_mappings = generate_crystal_copies_from_operations(ref_hierarchy=ref_hierarchy,
                                                                              crystal_symmetry=crystal_symmetry,
                                                                              sym_ops_mat=sym_ops_mat)

    return sym_ops_mat, sym_hierarchies, chain_mappings

def get_symmetry_operations_to_generate_crystal_contacts(ref_hierarchy, crystal_symmetry, buffer_thickness):
    """Use an alternate method to identify the symmetry operations required to generate crystal contacts"""

    # Extract the xray structure from the reference hierarchy
    ref_struc = ref_hierarchy.extract_xray_structure(crystal_symmetry=crystal_symmetry)
    ref_atoms = ref_hierarchy.atoms()

    # Extract the mappings that will tell us the adjacent symmetry copies
    asu_mappings = ref_struc.asu_mappings(buffer_thickness=buffer_thickness)
    uc = asu_mappings.unit_cell()
    # Symmetry operations for each atom
    mappings = asu_mappings.mappings()

    # There should be one mappings list per atom
    assert len(ref_struc.scatterers()) == len(mappings)

    # Get all atom pairs within distance_cutoff distance
    pair_generator = crystal.neighbors_fast_pair_generator(
        asu_mappings,
        distance_cutoff=buffer_thickness)
    sym_operations = []
    for pair in pair_generator:
      # obtain rt_mx_ji - symmetry operator that should be applied to j-th atom
      # to transfer it to i-th atom
      rt_mx_i = asu_mappings.get_rt_mx_i(pair)
      rt_mx_j = asu_mappings.get_rt_mx_j(pair)
      rt_mx_ji = rt_mx_i.inverse().multiply(rt_mx_j)
      # if it is not a unit matrix, that is symmetry related pair of atoms
      if not rt_mx_ji.is_unit_mx():
        if rt_mx_ji not in sym_operations:
          sym_operations.append(rt_mx_ji)

    return sym_operations

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
            chain_mappings[old_id].append(chain.id)

        # Add the hierarchy to the output dict, referenced by the symmetry operation
        sym_hierarchies.append(new_hierarchy)

    return sym_hierarchies, chain_mappings

def combine_hierarchies(list_of_hierarchies):
    """Combine a list of hierarchies into one hierarchy -- Requires all of the chain identifiers to be unique"""
    top_h = list_of_hierarchies[0].deep_copy()
    for next_h in list_of_hierarchies[1:]: top_h.transfer_chains_from_other(next_h.deep_copy())
    return top_h

if __name__=='__main__':

    input_file = './reference.pdb'

    inp = iotbx.pdb.input(input_file)
    hie = inp.construct_hierarchy()

    for method in [1,2]:
        sym_ops, contact_mappings, sym_hierarchies, chain_mappings = generate_adjacent_symmetry_copies(    ref_hierarchy=hie,
                                                                                                           crystal_symmetry=inp.crystal_symmetry(),
                                                                                                           buffer_thickness=50,
                                                                                                           method=method)

        print 'CHAIN MAPPINGS:'
        for ch in chain_mappings.keys():
            print '\tCHAIN {!s} maps to {!s}'.format(ch, chain_mappings[ch])
        print 'SYMMETRY HEIRARCHIES:'
        for x in sym_hierarchies:
            print '\t',x
        print 'SYMMETRY OPERATIONS:'
        for x in sorted(sym_ops, key=lambda m: str(m)):
            print '\t',x
        print '{!s} SYMMETRY COPIES GENERATED'.format(len(sym_hierarchies))

        combined_sym_hierarchy = combine_hierarchies(sym_hierarchies)

        output_file = input_file.replace('.pdb', '-contacts-method{!s}.pdb'.format(method))

        assert not os.path.exists(output_file)
        combined_sym_hierarchy.write_pdb_file(output_file)


import copy, time

import cctbx.uctbx, cctbx.sgtbx, cctbx.maptbx
import scitbx.sparse

from scitbx.array_family import flex

from giant.structure.align import *
from giant.xray.symmetry import get_crystal_contact_operators

def make_supercell(unit_cell, size=(3,3,3)):
    """Create an enlarged supercell composed of multiple unit cells"""
    assert isinstance(size, int) or (len(size) == 3), 'Size must be either single scale factor or 3 scale factors for (x,y,z)'
    if isinstance(size, int):   scales = [size]*3   + [1]*3
    else:                       scales = list(size) + [1]*3
    assert len(scales) == 6, 'Ooops, something seems to have gone wrong...: size {}, scales {}'.format(size, scales)
    old_params = unit_cell.parameters()
    new_params = [a*b for a,b in zip(old_params, scales)]
    return cctbx.uctbx.unit_cell(new_params)

def calculate_offset_to_centre_grid(grid_dimensions, centre_on):
    return [a-(b/2.0) for a,b in zip(centre_on, grid_dimensions)]

def get_subset_of_grid_points(gridding, grid_indices):
    """Use a set of indices to mask the grid - returns masked grid points"""
    import numpy
    mask_binary = numpy.zeros(gridding.n_grid_points(), dtype=bool)
    mask_binary.put(grid_indices, True)
    grid = flex.grid(gridding.n_real())
    for p in flex.nested_loop(gridding.n_real()):
        if mask_binary[grid(p)]: yield p

def create_native_map(native_crystal_symmetry, native_sites, native_hierarchy, alignment, reference_map, site_mask_radius=6, step=0.7, filename=None):
    """
    Transform the reference-aligned map back to the native crystallographic frame
    native_sites            - defines region that map will be masked around
    native_hierarchy        - full model for the crystal
    native_crystal_symmetry - crystal symmetry
    reference_map           - basic_map object of the map in the reference frame
    alignment               - Alignment object used to map between the reference frame and the native frame
    site_mask_radius        - Define mask radius around native_sites
    step                    - grid sampling step
    """

    start=time.time()

    native_unit_cell = native_crystal_symmetry.unit_cell()
    native_space_group = native_crystal_symmetry.space_group()

    # ===============================================================================>>>
    # Create a supercell containing the protein model at the centre
    # ===============================================================================>>>

    # create supercell in the native frame
    supercell = make_supercell(native_unit_cell)

    sc_gridding = cctbx.maptbx.crystal_gridding(unit_cell=supercell,        step=step)
    uc_gridding = cctbx.maptbx.crystal_gridding(unit_cell=native_unit_cell, step=step)

    # calculate the origin of the supercell (centred on the protein model)
    # adding origin translates "grid frame" to "crystallographic frame"
    origin = calculate_offset_to_centre_grid(grid_dimensions=supercell.parameters()[0:3], centre_on=native_sites.mean())

    # ===============================================================================>>>
    # Create a masked grid around the protein model
    # ===============================================================================>>>

    # sample the map points near to the protein (transform the structure to be centre of the grid)
    sites_shifted = native_sites - origin
    masked_points_indices = cctbx.maptbx.grid_indices_around_sites(unit_cell=supercell,
                                fft_n_real=sc_gridding.n_real(), fft_m_real=sc_gridding.n_real(),
                                sites_cart=sites_shifted, site_radii=flex.double(sites_shifted.size(),site_mask_radius))
    masked_points_grid_iter = get_subset_of_grid_points(gridding=sc_gridding, grid_indices=masked_points_indices)
    g2c = cctbx.maptbx.grid2cart(sc_gridding.n_real(), supercell.orthogonalization_matrix())
    masked_points_cart = flex.vec3_double(map(g2c, masked_points_grid_iter)) + origin

    from bamboo.pymol_utils.shapes import Sphere
    points = ['from pymol import cmd','from pymol.cgo import *']
    for i,p in enumerate(masked_points_cart):
        if i%100==0: points.append(Sphere(p, 0.2).as_cmd('steve'))
    points = '\n'.join(points)
    with open(filename+'.pml.py', 'w') as fh: fh.write(points)

    # ===============================================================================>>>
    # Sample masked points from the reference-aligned map
    # ===============================================================================>>>

    # transform points to the reference frame
    masked_points_transformed = alignment.nat2ref(masked_points_cart)
    # Sample the map at these points
    masked_values = reference_map.get_cart_values(masked_points_transformed)

    # ===============================================================================>>>
    # Create a native-aligned map and populate with the masked values
    # ===============================================================================>>>

    # Create a map of the density
    sc_map_data = scitbx.sparse.vector(sc_gridding.n_grid_points(), dict(zip(masked_points_indices, masked_values))).as_dense_vector()
    sc_map_data.reshape(flex.grid(sc_gridding.n_real()))
    # Transform the points back to the native frame (simple origin shift)
    sc_map_data = cctbx.maptbx.rotate_translate_map(unit_cell          = supercell,
                                                    map_data           = sc_map_data,
                                                    rotation_matrix    = scitbx.matrix.rec([1,0,0,0,1,0,0,0,1], (3,3)).elems,
                                                    translation_vector = scitbx.matrix.rec([-a for a in origin], (3,1)).elems    )

    # ===============================================================================>>>
    # Select the first (on origin) unit cell of the supercell
    # ===============================================================================>>>

    # Get the indices for the first unit cell
    supercell_grid = flex.grid(sc_gridding.n_real())
    supercell_mask = map(supercell_grid, flex.nested_loop(uc_gridding.n_real()))
    # Extract the map data for those values (and reshape to the right size of the unit cell)
    uc_map_data = sc_map_data.select(supercell_mask)
    uc_map_data.reshape(flex.grid(uc_gridding.n_real()))

    # Create a copy to be the output map
    combined_uc_map_data = copy.copy(uc_map_data)

    # ===============================================================================>>>
    # Apply symmetry operations to generate whole unit cell
    # ===============================================================================>>>

    # Get the symmetry operations for adjacent crystal copies
    sym_ops = get_crystal_contact_operators(hierarchy=native_hierarchy,
                                            crystal_symmetry=native_crystal_symmetry,
                                            distance_cutoff=5.0)

    for sym_op in sym_ops:
        # Get the transformation matrix
        rt_mx = sym_op.as_rational().as_float()
        # Translate the map
        rt_map_data = cctbx.maptbx.rotate_translate_map(unit_cell          = supercell,
                                                        map_data           = uc_map_data,
                                                        rotation_matrix    = rt_mx.r.elems,
                                                        translation_vector = rt_mx.t.elems    )
        # Set any values that are filled in combined_uc_map_data to 0
        rt_map_data.set_selected(combined_uc_map_data!=0, 0)
        # Add values to combined_uc_map_data
        combined_uc_map_data = combined_uc_map_data + rt_map_data

    # ===============================================================================>>>
    # Write output maps
    # ===============================================================================>>>

    if filename is not None:
        iotbx.ccp4_map.write_ccp4_map(  file_name   = filename,
                                        unit_cell   = native_unit_cell,
                                        space_group = native_space_group,
                                        map_data    = combined_uc_map_data,
                                        labels      = flex.std_string(['Map from pandda'])     )
#        iotbx.ccp4_map.write_ccp4_map(  file_name   = filename.replace('.ccp4','.supercell.ccp4'),
#                                        unit_cell   = supercell,
#                                        space_group = cctbx.sgtbx.space_group('P1'),
#                                        map_data    = sc_map_data,
#                                        labels      = flex.std_string(['Map from pandda'])     )

    return uc_map_data


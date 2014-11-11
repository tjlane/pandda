import mmtbx.f_model
from cctbx import maptbx
from libtbx.math_utils import ifloor, iceil
from scitbx.array_family import flex
from numpy import math

def calculate_sampling_distance(resolution, res_factor):
    return resolution*res_factor

def get_bounding_box_for_structure(structure):
    return (structure.sites_cart().min(), structure.sites_cart().max())

def calculate_grid_size(min_carts, max_carts, grid_spacing):
    """Calculate the number of points to be sampled for a box size and sampling distance. Returns the size of the box in A, and the number of points to be sampled along each axis."""

    # Extent of the box (dX, dY, dZ)
    cart_size = tuple([max_c-min_c for min_c, max_c in zip(min_carts, max_carts)])

    # Calculate the number of grid points from the grid sampling
    grid_size = tuple([iceil((1.0*c_size)/grid_spacing)+1 for c_size in cart_size])

    return cart_size, grid_size

def create_cartesian_grid(min_carts, max_carts, grid_spacing):
    """Create the grid of cartesian sites to be sampled based on a spacing between sampling points (in A)"""

    # Calculate the size of the grid (cartesian) and the number of points to be sampled along each axis
    box_size, grid_size = calculate_grid_size(min_carts, max_carts, grid_spacing)

    # Re-calculate the grid spacing ?

    # Grid point volume
    grid_point_volume = grid_spacing**3

    # Calculate the grid points in cartesian space
    grid_points = flex.vec3_double([tuple([min_carts[i]+grid_spacing*grid_point[i] for i in [0,1,2]]) for grid_point in flex.nested_loop(grid_size)])

    return box_size, grid_size, grid_points

def get_grid_points_within_distance_cutoff_of_origin(grid_spacing, distance_cutoff):
    """Find all points on isotropic grid within distance_cutoff of the origin"""

    # Normalise the cutoff to make the grid spacing grid-independent
    grid_index_cutoff = (1.0*distance_cutoff)/grid_spacing

    return get_grid_points_within_index_cutoff_of_origin(grid_index_cutoff=grid_index_cutoff)

def get_grid_points_within_index_cutoff_of_origin(grid_index_cutoff):
    """Find all points relative to the origin within a certain number of grid points"""

    # Round the max grid index down to the nearest int - outer bound on the x,y,z coords
    outer_bound_box = int(grid_index_cutoff)
    # Calculate r/sqrt(3) - inner bound on x,y,z coords
    inner_bound_box = grid_index_cutoff/math.sqrt(3)
    # Calculate r^2 - limiting sphere
    rad_sq = grid_index_cutoff**2
    # List of allowed grid indices
    grid_indices = []

    for x,y,z in flex.nested_loop(begin=(-outer_bound_box, -outer_bound_box, -outer_bound_box),end=(outer_bound_box+1, outer_bound_box+1, outer_bound_box+1)):
        if (abs(x) <= inner_bound_box) and (abs(y) <= inner_bound_box) and (abs(z) <= inner_bound_box):
            grid_indices.append((x, y, z))
        elif (x**2 + y**2 + z**2) <= rad_sq:
            grid_indices.append((x, y, z))

    return grid_indices

def combine_grid_point_and_grid_vectors(start_point, grid_vectors):
    """Take a list of grid vectors and add them to a particular grid point"""

    start_point = flex.int(start_point)
    grid_points = [tuple(start_point+flex.int(vec)) for vec in grid_vectors]
    return grid_points



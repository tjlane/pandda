from numpy import math
from scitbx.array_family import flex
from libtbx.math_utils import ifloor, iceil

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

def get_grid_points_within_distance_cutoff_of_cart_sites(cart_sites, grid_spacing, distance_cutoff):
    """Find all points on isotropic grid within distance cutoff of cartesian sites"""

    # Normalise the cutoff and the sites to make them grid-independent
    grid_index_cutoff = (1.0*distance_cutoff)/grid_spacing
    grid_sites = [tuple([1.0*c/grid_spacing for c in coords]) for coords in cart_sites]

    return get_grid_points_within_index_cutoff_of_grid_sites(grid_sites=grid_sites, grid_index_cutoff=grid_index_cutoff)

def get_grid_points_within_index_cutoff_of_grid_sites(grid_sites, grid_index_cutoff):
    """Find all points on a grid within a certain number of grid points of grid sites (not necessarily integer sites)"""

    # Calculate the size of the grid we need to check over
    min_x = ifloor(min([s[0] for s in grid_sites]) - grid_index_cutoff)
    if min_x < 0: min_x=0
    min_y = ifloor(min([s[1] for s in grid_sites]) - grid_index_cutoff)
    if min_y < 0: min_y=0
    min_z = ifloor(min([s[2] for s in grid_sites]) - grid_index_cutoff)
    if min_z < 0: min_z=0
    max_x = iceil(max([s[0] for s in grid_sites]) + grid_index_cutoff)
    max_y = iceil(max([s[1] for s in grid_sites]) + grid_index_cutoff)
    max_z = iceil(max([s[2] for s in grid_sites]) + grid_index_cutoff)
    # Grid Extremities
    min_grid = (min_x, min_y, min_z)
    max_grid = (max_x+1, max_y+1, max_z+1)

    # Round the max grid index down to the nearest int - outer bound on the dx,dy,dz values
    outer_bound_box = int(grid_index_cutoff)
    # Calculate r/sqrt(3) - inner bound on dx, dy, dz values
    inner_bound_box = grid_index_cutoff/math.sqrt(3)
    # Calculate r^2 - limiting sphere
    rad_sq = grid_index_cutoff**2

    # List of allowed grid indices
    grid_indices = []

    for gp in flex.nested_loop(min_grid, max_grid):
        for site in grid_sites:
            dx, dy, dz = [abs(p1-p2) for p1,p2 in zip(gp, site)]

            if (dx > outer_bound_box) or (dy > outer_bound_box) or (dz > outer_bound_box):
                continue
            elif (dx <= inner_bound_box) and (dy <= inner_bound_box) and (dz <= inner_bound_box):
                grid_indices.append(gp)
                break
            elif (dx**2 + dy**2 + dz**2) <= rad_sq:
                grid_indices.append(gp)
                break

    return grid_indices




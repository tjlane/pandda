from scipy import spatial
import numpy
from scitbx.array_family import flex
from libtbx.math_utils import ifloor, iceil

class grid_partition(object):
    def __init__(self, grid_size, grid_spacing, atomic_hierarchy):
        """Partition a grid based on the nearest neighbour calpha atom for each grid site"""

        # Save inputs
        self.hierarchy = atomic_hierarchy
        self.grid_size = grid_size
        self.grid_spacing = grid_spacing
        self.grid_indexer = flex.grid(grid_size)

        # Calculate partition variables
        atoms = [at for at in atomic_hierarchy.atoms_with_labels()]
        self.grid_sites = list(flex.nested_loop(grid_size))
        self.atom_sites_grid = numpy.array([a.xyz for a in atoms])/grid_spacing
        self.atom_sites_cart = numpy.array([a.xyz for a in atoms])
        self.atom_labels = [(a.chain_id, a.resid()) for a in atoms]

        # Distances from grid sites to nearest atom site
        self.nn_dists = None
        # Index of nearest atom to grid sites
        self.nn_groups = None
        # Label of nearest atom to grid sites
        self.nn_atom_labels = None

    def partition_grid(self):
        """Find the nearest neighbour for each grid point"""
        tree = spatial.KDTree(data=self.atom_sites_grid)
        self.nn_dists, self.nn_groups = tree.query(self.grid_sites)
        self.nn_atom_labels = [self.atom_labels[i] for i in self.nn_groups]
    def query_by_grid_indices(self, idxs):
        """Return the atom label for a grid site index"""
        if not self.nn_atom_labels: self.partition_grid()
        return [self.nn_atom_labels[i] for i in idxs]
    def query_by_grid_points(self, gps):
        """Return the atom label for a grid point"""
        if not self.nn_atom_labels: self.partition_grid()
        return [self.nn_atom_labels[self.grid_indexer(g)] for g in gps]

    def query_by_cart_points(self, sites_cart):
        """Dynamically calculate the nearest atom site to the input points"""
        tree = spatial.KDTree(data=self.atom_sites_cart)
        nn_dists, nn_groups = tree.query(self.grid_sites)
        return [self.atom_labels[i] for i in nn_groups]

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
    grid_points_cart = flex.vec3_double([tuple([min_carts[i]+grid_spacing*grid_point[i] for i in [0,1,2]]) for grid_point in flex.nested_loop(grid_size)])

    return box_size, grid_size, grid_points_cart

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
    inner_bound_box = grid_index_cutoff/numpy.math.sqrt(3)
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

def get_grid_points_within_distance_cutoff_of_cart_sites(cart_sites, grid_spacing, max_dist, min_dist=None):
    """Find all points on isotropic grid within distance cutoff of cartesian sites"""

    # Normalise the cutoff and the sites to make them grid-independent
    max_grid_dist = (1.0*max_dist)/grid_spacing
    if min_dist: min_grid_dist = (1.0*min_dist)/grid_spacing
    else:        min_grid_dist = None
    grid_sites = [tuple([1.0*c/grid_spacing for c in coords]) for coords in cart_sites]

    return get_grid_points_within_index_cutoff_of_grid_sites(grid_sites=grid_sites, max_grid_dist=max_grid_dist, min_grid_dist=min_grid_dist)

def get_grid_points_within_distance_cutoff_of_cart_sites_2(cart_sites, grid_spacing, max_dist, min_dist=None):
    """Find all points on isotropic grid within distance cutoff of cartesian sites"""

    # Normalise the cutoff and the sites to make them grid-independent
    max_grid_dist = (1.0*max_dist)/grid_spacing
    if min_dist: min_grid_dist = (1.0*min_dist)/grid_spacing
    else:        min_grid_dist = None
    grid_sites = [tuple([1.0*c/grid_spacing for c in coords]) for coords in cart_sites]

    return get_grid_points_within_index_cutoff_of_grid_sites_2(grid_sites=grid_sites, max_grid_dist=max_grid_dist, min_grid_dist=min_grid_dist)

def get_grid_points_within_index_cutoff_of_grid_sites(grid_sites, max_grid_dist, min_grid_dist=None):
    """Find all points on a grid within a certain number of grid points of grid sites (not necessarily integer sites)"""

    # Calculate the size of the grid we need to check over
    min_x = ifloor(min([s[0] for s in grid_sites]) - max_grid_dist)
    if min_x < 0: min_x=0
    min_y = ifloor(min([s[1] for s in grid_sites]) - max_grid_dist)
    if min_y < 0: min_y=0
    min_z = ifloor(min([s[2] for s in grid_sites]) - max_grid_dist)
    if min_z < 0: min_z=0
    max_x = iceil(max([s[0] for s in grid_sites]) + max_grid_dist)
    max_y = iceil(max([s[1] for s in grid_sites]) + max_grid_dist)
    max_z = iceil(max([s[2] for s in grid_sites]) + max_grid_dist)
    # Grid Extremities
    min_grid = (min_x, min_y, min_z)
    max_grid = (max_x+1, max_y+1, max_z+1)

    # Round the max grid distance down to the nearest int - outer bound on the dx,dy,dz values
    outer_bound_box = int(max_grid_dist)
    # Calculate r/sqrt(3) - inner bound on dx, dy, dz values
    inner_bound_box = max_grid_dist/numpy.math.sqrt(3)
    # Calculate r^2 - limiting sphere
    rad_sq = max_grid_dist**2

    # List of allowed grid indices
    outer_indices = []
    inner_indices = []

    # Iterate through and add valid points
    for gp in flex.nested_loop(min_grid, max_grid):
        for site in grid_sites:
            dx, dy, dz = [abs(p1-p2) for p1,p2 in zip(gp, site)]

            if (dx > outer_bound_box) or (dy > outer_bound_box) or (dz > outer_bound_box):
                continue
            elif (dx <= inner_bound_box) and (dy <= inner_bound_box) and (dz <= inner_bound_box):
                outer_indices.append(gp)
                break
            elif (dx**2 + dy**2 + dz**2) <= rad_sq:
                outer_indices.append(gp)
                break

    # Filter the grid points that are too close to the protein
    if min_grid_dist:
        # Round the min grid distance up to the nearest int - outer bound on the dx,dy,dz values
        outer_bound_box = int(min_grid_dist) + 1
        # Calculate r/sqrt(3) - inner bound on dx, dy, dz values
        inner_bound_box = min_grid_dist/numpy.math.sqrt(3)
        # Calculate r^2 - limiting sphere
        rad_sq = min_grid_dist**2

        # Iterate through and add valid points
        for gp in outer_indices:
            for site in grid_sites:
                dx, dy, dz = [abs(p1-p2) for p1,p2 in zip(gp, site)]

                if (dx > outer_bound_box) or (dy > outer_bound_box) or (dz > outer_bound_box):
                    continue
                elif (dx <= inner_bound_box) and (dy <= inner_bound_box) and (dz <= inner_bound_box):
                    inner_indices.append(gp)
                    break
                elif (dx**2 + dy**2 + dz**2) <= rad_sq:
                    inner_indices.append(gp)
                    break

    if min_grid_dist:
        total_indices = [gp for gp in outer_indices if gp not in inner_indices]
    else:
        total_indices = outer_indices

    return total_indices, outer_indices, inner_indices

def get_grid_points_within_index_cutoff_of_grid_sites_2(grid_sites, max_grid_dist, min_grid_dist=None):
    """Find all points on a grid within a certain number of grid points of grid sites (not necessarily integer sites)"""

    # Find the size of the grid that we'll need
    max_x = iceil(max([s[0] for s in grid_sites]) + max_grid_dist)
    max_y = iceil(max([s[1] for s in grid_sites]) + max_grid_dist)
    max_z = iceil(max([s[2] for s in grid_sites]) + max_grid_dist)
    max_grid = (max_x+2, max_y+2, max_z+2)
    print 'MAX GRID:', max_grid

    # Round the grid sites to the nearest grid point
    int_grid_sites = [tuple(map(int, site)) for site in grid_sites]

    # Grid objects
    grid_indexer = flex.grid(max_grid)
    grid_size = flex.product(flex.int(max_grid))

    # Masks
    outer_mask = numpy.zeros(grid_size, dtype=int)
    inner_mask = numpy.zeros(grid_size, dtype=int)

    # Find all of the grid vectors within max_dist of grid sites
    outer_grid_vectors = get_grid_points_within_index_cutoff_of_origin(grid_index_cutoff=max_grid_dist)
    outer_indices = []
    for site in int_grid_sites:
        outer_indices.extend(combine_grid_point_and_grid_vectors(site, grid_vectors=outer_grid_vectors))
    # Iterate through the grid and create binary masks
    print 'OUTER LENGTH:', len(outer_indices)
    a=[gp for gp in outer_indices]
    print 'MAX GP:', max(a)
    print 'MAX X:', max([x[0] for x in a])
    print 'MAX Y:', max([x[1] for x in a])
    print 'MAX Z:', max([x[2] for x in a])
    a=[grid_indexer(gp) for gp in outer_indices]
    print 'MAX INDEX:', max(a)
    [outer_mask.put(i, 1) for i,gp in zip(a,outer_indices)]

    if min_grid_dist:
        # Find all of the grid vectors within min_dist of grid sites
        inner_grid_vectors = get_grid_points_within_index_cutoff_of_origin(grid_index_cutoff=min_grid_dist)
        inner_indices = []
        for site in int_grid_sites:
            inner_indices.extend(combine_grid_point_and_grid_vectors(site, grid_vectors=inner_grid_vectors))
        # Iterate through the grid and create binary masks
        [inner_mask.put(grid_indexer(gp), 1) for gp in inner_indices]

    # Convert from the mask back to grid points
    outer_indices = [gp for i, gp in enumerate(flex.nested_loop(max_grid)) if outer_mask[i]==1]
    inner_indices = [gp for i, gp in enumerate(flex.nested_loop(max_grid)) if inner_mask[i]==1]
    total_indices = [gp for i, gp in enumerate(flex.nested_loop(max_grid)) if (inner_mask[i]==0 and outer_mask[i]==1)]

    return total_indices, outer_indices, inner_indices




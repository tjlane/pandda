import giant.logs as lg
logger = lg.getLogger(__name__)

import itertools, collections
import math, numpy
from scitbx.matrix import sym, diag
from scitbx.linalg import eigensystem_real_symmetric
from libtbx import adopt_init_args, group_args
from pandemic.adp import constants

from pandemic.adp.echt.analysis.cluster_tls_groups.utils import \
    make_selection_strings, tls_group_comparison_matrix, tls_group_adjacency_matrix
from pandemic.adp.echt.analysis.cluster_tls_groups.pymol import base_pymol_script


def uij_overlap_mass(a,b):
    """Calculate the common overlap between two symmetric tensors"""
    # Difference between two inputs
    d = tuple([av-bv for av,bv in zip(a,b)])
    # Diagonalise
    es = eigensystem_real_symmetric(d)
    vals = es.values()
    vecs = es.vectors().as_scitbx_matrix()
    # Set negative eigenvalues to zero and back transform
    vals.set_selected(vals<0.0, 0.0)
    vals = diag(diag_elems=vals)
    d_mod = vecs.transpose() * vals * vecs
    # Subtract difference from input to give overlap
    overlap_sym = tuple([av-dv for av,dv in zip(a,d_mod.as_sym_mat3())])
    overlap_iso = numpy.mean(overlap_sym[:3])
    return overlap_iso

def uij_overlap_mass_function(uijs_a, uijs_b):
    """Calculate common overlaps between two arrays of symmetric tensors"""
    u_dists = numpy.array([uij_overlap_mass(a,b) for a,b in zip(uijs_a,uijs_b)])
    return u_dists

def uij_overlap_mass_function_weighted_average(uijs_a, uijs_b):
    """Calculate weighted average common overlaps between two arrays of symmetric tensors (weights calculated as overlap/u_input)"""
    overlaps = uij_overlap_mass_function(uijs_a, uijs_b)
    if overlaps.sum() == 0.0: return 0.0
    u_sizes = numpy.array([numpy.mean(u[0:3]) for u in uijs_a])
    u_sizes[u_sizes == 0.0] = 1e-6 # avoid zerodivision error
    weights = overlaps / u_sizes
    o = numpy.average(overlaps, weights=weights)
    return o

def uij_overlap_mass_function_simple_average(uijs_a, uijs_b):
    """Calculate simple average common overlaps between two arrays of symmetric tensors"""
    overlaps = uij_overlap_mass_function(uijs_a, uijs_b)
    o = numpy.mean(overlaps)
    return o

def make_pymol_script(
    tls_groups,
    connectivity,
    similarity_matrix,
    indices_sets,
    model_structure,
    output_filename,
    group_colours,
    ):
    """Make pymol script to show results of overlap analysis"""

    from giant.pymol_utils import shapes

    minimum = 0.
    maximum = 1e-16 + similarity_matrix[~numpy.eye(similarity_matrix.shape[0], dtype=bool)].max()
    midpoint = (minimum + maximum) / 2.

    def map_value_to_colour(v, minimum=minimum, midpoint=midpoint, maximum=maximum):
        r = min(
            (maximum - v)/(maximum - midpoint),
            1.0,
            )
        g = min(
            (v - minimum)/(midpoint - minimum),
            1.0,
            )
        b = min(
            (maximum - v)/(maximum - midpoint),
            (v - minimum)/(midpoint - minimum),
            )
        return (r,g,b)

    def map_value_to_radius(v, scale=0.5/maximum):
        r = scale*v
        return r

    sc = base_pymol_script(
        tls_groups = tls_groups,
        connectivity_matrix = connectivity,
        model_structure = model_structure,
        output_filename = output_filename,
        )

    sc.add_comparison_matrix(
        comparison_matrix = similarity_matrix,
        colour_function = map_value_to_colour,
        radius_function = map_value_to_radius,
        obj_name = 'tls-group-overlaps',
        min_max_values = (minimum, maximum),
        )

    s = sc.script

    ##################################################

    # Iterate through different input levels
    group_counter = 1
    for i_l, groups_of_indices in enumerate(indices_sets):

        group_shapes = []

        for indices in groups_of_indices:
            indices = sorted(indices)

            group_colour = next(group_colours)[:3]

            for i in indices:
                for j in indices:
                    if i==j:
                        break
                    if not connectivity[i,j]:
                        continue

                    cyl = shapes.Cylinder(
                        start=sc.centroids[i],
                        end=sc.centroids[j],
                        colors=[group_colour, group_colour],
                        radius=map_value_to_radius(maximum),
                    )
                    group_shapes.append(cyl)

        s.add_shapes(group_shapes, obj='new_groups', state=i_l+1)
    #s.disable(obj='new_groups')

    sc.write()

    return


class ClusterTLSGroups_OverlapMass(object):


    _comparison_functions = {
        'simple_average' : uij_overlap_mass_function_simple_average,
        'weighted_average' : uij_overlap_mass_function_weighted_average,
    }

    def __init__(self,
        xyz_cutoff_distance = 4.0,
        comparison = 'weighted_average',
        ):
        adopt_init_args(self, locals())

        self.metric_type = 'similarity'

        if self.comparison not in self._comparison_functions.keys():
            raise Sorry(
                'Invalid comparison function ({}). Must be one of: {}'.format(
                    self.comparison,
                    ', '.join(self._comparison_functions.keys()),
                    )
                )

        from pandemic.adp.echt.analysis.cluster_tls_groups.identify_modules import IdentifyModules
        self.identify_modules = IdentifyModules(
            comparison_metric = 'similarity',
            )

    def __call__(self,
        tls_groups,
        ):

        uij_overlap_function = self._comparison_functions[self.comparison]

        tls_overlap = constants.EIGHTPISQ * tls_group_comparison_matrix(
            tls_groups = tls_groups,
            comparison_function = uij_overlap_function,
            make_symmetric = None,
            )

        xyz_distance = tls_group_adjacency_matrix(
            tls_groups = tls_groups,
            )

        connectivity = (xyz_distance < self.xyz_cutoff_distance)

        module_output = self.identify_modules(
            connectivity = connectivity,
            comparison_matrix = tls_overlap,
            )

        result = group_args(
            comparison_matrix = tls_overlap,
            connectivity = connectivity,
            module_info = module_output,
            cumulative_hierarchy = module_output.cumulative_hierarchy,
            reduced_hierarchy = module_output.reduced_hierarchy,
            )

        logger(self.summary(result))

        return result

    def description(self):
        text = """
        > Clustering Parameters
        Similarity metric used: Overlap Mass ({comparison})
        Connectivity cutoff: {xyz_threshold}

        > Clustering Metric Description
        The overlap between TLS Group 1 and TLS Group 2 is the amount of B-factor of TLS Group 2 that is reproduced by using the TLS matrices from TLS Group 1 on the coordinates of TLS Group 2.
        This can be thought of as the amount of B-factor of TLS Group 2 that would be the same if we replaced TLS Group 2 with TLS Group 1.
        NOTE! This metric is not symmetric: The overlap of Group1 -> Group2 is not the same as Group2 -> Group1.
        Using this metric, clustering is at different thresholds to identify groups that could be merged into a larger group.
        During clustering, connections are only considered to neighbouring groups with {xyz_threshold} Angstrom.""".format(
            comparison = self.comparison,
            xyz_threshold = self.xyz_cutoff_distance,
            )
        return text

    def summary(self, result):

        r = result

        off_diagonal_values = r.comparison_matrix[~numpy.eye(r.comparison_matrix.shape[0], dtype=bool)]
        thresholds = list(r.module_info.threshold_unique_modules.keys())

        s = ""
        s += "Overlap between groups: {:.3f} - {:.3f} A^2.\n".format(off_diagonal_values.min(), off_diagonal_values.max())
        s += "\n"
        if thresholds:
            s += "Modules identified at thresholds from {:.3f} - {:.3f} A^2.\n".format(min(thresholds),max(thresholds))
            s += "\n"
            s += "Thresholds & Modules:"
            for thresh, modules in r.module_info.threshold_unique_modules.items():
                s += "\n\n"
                s += "> Threshold: {}\n".format(thresh)
                s += "  New potential groupings: {}".format(", ".join(['+'.join(map(str,[i+1 for i in sorted(m)])) for m in modules]))
        else:
            s += "No modules identified"

        return s

    def write_output(self,
        result,
        tls_groups,
        model_structure,
        output_prefix,
        output_directory = './',
        ):

        output_files = collections.OrderedDict()

        y_width = self.identify_modules.threshold_delta

        ###############################################

        # What is the range in the thresholds?
        if result.module_info.threshold_unique_hierarchies:
            y_lim = (
                min(result.module_info.threshold_unique_hierarchies.keys()),
                max(result.module_info.threshold_unique_hierarchies.keys()),
                )
        else:
            y_lim = None

        # Sort hierarchies in ascending order
        threshold_unique_hierarchies = sorted(result.module_info.threshold_unique_hierarchies.items())

        ###############################################

        filename = str(output_directory / (output_prefix+'-new_levels.png'))
        output_files['modules_png'] = filename

        from pandemic.adp.echt.analysis.cluster_tls_groups.plots import LevelPlotAccumulator

        # Initialise plot object
        output_plot = LevelPlotAccumulator(
            n_groups = len(tls_groups),
            title = "New TLS groupings generated\nfrom {} input TLS groups".format(len(tls_groups)),
            x_label = 'Existing TLS Groups',
            y_label = 'B-factor overlap threshold ($\\AA^2$)',
            y_width = y_width,
            y_lim = y_lim,
            )

        # Add levels to plot
        last_threshold = None
        for i_block, (threshold, unique_hierarchy) in enumerate(threshold_unique_hierarchies):
            #if i_block > 0:
                #output_plot.add_line(y_value = threshold - y_width/2. , linestyle='--', color='r')
            #for i, indices_groups in enumerate(unique_hierarchy):
            label = str(threshold)
            output_plot.add_level(
                list_of_lists_of_indices_tuples=unique_hierarchy, #sorted(map(sorted,unique_hierarchy)),
                y_value = threshold,
                label = label,
                )
            last_threshold = threshold
        output_plot.simplify_yaxis_labels()
        output_plot.write(filename)

        ###############################################

        filename = str(output_directory / (output_prefix+'-pymol.py'))
        output_files['pymol_script'] = filename
        make_pymol_script(
            tls_groups = tls_groups,
            connectivity = result.connectivity,
            similarity_matrix = result.comparison_matrix,
            indices_sets = result.cumulative_hierarchy,
            output_filename = filename,
            model_structure = model_structure,
            group_colours = itertools.cycle(output_plot.used_colours),
            )

        ###############################################

        return output_files

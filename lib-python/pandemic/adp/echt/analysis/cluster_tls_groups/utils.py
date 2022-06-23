import numpy
from pandemic.adp.echt.tls import MultiDatasetTLSGroup
from scitbx.array_family import flex

def make_matrix_symmetric_min(m):
    return numpy.min([m, m.T], axis=0)
def make_matrix_symmetric_mean(m):
    return numpy.mean([m, m.T], axis=0)
def make_matrix_symmetric_max(m):
    return numpy.max([m, m.T], axis=0)

make_symmetric_hash = {
    'min'  : make_matrix_symmetric_min,
    'mean' : make_matrix_symmetric_mean,
    'max'  : make_matrix_symmetric_max,
}

def tls_group_comparison_matrix(tls_groups, comparison_function, make_symmetric=None, use_shortcut_for_zero_value_elements=True):
    """Returns n x n matrix where [i,j] represents comparison of how parameters of i reproduce values of j"""

    if make_symmetric is not None: 
        assert make_symmetric in make_symmetric_hash, 'make_symmetric={} is not implemented'.format(make_symmetric)
        make_symmetric_function = make_symmetric_hash[make_symmetric]
    else: 
        make_symmetric_function = None

    n_groups = len(tls_groups)

    tls_dist = numpy.zeros((n_groups, n_groups))

    # What happens if input is zero?
    if use_shortcut_for_zero_value_elements is True:
        zero_val = comparison_function(
            flex.sym_mat3_double([(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)]), 
            flex.sym_mat3_double([(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)]),
            )

    for ii, ig in enumerate(tls_groups):
        iu = ig.uijs().as_1d()
        # shortcut for zero-value elements
        if use_shortcut_for_zero_value_elements and iu.as_double().all_approx_equal(0.0):
            tls_dist[:,ii] == zero_val
            continue
        ic = ig.coordinates.as_1d()
        for jj, jg in enumerate(tls_groups):
            #if ii == jj: continue
            # Extract uijs for coordinates from i using parameters from j
            jg_mix = MultiDatasetTLSGroup(
                index=0, label='',
                tls_parameters = jg.tls_parameters, # original tls parameters
                coordinates = ig.coordinates,       # !!! use NEW coordinates !!!
                origins = jg.origins,               # original origins
                )
            ju = jg_mix.uijs().as_1d()
            # shortcut for zero-value elements
            if use_shortcut_for_zero_value_elements and ju.as_double().all_approx_equal(0.0):
                tls_dist[jj,ii] == zero_val
                continue
            # Differences between the uij values
            u_dist = comparison_function(iu, ju)
            # How the parameters of **j** reproduce **i**
            tls_dist[jj,ii] = u_dist

    if make_symmetric is not None: 
        # tls_dist is not symmetric yet: choose the minimum distance between the group
        tls_dist = make_symmetric_function(tls_dist) 

    return tls_dist

def tls_group_adjacency_matrix(tls_groups):

    n_groups = len(tls_groups)

    xyz_dist = numpy.zeros((n_groups, n_groups))

    for ii, ig in enumerate(tls_groups):
        iu = ig.uijs().as_1d()
        ic = ig.coordinates.as_1d()
        for jj, jg in enumerate(tls_groups):
            if ii >= jj: continue
            jc = jg.coordinates.as_1d()
            x_dist = ic.min_distance_between_any_pair(jc)
            xyz_dist[ii,jj] = x_dist
            xyz_dist[jj,ii] = x_dist

    return xyz_dist

def make_selection_strings(selections, clusters):

    sel_strings = []

    assert len(selections) > numpy.concatenate(clusters).max()

    for cluster_indices in clusters:
        sels = [selections[i] for i in cluster_indices]
        new_sel = '({})'.format(') or ('.join(sels))
        sel_strings.append(new_sel)
    return sel_strings

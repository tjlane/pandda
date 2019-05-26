from pytest import approx, raises

import numpy

tls_matrices = [
    {
        'T':    ( 1.326,     0.581,     0.419,    -0.759,    -0.408,     0.321),
        'L':    ( 6.689,     2.411,     2.637,    -3.229,    -3.728,     1.615),
        'S':    (-0.044,    -0.093,    -0.056,     0.008,    -0.333,    -0.125,    -0.529,     0.467,     0.377),
        }, 
    {
        'T':    ( 0.463,     0.321,     0.868,    -0.025,     0.295,    -0.008),
        'L':    (16.297,     4.092,     4.525,    -7.602,     4.047,    -1.021),
        'S':    (-0.081,     0.114,    -0.340,     0.153,    -0.043,     0.200,     0.371,    -0.147,     0.125),
        }, 
    {
        'T':    ( 0.601,     0.212,     0.270,    -0.005,     0.035,     0.063),
        'L':    ( 3.723,    26.197,    13.389,    -6.258,     3.596,    -3.643),
        'S':    (-0.516,    -0.114,     0.314,    -0.020,     0.379,     0.159,    -0.092,    -0.102,     0.137),
        },
]
tls_matrices_arr = numpy.zeros((3, 21))
for i, vals in enumerate(tls_matrices):
    tls_matrices_arr[i] = vals['T'] + vals['L'] + vals['S']

tls_amplitudes_arr = numpy.array([
    [3., 2.],
    [5., 6.],
    [1., 4.],
])

coordinates = numpy.array([
    [(0.,0.,0.), (1.,0.,0.), (0.,1.,0.), (0.,0.,1.)],
    [(1.,1.,1.), (2.,1.,1.), (1.,2.,1.), (1.,1.,2.)],
])

origins = numpy.array([
    (0.,0.,0.,),
    (1.,1.,1.,),
])

def test_multi_dataset_tls_group():

    from mmtbx.tls.utils import TLSMatricesAndAmplitudesList

    from scitbx.array_family import flex

    tls_matrices_flex = flex.double(tls_matrices_arr)
    assert tls_matrices_flex.all() == tls_matrices_arr.shape

    tls_amplitudes_flex = flex.double(tls_amplitudes_arr)
    assert tls_amplitudes_flex.all() == tls_amplitudes_arr.shape

    tls_parameters = TLSMatricesAndAmplitudesList(
        matrix_values = tls_matrices_flex,
        amplitude_values = tls_amplitudes_flex,
        )

    assert tls_parameters.size() == len(tls_matrices)

    from pandemic.adp.echt.tls import MultiDatasetTLSGroup

    mdg = MultiDatasetTLSGroup(
        index = 33,
        label = 'grouptest',
        tls_parameters = tls_parameters,
        coordinates = coordinates,
        origins = origins,
        )

    assert mdg.index == 33
    assert mdg.label == 'grouptest'

    assert mdg.n_modes == 3
    assert mdg.n_datasets == 2
    assert mdg.n_atoms == 4

    assert mdg.coordinates.nd() == 2
    assert mdg.origins.nd() == 1

    assert mdg.coordinates.all() == (2,4)
    assert mdg.origins.all() == (2,)

    assert numpy.array(mdg.coordinates) == approx(coordinates.reshape((8,3)))
    assert numpy.array(mdg.origins) == approx(origins.reshape((2,3)))

    uijs = mdg.uijs()

    assert uijs.all() == (mdg.n_datasets, mdg.n_atoms)

    assert numpy.array(uijs) == approx(numpy.array(tls_parameters.uijs(
        sites_carts = mdg.coordinates,
        origins = mdg.origins,
        )))

    uijs_by_mode = mdg.uijs_by_mode()

    assert [u.all() for u in uijs_by_mode] == [(mdg.n_datasets, mdg.n_atoms)]*mdg.n_modes

    uijs_by_mode_sum = uijs_by_mode[0]
    for u in uijs_by_mode[1:]:
        uijs_by_mode_sum = uijs_by_mode_sum + u
    assert numpy.array(uijs) == approx(numpy.array(uijs_by_mode_sum))

    uijs_unmultiplied = mdg.uijs_unmultiplied()

    # Check that unmultiplied components equivalent to full output
    for u_mult, (amps, u_unmult) in zip(uijs_by_mode, uijs_unmultiplied):

        assert amps.all() == (mdg.n_datasets,)
        assert [u.all() for u in u_unmult] == [(mdg.n_atoms,)]*mdg.n_datasets

        u_mult_c = [flex.sym_mat3_double(a*u.as_double()) for a, u in zip(amps, u_unmult)]
        assert numpy.array(u_mult).flatten() == approx(numpy.array(u_mult_c).flatten())

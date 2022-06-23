from pytest import approx, raises

def test_uij_overlap_mass():

    from pandemic.adp.echt.analysis.cluster_tls_groups.overlap import uij_overlap_mass, uij_overlap_mass_function, \
        uij_overlap_mass_function_weighted_average, uij_overlap_mass_function_simple_average

    a_b_o = [
        (( 1., 1., 1., 0., 0., 0.), ( 1., 1., 1., 0., 0., 0.), 1.),
        (( 2., 1., 1., 0., 0., 0.), ( 1., 1., 1., 0., 0., 0.), 1.),
        (( 1., 1., 1.,-1., 0., 0.), ( 5., 5., 5., 0., 0., 0.), 1.),
        (( 1., 1., 5., 0., 0., 0.), ( 1., 1., 0., 0., 0., 0.), 2./3.),
        (( 1., 1., 1.,-1., 0., 0.), ( 1., 1., 1., 1., 0., 0.), 1./3.),
        (( 1., 1., 1.,-1., 0., 0.), ( 1., 1., 0., 1., 0., 0.), 0.),
        (( 1., 1., 1., 0., 0., 0.), ( 0., 0., 0., 0., 0., 0.), 0.),
        (( 1., 0., 0., 0., 0., 0.), ( 0., 1., 0., 0., 0., 0.), 0.),
        ]

    for a,b,o in a_b_o:
        assert o == approx(uij_overlap_mass(a,b))

    a_b_o_z = list(zip(*a_b_o))

    assert a_b_o_z[2] == approx(uij_overlap_mass_function(a_b_o_z[0], a_b_o_z[1]).tolist())
    assert a_b_o_z[2] == approx(uij_overlap_mass_function(a_b_o_z[1], a_b_o_z[0]).tolist())

    assert 0.5 == approx(uij_overlap_mass_function_simple_average(a_b_o_z[0], a_b_o_z[1]))
    assert 0.5 == approx(uij_overlap_mass_function_simple_average(a_b_o_z[1], a_b_o_z[0]))

    assert 0.905771495878 == approx(uij_overlap_mass_function_weighted_average(a_b_o_z[0], a_b_o_z[1]))
    assert 0.842767295597 == approx(uij_overlap_mass_function_weighted_average(a_b_o_z[1], a_b_o_z[0]))


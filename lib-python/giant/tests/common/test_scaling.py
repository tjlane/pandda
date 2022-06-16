from pytest import approx
from scitbx.array_family import flex

def test_LinearScaling():

    v0 = 16.0
    v1 = 15.0

    t1_x_values = flex.double(list(range(1000)))

    t1_ref_values = v1*t1_x_values + v0
    t1_mov_values = flex.double(t1_x_values)

    from giant.common.scaling import LinearScaling

    ls = LinearScaling(
        x_values = t1_x_values,
        ref_values = t1_ref_values,
        scl_values = t1_mov_values,
        )

    assert v0 == approx(
        round(float(ls.optimised_values[0]),3)
        )
    assert v1 == approx(
        round(float(ls.optimised_values[1]),3)
        )

    # Check scaling is working correctly
    scaled = ls.transform(t1_mov_values)

    assert 0.0 == approx(
        flex.sum(t1_ref_values-scaled), abs=1e-6,
        )


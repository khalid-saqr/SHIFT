import numpy as np

from shift.mechanics import (
    cauchy_stress,
    conventional_lamb_vector,
    lamb_vector,
    newtonian_extra_stress,
    second_invariant_q,
    vorticity_from_gradient,
    wall_shear_stress,
)


def test_simple_shear_q_vorticity_and_lamb_sign_conventions() -> None:
    shear = 7.5
    gradient = np.zeros((1, 3, 3))
    gradient[0, 0, 1] = shear
    velocity = np.array([[2.0, 0.0, 0.0]])

    q = second_invariant_q(gradient)
    omega = vorticity_from_gradient(gradient)
    ell = lamb_vector(velocity, omega)
    conventional = conventional_lamb_vector(velocity, omega)

    assert np.array_equal(q, np.array([0.0]))
    assert np.array_equal(omega, np.array([[0.0, 0.0, -shear]]))
    assert np.array_equal(ell, np.array([[0.0, 2.0 * shear, 0.0]]))
    assert np.array_equal(conventional, -ell)


def test_wall_shear_is_tangential_projection_of_cauchy_traction() -> None:
    mu = 0.004
    gradient = np.zeros((3, 3))
    gradient[0, 1] = 300.0
    extra = newtonian_extra_stress(gradient, mu)
    sigma = cauchy_stress(12500.0, extra)
    normal = np.array([0.0, 1.0, 0.0])

    tau_w = wall_shear_stress(sigma, normal)

    assert np.allclose(tau_w, np.array([1.2, 0.0, 0.0]))

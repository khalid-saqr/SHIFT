import numpy as np

from shift.operators import (
    CartesianGrid,
    curl,
    divergence,
    first_derivative,
    scalar_gradient,
    scalar_laplacian,
    second_derivative,
    temporal_derivative,
    wall_normal_derivative,
)
from shift.uncertainty import ConvergenceSeries, rms_error


def test_nonuniform_polynomial_derivatives_are_exact_to_roundoff() -> None:
    x = np.array([0.0, 0.04, 0.13, 0.31, 0.58, 0.81, 1.0])
    first = first_derivative(x**2, x, axis=0)
    second = second_derivative(x**3, x, axis=0)
    np.testing.assert_allclose(first, 2.0 * x, rtol=0.0, atol=2e-13)
    np.testing.assert_allclose(second, 6.0 * x, rtol=0.0, atol=2e-11)


def test_cartesian_vector_calculus_on_polynomial_field() -> None:
    x = np.linspace(-0.5, 0.5, 7)
    y = np.linspace(-0.4, 0.6, 8)
    z = np.linspace(-0.3, 0.7, 9)
    grid = CartesianGrid(x, y, z)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

    scalar = X**2 + 2.0 * Y**2 + 3.0 * Z**2
    gradient = scalar_gradient(scalar, grid)
    expected_gradient = np.stack((2.0 * X, 4.0 * Y, 6.0 * Z), axis=-1)
    np.testing.assert_allclose(gradient, expected_gradient, atol=2e-12)
    np.testing.assert_allclose(scalar_laplacian(scalar, grid), 12.0, atol=3e-11)

    vector = np.stack((Y, Z, X), axis=-1)
    np.testing.assert_allclose(divergence(vector, grid), 0.0, atol=2e-13)
    expected_curl = -np.ones(vector.shape)
    np.testing.assert_allclose(curl(vector, grid), expected_curl, atol=2e-13)


def test_temporal_and_wall_normal_derivatives_use_explicit_coordinates() -> None:
    time = np.array([0.0, 0.02, 0.07, 0.15, 0.28])
    signal = time**2
    np.testing.assert_allclose(
        temporal_derivative(signal, time), 2.0 * time, atol=2e-14
    )

    distance = np.array([0.0, 0.00004, 0.00013, 0.00025])
    shear_rate = 300.0
    height = 0.004
    velocity = shear_rate * distance * (1.0 - distance / height)
    wall_gradient = wall_normal_derivative(velocity, distance, axis=0)
    assert np.isclose(wall_gradient, shear_rate, rtol=0.0, atol=1e-10)


def test_first_and_second_derivatives_show_second_order_convergence() -> None:
    first_errors = []
    second_errors = []
    spacing = []
    for count in (17, 33, 65, 129):
        x = np.linspace(0.0, 1.0, count)
        field = np.sin(2.0 * np.pi * x)
        exact_first = 2.0 * np.pi * np.cos(2.0 * np.pi * x)
        exact_second = -(2.0 * np.pi) ** 2 * np.sin(2.0 * np.pi * x)
        first_errors.append(rms_error(first_derivative(field, x, axis=0), exact_first))
        second_errors.append(
            rms_error(second_derivative(field, x, axis=0), exact_second)
        )
        spacing.append(float(np.max(np.diff(x))))

    first_series = ConvergenceSeries(
        np.array(spacing), np.array(first_errors), "first derivative", "s^-1"
    )
    second_series = ConvergenceSeries(
        np.array(spacing), np.array(second_errors), "second derivative", "s^-2"
    )
    assert first_series.observed_order(finest_levels=3) > 1.9
    assert second_series.observed_order(finest_levels=3) > 1.9

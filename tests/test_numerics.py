import numpy as np

from shift.mechanics import (
    cauchy_stress,
    newtonian_extra_stress,
    wall_shear_stress,
)
from shift.models.taylor_green import TaylorGreenVortex
from shift.numerics import SampledFlow, reconstruct_source04
from shift.operators import CartesianGrid, wall_normal_derivative


def _taylor_green_result(count: int):
    rho = 1050.0
    mu = 0.004
    model = TaylorGreenVortex(0.2, 200.0, rho, mu)
    length = 2.0 * np.pi / model.wave_number_m_1
    x = np.linspace(0.0, length, count)
    y = np.linspace(0.0, length, count)
    z = np.linspace(0.0, 0.003, 4)
    time = np.linspace(0.0, 0.4, count)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    velocity = []
    pressure = []
    for instant in time:
        velocity.append(model.velocity(X, Y, Z, instant))
        pressure.append(model.pressure_pa(X, Y, Z, instant))
    flow = SampledFlow(
        np.asarray(velocity),
        np.asarray(pressure),
        CartesianGrid(x, y, z),
        rho,
        mu,
        time_s=time,
    )
    return reconstruct_source04(flow)


def test_taylor_green_identity_residuals_converge_under_refinement() -> None:
    coarse = _taylor_green_result(9).residuals
    fine = _taylor_green_result(17).residuals
    assert fine.incompressibility_s_1 < 0.35 * coarse.incompressibility_s_1
    assert fine.charge_s_2 < 0.5 * coarse.charge_s_2
    assert fine.crocco_m_s2 < 0.25 * coarse.crocco_m_s2


def test_wall_normal_operator_recovers_exact_poiseuille_wss() -> None:
    viscosity = 0.004
    target_wss = 1.2
    shear_rate = target_wss / viscosity
    height = 0.004
    distance = np.array([0.0, 0.00005, 0.00017, 0.00031])
    velocity_x = shear_rate * distance * (1.0 - distance / height)
    du_dn = wall_normal_derivative(velocity_x, distance, axis=0)

    gradient = np.zeros((3, 3))
    gradient[0, 1] = du_dn
    stress = cauchy_stress(0.0, newtonian_extra_stress(gradient, viscosity))
    wss = wall_shear_stress(stress, np.array([0.0, 1.0, 0.0]))
    np.testing.assert_allclose(
        wss,
        np.array([target_wss, 0.0, 0.0]),
        atol=1e-12,
    )

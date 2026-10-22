"""Exact steady parallel-plate Poiseuille benchmark used by source04."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from ..contracts import CharacteristicScales
from ..mechanics import (
    CroccoComposition,
    LambUnsteadyDiagnostics,
    arealized_finite_volume_lamb,
    cauchy_stress,
    conventional_lamb_vector,
    crocco_composition,
    lamb_unsteady_diagnostics,
    lamb_vector,
    newtonian_extra_stress,
    second_invariant_q,
    vorticity_from_gradient,
    wall_shear_stress,
)
from ..metrics import trapezoid_measure_weights
from ..state import PrimaryMechanicalState

FloatArray = NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class ParallelPlatePoiseuilleState:
    """Sampled exact fields on the lower-wall near-wall layer."""

    height_m: float
    pressure_drop_gradient_pa_m: float
    y_m: FloatArray
    velocity_m_s: FloatArray
    velocity_gradient_s_1: FloatArray
    velocity_divergence_s_1: FloatArray
    vorticity_s_1: FloatArray
    q_invariant_s_2: FloatArray
    lamb_m_s2: FloatArray
    conventional_lamb_m_s2: FloatArray
    hydrodynamic_charge_s_2: FloatArray
    laplacian_kinetic_energy_s_2: FloatArray
    c_unsteady_m_s2: FloatArray
    c_head_m_s2: FloatArray
    c_dissipative_m_s2: FloatArray
    composition: CroccoComposition
    lamb_unsteady: LambUnsteadyDiagnostics
    lambda_l_pa: FloatArray
    lambda_l_exact_pa: FloatArray
    wall_shear_pa: FloatArray

    @property
    def primary_state(self) -> PrimaryMechanicalState:
        return PrimaryMechanicalState(
            self.q_invariant_s_2,
            self.hydrodynamic_charge_s_2,
            self.c_unsteady_m_s2,
            self.c_head_m_s2,
            self.c_dissipative_m_s2,
        )


def exact_parallel_plate_poiseuille(
    *,
    height_m: float,
    near_wall_thickness_m: float,
    density_kg_m3: float,
    viscosity_pa_s: float,
    wall_shear_pa: float,
    scales: CharacteristicScales,
    activity_fraction: float = 1e-3,
    sample_count: int = 4001,
) -> ParallelPlatePoiseuilleState:
    """Construct the exact lower-wall source04 benchmark without finite differences."""
    for name, value in (
        ("height_m", height_m),
        ("near_wall_thickness_m", near_wall_thickness_m),
        ("density_kg_m3", density_kg_m3),
        ("viscosity_pa_s", viscosity_pa_s),
        ("wall_shear_pa", wall_shear_pa),
    ):
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and strictly positive")
    if height_m <= 2.0 * near_wall_thickness_m:
        raise ValueError(
            "height_m must exceed twice near_wall_thickness_m so the declared "
            "near-wall layer remains in the lower half-channel"
        )
    if sample_count < 3:
        raise ValueError("sample_count must be at least 3")

    y = np.linspace(0.0, near_wall_thickness_m, sample_count)
    shear_rate = wall_shear_pa / viscosity_pa_s
    pressure_drop_gradient = 2.0 * wall_shear_pa / height_m
    u = shear_rate * y * (1.0 - y / height_m)
    du_dy = shear_rate * (1.0 - 2.0 * y / height_m)
    d2u_dy2 = np.full_like(y, -2.0 * shear_rate / height_m)

    velocity = np.zeros((sample_count, 3), dtype=float)
    velocity[:, 0] = u
    gradient = np.zeros((sample_count, 3, 3), dtype=float)
    gradient[:, 0, 1] = du_dy
    divergence = np.trace(gradient, axis1=-2, axis2=-1)
    omega = vorticity_from_gradient(gradient)
    q = second_invariant_q(gradient)
    ell = lamb_vector(velocity, omega)
    conventional_lamb = conventional_lamb_vector(velocity, omega)

    lap_k = du_dy**2 + u * d2u_dy2
    q_h = -(du_dy**2 + u * d2u_dy2)

    c_u = np.zeros_like(velocity)
    c_h = np.zeros_like(velocity)
    c_h[:, 0] = -pressure_drop_gradient / density_kg_m3
    c_h[:, 1] = ell[:, 1]
    c_d = np.zeros_like(velocity)
    c_d[:, 0] = pressure_drop_gradient / density_kg_m3

    composition = crocco_composition(
        ell,
        c_u,
        c_h,
        c_d,
        acceleration_scale_m_s2=scales.acceleration_m_s2,
        activity_fraction=activity_fraction,
    )
    lu = lamb_unsteady_diagnostics(
        ell,
        c_u,
        acceleration_scale_m_s2=scales.acceleration_m_s2,
        activity_fraction=activity_fraction,
    )

    line_weights = trapezoid_measure_weights(y)
    lambda_l = arealized_finite_volume_lamb(
        ell,
        density_kg_m3=density_kg_m3,
        volume_weights_m3=line_weights,
        wall_area_m2=1.0,
    )
    lambda_exact = np.array([0.0, 0.5 * density_kg_m3 * u[-1] ** 2, 0.0])

    wall_gradient = gradient[0]
    extra = newtonian_extra_stress(wall_gradient, viscosity_pa_s)
    sigma = cauchy_stress(0.0, extra)
    wall_normal_into_fluid = np.array([0.0, 1.0, 0.0])
    tau_w = wall_shear_stress(sigma, wall_normal_into_fluid)

    return ParallelPlatePoiseuilleState(
        height_m=float(height_m),
        pressure_drop_gradient_pa_m=float(pressure_drop_gradient),
        y_m=y,
        velocity_m_s=velocity,
        velocity_gradient_s_1=gradient,
        velocity_divergence_s_1=divergence,
        vorticity_s_1=omega,
        q_invariant_s_2=q,
        lamb_m_s2=ell,
        conventional_lamb_m_s2=conventional_lamb,
        hydrodynamic_charge_s_2=q_h,
        laplacian_kinetic_energy_s_2=lap_k,
        c_unsteady_m_s2=c_u,
        c_head_m_s2=c_h,
        c_dissipative_m_s2=c_d,
        composition=composition,
        lamb_unsteady=lu,
        lambda_l_pa=lambda_l,
        lambda_l_exact_pa=lambda_exact,
        wall_shear_pa=tau_w,
    )


def closed_form_d_q(
    *,
    height_a_m: float,
    height_b_m: float,
    near_wall_thickness_m: float,
) -> float:
    """Closed-form d_q for the locked scaling U_c=(tau/mu) delta, L_c=delta."""
    d1 = 1.0 / height_a_m - 1.0 / height_b_m
    d2 = 1.0 / height_a_m**2 - 1.0 / height_b_m**2
    delta = near_wall_thickness_m
    squared = 36.0 * (
        delta**2 * d1**2 / 3.0
        - delta**3 * d1 * d2 / 2.0
        + delta**4 * d2**2 / 5.0
    )
    return float(np.sqrt(squared))

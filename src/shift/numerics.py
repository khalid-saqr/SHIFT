"""Numerical reconstruction of source04 fields from sampled velocity and pressure."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .certify import IdentityResiduals, identity_residuals
from .mechanics import (
    conventional_lamb_vector,
    hydrodynamic_charge_residual,
    lamb_vector,
    second_invariant_q,
    vorticity_from_gradient,
)
from .operators import (
    CartesianGrid,
    divergence,
    scalar_gradient,
    scalar_laplacian,
    temporal_derivative,
    vector_gradient,
    vector_laplacian,
)
from .state import PrimaryMechanicalState

FloatArray = NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class SampledFlow:
    """Velocity/pressure samples with explicit spatial and temporal support."""

    velocity_m_s: FloatArray
    pressure_pa: FloatArray
    grid: CartesianGrid
    density_kg_m3: float
    viscosity_pa_s: float
    time_s: FloatArray | None = None
    steady: bool = False
    potential_m2_s2: FloatArray | None = None

    def __post_init__(self) -> None:
        velocity = np.asarray(self.velocity_m_s, dtype=float)
        pressure = np.asarray(self.pressure_pa, dtype=float)
        if not np.isfinite(self.density_kg_m3) or self.density_kg_m3 <= 0.0:
            raise ValueError("density_kg_m3 must be finite and > 0")
        if not np.isfinite(self.viscosity_pa_s) or self.viscosity_pa_s <= 0.0:
            raise ValueError("viscosity_pa_s must be finite and > 0")
        if self.steady:
            expected_velocity = self.grid.shape + (3,)
            expected_pressure = self.grid.shape
            if self.time_s is not None:
                raise ValueError("steady flow must not provide time_s")
        else:
            if self.time_s is None:
                raise ValueError("unsteady flow requires explicit time_s")
            time = np.asarray(self.time_s, dtype=float)
            if time.ndim != 1 or time.size < 3:
                raise ValueError("time_s must be a 1-D coordinate with at least 3 points")
            if not np.all(np.isfinite(time)) or not np.all(np.diff(time) > 0.0):
                raise ValueError("time_s must be finite and strictly increasing")
            expected_velocity = (time.size,) + self.grid.shape + (3,)
            expected_pressure = (time.size,) + self.grid.shape
            object.__setattr__(self, "time_s", time)
        if velocity.shape != expected_velocity or pressure.shape != expected_pressure:
            raise ValueError("velocity/pressure shapes do not match declared support")
        if not np.all(np.isfinite(velocity)) or not np.all(np.isfinite(pressure)):
            raise ValueError("velocity and pressure samples must be finite")
        object.__setattr__(self, "velocity_m_s", velocity)
        object.__setattr__(self, "pressure_pa", pressure)
        if self.potential_m2_s2 is not None:
            potential = np.asarray(self.potential_m2_s2, dtype=float)
            if potential.shape != expected_pressure or not np.all(np.isfinite(potential)):
                raise ValueError("potential must match scalar flow support and be finite")
            object.__setattr__(self, "potential_m2_s2", potential)


@dataclass(frozen=True, slots=True)
class NumericalSource04Fields:
    """Numerically reconstructed source04 fields and raw identity residuals."""

    primary_state: PrimaryMechanicalState
    velocity_gradient_s_1: FloatArray
    velocity_divergence_s_1: FloatArray
    vorticity_s_1: FloatArray
    lamb_m_s2: FloatArray
    conventional_lamb_m_s2: FloatArray
    laplacian_kinetic_energy_s_2: FloatArray
    hydrodynamic_charge_identity_s_2: FloatArray
    charge_identity_residual_s_2: FloatArray
    residuals: IdentityResiduals


def reconstruct_source04(flow: SampledFlow) -> NumericalSource04Fields:
    """Reconstruct the primary source04 state with independent charge routes."""
    velocity = flow.velocity_m_s
    gradient = vector_gradient(velocity, flow.grid)
    velocity_divergence = np.trace(gradient, axis1=-2, axis2=-1)
    q_invariant = second_invariant_q(gradient)
    vorticity = vorticity_from_gradient(gradient)
    ell = lamb_vector(velocity, vorticity)
    conventional_lamb = conventional_lamb_vector(velocity, vorticity)
    q_h = divergence(conventional_lamb, flow.grid)

    kinetic_energy = 0.5 * np.sum(velocity**2, axis=-1)
    laplacian_k = scalar_laplacian(kinetic_energy, flow.grid)
    q_h_identity = -2.0 * q_invariant - laplacian_k
    charge_residual = hydrodynamic_charge_residual(q_h, q_invariant, laplacian_k)

    if flow.steady:
        c_unsteady = np.zeros_like(velocity)
    else:
        c_unsteady = temporal_derivative(velocity, flow.time_s, axis=0)

    total_head = flow.pressure_pa / flow.density_kg_m3 + kinetic_energy
    if flow.potential_m2_s2 is not None:
        total_head = total_head + flow.potential_m2_s2
    c_head = scalar_gradient(total_head, flow.grid)
    kinematic_viscosity = flow.viscosity_pa_s / flow.density_kg_m3
    c_dissipative = -kinematic_viscosity * vector_laplacian(velocity, flow.grid)

    primary_state = PrimaryMechanicalState(
        q_invariant,
        q_h,
        c_unsteady,
        c_head,
        c_dissipative,
    )
    residuals = identity_residuals(
        velocity_divergence_s_1=velocity_divergence,
        hydrodynamic_charge_s_2=q_h,
        q_invariant_s_2=q_invariant,
        laplacian_kinetic_energy_s_2=laplacian_k,
        lamb_m_s2=ell,
        c_unsteady_m_s2=c_unsteady,
        c_head_m_s2=c_head,
        c_dissipative_m_s2=c_dissipative,
    )
    return NumericalSource04Fields(
        primary_state=primary_state,
        velocity_gradient_s_1=gradient,
        velocity_divergence_s_1=velocity_divergence,
        vorticity_s_1=vorticity,
        lamb_m_s2=ell,
        conventional_lamb_m_s2=conventional_lamb,
        laplacian_kinetic_energy_s_2=laplacian_k,
        hydrodynamic_charge_identity_s_2=q_h_identity,
        charge_identity_residual_s_2=charge_residual,
        residuals=residuals,
    )

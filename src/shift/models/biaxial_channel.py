"""Exact multiharmonic biaxial pulsatile parallel-plate Navier--Stokes flow."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from shift.certify import IdentityResiduals, identity_residuals
from shift.contracts import CharacteristicScales
from shift.mechanics import (
    CroccoComposition,
    LambUnsteadyDiagnostics,
    conventional_lamb_vector,
    crocco_composition,
    lamb_unsteady_diagnostics,
    lamb_vector,
    second_invariant_q,
    vorticity_from_gradient,
)
from shift.state import PrimaryMechanicalState

FloatArray = NDArray[np.float64]
ComplexArray = NDArray[np.complex128]


@dataclass(frozen=True, slots=True)
class BiaxialWallShearWaveform:
    """Protocol-form wall-shear waveform with axial and transverse harmonics."""

    shear_scale_pa: float
    axial_first_harmonic: float
    axial_second_harmonic: float
    transverse_first_harmonic: float

    def __post_init__(self) -> None:
        values = (
            self.shear_scale_pa,
            self.axial_first_harmonic,
            self.axial_second_harmonic,
            self.transverse_first_harmonic,
        )
        if not all(np.isfinite(value) for value in values):
            raise ValueError("wall-shear waveform parameters must be finite")
        if self.shear_scale_pa <= 0.0:
            raise ValueError("shear_scale_pa must be strictly positive")

    def evaluate(self, time_s: ArrayLike, *, period_s: float) -> FloatArray:
        """Evaluate the prescribed lower-wall vector WSS [Pa]."""
        time = np.asarray(time_s, dtype=float)
        if not np.isfinite(period_s) or period_s <= 0.0:
            raise ValueError("period_s must be finite and strictly positive")
        phase = 2.0 * np.pi * time / period_s
        result = np.zeros(time.shape + (3,), dtype=float)
        result[..., 0] = self.shear_scale_pa * (
            1.0
            + self.axial_first_harmonic * np.sin(phase)
            + self.axial_second_harmonic * np.sin(2.0 * phase)
        )
        result[..., 2] = (
            self.shear_scale_pa
            * self.transverse_first_harmonic
            * np.cos(phase)
        )
        return result

    def complex_harmonics_pa(self) -> ComplexArray:
        """Return complex WSS amplitudes for harmonics k=1 and k=2."""
        harmonics = np.zeros((2, 3), dtype=complex)
        harmonics[0, 0] = -1j * self.shear_scale_pa * self.axial_first_harmonic
        harmonics[1, 0] = -1j * self.shear_scale_pa * self.axial_second_harmonic
        harmonics[0, 2] = self.shear_scale_pa * self.transverse_first_harmonic
        return harmonics

    @property
    def mean_pa(self) -> FloatArray:
        """Mean vector wall shear [Pa]."""
        return np.array([self.shear_scale_pa, 0.0, 0.0], dtype=float)


def pressure_drop_harmonic_from_wall_shear(
    wall_shear_hat_pa: complex,
    *,
    angular_frequency_rad_s: float,
    height_m: float,
    kinematic_viscosity_m2_s: float,
) -> complex:
    """Invert the exact harmonic WSS transfer for the pressure-drop gradient."""
    if angular_frequency_rad_s <= 0.0 or not np.isfinite(angular_frequency_rad_s):
        raise ValueError("angular_frequency_rad_s must be finite and > 0")
    if height_m <= 0.0 or not np.isfinite(height_m):
        raise ValueError("height_m must be finite and > 0")
    if kinematic_viscosity_m2_s <= 0.0 or not np.isfinite(
        kinematic_viscosity_m2_s
    ):
        raise ValueError("kinematic_viscosity_m2_s must be finite and > 0")
    lam = np.sqrt(1j * angular_frequency_rad_s / kinematic_viscosity_m2_s)
    transfer = (
        kinematic_viscosity_m2_s
        * lam
        / (1j * angular_frequency_rad_s)
        * np.tanh(lam * height_m / 2.0)
    )
    return complex(wall_shear_hat_pa / transfer)


def _validate_coordinate(
    value: ArrayLike,
    name: str,
    *,
    lower: float,
    upper: float,
    minimum_points: int,
) -> FloatArray:
    coordinate = np.asarray(value, dtype=float)
    if coordinate.ndim != 1 or coordinate.size < minimum_points:
        raise ValueError(f"{name} must be 1-D with at least {minimum_points} points")
    if not np.all(np.isfinite(coordinate)) or not np.all(np.diff(coordinate) > 0.0):
        raise ValueError(f"{name} must be finite and strictly increasing")
    scale = max(1.0, abs(lower), abs(upper))
    tolerance = 1e-12 * scale
    if coordinate[0] < lower - tolerance or coordinate[-1] > upper + tolerance:
        raise ValueError(f"{name} lies outside the declared physical interval")
    return coordinate


def _harmonic_profiles(
    y_m: FloatArray,
    *,
    pressure_drop_hat_pa_m: complex,
    angular_frequency_rad_s: float,
    height_m: float,
    density_kg_m3: float,
    kinematic_viscosity_m2_s: float,
) -> tuple[ComplexArray, ComplexArray, ComplexArray]:
    lam = np.sqrt(1j * angular_frequency_rad_s / kinematic_viscosity_m2_s)
    denominator = np.cosh(lam * height_m / 2.0)
    shifted = lam * (y_m - height_m / 2.0)
    prefactor = pressure_drop_hat_pa_m / (
        1j * density_kg_m3 * angular_frequency_rad_s
    )
    velocity = prefactor * (1.0 - np.cosh(shifted) / denominator)
    first = prefactor * (-lam * np.sinh(shifted) / denominator)
    second = prefactor * (-lam**2 * np.cosh(shifted) / denominator)
    return velocity, first, second


def _evaluate_component(
    y_m: FloatArray,
    time_s: FloatArray,
    *,
    mean_wall_shear_pa: float,
    wall_shear_harmonics_pa: ComplexArray,
    period_s: float,
    height_m: float,
    density_kg_m3: float,
    dynamic_viscosity_pa_s: float,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray, ComplexArray]:
    """Evaluate one tangential velocity component and its exact derivatives."""
    kinematic_viscosity = dynamic_viscosity_pa_s / density_kg_m3
    velocity = np.zeros((time_s.size, y_m.size), dtype=float)
    first_y = np.zeros_like(velocity)
    second_y = np.zeros_like(velocity)
    first_t = np.zeros_like(velocity)
    pressure_drop = np.full(time_s.size, 2.0 * mean_wall_shear_pa / height_m)
    pressure_harmonics = np.zeros(2, dtype=complex)

    mean_gradient = 2.0 * mean_wall_shear_pa / height_m
    velocity += mean_gradient / (2.0 * dynamic_viscosity_pa_s) * (
        y_m * (height_m - y_m)
    )
    first_y += mean_gradient / (2.0 * dynamic_viscosity_pa_s) * (
        height_m - 2.0 * y_m
    )
    second_y += -mean_gradient / dynamic_viscosity_pa_s

    base_omega = 2.0 * np.pi / period_s
    for index, wall_shear_hat in enumerate(wall_shear_harmonics_pa, start=1):
        if wall_shear_hat == 0.0:
            continue
        angular_frequency = index * base_omega
        pressure_hat = pressure_drop_harmonic_from_wall_shear(
            complex(wall_shear_hat),
            angular_frequency_rad_s=angular_frequency,
            height_m=height_m,
            kinematic_viscosity_m2_s=kinematic_viscosity,
        )
        pressure_harmonics[index - 1] = pressure_hat
        u_hat, du_hat, d2u_hat = _harmonic_profiles(
            y_m,
            pressure_drop_hat_pa_m=pressure_hat,
            angular_frequency_rad_s=angular_frequency,
            height_m=height_m,
            density_kg_m3=density_kg_m3,
            kinematic_viscosity_m2_s=kinematic_viscosity,
        )
        temporal_phase = np.exp(1j * angular_frequency * time_s)
        phase = temporal_phase[:, None]
        velocity += np.real(phase * u_hat[None, :])
        first_y += np.real(phase * du_hat[None, :])
        second_y += np.real(phase * d2u_hat[None, :])
        first_t += np.real(
            1j * angular_frequency * phase * u_hat[None, :]
        )
        pressure_drop += np.real(pressure_hat * temporal_phase)

    return (
        velocity,
        first_y,
        second_y,
        first_t,
        pressure_drop,
        pressure_harmonics,
    )


@dataclass(frozen=True, slots=True)
class ExactBiaxialChannelState:
    """Exact periodic state and analytical source04 fields on supplied support."""

    height_m: float
    period_s: float
    near_wall_thickness_m: float
    density_kg_m3: float
    dynamic_viscosity_pa_s: float
    waveform: BiaxialWallShearWaveform
    y_m: FloatArray
    time_s: FloatArray
    velocity_m_s: FloatArray
    du_dy_s_1: FloatArray
    d2u_dy2_m_1_s: FloatArray
    du_dt_m_s2: FloatArray
    pressure_drop_gradient_pa_m: FloatArray
    pressure_gradient_pa_m: FloatArray
    velocity_gradient_s_1: FloatArray
    velocity_divergence_s_1: FloatArray
    vorticity_s_1: FloatArray
    kinetic_energy_m2_s2: FloatArray
    q_invariant_s_2: FloatArray
    lamb_m_s2: FloatArray
    conventional_lamb_m_s2: FloatArray
    hydrodynamic_charge_s_2: FloatArray
    laplacian_kinetic_energy_s_2: FloatArray
    c_unsteady_m_s2: FloatArray
    c_head_m_s2: FloatArray
    c_dissipative_m_s2: FloatArray
    wall_shear_target_pa: FloatArray
    wall_shear_exact_pa: FloatArray
    pressure_drop_harmonics_pa_m: ComplexArray
    composition: CroccoComposition
    lamb_unsteady: LambUnsteadyDiagnostics
    lambda_l_pa: FloatArray
    residuals: IdentityResiduals

    @property
    def primary_state(self) -> PrimaryMechanicalState:
        return PrimaryMechanicalState(
            self.q_invariant_s_2,
            self.hydrodynamic_charge_s_2,
            self.c_unsteady_m_s2,
            self.c_head_m_s2,
            self.c_dissipative_m_s2,
        )

    def pressure_pa(
        self,
        x_m: ArrayLike,
        z_m: ArrayLike,
        *,
        reference_pressure_pa: float = 0.0,
    ) -> FloatArray:
        """Evaluate p=p0-G_x(t)x-G_z(t)z on any broadcastable x-z support."""
        if not np.isfinite(reference_pressure_pa):
            raise ValueError("reference_pressure_pa must be finite")
        x, z = np.broadcast_arrays(
            np.asarray(x_m, dtype=float), np.asarray(z_m, dtype=float)
        )
        if not np.all(np.isfinite(x)) or not np.all(np.isfinite(z)):
            raise ValueError("pressure coordinates must be finite")
        reshape = (self.time_s.size,) + (1,) * x.ndim
        gx = self.pressure_drop_gradient_pa_m[:, 0].reshape(reshape)
        gz = self.pressure_drop_gradient_pa_m[:, 2].reshape(reshape)
        return reference_pressure_pa - gx * x - gz * z


def exact_biaxial_pulsatile_channel(
    *,
    height_m: float,
    period_s: float,
    near_wall_thickness_m: float,
    density_kg_m3: float,
    dynamic_viscosity_pa_s: float,
    waveform: BiaxialWallShearWaveform,
    scales: CharacteristicScales,
    activity_fraction: float,
    y_m: ArrayLike,
    time_s: ArrayLike,
) -> ExactBiaxialChannelState:
    """Generate the exact periodic biaxial channel state on explicit y-t support."""
    for name, value in (
        ("height_m", height_m),
        ("period_s", period_s),
        ("near_wall_thickness_m", near_wall_thickness_m),
        ("density_kg_m3", density_kg_m3),
        ("dynamic_viscosity_pa_s", dynamic_viscosity_pa_s),
    ):
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and strictly positive")
    if height_m <= 2.0 * near_wall_thickness_m:
        raise ValueError(
            "height_m must exceed twice near_wall_thickness_m for the locked "
            "lower-wall slab"
        )
    if activity_fraction < 0.0 or not np.isfinite(activity_fraction):
        raise ValueError("activity_fraction must be finite and non-negative")
    y = _validate_coordinate(
        y_m,
        "y_m",
        lower=0.0,
        upper=height_m,
        minimum_points=2,
    )
    time = _validate_coordinate(
        time_s,
        "time_s",
        lower=0.0,
        upper=period_s,
        minimum_points=2,
    )

    wall_harmonics = waveform.complex_harmonics_pa()
    ux = _evaluate_component(
        y,
        time,
        mean_wall_shear_pa=waveform.mean_pa[0],
        wall_shear_harmonics_pa=wall_harmonics[:, 0],
        period_s=period_s,
        height_m=height_m,
        density_kg_m3=density_kg_m3,
        dynamic_viscosity_pa_s=dynamic_viscosity_pa_s,
    )
    uz = _evaluate_component(
        y,
        time,
        mean_wall_shear_pa=waveform.mean_pa[2],
        wall_shear_harmonics_pa=wall_harmonics[:, 2],
        period_s=period_s,
        height_m=height_m,
        density_kg_m3=density_kg_m3,
        dynamic_viscosity_pa_s=dynamic_viscosity_pa_s,
    )

    velocity = np.zeros((time.size, y.size, 3), dtype=float)
    du_dy = np.zeros_like(velocity)
    d2u_dy2 = np.zeros_like(velocity)
    du_dt = np.zeros_like(velocity)
    pressure_drop = np.zeros((time.size, 3), dtype=float)
    pressure_harmonics = np.zeros((2, 3), dtype=complex)
    for component, values in ((0, ux), (2, uz)):
        velocity[..., component] = values[0]
        du_dy[..., component] = values[1]
        d2u_dy2[..., component] = values[2]
        du_dt[..., component] = values[3]
        pressure_drop[:, component] = values[4]
        pressure_harmonics[:, component] = values[5]

    gradient = np.zeros((time.size, y.size, 3, 3), dtype=float)
    gradient[..., 0, 1] = du_dy[..., 0]
    gradient[..., 2, 1] = du_dy[..., 2]
    divergence = np.trace(gradient, axis1=-2, axis2=-1)
    vorticity = vorticity_from_gradient(gradient)
    q_invariant = second_invariant_q(gradient)
    lamb = lamb_vector(velocity, vorticity)
    conventional_lamb = conventional_lamb_vector(velocity, vorticity)

    kinetic_energy = 0.5 * np.sum(velocity**2, axis=-1)
    kinetic_energy_y = np.sum(velocity * du_dy, axis=-1)
    kinetic_energy_yy = np.sum(du_dy**2 + velocity * d2u_dy2, axis=-1)
    laplacian_k = kinetic_energy_yy
    q_h = -kinetic_energy_yy

    c_unsteady = du_dt.copy()
    c_head = np.zeros_like(velocity)
    c_head[..., 0] = -pressure_drop[:, None, 0] / density_kg_m3
    c_head[..., 1] = kinetic_energy_y
    c_head[..., 2] = -pressure_drop[:, None, 2] / density_kg_m3
    kinematic_viscosity = dynamic_viscosity_pa_s / density_kg_m3
    c_dissipative = -kinematic_viscosity * d2u_dy2

    target_wss = waveform.evaluate(time, period_s=period_s)
    exact_wss = dynamic_viscosity_pa_s * du_dy[:, 0, :]
    pressure_gradient = -pressure_drop

    composition = crocco_composition(
        lamb,
        c_unsteady,
        c_head,
        c_dissipative,
        acceleration_scale_m_s2=scales.acceleration_m_s2,
        activity_fraction=activity_fraction,
    )
    lu = lamb_unsteady_diagnostics(
        lamb,
        c_unsteady,
        acceleration_scale_m_s2=scales.acceleration_m_s2,
        activity_fraction=activity_fraction,
    )

    delta = np.array([0.0, near_wall_thickness_m], dtype=float)
    ux_delta = _evaluate_component(
        delta,
        time,
        mean_wall_shear_pa=waveform.mean_pa[0],
        wall_shear_harmonics_pa=wall_harmonics[:, 0],
        period_s=period_s,
        height_m=height_m,
        density_kg_m3=density_kg_m3,
        dynamic_viscosity_pa_s=dynamic_viscosity_pa_s,
    )[0][:, 1]
    uz_delta = _evaluate_component(
        delta,
        time,
        mean_wall_shear_pa=waveform.mean_pa[2],
        wall_shear_harmonics_pa=wall_harmonics[:, 2],
        period_s=period_s,
        height_m=height_m,
        density_kg_m3=density_kg_m3,
        dynamic_viscosity_pa_s=dynamic_viscosity_pa_s,
    )[0][:, 1]
    kinetic_energy_delta = 0.5 * (ux_delta**2 + uz_delta**2)
    lambda_l = np.zeros((time.size, 3), dtype=float)
    lambda_l[:, 1] = density_kg_m3 * kinetic_energy_delta

    residuals = identity_residuals(
        velocity_divergence_s_1=divergence,
        hydrodynamic_charge_s_2=q_h,
        q_invariant_s_2=q_invariant,
        laplacian_kinetic_energy_s_2=laplacian_k,
        lamb_m_s2=lamb,
        c_unsteady_m_s2=c_unsteady,
        c_head_m_s2=c_head,
        c_dissipative_m_s2=c_dissipative,
    )

    return ExactBiaxialChannelState(
        height_m=float(height_m),
        period_s=float(period_s),
        near_wall_thickness_m=float(near_wall_thickness_m),
        density_kg_m3=float(density_kg_m3),
        dynamic_viscosity_pa_s=float(dynamic_viscosity_pa_s),
        waveform=waveform,
        y_m=y,
        time_s=time,
        velocity_m_s=velocity,
        du_dy_s_1=du_dy,
        d2u_dy2_m_1_s=d2u_dy2,
        du_dt_m_s2=du_dt,
        pressure_drop_gradient_pa_m=pressure_drop,
        pressure_gradient_pa_m=pressure_gradient,
        velocity_gradient_s_1=gradient,
        velocity_divergence_s_1=divergence,
        vorticity_s_1=vorticity,
        kinetic_energy_m2_s2=kinetic_energy,
        q_invariant_s_2=q_invariant,
        lamb_m_s2=lamb,
        conventional_lamb_m_s2=conventional_lamb,
        hydrodynamic_charge_s_2=q_h,
        laplacian_kinetic_energy_s_2=laplacian_k,
        c_unsteady_m_s2=c_unsteady,
        c_head_m_s2=c_head,
        c_dissipative_m_s2=c_dissipative,
        wall_shear_target_pa=target_wss,
        wall_shear_exact_pa=exact_wss,
        pressure_drop_harmonics_pa_m=pressure_harmonics,
        composition=composition,
        lamb_unsteady=lu,
        lambda_l_pa=lambda_l,
        residuals=residuals,
    )

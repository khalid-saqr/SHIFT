"""Exact dimensional Taylor-Green vortex used for Phase-2 verification."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class TaylorGreenVortex:
    """Two-dimensional incompressible Navier-Stokes solution embedded in 3-D."""

    velocity_scale_m_s: float
    wave_number_m_1: float
    density_kg_m3: float
    viscosity_pa_s: float

    def __post_init__(self) -> None:
        for name, value in (
            ("velocity_scale_m_s", self.velocity_scale_m_s),
            ("wave_number_m_1", self.wave_number_m_1),
            ("density_kg_m3", self.density_kg_m3),
            ("viscosity_pa_s", self.viscosity_pa_s),
        ):
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and strictly positive")

    @property
    def kinematic_viscosity_m2_s(self) -> float:
        return self.viscosity_pa_s / self.density_kg_m3

    def decay(self, time_s: ArrayLike) -> FloatArray:
        time = np.asarray(time_s, dtype=float)
        if not np.all(np.isfinite(time)) or np.any(time < 0.0):
            raise ValueError("time_s must be finite and non-negative")
        rate = 2.0 * self.kinematic_viscosity_m2_s * self.wave_number_m_1**2
        return np.exp(-rate * time)

    def velocity(
        self,
        x_m: ArrayLike,
        y_m: ArrayLike,
        z_m: ArrayLike,
        time_s: ArrayLike,
    ) -> FloatArray:
        x, y, z, time = np.broadcast_arrays(
            np.asarray(x_m, dtype=float),
            np.asarray(y_m, dtype=float),
            np.asarray(z_m, dtype=float),
            np.asarray(time_s, dtype=float),
        )
        del z
        amplitude = self.velocity_scale_m_s * self.decay(time)
        phase_x = self.wave_number_m_1 * x
        phase_y = self.wave_number_m_1 * y
        u = amplitude * np.sin(phase_x) * np.cos(phase_y)
        v = -amplitude * np.cos(phase_x) * np.sin(phase_y)
        w = np.zeros_like(u)
        return np.stack((u, v, w), axis=-1)

    def pressure_pa(
        self,
        x_m: ArrayLike,
        y_m: ArrayLike,
        z_m: ArrayLike,
        time_s: ArrayLike,
    ) -> FloatArray:
        x, y, z, time = np.broadcast_arrays(
            np.asarray(x_m, dtype=float),
            np.asarray(y_m, dtype=float),
            np.asarray(z_m, dtype=float),
            np.asarray(time_s, dtype=float),
        )
        del z
        amplitude_sq = (self.velocity_scale_m_s * self.decay(time)) ** 2
        return (
            0.25
            * self.density_kg_m3
            * amplitude_sq
            * (
                np.cos(2.0 * self.wave_number_m_1 * x)
                + np.cos(2.0 * self.wave_number_m_1 * y)
            )
        )

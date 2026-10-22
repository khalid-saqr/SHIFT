"""Locked Phase-3 state family and exact-state generation wrapper."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np

from .contracts import CharacteristicScales, StateWeights
from .models.biaxial_channel import (
    BiaxialWallShearWaveform,
    ExactBiaxialChannelState,
    exact_biaxial_pulsatile_channel,
)

PHASE3_RHO_KG_M3 = 1050.0
PHASE3_MU_PA_S = 0.004
PHASE3_DELTA_M = 0.0004
PHASE3_LC_M = 0.0004
PHASE3_UC_M_S = 0.12
PHASE3_TAU_C_PA = 1.2
PHASE3_ACTIVITY_FRACTION = 1e-3
PHASE3_SCALES = CharacteristicScales(
    PHASE3_UC_M_S,
    PHASE3_LC_M,
    PHASE3_RHO_KG_M3,
    PHASE3_TAU_C_PA,
)
PHASE3_STATE_WEIGHTS = StateWeights(0.2, 0.2, 0.2, 0.2, 0.2)

PHASE3_HEIGHTS_M = (0.002, 0.003, 0.004, 0.006, 0.008, 0.012, 0.016)
PHASE3_PERIODS_S = (0.6, 0.8, 1.0, 1.2)
PHASE3_A_VALUES = (0.0, 0.75, 1.5)
PHASE3_B_VALUES = (0.0, 0.3)
PHASE3_C_VALUES = (0.0, 0.6, 1.2)
PHASE3_NY_LEVELS = (17, 33, 65, 129, 257)
PHASE3_NT_LEVELS = (32, 64, 128, 256, 512)
PHASE3_PRODUCTION_NY = 257
PHASE3_PRODUCTION_NT = 512


def _locked_member(value: float, allowed: tuple[float, ...], name: str) -> float:
    numeric = float(value)
    if numeric not in allowed:
        raise ValueError(f"{name}={numeric!r} is outside the locked Phase-3 grid")
    return numeric


def _token(value: float) -> str:
    return f"{value:g}".replace(".", "p")


@dataclass(frozen=True, slots=True)
class Phase3Case:
    """One immutable member of the locked 504-state Phase-3 family."""

    height_m: float
    period_s: float
    a: float
    b: float
    c: float

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "height_m",
            _locked_member(self.height_m, PHASE3_HEIGHTS_M, "height_m"),
        )
        object.__setattr__(
            self,
            "period_s",
            _locked_member(self.period_s, PHASE3_PERIODS_S, "period_s"),
        )
        object.__setattr__(self, "a", _locked_member(self.a, PHASE3_A_VALUES, "a"))
        object.__setattr__(self, "b", _locked_member(self.b, PHASE3_B_VALUES, "b"))
        object.__setattr__(self, "c", _locked_member(self.c, PHASE3_C_VALUES, "c"))

    @property
    def case_id(self) -> str:
        """Stable human-readable ID derived only from locked parameters."""
        return (
            f"H{_token(self.height_m * 1e3)}mm_"
            f"T{_token(self.period_s)}s_"
            f"A{_token(self.a)}_B{_token(self.b)}_C{_token(self.c)}"
        )

    @property
    def waveform(self) -> BiaxialWallShearWaveform:
        return BiaxialWallShearWaveform(
            shear_scale_pa=PHASE3_TAU_C_PA,
            axial_first_harmonic=self.a,
            axial_second_harmonic=self.b,
            transverse_first_harmonic=self.c,
        )


def locked_phase3_cases() -> tuple[Phase3Case, ...]:
    """Return all 504 states in deterministic protocol order."""
    return tuple(
        Phase3Case(height, period, a, b, c)
        for period, a, b, c, height in product(
            PHASE3_PERIODS_S,
            PHASE3_A_VALUES,
            PHASE3_B_VALUES,
            PHASE3_C_VALUES,
            PHASE3_HEIGHTS_M,
        )
    )


def generate_locked_phase3_state(
    case: Phase3Case,
    *,
    wall_normal_samples: int = PHASE3_PRODUCTION_NY,
    temporal_samples: int = PHASE3_PRODUCTION_NT,
) -> ExactBiaxialChannelState:
    """Generate one exact Phase-3 state on an approved refinement support."""
    if not isinstance(case, Phase3Case):
        raise TypeError("case must be a locked Phase3Case")
    if wall_normal_samples not in PHASE3_NY_LEVELS:
        raise ValueError(
            "wall_normal_samples must use a locked Phase-3 refinement level"
        )
    if temporal_samples not in PHASE3_NT_LEVELS:
        raise ValueError("temporal_samples must use a locked Phase-3 refinement level")
    y = np.linspace(0.0, PHASE3_DELTA_M, wall_normal_samples)
    time = np.linspace(0.0, case.period_s, temporal_samples)
    return exact_biaxial_pulsatile_channel(
        height_m=case.height_m,
        period_s=case.period_s,
        near_wall_thickness_m=PHASE3_DELTA_M,
        density_kg_m3=PHASE3_RHO_KG_M3,
        dynamic_viscosity_pa_s=PHASE3_MU_PA_S,
        waveform=case.waveform,
        scales=PHASE3_SCALES,
        activity_fraction=PHASE3_ACTIVITY_FRACTION,
        y_m=y,
        time_s=time,
    )

"""Small immutable objects that enforce source04 experiment invariants."""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True, slots=True)
class CharacteristicScales:
    """Declared positive characteristic scales used by source04 normalization."""

    velocity_m_s: float
    length_m: float
    density_kg_m3: float
    shear_pa: float

    def __post_init__(self) -> None:
        for name, value in (
            ("velocity_m_s", self.velocity_m_s),
            ("length_m", self.length_m),
            ("density_kg_m3", self.density_kg_m3),
            ("shear_pa", self.shear_pa),
        ):
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and strictly positive")

    @property
    def acceleration_m_s2(self) -> float:
        return self.velocity_m_s**2 / self.length_m

    @property
    def charge_s_2(self) -> float:
        return self.velocity_m_s**2 / self.length_m**2

    @property
    def lambda_pa(self) -> float:
        return self.density_kg_m3 * self.velocity_m_s**2


@dataclass(frozen=True, slots=True)
class StateWeights:
    """Strictly positive weights for [Q*, q_H*, C_U*, C_H*, C_D*]."""

    q_invariant: float
    hydrodynamic_charge: float
    unsteady: float
    head_gradient: float
    dissipative: float

    def __post_init__(self) -> None:
        values = self.as_tuple()
        if any((not math.isfinite(v) or v <= 0.0) for v in values):
            raise ValueError("all source04 state weights must be finite and > 0")
        if not math.isclose(sum(values), 1.0, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError("source04 state weights must sum to 1")

    def as_tuple(self) -> tuple[float, float, float, float, float]:
        return (
            self.q_invariant,
            self.hydrodynamic_charge,
            self.unsteady,
            self.head_gradient,
            self.dissipative,
        )


@dataclass(frozen=True, slots=True)
class Threshold:
    """A non-negative threshold with mandatory independent provenance."""

    name: str
    value: float
    provenance: str

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("threshold name must be non-empty")
        if not math.isfinite(self.value) or self.value < 0.0:
            raise ValueError("threshold value must be finite and non-negative")
        if not self.provenance.strip():
            raise ValueError(
                "threshold provenance is required; arbitrary software defaults are invalid"
            )

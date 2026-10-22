"""Typed source04 primary mechanical-state containers."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .contracts import CharacteristicScales

FloatArray = NDArray[np.float64]


def _scalar_field(value: FloatArray, name: str) -> FloatArray:
    array = np.asarray(value, dtype=float)
    if array.ndim < 1:
        raise ValueError(f"{name} must be an array-valued field")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _vector_field(value: FloatArray, name: str) -> FloatArray:
    array = np.asarray(value, dtype=float)
    if array.ndim < 2 or array.shape[-1] != 3:
        raise ValueError(f"{name} must have trailing vector dimension 3")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


@dataclass(frozen=True, slots=True)
class PrimaryMechanicalState:
    """Dimensional source04 primary state [Q, q_H, C_U, C_H, C_D]."""

    q_invariant_s_2: FloatArray
    hydrodynamic_charge_s_2: FloatArray
    c_unsteady_m_s2: FloatArray
    c_head_m_s2: FloatArray
    c_dissipative_m_s2: FloatArray

    def __post_init__(self) -> None:
        q = _scalar_field(self.q_invariant_s_2, "q_invariant_s_2")
        qh = _scalar_field(self.hydrodynamic_charge_s_2, "hydrodynamic_charge_s_2")
        cu = _vector_field(self.c_unsteady_m_s2, "c_unsteady_m_s2")
        ch = _vector_field(self.c_head_m_s2, "c_head_m_s2")
        cd = _vector_field(self.c_dissipative_m_s2, "c_dissipative_m_s2")
        if q.shape != qh.shape:
            raise ValueError("Q and q_H fields must have identical shapes")
        for name, field in (("C_U", cu), ("C_H", ch), ("C_D", cd)):
            if field.shape[:-1] != q.shape:
                raise ValueError(f"{name} support must match scalar state support")

    def dimensionless(self, scales: CharacteristicScales) -> "DimensionlessState":
        return DimensionlessState(
            q_invariant=np.asarray(self.q_invariant_s_2) / scales.charge_s_2,
            hydrodynamic_charge=(
                np.asarray(self.hydrodynamic_charge_s_2) / scales.charge_s_2
            ),
            c_unsteady=np.asarray(self.c_unsteady_m_s2) / scales.acceleration_m_s2,
            c_head=np.asarray(self.c_head_m_s2) / scales.acceleration_m_s2,
            c_dissipative=(
                np.asarray(self.c_dissipative_m_s2) / scales.acceleration_m_s2
            ),
        )


@dataclass(frozen=True, slots=True)
class DimensionlessState:
    """Dimensionless source04 primary state [Q*, q_H*, C_U*, C_H*, C_D*]."""

    q_invariant: FloatArray
    hydrodynamic_charge: FloatArray
    c_unsteady: FloatArray
    c_head: FloatArray
    c_dissipative: FloatArray

    def __post_init__(self) -> None:
        q = _scalar_field(self.q_invariant, "q_invariant")
        qh = _scalar_field(self.hydrodynamic_charge, "hydrodynamic_charge")
        cu = _vector_field(self.c_unsteady, "c_unsteady")
        ch = _vector_field(self.c_head, "c_head")
        cd = _vector_field(self.c_dissipative, "c_dissipative")
        if q.shape != qh.shape:
            raise ValueError("Q* and q_H* fields must have identical shapes")
        for name, field in (("C_U*", cu), ("C_H*", ch), ("C_D*", cd)):
            if field.shape[:-1] != q.shape:
                raise ValueError(f"{name} support must match scalar state support")

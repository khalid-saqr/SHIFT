"""Source04 closure distances on already-corresponded samples."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .contracts import StateWeights
from .state import DimensionlessState

FloatArray = NDArray[np.float64]


def trapezoid_measure_weights(coordinate: ArrayLike) -> FloatArray:
    """Positive 1-D trapezoid quadrature weights for a strictly increasing grid."""
    x = np.asarray(coordinate, dtype=float)
    if x.ndim != 1 or x.size < 2:
        raise ValueError("coordinate must be a 1-D array with at least two points")
    dx = np.diff(x)
    if not np.all(np.isfinite(x)) or not np.all(dx > 0.0):
        raise ValueError("coordinate must be finite and strictly increasing")
    weights = np.empty_like(x)
    weights[0] = 0.5 * dx[0]
    weights[-1] = 0.5 * dx[-1]
    if x.size > 2:
        weights[1:-1] = 0.5 * (x[2:] - x[:-2])
    return weights


def _validate_weights(weights: ArrayLike, support_shape: tuple[int, ...]) -> FloatArray:
    result = np.asarray(weights, dtype=float)
    if result.shape != support_shape:
        raise ValueError(f"measure weights must have shape {support_shape}")
    if not np.all(np.isfinite(result)) or np.any(result < 0.0):
        raise ValueError("measure weights must be finite and non-negative")
    if float(np.sum(result)) <= 0.0:
        raise ValueError("measure weights must have positive total measure")
    return result


def _weighted_mean(value: FloatArray, weights: FloatArray) -> float:
    return float(np.sum(value * weights) / np.sum(weights))


def scalar_distance(delta: ArrayLike, measure_weights: ArrayLike) -> float:
    """Root-mean-square distance for a scalar field."""
    difference = np.asarray(delta, dtype=float)
    weights = _validate_weights(measure_weights, difference.shape)
    return float(np.sqrt(_weighted_mean(difference**2, weights)))


def vector_distance(delta: ArrayLike, measure_weights: ArrayLike) -> float:
    """Root-mean-square distance for a vector field."""
    difference = np.asarray(delta, dtype=float)
    if difference.ndim < 1 or difference.shape[-1] != 3:
        raise ValueError("vector difference must have trailing dimension 3")
    weights = _validate_weights(measure_weights, difference.shape[:-1])
    squared_norm = np.sum(difference**2, axis=-1)
    return float(np.sqrt(_weighted_mean(squared_norm, weights)))


def wss_distance(
    wss_a_pa: ArrayLike,
    wss_b_pa: ArrayLike,
    *,
    shear_scale_pa: float,
    measure_weights: ArrayLike | None = None,
) -> float:
    """Dimensionless source04 distance for the complete vector WSS field."""
    if shear_scale_pa <= 0.0 or not np.isfinite(shear_scale_pa):
        raise ValueError("shear_scale_pa must be finite and > 0")
    delta = np.asarray(wss_a_pa, dtype=float) - np.asarray(wss_b_pa, dtype=float)
    if delta.shape[-1] != 3:
        raise ValueError("WSS fields must have trailing dimension 3")
    if delta.ndim == 1:
        return float(np.linalg.norm(delta) / shear_scale_pa)
    if measure_weights is None:
        raise ValueError("measure_weights are required for sampled WSS fields")
    return vector_distance(delta / shear_scale_pa, measure_weights)


def state_distance(
    state_a: DimensionlessState,
    state_b: DimensionlessState,
    *,
    state_weights: StateWeights,
    measure_weights: ArrayLike,
) -> float:
    """Primary source04 operational distance d_M."""
    dq = np.asarray(state_a.q_invariant) - np.asarray(state_b.q_invariant)
    dqh = np.asarray(state_a.hydrodynamic_charge) - np.asarray(
        state_b.hydrodynamic_charge
    )
    dcu = np.asarray(state_a.c_unsteady) - np.asarray(state_b.c_unsteady)
    dch = np.asarray(state_a.c_head) - np.asarray(state_b.c_head)
    dcd = np.asarray(state_a.c_dissipative) - np.asarray(state_b.c_dissipative)
    support = dq.shape
    weights = _validate_weights(measure_weights, support)
    if dqh.shape != support:
        raise ValueError("state scalar supports differ")
    for vector in (dcu, dch, dcd):
        if vector.shape != support + (3,):
            raise ValueError("state vector supports differ")
    wq, wh, wu, whead, wd = state_weights.as_tuple()
    integrand = (
        wq * dq**2
        + wh * dqh**2
        + wu * np.sum(dcu**2, axis=-1)
        + whead * np.sum(dch**2, axis=-1)
        + wd * np.sum(dcd**2, axis=-1)
    )
    return float(np.sqrt(_weighted_mean(integrand, weights)))


def lamb_magnitude_distance(
    lamb_a_star: ArrayLike, lamb_b_star: ArrayLike, measure_weights: ArrayLike
) -> float:
    """Nested d_ell based on dimensionless Lamb-vector magnitude."""
    a = np.asarray(lamb_a_star, dtype=float)
    b = np.asarray(lamb_b_star, dtype=float)
    if a.shape != b.shape or a.shape[-1] != 3:
        raise ValueError("Lamb fields must have identical trailing vector dimension 3")
    return scalar_distance(
        np.linalg.norm(a, axis=-1) - np.linalg.norm(b, axis=-1),
        measure_weights,
    )


def composition_distance(
    composition_a: tuple[ArrayLike, ArrayLike, ArrayLike],
    composition_b: tuple[ArrayLike, ArrayLike, ArrayLike],
    *,
    active_a: ArrayLike,
    active_b: ArrayLike,
    measure_weights: ArrayLike,
) -> float | None:
    """Nested d_comp on the common Crocco-active set; None if empty."""
    a = [np.asarray(item, dtype=float) for item in composition_a]
    b = [np.asarray(item, dtype=float) for item in composition_b]
    support = a[0].shape
    if any(item.shape != support for item in (*a, *b)):
        raise ValueError("composition fields must have identical support")
    mask = np.asarray(active_a, dtype=bool) & np.asarray(active_b, dtype=bool)
    if mask.shape != support:
        raise ValueError("activity masks must match composition support")
    if not np.any(mask):
        return None
    weights = _validate_weights(measure_weights, support)[mask]
    squared = sum((ai[mask] - bi[mask]) ** 2 for ai, bi in zip(a, b, strict=True))
    return float(np.sqrt(np.sum(weights * squared) / (3.0 * np.sum(weights))))


def finite_volume_lamb_distance(
    lambda_a_star: ArrayLike,
    lambda_b_star: ArrayLike,
    *,
    time_weights: ArrayLike | None = None,
) -> float:
    """Nested d_Lambda for dimensionless vector finite-volume descriptors."""
    delta = np.asarray(lambda_a_star, dtype=float) - np.asarray(
        lambda_b_star, dtype=float
    )
    if delta.shape[-1] != 3:
        raise ValueError("Lambda descriptors must have trailing dimension 3")
    if delta.ndim == 1:
        return float(np.linalg.norm(delta))
    if time_weights is None:
        raise ValueError("time_weights are required for time-resolved Lambda")
    return vector_distance(delta, time_weights)


def lamb_unsteady_distance(
    balance_a: ArrayLike,
    alignment_a: ArrayLike,
    balance_b: ArrayLike,
    alignment_b: ArrayLike,
    *,
    balance_weight: float,
    alignment_weight: float,
    measure_weights: ArrayLike,
) -> float | None:
    """Nested d_LU on the common set where both B_LU and A_LU are defined."""
    if balance_weight <= 0.0 or alignment_weight <= 0.0:
        raise ValueError("Lamb-unsteady weights must be strictly positive")
    if not np.isclose(balance_weight + alignment_weight, 1.0, atol=1e-12):
        raise ValueError("Lamb-unsteady weights must sum to 1")
    ba = np.asarray(balance_a, dtype=float)
    aa = np.asarray(alignment_a, dtype=float)
    bb = np.asarray(balance_b, dtype=float)
    ab = np.asarray(alignment_b, dtype=float)
    if not (ba.shape == aa.shape == bb.shape == ab.shape):
        raise ValueError("Lamb-unsteady diagnostic supports must match")
    common = np.isfinite(ba) & np.isfinite(aa) & np.isfinite(bb) & np.isfinite(ab)
    if not np.any(common):
        return None
    weights = _validate_weights(measure_weights, ba.shape)[common]
    squared = (
        balance_weight * (ba[common] - bb[common]) ** 2
        + alignment_weight * (aa[common] - ab[common]) ** 2
    )
    return float(np.sqrt(np.sum(weights * squared) / np.sum(weights)))


@dataclass(frozen=True, slots=True)
class DistanceBundle:
    """Primary and nested source04 distances for one corresponded pair."""

    d_tau: float
    d_m: float
    d_q_invariant: float
    d_hydrodynamic_charge: float
    d_lamb: float
    d_composition: float | None
    d_lambda: float
    d_lu: float | None

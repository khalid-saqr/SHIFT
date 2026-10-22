"""Source04 mechanical definitions independent of any flow solver."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]


def _array(value: ArrayLike) -> FloatArray:
    return np.asarray(value, dtype=float)


def vector_magnitude(vector: ArrayLike) -> FloatArray:
    """Euclidean magnitude of a vector field with trailing dimension 3."""
    array = _array(vector)
    if array.shape[-1] != 3:
        raise ValueError("vector field must have trailing dimension 3")
    return np.linalg.norm(array, axis=-1)


def strain_rotation(velocity_gradient_s_1: ArrayLike) -> tuple[FloatArray, FloatArray]:
    """Return symmetric strain and antisymmetric rotation tensors."""
    gradient = _array(velocity_gradient_s_1)
    if gradient.shape[-2:] != (3, 3):
        raise ValueError("velocity gradient must have trailing shape (3, 3)")
    transpose = np.swapaxes(gradient, -1, -2)
    return 0.5 * (gradient + transpose), 0.5 * (gradient - transpose)


def second_invariant_q(velocity_gradient_s_1: ArrayLike) -> FloatArray:
    """Source04 second-invariant strain-rotation descriptor Q [s^-2]."""
    strain, rotation = strain_rotation(velocity_gradient_s_1)
    strain_sq = np.sum(strain * strain, axis=(-2, -1))
    rotation_sq = np.sum(rotation * rotation, axis=(-2, -1))
    return 0.5 * (rotation_sq - strain_sq)


def vorticity_from_gradient(velocity_gradient_s_1: ArrayLike) -> FloatArray:
    """Curl of velocity from A_ij = partial u_i / partial x_j."""
    gradient = _array(velocity_gradient_s_1)
    if gradient.shape[-2:] != (3, 3):
        raise ValueError("velocity gradient must have trailing shape (3, 3)")
    return np.stack(
        (
            gradient[..., 2, 1] - gradient[..., 1, 2],
            gradient[..., 0, 2] - gradient[..., 2, 0],
            gradient[..., 1, 0] - gradient[..., 0, 1],
        ),
        axis=-1,
    )


def lamb_vector(velocity_m_s: ArrayLike, vorticity_s_1: ArrayLike) -> FloatArray:
    """Present-paper Lamb convention ell = u x omega [m/s^2]."""
    return np.cross(_array(velocity_m_s), _array(vorticity_s_1))


def conventional_lamb_vector(
    velocity_m_s: ArrayLike, vorticity_s_1: ArrayLike
) -> FloatArray:
    """Conventional Lamb vector L = omega x u = -ell [m/s^2]."""
    return np.cross(_array(vorticity_s_1), _array(velocity_m_s))


def hydrodynamic_charge_from_identity(
    q_invariant_s_2: ArrayLike, laplacian_kinetic_energy_s_2: ArrayLike
) -> FloatArray:
    """Exact incompressible identity q_H = -2Q - laplacian(K)."""
    return -2.0 * _array(q_invariant_s_2) - _array(laplacian_kinetic_energy_s_2)


def hydrodynamic_charge_residual(
    hydrodynamic_charge_s_2: ArrayLike,
    q_invariant_s_2: ArrayLike,
    laplacian_kinetic_energy_s_2: ArrayLike,
) -> FloatArray:
    """Residual R_q = q_H + 2Q + laplacian(K)."""
    return (
        _array(hydrodynamic_charge_s_2)
        + 2.0 * _array(q_invariant_s_2)
        + _array(laplacian_kinetic_energy_s_2)
    )


def newtonian_extra_stress(
    velocity_gradient_s_1: ArrayLike, viscosity_pa_s: float
) -> FloatArray:
    """Newtonian extra stress tau = 2 mu S [Pa]."""
    if not np.isfinite(viscosity_pa_s) or viscosity_pa_s <= 0.0:
        raise ValueError("viscosity_pa_s must be finite and strictly positive")
    strain, _ = strain_rotation(velocity_gradient_s_1)
    return 2.0 * viscosity_pa_s * strain


def cauchy_stress(pressure_pa: ArrayLike, extra_stress_pa: ArrayLike) -> FloatArray:
    """Cauchy stress sigma = -p I + tau [Pa]."""
    extra = _array(extra_stress_pa)
    if extra.shape[-2:] != (3, 3):
        raise ValueError("extra stress must have trailing shape (3, 3)")
    pressure = _array(pressure_pa)
    if pressure.shape != extra.shape[:-2]:
        try:
            pressure = np.broadcast_to(pressure, extra.shape[:-2])
        except ValueError as exc:
            raise ValueError(
                "pressure support is incompatible with stress support"
            ) from exc
    return extra - pressure[..., None, None] * np.eye(3)


def _validate_unit_normals(normal: ArrayLike, atol: float = 1e-12) -> FloatArray:
    normals = _array(normal)
    if normals.shape[-1] != 3:
        raise ValueError("normal must have trailing dimension 3")
    norms = np.linalg.norm(normals, axis=-1)
    if not np.allclose(norms, 1.0, rtol=0.0, atol=atol):
        raise ValueError("wall normals must be unit vectors")
    return normals


def wall_traction(cauchy_stress_pa: ArrayLike, normal: ArrayLike) -> FloatArray:
    """Full traction t_w = sigma n [Pa]."""
    stress = _array(cauchy_stress_pa)
    if stress.shape[-2:] != (3, 3):
        raise ValueError("Cauchy stress must have trailing shape (3, 3)")
    normals = _validate_unit_normals(normal)
    return np.einsum("...ij,...j->...i", stress, normals)


def wall_shear_stress(cauchy_stress_pa: ArrayLike, normal: ArrayLike) -> FloatArray:
    """Tangential wall traction (I - nn^T) sigma n [Pa]."""
    normals = _validate_unit_normals(normal)
    traction = wall_traction(cauchy_stress_pa, normals)
    normal_component = np.sum(traction * normals, axis=-1, keepdims=True) * normals
    return traction - normal_component


def crocco_residual(
    lamb_m_s2: ArrayLike,
    c_unsteady_m_s2: ArrayLike,
    c_head_m_s2: ArrayLike,
    c_dissipative_m_s2: ArrayLike,
) -> FloatArray:
    """R_CV = ell - C_U - C_H - C_D [m/s^2]."""
    return (
        _array(lamb_m_s2)
        - _array(c_unsteady_m_s2)
        - _array(c_head_m_s2)
        - _array(c_dissipative_m_s2)
    )


def normalized_charge_residual(
    hydrodynamic_charge_s_2: ArrayLike,
    q_invariant_s_2: ArrayLike,
    laplacian_kinetic_energy_s_2: ArrayLike,
) -> FloatArray:
    """Source04 normalized charge residual; undefined where denominator is zero."""
    qh = _array(hydrodynamic_charge_s_2)
    q = _array(q_invariant_s_2)
    lap_k = _array(laplacian_kinetic_energy_s_2)
    numerator = np.abs(qh + 2.0 * q + lap_k)
    denominator = np.abs(qh) + 2.0 * np.abs(q) + np.abs(lap_k)
    result = np.full(np.broadcast_shapes(numerator.shape, denominator.shape), np.nan)
    np.divide(numerator, denominator, out=result, where=denominator > 0.0)
    return result


def normalized_crocco_residual(
    lamb_m_s2: ArrayLike,
    c_unsteady_m_s2: ArrayLike,
    c_head_m_s2: ArrayLike,
    c_dissipative_m_s2: ArrayLike,
) -> FloatArray:
    """Source04 epsilon_CV; undefined where all participating vectors vanish."""
    ell = _array(lamb_m_s2)
    cu = _array(c_unsteady_m_s2)
    ch = _array(c_head_m_s2)
    cd = _array(c_dissipative_m_s2)
    residual = vector_magnitude(crocco_residual(ell, cu, ch, cd))
    denominator = (
        vector_magnitude(ell)
        + vector_magnitude(cu)
        + vector_magnitude(ch)
        + vector_magnitude(cd)
    )
    result = np.full(np.broadcast_shapes(residual.shape, denominator.shape), np.nan)
    np.divide(residual, denominator, out=result, where=denominator > 0.0)
    return result


@dataclass(frozen=True, slots=True)
class CroccoComposition:
    """Magnitude fractions and cancellation on the declared Crocco-active set."""

    d_unsteady: FloatArray
    d_head: FloatArray
    d_dissipative: FloatArray
    cancellation: FloatArray
    active: NDArray[np.bool_]


def crocco_composition(
    lamb_m_s2: ArrayLike,
    c_unsteady_m_s2: ArrayLike,
    c_head_m_s2: ArrayLike,
    c_dissipative_m_s2: ArrayLike,
    *,
    acceleration_scale_m_s2: float,
    activity_fraction: float,
) -> CroccoComposition:
    """Compute D_U, D_H, D_D and chi_C on the source04 active set."""
    if acceleration_scale_m_s2 <= 0.0 or not np.isfinite(acceleration_scale_m_s2):
        raise ValueError("acceleration_scale_m_s2 must be finite and > 0")
    if activity_fraction < 0.0 or not np.isfinite(activity_fraction):
        raise ValueError("activity_fraction must be finite and >= 0")
    ell = _array(lamb_m_s2)
    cu = _array(c_unsteady_m_s2)
    ch = _array(c_head_m_s2)
    cd = _array(c_dissipative_m_s2)
    mu = vector_magnitude(cu)
    mh = vector_magnitude(ch)
    md = vector_magnitude(cd)
    denominator = mu + mh + md
    active = denominator > activity_fraction * acceleration_scale_m_s2
    shape = denominator.shape
    du = np.full(shape, np.nan)
    dh = np.full(shape, np.nan)
    dd = np.full(shape, np.nan)
    cancellation = np.full(shape, np.nan)
    np.divide(mu, denominator, out=du, where=active)
    np.divide(mh, denominator, out=dh, where=active)
    np.divide(md, denominator, out=dd, where=active)
    cancellation[active] = 1.0 - vector_magnitude(ell)[active] / denominator[active]
    return CroccoComposition(du, dh, dd, cancellation, active)


@dataclass(frozen=True, slots=True)
class LambUnsteadyDiagnostics:
    """Source04 nonlinear-unsteady magnitude balance and directional alignment."""

    magnitude_balance: FloatArray
    alignment: FloatArray
    magnitude_active: NDArray[np.bool_]
    alignment_active: NDArray[np.bool_]


def lamb_unsteady_diagnostics(
    lamb_m_s2: ArrayLike,
    c_unsteady_m_s2: ArrayLike,
    *,
    acceleration_scale_m_s2: float,
    activity_fraction: float,
) -> LambUnsteadyDiagnostics:
    """Compute B_LU and A_LU with source04 activity rules."""
    if acceleration_scale_m_s2 <= 0.0 or not np.isfinite(acceleration_scale_m_s2):
        raise ValueError("acceleration_scale_m_s2 must be finite and > 0")
    if activity_fraction < 0.0 or not np.isfinite(activity_fraction):
        raise ValueError("activity_fraction must be finite and >= 0")
    ell = _array(lamb_m_s2)
    cu = _array(c_unsteady_m_s2)
    ell_mag = vector_magnitude(ell)
    cu_mag = vector_magnitude(cu)
    threshold = activity_fraction * acceleration_scale_m_s2
    total = ell_mag + cu_mag
    magnitude_active = total > threshold
    balance = np.full(total.shape, np.nan)
    np.divide(ell_mag - cu_mag, total, out=balance, where=magnitude_active)
    alignment_active = (ell_mag > threshold) & (cu_mag > threshold)
    alignment = np.full(total.shape, np.nan)
    dot = np.sum(ell * cu, axis=-1)
    np.divide(
        dot,
        ell_mag * cu_mag,
        out=alignment,
        where=alignment_active,
    )
    return LambUnsteadyDiagnostics(
        balance, alignment, magnitude_active, alignment_active
    )


def arealized_finite_volume_lamb(
    lamb_m_s2: ArrayLike,
    *,
    density_kg_m3: float,
    volume_weights_m3: ArrayLike,
    wall_area_m2: float,
) -> FloatArray:
    """Lambda_L^(delta) = A^-1 integral rho ell dV [Pa]."""
    if density_kg_m3 <= 0.0 or not np.isfinite(density_kg_m3):
        raise ValueError("density_kg_m3 must be finite and > 0")
    if wall_area_m2 <= 0.0 or not np.isfinite(wall_area_m2):
        raise ValueError("wall_area_m2 must be finite and > 0")
    lamb = _array(lamb_m_s2)
    weights = _array(volume_weights_m3)
    if lamb.shape[-1] != 3 or weights.shape != lamb.shape[:-1]:
        raise ValueError("volume weights must match Lamb-vector sample support")
    integral = np.sum(
        lamb * weights[..., None], axis=tuple(range(weights.ndim))
    )
    return density_kg_m3 * integral / wall_area_m2

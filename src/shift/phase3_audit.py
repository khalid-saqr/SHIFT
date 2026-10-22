"""Predeclared Phase-3 numerical uncertainty audit.

This module executes the numerical audit frozen in ``docs/PHASE3_PROTOCOL.md``.
It does not enumerate WSS pairs or evaluate closure outcomes.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import cache

import numpy as np
from numpy.typing import NDArray

from .certify import IdentityResiduals, identity_residuals
from .mechanics import (
    conventional_lamb_vector,
    lamb_vector,
    second_invariant_q,
    vorticity_from_gradient,
)
from .metrics import (
    finite_volume_lamb_distance,
    lamb_magnitude_distance,
    scalar_distance,
    state_distance,
    trapezoid_measure_weights,
    vector_distance,
    wss_distance,
)
from .operators import differentiate
from .phase3 import (
    PHASE3_NT_LEVELS,
    PHASE3_NY_LEVELS,
    PHASE3_PERIODS_S,
    PHASE3_PRODUCTION_NT,
    PHASE3_PRODUCTION_NY,
    PHASE3_SCALES,
    PHASE3_STATE_WEIGHTS,
    Phase3Case,
    generate_locked_phase3_state,
    locked_phase3_cases,
)
from .state import PrimaryMechanicalState

FloatArray = NDArray[np.float64]

PHASE3_CONVERGENCE_GATE = 1.9

# A convergence order is fitted only against a numerical axis that changes the
# corresponding reconstructed quantity. This separates spatial and temporal
# operator convergence without treating the fixed orthogonal-axis error floor as
# an order-zero failure.
SPATIAL_GATE_QUANTITIES = ("wss", "qH", "CH", "CD", "ell", "Lambda")
TEMPORAL_GATE_QUANTITIES = ("CU",)
STRUCTURAL_ZERO_QUANTITIES = ("Q",)

ERROR_QUANTITIES = (
    "wss",
    "Q",
    "qH",
    "CU",
    "CH",
    "CD",
    "M",
    "ell",
    "Lambda",
    "charge_residual",
    "crocco_residual",
    "incompressibility_residual",
)


@dataclass(frozen=True, slots=True)
class Phase3AuditLevel:
    """Global maximum errors for one locked refinement level."""

    samples: int
    spacing: float
    maxima: dict[str, float]
    worst_case_ids: dict[str, str | None]


@dataclass(frozen=True, slots=True)
class Phase3ConvergenceGate:
    """Observed finest-three-level orders and the predeclared gate decision."""

    spatial_orders: dict[str, float | None]
    temporal_orders: dict[str, float | None]
    spatial_passed: dict[str, bool]
    temporal_passed: dict[str, bool]
    structural_zero_passed: dict[str, bool]

    @property
    def passed(self) -> bool:
        return bool(
            all(self.spatial_passed.values())
            and all(self.temporal_passed.values())
            and all(self.structural_zero_passed.values())
        )


@dataclass(frozen=True, slots=True)
class Phase3AuditResult:
    """Complete pre-pair Phase-3 numerical audit result."""

    state_count: int
    spatial_levels: tuple[Phase3AuditLevel, ...]
    temporal_levels: tuple[Phase3AuditLevel, ...]
    convergence: Phase3ConvergenceGate
    production_one_state_bounds: dict[str, float]
    production_worst_case_ids: dict[str, str | None]

    @property
    def gate_passed(self) -> bool:
        return self.convergence.passed

    def pair_tolerances(self) -> dict[str, float]:
        """Return conservative two-state bounds only after the gate passes."""
        if not self.gate_passed:
            raise RuntimeError(
                "Phase-3 convergence gate failed; pair tolerances must not be frozen"
            )
        return {
            "epsilon_tau": 2.0 * self.production_one_state_bounds["wss"],
            "epsilon_M": 2.0 * self.production_one_state_bounds["M"],
            "epsilon_Q": 2.0 * self.production_one_state_bounds["Q"],
            "epsilon_q": 2.0 * self.production_one_state_bounds["qH"],
            "epsilon_ell": 2.0 * self.production_one_state_bounds["ell"],
            "epsilon_Lambda": 2.0 * self.production_one_state_bounds["Lambda"],
        }

    def as_dict(self) -> dict[str, object]:
        data = asdict(self)
        data["gate_passed"] = self.gate_passed
        data["pair_tolerances_frozen"] = self.gate_passed
        data["pair_tolerances"] = (
            self.pair_tolerances() if self.gate_passed else None
        )
        return data


@dataclass(frozen=True, slots=True)
class _NumericalChannelFields:
    primary_state: PrimaryMechanicalState
    lamb_m_s2: FloatArray
    wall_shear_pa: FloatArray
    lambda_l_pa: FloatArray
    residuals: IdentityResiduals


def _apply_matrix(values: FloatArray, matrix: FloatArray, *, axis: int) -> FloatArray:
    moved = np.moveaxis(np.asarray(values, dtype=float), axis, 0)
    differentiated = np.tensordot(matrix, moved, axes=(1, 0))
    return np.moveaxis(differentiated, 0, axis)


@cache
def _operator_matrix(samples: int, extent: float, derivative_order: int) -> FloatArray:
    coordinate = np.linspace(0.0, extent, samples)
    identity = np.eye(samples, dtype=float)
    return differentiate(
        identity,
        coordinate,
        axis=0,
        derivative_order=derivative_order,
    )


def _reconstruct_numerically(state) -> _NumericalChannelFields:
    """Reconstruct source04 fields using the locked Phase-2 operators."""
    y = state.y_m
    time = state.time_s
    dy = _operator_matrix(y.size, float(y[-1] - y[0]), 1)
    d2y = _operator_matrix(y.size, float(y[-1] - y[0]), 2)
    dt = _operator_matrix(time.size, float(time[-1] - time[0]), 1)

    velocity = state.velocity_m_s
    du_dy = _apply_matrix(velocity, dy, axis=1)
    d2u_dy2 = _apply_matrix(velocity, d2y, axis=1)
    du_dt = _apply_matrix(velocity, dt, axis=0)

    gradient = np.zeros(velocity.shape[:-1] + (3, 3), dtype=float)
    gradient[..., 0, 1] = du_dy[..., 0]
    gradient[..., 2, 1] = du_dy[..., 2]
    divergence = np.trace(gradient, axis1=-2, axis2=-1)
    q_invariant = second_invariant_q(gradient)
    vorticity = vorticity_from_gradient(gradient)
    ell = lamb_vector(velocity, vorticity)
    conventional_lamb = conventional_lamb_vector(velocity, vorticity)

    q_h = _apply_matrix(conventional_lamb[..., 1], dy, axis=1)
    kinetic_energy = 0.5 * np.sum(velocity**2, axis=-1)
    laplacian_k = _apply_matrix(kinetic_energy, d2y, axis=1)
    grad_k_y = _apply_matrix(kinetic_energy, dy, axis=1)

    c_unsteady = du_dt
    c_head = np.zeros_like(velocity)
    pressure_gradient = state.pressure_gradient_pa_m / state.density_kg_m3
    c_head[..., 0] = pressure_gradient[:, None, 0]
    c_head[..., 1] = grad_k_y
    c_head[..., 2] = pressure_gradient[:, None, 2]
    nu = state.dynamic_viscosity_pa_s / state.density_kg_m3
    c_dissipative = -nu * d2u_dy2

    wall_shear = state.dynamic_viscosity_pa_s * du_dy[:, 0, :]
    y_weights = trapezoid_measure_weights(y)
    lambda_l = state.density_kg_m3 * np.tensordot(
        ell,
        y_weights,
        axes=(1, 0),
    )

    primary = PrimaryMechanicalState(
        q_invariant,
        q_h,
        c_unsteady,
        c_head,
        c_dissipative,
    )
    residuals = identity_residuals(
        velocity_divergence_s_1=divergence,
        hydrodynamic_charge_s_2=q_h,
        q_invariant_s_2=q_invariant,
        laplacian_kinetic_energy_s_2=laplacian_k,
        lamb_m_s2=ell,
        c_unsteady_m_s2=c_unsteady,
        c_head_m_s2=c_head,
        c_dissipative_m_s2=c_dissipative,
    )
    return _NumericalChannelFields(
        primary_state=primary,
        lamb_m_s2=ell,
        wall_shear_pa=wall_shear,
        lambda_l_pa=lambda_l,
        residuals=residuals,
    )


def _one_state_errors(state) -> dict[str, float]:
    numerical = _reconstruct_numerically(state)
    y_weights = trapezoid_measure_weights(state.y_m)
    time_weights = trapezoid_measure_weights(state.time_s)
    measure = time_weights[:, None] * y_weights[None, :]
    exact_star = state.primary_state.dimensionless(PHASE3_SCALES)
    numerical_star = numerical.primary_state.dimensionless(PHASE3_SCALES)

    errors = {
        "wss": wss_distance(
            numerical.wall_shear_pa,
            state.wall_shear_exact_pa,
            shear_scale_pa=PHASE3_SCALES.shear_pa,
            measure_weights=time_weights,
        ),
        "Q": scalar_distance(
            numerical_star.q_invariant - exact_star.q_invariant,
            measure,
        ),
        "qH": scalar_distance(
            numerical_star.hydrodynamic_charge
            - exact_star.hydrodynamic_charge,
            measure,
        ),
        "CU": vector_distance(
            numerical_star.c_unsteady - exact_star.c_unsteady,
            measure,
        ),
        "CH": vector_distance(
            numerical_star.c_head - exact_star.c_head,
            measure,
        ),
        "CD": vector_distance(
            numerical_star.c_dissipative - exact_star.c_dissipative,
            measure,
        ),
        "M": state_distance(
            numerical_star,
            exact_star,
            state_weights=PHASE3_STATE_WEIGHTS,
            measure_weights=measure,
        ),
        "ell": lamb_magnitude_distance(
            numerical.lamb_m_s2 / PHASE3_SCALES.acceleration_m_s2,
            state.lamb_m_s2 / PHASE3_SCALES.acceleration_m_s2,
            measure,
        ),
        "Lambda": finite_volume_lamb_distance(
            numerical.lambda_l_pa / PHASE3_SCALES.lambda_pa,
            state.lambda_l_pa / PHASE3_SCALES.lambda_pa,
            time_weights=time_weights,
        ),
        "charge_residual": numerical.residuals.charge_s_2,
        "crocco_residual": numerical.residuals.crocco_m_s2,
        "incompressibility_residual": numerical.residuals.incompressibility_s_1,
    }
    return {name: float(value) for name, value in errors.items()}


def _empty_maxima() -> dict[str, float]:
    return {name: 0.0 for name in ERROR_QUANTITIES}


def _empty_worst() -> dict[str, str | None]:
    return {name: None for name in ERROR_QUANTITIES}


def _update_global(
    maxima: dict[str, float],
    worst: dict[str, str | None],
    case: Phase3Case,
    errors: dict[str, float],
) -> None:
    for name, value in errors.items():
        if value > maxima[name]:
            maxima[name] = value
            worst[name] = case.case_id


def _observed_order(
    levels: tuple[Phase3AuditLevel, ...],
    quantity: str,
) -> float | None:
    finest = levels[-3:]
    errors = np.array([level.maxima[quantity] for level in finest], dtype=float)
    if np.any(errors <= 0.0):
        return None
    spacing = np.array([level.spacing for level in finest], dtype=float)
    slope, _ = np.polyfit(np.log(spacing), np.log(errors), 1)
    return float(slope)


def _spatial_levels(
    cases: tuple[Phase3Case, ...],
) -> tuple[Phase3AuditLevel, ...]:
    levels = []
    for samples in PHASE3_NY_LEVELS:
        maxima = _empty_maxima()
        worst = _empty_worst()
        for case in cases:
            state = generate_locked_phase3_state(
                case,
                wall_normal_samples=samples,
                temporal_samples=PHASE3_PRODUCTION_NT,
            )
            _update_global(maxima, worst, case, _one_state_errors(state))
        levels.append(
            Phase3AuditLevel(
                samples=samples,
                spacing=state.near_wall_thickness_m / (samples - 1),
                maxima=maxima,
                worst_case_ids=worst,
            )
        )
    return tuple(levels)


def _temporal_levels(
    cases: tuple[Phase3Case, ...],
    production_level: Phase3AuditLevel,
) -> tuple[Phase3AuditLevel, ...]:
    levels = []
    max_period = max(PHASE3_PERIODS_S)
    for samples in PHASE3_NT_LEVELS:
        if samples == PHASE3_PRODUCTION_NT:
            levels.append(
                Phase3AuditLevel(
                    samples=samples,
                    spacing=max_period / (samples - 1),
                    maxima=dict(production_level.maxima),
                    worst_case_ids=dict(production_level.worst_case_ids),
                )
            )
            continue

        maxima = _empty_maxima()
        worst = _empty_worst()
        for case in cases:
            state = generate_locked_phase3_state(
                case,
                wall_normal_samples=PHASE3_PRODUCTION_NY,
                temporal_samples=samples,
            )
            _update_global(maxima, worst, case, _one_state_errors(state))
        levels.append(
            Phase3AuditLevel(
                samples=samples,
                spacing=max_period / (samples - 1),
                maxima=maxima,
                worst_case_ids=worst,
            )
        )
    return tuple(levels)


def run_phase3_uncertainty_audit() -> Phase3AuditResult:
    """Run the complete locked 504-state Phase-3 pre-pair uncertainty audit."""
    cases = locked_phase3_cases()
    spatial = _spatial_levels(cases)
    production = spatial[-1]
    temporal = _temporal_levels(cases, production)

    spatial_orders = {
        name: _observed_order(spatial, name) for name in SPATIAL_GATE_QUANTITIES
    }
    temporal_orders = {
        name: _observed_order(temporal, name) for name in TEMPORAL_GATE_QUANTITIES
    }
    spatial_passed = {
        name: order is not None and order >= PHASE3_CONVERGENCE_GATE
        for name, order in spatial_orders.items()
    }
    temporal_passed = {
        name: order is not None and order >= PHASE3_CONVERGENCE_GATE
        for name, order in temporal_orders.items()
    }
    structural_zero_passed = {
        name: production.maxima[name] == 0.0 for name in STRUCTURAL_ZERO_QUANTITIES
    }
    convergence = Phase3ConvergenceGate(
        spatial_orders=spatial_orders,
        temporal_orders=temporal_orders,
        spatial_passed=spatial_passed,
        temporal_passed=temporal_passed,
        structural_zero_passed=structural_zero_passed,
    )
    return Phase3AuditResult(
        state_count=len(cases),
        spatial_levels=spatial,
        temporal_levels=temporal,
        convergence=convergence,
        production_one_state_bounds=dict(production.maxima),
        production_worst_case_ids=dict(production.worst_case_ids),
    )

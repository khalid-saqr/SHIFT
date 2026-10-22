"""Locked analytical regression benchmarks for the source04 theory kernel."""

from dataclasses import asdict, dataclass

import numpy as np

from .certify import identity_residuals
from .contracts import CharacteristicScales, StateWeights
from .metrics import (
    DistanceBundle,
    composition_distance,
    finite_volume_lamb_distance,
    lamb_magnitude_distance,
    lamb_unsteady_distance,
    scalar_distance,
    state_distance,
    trapezoid_measure_weights,
    wss_distance,
)
from .models.poiseuille import exact_parallel_plate_poiseuille

RHO_KG_M3 = 1050.0
MU_PA_S = 0.004
TAU0_PA = 1.20
HEIGHT_A_M = 0.004
DELTA_M = 0.0004
SHEAR_RATE_S_1 = TAU0_PA / MU_PA_S
U_C_M_S = SHEAR_RATE_S_1 * DELTA_M
SCALES = CharacteristicScales(U_C_M_S, DELTA_M, RHO_KG_M3, TAU0_PA)
STATE_WEIGHTS = StateWeights(0.2, 0.2, 0.2, 0.2, 0.2)
ACTIVITY_FRACTION = 1e-3


@dataclass(frozen=True, slots=True)
class PoiseuillePairResult:
    """Locked representative exact-pair result."""

    height_a_m: float
    height_b_m: float
    distances: DistanceBundle
    max_incompressibility_residual_s_1: float
    max_charge_residual_s_2: float
    max_crocco_residual_m_s2: float
    max_lambda_quadrature_error_pa: float

    def as_dict(self) -> dict[str, object]:
        data = asdict(self)
        data["distances"] = asdict(self.distances)
        return data


def compare_poiseuille_heights(
    height_a_m: float,
    height_b_m: float,
    *,
    sample_count: int = 4001,
) -> PoiseuillePairResult:
    """Evaluate the complete source04 hierarchy for two exact same-WSS states."""
    common = dict(
        near_wall_thickness_m=DELTA_M,
        density_kg_m3=RHO_KG_M3,
        viscosity_pa_s=MU_PA_S,
        wall_shear_pa=TAU0_PA,
        scales=SCALES,
        activity_fraction=ACTIVITY_FRACTION,
        sample_count=sample_count,
    )
    a = exact_parallel_plate_poiseuille(height_m=height_a_m, **common)
    b = exact_parallel_plate_poiseuille(height_m=height_b_m, **common)
    if not np.array_equal(a.y_m, b.y_m):
        raise RuntimeError("exact benchmark states do not share the same support")
    measure = trapezoid_measure_weights(a.y_m)
    a_star = a.primary_state.dimensionless(SCALES)
    b_star = b.primary_state.dimensionless(SCALES)

    d_tau = wss_distance(
        a.wall_shear_pa,
        b.wall_shear_pa,
        shear_scale_pa=SCALES.shear_pa,
    )
    d_m = state_distance(
        a_star,
        b_star,
        state_weights=STATE_WEIGHTS,
        measure_weights=measure,
    )
    d_q_invariant = scalar_distance(
        a_star.q_invariant - b_star.q_invariant, measure
    )
    d_hydrodynamic_charge = scalar_distance(
        a_star.hydrodynamic_charge - b_star.hydrodynamic_charge, measure
    )
    d_lamb = lamb_magnitude_distance(
        a.lamb_m_s2 / SCALES.acceleration_m_s2,
        b.lamb_m_s2 / SCALES.acceleration_m_s2,
        measure,
    )
    d_composition = composition_distance(
        (
            a.composition.d_unsteady,
            a.composition.d_head,
            a.composition.d_dissipative,
        ),
        (
            b.composition.d_unsteady,
            b.composition.d_head,
            b.composition.d_dissipative,
        ),
        active_a=a.composition.active,
        active_b=b.composition.active,
        measure_weights=measure,
    )
    d_lambda = finite_volume_lamb_distance(
        a.lambda_l_pa / SCALES.lambda_pa,
        b.lambda_l_pa / SCALES.lambda_pa,
    )
    d_lu = lamb_unsteady_distance(
        a.lamb_unsteady.magnitude_balance,
        a.lamb_unsteady.alignment,
        b.lamb_unsteady.magnitude_balance,
        b.lamb_unsteady.alignment,
        balance_weight=0.5,
        alignment_weight=0.5,
        measure_weights=measure,
    )

    residual_a = identity_residuals(
        velocity_divergence_s_1=a.velocity_divergence_s_1,
        hydrodynamic_charge_s_2=a.hydrodynamic_charge_s_2,
        q_invariant_s_2=a.q_invariant_s_2,
        laplacian_kinetic_energy_s_2=a.laplacian_kinetic_energy_s_2,
        lamb_m_s2=a.lamb_m_s2,
        c_unsteady_m_s2=a.c_unsteady_m_s2,
        c_head_m_s2=a.c_head_m_s2,
        c_dissipative_m_s2=a.c_dissipative_m_s2,
    )
    residual_b = identity_residuals(
        velocity_divergence_s_1=b.velocity_divergence_s_1,
        hydrodynamic_charge_s_2=b.hydrodynamic_charge_s_2,
        q_invariant_s_2=b.q_invariant_s_2,
        laplacian_kinetic_energy_s_2=b.laplacian_kinetic_energy_s_2,
        lamb_m_s2=b.lamb_m_s2,
        c_unsteady_m_s2=b.c_unsteady_m_s2,
        c_head_m_s2=b.c_head_m_s2,
        c_dissipative_m_s2=b.c_dissipative_m_s2,
    )
    lambda_error = max(
        float(np.max(np.abs(a.lambda_l_pa - a.lambda_l_exact_pa))),
        float(np.max(np.abs(b.lambda_l_pa - b.lambda_l_exact_pa))),
    )

    return PoiseuillePairResult(
        height_a_m=float(height_a_m),
        height_b_m=float(height_b_m),
        distances=DistanceBundle(
            d_tau=d_tau,
            d_m=d_m,
            d_q_invariant=d_q_invariant,
            d_hydrodynamic_charge=d_hydrodynamic_charge,
            d_lamb=d_lamb,
            d_composition=d_composition,
            d_lambda=d_lambda,
            d_lu=d_lu,
        ),
        max_incompressibility_residual_s_1=max(
            residual_a.incompressibility_s_1,
            residual_b.incompressibility_s_1,
        ),
        max_charge_residual_s_2=max(residual_a.charge_s_2, residual_b.charge_s_2),
        max_crocco_residual_m_s2=max(
            residual_a.crocco_m_s2, residual_b.crocco_m_s2
        ),
        max_lambda_quadrature_error_pa=lambda_error,
    )


def representative_poiseuille_pair() -> PoiseuillePairResult:
    """Return the locked H_B = 2 H_A source04 analytical regression pair."""
    return compare_poiseuille_heights(HEIGHT_A_M, 2.0 * HEIGHT_A_M)


def poiseuille_sweep(
    *, gamma_min: float = 0.5, gamma_max: float = 3.0, cases: int = 101
) -> list[PoiseuillePairResult]:
    """Locked continuous same-WSS height-ratio family."""
    if gamma_min <= 0.0 or gamma_max < gamma_min or cases < 2:
        raise ValueError("invalid gamma sweep definition")
    return [
        compare_poiseuille_heights(HEIGHT_A_M, gamma * HEIGHT_A_M)
        for gamma in np.linspace(gamma_min, gamma_max, cases)
    ]

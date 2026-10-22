import numpy as np
import pytest

from shift.models.biaxial_channel import (
    BiaxialWallShearWaveform,
    exact_biaxial_pulsatile_channel,
)
from shift.models.poiseuille import exact_parallel_plate_poiseuille
from shift.phase3 import (
    PHASE3_ACTIVITY_FRACTION,
    PHASE3_DELTA_M,
    PHASE3_MU_PA_S,
    PHASE3_RHO_KG_M3,
    PHASE3_SCALES,
    PHASE3_TAU_C_PA,
    Phase3Case,
    generate_locked_phase3_state,
    locked_phase3_cases,
)


def _unsteady_state():
    height = 0.008
    period = 0.8
    y = np.unique(
        np.concatenate(
            (np.linspace(0.0, height, 65), np.array([PHASE3_DELTA_M]))
        )
    )
    time = np.linspace(0.0, period, 129)
    waveform = BiaxialWallShearWaveform(PHASE3_TAU_C_PA, 0.75, 0.3, 0.6)
    return exact_biaxial_pulsatile_channel(
        height_m=height,
        period_s=period,
        near_wall_thickness_m=PHASE3_DELTA_M,
        density_kg_m3=PHASE3_RHO_KG_M3,
        dynamic_viscosity_pa_s=PHASE3_MU_PA_S,
        waveform=waveform,
        scales=PHASE3_SCALES,
        activity_fraction=PHASE3_ACTIVITY_FRACTION,
        y_m=y,
        time_s=time,
    )


def test_exact_unsteady_state_satisfies_boundaries_periodicity_and_wss() -> None:
    state = _unsteady_state()

    assert np.max(np.abs(state.velocity_m_s[:, 0, :])) <= 1e-13
    assert np.max(np.abs(state.velocity_m_s[:, -1, :])) <= 1e-13
    assert np.max(np.abs(state.velocity_m_s[0] - state.velocity_m_s[-1])) <= 1e-13
    pressure_cycle_gap = (
        state.pressure_drop_gradient_pa_m[0]
        - state.pressure_drop_gradient_pa_m[-1]
    )
    assert np.max(np.abs(pressure_cycle_gap)) <= 1e-10
    assert np.max(
        np.abs(state.wall_shear_exact_pa - state.wall_shear_target_pa)
    ) <= 2e-12


def test_exact_unsteady_state_closes_momentum_and_source04_identities() -> None:
    state = _unsteady_state()
    nu = PHASE3_MU_PA_S / PHASE3_RHO_KG_M3
    rhs = (
        state.pressure_drop_gradient_pa_m[:, None, :] / PHASE3_RHO_KG_M3
        + nu * state.d2u_dy2_m_1_s
    )

    assert np.max(np.abs(state.du_dt_m_s2 - rhs)) <= 2e-12
    assert np.max(np.abs(state.q_invariant_s_2)) <= 1e-12
    assert state.residuals.incompressibility_s_1 <= 1e-12
    assert state.residuals.charge_s_2 <= 1e-12
    assert state.residuals.crocco_m_s2 <= 2e-12


def test_exact_finite_volume_lamb_reference_is_rho_k_at_delta() -> None:
    state = _unsteady_state()
    delta_index = int(np.flatnonzero(np.isclose(state.y_m, PHASE3_DELTA_M))[0])
    expected = np.zeros_like(state.lambda_l_pa)
    expected[:, 1] = (
        PHASE3_RHO_KG_M3 * state.kinetic_energy_m2_s2[:, delta_index]
    )
    assert np.allclose(state.lambda_l_pa, expected, rtol=0.0, atol=1e-12)


def test_pressure_evaluator_has_the_declared_exact_gradient() -> None:
    state = _unsteady_state()
    dx = 1e-5
    dz = 2e-5
    p0 = state.pressure_pa(0.0, 0.0)
    px = state.pressure_pa(dx, 0.0)
    pz = state.pressure_pa(0.0, dz)

    dpdx = (px - p0) / dx
    dpdz = (pz - p0) / dz
    assert np.allclose(
        dpdx, state.pressure_gradient_pa_m[:, 0], rtol=1e-12, atol=1e-10
    )
    assert np.allclose(
        dpdz, state.pressure_gradient_pa_m[:, 2], rtol=1e-12, atol=1e-10
    )


def test_exact_model_requires_explicit_lower_wall_support() -> None:
    waveform = BiaxialWallShearWaveform(PHASE3_TAU_C_PA, 0.75, 0.3, 0.6)
    with pytest.raises(ValueError, match="lower wall"):
        exact_biaxial_pulsatile_channel(
            height_m=0.008,
            period_s=0.8,
            near_wall_thickness_m=PHASE3_DELTA_M,
            density_kg_m3=PHASE3_RHO_KG_M3,
            dynamic_viscosity_pa_s=PHASE3_MU_PA_S,
            waveform=waveform,
            scales=PHASE3_SCALES,
            activity_fraction=PHASE3_ACTIVITY_FRACTION,
            y_m=np.array([1e-4, 2e-4]),
            time_s=np.array([0.0, 0.8]),
        )


def test_steady_member_reduces_to_locked_phase1_poiseuille_solution() -> None:
    height = 0.008
    y = np.linspace(0.0, PHASE3_DELTA_M, 33)
    time = np.array([0.0, 0.8])
    steady = exact_biaxial_pulsatile_channel(
        height_m=height,
        period_s=0.8,
        near_wall_thickness_m=PHASE3_DELTA_M,
        density_kg_m3=PHASE3_RHO_KG_M3,
        dynamic_viscosity_pa_s=PHASE3_MU_PA_S,
        waveform=BiaxialWallShearWaveform(PHASE3_TAU_C_PA, 0.0, 0.0, 0.0),
        scales=PHASE3_SCALES,
        activity_fraction=PHASE3_ACTIVITY_FRACTION,
        y_m=y,
        time_s=time,
    )
    phase1 = exact_parallel_plate_poiseuille(
        height_m=height,
        near_wall_thickness_m=PHASE3_DELTA_M,
        density_kg_m3=PHASE3_RHO_KG_M3,
        viscosity_pa_s=PHASE3_MU_PA_S,
        wall_shear_pa=PHASE3_TAU_C_PA,
        scales=PHASE3_SCALES,
        activity_fraction=PHASE3_ACTIVITY_FRACTION,
        sample_count=33,
    )

    assert np.allclose(steady.velocity_m_s[0], phase1.velocity_m_s, atol=1e-14)
    assert np.allclose(
        steady.velocity_gradient_s_1[0], phase1.velocity_gradient_s_1, atol=1e-12
    )
    assert np.allclose(
        steady.hydrodynamic_charge_s_2[0],
        phase1.hydrodynamic_charge_s_2,
        atol=1e-10,
    )
    assert np.allclose(steady.c_head_m_s2[0], phase1.c_head_m_s2, atol=1e-12)
    assert np.allclose(
        steady.c_dissipative_m_s2[0], phase1.c_dissipative_m_s2, atol=1e-12
    )
    assert np.allclose(steady.lambda_l_pa[0], phase1.lambda_l_exact_pa, atol=1e-12)
    assert np.allclose(steady.wall_shear_exact_pa[0], phase1.wall_shear_pa, atol=1e-12)


def test_locked_phase3_family_contains_exactly_504_unique_states() -> None:
    cases = locked_phase3_cases()
    assert len(cases) == 504
    assert len({case.case_id for case in cases}) == 504
    assert len({case.height_m for case in cases}) == 7
    assert len({case.period_s for case in cases}) == 4
    assert len({(case.a, case.b, case.c) for case in cases}) == 18


def test_locked_generator_uses_only_protocol_grid_and_refinement_levels() -> None:
    case = Phase3Case(0.016, 0.6, 1.5, 0.3, 1.2)
    state = generate_locked_phase3_state(
        case,
        wall_normal_samples=17,
        temporal_samples=32,
    )
    assert state.velocity_m_s.shape == (32, 17, 3)
    assert np.isclose(state.y_m[-1], PHASE3_DELTA_M)
    assert np.isclose(state.time_s[-1], case.period_s)
    assert np.max(
        np.abs(state.wall_shear_exact_pa - state.wall_shear_target_pa)
    ) <= 2e-12

    with pytest.raises(ValueError):
        Phase3Case(0.005, 0.6, 1.5, 0.3, 1.2)
    with pytest.raises(ValueError):
        generate_locked_phase3_state(case, wall_normal_samples=31, temporal_samples=32)
    with pytest.raises(ValueError):
        generate_locked_phase3_state(case, wall_normal_samples=17, temporal_samples=100)

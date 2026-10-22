import numpy as np

from shift.benchmarks import (
    HEIGHT_A_M,
    poiseuille_sweep,
    representative_poiseuille_pair,
)
from shift.models.poiseuille import closed_form_d_q


def test_representative_pair_matches_locked_source04_values() -> None:
    result = representative_poiseuille_pair()
    d = result.distances

    assert d.d_tau == 0.0
    assert d.d_q_invariant == 0.0
    assert np.isclose(
        d.d_hydrodynamic_charge, 0.15380182234880216, rtol=2e-8, atol=1e-12
    )
    assert np.isclose(d.d_lamb, 0.061499131829129146, rtol=2e-8, atol=1e-12)
    assert np.isclose(
        d.d_composition, 0.030259329233964294, rtol=2e-8, atol=1e-12
    )
    assert np.isclose(d.d_lambda, 0.0462500013281249, rtol=2e-8, atol=1e-12)
    assert np.isclose(d.d_m, 0.07424704720054744, rtol=2e-8, atol=1e-12)
    assert d.d_lu is None

    assert result.max_incompressibility_residual_s_1 <= 1e-12
    assert result.max_charge_residual_s_2 <= 1e-12
    assert result.max_crocco_residual_m_s2 <= 1e-12
    assert result.max_lambda_quadrature_error_pa <= 1e-7


def test_dq_matches_independent_closed_form_expression() -> None:
    result = representative_poiseuille_pair()
    closed = closed_form_d_q(
        height_a_m=HEIGHT_A_M,
        height_b_m=2.0 * HEIGHT_A_M,
        near_wall_thickness_m=0.0004,
    )
    assert np.isclose(
        result.distances.d_hydrodynamic_charge,
        closed,
        rtol=2e-8,
        atol=1e-12,
    )


def test_locked_height_ratio_family_is_exactly_wss_equivalent() -> None:
    sweep = poiseuille_sweep()
    assert len(sweep) == 101
    assert all(item.distances.d_tau == 0.0 for item in sweep)
    assert all(item.max_incompressibility_residual_s_1 <= 1e-12 for item in sweep)
    assert all(item.max_charge_residual_s_2 <= 1e-12 for item in sweep)
    assert all(item.max_crocco_residual_m_s2 <= 1e-12 for item in sweep)

    gamma_one = sweep[20]
    assert np.isclose(gamma_one.height_b_m / gamma_one.height_a_m, 1.0)
    assert gamma_one.distances.d_m < 1e-13
    assert gamma_one.distances.d_hydrodynamic_charge < 1e-13
    assert gamma_one.distances.d_lambda < 1e-13

    nontrivial = [item for i, item in enumerate(sweep) if i != 20]
    assert all(item.distances.d_m > 0.0 for item in nontrivial)


def test_representative_pair_wss_and_state_separation_are_not_threshold_tuned() -> None:
    result = representative_poiseuille_pair()
    assert result.distances.d_tau == 0.0
    assert result.distances.d_m > 0.07

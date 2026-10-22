import numpy as np
import pytest

from shift.phase3_audit import run_phase3_uncertainty_audit


def test_full_locked_phase3_audit_reproduces_pre_pair_hard_stop() -> None:
    result = run_phase3_uncertainty_audit()

    assert result.state_count == 504
    assert not result.gate_passed
    assert result.convergence.spatial_passed["wss"]
    assert not result.convergence.spatial_passed["qH"]
    assert result.convergence.spatial_passed["CH"]
    assert result.convergence.spatial_passed["CD"]
    assert result.convergence.spatial_passed["ell"]
    assert result.convergence.spatial_passed["Lambda"]
    assert result.convergence.temporal_passed["CU"]
    assert result.convergence.structural_zero_passed["Q"]

    assert np.isclose(
        result.convergence.spatial_orders["qH"],
        1.5339639140035293,
        rtol=2e-5,
        atol=0.0,
    )
    assert np.isclose(
        result.convergence.spatial_orders["wss"],
        1.9971243794196796,
        rtol=2e-5,
        atol=0.0,
    )
    assert np.isclose(
        result.convergence.temporal_orders["CU"],
        2.0187018765118663,
        rtol=2e-5,
        atol=0.0,
    )
    assert np.isclose(
        result.production_one_state_bounds["qH"],
        5.4964502277711506e-05,
        rtol=2e-5,
        atol=0.0,
    )
    with pytest.raises(RuntimeError, match="must not be frozen"):
        result.pair_tolerances()

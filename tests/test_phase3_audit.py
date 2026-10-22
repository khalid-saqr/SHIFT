import numpy as np
import pytest

from shift.phase3 import Phase3Case, generate_locked_phase3_state
from shift.phase3_audit import (
    PHASE3_CONVERGENCE_GATE,
    Phase3AuditLevel,
    Phase3AuditResult,
    Phase3ConvergenceGate,
    _one_state_errors,
)


def test_phase3_case_audit_uses_numerical_reconstruction() -> None:
    case = Phase3Case(0.008, 0.8, 0.75, 0.3, 0.6)
    state = generate_locked_phase3_state(
        case,
        wall_normal_samples=65,
        temporal_samples=128,
    )
    errors = _one_state_errors(state)

    assert errors["Q"] == 0.0
    assert errors["wss"] > 0.0
    assert errors["qH"] > 0.0
    assert errors["CU"] > 0.0
    assert errors["CH"] > 0.0
    assert errors["CD"] > 0.0
    assert errors["M"] > 0.0
    assert errors["charge_residual"] > 0.0
    assert errors["crocco_residual"] > 0.0
    assert errors["incompressibility_residual"] == 0.0


def test_locked_qh_defining_route_does_not_meet_phase3_spatial_gate() -> None:
    case = Phase3Case(0.016, 0.6, 1.5, 0.3, 0.0)
    samples = np.array([65, 129, 257], dtype=int)
    spacing = 0.0004 / (samples - 1)
    errors = []
    for count in samples:
        state = generate_locked_phase3_state(
            case,
            wall_normal_samples=int(count),
            temporal_samples=512,
        )
        errors.append(_one_state_errors(state)["qH"])

    order, _ = np.polyfit(np.log(spacing), np.log(errors), 1)
    assert order < PHASE3_CONVERGENCE_GATE
    assert 1.45 < order < 1.65


def test_failed_gate_refuses_to_issue_pair_tolerances() -> None:
    level = Phase3AuditLevel(
        samples=257,
        spacing=0.0004 / 256,
        maxima={
            "wss": 1e-6,
            "Q": 0.0,
            "qH": 1e-5,
            "CU": 1e-6,
            "CH": 1e-6,
            "CD": 1e-8,
            "M": 1e-5,
            "ell": 1e-6,
            "Lambda": 1e-6,
            "charge_residual": 1.0,
            "crocco_residual": 1e-3,
            "incompressibility_residual": 0.0,
        },
        worst_case_ids={},
    )
    convergence = Phase3ConvergenceGate(
        spatial_orders={"qH": 1.5},
        temporal_orders={"CU": 2.0},
        spatial_passed={"qH": False},
        temporal_passed={"CU": True},
        structural_zero_passed={"Q": True},
    )
    result = Phase3AuditResult(
        state_count=504,
        spatial_levels=(level,),
        temporal_levels=(level,),
        convergence=convergence,
        production_one_state_bounds=level.maxima,
        production_worst_case_ids=level.worst_case_ids,
    )

    assert not result.gate_passed
    with pytest.raises(RuntimeError, match="must not be frozen"):
        result.pair_tolerances()
    assert result.as_dict()["pair_tolerances"] is None

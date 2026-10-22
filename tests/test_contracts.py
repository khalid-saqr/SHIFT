import math

import pytest

from shift.contracts import CharacteristicScales, StateWeights, Threshold


def test_characteristic_scales_match_source04_definitions() -> None:
    scales = CharacteristicScales(
        velocity_m_s=0.12,
        length_m=4.0e-4,
        density_kg_m3=1050.0,
        shear_pa=1.2,
    )
    assert math.isclose(scales.acceleration_m_s2, 36.0)
    assert math.isclose(scales.charge_s_2, 90000.0)
    assert math.isclose(scales.lambda_pa, 15.12)


def test_characteristic_scales_reject_nonpositive_values() -> None:
    with pytest.raises(ValueError):
        CharacteristicScales(0.0, 1.0, 1.0, 1.0)


def test_state_weights_require_strict_positivity_and_unit_sum() -> None:
    weights = StateWeights(0.2, 0.2, 0.2, 0.2, 0.2)
    assert math.isclose(sum(weights.as_tuple()), 1.0)

    with pytest.raises(ValueError):
        StateWeights(0.0, 0.25, 0.25, 0.25, 0.25)

    with pytest.raises(ValueError):
        StateWeights(0.1, 0.1, 0.1, 0.1, 0.1)


def test_threshold_requires_independent_provenance() -> None:
    threshold = Threshold("epsilon_tau", 0.01, "grid-convergence study v1")
    assert threshold.value == 0.01

    with pytest.raises(ValueError):
        Threshold("epsilon_tau", 0.01, "")

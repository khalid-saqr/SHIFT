import numpy as np

from shift.contracts import StateWeights
from shift.metrics import (
    composition_distance,
    lamb_unsteady_distance,
    scalar_distance,
    state_distance,
    trapezoid_measure_weights,
)
from shift.state import DimensionlessState


def test_trapezoid_weights_integrate_linear_function_exactly() -> None:
    x = np.array([0.0, 0.1, 0.4, 1.0])
    weights = trapezoid_measure_weights(x)
    integral = float(np.sum(weights * x))
    assert np.isclose(integral, 0.5)


def test_state_distance_uses_declared_positive_weights() -> None:
    support = (3,)
    zeros = np.zeros(support)
    zero_vec = np.zeros(support + (3,))
    unit = np.ones(support)
    unit_vec = np.zeros(support + (3,))
    unit_vec[:, 0] = 1.0
    a = DimensionlessState(zeros, zeros, zero_vec, zero_vec, zero_vec)
    b = DimensionlessState(unit, unit, unit_vec, unit_vec, unit_vec)
    weights = StateWeights(0.2, 0.2, 0.2, 0.2, 0.2)
    measure = np.ones(3)

    distance = state_distance(
        a, b, state_weights=weights, measure_weights=measure
    )
    assert np.isclose(distance, 1.0)
    assert np.isclose(scalar_distance(unit, measure), 1.0)


def test_empty_activity_sets_remain_undefined() -> None:
    x = np.array([np.nan, np.nan])
    measure = np.ones(2)
    assert (
        lamb_unsteady_distance(
            x,
            x,
            x,
            x,
            balance_weight=0.5,
            alignment_weight=0.5,
            measure_weights=measure,
        )
        is None
    )

    zeros = np.zeros(2)
    inactive = np.zeros(2, dtype=bool)
    assert (
        composition_distance(
            (zeros, zeros, zeros),
            (zeros, zeros, zeros),
            active_a=inactive,
            active_b=inactive,
            measure_weights=measure,
        )
        is None
    )

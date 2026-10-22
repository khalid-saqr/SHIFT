import numpy as np

from shift.certify import (
    IdentityResiduals,
    IdentityThresholds,
    qualify_identity_residuals,
)
from shift.uncertainty import (
    ConvergenceSeries,
    UncertaintyBudget,
    UncertaintyComponent,
    richardson_uncertainty,
    threshold_from_budget,
)


def test_convergence_series_recovers_known_order() -> None:
    h = np.array([0.2, 0.1, 0.05, 0.025])
    error = 3.0 * h**2
    series = ConvergenceSeries(h, error, "manufactured derivative", "s^-1")
    assert np.isclose(series.observed_order(), 2.0, atol=1e-12)


def test_richardson_uncertainty_recovers_fine_grid_error() -> None:
    truth = 7.0
    coefficient = 4.0
    coarse_h = 0.1
    fine_h = 0.05
    coarse = truth + coefficient * coarse_h**2
    fine = truth + coefficient * fine_h**2
    estimate = richardson_uncertainty(
        coarse_value=coarse,
        fine_value=fine,
        coarse_spacing=coarse_h,
        fine_spacing=fine_h,
        verified_order=2.0,
        name="spatial discretization",
        units="s^-1",
        provenance="second-order operator verified on manufactured fields",
    )
    assert np.isclose(estimate.uncertainty.bound, coefficient * fine_h**2)


def test_budget_generates_threshold_with_complete_provenance() -> None:
    budget = UncertaintyBudget(
        (
            UncertaintyComponent("spatial", 0.03, "s^-1", "grid refinement A"),
            UncertaintyComponent("temporal", 0.01, "s^-1", "time refinement B"),
        )
    )
    threshold = threshold_from_budget("continuity", budget)
    assert threshold.value == 0.04
    assert "spatial" in threshold.provenance
    assert "temporal" in threshold.provenance


def test_qualification_uses_supplied_uncertainty_thresholds_exactly() -> None:
    continuity = threshold_from_budget(
        "continuity",
        UncertaintyBudget((UncertaintyComponent("grid", 0.1, "s^-1", "study"),)),
    )
    charge = threshold_from_budget(
        "charge",
        UncertaintyBudget((UncertaintyComponent("grid", 2.0, "s^-2", "study"),)),
    )
    crocco = threshold_from_budget(
        "crocco",
        UncertaintyBudget((UncertaintyComponent("grid", 0.5, "m/s^2", "study"),)),
    )
    thresholds = IdentityThresholds(continuity, charge, crocco)
    certificate = qualify_identity_residuals(
        IdentityResiduals(0.1, 2.0001, 0.49), thresholds
    )
    assert certificate.incompressibility_passed
    assert not certificate.charge_passed
    assert certificate.crocco_passed
    assert not certificate.passed

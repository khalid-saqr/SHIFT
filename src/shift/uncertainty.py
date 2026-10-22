"""Outcome-independent convergence and numerical-uncertainty utilities."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .contracts import Threshold

FloatArray = NDArray[np.float64]


def rms_error(numerical: ArrayLike, reference: ArrayLike) -> float:
    """Root-mean-square absolute error on identical sample support."""
    numerical_array = np.asarray(numerical, dtype=float)
    reference_array = np.asarray(reference, dtype=float)
    if numerical_array.shape != reference_array.shape:
        raise ValueError("numerical and reference arrays must have identical shapes")
    difference = numerical_array - reference_array
    if not np.all(np.isfinite(difference)):
        raise ValueError("error inputs must be finite")
    return float(np.sqrt(np.mean(difference**2)))


def max_error(numerical: ArrayLike, reference: ArrayLike) -> float:
    """Maximum absolute error on identical sample support."""
    numerical_array = np.asarray(numerical, dtype=float)
    reference_array = np.asarray(reference, dtype=float)
    if numerical_array.shape != reference_array.shape:
        raise ValueError("numerical and reference arrays must have identical shapes")
    difference = numerical_array - reference_array
    if not np.all(np.isfinite(difference)):
        raise ValueError("error inputs must be finite")
    return float(np.max(np.abs(difference)))


@dataclass(frozen=True, slots=True)
class ConvergenceSeries:
    """Errors from a systematically refined sequence with known reference truth."""

    spacing: FloatArray
    error: FloatArray
    quantity: str
    units: str

    def __post_init__(self) -> None:
        spacing = np.asarray(self.spacing, dtype=float)
        error = np.asarray(self.error, dtype=float)
        if spacing.ndim != 1 or error.ndim != 1 or spacing.size != error.size:
            raise ValueError("spacing and error must be equal-length 1-D arrays")
        if spacing.size < 3:
            raise ValueError("at least three refinement levels are required")
        if not np.all(np.isfinite(spacing)) or np.any(spacing <= 0.0):
            raise ValueError("spacing values must be finite and strictly positive")
        if not np.all(np.diff(spacing) < 0.0):
            raise ValueError("spacing must be strictly decreasing from coarse to fine")
        if not np.all(np.isfinite(error)) or np.any(error < 0.0):
            raise ValueError("errors must be finite and non-negative")
        if not self.quantity.strip() or not self.units.strip():
            raise ValueError("quantity and units must be non-empty")
        object.__setattr__(self, "spacing", spacing)
        object.__setattr__(self, "error", error)

    def observed_order(self, *, finest_levels: int | None = None) -> float:
        """Log-log convergence slope p for error proportional to h**p."""
        h = self.spacing
        e = self.error
        if finest_levels is not None:
            if finest_levels < 2 or finest_levels > h.size:
                raise ValueError("finest_levels must be between 2 and series length")
            h = h[-finest_levels:]
            e = e[-finest_levels:]
        positive = e > 0.0
        if np.count_nonzero(positive) < 2:
            raise ValueError("at least two strictly positive errors are required")
        slope, _ = np.polyfit(np.log(h[positive]), np.log(e[positive]), 1)
        return float(slope)

    @property
    def finest_error(self) -> float:
        return float(self.error[-1])


@dataclass(frozen=True, slots=True)
class UncertaintyComponent:
    """One non-negative absolute uncertainty bound with explicit provenance."""

    name: str
    bound: float
    units: str
    provenance: str

    def __post_init__(self) -> None:
        if not self.name.strip() or not self.units.strip() or not self.provenance.strip():
            raise ValueError("name, units, and provenance must be non-empty")
        if not np.isfinite(self.bound) or self.bound < 0.0:
            raise ValueError("uncertainty bound must be finite and non-negative")


@dataclass(frozen=True, slots=True)
class RichardsonEstimate:
    """Fine-grid discretization uncertainty estimated from verified order."""

    uncertainty: UncertaintyComponent
    refinement_ratio: float
    verified_order: float


def richardson_uncertainty(
    *,
    coarse_value: float,
    fine_value: float,
    coarse_spacing: float,
    fine_spacing: float,
    verified_order: float,
    name: str,
    units: str,
    provenance: str,
) -> RichardsonEstimate:
    """Estimate fine-grid error from a coarse/fine pair and verified order.

    The estimate is |fine - coarse| / (r**p - 1), where r is the ratio of
    characteristic spacings and p is an independently verified convergence order.
    """
    for label, value in (
        ("coarse_value", coarse_value),
        ("fine_value", fine_value),
        ("coarse_spacing", coarse_spacing),
        ("fine_spacing", fine_spacing),
        ("verified_order", verified_order),
    ):
        if not np.isfinite(value):
            raise ValueError(f"{label} must be finite")
    if coarse_spacing <= fine_spacing or fine_spacing <= 0.0:
        raise ValueError("coarse_spacing must be greater than positive fine_spacing")
    if verified_order <= 0.0:
        raise ValueError("verified_order must be strictly positive")
    ratio = coarse_spacing / fine_spacing
    denominator = ratio**verified_order - 1.0
    bound = abs(fine_value - coarse_value) / denominator
    component = UncertaintyComponent(
        name=name,
        bound=float(bound),
        units=units,
        provenance=provenance,
    )
    return RichardsonEstimate(component, float(ratio), float(verified_order))


@dataclass(frozen=True, slots=True)
class UncertaintyBudget:
    """Conservative additive budget of absolute uncertainty bounds."""

    components: tuple[UncertaintyComponent, ...]

    def __post_init__(self) -> None:
        if not self.components:
            raise ValueError("uncertainty budget requires at least one component")
        units = {component.units for component in self.components}
        if len(units) != 1:
            raise ValueError("all uncertainty components must use identical units")

    @property
    def units(self) -> str:
        return self.components[0].units

    @property
    def total_bound(self) -> float:
        return float(sum(component.bound for component in self.components))

    @property
    def provenance(self) -> str:
        details = "; ".join(
            f"{component.name}={component.bound:g} {component.units}"
            f" [{component.provenance}]"
            for component in self.components
        )
        return f"conservative additive uncertainty budget: {details}"


def threshold_from_budget(name: str, budget: UncertaintyBudget) -> Threshold:
    """Create a provenance-bearing threshold exactly equal to a declared budget."""
    if not name.strip():
        raise ValueError("threshold name must be non-empty")
    return Threshold(
        name=name,
        value=budget.total_bound,
        provenance=f"{budget.provenance}; units={budget.units}",
    )

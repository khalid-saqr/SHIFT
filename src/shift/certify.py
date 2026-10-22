"""Internal mathematical residuals required before source04 interpretation."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from .mechanics import crocco_residual, hydrodynamic_charge_residual, vector_magnitude


@dataclass(frozen=True, slots=True)
class IdentityResiduals:
    """Maximum absolute residuals; no application-specific pass threshold is implied."""

    incompressibility_s_1: float
    charge_s_2: float
    crocco_m_s2: float


def identity_residuals(
    *,
    velocity_divergence_s_1: ArrayLike,
    hydrodynamic_charge_s_2: ArrayLike,
    q_invariant_s_2: ArrayLike,
    laplacian_kinetic_energy_s_2: ArrayLike,
    lamb_m_s2: ArrayLike,
    c_unsteady_m_s2: ArrayLike,
    c_head_m_s2: ArrayLike,
    c_dissipative_m_s2: ArrayLike,
) -> IdentityResiduals:
    """Compute source04 continuity, charge, and Crocco residual maxima."""
    divergence = np.asarray(velocity_divergence_s_1, dtype=float)
    charge = hydrodynamic_charge_residual(
        hydrodynamic_charge_s_2,
        q_invariant_s_2,
        laplacian_kinetic_energy_s_2,
    )
    crocco = crocco_residual(
        lamb_m_s2,
        c_unsteady_m_s2,
        c_head_m_s2,
        c_dissipative_m_s2,
    )
    return IdentityResiduals(
        incompressibility_s_1=float(np.max(np.abs(divergence))),
        charge_s_2=float(np.max(np.abs(charge))),
        crocco_m_s2=float(np.max(vector_magnitude(crocco))),
    )

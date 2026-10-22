# SHIFT architecture

SHIFT should remain a small scientific core with clearly bounded extensions. Architecture follows the experiment, not the other way around.

## Architectural rules

1. **Theory is independent of solver.** State definitions and metrics do not depend on a particular CFD package.
2. **Dimensional fields first.** Normalize only at explicit theory boundaries.
3. **Coordinates are explicit.** Arrays must carry enough metadata to identify component order, spatial coordinates, time support, units, and wall orientation.
4. **Operators are auditable.** Gradients, divergence, curl, Laplacians, temporal derivatives, and wall derivatives expose their scheme and support.
5. **No silent reconstruction.** Pressure, missing velocity components, wall normals, or correspondences are never inferred without provenance.
6. **Qualification is separate from outcome evaluation.** Numerical certificates cannot depend on whether closure succeeds.
7. **Small public API.** Internal implementation can evolve; scientific result objects and CLI behavior remain deliberate.

## Planned modules

Modules are added only when implementation starts; empty placeholder packages are avoided.

- `shift.contracts` — immutable scales, weights, thresholds, and validation.
- `shift.fields` — typed dimensional velocity/pressure/wall field containers.
- `shift.geometry` — wall meshes, normals, correspondence, and $V_\delta$.
- `shift.operators` — spatial/temporal differential operators and convergence metadata.
- `shift.mechanics` — traction/WSS, $Q$, vorticity, Lamb vector, $q_H$, Crocco terms, finite-volume descriptors.
- `shift.certify` — continuity and identity residuals plus uncertainty/convergence qualification.
- `shift.metrics` — $d_\tau,d_{\mathcal M}$ and nested source04 distances.
- `shift.models` — parametric exact and numerical flow families used by the experiment.
- `shift.search` — WSS-only candidate matching and frozen-pair manifests.
- `shift.io` — explicit adapters for HDF5/VTK/NPZ and later experimental/CFD sources.
- `shift.report` — machine-readable manifests, source tables, and publication graphics.
- `shift.cli` — stable user-facing commands.

## Data-source adapters

An adapter declares an observability receipt for every requested primitive:

- `OBSERVED` — directly measured/stored;
- `DERIVED` — obtained by a declared deterministic operator from observed fields;
- `RECONSTRUCTED` — obtained from an inverse/model problem;
- `PROJECTED` — resolved or filtered representation of a higher-dimensional quantity;
- `NOT_AVAILABLE` — unsupported by the source.

The core must be able to refuse an analysis when the requested source04 coordinate is not supportable.

## Solver strategy

SHIFT should not implement a general-purpose CFD solver unless the experiment demonstrates that this is necessary. Prefer a narrow, verified state-generation layer that can use exact solutions or a mature external numerical backend while keeping all source04 diagnostics and certification inside SHIFT.

## Performance strategy

Correctness and auditability come first, but avoid obvious waste:

- vectorize pairwise metrics;
- cache immutable state products by content hash;
- separate cheap WSS screening from expensive full-state evaluation;
- operate locally on large files rather than through cloud mounts;
- benchmark the expensive path before population-scale runs.

Optimization must preserve identical scientific outputs within declared numerical tolerance.

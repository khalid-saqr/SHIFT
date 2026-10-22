# SHIFT

**Shear-equivalent Hemodynamic Invariance Falsification Test**

SHIFT is the reference computational framework for testing whether equivalence of the full space–time wall-shear-stress (WSS) vector field closes a declared near-wall hemodynamic state.

The project implements the mechanical-closure theory developed for arterial mechanobiology and is being built as a public, reproducible companion to the associated paper.

## Scientific question

SHIFT tests the operational closure hypothesis

$$
d_\tau(A,B) \le \varepsilon_\tau
\;\Longrightarrow\;
d_{\mathcal M}(A,B) \le \varepsilon_{\mathcal M}.
$$

A reproducible pair satisfying

$$
d_\tau(A,B) \le \varepsilon_\tau,
\qquad
d_{\mathcal M}(A,B) > \varepsilon_{\mathcal M}
$$

is a **WSS-equivalent but diagnostically distinct near-wall state** for the declared state representation.

The primary dimensionless state is

$$
\mathbf m^*=
[Q^*,q_H^*,\mathbf C_U^*,\mathbf C_H^*,\mathbf C_D^*].
$$

SHIFT treats the diagnostic hierarchy as a closure-testing sequence, not as a biological causal hierarchy.

## Design principles

SHIFT is intentionally strict about scientific provenance:

- the complete WSS vector field is the primary shear observable; scalar summaries do not replace it in the closure test;
- internal mathematical identities are verified before derived state classifications are interpreted;
- tolerances are derived from independently quantified numerical or experimental uncertainty, not imported as arbitrary software defaults;
- state weights, characteristic scales, control-volume thickness, mappings, tolerances, and numerical operators are declared before closure outcomes are inspected;
- WSS-equivalent candidate pairs are frozen from WSS information before state-distance outcomes are evaluated;
- resolution and operator-sensitivity studies are first-class outputs;
- undefined diagnostics remain undefined rather than being zero-filled;
- every public result is traceable to a configuration, software version, and machine-readable manifest.

The locked scientific contract is documented in [`docs/SCIENTIFIC_CONTRACT.md`](docs/SCIENTIFIC_CONTRACT.md).

## Current status

**Phase 1 — exact theory kernel** and **Phase 2 — numerical operators and uncertainty** are implemented. The **Phase-3 parametric closure experiment is scientifically locked** in [`docs/PHASE3_PROTOCOL.md`](docs/PHASE3_PROTOCOL.md), and its exact state generator is now implemented.

The locked Phase-3 family is an exact multiharmonic biaxial pulsatile parallel-plate Navier–Stokes family with 504 predeclared states and 1512 predeclared same-WSS candidate pairs. Its parameter grid, refinement levels, uncertainty derivation, WSS-only pair-freezing rule, primary closure decision, nested outcomes, and stopping rule remain fixed.

The package currently provides:

- solver-independent source04 mechanics and the complete primary/nested distance hierarchy;
- the locked exact parallel-plate Poiseuille counterexample;
- auditable first- and second-derivative operators on explicit nonuniform Cartesian coordinates;
- gradient, divergence, curl, scalar/vector Laplacian, temporal derivative, and wall-normal derivative operators;
- numerical reconstruction of $Q$, $q_H$, Lamb fields, $C_U$, $C_H$, and $C_D$ from sampled velocity and pressure;
- independent hydrodynamic-charge and Crocco residual routes;
- exact Taylor-Green Navier-Stokes verification;
- convergence-series, Richardson uncertainty, conservative uncertainty budgets, and provenance-bearing qualification thresholds;
- qualification certificates that compare residuals with already-derived uncertainty bounds without changing either quantity;
- exact multiharmonic biaxial pulsatile channel fields with analytical velocity, pressure-gradient, WSS, source04 state, Crocco, and finite-volume Lamb references;
- deterministic enumeration of the locked 504-state Phase-3 family and refusal of off-protocol parameter/refinement values.

The Phase-3 uncertainty audit, tolerance freeze, WSS-pair manifest, and closure outcomes are not implemented yet. Curved endovascular geometry, CFD solving, and experimental-data adapters remain later phases.

Run the locked analytical benchmark with:

```bash
shift benchmark poiseuille
```

For the representative pair $H_B=2H_A$, the regression target includes $d_\tau=0$, $d_Q=0$, $d_q=0.1538018223$, and $d_{\mathcal M}=0.0742470472$, with $d_{LU}$ remaining undefined for the steady solution. The complete derivation and reference values are documented in [`docs/EXACT_POISEUILLE_BENCHMARK.md`](docs/EXACT_POISEUILLE_BENCHMARK.md).

The numerical verification and uncertainty contract is documented in [`docs/NUMERICAL_VERIFICATION.md`](docs/NUMERICAL_VERIFICATION.md).

A locked Phase-3 exact state can be generated directly from the protocol table:

```python
from shift.phase3 import generate_locked_phase3_state, locked_phase3_cases

case = locked_phase3_cases()[0]
state = generate_locked_phase3_state(
    case,
    wall_normal_samples=17,
    temporal_samples=32,
)
```

The default generator support is the locked production representation $N_y=257$, $N_t=512$; smaller values above are also protocol-approved refinement levels.

## Intended capabilities

The codebase grows only where required by the experiment. The staged target includes:

1. generation or ingestion of incompressible hemodynamic states;
2. vascular wall geometry, wall correspondence, and near-wall control-volume construction;
3. WSS evaluation from the full Cauchy traction;
4. verified three-dimensional differential operators with convergence auditing;
5. computation of $Q$, $q_H$, Lamb-vector fields, generalized Crocco terms, nonlinear–unsteady diagnostics, and finite-volume Lamb descriptors;
6. independent numerical certification of incompressibility, hydrodynamic-charge identity, and Crocco closure;
7. uncertainty-derived tolerances and reproducible WSS-equivalence search;
8. the complete source04 distance hierarchy;
9. machine-readable experiment manifests and publication figures/tables;
10. later adapters for experimentally measured or externally simulated vascular fields when their observability supports the requested quantity.

## Install

```bash
python -m pip install -e .
shift info
shift doctor
shift benchmark poiseuille
```

For development:

```bash
python -m pip install -e '.[dev]'
pytest
ruff check .
```

Python 3.11 or newer is required. NumPy is the only runtime dependency.

## Repository map

```text
src/shift/                  installable Python package
  contracts.py              immutable scales, weights, thresholds
  state.py                  dimensional/dimensionless source04 primary states
  mechanics.py              Q, vorticity, Lamb, WSS, Crocco, finite-volume tools
  operators.py              auditable spatial, temporal, wall-normal derivatives
  numerics.py               sampled-field source04 reconstruction
  uncertainty.py            convergence, Richardson estimates, uncertainty budgets
  certify.py                raw residuals and uncertainty-derived qualification
  metrics.py                d_tau, d_M, and nested source04 distances
  benchmarks.py             locked analytical regression entry points
  phase3.py                 locked Phase-3 constants, 504 cases, exact wrapper
  models/poiseuille.py      exact parallel-plate Poiseuille family
  models/taylor_green.py    exact Navier-Stokes numerical-verification family
  models/biaxial_channel.py exact multiharmonic biaxial pulsatile channel
  cli.py                    public command-line entry point
tests/                      contracts, mechanics, numerics, uncertainty, regressions
docs/SCIENTIFIC_CONTRACT.md theory-to-code invariants
docs/EXPERIMENT_PROTOCOL.md global predeclared experiment ordering
docs/EXACT_POISEUILLE_BENCHMARK.md exact Phase-1 benchmark derivation
docs/NUMERICAL_VERIFICATION.md Phase-2 operator and uncertainty contract
docs/PHASE3_PROTOCOL.md      locked Phase-3 parametric closure experiment
docs/ARCHITECTURE.md        bounded module architecture
docs/ROADMAP.md             staged implementation plan
examples/                   publication-oriented examples as they are added
```

## Reproducibility policy

A SHIFT result is intended to be publishable only when its run manifest identifies the code revision, configuration, numerical resolution, derivative scheme, characteristic scales, state weights, uncertainty procedure, thresholds, geometry mapping, and resulting qualification status. The implementation must preserve enough intermediate information to audit every closure classification.

## Citation and license

Citation metadata are provided in `CITATION.cff` and will be updated when the companion paper metadata are final. A release license must be selected before the first archival software release.

## Development rule

Scientific definitions in `docs/SCIENTIFIC_CONTRACT.md` are normative. Changes to those definitions require an explicit theory-version change and must never be introduced merely because a dataset, solver, or result is inconvenient.

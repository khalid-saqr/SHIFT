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

**Phase 1 — exact theory kernel** is implemented. The current package provides solver-independent source04 mechanics, dimensionless primary-state containers, the primary and nested closure distances, internal identity residuals, and an exact parallel-plate Poiseuille regression family. Numerical derivative operators, uncertainty-derived application thresholds, CFD backends, and vascular-data adapters belong to later phases and are not claimed by the current release.

Run the locked analytical benchmark with:

```bash
shift benchmark poiseuille
```

For the representative pair $H_B=2H_A$, the regression target includes $d_\tau=0$, $d_Q=0$, $d_q=0.1538018223$, and $d_{\mathcal M}=0.0742470472$, with $d_{LU}$ remaining undefined for the steady solution. The complete derivation and reference values are documented in [`docs/EXACT_POISEUILLE_BENCHMARK.md`](docs/EXACT_POISEUILLE_BENCHMARK.md).

## Intended capabilities

The codebase will grow only where required by the experiment. Planned capabilities include:

1. generation or ingestion of incompressible hemodynamic states;
2. vascular wall geometry, wall correspondence, and near-wall control-volume construction;
3. WSS evaluation from the full Cauchy traction;
4. three-dimensional differential operators with convergence auditing;
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

Python 3.11 or newer is required. NumPy is the only runtime dependency in the exact theory kernel.

## Repository map

```text
src/shift/                  installable Python package
  contracts.py              immutable scales, weights, thresholds
  state.py                  dimensional/dimensionless source04 primary states
  mechanics.py              Q, vorticity, Lamb, WSS, Crocco, finite-volume tools
  certify.py                continuity, charge, and Crocco residuals
  metrics.py                d_tau, d_M, and nested source04 distances
  benchmarks.py             locked analytical regression entry points
  models/poiseuille.py      exact parallel-plate Poiseuille family
  cli.py                    public command-line entry point
tests/                      contract, mechanics, metric, CLI, and regression tests
docs/SCIENTIFIC_CONTRACT.md theory-to-code invariants
docs/EXPERIMENT_PROTOCOL.md predeclared experiment ordering
docs/EXACT_POISEUILLE_BENCHMARK.md exact Phase-1 benchmark derivation
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

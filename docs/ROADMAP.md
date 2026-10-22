# SHIFT implementation roadmap

The roadmap is intentionally staged. A later phase starts only after the earlier scientific interfaces and tests are stable.

## Phase 0 — foundation

- [x] repository scientific charter
- [x] source04 contract
- [x] experiment ordering and anti-outcome-leakage rules
- [x] package/CLI skeleton
- [x] CI and contract tests
- [ ] select archival software license
- [ ] finalize companion-paper citation metadata

## Phase 1 — exact theory kernel

- [ ] dimensional field containers
- [ ] analytical differential operators for exact benchmark fields
- [ ] WSS from full traction
- [ ] \(Q\), vorticity, Lamb vector, \(q_H\)
- [ ] generalized Crocco terms
- [ ] finite-volume Lamb descriptor
- [ ] primary and nested distances
- [ ] reproduce the source04 parallel-plate counterexample exactly

Exit criterion: all closed-form benchmark values reproduce independently derived reference values and identity residuals are at floating-point level where analytically expected.

## Phase 2 — numerical operator and uncertainty kernel

- [ ] structured-grid spatial derivatives
- [ ] temporal derivatives
- [ ] wall-normal derivatives
- [ ] convergence studies
- [ ] uncertainty-derived thresholds
- [ ] qualification certificates

Exit criterion: known manufactured solutions demonstrate declared convergence behavior and tolerance derivation is outcome-independent.

## Phase 3 — parametric closure experiment

- [ ] choose bounded hemodynamic parameter family
- [ ] generate admissible states
- [ ] freeze WSS-equivalent pair search
- [ ] evaluate the complete source04 hierarchy
- [ ] resolution/robustness study
- [ ] publication result bundle

Exit criterion: the full experiment can be reproduced from one versioned configuration and produces machine-readable evidence for every reported claim.

## Phase 4 — endovascular geometry

- [ ] curved wall representation
- [ ] one-to-one normalized wall correspondence
- [ ] near-wall volume generation
- [ ] physiologic waveform utilities
- [ ] idealized stenosis/aneurysm/bifurcation benchmark geometries where scientifically justified

Exit criterion: geometry-aware results pass the same qualification and convergence requirements as simpler benchmarks.

## Phase 5 — external-data adapters

- [ ] CFD field adapter
- [ ] volumetric PIV/PTV adapter
- [ ] 4D-flow MRI adapter where observability permits
- [ ] explicit pressure-reconstruction provenance
- [ ] uncertainty propagation and refusal rules

Experimental adapters are challenges to the theory implementation, not replacements for the controlled parametric experiment.

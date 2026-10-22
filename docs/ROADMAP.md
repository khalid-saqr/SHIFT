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

- [x] dimensional exact-benchmark field and primary-state containers
- [x] analytical derivatives for the exact benchmark family
- [x] WSS from full traction
- [x] $Q$, vorticity, Lamb vector, $q_H$
- [x] generalized Crocco terms
- [x] finite-volume Lamb descriptor
- [x] primary and nested distances
- [x] reproduce the source04 parallel-plate counterexample exactly

Exit criterion: all closed-form benchmark values reproduce independently derived reference values and identity residuals are at floating-point level where analytically expected. **Satisfied by the locked Poiseuille regression suite.**

## Phase 2 — numerical operator and uncertainty kernel

- [x] structured-grid spatial derivatives
- [x] temporal derivatives
- [x] wall-normal derivatives
- [x] convergence studies
- [x] uncertainty-derived thresholds
- [x] qualification certificates
- [x] numerical source04 reconstruction from sampled velocity/pressure fields
- [x] exact Taylor-Green Navier-Stokes verification model
- [x] Poiseuille wall-gradient/WSS numerical regression

Exit criterion: known manufactured solutions demonstrate declared convergence behavior and tolerance derivation is outcome-independent. **Satisfied by the Phase-2 operator, Taylor-Green, uncertainty, and wall-WSS regression suite.**

## Phase 3 — parametric closure experiment

- [x] choose and lock bounded hemodynamic parameter family
- [x] freeze parameter grid, numerical audit, tolerance derivation, and pair-selection protocol
- [x] implement the exact biaxial pulsatile channel state generator
- [x] implement deterministic enumeration of the locked 504-state family
- [ ] run the predeclared numerical uncertainty audit
- [ ] derive and freeze Phase-3 tolerances
- [ ] compute WSS only and freeze the qualified pair manifest
- [ ] evaluate the complete source04 hierarchy for frozen pairs
- [ ] analytical-versus-numerical classification audit
- [ ] publication result bundle

The locked Phase-3 design is documented in [`PHASE3_PROTOCOL.md`](PHASE3_PROTOCOL.md). It contains 504 exact periodic states and 1512 predeclared same-WSS candidate pairs. The exact state generator implements this frozen family without changing its parameter grid, refinement levels, or outcome firewall.

Exit criterion: one versioned command reproduces the full state family, uncertainty audit, frozen-pair manifest, complete source04 hierarchy, and closure decision without post-outcome protocol changes.

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

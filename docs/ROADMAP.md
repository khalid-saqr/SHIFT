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
- [x] run the predeclared numerical uncertainty audit
- [ ] derive and freeze Phase-3 tolerances — **not permitted after the failed convergence gate**
- [ ] compute WSS only and freeze the qualified pair manifest — **not permitted after the failed convergence gate**
- [ ] evaluate the complete source04 hierarchy for frozen pairs — **not reached**
- [ ] analytical-versus-numerical classification audit — **not reached**
- [ ] publication result bundle — only the pre-pair feasibility/audit result is currently available

The locked Phase-3 design is documented in [`PHASE3_PROTOCOL.md`](PHASE3_PROTOCOL.md). The complete 504-state numerical audit is documented in [`PHASE3_UNCERTAINTY_AUDIT.md`](PHASE3_UNCERTAINTY_AUDIT.md).

The predeclared convergence requirement was $p\ge1.9$. The defining hydrodynamic-charge route produced

$$
p_{q_H}=1.53396,
$$

so the locked hard-stop rule was triggered before tolerance freezing or WSS-pair inspection. No Phase-3 closure outcome has been computed.

Phase-3 execution status:

$$
\boxed{\text{STOPPED AT THE PRE-PAIR NUMERICAL CONVERGENCE GATE}}
$$

The protocol must not be relaxed or silently revised to continue this experiment. Any future numerical redesign is a new protocol/version, not a continuation of this locked Phase-3 run.

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

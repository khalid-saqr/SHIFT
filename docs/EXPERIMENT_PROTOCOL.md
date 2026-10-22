# SHIFT experiment protocol

This protocol separates numerical qualification, WSS matching, and state-outcome evaluation so that the experiment cannot tune itself to the result.

## Stage 0 — lock the run definition

Write and hash a machine-readable configuration containing:

- theory/protocol version;
- geometry family and parameter domain;
- fluid properties and body-force model;
- temporal forcing;
- spatial and temporal discretization plan;
- derivative and interpolation operators;
- wall mapping and near-wall thickness $\delta$;
- characteristic scales;
- strictly positive state weights;
- uncertainty/convergence procedure;
- pair-selection rule.

No closure classification is allowed without this manifest.

## Stage 1 — construct admissible fluid states

Generate or ingest velocity and pressure states with explicit units, coordinates, time support, wall geometry, and provenance. The first parametric benchmark family must include exact or independently verifiable solutions before more complex endovascular geometries are admitted.

## Stage 2 — numerical qualification only

Before pair matching, quantify:

- continuity error;
- charge-identity residual;
- Crocco residual;
- spatial convergence;
- temporal convergence;
- derivative-stencil sensitivity;
- wall-gradient sensitivity;
- pressure-gradient sensitivity where applicable;
- control-volume sensitivity.

This stage may reject a state as numerically inadmissible. It must not use $d_{\mathcal M}$ to choose thresholds.

## Stage 3 — derive and freeze uncertainty thresholds

Derive tolerances from the predeclared convergence, uncertainty, or repeatability procedure. Store the threshold values, derivation inputs, and hashes. Thresholds become immutable for that protocol version.

## Stage 4 — WSS-only equivalence search

Compute the complete space-time vector WSS representation and $d_\tau$. Select candidate pairs using only WSS information and the frozen $\varepsilon_\tau$. Freeze the candidate-pair list and its hash before state-distance outcomes are inspected.

## Stage 5 — state-space evaluation

For the frozen WSS-equivalent pairs, compute:

- $d_Q$;
- $d_q$;
- $d_\ell$;
- $d_{\rm comp}$;
- $d_\Lambda$;
- $d_{LU}$;
- primary $d_{\mathcal M}$.

Undefined activity-set metrics remain undefined.

## Stage 6 — closure classification

Apply the already frozen thresholds. A pair with

$$
d_\tau\le\varepsilon_\tau,
\qquad
d_{\mathcal M}>\varepsilon_{\mathcal M}
$$

is a source04 counterexample for the tested flow class and declared state representation.

## Stage 7 — robustness

Repeat the declared resolution and sensitivity studies without changing the decision rule. Report whether each state distinction survives them.

## Stage 8 — publication artifact

Every reported run exports:

- immutable configuration;
- code revision;
- environment summary;
- state qualification table;
- threshold derivation table;
- frozen pair list;
- complete distance table;
- internal residuals;
- convergence results;
- figure/table source data;
- checksums and a final manifest.

A figure without its underlying source table is not considered a complete SHIFT result.

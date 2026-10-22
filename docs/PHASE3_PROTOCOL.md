# Phase 3 — locked parametric closure experiment

This document freezes the Phase-3 scientific design before any Phase-3 state-generator, pair-search, or closure-outcome code is written.

The experiment tests the source04 wall-shear closure hypothesis in a bounded exact Navier–Stokes family while preserving the repository ordering

$$
\text{model definition}
\rightarrow
\text{numerical uncertainty}
\rightarrow
\text{freeze tolerances}
\rightarrow
\text{WSS-only pair qualification}
\rightarrow
\text{freeze pairs}
\rightarrow
\text{mechanical-state evaluation}.
$$

No later implementation decision may change the family, parameter grid, candidate-pair population, uncertainty procedure, or closure definitions merely because a result is inconvenient.

## 1. Selected Navier–Stokes family

Phase 3 uses an exact **multiharmonic biaxial pulsatile parallel-plate flow**.

Let the stationary walls be located at $y=0$ and $y=H$. The velocity is

$$
\mathbf u(y,t)
=
\bigl(u_x(y,t),0,u_z(y,t)\bigr),
$$

and the pressure is

$$
p(x,z,t)
=
p_0-G_x(t)x-G_z(t)z.
$$

The flow is incompressible, Newtonian, periodic in time, and free of nonconservative body forces. The two tangential velocity components depend only on the wall-normal coordinate, so the convective acceleration vanishes identically while the unsteady, pressure-gradient, viscous, Lamb, charge, and finite-volume terms remain nontrivial.

For $j\in\{x,z\}$,

$$
\frac{\partial u_j}{\partial t}
=
\frac{G_j(t)}{\rho}
+
\nu\frac{\partial^2u_j}{\partial y^2},
\qquad
u=\frac{\mu}{\rho},
$$

with

$$
u_j(0,t)=u_j(H,t)=0.
$$

The exact periodic solution is used. Transient start-up states are outside Phase 3.

## 2. Why this family is selected

This family is the strongest pre-geometry experiment satisfying all of the source04 operational requirements at once.

1. **Exact Navier–Stokes consistency.** Velocity, pressure, WSS, and all source04 fields have analytical reference values.
2. **Full-vector WSS.** The lower-wall shear has two tangential components and can change direction through the cycle.
3. **Active unsteady mechanics.** $\mathbf C_U$ is nonzero for pulsatile members, extending the steady Phase-1 test.
4. **Exact WSS control.** Each Fourier component of the prescribed wall-shear waveform can be inverted to the pressure-gradient harmonic required for any channel height.
5. **Common near-wall measure.** Every comparison uses the same planar lower-wall patch and the same slab $0\le y\le\delta$, avoiding a geometry-dependent volume Jacobian or wall-correspondence ambiguity.
6. **Direct Phase-1 continuity.** The steady uniaxial member reduces to the exact parallel-plate Poiseuille family already locked in Phase 1.
7. **Clean staging.** Curvature, stenosis, aneurysm, bifurcation, and patient-specific geometry remain later challenges after the closure machinery has been tested in this exact time-dependent vector-shear family.

Multidirectional WSS is retained deliberately because disturbed arterial flows can exhibit substantial transverse shear and changing WSS direction. Phase 3 therefore tests equivalence of the complete WSS vector rather than a scalar shear magnitude.

## 3. Exact harmonic transfer

Let

$$
\omega=\frac{2\pi}{T}.
$$

For harmonic $k\ge1$ with $\Omega_k=k\omega$, define

$$
\lambda_k
=
\sqrt{\frac{i\Omega_k}{\nu}}.
$$

If $\widehat G_{j,k}$ is the complex pressure-gradient amplitude, then

$$
\widehat u_{j,k}(y)
=
\frac{\widehat G_{j,k}}{i\rho\Omega_k}
\left[
1-
\frac{\cosh\left(\lambda_k(y-H/2)\right)}
{\cosh\left(\lambda_kH/2\right)}
\right].
$$

The lower-wall shear harmonic is

$$
\widehat\tau_{j,k}
=
\frac{\nu\lambda_k}{i\Omega_k}
\tanh\left(\frac{\lambda_kH}{2}\right)
\widehat G_{j,k}.
$$

Therefore the pressure-gradient harmonic required to generate a prescribed wall-shear harmonic is fixed uniquely by

$$
\widehat G_{j,k}
=
\widehat\tau_{j,k}
\frac{i\Omega_k}
{\nu\lambda_k\tanh(\lambda_kH/2)}.
$$

For the mean component,

$$
\tau_{j,0}=\frac{H}{2}G_{j,0},
\qquad
G_{j,0}=\frac{2\tau_{j,0}}{H}.
$$

This inversion is performed before any state-distance calculation. It constructs states that share a target wall-shear waveform without using $d_{\mathcal M}$ or any nested source04 outcome.

## 4. Locked physical constants and source04 scales

The Phase-3 fluid and normalization constants are inherited unchanged from the locked source04 analytical experiment:

| Quantity | Value |
| --- | ---: |
| Density $\rho$ | $1050\ \mathrm{kg\,m^{-3}}$ |
| Dynamic viscosity $\mu$ | $0.004\ \mathrm{Pa\,s}$ |
| Near-wall thickness $\delta$ | $0.4\ \mathrm{mm}$ |
| Length scale $L_c$ | $0.4\ \mathrm{mm}$ |
| Velocity scale $U_c$ | $0.12\ \mathrm{m\,s^{-1}}$ |
| Shear scale $\tau_c$ | $1.2\ \mathrm{Pa}$ |
| State weights | $w_Q=w_q=w_U=w_H=w_D=0.2$ |
| Crocco activity fraction $\eta_a$ | $10^{-3}$ |

No body-force term is introduced.

## 5. Locked parameter grid

### 5.1 Channel height

$$
H\in\{2,3,4,6,8,12,16\}\ \mathrm{mm}.
$$

There are seven height levels. The Phase-1 representative $4$ mm and $8$ mm pair is retained inside this grid.

### 5.2 Cardiac period

$$
T\in\{0.6,0.8,1.0,1.2\}\ \mathrm{s}.
$$

Pairs are formed only within the same period stratum. No phase-rescaling is used to compare different physical periods.

### 5.3 Prescribed full-vector WSS waveforms

Let

$$
\phi=\frac{2\pi t}{T}.
$$

The lower-wall target is

$$
\frac{\boldsymbol\tau_w(t)}{\tau_c}
=
\left[1+A\sin\phi+B\sin(2\phi)\right]\mathbf e_x
+
C\cos\phi\,\mathbf e_z,
$$

with the deterministic coefficient grid

$$
A\in\{0,0.75,1.5\},
$$

$$
B\in\{0,0.3\},
$$

and

$$
C\in\{0,0.6,1.2\}.
$$

The 18 coefficient combinations are all retained.

The parameters have distinct roles:

- $A$ changes primary pulsatility and introduces axial reversal at the largest level;
- $B$ adds a second harmonic and waveform asymmetry;
- $C$ introduces transverse shear and time-varying WSS direction.

The steady regression waveform is $(A,B,C)=(0,0,0)$.

## 6. State and candidate-pair counts

The complete state family contains

$$
7\times4\times18=504
$$

states.

Within each fixed $(T,A,B,C)$ stratum, all distinct-height pairs are candidate same-WSS comparisons. Seven heights produce

$$
\binom{7}{2}=21
$$

candidate pairs per stratum.

The complete predeclared candidate population is therefore

$$
4\times18\times21=1512
$$

pairs.

No additional pair is introduced later, and no candidate pair is removed because of a mechanical-state outcome.

## 7. Near-wall comparison domain

All comparisons use the lower wall and the physical slab

$$
0\le y\le\delta.
$$

The wall patch can be assigned any finite common area $A_\Sigma$ because all Phase-3 fields are homogeneous in $x$ and $z$; the area cancels from the normalized space-time averages.

The common wall-normal coordinate is the physical distance $y$ from the lower wall. Because the compared domains are planar and use the same $\delta$, no geometry mapping or Jacobian convention is needed.

The observation interval is exactly one periodic cycle $0\le t\le T$.

## 8. Analytical source04 structure of the family

For

$$
\mathbf u=(u_x,0,u_z),
$$

the second-invariant descriptor satisfies

$$
Q=0
$$

analytically for every Phase-3 state.

Let

$$
K=\frac12\left(u_x^2+u_z^2\right).
$$

Then

$$
\boldsymbol\ell
=
\mathbf u\times\boldsymbol\omega
=
\left(0,\frac{\partial K}{\partial y},0\right),
$$

and

$$
q_H
=
-\frac{\partial^2K}{\partial y^2}.
$$

The generalized Crocco terms are

$$
\mathbf C_U
=
\left(
\frac{\partial u_x}{\partial t},
0,
\frac{\partial u_z}{\partial t}
\right),
$$

$$
\mathbf C_H
=
\left(
-\frac{G_x}{\rho},
\frac{\partial K}{\partial y},
-\frac{G_z}{\rho}
\right),
$$

and

$$
\mathbf C_D
=
-\nu
\left(
\frac{\partial^2u_x}{\partial y^2},
0,
\frac{\partial^2u_z}{\partial y^2}
\right).
$$

The exact identity

$$
\boldsymbol\ell
=
\mathbf C_U+\mathbf C_H+\mathbf C_D
$$

holds throughout the family.

The arealized finite-volume Lamb descriptor over the lower-wall slab has the exact reference

$$
\boldsymbol\Lambda_L^{(\delta)}
=
\rho\int_0^\delta\boldsymbol\ell\,\mathrm dy
=
\left(0,\rho K(\delta,t),0\right),
$$

because no slip gives $K(0,t)=0$.

The family therefore supplies exact truth for every primary state coordinate and every internal identity required by source04.

## 9. Numerical uncertainty audit

The numerical experiment uses sampled exact fields and reconstructs source04 quantities through SHIFT numerical operators. Analytical truth is retained only as the independent error reference.

### 9.1 Fixed refinement levels

Wall-normal sample counts are

$$
N_y\in\{17,33,65,129,257\}.
$$

Temporal sample counts are

$$
N_t\in\{32,64,128,256,512\}.
$$

The wall-normal grid spans $0\le y\le\delta$. The time grid spans one complete period including the endpoints.

Spatial convergence is assessed by varying $N_y$ while holding $N_t=512$. Temporal convergence is assessed by varying $N_t$ while holding $N_y=257$.

The final production representation is $N_y=257$, $N_t=512$, conditional on the convergence gate below.

### 9.2 Convergence gate

For every non-roundoff-limited quantity used to construct WSS or the primary source04 state, the maximum analytical reconstruction error across the complete 504-state family is tracked at every refinement level.

The finest three levels must show observed order

$$
p\ge1.9.
$$

Quantities whose exact value is structurally zero and whose numerical error is already roundoff-limited are assessed by absolute error rather than by fitting an unstable log-log slope.

If this convergence gate fails, Phase 3 stops before WSS pair freezing. Resolution levels, operators, and the gate are not modified after pair-state outcomes are available.

## 10. Derivation of Phase-3 tolerances

Because exact truth is available, Phase 3 uses the finest-grid analytical reconstruction error as its primary uncertainty source. Richardson estimates from the verified convergence series are retained as an independent cross-check.

For a quantity $j$, define the one-state numerical error on state $s$ as

$$
e_j(s)
=
d_j\left(j_{\mathrm{num}}(s),j_{\mathrm{exact}}(s)\right),
$$

using the same normalization and measure as the corresponding source04 comparison.

The global one-state uncertainty bound is

$$
U_j
=
\max_{s\in\mathcal S}e_j(s),
$$

where $\mathcal S$ is the complete 504-state family.

For ordinary pair distances, the predeclared pair-equivalence tolerance is the conservative two-state bound

$$
\boxed{
\varepsilon_j=2U_j.
}
$$

In particular,

$$
\varepsilon_\tau=2U_\tau,
\qquad
\varepsilon_{\mathcal M}=2U_{\mathcal M}.
$$

The same construction is used for $\varepsilon_Q$, $\varepsilon_q$, $\varepsilon_\ell$, and $\varepsilon_\Lambda$.

For normalized diagnostics whose definition depends on an activity set, numerical/exact activity-set disagreement is reported as an uncertainty-limited nested diagnostic. It does not alter the primary $d_{\mathcal M}$ closure result and is never converted to a zero distance.

No percentage tolerance, hidden safety multiplier, percentile replacement, or outcome-derived correction is permitted.

## 11. WSS-only pair qualification and freeze

The pair workflow is fixed as follows.

1. Finish the complete numerical convergence and uncertainty audit.
2. Write all uncertainty bounds and tolerances to a versioned manifest.
3. Compute WSS only for the 1512 predeclared candidate pairs.
4. Evaluate

$$
d_\tau
=
\frac{
\left\|
\boldsymbol\tau_{w,A}-\boldsymbol\tau_{w,B}
\right\|_{\Sigma,T}
}{\tau_c}.
$$

5. Retain every candidate satisfying

$$
d_\tau\le\varepsilon_\tau.
$$

6. Write the retained pair IDs, parameters, $d_\tau$, tolerance provenance, configuration hash, and code revision to an immutable frozen-pair manifest.
7. Only after the manifest hash exists may any code path evaluate $d_{\mathcal M}$ or the nested state distances for those pairs.

The pair set cannot be changed using $d_{\mathcal M}$, $d_q$, Lamb-vector quantities, Crocco composition, finite-volume descriptors, or any other mechanical-state result.

## 12. Primary closure decision

The primary hypothesis is unchanged:

$$
H_{\mathrm{shear}}:
\quad
d_\tau\le\varepsilon_\tau
\Longrightarrow
d_{\mathcal M}\le\varepsilon_{\mathcal M}.
$$

For each frozen WSS-equivalent pair, SHIFT evaluates the primary source04 distance with the locked equal weights.

A pair is a numerical counterexample when

$$
d_\tau\le\varepsilon_\tau,
\qquad
d_{\mathcal M}>\varepsilon_{\mathcal M}.
$$

Because this Phase-3 family has analytical truth, a reported counterexample must also have the same closure classification under the analytical reference evaluation. Disagreement between numerical and analytical classification is treated as numerical failure, not as evidence for closure or non-closure.

One reproducible counterexample rejects wall-shear closure for this declared family and state representation. The full distribution of margins

$$
d_{\mathcal M}-\varepsilon_{\mathcal M}
$$

is reported descriptively for all frozen pairs.

## 13. Nested source04 outcomes

After pair freezing, every frozen pair is evaluated for

$$
d_Q,
\quad
d_q,
\quad
d_\ell,
\quad
d_{\mathrm{comp}},
\quad
d_\Lambda,
\quad
d_{LU}.
$$

The family has $Q=0$ analytically for every state. This feature is intentional: Phase 3 directly tests whether matching the full WSS vector and matching the strain–rotation descriptor can coexist with reproducible differences in hydrodynamic charge and momentum-balance structure.

Unlike the steady Phase-1 family, pulsatile members have active $\mathbf C_U$. The nonlinear–unsteady balance diagnostic therefore becomes available whenever its source04 activity conditions are satisfied.

Undefined nested diagnostics remain undefined.

## 14. Controls

The following controls are fixed before implementation.

- **Phase-1 regression control:** $(A,B,C)=(0,0,0)$ includes the steady same-WSS Poiseuille family and the $H=4$ mm versus $H=8$ mm representative pair.
- **Analytical identity control:** incompressibility, hydrodynamic-charge identity, and generalized Crocco balance are evaluated against their exact values.
- **Numerical self-consistency control:** final numerical classifications are compared with exact analytical classifications without using the latter to retune thresholds.
- **Pair-selection firewall:** the frozen-pair manifest is produced before state-distance evaluation and is hash-addressed in the result bundle.

## 15. Primary outputs

The publication-grade Phase-3 result bundle must contain at least:

- the complete protocol/configuration manifest;
- the 504-state parameter table;
- spatial and temporal convergence tables;
- uncertainty budgets and all derived tolerances;
- the 1512-candidate WSS table;
- the frozen-pair manifest and hash;
- all primary and nested distances for frozen pairs;
- analytical-versus-numerical classification agreement;
- closure decision for $H_{\mathrm{shear}}$;
- nested sufficiency decisions where the diagnostics are defined;
- machine-readable tables used to generate every publication figure.

## 16. Out of scope for Phase 3

Phase 3 does not add:

- curved-vessel geometry;
- stenosis, aneurysm, or bifurcation morphology;
- a general CFD solver;
- mesh-generation infrastructure;
- experimental PIV/PTV or 4D-flow adapters;
- post-outcome tuning of waveform coefficients, heights, periods, resolution, tolerances, or pair membership.

Those extensions are later challenges to a closure framework whose exact time-dependent parametric experiment has already been frozen.

## 17. Exit criterion

Phase 3 is complete only when one versioned command can regenerate the full 504-state family, reproduce the numerical uncertainty audit, construct and freeze the WSS-qualified pair manifest, evaluate every frozen pair with the complete source04 hierarchy, and emit the final closure decision without any post-outcome protocol change.

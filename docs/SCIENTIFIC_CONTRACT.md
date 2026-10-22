# SHIFT scientific contract

This document is the normative theory-to-code contract for the source04 implementation of SHIFT. Code may add numerically necessary metadata, but it must not silently alter these definitions.

## 1. Admissible flow class

Primary fields are velocity \(\mathbf u(\mathbf x,t)\) and pressure \(p(\mathbf x,t)\), with constant density \(\rho\) and

\[
\nabla\cdot\mathbf u=0.
\]

The momentum equation is represented as

\[
\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u
=-\nabla\left(\frac p\rho\right)+\mathbf D+\mathbf b,
\qquad
\mathbf D=\frac1\rho\nabla\cdot\boldsymbol\tau.
\]

SHIFT source04 applications are restricted to absent or conservative body forces. For \(\mathbf b=-\nabla\Phi\), define

\[
\mathcal H=\frac p\rho+K+\Phi,
\qquad K=\frac12|\mathbf u|^2.
\]

## 2. Wall shear stress is the full vector boundary observable

With Cauchy stress \(\boldsymbol\sigma=-p\mathbf I+\boldsymbol\tau\), wall normal \(\mathbf n\), and wall traction \(\mathbf t_w=\boldsymbol\sigma\mathbf n\), WSS is the tangential projection

\[
\boldsymbol\tau_w=(\mathbf I-\mathbf n\mathbf n^T)\mathbf t_w.
\]

The primary closure test compares the full space-time WSS vector field. Scalar summaries, OSI-like metrics, spectra, or temporal averages may be derived diagnostics but do not replace that field in \(d_\tau\).

For corresponding wall patches,

\[
\|\mathbf g\|_{\Sigma,T}
=\left[\frac1{A_\Sigma T}\int_0^T\!\int_\Sigma\|\mathbf g\|^2\,dA\,dt\right]^{1/2},
\]

and

\[
d_\tau(A,B)=\frac{\|\boldsymbol\tau_{w,A}-\boldsymbol\tau_{w,B}\|_{\Sigma,T}}{\tau_c}.
\]

If geometries differ, a one-to-one wall correspondence must be declared before evaluating \(d_\tau\).

## 3. Three-dimensional state descriptors

Let \(\mathbf A=\nabla\mathbf u=\mathbf S+\boldsymbol\Omega\). For incompressible flow,

\[
Q=\frac12(\boldsymbol\Omega:\boldsymbol\Omega-\mathbf S:\mathbf S)
=-\frac12\operatorname{tr}(\mathbf A^2).
\]

SHIFT uses \(Q\) as a second-invariant strain-rotation descriptor, not as a complete 3-D topology classifier.

Use the Lamb convention

\[
\boldsymbol\omega=\nabla\times\mathbf u,
\qquad
\boldsymbol\ell=\mathbf u\times\boldsymbol\omega.
\]

The exact incompressible kinematic identity is

\[
q_H=\nabla\cdot(\boldsymbol\omega\times\mathbf u)
=-2Q-\nabla^2K.
\]

The generalized Crocco-Vazsonyi decomposition is

\[
\boldsymbol\ell
=\mathbf C_U+\mathbf C_H+\mathbf C_D,
\]

with

\[
\mathbf C_U=\frac{\partial\mathbf u}{\partial t},
\qquad
\mathbf C_H=\nabla\mathcal H,
\qquad
\mathbf C_D=-\mathbf D.
\]

These are algebraic components of one momentum identity, not independent biological stimuli or separate forces on the endothelium.

## 4. Finite-volume Lamb descriptor

For a declared wall-adjacent control volume \(V_\delta\) associated with wall area \(A_\Sigma\),

\[
\boldsymbol\Lambda_L^{(\delta)}
=\frac1{A_\Sigma}\int_{V_\delta}\rho\boldsymbol\ell\,dV.
\]

This is a coarse-grained control-volume descriptor. Although it has stress dimensions, it is not wall traction and is not identified as a force transmitted to an endothelial cell.

## 5. Characteristic scales and primary state

For declared \(U_c>0\), \(L_c>0\), \(\rho>0\), and \(\tau_c>0\),

\[
a_c=\frac{U_c^2}{L_c},
\qquad
q_c=\frac{U_c^2}{L_c^2},
\qquad
\Lambda_c=\rho U_c^2.
\]

The primary dimensionless operational state is

\[
\mathbf m^*=[Q^*,q_H^*,\mathbf C_U^*,\mathbf C_H^*,\mathbf C_D^*].
\]

The state weights obey

\[
w_Q+w_q+w_U+w_H+w_D=1,
\qquad w_j>0.
\]

Weights are fixed a priori and are not optimized from the closure outcome.

## 6. Primary state distance

For corresponding near-wall domains and a common observation interval,

\[
\begin{aligned}
d_{\mathcal M}^2(A,B)=\frac1{V_\delta T}\int_0^T\!\int_{V_\delta}
[&w_Q(\Delta Q^*)^2+w_q(\Delta q_H^*)^2\\
&+w_U\|\Delta\mathbf C_U^*\|^2
+w_H\|\Delta\mathbf C_H^*\|^2
+w_D\|\Delta\mathbf C_D^*\|^2]\,dV\,dt.
\end{aligned}
\]

The metric, weights, characteristic scales, state coordinates, control volume, and tolerance together define the operational comparison; they are not claimed to be a unique metric on the complete fluid state.

## 7. Closure hypothesis

WSS equivalence is

\[
d_\tau(A,B)\le\varepsilon_\tau.
\]

State equivalence is

\[
d_{\mathcal M}(A,B)\le\varepsilon_{\mathcal M}.
\]

The principal hypothesis is

\[
H_{\rm shear}:
\quad d_\tau\le\varepsilon_\tau
\Longrightarrow
d_{\mathcal M}\le\varepsilon_{\mathcal M}.
\]

A counterexample satisfies

\[
d_\tau\le\varepsilon_\tau,
\qquad
d_{\mathcal M}>\varepsilon_{\mathcal M}.
\]

It is reported as a WSS-equivalent but diagnostically distinct near-wall state pair for the declared representation and tested flow class.

## 8. Mandatory internal verification

Before state-degeneracy results are interpreted, SHIFT must verify:

1. incompressibility: \(\nabla\cdot\mathbf u\simeq0\);
2. hydrodynamic charge: \(R_q=q_H+2Q+\nabla^2K\simeq0\);
3. Crocco closure: \(\mathbf R_{CV}=\boldsymbol\ell-\mathbf C_U-\mathbf C_H-\mathbf C_D\simeq\mathbf0\).

A normalized charge residual may be reported as

\[
\epsilon_q=\frac{|q_H+2Q+\nabla^2K|}{|q_H|+2|Q|+|\nabla^2K|}.
\]

Crocco residual normalization follows the theory definition. Numerical acceptance levels are not universal physical constants; they must be justified from independent uncertainty, repeatability, or convergence analysis for the application.

## 9. Nested diagnostics

SHIFT must ultimately support the source04 hierarchy \(d_Q,d_q,d_\ell,d_{\rm comp},d_\Lambda,d_{LU}\), including their activity-set rules. If a required common active set is empty, the corresponding diagnostic is undefined and is not reported as zero.

## 10. Anti-goalpost rules

All of the following are fixed before the closure outcome is inspected:

- state coordinates and sign conventions;
- characteristic scales;
- positive state weights;
- wall and near-wall correspondence;
- \(\delta\) and control-volume construction;
- numerical operators and filtering policy;
- activity rules;
- uncertainty procedure;
- \(\varepsilon_\tau,\varepsilon_{\mathcal M}\) and nested tolerances;
- WSS-only candidate-pair selection rule.

Any later change creates a new protocol version and is reported as such. A result is reproducibly supported only when it survives the declared resolution and uncertainty tests.

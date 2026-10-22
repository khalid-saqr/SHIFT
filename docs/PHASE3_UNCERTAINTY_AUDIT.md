# Phase 3 numerical uncertainty audit

This document records the pre-pair numerical uncertainty audit required by [`PHASE3_PROTOCOL.md`](PHASE3_PROTOCOL.md). The audit was executed over the complete locked family before any WSS-pair qualification or mechanical-state outcome was implemented.

## Audit population

The audit covers all

$$
7\times4\times18=504
$$

locked Phase-3 states.

Spatial refinement uses

$$
N_y\in\{17,33,65,129,257\}
$$

at fixed $N_t=512$.

Temporal refinement uses

$$
N_t\in\{32,64,128,256,512\}
$$

at fixed $N_y=257$.

The finest-three-level convergence requirement remains

$$
p\ge1.9.
$$

For each refinement axis, an order is fitted only for quantities changed by a numerical operator on that axis. Thus WSS, $q_H$, $C_H$, $C_D$, Lamb magnitude, and $\Lambda_L$ are spatial gates; $C_U$ is the temporal gate; $Q=0$ is treated as the locked structural-zero case. The combined primary-state reconstruction error $U_{\mathcal M}$ and identity residuals are reported diagnostically but are not used to redefine the already-locked gate.

The pressure-gradient contribution to $C_H$ is the exact prescribed forcing from the locked family. Phase 3 defines no numerical $x$ or $z$ grid, so no unplanned pressure-gradient differentiation is introduced during this audit.

## Global spatial convergence

The errors below are global maxima over all 504 states at each wall-normal resolution. WSS and source04-distance errors are dimensionless. The charge residual is in $\mathrm{s^{-2}}$ and the Crocco residual in $\mathrm{m\,s^{-2}}$.

| Quantity | $N_y=17$ | $33$ | $65$ | $129$ | $257$ | Finest-three $p$ | Gate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| WSS | $7.97818\times10^{-4}$ | $2.01583\times10^{-4}$ | $5.06641\times10^{-5}$ | $1.26997\times10^{-5}$ | $3.17916\times10^{-6}$ | $1.9971$ | pass |
| $Q$ | $0$ | $0$ | $0$ | $0$ | $0$ | structural zero | pass |
| $q_H$ | $4.41792\times10^{-3}$ | $1.38271\times10^{-3}$ | $4.60915\times10^{-4}$ | $1.57969\times10^{-4}$ | $5.49645\times10^{-5}$ | $1.5340$ | **fail** |
| $C_H$ | $2.65102\times10^{-3}$ | $6.34410\times10^{-4}$ | $1.54803\times10^{-4}$ | $3.82080\times10^{-5}$ | $9.48919\times10^{-6}$ | $2.0140$ | pass |
| $C_D$ | $2.65733\times10^{-5}$ | $4.98380\times10^{-6}$ | $9.67634\times10^{-7}$ | $1.97605\times10^{-7}$ | $4.27596\times10^{-8}$ | $2.2501$ | pass |
| Lamb magnitude | $1.81088\times10^{-4}$ | $4.31304\times10^{-5}$ | $1.05138\times10^{-5}$ | $2.59475\times10^{-6}$ | $6.44469\times10^{-7}$ | $2.0140$ | pass |
| $\Lambda_L$ | $1.14635\times10^{-3}$ | $2.85909\times10^{-4}$ | $7.13924\times10^{-5}$ | $1.78375\times10^{-5}$ | $4.45804\times10^{-6}$ | $2.0006$ | pass |

The defining hydrodynamic-charge reconstruction therefore fails the locked spatial convergence requirement.

## Temporal convergence

$C_U$ is the primary state coordinate changed by the temporal derivative operator. Its global error sequence is

| $N_t$ | $32$ | $64$ | $128$ | $256$ | $512$ |
| --- | ---: | ---: | ---: | ---: | ---: |
| $U_{C_U}$ | $2.74039\times10^{-4}$ | $6.31182\times10^{-5}$ | $1.50497\times10^{-5}$ | $3.66948\times10^{-6}$ | $9.05704\times10^{-7}$ |

The finest-three observed order is

$$
p_{C_U}=2.0187,
$$

which passes the locked temporal gate.

## Production-grid diagnostic bounds

At $N_y=257$ and $N_t=512$, the global one-state numerical reconstruction bounds are:

| Quantity | Maximum one-state error | Worst locked state |
| --- | ---: | --- |
| WSS | $3.17916\times10^{-6}$ | `H3mm_T0p6s_A1p5_B0p3_C1p2` |
| $Q$ | $0$ | structural zero |
| $q_H$ | $5.49645\times10^{-5}$ | `H16mm_T0p6s_A1p5_B0p3_C0` |
| $C_U$ | $9.05704\times10^{-7}$ | `H2mm_T0p6s_A1p5_B0p3_C1p2` |
| $C_H$ | $9.48919\times10^{-6}$ | `H2mm_T0p6s_A1p5_B0p3_C1p2` |
| $C_D$ | $4.27596\times10^{-8}$ | `H6mm_T0p6s_A1p5_B0p3_C1p2` |
| $d_{\mathcal M}(\mathrm{num},\mathrm{exact})$ | $2.47500\times10^{-5}$ | `H16mm_T0p6s_A1p5_B0p3_C0` |
| Lamb magnitude | $6.44469\times10^{-7}$ | `H16mm_T0p6s_A1p5_B0p3_C0` |
| $\Lambda_L$ | $4.45804\times10^{-6}$ | `H2mm_T0p6s_A1p5_B0p3_C1p2` |

The maximum production-grid internal residuals are

$$
R_q^{\max}=225.6534632\ \mathrm{s^{-2}},
$$

$$
R_{CV}^{\max}=1.9164496\times10^{-3}\ \mathrm{m\,s^{-2}},
$$

with incompressibility residual exactly zero for the declared channel kinematics.

These production-grid values are diagnostic uncertainty measurements only. They do not become pair-equivalence tolerances because the convergence gate failed.

## Predeclared decision

The Phase-3 protocol states that failure of the convergence gate stops the experiment before WSS-pair freezing. The observed value

$$
\boxed{p_{q_H}=1.53396<1.9}
$$

therefore triggers that stop rule.

Accordingly:

$$
\boxed{\text{Phase 3 stops before tolerance freeze and before pair search.}}
$$

No $\varepsilon_\tau$, $\varepsilon_{\mathcal M}$, $\varepsilon_q$, or other pair-equivalence tolerance is frozen from this audit. No candidate pair is inspected. No closure outcome is computed.

The implementation enforces this decision: `Phase3AuditResult.pair_tolerances()` raises when the convergence gate has failed, and the command

```bash
shift phase3 audit
```

returns the machine-readable audit with `gate_passed: false` and `pair_tolerances: null`.

## Numerical interpretation

The observed $q_H$ behavior is consistent with the defining numerical route used by SHIFT Phase 2:

$$
q_H=\nabla\cdot(\boldsymbol\omega\times\mathbf u).
$$

In this family, the route composes one numerical wall-normal derivative inside vorticity with a second numerical wall-normal derivative of the resulting conventional Lamb field. The global error is boundary-sensitive and exhibits approximately $3/2$-order behavior in the locked $L^2$-type state measure. This interpretation explains the observed convergence pattern; it does not alter the frozen gate or replace the defining $q_H$ route with the independent identity route.

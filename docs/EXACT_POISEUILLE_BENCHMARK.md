# Exact parallel-plate Poiseuille regression benchmark

This benchmark is the first executable regression target for the source04
mechanical-closure theory. It uses an exact steady incompressible
parallel-plate Poiseuille family and contains no CFD solver, interpolation,
finite-difference derivative, or experimental reconstruction.

## Locked physical definition

The lower wall is shared by every state. The axial velocity is

$$
u(y)=\frac{\tau_0}{\mu}y\left(1-\frac{y}{H}\right).
$$

The imposed pressure-drop gradient satisfies

$$
-\frac{\mathrm dp}{\mathrm dx}=G=\frac{2\tau_0}{H}.
$$

The common physical near-wall domain is $0\le y\le\delta$. The benchmark
uses

- $\rho=1050\ \mathrm{kg\,m^{-3}}$;
- $\mu=0.004\ \mathrm{Pa\,s}$;
- $\tau_0=1.2\ \mathrm{Pa}$;
- $H_A=4\ \mathrm{mm}$;
- $\delta=0.4\ \mathrm{mm}$;
- $U_c=(\tau_0/\mu)\delta=0.12\ \mathrm{m\,s^{-1}}$;
- $L_c=\delta$;
- equal strictly positive primary-state weights, $w_j=0.2$.

Every member of the family has the same complete steady WSS vector on the
shared lower wall.

## Exact source04 fields

For $\mathbf u=(u(y),0,0)$,

$$
\boldsymbol\omega=(0,0,-u'),
\qquad
Q=0,
$$

and the present-paper Lamb convention gives

$$
\boldsymbol\ell=\mathbf u\times\boldsymbol\omega=(0,uu',0).
$$

The hydrodynamic charge is evaluated from the defining one-dimensional
Lamb-divergence route,

$$
q_H=-\left[(u')^2+uu''\right],
$$

while the independent kinematic identity uses

$$
q_H=-2Q-\nabla^2K,
\qquad
\nabla^2K=(u')^2+uu''.
$$

For this steady Newtonian solution,

$$
\mathbf C_U=\mathbf 0,
$$

$$
\mathbf C_H=
\left(-\frac{G}{\rho},uu',0\right),
$$

and

$$
\mathbf C_D=
\left(\frac{G}{\rho},0,0\right).
$$

Therefore $\boldsymbol\ell=\mathbf C_U+\mathbf C_H+\mathbf C_D$ exactly.

The arealized finite-volume Lamb descriptor reduces to

$$
\boldsymbol\Lambda_L^{(\delta)}
=
\int_0^\delta \rho\boldsymbol\ell\,\mathrm dy
=
\left(0,\frac{\rho}{2}u(\delta)^2,0\right).
$$

## Locked representative regression pair

For $H_B=2H_A=8\ \mathrm{mm}$, the implementation must reproduce the locked
source04 analytical results within floating-point/quadrature regression
precision:

| quantity | reference |
| --- | ---: |
| $d_\tau$ | $0$ |
| $d_Q$ | $0$ |
| $d_q$ | $0.15380182234880216$ |
| $d_\ell$ | $0.061499131829129146$ |
| $d_{\mathrm{comp}}$ | $0.030259329233964294$ |
| $d_\Lambda$ | $0.0462500013281249$ |
| $d_{\mathcal M}$ | $0.07424704720054744$ |
| $d_{LU}$ | undefined |

$d_{LU}$ is undefined because the exact benchmark is steady and
$\mathbf C_U=\mathbf 0$, so the directional alignment $A_{LU}$ has no active
set.

The suite also checks the 101-state family
$\gamma=H_B/H_A\in[0.5,3.0]$. All members retain exact WSS equivalence, while
$\gamma=1$ is the unique state-identical member of this family.

## Independent closed-form cross-check

The source04 $d_q$ value is additionally checked against the closed-form
expression

$$
d_q^2
=
36\left[
\frac{\delta^2\Delta_1^2}{3}
-\frac{\delta^3\Delta_1\Delta_2}{2}
+\frac{\delta^4\Delta_2^2}{5}
\right],
$$

where

$$
\Delta_1=\frac{1}{H_A}-\frac{1}{H_B},
\qquad
\Delta_2=\frac{1}{H_A^2}-\frac{1}{H_B^2}.
$$

This cross-check is independent of the numerical quadrature used by the
generic metric functions.

## Run it

After installation:

```bash
shift benchmark poiseuille
```

The command emits machine-readable JSON containing the complete representative
pair distances and internal residuals.

The numerical assertions used by this regression suite are software-regression
precision checks for an exact manufactured solution. They are not physiological
or experimental admissibility tolerances and do not define later SHIFT closure
thresholds.

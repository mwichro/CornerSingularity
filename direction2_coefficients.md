# Direction 2 — analytic computation of the bifurcation coefficients $c_2$, $b$

Goal: turn Assumption 6 (`ass:subcrit`) from *assumed* into *computed*, by evaluating the
quadratic ($c_2$) and cubic ($b$) Lyapunov–Schmidt coefficients for the neo-Hookean energy
with $f_{\log}$. All formulas below are CAS-verified (sympy).

## Setup

Reduced energy along the critical mode, $\mathcal E_{\mathrm{red}}(w)=\int_{B_\rho}\Psi(\bm F_0+w\,\Grad\bm\phi)\,dV+(\text{load, linear in }w)$.
Writing $E_k:=\partial_w^k\mathcal E_{\mathrm{red}}|_0$, the normal-form coefficients are
$$ c_2=\tfrac12 E_3,\qquad b=\tfrac16 E_4 $$
(matching the paper's $b=\tfrac16\langle\bm\phi,\partial^3_{\bm u}G[\bm\phi^3]\rangle$, since $G=\partial_{\bm u}\mathcal E$).

## Key structural fact: the deviatoric part drops out

$\Psi=\frac{\mu_e}{2}(\tr\bm F^T\bm F-d)+f(J)$. The deviatoric term is **quadratic in $\bm F$**, so its
third and all higher Fréchet derivatives vanish. Therefore

> **$c_2$ and $b$ are purely volumetric — they come entirely from $f(J)$.**

This is clean and worth stating in the paper on its own: the deviatoric (Hookean) response never
contributes to the nonlinearity of the bifurcation; only the barrier $f$ does.

## Explicit formulas

Let $J(w)=\det(\bm F_0+w\Grad\bm\phi)=\sum_k D_k w^k$ with $D_k$ the mixed discriminants already in
the paper: $D_0=J_0$, $D_1=\operatorname{cof}\bm F_0:\Grad\bm\phi$, $D_2=d_\phi$ (and $D_k=0$ for $k>d$).
Faà di Bruno (CAS-verified) gives, with $f^{(k)}=f^{(k)}(J_0)$,

$$ c_2=\tfrac12\!\int_{B_\rho}\!\big[f''' D_1^3+6f'' D_1 D_2+6f' D_3\big]\,dV $$
$$ b=\tfrac16\!\int_{B_\rho}\!\big[f'''' D_1^4+12f''' D_1^2 D_2+12f'' D_2^2+24f'' D_1 D_3+24f' D_4\big]\,dV $$

For $d=2$ ($D_3=D_4=0$):
$$ c_2=\tfrac12\!\int D_1\big[f''' D_1^2+6f'' D_2\big]dV,\qquad
   b=\tfrac16\!\int\big[f'''' D_1^4+12f''' D_1^2 D_2+12f'' D_2^2\big]dV. $$

The $f_{\log}$ derivatives (CAS-verified):
$f''J^2=\mu_e+2\lambda_e-2\lambda_e\log J$, $\ f'''J^3=4\lambda_e\log J-6\lambda_e-2\mu_e$, $\ f''''J^4=6\mu_e+22\lambda_e-12\lambda_e\log J$.

## Result 1 — the symmetry/vanishing of $c_2$ is now *derived*, not assumed

Every term of $c_2$ carries an **odd** power of $D_1$ ($D_1^3$ and $D_1 D_2$; $D_3$-term absent in 2D).
$D_1=\operatorname{cof}\bm F_0:\Grad\bm\phi$ is **odd** under the lateral reflection $R$ when $\bm\phi$ is
$R$-antisymmetric (and $\bm F_0$ is $R$-symmetric), while $f^{(k)}(J_0)$ and $D_2$ are $R$-even. Hence the
integrand is odd over the $R$-symmetric domain and

> **$c_2=0$ exactly for an $R$-antisymmetric critical mode** — a rigorous derivation of the §`sub:symmetry`
> claim, with the symmetry mechanism made explicit (it kills the *odd-in-$D_1$* structure). For a symmetric or
> non-symmetric mode $c_2\neq0$ generically, confirming the reviewer's transcritical concern in that case.

## Result 2 — the sign of $b$ is a geometry-dependent $O(1)$ competition (NOT model-forced)

**Regime caveat (important).** The bifurcation is analysed at $t_c$, and by the paper's own geometric
condition *buckling precedes collapse* (`rem:framework`, Lemma `lem:attain`:
$\inf_{B_\rho}\lambda_{\min}(\bm F_0(t_c))\ge c>0$), the base state at $t_c$ has $J_0\ge c^d>0$ **everywhere** in
$B_\rho$ — the barrier does **not** blow up at the bifurcation point. So the relevant evaluation is in the
$J_0=\mathcal O(1)$ regime, **not** $J_0\to0$. (An earlier version of this note sized the terms as $J_0\to0$, the
*post-$t_c$ collapse* regime; that is the wrong point for the LS coefficients.)

At $t_c$, with $J_0=\mathcal O(1)$ (and $J_0<1<J_*$, so $f''(J_0)>0$), the direct (frozen-$\bm\psi$) coefficient is
$$ b_{\mathrm{direct}} = 2\!\int_{B_\rho} f''(J_0)\,d_\phi^2\,dV \;=\; \mathcal O(1)\;>\;0, $$
an $\mathcal O(1)$ **positive** (supercritical-leaning) quantity — energetically the relief acting at fourth order,
$\tfrac12 f''(J_0)d_\phi^2 w^4>0$, stabilises. The full LS cubic carries the standard feedback correction
$$ b = b_{\mathrm{direct}} - \big\langle \partial^2_{\bm u}G[\bm\phi,\bm\phi],\,(Q\mathcal L_c Q)^{-1}\,\partial^2_{\bm u}G[\bm\phi,\bm\phi]\big\rangle, $$
the second term **negative** and also $\mathcal O(1)$. So

> **the sign of $b$ is a genuine $\mathcal O(1)$-vs-$\mathcal O(1)$ competition, decided by the corner geometry and
> loading — not by the material model alone.**

Everything in the competition except $f''$ carries geometry: $J_0(\bm x)$ (base state = geometry × loading),
$d_\phi=\det\Grad\bm\phi$ and the operator $Q\mathcal L_c Q$ (corner angle, D–N configuration, the mode). Hence the
same material $f_{\log}$ can be subcritical at a flat creasing surface (where $f''(J_0)=\mathcal O(1)$ at the Biot
threshold) yet have a *different* sign at the $90^\circ$ D–N corner. This is exactly the reviewer's Direction-1 point
that "$b$ depends heavily on the asymmetry of the boundary conditions," and it is why flat-surface creasing
subcriticality cannot be transplanted to the corner.

## Where this leaves the assumption

$b<0$ (subcriticality) is neither forced nor excluded analytically: it requires the geometry-dependent $\bm\psi$
feedback to outweigh the $\mathcal O(1)$ positive direct term. The direct term *leans the other way*
(supercritical), so there is no soft argument for $b<0$ — which justifies keeping it an honest assumption, and
strengthens the case for the geometry-robust `rem:transcritical` (the main result should not rest on a sign the
geometry decides).

## Take-aways for the paper

1. **State the deviatoric-drop-out** (clean, free, strengthens the structure of §4): $c_2,b$ are purely volumetric.
2. **Upgrade §`sub:symmetry`**: $c_2=0$ for antisymmetric modes is now *derived* from the odd-in-$D_1$ structure of
   the explicit formula, not merely asserted from abstract sector parity.
3. **The sign of $b$ is geometry-dependent, not model-only.** At $t_c$ both competing $\mathcal O(1)$ terms are
   geometry-built; the direct term leans supercritical. So $b<0$ is an honest assumption, not a soft theorem, and
   `rem:transcritical` (robust to the sign/type) is the right hedge. For Direction 3 (FEM) the quantity to measure is
   the *full* $b=\tfrac16\langle\bm\phi,\partial^3_{\bm u}G[\bm\phi^3]\rangle$ **with** the reduced $\bm\psi(w)$ solved
   (not the frozen-$\bm\psi$ integral), since the sign lives in the $\bm\psi$-feedback vs. direct-term balance —
   and that balance is set by the specific corner geometry and load.

## What I did *not* put in the paper

These findings are exploratory and one of them (Result 2) is in tension with the current $b<0$ narrative, so I left
the manuscript edits to the two safe items (bisector wording fix + robustness remark `rem:transcritical`). Whether/how
to fold Results 1–2 into §4 / §`sub:symmetry` is your call.

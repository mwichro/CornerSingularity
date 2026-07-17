# Gap review — what is not proven and what blocks the paper

Date: 2026-07-17. Scope: mathematical gaps and open problems only, ordered by severity.
Small writing/typo issues deliberately omitted. Items marked **[fixed in this pass]** were
edited directly in the `.tex`; everything else needs real work.

---

## A. Gaps in proved statements (claims presented as theorems that are not fully proved)

### A1. Theorem `thm:A` case (A): the pencil-stability step is not covered by the cited theory
The proof freezes the linear-elasticity tensor at the vertex and treats
`R = L(F) − L_lin` as a "subordinate perturbation of relative bound 0" because
`‖R‖ = o(r^{−δ})` for every `δ`, concluding the corner pencil — hence the exponents `α`
and profiles `Φ` — is unchanged. The problem: `o(r^{−δ})` for all `δ` includes
perturbations that are **O(1) or O(log 1/r) and do not vanish at the vertex**, and that is
exactly the situation here. Under the tensile kinematics `F = I + O(K r^{α−1})`, `F` does
*not* tend to `I` as `r → 0`, so `L(F(x))` does not converge to `L_lin`; along the
self-similar field it has an angular, `K`-dependent limit (and for `f_log` the volumetric
coefficient even grows like `log(1/r)`). Kondratiev/KMR stability theorems for the pencil
require the principal-part perturbation to vanish at the vertex (o(1), usually with a Dini
or Hölder rate); a non-vanishing bounded perturbation of the principal symbol *does* shift
pencil eigenvalues in general. As written, "the exponents coincide with linearised
elasticity" is proved only in the mild-loading limit `|K̂| → 0`.

**Propagates to:** `prop:blowup`(ii) (which cites this step for "common vertex pencil"),
`§sub:edgeexp` ("the nonlinear tangent inherits exactly this α"), and every downstream use
of the universal `(ν, K̂)` phase diagram.

**Possible repairs:** (i) restate the pencil claim as perturbative in `K̂` (small
stress-intensity), which is likely all the numerics need; or (ii) freeze at the genuine
angular limit `L_∞(θ)` of `L` along the self-similar field and study *that* pencil (its
eigenvalues will be `K̂`-dependent — which would itself be a finding); or (iii) prove a
log-perturbation stability theorem, which does not exist off the shelf.

**STATUS (2026-07-17): closed in principle via route (ii) — see `thmA_pencil_repair.md`.**
The genuine angular limit turns out to be `μ_e·Id` (the decoupled vector Laplacian), not a
`K̂`-dependent tensor: every nonlinear term of `eq:Lnh` carries two factors of `F^{-T}`,
and under the tensor Williams ansatz with nondegenerate angular determinant
`μ(θ) = α Φ×Φ' > 0`, both stretches diverge for generic θ and `F^{-1} → 0` like
`r^{1−α}` ("screening"). KMR then applies with the Laplacian pencil (exponents 1,3,5,… —
no singularity of the increments at the vertex); the Williams exponents `α(ν)` survive as
intermediate asymptotics on `r_K ≪ r ≪ ρ`, which is what the numerics and the (ν,K̂)
universality actually use. Implemented in the paper: new hypothesis (W) + Lemma
`lem:screen`, rewritten case (A) statement and proof, updated `rem:core`,
`prop:blowup(ii)`, `§sub:edgeexp`, conclusion. Remaining (bounded) debt, flagged in the
text: (a) absorption of the degenerate-ray layers (log-sized perturbation on angular sets
of algebraically vanishing measure — Meyers-type argument owed); (b) transmission of the
Dirichlet condition through the collapsing elastic boundary layer at the clamped face;
(c) hypothesis (W) itself should be verified numerically from the computed `Φ(θ; ν)`
(one-line check), and self-consistency of the ansatz below `r_K` remains as in A4. This
also partially resolves A7 (the one-large-stretch regime is now *derived* to be the
boundary-layer regime, not the generic one) and sharpens A2 (the finite-extensibility-cap
exemplar was not merely unexhibited but inconsistent: (NC) + `λ_max → ∞` forces `J → ∞`,
so a capped law cannot carry the assumed kinematics — now stated in `rem:core`; the pure
log barrier `f = −μ_e log J` is named in the theorem as the model exemplar).

### A2. Theorem `thm:A` case (A) has no exhibited energy satisfying its hypotheses
Case (A) needs simultaneously: barrier (`ass:reg`), inverse-channel growth (`ass:growth`),
mild growth `Θ = O(log^p J)` (`ass:mild`), and `f'' ≥ 0` globally. Of the models in the
paper, `f_quad` fails `ass:reg` and `ass:mild`; `f_log` fails `f'' ≥ 0` (the dilation
core). The "finite-extensibility cap `J_max ≤ J*`" invoked in the theorem is not an energy
in the framework — a constraint `J < J_max` is not a `C²` function on `GL⁺(d)`, and its
`Θ` is never computed. So the "genuine corner result" (A) is currently a theorem about an
empty (or at least unexhibited) class. Note `f = −μ_e log J` (i.e. `λ_e = 0`) does satisfy
everything (`f'' = μ_e/J² > 0` globally, `f''J² + f'J = 0`, `Θ ≡ μ_e` bounded), but it
removes the second Lamé constant; whether an `f` with barrier + global convexity +
logarithmic `Θ` + independent volumetric stiffness exists is a nontrivial question the
paper should answer or explicitly concede (candidate obstruction: `f'' ≥ 0` with
`f'(J)J = O(log J)` forces `f'(J) → 0⁻`, i.e. asymptotically pressure-free bulk response).
Either exhibit a concrete admissible `f` in a remark, or downgrade (A) to a conditional
statement.

**STATUS (2026-07-17): closed — the class is nonempty with both moduli free, and the
obstruction sharpened into a characterisation.** In log variables `t = log J`,
`h(t) := f'(J)J`, convexity reads `h' ≥ h` and `Θ = |h'−h| + |h|`. (i) *Forced sign
(proved, now `rem:tempered`(i)):* convexity + mild growth force `f' ≤ 0 everywhere` —
admissible materials are never volumetrically tensile; `f_log` exits the class precisely
because its volumetric stress turns positive past `J_p = e^{μ_e/2λ_e}`. (ii) *Explicit
exemplar (`eq:ftempered`):* saturate `f_log`'s linear `h` with
`h_T = −μ_e + 2λ_e T tanh(t/T)`, `0 < T < μ_e/(2λ_e)`, i.e.
`f_T(J) = −μ_e log J + 2λ_e T² log cosh(log J / T)`. Verified (analytically and
numerically): `C^∞`, `f_T(1)=0`, `f_T'(1)=−μ_e`, `f_T''(1)=μ_e+2λ_e` (identical linearised
moduli and Poisson ratio to `f_log`; `f_T = f_log + O(log⁴J)` at `J=1`), global convexity
floor `f''J² ≥ μ_e − 2λ_e T > 0`, `Θ ≤ 3μ_e + 2λ_e` bounded (`p=0`), log barrier at
`J→0⁺`. Interpolates pure log (`T→0`) ↔ `f_log` (inadmissible limit `T→∞`).
(iii) *One-sided temper:* gluing `f_log` on `J ≤ 1` to `f_T` on `J ≥ 1` is `C³` and leaves
every compressive statement of the paper pointwise unchanged while removing the dilation
core (`r_* = 0`) — the thm:A case split is an artefact of the untempered `log²` tail and
can be deleted for the concrete model at no cost. Implemented: new `rem:tempered` after
`ass:mild`, exemplar named in `thm:A`(A), table row in `§sub:growth`, finite-extensibility
caution in `rem:core` now points to the tempered penalty. Open residue: whether to
actually *adopt* the one-sided temper in the concrete realization (would make the tensile
side unconditionally case (A)) is a modelling decision left to the authors; numerics for
compressive states are unaffected by construction.

### A3. No existence/regularity theory for the compressive fundamental path — partial circularity
The whole of §4 and §5 presupposes a fundamental branch `u_0(t)` that (a) exists and is
continuous in `t`, (b) has second-variation form `Q_t` bounded (needed for `μ_1(t)` to be
finite and continuous in `lem:tc`), (c) admits the *linear* Kondratiev expansion with
pointwise gradient control of the remainder (used quantitatively in `prop:tipdisc`), and
(d) has a monotone load-to-stretch map (`lem:bcwedge`). None of this is proved for the
compressive corner. Worse, boundedness of the coefficients of `Q_t` for `0 < t < t_c`
already requires `λ_min(F_0(t))` bounded below on `B_ρ` — i.e. a (BC)-type control — while
(BC) is introduced later and only *at* `t ≤ t_c`. The text itself admits (rem:bc) that the
linear prediction collapses the tip at *any* load, so it is not even clear `Q_t` is a
bounded form for small `t` without an argument that the nonlinear state regularises the
tip. This should become an explicit standing hypothesis on the fundamental branch (with
the circularity acknowledged and localised), or a small-load existence lemma (implicit
function theorem off `t = 0` in a weighted space) should be added — the latter seems
provable with the tools already in the paper and would anchor the whole construction.

**STATUS (2026-07-17): largely closed — existence proved, circularity dissolved,
irreducible hypotheses named; A5 closed as a by-product, A4 folded into (FP).**
Findings and implementation:
1. *The IFT route proposed above does not work* — for any `t > 0` the linear Williams
   field predicts `J < 0` inside the tip disc, so the linear solution is not even in the
   energy's domain and no contraction around it is available. The correct anchor is
   variational: new `lem:fundexist` proves existence of a global minimiser for every `t`
   by the direct method, using (a) the 2D affine-cofactor identity ⇒ `det∇u` weakly
   continuous along bounded-`H¹` sequences with no singular part (the distributional
   limit is pinned to the `L¹` function `det∇u`), (b) Goffman–Serrin/affine-minorant
   lower semicontinuity for the convex penalty (no superlinearity needed), (c) barrier ⇒
   `J > 0` a.e. Also: `E_min(t) ≤ Ct²` and `E_min` continuous. Requires the penalty
   globally convex — a *second dividend of the A2 tempered class* (`f_log` itself is
   obstructed by the same large-`J` nonconvexity as its dilation core; stated honestly).
2. *What remains unprovable is named once, up front:* new standing Assumption (FP)
   (`ass:fp`), before its first use: (i) pre-collapse `λ_min` floor (compressive analogue
   of (NC); no screening under compression, cf. A1), (ii) the corner expansion with
   gradient envelope (the former `prop:tipdisc` hypothesis — A4 folded in), (iii) the
   monotone load-to-stretch map (M). New `rem:fp` states why each resists proof
   (minimisers of polyconvex energies have no regularity theory; selection continuity not
   automatic).
3. *The circularity is dissolved, not just acknowledged:* `lem:tc` rewired — it needs the
   floor only on `[0, t_*]` (the Biot load, where `λ_min ≈ λ_★ = O(1)`), and under
   `t_* < t_coll` it now *concludes* `t_c ≤ t_* < t_coll`: the (BC) ordering becomes a
   theorem given one quantitative input, a tip floor at the Biot load (checkable via
   `eq:scsafe`), rather than a hypothesis at `t_c`. Stated in `lem:tc`, `rem:fp`, and
   appended to `rem:bc`.
4. *A5 closed:* `lem:measrep` restated for any branch with `E ≤ E(t)` (minimality never
   used, only the energy bound); the minimiser branch qualifies by construction, and
   unstable continuations qualify whenever the monitored energy stays below the bound.
Remaining open (as expected, now sharply localised): the single lower bound
`t_coll > t_*` (item B2 of this review), and the (FP)(ii) expansion for the nonlinear
compressive state (shares its fate with `rem:core`'s inner problem).

### A4. `prop:tipdisc` applies linear corner asymptotics to the nonlinear state
The hypothesis "`|∇u_0^♯| ≤ C_♯ (r/ρ)^{α'−1}` \cite{KMR2001,MazyaRossmann2010}" cites
linear theory for the expansion of the *nonlinear* fundamental state in the compressive
regime — precisely where Theorem `thm:A` (tension only) gives nothing. `M(t)` bounded
along the path is likewise assumed. The proposition is fine as a conditional statement,
but it is currently listed in the ledger as an unconditional estimate. Mark its hypothesis
as part of the same standing assumption as A3, or restrict it to the regime where a
nonlinear expansion is actually available.

**STATUS (2026-07-17): folded into (FP)(ii) — see A3 status.** `prop:tipdisc` now cites
the expansion as hypothesis (FP)(ii) explicitly ("a genuine hypothesis on the nonlinear
state, not a citation"), with the no-screening-under-compression reason recorded in
`rem:fp`.

### A5. `lem:measrep` assumes the fundamental state is an energy minimiser
The Chebyshev argument needs `E(u_0(t)) ≤ E(trial)`, i.e. global (or at least
competitor-beating) minimality under displacement control. Before `t_c` the fundamental
branch is only a *local* minimiser in general (and past `t_c`, where collapse is actually
threatened, it is unstable, so the lemma does not apply there at all). Either weaken to
"any branch along which the energy stays below the explicit trial bound" (which is what is
really used, and is checkable), or restrict the statement to `t < t_c` with a stability
argument.

**STATUS (2026-07-17): closed — see A3 status, item 4.** `lem:measrep` restated for any
energy-bounded branch; the global-minimiser branch of the new `lem:fundexist` supplies
the bound by construction, and the past-`t_c` regime is covered by monitoring the energy.

### A6. Genericity invoked without a parameter
Two places rely on "generic" with nothing being perturbed: simplicity of the ground state
within the symmetry sector (`§sub:symmetry`, feeding `ass:bif`), and transversal vanishing
of `λ_min` via "Sard" (`lem:sval` — Sard needs a family; a fixed field's zero set has no
genericity by itself). For a *fixed* material, geometry and load path these are just
assumptions. Cheap repair: treat `(ν, K̂)` as the parameter family (the paper already
works on that rectangle) and phrase both as "for a.e. `(ν, K̂)`", or fold them into the
numerically-checked column of the ledger (the computed spectral gap in §6.3 is exactly a
simplicity check — say so).

### A7. The "one-large-stretch" kinematics is asserted, not derived
Scenario A's `J ~ r^{α−1}` (one diverging stretch) and Scenario B's `J ~ λ_min` are taken
from a verbal characterisation of Dirichlet–Neumann openings. The actual singular-value
structure of `∇u = K r^{α−1}(αΦ⊗e_r + Φ'⊗e_θ)` is computable from the Williams profile
`Φ(θ)` that the paper already computes to ten digits — check it once (does exactly one
singular value of the angular matrix dominate, for all θ, over the physical ν range?) and
cite the computation. Same for the self-consistency of (NC): currently a plausibility
remark, never verified on any actual tensile solution.

---

## B. Open hypotheses the architecture still rests on (correctly labelled, but load-bearing)

Ranked by how much falls if they fail.

1. **Quintic sign `d > 0`** — the *least* supported clause in the paper: the frozen part is
   provably adverse (`d_dir < 0`), the feedback `d_ψ` is unsigned by the `lem:fbsign`
   argument, and — unlike `b` — there is **no numerical computation of `d` anywhere**
   (§6.3 is silent; now stated explicitly after this pass). If `d ≤ 0` the fold structure
   `(t_L, w_L)`, the O(1) creased branch, coexistence/hysteresis and `prop:physical`'s
   post-snap story all change. Note what survives without `d`: `eq:wscale`
   (`w ~ ε^{1/3}` on the unstable continuation) and hence `thm:reg`/`prop:narrow` need
   only `b < 0`; the paper could restructure to make the regularisation results
   `d`-independent and quarantine `d > 0` to the physical-branch narrative. Either compute
   `d` in the wedge solve or do that restructuring.
2. **(BC) pointwise floor on the tip disc** — reduced (nicely) to one shrinking disc, but
   the lower bound on the barrier-regularised collapse load `t_coll` is open, and it is
   the hinge for attainment (`lem:attain`) and thus for the entire bifurcation apparatus.
3. **`eq:dphisign` (`d_φ ≥ c_d > 0` near `Σ_c`)** — needed only in the collapse limit, but
   the accessible numerics *lean against it* (deep-skin ratio `ω²/|det E| ≈ 0.6 < 1`,
   drifting down). `lem:orient`(iii)/`rem:orient` reduce it to a closed-form
   envelope-extremum question on the Stroh data — that question is answerable now with
   `surface_instability.py` outputs and should be answered before the sign is leaned on
   any further. **[floor made uniform in this pass — see C3]**
4. **(AR) rate `δJ = O(J_0/√|log J_0|)`** — only qualitatively supported (base states reach
   `J_0 ≈ 0.56`, nowhere near the limit). `thm:reg` invokes the precise rate.
5. **`|b_ψ| > b_dir` on the finite rectangle** — evidence currently from a surrogate cell
   with unsettled convergence (see C1); the inversion-free `Γ` certificate of `cor:cert` /
   `prop:certdeg` is *recommended by the paper but has never been computed*. This is the
   cheapest high-value computation available: one trilinear integral + one Rayleigh
   quotient on the existing discretisation, immune to the `cond(B)` blow-up, and it
   simultaneously tests hypothesis (b) of `prop:certdeg` (does the coupling survive
   `K̂ → 0`?).
6. **(ND) post-snap non-degeneracy** — posited on physical grounds; only a full nonlinear
   continuation (FEniCS-type path following with deflation) can check it. Honestly
   labelled; just noting that `prop:physical`'s headline "no tangent blow-up on the
   physical branch" is conditional on it.

---

## C. Numerics vs. text (the empirical weak link)

1. **The prescribed computation has not been done.** §6.3 prescribes a
   Mellin×Chebyshev solve on the wedge with the physical `r^{α−1}` grading;
   `graded_buckling.py` actually solves a flat half-space cell with an ad-hoc *linear*
   ramp. The text presented this proxy as "the model wedge". Its own refinement data
   undermine the headline: the discrete spectral gap decays monotonically
   (5.8→2.8×10⁻²) — consistent with the flat resonance re-emerging in the continuum
   limit, i.e. the surrogate may have no well-defined `b` at all — while `r` oscillates
   (−2.09/−2.45/−2.43/−2.21) and `cond(B)` grows to ~4×10⁸.
   **[fixed in this pass: §6.3 and the conclusion now state the surrogate honestly, drop
   "robust/settles/O(1) gap", flag the gap decay, and state that `d` is untested]** The
   underlying gap remains: build the real wedge solver (or at minimum run the `Γ`
   certificate, item B5, which sidesteps most of the conditioning pathology).
2. **`power`-profile margins.** The closer-to-physical `power` grading profile in
   `graded_buckling.py` gives thinner margins (r ≈ −1.53 at ν = 0.417) than the linear
   ramp reported in the tables (−1.86 to −2.19 at that ν). The reported rectangle sweep
   uses the more favourable profile; if the tables stay, the thinner `power` numbers
   should appear next to them, or the tables should be replaced by `Γ` values.
3. **`lem:relief` division by `d_φ`** — Young's inequality in case (ii) uses `1/(2 d_φ)`;
   a.e. positivity without a uniform floor makes `w_0` and the constants in `eq:Jlower`
   non-uniform on the collapse tube. **[fixed in this pass: `eq:dphisign` and
   `lem:relief` now require a uniform floor `d_φ ≥ c_d > 0`]**

---

## D. Manuscript-level holes (structural, not mathematical)

1. `splitted/introduction.tex` is **empty** (one comment line). The literature discussion
   that exists (§1.1–1.2) has no actual citations for Knowles–Sternberg, Ball, Kondratiev
   — the claims are cited by name-and-decade in prose.
2. `splitted/appendix.tex` is empty; several proofs (the Gårding dyadic-freezing argument
   in `thm:A`, the compressible Biot secular analysis behind `lem:biot`/`rem:threshold`)
   promise detail that would naturally live there.
3. `lem:biot` cites \cite{Biot1963,Hong2009} for existence, exponential decay, and strict
   negativity `q(λ_∥) < 0` of the surface mode for the *compressible* `f_log` material;
   those references treat the incompressible/gel cases. The compressible secular analysis
   is exactly what `surface_instability.py` implements — either cite a compressible
   reference or state the secular condition and its root in the appendix.
4. Stale agent file: `agents/finite-fem-paper.agent.md` describes a different project
   (PreTensorFEM GPU-FEM paper) and contradicts `AGENTS.md`'s "no references to the
   codebase" rule; it should be deleted or rewritten for this repo.

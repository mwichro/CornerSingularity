# Gap review — current state

First pass 2026-07-17; updated 2026-07-18 after the repair passes and the whole-paper
review. Margin notes in the PDF (`\claude{...}`, olive) mark in-place review comments;
this file keeps only what still needs work, plus a one-line ledger of what was closed.
Details of the closed items live in the paper itself, in `thmA_pencil_repair.md`, and in
`handoff.md` §5.0.

---

## OPEN — ranked by how much falls if it fails

### 1. Quintic sign `d > 0` (least supported claim in the paper)
Frozen part provably adverse (`d_dir < 0`), feedback `d_ψ` unsigned by the `lem:fbsign`
argument, and **no numerical computation of `d` exists anywhere** (now stated in
`rem:fbstatus`, `sub:core`, and the ledger row). If `d ≤ 0`: the fold `(t_L, w_L)`, the
O(1) creased branch, coexistence/hysteresis, and `prop:physical`'s post-snap story all
change. What survives without it: `eq:wscale` and hence `thm:reg`/`prop:narrow` need only
`b < 0`. Either compute `d` in the wedge solve or restructure to quarantine `d > 0` to
the physical-branch narrative.

### 2. Tip floor `t_coll > t_*` (the single remaining hinge of (BC))
After the A3 rewiring, (BC) is *derived* from `lem:tc` given (FP)(i) and this one bound:
the barrier-regularised collapse load must exceed the Biot load, i.e. the tip disc of
`prop:tipdisc` must not collapse by the time the free-face stretch reaches `λ_★`
(checkable via `eq:scsafe`). This powers attainment (`lem:attain`) and hence the whole
bifurcation apparatus. The paper calls it its "irreducible analytic core" — correctly.

### 3. `eq:dphisign` (`d_φ ≥ c_d > 0` near `Σ_c`) — evidence currently *leans against*
Deep-skin ratio `ω²/|det E| ≈ 0.6 < 1` and drifting down under deepening; the deepest-J₀
decile has `d_φ > 0` essentially nowhere. `lem:orient`(iii)/`rem:orient` reduce it to a
closed-form envelope-extremum question on the Stroh data of the compressible surface
eigenmode — **answerable now with `surface_instability.py` outputs, no wedge solver
needed** (margin note placed at `rem:orient`). Should be answered before the relief floor
is leaned on further; `thm:reg` now lists the floor explicitly among its hypotheses.

### 4. (AR) rate `δJ = O(J₀/√|log J₀|)`
Qualitative content confirmed (rotation defect 0.08–0.16, co-energy bounded by
cancellation), but base states reach only `J₀ ≈ 0.56` — the asymptotic rate `thm:reg`
invokes is untested and untestable with current tooling.

### 5. The true wedge solve (continuum limit of the surrogate)
`graded_buckling.py` is a flat graded half-space cell with monotone ramp, not the
prescribed Mellin×Chebyshev wedge with `r^{α−1}` grading. Its discrete gap decays under
refinement (5.8→2.8e-2), so the surrogate may lack a well-defined continuum `b`. The Γ
certificate (see ledger) certified the *discrete* sign at every resolution and
cross-validated the bordered solve, and the thinner power-profile margins are now
reported in §6 — but the continuum question and the physical-profile computation stand.
The wedge solver's job list: Γ on the rectangle with physical grading, the tip-disc (BC)
floor, and the quintic `d`.

### 6. (ND) post-snap non-degeneracy
Only a full nonlinear continuation (path following + deflation, FEniCS-class) can check
it; `prop:physical`'s "no tangent blow-up on the physical branch" is conditional on it.

### 7. Theorem A residual debts (bounded, flagged in the text)
(a) layer-absorption lemma (log-sized perturbation on angular sets of algebraically
vanishing measure — Meyers-type argument owed); (b) transmission of the Dirichlet
condition through the collapsing clamped-face elastic layer; (c) for `ν < ν_W ≈ 0.404`
the profile determinant is negative on a clamped-side sector (hypothesis (W) fails
there — `§sub:numW`), so the screened-vertex statement holds only on the free-face
sector and the self-consistent nonlinear inner state on the negative sector is open
(same class as `rem:core`'s compressive inner problem); (d) self-consistency of the
Williams ansatz below the crossover scale `r_K` in general, and the (FP)(ii) expansion
for the compressive state.

### 8. Modelling decision (authors'): adopt the one-sided tempered penalty?
`rem:tempered`(iii): gluing `f_log` (compression) to the tempered tail (dilation) is C³,
leaves every compressive statement pointwise unchanged, removes the dilation core
(`r_* = 0`, tensile side unconditionally case (A)), and makes `lem:fundexist`'s
existence theorem apply. Costs nothing mathematically; it is a modelling choice.

### 9. Manuscript-level
- `introduction.tex` still empty — an inline todo note in the PDF now specifies the
  suggested five-paragraph structure, including where `d = 2` is essential.
- §1.1–1.2 cite Knowles–Sternberg, Ball, Kondratiev by name-and-decade only; real
  citations needed.
- `appendix.tex` empty and not even `\input` in `main.tex`; the compressible Biot
  secular analysis (used by `lem:biot` Step 1 — margin note placed there: the cited
  refs cover the incompressible/gel cases only) and the Gårding dyadic-freezing detail
  belong there.
- Stale agent file `agents/finite-fem-paper.agent.md` (different project); delete.

---

## RESOLVED — one-line ledger (details in the paper / thmA_pencil_repair.md / handoff §5.0)

- **A1 pencil-stability step of thm:A** — replaced by the screening route: vertex limit
  is `μ_e·Id` (two-scale statement; hypothesis (W) + `lem:screen` + `rem:Wstatus`);
  Williams exponents are intermediate asymptotics. Residuals → OPEN 7.
- **A2 empty case-(A) class** — nonempty with both moduli free: forced-sign
  characterisation + tempered penalty `eq:ftempered` (`rem:tempered`); cap exemplar was
  inconsistent with (NC), removed.
- **A3 fundamental path** — existence proved for convex penalties (`lem:fundexist`, 2D
  direct method via affine cofactor); irreducibles named once in (FP)/`rem:fp`;
  (BC)-circularity dissolved (`lem:tc` concludes `t_c ≤ t_* < t_coll`; noted at
  `ass:bc`). Residual → OPEN 2.
- **A4 tipdisc linear asymptotics** — folded into (FP)(ii), cited as hypothesis.
- **A5 measrep minimality** — restated for any energy-bounded branch.
- **A6 genericity without a parameter** — (T) named in `lem:sval` + `rem:transgen`
  (Sard misuse removed; unconditional log-tube fallback wired in);
  simplicity grounded in von Neumann–Wigner over `(ν, K̂)` (`§sub:symmetry`).
- **A7 one-large-stretch kinematics** — derived (layer regime) and checked
  (`williams_profile_W.py`, `§sub:numW`): `tr G ≠ 0` at all zeros; (W) holds iff
  `ν > ν_W ≈ 0.4042` — partial falsification recorded honestly. Residual → OPEN 7(c).
- **B5 Γ certificate** — computed (`gamma_certificate.py`): `b < 0` certified at all 12
  rectangle points via monotone CG-Krylov test fields (k₁ = 66–184); deep CG reproduces
  the bordered ratio (2.12 vs 2.17), retiring the cond(B) worry; `prop:certdeg`'s
  hypotheses not operative on the monotone-graded cell (honest negative, in
  `rem:certnum`). Residual → OPEN 5.
- **C2 power-profile margins** — now reported in §6.3 next to the rectangle table.
- **C3 `d_φ` floor** — `lem:relief`/`eq:dphisign` carry a uniform floor `c_d`.
- **Abstract overclaims** — "preserving Kondratiev exponents" → two-scale statement;
  "we prove creasing preempts collapse" → reduced-to-tip-disc-criterion phrasing.
- **Whole-paper consistency pass (2026-07-18)** — `prop:physical` no longer claims an
  "unconditional" bound (cites (FP)(i) + `lem:tc`); `thm:reg` lists `eq:dphisign` among
  its hypotheses; `lem:biot`'s conclusion routed through the rewired `lem:tc`;
  self-consistency-of-(NC) remark cross-linked to `lem:screen`/`§sub:numW`; ledger row
  for subcriticality split into cubic (certified) vs quintic (open); (AR) "confirmed"
  softened to "qualitative content confirmed".

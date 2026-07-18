# Handoff — CornerSingularity checkpoint

Date: 2026-07-02 (supersedes the 2026-06-22 handoff; that session's d_φ story is
condensed in §4 below and fully recorded in memory).
**Standing constraint: do NOT read review3–review8.md (independent review); reading
python scripts is fine.**

## 0. Where things stand (one paragraph)

The paper's proof skeleton is sound and now substantially stronger: this session
(2026-07-02) added four theorem-level results (B1–B4 below) that convert two of the
"numerically checked" residues into part-theorems. The manuscript compiles clean
(`pdflatex` ×2 + `bibtex`, exit 0, no undefined refs/citations). The numerical weak
link is unchanged and known: `graded_buckling.py` uses a flat half-space cell with an
ad-hoc **linear** ramp instead of the wedge with the physical r^{α−1} grading, and its
own convergence data undermine the headline (gap shrinks under refinement
5.8→2.8e-2, cond(B) grows to ~4e8, ratio r drifts −2.09/−2.45/−2.43/−2.21; output in
`graded_buckling.output`, do not re-run to "verify"). The new theorems were designed
so the eventual definitive computation needs far less from the numerics.

## 1. NEW this session — four proof additions (all in the .tex, compiling)

**B1. Feedback sign is a theorem.** `lem:fbsign` + `cor:cert` + `rem:fbstatus`,
new §`sub:fbsign` in `splitted/main_text.tex` (right after prop:LS). At the first
crossing L_c ≥ 0 with simple isolated kernel ⇒ the LS feedback is
b_ψ = −2⟨R,(QL_cQ)⁻¹R⟩ ≤ 0 with R := ½E₃[φ,φ,·] on V^⊥. Verified consistent with
the code (`b = 4·E4 + 2·Rchi`, `χ = −Q⁻¹R` in graded_buckling.py). Consequences:
- ass:subcrit's cubic clause is now ONLY the magnitude |b_ψ| > b_dir (i.e. r < −1);
  r < 0 in the tables is a theorem, not evidence.
- `cor:cert`: any test field v with **E₃[φ,φ,v]² > 2 b_dir ⟨v,L_c v⟩ certifies b<0**
  — one trilinear integral + one Rayleigh quotient, NO bordered solve, immune to the
  cond(B) blow-up. This is the recommended statistic for the definitive computation.
- The quintic feedback d_ψ is NOT signed by this argument (says so in rem:fbstatus).

**B2. Small-K̂ subcriticality, conditional.** `prop:certdeg` + `rem:certnum` in
`splitted/concrete_realization.tex` (after rem:flatdegenerate). As K̂→0 the second
harmonic goes neutral (flat resonance), so the certification ratio Γ (eq:gammadef)
→ ∞ and b(ν,K̂) < 0 for all small K̂ — conditional on hypothesis (b): the
φ²→second-harmonic coupling E₃[φ,φ,v]²/(b_dir‖v‖²) stays bounded below. That
coupling non-degeneracy is the checkable residue.

**B3. Orientation identities.** `lem:orient` + `rem:orient` in §sub:reliefd
(concrete_realization.tex). (i) d_φ is a null Lagrangian (∫_ω d_φ = ∮ φ₁∂_τφ₂) —
global sign statistics are structurally meaningless; (ii) surface wave
φ=Re[U(n)e^{ikx}]: x-mean is exact, d̄_φ = −(k/2)∂_n Im(U₁Ū₂); (iii) quadrature mode
(φ₁=−a sin kx, φ₂=b cos kx): **pointwise** d_φ = −k[ab′cos² + a′b sin²], so
one-signed monotone-decaying envelopes ⇒ d_φ > 0 everywhere, no genericity. Residue
of eq:dphisign reduced to one closed-form question on the Stroh data: does the
envelope extremum (the negative skin) detach from the collapse tube as J₀→0.

**B4. (BC) quantified.** `lem:measrep` + `prop:tipdisc` in concrete_realization.tex
(after rem:bc). Ψ ≥ λ_e log²J pointwise ⇒ |{J₀≤δ}| ≤ E(t)/(λ_e log²δ) under
displacement control (measure repulsion — NOTE: does NOT exclude algebraic pointwise
tip collapse at finite energy; don't overclaim). Kondratiev floor
λ_min ≥ 1 − M(t)(r/ρ)^{α−1}, M = c_Φ μ_e|K̂| + C_♯ ⇒ any collapse is confined to the
tip disc r₀ = ρ(2M)^{1/(1−α)}; the (SC) point P is provably safe iff
M(t*) < ½(r_P/ρ)^{1−α} (eq:scsafe). (BC)'s open residue = pointwise floor on that
one shrinking disc.

**Consistency edits:** rem:bsign (main_text), §sub:core items 1–2 + closing prose,
tab:status rows for (BC) and subcriticality, numerical.tex §sub:numb framing (tables
test magnitude only; Γ recommended). All uncommitted, together with the pre-existing
working-tree changes.

## 2. Current hypothesis ledger (post-B1–B4)

- **Proved, no numerics:** material Tier 1; Θ_log vs Θ_quad dichotomy; thm:A;
  prop:blowup; lem:biot (under (SC)); lem:volcoeff (c₂=0, b_dir>0, d_dir<0);
  prop:bcount (α>½ ⇒ b=b(ν,K̂)); **lem:fbsign (b_ψ≤0)**; **lem:orient identities**;
  **lem:measrep**; **prop:tipdisc**; lem:attain (under (BC)); lem:rotreg (at t_c).
- **Conditional theorems:** prop:certdeg (b<0 for small K̂, given coupling
  non-degeneracy); (SC)-point safety given eq:scsafe.
- **Still open / numerical:** magnitude |b_ψ|>b_dir on the finite (ν,K̂) rectangle
  (via Γ>1, cor:cert); quintic d>0 (unsigned); (AR) rate eq:rotrate as J₀→0;
  eq:dphisign on the skin (= envelope/Stroh question of lem:orient); pointwise (BC)
  floor on the tip disc; (ND) post-snap non-degeneracy.

## 3. Numerics status (unchanged; do not re-run to "verify")

`graded_buckling.py` runs section-wise (`python3 graded_buckling.py C E J`, ~600s
timeout if run all); output snapshot in `graded_buckling.output`. Known defects:
(a) wrong object — flat cell + linear ramp, not the Mellin×Chebyshev wedge with
r^{α−1} grading that numerical.tex §sub:numb prescribes; the `power` profile
(section I) is the closer proxy and gives thinner margins (r≈−1.53 at ν=0.417);
(b) self-undermining convergence (gap→0 under refinement and under Lx↑, cond(B)→4e8,
r not settled). numerical.tex §sub:numb prose still overstates the tables in places
("settles near −2.3") — softening it is pending, or moot once the real wedge solver
exists. Validated and fine: `williams_clampedfree.py`, `mellin_chebyshev_pencil.py`,
`surface_instability.py`. `dominant_balance_b.py` has a syntax error (unused by
paper).

## 4. Condensed history (details in memory)

- d_φ>0 was globally false; rescoped to Σ_c with a sign-split lem:relief proof
  (case (i) needs no sign — proved); accessible evidence mildly disconfirms the sign
  at Σ_c itself (skin ratio ω²/|det E|≈0.6, drifting down). B3 now explains this
  skin structure exactly (envelope extremum) and proves positivity on the quadrature
  bulk. See memory: relief-dphi-localized-to-collapse-locus.
- b<0 story: c₂=0 proved by parity; frozen terms adverse; b_ψ≤0 now proved (B1);
  magnitude is the open residue. See memory:
  bifurcation-coefficient-sign-geometry-dependent (updated 2026-07-02).
- Geometry note: a quarter-disc body changes nothing qualitative (prop:blowup /
  rem:cleanmodel — same tangent wedge, only K̂ differs); it is merely a cleaner
  exemplar (no polyhedral vertices).

## 5. Next steps, in recommended order

0. **ADDENDUM 2026-07-17/18 (later session): items 1 and 5 are DONE, plus a large
   proof-repair pass.** See `gap_review.md` (status blocks under A1–A7, B5) and
   `thmA_pencil_repair.md` for details. Headlines: thm:A's pencil step replaced by the
   screening lemma (vertex limit is μ_e·Id, not L_lin; new ass (W) + lem:screen;
   hypothesis (W) checked numerically — holds iff ν > ν_W ≈ 0.404,
   `williams_profile_W.py`); case-(A) exemplar class populated (tempered log penalty
   `eq:ftempered`, rem:tempered); fundamental path: existence proved for convex
   penalty (lem:fundexist), standing assumption (FP) named, (BC)-circularity dissolved
   (lem:tc now concludes t_c ≤ t_* < t_coll); lem:measrep freed of minimality;
   genericity grounded in the (ν,K̂) family (hypothesis (T) + rem:transgen,
   von Neumann–Wigner for simplicity). **Γ certificate run** (`gamma_certificate.py`,
   output committed): certified b<0 at all 12 rectangle points via monotone CG-Krylov
   test fields (k₁ = 66–184; single cheap fields insufficient, Γ(R) ≈ 0.1); deep CG
   reproduces the bordered ratio (2.12 vs −r = 2.17) — cond(B) worry retired;
   prop:certdeg's hypotheses (a)/(b) NOT operative on the monotone-graded cell (no
   neutrality, erratic coupling) — that test still needs the periodic/wedge geometry.
   numerical.tex §sub:numb softening also done (this covered old item 5).

1. **Γ-based certification run** (cheapest, highest value): on the existing (or any)
   discretisation, compute Γ = E₃[φ,φ,v]²/(2 b_dir⟨v,L_c v⟩) with v = the
   second-harmonic test field, over the (ν,K̂) rectangle. Γ>1 certifies b<0
   (cor:cert) without the ill-conditioned bordered solve. Also directly tests
   prop:certdeg's hypothesis (b) (does the coupling die as K̂→0?).
   **[DONE — see addendum 0.]**
2. **Stroh/envelope evaluation for lem:orient(iii)**: closed-form envelopes a(n),
   b(n) of the compressible surface eigenmode from `surface_instability.py` data;
   locate the envelope extremum vs the collapse depth — settles the eq:dphisign
   skin question analytically or semi-analytically.
3. **The real wedge solver** (what numerical.tex actually prescribes):
   Mellin×Chebyshev on L_W with the physical r^{α−1} grading; its job is now only
   (a) Γ on the rectangle, (b) the tip-disc (BC) floor of prop:tipdisc, (c) the
   quintic sign d>0.
4. **FEniCS direct continuation** (only route to the genuine corner fundamental
   state and (ND)): monitor min λ_min along the path (BC), arc-length through the
   bifurcation, continue post-snap. Hard: needs imperfection seeding + deflation +
   graded mesh. User was considering installing FEniCS; not started.
5. Soften remaining over-claims in numerical.tex §sub:numb if the definitive
   computation is still far off. **[DONE — see addendum 0.]**

## 6. Files & memory

- Edited 2026-07-02: `splitted/main_text.tex` (§sub:fbsign, rem:bsign),
  `splitted/concrete_realization.tex` (prop:certdeg, lem:orient, lem:measrep,
  prop:tipdisc, sub:core, tab:status), `splitted/numerical.tex` (framing).
- Memory: `proof-strengthening-2026-07.md` (this session),
  `bifurcation-coefficient-sign-geometry-dependent.md` (updated),
  `relief-dphi-localized-to-collapse-locus.md` (updated),
  `review4-proof-validity-fixes.md`, `blowup-localization-and-alpha-homogeneity.md`.

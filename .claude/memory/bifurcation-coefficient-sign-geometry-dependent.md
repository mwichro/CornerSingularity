---
name: bifurcation-coefficient-sign-geometry-dependent
description: "CornerSingularity paper — sign of the cubic LS coefficient b (subcritical vs supercritical) is geometry-dependent, not model-only; deviatoric drops out; c2=0 derived for antisymmetric modes"
metadata: 
  node_type: memory
  type: project
  originSessionId: 622a60ea-580f-4e49-aadb-f94b2bc91a69
---

In the CornerSingularity manuscript (`splitted/main_text.tex` §4 imperfect-bifurcation analysis,
Assumption `ass:subcrit`), reviewer 2 (`review2.md`) challenged the subcritical-pitchfork assumption.
Key analytic findings (Direction 2, derivation in `direction2_coefficients.md`, scripts
`symbolic.py` = Faà di Bruno / f_log derivatives, `dominant_balance_b.py` = term sizing; both CAS/numerics-verified):

- **Deviatoric part drops out**: Ψ = (μ_e/2)(tr FᵀF − d) + f(J); the deviatoric term is quadratic in F, so
  3rd+ Fréchet derivatives vanish ⇒ the quadratic coeff c2 and cubic coeff b are **purely volumetric** (only f(J)).
  This is model-structural, geometry-independent.
- **c2 = 0 is derived, not assumed**: c2 = ½∫[f''' D1³ + 6 f'' D1 D2] (d=2), every term odd in
  D1 = cof F0 : Grad φ. For an R-antisymmetric mode (lateral reflection R of `sub:symmetry`, NOT the bisector
  reflection which swaps Γ_D↔Γ_N) the integrand is odd ⇒ c2=0 exactly. Generic/symmetric mode ⇒ c2≠0 (transcritical).
- **Sign of b is GEOMETRY-dependent, not model-only** (user's key point). At t_c, the paper's own
  "buckling precedes collapse" condition forces J0 ≥ c^d > 0 everywhere (no barrier blowup at the bifurcation point).
  So b = b_direct − (ψ-feedback correction) is an O(1)-vs-O(1) competition: b_direct = 2∫f''(J0) d_φ² > 0
  (supercritical-leaning), ψ-correction < 0. Both terms except f'' are geometry-built (J0, d_φ=det∇φ, QL_cQ).
  Same material f_log can crease subcritically on a flat surface yet differ at the 90° D–N corner.
  ⇒ b<0 is an honest assumption, not a soft theorem.

**Correction I made**: my first `direction2_coefficients.md` sized b_direct ~ |log J0|/J0² as J0→0 — WRONG regime
(that's the post-t_c collapse, not the bifurcation point where J0=O(1)). Fixed.

**Manuscript edits made** (all now IN the paper, compiles clean):
(1) fixed "corner bisector" wording in `ass:subcrit` (bisector swaps D↔N, is NOT a symmetry);
(2) added `rem:transcritical` — main result robust to pitchfork (Koiter 2/3, ε^{-4/3}) vs transcritical
(Koiter 1/2, ε^{-2}log) type;
(3) FULL REFRAME of `ass:subcrit` (user chose this): added Lemma `lem:volcoeff` in concrete §sub:symmetry
(deviatoric drop-out; explicit c2,b formulas; c2=0 proved via odd-in-D1, no ψ-correction since c2 is 3rd-order;
frozen b_dir=2∫f''d_φ²>0, d_dir=∫f'''d_φ³<0 in compression — ADVERSE to subcriticality). Reframed `ass:subcrit` as
a residual geometry-dependent sign hypothesis (b<0<d), added `rem:bsign` (signs decided by ψ-feedback = geometry,
verified numerically via cut-FEM; if b>0 supercritical is even more benign). Updated the hypothesis ledger table +
item 2 in concrete §sub:core. Quintic E5/E6 added to `symbolic.py`.

**UPDATE 2026-06-15 — b is now COMPUTED, not just assumed.** `graded_buckling.py` (Stage 3, consistent weak
Galerkin on the monotonically-graded model wedge; replaces the divergent Stage-2 of `creasing_buckling.py`) returns
**b<0 (subcritical) uniformly across the (ν,K̂) rectangle**. Mechanism: flat (K̂=0) cell is already indefinite at
λ★ (no discrete mode, the rem:flatdegenerate degeneracy); monotone grading lifts it → discrete localised kernel with
O(1) gap (≈4–5e-2). Sign read off the normalisation-invariant ratio r=b_feedback/b_E4≈−2.3 (<−1 ⇒ b<0); feedback
overwhelms the adverse frozen quartic. So `ass:subcrit` now HOLDS by computation. Written up in `numerical.tex`
§sub:numb. Note: this uniform b<0 is the settled numerical answer; the "geometry-dependent sign" above was the prior
theoretical expectation (the competition is real but resolves to b<0 over the physical range).

Related: [[cornersingularity-paper-overview]], [[review4-proof-validity-fixes]]

**UPDATE 2026-07-02 — feedback sign is now a THEOREM.** Added Lemma `lem:fbsign` (main_text §sub:fbsign, after
prop:LS): since G is a gradient and L_c ≥ 0 at the first crossing with simple isolated kernel, the LS feedback is
b_ψ = −2⟨R,(QL_cQ)⁻¹R⟩ ≤ 0 always, with R := ½E₃[φ,φ,·] on V^⊥ (verified consistent with graded_buckling.py's
b = 4·E4 + 2·Rchi, χ = −Q⁻¹R). So ass:subcrit's cubic clause is now ONLY the magnitude |b_ψ|>b_dir (r<−1);
r<0 in the tables is no longer evidence of anything. Companion results: Cor `cor:cert` (variational certification:
any test field v with E₃[φ,φ,v]² > 2 b_dir ⟨v,L_c v⟩ certifies b<0 — no bordered solve, immune to cond(B) blowup)
and Prop `prop:certdeg` (concrete_realization, after rem:flatdegenerate): b<0 proved for small K̂ provided the
φ²→second-harmonic coupling doesn't vanish with the gap (non-degeneracy hypothesis (b)). Quintic feedback d_ψ NOT
signed by this argument. Γ of eq:gammadef is the recommended numerical statistic going forward.

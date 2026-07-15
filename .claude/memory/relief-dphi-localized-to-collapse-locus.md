---
name: relief-dphi-localized-to-collapse-locus
description: "CornerSingularity: d_phi>0 relief hypothesis was numerically false globally; rescoped to a neighbourhood of Sigma_c + sign-split proof"
metadata: 
  node_type: memory
  type: project
  originSessionId: 00579acd-1f07-4299-b9f1-21ef01dfd412
---

The relief-orientation condition `d_φ>0` (eq:dphisign) in lem:relief was originally
asserted "a.e. on B_ρ" and "observed" on the computed mode. The localized diagnostic
in graded_buckling.py landau_b (section (J): relief_J0_10/25/50, relief_modewt,
relief_sevwt) showed this is **globally false and the wrong statistic**: det∇φ
oscillates for a wavy mode (global ~25–36% positive), is ~70–80% on the mode-active
mid-compression annulus, but **reverses (≈0% positive) in a thin surface layer at the
deepest J0**. And J0min≈0.7 there, so the cell never reaches collapse — the sign ON
Σ_c (J0→0) is not probed at all.

**Fix (2026-06-22, uncommitted):** rescoped, not dropped.
- eq:dphisign restated a.e. on a neighbourhood of Σ_c (the only region the relief
  floor is invoked), not on B_ρ.
- lem:relief proof case (i) made a **sign-split**: off Σ_c, J0 bounded below, choose
  w0≤√(J0/4|d_φ|) so a negative w²d_φ is absorbed into the J0 margin → J≥½J0 with the
  load floor C1·λc alone. So thm:reg's ε^{-4/3} no longer needs the d_φ sign anywhere
  J0 is bounded below.
- Case (ii) (near collapse) is the only place d_φ>0 is genuinely used (Young
  denominator), now flagged — like the (AR) rate eq:rotrate — a genuine hypothesis
  in the collapse limit (parallels rem:rotnum honesty).

**Verification attempt (same day, eq:dphisplit added):** exact 2D identity
`d_φ = det E + ω²` (strain/spin split; cross term vanishes since E symmetric;
verified to 1e-17). So d_φ>0 ⟺ ω²>|det E| (spin dominates strain determinant —
same rotation-dominance as AR). Diagnostics on the mode:
- frac(d_φ>0) = frac(ω²>|det E|) exactly.
- Where the mode concentrates (|∇φ|²-weighted) spin dominates ~2.4:1 → d_φ>0 (~70-80%).
- Mode PEAKS at the 8th J0-percentile (i.e. AT the collapse skin), but the deepest
  ~10% J0 decile is shear-dominated: det E<0, ω²/|det E|≈0.6<1 → d_φ<0 there.
- Mode-following depth sweep (lower l1bar past critical, track eigenvector by overlap)
  reaches only J0min≈0.5 (barrier f''~|logJ0|/J0² kills conditioning below ~0.3);
  the deep-skin ratio DRIFTS DOWN 0.65→0.61 as J0 deepens — leans AGAINST d_φ>0 at Σ_c.
**Upshot:** could NOT verify; accessible evidence mildly disconfirms d_φ>0 at the
literal collapse locus. BUT in the whole computable range (J0min≳0.5) λc is bounded
below so case-(i) sign-split covers it with NO sign hypothesis (proved). The residual
is purely the unreachable J0→0 limit. sub:reliefd wording corrected from "supported by
the accessible range" to "neither verified nor favoured there."
- prop:bcount(ii): d_φ³ noted sign-indefinite; "frozen quintic adverse" reworded
  "where relief acts"; conclusion (feedback decides subcriticality) unchanged.

This is the same honesty pattern as [[review4-proof-validity-fixes]] and the b<0
caveat in [[bifurcation-coefficient-sign-geometry-dependent]]: a numerically-decided
residue rescoped to exactly where the proof uses it. See also [[blowup-localization-and-alpha-homogeneity]].

**UPDATE 2026-07-02 — orientation now has exact identities (Lemma `lem:orient`, sub:reliefd).**
(i) d_φ is a null Lagrangian: ∫_ω d_φ = ∮ φ₁∂_τφ₂ — bulk carries NO net orientation, global fraction statistic
structurally meaningless (rescoping to Σ_c is forced, not just prudent). (ii) For a surface wave φ=Re[U(n)e^{ikx}],
the x-mean is exact: d̄_φ(n) = −(k/2) d/dn Im(U₁Ū₂); depth integral = surface phase term. (iii) Quadrature mode
φ₁=−a(n)sin kx, φ₂=b(n)cos kx ⇒ d_φ = −k[ab′cos² + a′b sin²] POINTWISE ⇒ d_φ>0 everywhere iff ab′<0 and a′b<0;
one-signed monotone-decaying envelopes ⇒ eq:dphisign holds pointwise, no genericity. Open residue reduced to: does
the envelope extremum (the negative skin) detach from or ride on the collapse tube as J0→0 — a closed-form question
on the Stroh data (two exponents + amplitude ratios of the compressible surface eigenmode), left to the wedge
computation. Explains numerics: ~70–80% weighted fractions = quadrature bulk; negative decile = envelope-extremum skin.

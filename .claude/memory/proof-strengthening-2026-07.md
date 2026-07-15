---
name: proof-strengthening-2026-07
description: "CornerSingularity: 2026-07-02 batch of new theorems (B1–B4) — feedback sign lemma, small-K̂ certified subcriticality, orientation identities, barrier measure-repulsion + tip-disc confinement of (BC)"
metadata: 
  node_type: memory
  type: project
  originSessionId: 4679a47b-b078-49b2-b50e-83942e97fd86
---

Four proof-side additions made 2026-07-02 (all compile clean, bibtex resolves; uncommitted with the rest of the
working tree). Motivation: graded_buckling.py numerics are shaky (gap shrinks under refinement, cond(B) grows,
ratio r doesn't settle — the linear-ramp grading is a flat-resonance artifact), so shift weight onto theorems.

- **B1** `lem:fbsign` + `cor:cert` + `rem:fbstatus` (main_text.tex, new §sub:fbsign after prop:LS):
  b_ψ = −2⟨R,Q⁻¹R⟩ ≤ 0 always; certification b<0 ⇐ E₃[φ,φ,v]² > 2 b_dir⟨v,L_c v⟩ for a single test field
  (no bordered solve). Details in [[bifurcation-coefficient-sign-geometry-dependent]].
- **B2** `prop:certdeg` + `rem:certnum` (concrete_realization.tex, after rem:flatdegenerate): b(ν,K̂)<0 for all
  small K̂, conditional on second-harmonic coupling non-degeneracy (hyp (b)); Γ (eq:gammadef) → ∞ as gap → 0.
- **B3** `lem:orient` + `rem:orient` (sub:reliefd): null-Lagrangian + surface-wave mean + quadrature pointwise
  sign identities for d_φ. Details in [[relief-dphi-localized-to-collapse-locus]].
- **B4** `lem:measrep` + `prop:tipdisc` (after rem:bc): Ψ ≥ λ_e log²J pointwise (deviatoric–volumetric pairing
  λᵢ²−1−log λᵢ² ≥ 0) ⇒ |{J₀≤δ}| ≤ E(t)/(λ_e log²δ) under displacement control; Kondratiev floor
  λ_min ≥ 1 − M(t)(r/ρ)^{α−1}, M = c_Φ μ_e|K̂| + C_♯ ⇒ collapse confined to tip disc r₀ = ρ(2M)^{1/(1−α)};
  (SC) point P safe iff M(t*) < ½(r_P/ρ)^{1−α} (eq:scsafe). (BC) residue = pointwise floor on that one disc.
  NOTE: log² coercivity does NOT exclude algebraic pointwise tip collapse (finite energy) — measure bound only;
  don't overclaim.

Consistency edits: rem:bsign, §sub:core items 1–2 + closing prose, tab:status rows (BC & subcriticality),
numerical.tex §sub:numb framing (tables now test magnitude only; Γ recommended as primary statistic).

Cross-checks done: lem:fbsign formulas match graded_buckling.py exactly (b = 4E4 + 2Rchi, R = ½E₃ source,
χ = −Q⁻¹R); lem:orient (iii) hand-verified twice against (ii) with U₁ = ia, U₂ = b.

Related: [[review4-proof-validity-fixes]], [[blowup-localization-and-alpha-homogeneity]]

---
name: review4-proof-validity-fixes
description: "CornerSingularity: how review4.md's 7 proof-validity findings were resolved (rescoping, not new proofs); #4 rotation condition handled via lem:rotreg + Assumption ass:rot + numerical rem:rotnum"
metadata: 
  node_type: memory
  type: project
  originSessionId: 48dae5ae-7c09-4c16-94b2-2ef495c80371
---

`review4.md` raised 7 findings that central claims were stronger than the displayed proofs. All addressed
2026-06-15 (paper compiles clean, 0 errors). Six were honesty/rescoping edits to match existing proofs:

- **#1 `thm:A`** (well-posedness, tensile corner): split into **(A)** globally-convex barrier (f''≥0, r_*=0, vertex
  INCLUDED → genuine Kondratiev corner result) and **(B)** f_log (r_*>0, vertex in non-elliptic core → only a regular
  uniformly-elliptic statement on the vertex-free annulus A_*, NO corner-pencil asymptotics across r=r_*). The flagged
  "annulus excludes the vertex" comparability step is valid only for (B); (A) gets the correct justification (log growth
  is a perturbation subordinate to the principal elasticity symbol, rel. bound 0 → Kondratiev–Maz'ya keeps the weight
  interval fixed).
- **#2 `lem:nh`**: ellipticity scoped to {f''(J)≥0} (was falsely claimed for all of f_log).
- **#3 `lem:Aellip`**: "can fail" where sharper cond |f''|J²‖F⁻¹‖²>μ_e binds, not "lost throughout core".
- **#5 `prop:physical`**: post-snap O(1) λ_min bound now from the energy barrier (f_log→∞ ⇒ J≥J_floor>0), NOT from
  lem:relief (which is only |w|≤w_0).
- **#6 `lem:bcwedge`/numerical**: buckling side t_c≤t_* unconditional; full ordering t_c<t_coll conditional on the open
  tip lower bound on t_coll.
- **#7**: compressive corner stays strongly elliptic (only tangent NORM diverges) — fixed the "ellipticity is lost"
  contradiction.

**#4 (the only one needing real work) — `lem:rot` rotation condition.** Old proof isolated f''(J₀)(δJ)² and called
the rest bounded. FALSE: the exact 2D second variation is
`⟨φ,Lφ⟩ = μ|∇φ|² + f''(J₀)(δJ)² + 2f'(J₀)det(∇φ)` (eq:secondvar, matches graded_buckling.py stiffness exactly), and
the dropped det term has prefactor 2f'~|log J₀|/J₀ — MORE singular than f''~|log J₀|. Resolution (NOT a single
theorem):
1. **`lem:rotreg` (rigorous)**: at the bifurcation J₀(t_c)>0, eq:secondvar is a regular form, neutral mode has
   ⟨φ,Lφ⟩=0, δJ & det∇φ finite — no rotation condition needed; lem:relief valid in its |w|≤w_0 range.
2. **Assumption `ass:rot` (AR)** for the near-collapse continuation: δJ=O(J₀/√|log J₀|), held by CANCELLATION (shear
   balances the singular det term), explicitly NOT by isolating f''.
3. **`rem:rotnum` (numerical, graded_buckling.py sections E/F)**: same critical mode gives co-energy≈1e-8 (confirms
   lem:rotreg), rotation defect ρ_rot=‖δJ‖/‖∇φ‖≈0.08–0.16 (stable ~0.07–0.09 under deeper compression, but J₀ only
   reaches ~0.56, NOT →0, so the asymptotic RATE is supported not proven), and T_det/T_vol2≈−6 to −10 (the dropped
   term dominates 6–10× — the reviewer was right; finiteness is by cancellation). `thm:reg` hypothesis list now also
   cites ass:rot.

Posture throughout mirrors how `ass:subcrit` is handled: rigorous where possible, else an assumption verified by the
spectral solver.

**review5.md (2nd pass, all 5 fixed 2026-06-15, compiles clean).** Reviewer accepted the thm:A split; raised follow-ups
on the fixes: (1) thm:A case (A) wrongly cited f_quad as example — f_quad fails ass:reg (no barrier) AND ass:mild
(Θ_quad algebraic); replaced with "f_log with finite-extensibility cap J_max≤J_*". (2) rem:rotnum/numerical overstated
(AR): the data verify qualitative volume-preservation + finiteness, NOT the J₀-dependent rate δJ=O(J₀/√|logJ₀|) (J₀ only
reaches ~0.56); added explicit "rate remains an assumption" caveat. (3) prop:physical's "∫f(J)<∞ ⇒ pointwise J≥J_floor"
is FALSE (integral bound = measure control only); restated post-snap bound as a non-degeneracy property (ND) of the
stable creased branch, numerically supported — compressed-branch part stays rigorous. (4) summary prose slips ("barrier-
forced theorem", "barrier-forced kinematic relief") fixed to "assumption (AR), numerically supported". (5) lem:rotreg
claimed lem:relief "applies in proven range" at t_c but lem:relief's proof used the small-J₀ step; added an explicit
finite-J₀ case (i) to lem:relief (bound holds for w small relative to J₀, NOT needing D₁ small) + near-collapse case (ii)
via AR. lem:rotreg now cites case (i).

Related: [[bifurcation-coefficient-sign-geometry-dependent]], [[blowup-localization-and-alpha-homogeneity]]

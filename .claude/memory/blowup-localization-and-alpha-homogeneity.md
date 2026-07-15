---
name: blowup-localization-and-alpha-homogeneity
description: "CornerSingularity §sub:blowup/§sub:bcount — blow-up reduction to a universal model wedge; alpha>1/2 homogeneity count means (BC) localizes (nu-only) but cubic b stays b(nu,Khat); script williams_clampedfree.py"
metadata: 
  node_type: memory
  type: project
  originSessionId: 80bc71de-57df-4594-acdd-d48a0bf0d7ed
---

CornerSingularity paper. Implemented (2026-06-15) the "localized to the corner/edge" formalization the user asked for.
Three additions to `splitted/concrete_realization.tex` (+ a tweak to `main_text.tex` rem:bsign), all compile clean (32 pp).

- **Blow-up reduction (Prop `prop:blowup`, §`sub:blowup`)**: rescale Ω_λ=λ⁻¹(Ω−x₀)→ tangent wedge W. Rescaled
  fundamental field û_λ=λ^{-α}(u₀(x₀+λξ)−…) → K r^α Φ(θ) in V^{1,2}_{β⋆}, rate O(λ^g), g=α'−α = spectral gap of the
  pencil. Operator/spectrum converge to model-wedge L_W; only the dimensionless SIF K̂=K ρ^{α-1}/μ_e survives from the
  body. 3D→2D: Fourier along smooth edge Σ decomposes L_W into planar quarter-plane fibres, leading exponent at zero
  along-edge wavenumber. Rem `rem:sifparam` = "one surviving global datum", everything else is a (ν,K̂) phase diagram.

- **α(ν) computed** in new script `williams_clampedfree.py` (Kolosov–Muskhelishvili φ=Cz^α, ψ=Ez^α; 4×4 BC determinant:
  clamped u_r=u_θ=0 at θ=0, traction-free σ_θθ=σ_rθ=0 at θ=π/2; κ=3−4ν plane strain). Smallest real root in (0,1):
  **α(ν)∈(0.60,0.87) for ν∈(0,½)**, e.g. α(0.3)≈0.711 (matches classical 90° fixed–free), α(0.45)≈0.622,
  α(0.499)≈0.595. So **α>½ throughout the physical range.**

- **Homogeneity dichotomy (Prop `prop:bcount`, §`sub:bcount`)**: in 2D the single-mode cubic integrand d_φ²~r^{4(α-1)},
  so ∫ d_φ² dV ~ ∫ r^{4α-3} dr converges at the tip iff α>½. Since α>½ here, b_dir is **tip-convergent but NOT
  tip-dominated** → O(1) residue from the matching annulus = the K̂ dependence. Therefore **b does NOT collapse to a
  bare-wedge number; it stays b(ν,K̂)** (graded-substrate creasing, Rem `rem:twoscales`: crease scale 1/k vs distance
  r_P). This is the analytic origin of the "O(1) geometry competition" of [[bifurcation-coefficient-sign-geometry-dependent]].
  Complementary α<½ (not realized here) would make b tip-localized & K̂-free.

- **(BC) localizes (Lem `lem:bcwedge`)**: t_c<t_coll is a comparison of two *scale-invariant* wedge thresholds (creasing
  λ_⋆, collapse), amplitude-free ⇒ depends on ν only, independent of body and K̂. So of the two open conditions: (BC)
  is the easy/localizing half (ν-only computation), the cubic sign is the hard/non-localizing half (ν,K̂ phase diagram).

Status table in §`sub:core` updated with both rows + a new "Localisation to model wedge: proved" row.

**Why this matters**: upgrades the "numerical part of the proof" from spot-checks on two bodies (H half-ball, C box) to
computations on a *universal* object (curve / phase diagram), and makes 3D reduce to the same 2D problem.

**Spectral solver prototype (no FEM needed).** User's idea: the cube/box is suited to spectral methods. Key correction:
the corner is *scale*-structured not *translation*-structured, so the natural transform is **Mellin** (= Fourier in
log r), not a Cartesian FFT; plain Fourier applies only along the smooth edge (the 3D->2D reduction) and possibly the
crease-along-edge wavenumber. Built and validated (2026-06-15):
- `mellin_chebyshev_pencil.py`: Kondratiev/Mellin elasticity pencil at the wedge via Chebyshev collocation in theta;
  separable u=r^lambda(a(theta),b(theta)) in plane-strain Navier -> quadratic eigenproblem in lambda. Reproduces
  alpha(nu) to ~1e-10 vs the secular determinant, and converges EXPONENTIALLY (N=12 angular modes -> machine precision;
  FEM is algebraic) -- the concrete "beats FEM" demonstration.
- `surface_instability.py`: Biot/creasing threshold of the compressed compressible neo-Hookean half-space (the buckling
  building block of Lemma biot). Stroh/characteristic-root secular method; det is purely imaginary, sign change in
  Im gives lambda_star. Validated: Lam/mu->inf gives lambda_star->0.5437 (classical incompressible Biot). Gives the
  (BC) ingredient lambda_star(nu), e.g. 0.514 at Lam/mu=1, 0.543 at Lam/mu=50.
Verdict on FEM: (BC) and alpha are essentially analytic/semi-analytic; the only genuine remaining solve is b(nu,Khat),
a 2D transverse buckling+Lyapunov-Schmidt problem best done Mellin x Chebyshev (+ edge-Fourier in 3D), NOT plain FFT
(the crease wavelength competes with the pre-stress grading in the SAME direction in 2D -- Rem twoscales).

Solver status (2026-06-15):
- `creasing_buckling.py` STAGE 1 (linear buckling, Chebyshev-in-depth half-space) VALIDATED: critical stretch matches
  the independent Stroh solver to ~1e-12. b_dir = 2 INT f''(J0) d_phi^2 computed, stable, >0 (confirms paper's adverse
  frozen sign).
- STAGE 2 (weakly-nonlinear Landau coefficient b = b_dir + b_psi via LS) IMPLEMENTED but NOT CONVERGED: feedback term
  diverges with N. Root cause diagnosed: phi taken from the STRONG-form (collocation) null space is fed into a WEAK-form
  (Galerkin) nonlinear assembly -> inconsistent spurious-mode content pollutes the 2nd-order source. FIX (next session):
  one consistent weak formulation -- phi from the Galerkin Q_1 generalised eigenproblem, m=2 harmonic in REAL (cos/sin)
  DOFs to kill the conjugation ambiguity. Do NOT trust the printed b sign until this is fixed.
Section sec:numerical / sub:numb in the paper correctly leaves b(nu,Khat) as the open computation; do not insert a b
sign until the weak-form solver converges.

RESOLVED (2026-06-15, `graded_buckling.py`): the open b(ν,K̂) solve is DONE and written into §sub:numb.
Consistent weak (Galerkin) scheme on the model wedge: Chebyshev along free face X with a **monotone** grading
λ∥(X)=λ̄(1+K̂ X/L) through the Biot threshold (turning point), mapped Chebyshev in depth. Critical mode φ = kernel of
the SAME weak operator inverted in the LS step (fixes the Stage-2 strong/weak inconsistency); 2nd harmonic in real
nodal DOFs. Results, each robust under refinement: (i) flat K̂=0 cell degenerate (indefinite at λ⋆, no discrete mode);
(ii) monotone grading → single LOCALISED neutral mode, discrete kernel with O(1) spectral gap (~0.04–0.05) → LS
non-resonant, ψ well defined; (iii) sign via the NORMALISATION-INVARIANT ratio b_fb/b_E4 (both ~‖φ‖⁴) which converges
to ≈−2.3 < −1 ⇒ feedback overwhelms the adverse frozen quartic ⇒ **b<0 SUBCRITICAL throughout the (ν,K̂) rectangle**.
CAVEAT: |b| itself is mode-amplitude-normalisation dependent (drifts with N under max|∇φ|=1); only the SIGN (the ratio)
is the converged physical output. CAVEAT 2: a PERIODIC (cosine) grading does NOT lift the degeneracy — it only splits
the scale-free Biot spectrum into near-degenerate clusters; the MONOTONE turning-point grading (the corner/Mellin
grading) is essential. Updates [[bifurcation-coefficient-sign-geometry-dependent]] (sign now pinned: b<0).

KEY FINDING (2026-06-15): the flat half-space single-mode b is NOT just mis-discretized -- it is genuinely ILL-DEFINED.
Verified with the validated linear solver: at l_star, smallest_sv ~ 1e-12 for ALL k=1,2,3,4 (scale-free Biot threshold
is wavenumber-independent). Hence the 2nd-harmonic (m=2) LS response is RESONANT (Q_{k=2} singular at l_star) and chi_2
diverges for any solver, weak or strong. Consequences: (1) the flat single-mode Landau coefficient does not exist -- this
is the precise reason flat creasing is non-classical/localized (Hong, Cao-Hutchinson handle it via localized modes);
(2) the CORNER is REQUIRED for b to be defined: the pre-stress grading breaks scale-invariance, lifts the all-k
degeneracy, selects a discrete mode (Lemma attain under BC), making b(nu,Khat) well-defined. So b must be computed on the
GRADED substrate / actual wedge (Mellin x Chebyshev, variable coeff), not the flat half-space. This confirms prop:bcount
+ rem:twoscales: b cannot be a bare-wedge (Khat-free) number. Candidate paper addition: a remark that the flat
single-mode reduction is degenerate and the corner regularizes it (strengthens why b is irreducibly geometric).

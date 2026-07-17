# Closing gap A1: the pencil-stability step of Theorem A

**Verdict up front.** The step as written cannot be repaired: the class
`‖R‖ = o(r^{−δ}) ∀δ` is *not* an admissible perturbation class for the corner
pencil, and the vertex limit of the nonlinear tangent is *not* the
linear-elasticity tensor — for any stress intensity `K ≠ 0`, not just large
ones. **But the theorem is repairable in a stronger form.** The missing
assumption is kinematic: once the Williams ansatz is imposed as a genuine
tensor ansatz (not just a norm bound), the 2D tangent formula `eq:Lnh` forces
the coefficient field to stabilise at the vertex to `μ_e 𝕀` — the *decoupled
vector Laplacian* — with an *algebraically* decaying remainder off explicit
thin angular layers. KMR theory then applies verbatim with the Laplacian
pencil, whose mixed Dirichlet–Neumann exponents on the quarter-plane are
`s = 1, 3, 5, …`: **the Newton-step increments are regular (Lipschitz) at the
vertex**, and the linear-elasticity Williams exponents `α(ν)` survive exactly
as *intermediate asymptotics* on the annulus `r_K ≪ r ≪ ρ`,
`r_K := |K|^{−1/(1−α)}`. This is both stronger and physically cleaner than the
original claim.

---

## 1. Why the original step is unfixable

The proof of `thm:A` case (A) writes `𝕃 = 𝕃_lin + ℝ` with `𝕃_lin = 𝕃(I)` and
argues that `‖ℝ‖ = o(r^{−δ})` for every `δ > 0` makes ℝ a "subordinate
perturbation" that leaves the pencil unmoved. Two independent failures:

1. **The perturbation class is wrong.** `o(r^{−δ}) ∀δ` contains every bounded
   perturbation and every `O(log 1/r)` perturbation. A *bounded, non-vanishing*
   perturbation of the principal symbol replaces the pencil wholesale (the
   pencil is determined by the limit coefficients at the vertex, and KMR
   stability requires `‖R(x)‖ → 0` as `r → 0` — continuity/stabilisation of
   coefficients, not mere subordination to algebraic weights). A log-growing
   coefficient is worse: in Mellin variables `t = log(1/r)` it is a linearly
   growing coefficient, outside the KMR class entirely.

2. **ℝ genuinely does not vanish relative to `𝕃_lin` at the vertex.** For any
   `K ≠ 0`, `F̄ = I + O(K r^{α−1})` leaves every neighbourhood of `I` as
   `r → 0`, so `𝕃(F̄) ↛ 𝕃(I)`. Smallness of `K̂` does not help pointwise: it
   only shrinks the radius `r_K = |K|^{−1/(1−α)}` below which the field is
   large. The claim "exponents coincide with linear elasticity" is therefore
   false *at the vertex* for every `K ≠ 0`; it is correct only on
   `r ≫ r_K` (intermediate asymptotics).

Also note an internal inconsistency in the case-(A) exemplar as stated: under
(NC) (`λ_min ≥ c`) and the tensile blow-up `λ_max ~ |K| r^{α−1}`, one has
`J ≥ λ_min λ_max ≥ c |K| r^{α−1} → ∞`, so a finite-extensibility **cap
`J ≤ J_max` is incompatible with the assumed corner kinematics**: a capped law
does not regularise the singular field, it forbids it (the state would have to
reorganise, invalidating the ansatz). The consistent exemplar of the case-(A)
hypothesis class is the pure logarithmic barrier `f(J) = −μ_e log J`
(`f''J² = μ_e > 0` globally, `Θ ≡ μ_e`, barrier at `J = 0⁺`, `f'(1) = −μ_e`).

## 2. The missing assumption: the tensor form of the Williams ansatz

Impose the ansatz as the paper already uses it implicitly:

```
F̄(r,θ) = I + K r^{α−1} ( G(θ) + o(1) ),   G(θ) = α Φ(θ)⊗ê_r + Φ'(θ)⊗ê_θ,
```

with `Φ` the Williams angular profile. In 2D, `det(a⊗ê_r + b⊗ê_θ) = a×b`, so

```
μ(θ) := det G(θ) = α ( Φ(θ) × Φ'(θ) ).
```

Two angular regimes, both *forced* by the ansatz (this resolves gap A7 — the
"one-large-stretch regime" is not a global property but a boundary-layer one):

- **Interior (generic θ): both stretches diverge.** Where `μ(θ) ≠ 0`,
  `det F̄ ≈ K² r^{2(α−1)} μ(θ) → ∞` and *both* singular values of `F̄` blow up
  like `r^{α−1}`. Since 2D cofactor is affine (`cof(I+M) = I + cof M`),

  ```
  ‖F̄^{-1}‖ = ‖cof F̄‖ / J  ≈  ‖G‖ (r/r_K)^{1−α} / μ(θ)  → 0   algebraically.
  ```

  (Checked numerically: `|F^{-1}|` tracks `r^{1−α}‖G‖/(K μ)` to 5% by
  `r = 10^{-4}`.) Orientation `J > 0` requires `μ(θ) ≥ 0`; this is a
  self-consistency constraint on the ansatz (cf. gap A4).

- **Degenerate rays: one-large-stretch layers.** `μ` vanishes at the Dirichlet
  face (`Φ(0) = 0` forces a second-order zero: `Φ×Φ' = (θ²/2)(a×b) + O(θ³)`)
  and possibly at finitely many interior angles. Near such a ray `G` is
  (nearly) rank-one, `F̄ = I + (big rank-one)`, and — provided
  `tr G ≠ 0` there, i.e. `J ~ K r^{α−1} tr G` — `λ_min` stays bounded below
  (checked: `λ_min → trG/‖G‖`-type constant). This is exactly the paper's
  one-large-stretch regime, now localised to layers

  ```
  Λ(r) = { θ : |μ(θ)| ≤ 2‖G‖_∞ (r/r_K)^{1−α} },     |Λ(r)| ≤ C (r/r_K)^{(1−α)/q₀}
  ```

  (`q₀` = maximal order of the zeros of μ; `q₀ = 2` at the clamped face).

**Nondegeneracy hypothesis (W):** μ ≥ 0 on `[0, π/2]`, with finitely many
zeros of finite order, and `tr G ≠ 0` at each zero. Generic for the
clamped–free Williams mode; checkable from the computed `Φ(θ; ν)`.

## 3. The screening lemma (the correct vertex limit)

From `eq:Lnh`, **the entire nonlinear part of 𝕃 carries two factors of
`F^{-T}`**:

```
𝕃(F̄) − μ_e 𝕀 = (f''J² + f'J) F̄^{-T}⊗F̄^{-T} − (f'J) F̄^{-T} ⊠ F̄^{-T},
```

so `‖𝕃(F̄) − μ_e 𝕀‖ ≤ Θ(J) ‖F̄^{-1}‖²`. Under mild growth
(`Θ(J) = O(log^p J)`) and hypothesis (W):

- **Off the layers** (`θ ∉ Λ(r)`):
  `‖𝕃(F̄) − μ_e 𝕀‖ ≤ C (log 1/r)^p (r/r_K)^{2(1−α)} μ(θ)^{−2} → 0`
  at an algebraic rate. The tangent *stabilises to the vector Laplacian*.
  For `f_log` the exact cancellation `f''J² + f'J = 2λ_e` (verified
  numerically) makes the rank-one term's coefficient *constant*; only the
  transposition term carries the log, and it dies under the `‖F^{-1}‖²`
  factor regardless.
- **On the layers:** (NC) holds with a computable `c`, and
  `‖𝕃‖ ≤ μ_e + C (log 1/r)^p` (the paper's `lem:Abound`, now valid exactly
  where its kinematics holds), on an angular set of algebraically vanishing
  width.
- Consequently `𝕃(F̄(r,·)) → μ_e 𝕀` in `L^q(0, π/2)` for every `q < ∞`, at an
  algebraic rate, though not in `L^∞`.

**Strong ellipticity bonus (affects case (B) and `rem:core`).** The
Legendre–Hadamard form is `μ_e |m|²|N|² + f''J² ξ²` with
`|ξ| ≤ |m||N|‖F̄^{-1}‖`. Off the layers, `|f''J²| ‖F̄^{-1}‖² ~
(log 1/r)(r/r_K)^{2(1−α)} → 0` even where `f''< 0`. So for `f_log` with
`λ_e > 0` **the dilation core is not a disc `{r ≤ r_*}`**: ellipticity can
fail only on the union of (i) the degenerate-ray layers and (ii) a bounded
transition collar `r ~ r_K` where the field is O(1)-nonlinear — a set of
vanishing angular measure near the vertex, not a full neighbourhood of it.
The "open inner problem" of `rem:core` shrinks accordingly.

## 4. The repaired theorem

Under Assumptions reg/growth/mild + (NC) + (W):

1. **Vertex pencil = mixed D–N vector Laplacian.** Off the layers the
   coefficient perturbation `𝕃 − μ_e 𝕀` decays algebraically, a textbook
   KMR-admissible perturbation; the layer contribution is `O(log 1/r)` on
   angular sets of width `O((r/r_K)^{(1−α)/q₀})` and is absorbable by a
   Meyers/higher-integrability perturbation argument (smallness in `L^q_θ`
   per dyadic annulus — the one remaining lemma-level technical step, see §5).
   The pencil of the limit operator `μ_e Δ` (components decouple; traction of
   the limit is `μ_e ∂_N`) on the quarter-plane with Dirichlet at `θ=0`,
   Neumann at `θ=π/2` has spectrum `s = 2k+1, k ≥ 0`: **no exponent in
   (0,1)** — the Newton-step corrections are Lipschitz at the vertex.
2. **Intermediate asymptotics = linear elasticity.** For `r_K ≪ r ≪ ρ`,
   `F̄ ≈ I` and `𝕃 ≈ 𝕃_lin`: the Williams exponents `α(ν) ∈ (0,1)` govern this
   range — which is the mesh- and matching-relevant one, and, since the inner
   exponents all exceed 1 > α, the `r^α` field remains the *dominant visible*
   asymptotics seen from any fixed scale. This is the precise surviving
   content of "passes like linear elasticity"; at the vertex itself the
   corner passes *better* than linear elasticity (extreme stretch screens the
   volumetric coupling — all `F`-dependence of the neo-Hookean tangent enters
   through `F^{-1} → 0`).
3. **Two-scale structure.** The Newton operator is a corner operator with a
   crossover at `r_K`: outer linear-elastic pencil, inner Laplacian pencil.
   Fredholmness in `V^{1,2}_β` holds for β admissible for *both* pencils; the
   energy weight `β = 0` is admissible for both (Laplacian interval is wide:
   `(−1, 1)`-type; elastic interval contains 0 since `α > 0`).

**Downstream (prop:blowup(ii)).** The wedge operator `𝓛_W` inherits the same
two-scale structure with `K̂` in place of `K`. Universality of `(α, Φ)` in
`(ν, K̂)` survives *as intermediate-range data*, which is what the numerics
resolve; the claim that all three operators "share the same leading vertex
pencil" should be dropped. Note also that at the *compressive* corner (where
`prop:blowup` is applied) `lem:Abound` is not available at the tip at all —
(NC) fails there — so the vertex claim was doubly unsupported; under
compression there is no screening (`J → 0` drives `f'' ~ J^{−2}` up through
the inverse channel) and the inner theory stays open (`rem:core`), as the
paper already concedes at the end of `prop:blowup(ii)`.

## 5. What remains open after the repair (much smaller)

1. **Layer absorption lemma.** Pencil/Fredholm stability under a perturbation
   bounded by `ε(r) + (log 1/r)^p 1_{Λ(r)}(θ)` with `|Λ(r)| → 0`
   algebraically. Standard route: Meyers higher integrability for the
   Laplacian-dominated operator on dyadic annuli + Hölder on the layer;
   should be a page, but it is genuine work not covered by a citation.
2. **Dirichlet-face transmission.** The clamped-face layer is not a small
   perturbation *of the boundary condition*: inside `θ ≲ (r/r_K)^{(1−α)/2}`
   the field is moderate and the operator is elastic-like. The claim that the
   collapsing elastic layer transmits the Dirichlet condition to the Laplacian
   bulk problem (so that the inner pencil is D–N as asserted) is a corner
   singular-perturbation statement. Plausible (the layer is asymptotically
   thin and uniformly elliptic), but it is the second place real analysis is
   owed.
3. **Self-consistency of the ansatz (gap A4 sharpened).** The state solves the
   nonlinear problem; with inner exponents `{1,3,…}` the self-consistent inner
   expansion of the *state* need not be the linear Williams field on
   `r ≪ r_K`. Hypothesis (W) should therefore be stated as a hypothesis on the
   state, with the linear Williams field justifying it only on `r ≳ r_K`.
   (Compare nonlinear crack-tip literature: Knowles–Sternberg-type inner
   fields differ from the linear ones — same phenomenon.)
4. **(W) verification.** `μ(θ) = α Φ×Φ' ≥ 0` with finite-order zeros and
   `tr G ≠ 0` at them is a one-line numerical check on the already-computed
   Williams profile `Φ(θ; ν)` — worth adding to the numerics section.

## 6. Status of the paper edits

Implemented in this pass (see `main_text.tex`, `concrete_realization.tex`):
new screening lemma + hypothesis (W); case (A) statement and proof rewritten
to the two-scale form (vertex pencil = Laplacian, linear-elastic exponents as
intermediate asymptotics); cap exemplar replaced by the pure log barrier and
the cap/(NC) incompatibility noted; `rem:core` updated (core = layers +
transition collar, not a disc); `prop:blowup(ii)` reworded to
intermediate-range universality. The two owed lemmas (layer absorption,
Dirichlet-face transmission) are flagged explicitly in the text as open
technical steps — but they are *bounded* debts, unlike the original step,
which was circular for every `K ≠ 0`.

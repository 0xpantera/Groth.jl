# Groth.jl — Package Reference & Notes

This is the canonical reference for the repository. Each package section lists
its purpose, key modules, noteworthy implementation details (original notes are
annotated rather than discarded), and follow-ups.

## GrothAlgebra

**Purpose**
- Prime-field arithmetic (BN254, secp256k1, and `GaloisField{p}`) plus generic
  polynomial and group utilities consumed by higher-level packages.

**Key modules**
- `src/FiniteFields.jl` — prime-field semantics plus the current default
  `BigInt` backend, with normalization, `invmod`, `powermod`, canonical integer
  conversion, display helpers, and `GaloisField{p}` for prime `p`.
- `src/Polynomial.jl` — Degree/leading-coefficient management, Horner
  evaluation, Lagrange interpolation, derivative, FFT scaffolding.
- `src/Group.jl` — Generic group/curve interface with scalar multiplication,
  w-NAF utilities, a Pippenger-backed variable-base MSM, and fixed-base helpers.

**Implementation notes**
- (original) FFT helpers were stubs; dense interpolation handled `h(x)`.
- (2025-09-29) Added `interpolate_prefix_points` so subset domains recover
  coefficients before padding for coset FFTs.

**Follow-ups**
- Evaluate FFT twiddle caching / mixed-radix support once the QAP domain aligns
  with arkworks.

## GrothCurves

**Purpose**
- BN254 tower fields, curve arithmetic, and pairing engine integration.

**Key modules**
- `src/BN254Fp2.jl`, `src/BN254Fp6_3over2.jl`, `src/BN254Fp12_2over6.jl` —
  Extension tower with optimized multiply/square/inverse.
- `src/BN254Curve.jl` — G1/G2 Jacobian ops, generators, checks.
- `src/BN254MillerLoop.jl`, `src/BN254FinalExp.jl`, `src/BN254Pairing.jl` —
  Optimal ate pairing pipeline and pairing engine wrapper.

**Implementation notes**
- Sparse Fp12 placement and Frobenius corrections follow standard BN254 refs.
- Projective point abstraction (`ProjectivePoint{Curve,F}`) enables additional
  curve engines in future work.
- BN254 G1/G2 scalar multiplication now dispatches to tuned w-NAF windows by
  default; the generic binary fallback remains available for other group types.

**Follow-ups**
- Prototype a second curve/engine (e.g., BLS12-381) to validate abstractions.
- Investigate GT cyclotomic optimisations if benchmarks surface bottlenecks.

## GrothProofs

**Purpose**
- R1CS/QAP conversion, Groth16 setup/prove/verify, and supporting tests.

**Key modules**
- `src/R1CS.jl`, `src/QAP.jl`, `src/Groth16.jl`.
- Tests in `test/runtests.jl` and `test/random_circuits.jl`.

**Implementation notes**
- (original) QAP interpolation used dense Lagrange polys over `[1..n]`, with a
  future FFT path planned.
- (2025-09-29) Coset path is default; dense vs coset equality asserted. Subset
  domains recover coefficients via barycentric interpolation before FFT (dense
  fallback retained for tests).
- (TODO) Align evaluation domain population with arkworks (fill entire domain,
  drop barycentric shim).

**Follow-ups**
- Optional proof aggregation (arkworks-style) if bandwidth savings are needed.
- Execute the domain-alignment work tracked in `docs/ROADMAP.md`.

## GrothExamples

**Purpose**
- Educational Pluto notebooks matching the tutorial narrative.

**Key modules**
- `src/r1cs_qap_pluto.jl` — AbstractAlgebra-first toy R1CS → QAP walkthrough.
- `src/r1cs_qap_groth_pluto.jl` — GrothAlgebra/GrothProofs version of the same flow.

**Notes**
- Notebooks use the shared-environment Pluto pattern (`Pkg.activate(joinpath(@__DIR__, ".."))`).
- Keeping both notebooks allows side-by-side comparison of pedagogical and package-native derivations.

## Benchmarks

**Purpose**
- `run.jl` / `plot.jl` microbenchmarks with JSON/PNG artefacts.

**Notes**
- Latest run: `results_2025-09-29_121914.json` (coset path default). README
  records medians; PNGs regenerated from this run.

## Docs

- `docs/RareSkills_Groth16_Map.md` — RareSkills mapping (updated Sep 2025).
- `docs/Implementation_vs_Arkworks.md` — Compare implementation choices with
  arkworks.

## Roadmap snapshot (see `docs/ROADMAP.md`)

- Immediate: align QAP domain population with arkworks.
- Stretch: proof aggregation, pairing optimisations, second curve prototype.

## Backend Migration Notes

- The current BN254 field layer still uses a `BigInt` default backend, but
  Stage 1 of the Montgomery roadmap extracts the backend boundary so future
  fixed-width field implementations can slot in without tying higher-level
  semantics to `.value::BigInt`.

---

# Groth.jl — Repository Status and Pairings-Oriented Overview

This document summarizes what’s implemented and well-documented in the repository, clarifies the mathematics versus representation and algorithm choices (with an emphasis on BN254 extension fields and pairings), and outlines concrete next steps to complete a functional Groth16 stack. Inline references point to files in this repo.

## Executive Summary

- Solid foundation in finite fields, basic algebra, and polynomials: clear, documented implementations with tests in `GrothAlgebra`.
- BN254 curve, extension fields, and the BN254 optimal ate pairing are implemented end-to-end in `GrothCurves`:
  - Fp2/Fp6/Fp12 extension tower with explicit formulas and optimized squarings.
  - G1/G2 in Jacobian coordinates with addition/doubling.
  - Miller loop with D-twist line placement and Frobenius correction steps.
  - Final exponentiation (easy + hard parts) with Frobenius maps and precomputed twist constants.
  - A large test suite exists under `GrothCurves/test` (good sign of maturity).
- R1CS and QAP conversion are present and readable in `GrothProofs`, but QAP uses a simplified evaluation domain and the Groth16 is a demonstrative, simplified version (not the full scheme yet).

Overall, the pairing stack looks coherent and well structured; the proof system side needs domain/FFT upgrades and a production-level Groth16 (setup/prove/verify) to align with the mathematical design.

## What’s Well Implemented and Documented

- Finite fields and polynomials (GrothAlgebra)
  - `GrothAlgebra/src/FiniteFields.jl`: BigInt-backed prime fields (BN254 prime, secp256k1, and `GaloisField{p}`) with normalized arithmetic, inversion via `invmod`, power via `powermod`, plus helpful display and utilities. Clean docstrings and straightforward APIs.
  - `GrothAlgebra/src/Polynomial.jl`: Polynomials over a field with degree/leading coefficient handling, Horner evaluation, Lagrange interpolation, derivative, and a placeholder for FFT multiply. Well commented; operations are correct and readable.
- `GrothAlgebra/src/Group.jl`: A generic `GroupElem` interface with scalar multiplication, w-NAF utilities, a Pippenger-backed variable-base MSM with Straus fallback, and helpers. Clear documentation of expectations from concrete curve types.

- BN254 extension tower and curve arithmetic (GrothCurves)
  - `GrothCurves/src/BN254Fp2.jl`: Fp2 with u² = −1; real/imag accessors, conjugate, norm, inverse, Frobenius as conjugation. Docstrings explain representation and formulas.
  - `GrothCurves/src/BN254Fp6_3over2.jl`: Fp6 over Fp2 via v³ = ξ with ξ = 9 + u (D-twist nonresidue). Uses Karatsuba-like multiplication and optimized squaring. Includes a clear `inv` via norm decomposition.
  - `GrothCurves/src/BN254Fp12_2over6.jl`: Fp12 over Fp6 via w² = v. Optimized squaring and inversion, conjugation, and a baseline Frobenius (later enhanced by final exponentiation module).
  - `GrothCurves/src/BN254Curve.jl`: G1 (over Fp) and G2 (over Fp2) in Jacobian coordinates with addition/doubling, affine conversion, curve membership checks, and standard generators. Nicely separated and easy to follow.
  - `GrothCurves/src/BN254MillerLoop.jl`: Miller loop for the optimal ate pairing on BN254.
    - D-twist line placement with explicit coefficient ordering, mixed addition, and careful mapping to sparse Fp12 coordinates in `evaluate_line`.
    - Frobenius-based endomorphism on G2 using constants (`P_POWER_ENDOMORPHISM_COEFF_*`) from known references, and the NAF loop count for u.
  - `GrothCurves/src/BN254FinalExp.jl`: Final exponentiation split into easy and hard parts.
    - Precomputes sextic-twist Frobenius `γ` constants (GAMMA1/2/3) and provides `frobenius_p1/p2/p3` for Fp12.
    - Implements the standard BN254 exponentiation ladder (`exp_by_u`, y₀..y₆ combination) consistent with literature.
  - `GrothCurves/src/BN254Pairing.jl`: Composes Miller loop + final exponentiation into `optimal_ate_pairing`, plus batch pairing.

- R1CS and QAP (GrothProofs)
  - `GrothProofs/src/R1CS.jl`: Clean R1CS representation with constraint checking and an example circuit (r = x·y·z·u). Well-documented and approachable.
  - `GrothProofs/src/QAP.jl`: Readable R1CS→QAP conversion via Lagrange interpolation and direct polynomial division for h(x). The domain is simplified to integers `[1..n]` (works for pedagogy; needs true roots-of-unity domain for FFT acceleration later).
  - `GrothProofs/src/Groth16.jl`: A simplified, educational Groth16 (CRS, prove, verify) that demonstrates the flow but is intentionally not the full construction.

[... original deep-dive continues ...]

## BN254 Field Extensions, Tower, and Pairing: A Focused Overview

This section distinguishes between (a) the math we need, (b) representation options, and (c) algorithmic choices taken here, with references.

### The Math We Need

- Fields and Tower
  - Base field: Fp with BN254 prime p. See `BN254_PRIME` in `GrothCurves/src/BN254Fp2.jl` and field arithmetic in `GrothAlgebra/src/FiniteFields.jl`.
  - Quadratic extension: Fp2 = Fp[u]/(u² + 1). Elements are a + b·u.
    - Multiplication: (a₀ + a₁u)(b₀ + b₁u) = (a₀b₀ − a₁b₁) + (a₀b₁ + a₁b₀)u.
    - Conjugate: \(\overline{a + bu} = a - bu\).
    - Norm: \(N(a + bu) = a^2 + b^2\).
    - Inverse: \((a + bu)^{-1} = (a - bu)/(a^2 + b^2)\).
  - Fp6 = Fp2[v]/(v³ - ξ) with ξ = 9 + u; split into cubic extension over Fp2.
  - Fp12 = Fp6[w]/(w² - v); final extension needed for GT.

- Curve Groups
  - G1: E(Fp) defined by y² = x³ + b with b = 3.
  - G2: defined over Fp2, using the sextic twist.

- Pairing
  - Optimal ate pairing e : G1 × G2 → GT with embedding degree 12.
  - Follows standard loop count (six bits of BN parameter u in NAF form) with Frobenius corrections.

### Representation Options (semantics-preserving)

- Field elements stored as small static tuples (`SVector`) to keep heap allocations down.
- G1/G2 stored in Jacobian coordinates (X:Y:Z) with the curve parameter a = 0 simplifying formulas.
- Fp12 stored as `Fp6_c0 + Fp6_c1 · w` with Fp6 stored as three Fp2 elements.

### Algorithmic Choices (performance / implementation detail)

- Fp2/Fp6/Fp12 multiplication uses Karatsuba-like decompositions and leverages nonresidue relationships to cut multiplications.
- Miller loop uses sparse line evaluation and explicit negation steps to minimise Fp12 multiplies.
- Final exponentiation uses Frobenius maps for the easy part and a standard BN ladder for the hard part (`exp_by_u`, `frobenius_p1/p2/p3`).
- MSM and fixed-base tables currently live in `GrothAlgebra/src/Group.jl`; variable-base MSM now uses a Pippenger backend with a small-input Straus fallback, while fixed-base follow-up work remains on the roadmap.

## Where We Are vs. What’s Left for Groth16

- Groth16 pipeline (CRS, prove, verify) is educationally complete but needs FFT-backed domain alignment to match arkworks (work in progress).
- Coset FFT path is implemented and default; dense path survives only for assertions/tests.
- Prepared verifier matches arkworks’ prepared path (batched pairing with a single final exponentiation).
- Proof aggregation is still outstanding.

## Implementation-Focused Math Notes (with Repo Mapping)

- Fp2 arithmetic (BN254, u² = −1) — `BN254Fp2.jl`
  - \((a_0 + a_1 u)(b_0 + b_1 u) = (a_0 b_0 - a_1 b_1) + (a_0 b_1 + a_1 b_0) u)\).
  - Conjugate \(\overline{a_0 + a_1 u} = a_0 - a_1 u\) and norm \(N(a) = a_0^2 + a_1^2\).
- Fp6 over Fp2 via v³ = ξ (ξ = 9 + u) — `BN254Fp6_3over2.jl`
  - Elements \(c_0 + c_1 v + c_2 v^2\) with Karatsuba cross-terms folded via ξ.
- Fp12 over Fp6 via w² = v — `BN254Fp12_2over6.jl`
  - Multiplication reduces `(d0 + d1w)(e0 + e1w)` using the relation w² = v.
- G1/G2 Jacobian formulas — `BN254Curve.jl`
  - Doubling: `M = 3X^2`, `S = 4XY^2`, etc. Mixed addition leverages affine second operand.
- Miller loop — `BN254MillerLoop.jl`
  - D-twist line evaluation, sparse Fp12 multiplication, Frobenius corrections.
- Final exponentiation — `BN254FinalExp.jl`
  - Easy part (`(f^p^6 / f)` etc.) and hard part via Frobenius + `exp_by_u` combos.

## Practical Distinctions: Math vs. Representation vs. Algorithms

- **Math invariants** — group laws, pairing bilinearity, field identities.
- **Representation** — coordinate systems, tower shapes, data layout; doesn’t change correctness but impacts performance.
- **Algorithms** — Karatsuba, sparse Fp12 ops, NAF loops, w-NAF MSM, etc.

## Glossary of Notation

- Fp, Fp2, Fp6, Fp12: Base and extension fields.
- G1, G2: Curve groups used in Groth16.
- GT: Multiplicative subgroup of Fp12 where the pairing lands.

## Concrete Next Steps Checklist

- [x] Coset FFT path default with dense assertion.
- [ ] Align QAP domain population with arkworks.
- [x] Implement variable-base Pippenger MSM and route prover query MSMs through it.
- [ ] Extend MSM work with fixed-base Pippenger / further setup-side tuning if benchmarks justify it.
- [ ] Prototype second pairing engine (BLS12-381).
- [ ] Port proof aggregation.

## Pointers Into This Repo

- `GrothAlgebra/src/*.jl`
- `GrothCurves/src/BN254*.jl`
- `GrothProofs/src/*.jl`
- `GrothExamples/`
- `benchmarks/`

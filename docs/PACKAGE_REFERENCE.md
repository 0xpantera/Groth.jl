# Groth.jl — Package Reference & Notes

This is the canonical reference for the repository. Each package section lists
its purpose, key modules, noteworthy implementation details (original notes are
annotated rather than discarded), and follow-ups.

## GrothAlgebra

**Purpose**
- Prime-field arithmetic (BN254, secp256k1, and `GaloisField{p}`) plus generic
  polynomial and group utilities consumed by higher-level packages.

**Key modules**
- `src/FiniteFields.jl` — prime-field semantics with a fixed-width Montgomery
  backend for `BN254Fq` / `BN254Fr`, plus the existing `BigInt`-backed generic
  fields (`GaloisField{p}`, secp256k1), canonical integer conversion, and
  display helpers.
- `src/Polynomial.jl` — Degree/leading-coefficient management, Horner
  evaluation, Lagrange interpolation, derivative, cached evaluation domains,
  and FFT / inverse FFT helpers.
- `src/Group.jl` — Generic group/curve interface with scalar multiplication,
  w-NAF utilities, a Pippenger-backed variable-base MSM, and fixed-base helpers.

**Implementation notes**
- (original) FFT helpers were stubs; dense interpolation handled `h(x)`.
- (2025-09-29) Added `interpolate_prefix_points` so subset domains recover
  coefficients before padding for coset FFTs.
- (2026-04-01) Stage 2 of the Montgomery roadmap moved `BN254Fq` and
  `BN254Fr` to a 4-limb Montgomery representation behind the Stage 1 backend
  hooks. The higher layers still see canonical field semantics; compatibility
  shims remain for callers that still inspect `.value`.
- (2026-04-01) Stage 3 cached per-domain twiddle/bitreverse metadata and added
  internal in-place `fft!` / `ifft!` helpers so BN254 `Fr` polynomial and QAP
  coset paths reuse the new Montgomery backend without extra setup/copy work.

**Follow-ups**
- Replace the remaining `BigInt`-based inversion path for BN254.
- Keep prover-shaped MSM follow-through tied to the deterministic benchmark
  fixtures.

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
- (2026-04-01) Stage 4 of the Montgomery roadmap replaced `SVector`-backed
  `Fp2`/`Fp6`/`Fp12` storage with concrete field members and specialized
  nonresidue helpers for `xi = 9 + u` and `v = (0, 1, 0)`. The first Stage 4
  artifact (`2026-04-01_145350`) cut `Fp6 mul` from `2.324 us` to `0.533 us`,
  `Fp12 mul` from `19.527 us` to `1.721 us`, and full pairing from `10.409 ms`
  to `3.843 ms` relative to the Stage 3 profile baseline.
- (2026-04-01) Stage 5 rewrote `BN254Curve.jl` around shared Jacobian helpers,
  specialized `Fp2` squaring in the curve hot path, and a lower-allocation
  `batch_to_affine!` implementation for both G1 and G2.

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
- (2026-04-02) The current follow-through work has shifted from broad backend
  migration to the remaining prover and pairing specialization steps surfaced
  by the Stage 8A baseline.

**Follow-ups**
- Optional proof aggregation (arkworks-style) if bandwidth savings are needed.
- Execute the remaining targeted specialization work tracked in `ROADMAP.md`.

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
- The benchmark harness now includes a dedicated `bn254_curve_kernels` family
  for direct G1/G2 add, double, and `to_affine` timings, alongside the Stage 5
  `batch_to_affine!` normalization families.
- The focused backend profiles are now staged through `stage8a`, with
  deterministic prover fixtures and `_semantic` proof checks wired into the
  `prove_full` path.
- Latest Stage 8A artifact: `benchmarks/artifacts/2026-04-01_223859`, with
  `generated_24_constraints prove_full` at `28.873 ms`, `msm_b_g2` at
  `4.334 ms`, `h_msm` at `7.036 ms`, and `final_c` at `2.804 ms`.

## Docs

- `docs/RareSkills_Groth16_Map.md` — RareSkills mapping (updated Sep 2025).
- `docs/Implementation_vs_Arkworks.md` — Compare implementation choices with
  arkworks.

## Roadmap snapshot (see `ROADMAP.md`)

- Immediate: limb-native inversion, final-exponentiation specialization,
  extension-field hot paths, safe G2 GLV exposure, and prover-shaped MSM
  specialization.
- Stretch: proof aggregation, second curve prototype, then Stage 9
  parallelism / accelerators once the single-thread backend is tighter.

## Backend Migration Notes

- Stage 1 extracted the backend boundary so field semantics are no longer tied
  directly to `.value::BigInt`.
- Stage 2 now runs `BN254Fq` and `BN254Fr` on a 4-limb Montgomery backend.
  Generic `GaloisField{p}` and `Secp256k1Field` still use the default `BigInt`
  path.
- Canonical integer conversion remains stable via `convert(BigInt, x)`.

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

# Implementation vs. Arkworks

This reference maps each major subsystem in Groth.jl to the analogous pieces in
arkworks (`ark-ff`/`ark-poly`/`ark-groth16`) and calls out where we mirror or
intentionally diverge from their design.

## Evaluation Domains & FFTs

| Topic | Arkworks | Groth.jl |
| --- | --- | --- |
| Domain sizing | `EvaluationDomain::new(num_constraints + num_inputs)` ensures every slot is populated before an IFFT | Matches arkworks' shape: constraints first, public-input selector rows next, then zero padding to the next power of two |
| Coset shift | `EvaluationDomain::get_coset(F::GENERATOR)` with tracked offset/inverses | `get_coset(domain, default_coset_offset)` tracks offset, inverse, and `offset_pow_size` |
| Vanishing on coset | `domain.evaluate_vanishing_polynomial(g)` with cached powers | Uses the full-domain closed form `g^n - 1`, cached as a constant inverse over the shifted coset |

Aug 2025 refactor summary:
- `prove_full` uses the coset-only quotient path; `compute_h_polynomial` remains
  the checked dense/coset helper for tests and debugging.
- QAP conversion feeds full-domain evaluation vectors to IFFT directly, so no
  subset coefficient recovery helper is needed for the Groth16 path.

## Multi-Scalar Multiplication (MSM) & Precomputation

| Topic | Arkworks | Groth.jl |
| --- | --- | --- |
| MSM backend | `VariableBaseMSM` (Straus + endomorphism on BLS curves) | `GrothAlgebra.multi_scalar_mul` now uses a Pippenger-style variable-base backend with a small-input Straus fallback; w-NAF helpers remain shared. Endomorphism optimisations TBD |
| Fixed-base tables | `FixedBaseMSM` with windowing | `FixedBaseTable` in `GrothAlgebra/src/Group.jl`; benchmarks cover table build & batch multiply |

Planned work: evaluate window sizes and endomorphisms once we profile the coset
prover path.

## Curve Abstractions & Pairing Engine

| Topic | Arkworks | Groth.jl |
| --- | --- | --- |
| Curve traits | `AffineCurve`, `ProjectiveCurve`, generic engines | `ProjectivePoint{Curve,F}` abstraction; pairing engine (`BN254Engine`) mirrors arkworks’ API |
| Pairing pipeline | Optimal ate with sparse line placement, Frobenius corrections | Same formulas; files named `BN254MillerLoop.jl`, `BN254FinalExp.jl`, `BN254Pairing.jl` |
| Engine extension | Additional curves (BLS12-381, etc.) provided via feature flags | Planned; abstractions ready but only BN254 implemented today |

## Groth16 Pipeline

| Topic | Arkworks | Groth.jl |
| --- | --- | --- |
| R1CS → QAP | Domain filled completely, IFFT → FFT on coset | Same structure; constraints, public-input selector slots, and zero padding are all explicit before IFFT |
| Prover | Coset FFT path by default, dense path available for debugging | `prove_full` uses coset-only H computation; dense/coset parity is kept in explicit debug/test helpers |
| Prepared verifier | `PreparedVerifyingKey` with batched pairing | `prepare_verifying_key`, `prepare_inputs`, `verify_with_prepared` mirror arkworks and use the pairing engine |
| Aggregation | `groth16::aggregate_proofs` available | Not yet ported; on roadmap |

## Documentation Cross-links

- RareSkills mapping: see `docs/RareSkills_Groth16_Map.md` for how textbook
  chapters map to Groth.jl modules.
- Package reference: `docs/PACKAGE_REFERENCE.md` summarises modules and
  implementation notes.

## Upcoming Alignment Tasks

1. Mirror more of arkworks' domain helper caching where it materially helps
   benchmarked prover paths.
2. Continue prover-side MSM and fixed-base tuning against `prove_full`.
3. Re-run benchmarks after each high-leverage domain/prover change.
4. Port proof aggregation once the core prover path is stable.

Feel free to extend this comparison as additional features (e.g., aggregation,
poly-commitments) land in the Julia stack.

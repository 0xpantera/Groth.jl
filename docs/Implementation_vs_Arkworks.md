# Implementation vs. Arkworks

This reference maps each major subsystem in Groth.jl to the analogous pieces in
arkworks (`ark-ff`/`ark-poly`/`ark-groth16`) and calls out where we mirror or
intentionally diverge from their design.

## Evaluation Domains & FFTs

| Topic | Arkworks | Groth.jl |
| --- | --- | --- |
| Domain sizing | `EvaluationDomain::new(num_constraints + num_inputs)` ensures every slot is populated before an IFFT | Currently: minimal power-of-two ≥ `num_constraints`; subset circuits pad coefficients after barycentric interpolation. **Next step:** match arkworks by filling every slot so we can drop the barycentric shim |
| Coset shift | `EvaluationDomain::get_coset(F::GENERATOR)` with tracked offset/inverses | `get_coset(domain, default_coset_offset)` tracks offset, inverse, and `offset_pow_size` |
| Vanishing on coset | `domain.evaluate_vanishing_polynomial(g)` with cached powers | We FFT `t(x)` on the coset; for full domains we use closed form `g^n - 1`. After we align domains we can mirror arkworks’ cached evaluation |

Aug 2025 refactor summary:
- Coset path is default (`compute_h_polynomial` asserts coset/dense equality).
- Subset domains temporarily recover coefficients via barycentric interpolation
  before padding in coefficient space.
- Planned alignment: populate the entire domain as arkworks does so the FFT/IFFT
  pair is used directly with no coefficient recovery helper.

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
| R1CS → QAP | Domain filled completely, IFFT → FFT on coset | Same structure; subset handling currently requires barycentric recovery (to be removed when we align domains) |
| Prover | Coset FFT path by default, dense path available for debugging | Coset path now default, dense used only for assertions; matches arkworks once domain alignment lands |
| Prepared verifier | `PreparedVerifyingKey` with batched pairing | `prepare_verifying_key`, `prepare_inputs`, `verify_with_prepared` mirror arkworks and use the pairing engine |
| Aggregation | `groth16::aggregate_proofs` available | Not yet ported; on roadmap |

## Documentation Cross-links

- RareSkills mapping: see `docs/RareSkills_Groth16_Map.md` for how textbook
  chapters map to Groth.jl modules.
- Package reference: `docs/PACKAGE_REFERENCE.md` summarises modules and
  implementation notes.

## Upcoming Alignment Tasks

1. Populate the entire evaluation domain (constraints + inputs) before calling
   IFFT, eliminating the barycentric interpolation helper.
2. Mirror arkworks’ domain helpers (vanishing evaluations, cached twiddles) once
   the layout matches.
3. Re-run benchmarks after domain alignment to ensure performance remains stable.
4. Port proof aggregation once the core prover alignment is complete.

Feel free to extend this comparison as additional features (e.g., aggregation,
poly-commitments) land in the Julia stack.

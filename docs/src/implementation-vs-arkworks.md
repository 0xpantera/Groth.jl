# [Implementation vs Arkworks](@id implementation-vs-arkworks)

This page maps each Groth.jl subsystem to the analogous arkworks component and
highlights where we match or intentionally diverge.

```@contents
Pages = ["implementation-vs-arkworks.md"]
Depth = 2
```

## Evaluation Domains & FFTs

| Topic | Arkworks | Groth.jl |
| --- | --- | --- |
| Domain sizing | `EvaluationDomain::new(num_constraints + num_inputs)` fills every slot before the IFFT. | Matches arkworks' shape: constraint rows first, public-input selector rows next, then zero padding to the next power of two. |
| Coset shift | `EvaluationDomain::get_coset(F::GENERATOR)` tracks offsets and inverses. | `get_coset(domain, default_coset_offset)` stores the offset, its inverse, and `offset_pow_size`. |
| Vanishing polynomial on a coset | Cached evaluation via `domain.evaluate_vanishing_polynomial(g)`. | Uses the full-domain closed form `g^n - 1`, cached as a constant inverse over the shifted coset. |

**Refactor snapshot (Sep 2025):**

- `prove_full` uses the coset-only quotient path; `compute_h_polynomial` remains
  the checked dense/coset helper for tests and debugging.
- QAP conversion feeds full-domain evaluation vectors to IFFT directly, so the
  Groth16 path no longer needs subset coefficient recovery.

## Multi-Scalar Multiplication (MSM)

| Topic | Arkworks | Groth.jl |
| --- | --- | --- |
| Variable-base MSM | Straus / Pippenger variants with endomorphism support where safe. | `GrothAlgebra.multi_scalar_mul` uses a Pippenger-style backend with a small-input Straus fallback; BN254 G1 has measured single-scalar GLV and an explicit subgroup-owned GLV-MSM path for the prover H/L query, while G2 exposes scalar GLV through an explicit subgroup-only helper rather than the arbitrary-point default. |
| Fixed-base tables | `FixedBaseMSM` window tables. | `FixedBaseTable` plus `build_fixed_table`, `mul_fixed`, and `batch_mul` mirror the workflow; setup uses measured BN254 scalar/GLV dispatch for G1 query generation and a fixed-window batch path for the G2 query. |

**Next steps:** keep MSM work tied to the real `prove_full` fixtures and avoid
treating synthetic MSM sweeps as the sole tuning signal.

## Curve Abstractions & Pairing Engine

| Topic | Arkworks | Groth.jl |
| --- | --- | --- |
| Curve traits | `AffineCurve`, `ProjectiveCurve`, generic engines. | `ProjectivePoint{Curve,F}` underpins `BN254Engine`; exposes the same pairing interface. |
| Pairing pipeline | Optimal ate with sparse line placement and Frobenius corrections. | Identical formulas, split across `BN254MillerLoop.jl`, `BN254FinalExp.jl`, `BN254Pairing.jl`. |
| Extensibility | Feature flags enable extra curves (e.g., BLS12-381). | Abstractions are in place; BN254 is live today and a second engine is on the roadmap. |

## Groth16 Pipeline

| Topic | Arkworks | Groth.jl |
| --- | --- | --- |
| R1CS → QAP | Domain fully populated, IFFT then FFT on the coset. | Same structure; constraints, public-input selector slots, and zero padding are all explicit before IFFT. |
| Prover/setup | Coset FFT path, dense available for debugging. | `prove_full` uses coset-only H computation and combines H/L into one subgroup-owned G1 GLV-MSM for `C`; `setup_full` uses measured query-generation dispatch; subgroup-owned fixed G2 elements use explicit GLV; dense/coset parity is kept in explicit debug/test helpers. |
| Prepared verifier | `PreparedVerifyingKey` batches pairings. | `prepare_verifying_key`, `prepare_inputs`, and `verify_with_prepared` mirror the API. |
| Aggregation | Optional `groth16::aggregate_proofs`. | Not yet ported; tracked on the roadmap. |

## Current Follow-Through Tasks

1. Replace the remaining `BigInt`-based inversion path in the BN254 Montgomery backend.
2. Specialize final exponentiation further around cyclotomic operations and measured BN254 structure.
3. Reduce value churn in extension-field hot paths where the current pairing code is still less in-place than arkworks.
4. Rebaseline `prove_full` after each high-leverage specialization stage instead of batching changes together.

Use this page as a checklist whenever behaviour changes—update the tables above as new features land.

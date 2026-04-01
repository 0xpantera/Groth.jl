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
| Domain sizing | `EvaluationDomain::new(num_constraints + num_inputs)` fills every slot before the IFFT. | Currently we pick the smallest power-of-two ≥ `num_constraints`. When the circuit does not span the full domain we recover the coefficients via barycentric interpolation before padding. **Next step:** populate every slot so the coset FFT can consume the coefficients directly. |
| Coset shift | `EvaluationDomain::get_coset(F::GENERATOR)` tracks offsets and inverses. | `get_coset(domain, default_coset_offset)` stores the offset, its inverse, and `offset_pow_size`. |
| Vanishing polynomial on a coset | Cached evaluation via `domain.evaluate_vanishing_polynomial(g)`. | For full domains we use the closed form `g^n - 1`. Subset domains FFT `t(x)` on the coset; once the padding matches arkworks we can reuse cached evaluations. |

**Refactor snapshot (Sep 2025):**

- Coset path is the default (`compute_h_polynomial` asserts coset vs dense parity).
- Subset domains temporarily rebuild coefficients via barycentric interpolation before the FFT padding step.
- Alignment plan: populate the full evaluation domain so the FFT/IFFT pair operates without the interpolation shim.

## Multi-Scalar Multiplication (MSM)

| Topic | Arkworks | Groth.jl |
| --- | --- | --- |
| Variable-base MSM | Straus / Pippenger variants with endomorphism support where safe. | `GrothAlgebra.multi_scalar_mul` uses a Pippenger-style backend with a small-input Straus fallback; BN254 G1 now has a measured hybrid GLV path, and G2 has an internal GLV path that is not yet the unconditional public default. |
| Fixed-base tables | `FixedBaseMSM` window tables. | `FixedBaseTable` plus `build_fixed_table`, `mul_fixed`, and `batch_mul` mirror the workflow; benchmarks track table build vs reuse. |

**Next steps:** keep MSM work tied to the real `prove_full` fixtures, pursue
safe G2 GLV exposure, and avoid treating synthetic MSM sweeps as the sole
tuning signal.

## Curve Abstractions & Pairing Engine

| Topic | Arkworks | Groth.jl |
| --- | --- | --- |
| Curve traits | `AffineCurve`, `ProjectiveCurve`, generic engines. | `ProjectivePoint{Curve,F}` underpins `BN254Engine`; exposes the same pairing interface. |
| Pairing pipeline | Optimal ate with sparse line placement and Frobenius corrections. | Identical formulas, split across `BN254MillerLoop.jl`, `BN254FinalExp.jl`, `BN254Pairing.jl`. |
| Extensibility | Feature flags enable extra curves (e.g., BLS12-381). | Abstractions are in place; BN254 is live today and a second engine is on the roadmap. |

## Groth16 Pipeline

| Topic | Arkworks | Groth.jl |
| --- | --- | --- |
| R1CS → QAP | Domain fully populated, IFFT then FFT on the coset. | Same structure; subset handling currently performs barycentric recovery before padding. |
| Prover | Coset FFT path, dense available for debugging. | Coset path is default; dense exists only for assertions. Will align fully once domains match. |
| Prepared verifier | `PreparedVerifyingKey` batches pairings. | `prepare_verifying_key`, `prepare_inputs`, and `verify_with_prepared` mirror the API. |
| Aggregation | Optional `groth16::aggregate_proofs`. | Not yet ported; tracked on the roadmap. |

## Current Follow-Through Tasks

1. Replace the remaining `BigInt`-based inversion path in the BN254 Montgomery backend.
2. Specialize final exponentiation further around cyclotomic operations and measured BN254 structure.
3. Reduce value churn in extension-field hot paths where the current pairing code is still less in-place than arkworks.
4. Find a safe way to expose G2 GLV acceleration on subgroup-owned points without weakening public semantics.
5. Rebaseline `prove_full` after each high-leverage specialization stage instead of batching changes together.

Use this page as a checklist whenever behaviour changes—update the tables above as new features land.

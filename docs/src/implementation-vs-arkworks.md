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
| Variable-base MSM | Straus with optional endomorphism on BLS curves. | `GrothAlgebra.multi_scalar_mul` implements Straus-style MSM; w-NAF helpers are shared with arkworks. Endomorphism optimisations are still TODO. |
| Fixed-base tables | `FixedBaseMSM` window tables. | `FixedBaseTable` plus `build_fixed_table`, `mul_fixed`, and `batch_mul` mirror the workflow; benchmarks track table build vs reuse. |

**Next steps:** profile the prover with window heuristics and consider endomorphisms once the coset path stabilises.

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

## Upcoming Alignment Tasks

1. Populate every evaluation domain slot (`num_constraints + num_inputs`) before the IFFT so we can delete the barycentric interpolation helper.
2. Mirror arkworks’ domain utilities (cached vanishing evaluations, twiddle reuse) once the layout matches.
3. Rebaseline benchmarks after the domain alignment to ensure prover performance stays on track.
4. Port proof aggregation only after the core prover is arkworks-aligned.

Use this page as a checklist whenever behaviour changes—update the tables above as new features land.

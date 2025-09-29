# Groth.jl Roadmap

This roadmap consolidates the previous planning notes (CODEX roadmap, polishing
plan, coset status, TODO) into a single reference. It tracks our progress toward
an arkworks-aligned, production-friendly Groth16 stack.

## Snapshot (September 2025)

- **GrothAlgebra** — Prime fields, polynomials, MSM helpers in place. Coset path
  now recovers coefficients via barycentric interpolation; FFT guard prevents
  truncation. Pending: align evaluation domain population with arkworks and then
  revisit FFT optimisations.
- **GrothCurves** — BN254 tower/curve/pairing fully implemented with pairing
  engine abstraction. Pending: second-curve prototype (e.g., BLS12-381) and
  optional GT optimisations.
- **GrothProofs** — R1CS/QAP/Groth16 path mirrors arkworks; coset FFT is the
  default. Dense path retained for assertions. Pending: domain alignment, optional
  proof aggregation, polishing docs.
- **GrothExamples** — Scripts and README now expose coset vs dense parity so the
  tutorial output matches implementation.
- **Benchmarks** — Run/plot workflow produces JSON + PNG artefacts; latest
  baseline `results_2025-09-29_121914.json` captures coset-default timings.

## Near-term Focus

1. **Align QAP domain with arkworks**
   - Populate all domain slots (`num_constraints + num_inputs`, next power of
     two) before IFFT.
   - Remove the barycentric interpolation helper, rely on native IFFT/FFT.
   - Update tests/examples/benchmarks accordingly.

2. **Tooling & verification hardening**
   - Run JET sweeps after major refactors, keep coset/dense parity assertions.
   - Extend documentation (package reference, RareSkills map, implementation vs
     arkworks) as behaviour changes.

## Upcoming Work

| Track | Items |
| --- | --- |
| Performance | Evaluate FFT twiddle caching, implement Pippenger-style variable/fixed-base MSM (with thresholds and w-NAF tables), keep batch normalisation fast, explore GT optimisations |

### MSM & Batch Optimisation Notes

- Variable-base MSM: introduce a Pippenger backend with window heuristics/thresholds and optional threading.
- Fixed-base MSM: build w-NAF tables for setup query generation (`FixedBaseTable`) and reuse in benchmarks.
- Batch normalisation: keep `batch_to_affine!` hot and normalise CRS vectors post-setup.
- Fallback strategy: keep naive loops for very small N and gate optimised paths behind thresholds.


| Features | Arkworks-style proof aggregation; second curve/pairing engine prototype |
| Tooling | Unified docs site (Documenter), CI matrix, contribution guide |
| Packaging | Compat bounds, CHANGELOGs, eventual registration |

## Completed Milestones

- Coset FFT path aligned (dense fallback removed, assertion on parity).
- Prepared verifier matches arkworks (batched pairing path).
- Benchmarks expanded (pairing & Groth16 hot paths with JSON/PNG artefacts).
- Groth16 tests cover multiple circuits, randomized seeds, prepared-path
  negatives, and randomized R1CS fixtures.

## References

- Package reference: `CODEX_ANALYSIS.md`
- Implementation vs arkworks: `docs/Implementation_vs_Arkworks.md`
- RareSkills mapping: `docs/RareSkills_Groth16_Map.md`

This roadmap should be updated as milestones land or priorities shift.

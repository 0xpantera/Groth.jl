# Groth.jl

Groth.jl is a Julia research platform for **BN254 algebra, pairings, and
Groth16**.

It combines:

- inspectable finite-field, extension-field, curve, and pairing code
- an end-to-end Groth16 setup / prove / verify pipeline
- benchmark and profiling infrastructure for prover hot paths
- Pluto-first educational material alongside performance-oriented engineering

This repository is **research-grade, not production software**. The goals are
clarity, mathematical correctness, REPL ergonomics, and serious performance
work without hiding the underlying algebra.

## What Is Implemented

- **GrothAlgebra**
  - BN254 `Fq` / `Fr` on a fixed-width Montgomery backend
  - generic polynomial utilities, cached evaluation domains, FFT / inverse FFT
  - variable-base MSM, fixed-base tables, and scalar-multiplication helpers
- **GrothCurves**
  - BN254 `Fp2` / `Fp6` / `Fp12`
  - G1 / G2 Jacobian arithmetic and batch normalization
  - optimal ate pairing with Miller loop and final exponentiation
- **GrothProofs**
  - R1CS and QAP conversion
  - Groth16 setup, proving, verification, and prepared-verifier flow
  - deterministic benchmark fixtures for `prove_full`
- **Benchmarks and docs**
  - reproducible JSON/PNG benchmark artifacts
  - primitive comparisons against `py_ecc` and local arkworks harnesses
  - roadmap and implementation notes tied to measured performance work

## Current Status

The project has moved well beyond a minimal Groth16 demo.

- BN254 primitives no longer run on the original `BigInt` hot path; the main
  backend is now Montgomery-based.
- Primitive benchmarks currently beat `py_ecc` across the tracked BN254 suite.
- The gap to arkworks has narrowed substantially, but arkworks is still ahead.
- QAP conversion now follows the arkworks domain shape: constraints first,
  public-input selector rows next, and zero padding to the next power of two.
- The current larger deterministic `prove_full` baseline after QAP domain
  alignment, coset-only H proving, and H/L MSM fusion is `26.643 ms` in the
  tracked summary
  [docs/src/assets/prove_full_msm_tuning_2026_05_11.json](./docs/src/assets/prove_full_msm_tuning_2026_05_11.json),
  down from the original `136.187 ms` baseline captured at the start of the
  performance investigation.
- The current larger deterministic `setup_full` fixture is `116.918 ms` in
  [docs/src/assets/setup_full_tuning_2026_05_11.json](./docs/src/assets/setup_full_tuning_2026_05_11.json),
  down from the `142.715 ms` pre-change baseline on the same fixture.
- Groth16 setup/proving now use an explicit G2 subgroup GLV helper for
  construction-owned key points, while generic G2 scalar multiplication remains
  the safe path for arbitrary verifier input.
- The active roadmap has shifted from broad backend replacement to targeted
  specialization: limb-native inversion, final-exponentiation work,
  prover-shaped MSM tuning, and then a fresh prover re-baseline.

See [ROADMAP.md](./ROADMAP.md) for the staged backend history and remaining
specialization work.

## Repository Layout

```text
Groth.jl/
├── GrothAlgebra/   # finite fields, polynomials, group utilities
├── GrothCurves/    # BN254 tower fields, curve arithmetic, pairing engine
├── GrothProofs/    # R1CS, QAP, Groth16 prover / verifier
├── GrothExamples/  # Pluto notebooks and walkthroughs
├── GrothCrypto/    # higher-level protocol space
├── benchmarks/     # BenchmarkTools environment, plots, profiling scripts
└── docs/           # reference docs, benchmarks, implementation notes
```

The sibling repositories in the workspace, such as `ark-works/`, `py_ecc/`,
and `zk-book/`, are reference checkouts. Active development happens in
`Groth.jl/`.

## Quick Start

```bash
# canonical workspace setup
julia --project=. -e 'using Pkg; Pkg.instantiate(workspace=true)'

# canonical full validation
julia --project=. scripts/test_all.jl

# package-scoped validation when intentionally narrowed
julia --project=GrothAlgebra -e 'using Pkg; Pkg.test()'
julia --project=GrothCurves -e 'using Pkg; Pkg.test()'
julia --project=GrothProofs -e 'using Pkg; Pkg.test()'

# benchmark harness
julia --project=. benchmarks/run.jl --list-profiles
julia --project=. benchmarks/run.jl --profile=quick
julia --project=. benchmarks/plot.jl

# docs
julia --project=docs docs/make.jl
```

Key notebooks live in `GrothExamples/`, starting with:

- `src/r1cs_qap_pluto.jl`
- `src/r1cs_qap_groth_pluto.jl`

## Documentation Map

- [ROADMAP.md](./ROADMAP.md) — staged BN254 backend roadmap and remaining work
- [benchmarks/README.md](./benchmarks/README.md) — benchmark methodology and
  current artifacts
- [docs/src/benchmarks.md](./docs/src/benchmarks.md) — docs-site benchmark page
- [docs/PACKAGE_REFERENCE.md](./docs/PACKAGE_REFERENCE.md) — package-level
  reference and repository notes
- [docs/Implementation_vs_Arkworks.md](./docs/Implementation_vs_Arkworks.md) —
  structural comparison with arkworks
- [docs/RareSkills_Groth16_Map.md](./docs/RareSkills_Groth16_Map.md) —
  textbook-to-code mapping

## Benchmark summary

The benchmark harness writes raw JSON/PNG artifacts under
`benchmarks/artifacts/`, which are intentionally ignored by git. The durable
external-comparison summary is tracked at
`docs/src/assets/external_benchmark_summary.json`, with narrative context in
`docs/src/benchmarks.md`.

Latest preserved `py_ecc` primitive comparison (`2026-04-01_174825`):

| Workload | Groth.jl | `py_ecc` | Result |
| --- | ---: | ---: | ---: |
| G1 scalar multiplication | `0.126 ms` | `0.264 ms` | Groth.jl `2.09x` faster |
| G2 scalar multiplication | `0.255 ms` | `1.492 ms` | Groth.jl `5.85x` faster |
| G1 naive accumulation, N=32 | `4.577 ms` | `12.452 ms` | Groth.jl `2.72x` faster |
| G2 naive accumulation, N=32 | `11.512 ms` | `69.238 ms` | Groth.jl `6.01x` faster |
| Single pairing | `3.140 ms` | `149.689 ms` | Groth.jl `47.67x` faster |

Tracked copies of the plots generated for that preserved artifact:

| Scalar multiplication | Pairing |
| --- | --- |
| ![Groth.jl vs py_ecc scalar multiplication](docs/src/assets/py_ecc_scalar_2026_04_01_174825.png) | ![Groth.jl vs py_ecc pairing](docs/src/assets/py_ecc_pairing_2026_04_01_174825.png) |
| G1 naive accumulation | G2 naive accumulation |
| ![Groth.jl vs py_ecc G1 naive accumulation](docs/src/assets/py_ecc_naive_accum_g1_2026_04_01_174825.png) | ![Groth.jl vs py_ecc G2 naive accumulation](docs/src/assets/py_ecc_naive_accum_g2_2026_04_01_174825.png) |

Refreshed local arkworks primitive comparison
(`2026-05-11_arkworks_bn254_refresh`):

| Workload | Groth.jl | arkworks | Result |
| --- | ---: | ---: | ---: |
| G1 scalar multiplication | `0.115 ms` | `0.00654 ms` | arkworks `17.58x` faster |
| G2 scalar multiplication | `0.223 ms` | `0.0170 ms` | arkworks `13.11x` faster |
| G1 naive accumulation, N=32 | `3.462 ms` | `0.182 ms` | arkworks `18.98x` faster |
| G2 naive accumulation, N=32 | `7.324 ms` | `0.481 ms` | arkworks `15.22x` faster |
| Single pairing | `2.960 ms` | `0.415 ms` | arkworks `7.14x` faster |

These are primitive-level measurements, not end-to-end Groth16 prover
comparisons.

Latest tracked `prove_full` prover fixture summary
(`2026-05-11_165756`, Stage 8 profile):

| Fixture | Domain | `prove_full` | Key prover phases |
| --- | ---: | ---: | --- |
| `sum_of_products_small` | `16` | `10.943 ms` | H+L MSM `5.619 ms`, final C `2.133 ms` |
| `generated_24_constraints` | `32` | `26.643 ms` | H+L MSM `11.029 ms`, final C `2.140 ms` |

The generated fixture improved from `28.636 ms` to `26.643 ms` versus the
previous coset-only H baseline by combining the separate H and L G1 MSMs into
one MSM for the `C` proof element.

Latest tracked `setup_full` fixture summary
(`2026-05-11_175228`, setup profile):

| Fixture | Domain | Baseline | Current |
| --- | ---: | ---: | ---: |
| `sum_of_products_small` | `16` | `47.910 ms` | `46.007 ms` |
| `generated_24_constraints` | `32` | `142.715 ms` | `116.918 ms` |

Setup now uses the BN254 G1 scalar dispatcher for G1 queries because the GLV
path beats fixed-base w-NAF on the measured full-width setup scalars; the G2
query keeps a fixed-window batch path. Fixed G2 key elements and the prover's
`delta_g2` randomizer term use an explicit subgroup-only GLV helper, preserving
the generic G2 scalar path for arbitrary on-curve points and verifier subgroup
checks.

## Performance Snapshot

Groth.jl now has two useful external reference points:

- **vs `py_ecc`**: current primitive benchmarks put Groth.jl ahead across the
  tracked BN254 scalar, accumulation, and pairing suite
- **vs arkworks**: Groth.jl has narrowed the gap sharply since the initial
  `BigInt` backend, but arkworks remains the stronger performance target

This repo is therefore best understood as:

- a serious Julia implementation of BN254 algebra, pairings, and Groth16
- a research and optimization platform
- not yet a drop-in replacement for a production Rust stack

## Development Notes

- Follow the repo and workspace `AGENTS.md` files.
- Use `execplans/` for non-trivial work.
- Keep benchmarks and docs in sync with user-visible behavior and measured
  performance changes.
- Prefer measured claims over aspirational ones.

## References

- [arkworks](https://github.com/arkworks-rs)
- [RareSkills Zero Knowledge Book](https://github.com/zkCollective/zk-book)

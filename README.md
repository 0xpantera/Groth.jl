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
- The current larger deterministic `prove_full` baseline is `28.873 ms` in
  [benchmarks/artifacts/2026-04-01_223859](./benchmarks/artifacts/2026-04-01_223859),
  down from the original `136.187 ms` baseline captured at the start of the
  performance investigation.
- The active roadmap has shifted from broad backend replacement to targeted
  specialization: limb-native inversion, final-exponentiation work, safer G2
  GLV exposure, prover-shaped MSM tuning, and then a fresh prover re-baseline.

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

```
                    ╔═══════════════════════════════════════════════════╗
                    ║                                                   ║
                    ║    ▄████  ██▀███   ▒█████  ▄▄▄█████▓ ██░ ██     ║
                    ║   ██▒ ▀█▒▓██ ▒ ██▒▒██▒  ██▒▓  ██▒ ▓▒▓██░ ██▒    ║
                    ║  ▒██░▄▄▄░▓██ ░▄█ ▒▒██░  ██▒▒ ▓██░ ▒░▒██▀▀██░    ║
                    ║  ░▓█  ██▓▒██▀▀█▄  ▒██   ██░░ ▓██▓ ░ ░▓█ ░██     ║
                    ║  ░▒▓███▀▒░██▓ ▒██▒░ ████▓▒░  ▒██▒ ░ ░▓█▒░██▓    ║
                    ║   ░▒   ▒ ░ ▒▓ ░▒▓░░ ▒░▒░▒░   ▒ ░░    ▒ ░░▒░▒    ║
                    ║    ░   ░   ░▒ ░ ▒░  ░ ▒ ▒░     ░     ▒ ░▒░ ░    ║
                    ║  ░ ░   ░   ░░   ░ ░ ░ ░ ▒    ░       ░  ░░ ░    ║
                    ║        ░    ░         ░ ░             ░  ░  ░    ║
                    ╚═══════════════════════════════════════════════════╝
                          ∴‥∵‥∴ zkSNARK SUPREMACY ∴‥∵‥∴

                    ░▒▓█►  ρяσνє єνєяутнιηg, яєνєαℓ ησтнιηg  ◄█▓▒░

         ╔═════════════════════════════════════════════════════════════╗
         ║  "In cryptography we trust, in zero-knowledge we thrive"   ║
         ║              ∞ milady privacy maximalism ∞                 ║
         ╚═════════════════════════════════════════════════════════════╝

                        ▓▓▓▓▓▓▓   ▓▓▓▓▓▓▓
                       ▓       ▓ ▓       ▓    ┌─────────────┐
                       ▓  ╳ ╳  ▓ ▓  ╳ ╳  ▓    │ COMMITMENT  │
                       ▓       ▓ ▓       ▓    │   HIDDEN    │
                        ▓     ▓   ▓     ▓     │  KNOWLEDGE  │
                         ▓▓▓▓▓     ▓▓▓▓▓      └─────────────┘
                           ║         ║
                        ═══╬═════════╬═══
                           ║         ║
                     ╔═════╩═════════╩═════╗
                     ║  TRUSTED SETUP CULT  ║
                     ╚═════════════════════╝
```

# Groth.jl

Research-focused, modular implementation of the Groth16 zero-knowledge proof system in Julia.
This repository explores cryptographic primitives (finite fields, curves, pairings) and the Groth16 protocol end-to-end in Julia. It is a research codebase and **not** production software.

## Repository layout

```
Groth.jl/
 ├── GrothAlgebra/   # finite fields, polynomials, group utilities
 ├── GrothCurves/    # BN254 curve + pairing engine
 ├── GrothProofs/    # R1CS, QAP, Groth16 prover/verifier
 ├── GrothExamples/  # Pluto notebooks and walkthroughs
 ├── GrothCrypto/    # placeholder for higher-level protocols
 ├── benchmarks/     # BenchmarkTools environment + plots
 └── docs/           # Roadmap, package reference, arkworks mapping
```

- `docs/PACKAGE_REFERENCE.md` — per-package summary, key modules, follow-ups
- `docs/Implementation_vs_Arkworks.md` — how our implementation compares to
  arkworks (domains, MSM, pairing, Groth16 pipeline)
- `docs/ROADMAP.md` — current priorities and completed milestones
- `docs/RareSkills_Groth16_Map.md` — link between the RareSkills ZK book and
  the Julia code

The sibling repositories `ark-works/` and `zk-book/` are reference checkouts
only; all active work lives in `Groth.jl/`.

## Getting started

```
# workspace setup (canonical)
julia --project=. -e 'using Pkg; Pkg.instantiate(workspace=true)'

# run all package tests (canonical)
julia --project=. scripts/test_all.jl

# package-scoped alternatives (when needed)
julia --project=GrothAlgebra -e 'using Pkg; Pkg.test()'
julia --project=GrothProofs -e 'using Pkg; Pkg.test()'

# benchmarks
julia --project=. benchmarks/run.jl
julia --project=. benchmarks/plot.jl
```

Key tutorials live in `GrothExamples/` as Pluto notebooks (starting with
`src/r1cs_qap_pluto.jl` and `src/r1cs_qap_groth_pluto.jl`).

## Project status

- **GrothAlgebra** — prime fields, polynomial utilities, MSM helpers. Coset FFT
  path uses coefficient-first padding; domain alignment with arkworks is the
  next major task.
- **GrothCurves** — BN254 tower, Jacobian curve ops, optimal ate pairing, and
  a pairing-engine abstraction ready for additional curves.
- **GrothProofs** — R1CS/QAP conversion, Groth16 setup/prove/verify mirroring
  arkworks. Coset FFT path is default; dense path retained only for assertions.
- **Benchmarks** — `benchmarks/results_2025-09-29_121914.json` captures the
  current coset-enabled timings; plots regenerated via `benchmarks/plot.jl`.

See `docs/ROADMAP.md` for detailed milestones and follow-up work (domain
alignment, MSM optimisations, second curve prototype, proof aggregation).

## Development guidelines

- Follow Julia style conventions: 4-space indent, `lowercase_with_underscores`
  for functions, docstrings on exported methods, and dispatch-friendly
  signatures that match existing APIs.
- Run the relevant `Pkg.test()` suites and update benchmarks when performance
  changes.
- Use concise imperative commit messages (e.g., `groth16: align coset domain`).
- Keep tutorial/docs in sync when user-visible output changes; the docs listed
  above are the canonical references.

## References

- [arkworks](https://github.com/arkworks-rs) (Groth16 and algebra implementations)
- [RareSkills Zero Knowledge Book](https://github.com/zkCollective/zk-book)

This repo explores cryptographic primitives (finite fields, curves, pairings) and the Groth16 protocol end-to-end in Julia.

Note: This project is for research and experimentation. It is not intended for production use.

## Project Structure

This is a monorepo containing several interconnected Julia packages:

### Core Packages

- **[GrothAlgebra](./GrothAlgebra)** - Foundation package with finite field arithmetic, polynomial operations, and group theory primitives
- **[GrothCurves](./GrothCurves)** - Elliptic curve implementations, focusing on BN254 (alt-bn128) with pairing support
- **[GrothProofs](./GrothProofs)** - Zero-knowledge proof systems including R1CS, QAP conversion, and Groth16 implementation
- **[GrothCrypto](./GrothCrypto)** - High-level cryptographic protocols built on top of the primitives
- **[GrothExamples](./GrothExamples)** - Educational Pluto notebooks and demonstrations

### Documentation & Tools

- **[docs/](./docs)** - Project documentation including the RareSkills→Groth16 map and polishing roadmap
- **[benchmarks/](./benchmarks)** - Performance benchmarks across packages

## Status Overview

- Algebra (GrothAlgebra)
  - Prime fields (Fp) and scalar field (Fr) implementations; polynomial arithmetic with interpolation and evaluation.
  - FFT/NTT and roots-of-unity domain planned (not yet implemented).
- Curves & Pairing (GrothCurves)
  - BN254 G1/G2 with Fp2/Fp6/Fp12 tower, optimal ate Miller loop, and final exponentiation.
  - Pairing engine abstraction (`AbstractPairingEngine`, `BN254Engine`) to support future curves while keeping BN254 optimized.
- Proofs (GrothProofs)
  - R1CS and QAP conversion in Fr; end-to-end Groth16 (CRS, prover, verifier) aligned with arkworks structure and equations.
  - Verifier enforces on-curve and subgroup checks; prepared verifier path batches pairings through the engine interface.
- Examples (GrothExamples)
  - Notebook-based walkthroughs (starting with toy R1CS -> QAP in Pluto).

## Roadmap Highlights

- Grow curve/generalization support: extend the pairing-engine interface to additional curves and trait-style abstractions.
- FFT/NTT + roots-of-unity domain for QAP and faster polynomial ops.
- MSM, fixed-base tables, and threaded hot paths for the prover plus batch normalization of CRS elements.
- Documentation polish: Documenter.jl site, Pluto notebooks, and contributor guides for multiple dispatch patterns.

## Pairing Engine Interface

- Depend on `AbstractPairingEngine` (re-exported by `GrothCurves`) when writing code that needs pairings. The interface consists of `miller_loop`, `final_exponentiation`, `pairing`, and `pairing_batch` methods specialised on an engine type parameterised by the curve family.
- `BN254Engine` is the default zero-sized backend; use the explicit form `pairing(BN254_ENGINE, P, Q)` or pass `engine=BN254_ENGINE` into helpers such as `GrothProofs.setup_full`. Convenience overloads without the engine argument remain available for interactive use.
- The engine test harness lives in `GrothCurves/test/test_pairing_engine_interface.jl` and checks zero-handling, bilinearity, and batch pairing consistency. Mirror those checks when introducing additional engines.

## Development

- This is a research project: APIs may evolve. Not for production use.
- Each package is a standalone Julia project with its own tests and docs.

## References

- arkworks-rs (groth16, algebra, relations)
- RareSkills ZK Book

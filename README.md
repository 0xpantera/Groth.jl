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
This repo explores cryptographic primitives (finite fields, curves, pairings) and the Groth16 protocol end-to-end in Julia.

Note: This project is for research and experimentation. It is not intended for production use.

## Project Structure

This is a monorepo containing several interconnected Julia packages:

### Core Packages

- **[GrothAlgebra](./GrothAlgebra)** - Foundation package with finite field arithmetic, polynomial operations, and group theory primitives
- **[GrothCurves](./GrothCurves)** - Elliptic curve implementations, focusing on BN254 (alt-bn128) with pairing support
- **[GrothProofs](./GrothProofs)** - Zero-knowledge proof systems including R1CS, QAP conversion, and Groth16 implementation
- **[GrothCrypto](./GrothCrypto)** - High-level cryptographic protocols built on top of the primitives
- **[GrothExamples](./GrothExamples)** - Educational examples and demonstrations

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
  - End-to-end demonstration of the Groth16 pipeline.

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

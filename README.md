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

- **[docs/](./docs)** - Overall project documentation and theory explanations
- **[benchmarks/](./benchmarks)** - Performance benchmarks across packages

## Status Overview

- Algebra (GrothAlgebra)
  - Prime fields (Fp) and scalar field (Fr) implementations; polynomial arithmetic with interpolation and evaluation.
  - FFT/NTT and roots-of-unity domain planned (not yet implemented).
- Curves & Pairing (GrothCurves)
  - BN254 G1/G2 with Fp2/Fp6/Fp12 tower, optimal ate Miller loop, and final exponentiation.
  - Comprehensive tests for field extensions and pairings.
- Proofs (GrothProofs)
  - R1CS and QAP conversion in Fr; end-to-end Groth16 (CRS, prover, verifier) aligned with arkworks structure and equations.
  - Verifier enforces on-curve and subgroup checks; tests include positive/negative cases for a sample circuit.
- Examples (GrothExamples)
  - End-to-end demonstration of the Groth16 pipeline.

## Roadmap Highlights

- Broaden Groth16 tests (multiple circuits, randomized r,s, multi-input IC).
- Prepared verifier (cache e(α,β); aggregate pairings via one multi-Miller loop + single final exp).
- FFT/NTT + roots-of-unity domain for QAP and faster polynomial ops.
- MSM and windowed scalar multiplication for the prover; batch normalization.

## Development

- This is a research project: APIs may evolve. Not for production use.
- Each package is a standalone Julia project with its own tests and docs.

## References

- arkworks-rs (groth16, algebra, relations)
- RareSkills ZK Book

# GrothCurves

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://0xpantera.github.io/GrothCurves.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://0xpantera.github.io/GrothCurves.jl/dev/)
[![Build Status](https://github.com/0xpantera/GrothCurves.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/0xpantera/GrothCurves.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Pairing engines

- `GrothCurves` exposes `AbstractPairingEngine{Curve}` plus the default `BN254Engine` for optimal ate pairings on BN254. Call `pairing(BN254_ENGINE, P, Q)` or reuse the exported convenience wrappers that omit the engine argument.
- Any new curve backend should implement `miller_loop`, `final_exponentiation`, and (optionally) `pairing_batch` for its engine type. See `src/BN254MillerLoop.jl` and `src/BN254Pairing.jl` for the reference implementation.
- The reusable interface tests live in `test/test_pairing_engine_interface.jl`; extend them when adding new engines to ensure bilinearity and batch consistency hold.

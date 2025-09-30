# Groth Curves

```@meta
CurrentModule = GrothCurves
DocTestSetup = :(using GrothCurves, Random)
```

```@contents
Depth = 2
```

## Overview

- BN254 tower fields, curve arithmetic, and the optimal ate pairing engine.
- Shared abstractions (`ProjectivePoint`, `BN254Engine`) prepare the ground for additional curves.

## Key Modules

- `src/BN254Fp2.jl`, `src/BN254Fp6_3over2.jl`, `src/BN254Fp12_2over6.jl` – extension tower with optimised multiply/square/inverse routines.
- `src/BN254Curve.jl` – G1/G2 Jacobian operations, generators, and curve membership checks.
- `src/BN254MillerLoop.jl`, `src/BN254FinalExp.jl`, `src/BN254Pairing.jl` – optimal ate Miller loop, final exponentiation, and pairing wrapper.

## Implementation Notes

- Sparse Fp12 placement and Frobenius corrections follow standard BN254 references.
- `ProjectivePoint{Curve,F}` keeps the pairing engine generic so future curves can reuse the pipeline.

## Follow-ups

- Prototype a second pairing engine (e.g., BLS12-381) to validate the abstractions.
- Investigate GT cyclotomic optimisations as benchmarks surface hotspots.

## Extension Fields

```@docs
Fp2Element
norm
conjugate
Fp6Element
Fp12Element
```

## Curve Operations

```@docs
ProjectivePoint
GrothCurves.doubling_step
GrothCurves.addition_step
GrothCurves.evaluate_line
```

### Example

```@example
using GrothCurves

P = g2_generator()
P2, coeffs_double = doubling_step(P)
Q = g2_generator()
PQ, coeffs_add = addition_step(P, Q)
(is_on_curve(P2), is_on_curve(PQ))
```

## Pairing Engine

```@docs
BN254Engine
GrothCurves.miller_loop
GrothCurves.final_exponentiation
GrothCurves.optimal_ate_pairing
GrothCurves.pairing
GrothCurves.pairing_batch
```

### Example

```@example
using GrothCurves

P = g1_generator()
Q = g2_generator()
pair_val = optimal_ate_pairing(P, Q)
batch_val = pairing_batch(BN254_ENGINE, [P], [Q])
(pair_val == batch_val)
```

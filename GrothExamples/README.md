# GrothExamples

Educational examples and demonstrations for the Groth zkSNARK ecosystem.

## Overview

This package contains practical examples demonstrating the use of the Groth16 implementation components:

- **Finite Field Arithmetic** - Working with BN254 and other fields
- **Polynomial Operations** - Interpolation, evaluation, and division
- **R1CS Construction** - Building constraint systems
- **QAP Conversion** - Converting R1CS to Quadratic Arithmetic Programs
- **Proof Generation** - Creating zero-knowledge proofs
- **Verification** - Verifying proofs using pairings

## Examples

### R1CS to QAP Demonstration
```julia
using GrothExamples
demonstrate_r1cs_qap()
```

### Multiplication Proof
```julia
using GrothExamples
multiplication_proof_example()
```

## Learning Path

1. Start with `test_r1cs_qap.jl` to understand the R1CS to QAP conversion
2. Study `multiplication_proof.jl` for a complete proof generation example
3. Explore the test files in each package for more detailed usage

## Dependencies

- GrothAlgebra - Core algebraic structures
- GrothCurves - Elliptic curve implementations
- GrothProofs - Proof system components
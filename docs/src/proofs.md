# Groth Proofs

```@meta
CurrentModule = GrothProofs
DocTestSetup = :(using GrothProofs, GrothAlgebra, GrothCurves, Random)
```

```@contents
Pages = ["proofs.md"]
Depth = 2
```

## Overview

- R1CS/QAP conversion plus Groth16 setup, proving, and verification.
- A seeded test suite exercises multiple circuit families and negative paths.

## Key Modules

- `src/R1CS.jl` – rank-1 constraint systems and reusable fixtures.
- `src/QAP.jl` – interpolation, evaluation domains, and the coset quotient helpers.
- `src/Groth16.jl` – setup/prove/verify APIs and prepared verifier flow.

## Implementation Notes

- Coset FFT path is the default; dense vs coset parity is asserted to guard regressions.
- Subset domains recover coefficients via barycentric interpolation before the FFT padding step until the domain alignment work completes.

## Follow-ups

- Align QAP domain population with arkworks, removing the interpolation shim.
- Investigate arkworks-style proof aggregation if bandwidth becomes a constraint.

## R1CS and QAP

```@docs
R1CS
Witness
r1cs_to_qap
QAP
evaluate_qap
compute_h_polynomial
polynomial_division
```

## Groth16 Pipeline

```@docs
Groth16Proof
ProvingKey
VerificationKey
PreparedVerificationKey
setup_full
prove_full
verify_full
prepare_verifying_key
prepare_inputs
verify_with_prepared
validate_witness_shape
public_inputs_from_witness
setup
prove
process_vk
verify
verify_prepared
```

### Example

```@example
using GrothProofs
using GrothAlgebra
using Random

r1cs = create_r1cs_example_multiplication()
witness = create_witness_multiplication(3, 5, 7, 11)
qap = r1cs_to_qap(r1cs)
# NOTE: For real deployments, use a CSPRNG instead of Random.GLOBAL_RNG.
keypair = setup_full(qap; rng = Random.MersenneTwister(42))
proof = prove_full(keypair.pk, qap, witness; rng = Random.MersenneTwister(1337))
public_inputs = r1cs.num_public > 1 ? witness.values[2:r1cs.num_public] : eltype(witness.values)[]
verify_full(keypair.vk, proof, public_inputs)
```

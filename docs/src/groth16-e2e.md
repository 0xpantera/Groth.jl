# Groth16 End-to-End

```@meta
CurrentModule = GrothProofs
DocTestSetup = :(using GrothProofs, GrothAlgebra, GrothCurves, Random)
```

This tutorial shows the high-level Groth16 flow in `GrothProofs`:

1. Build an R1CS and satisfying witness.
2. Run setup.
3. Generate a proof.
4. Verify with both raw and prepared verification keys.
5. Confirm failure on mutated public inputs.

## Conventions

- Witness layout is `[1, public..., private...]`.
- `public_inputs` excludes the leading `1` (arkworks-style).
- `validate_witness_shape` and `public_inputs_from_witness` enforce this contract.
- Current setup/prover randomness samplers are BN254Fr-oriented.
- For production usage, provide a CSPRNG-backed RNG.

## Complete Example

```@example
using GrothProofs
using Random

# 1) Build a circuit and witness.
r1cs = create_r1cs_example_sum_of_products()           # r = x*y + z*u
witness = create_witness_sum_of_products(3, 5, 7, 11)
validate_witness_shape(r1cs, witness)

# 2) Setup keys and QAP from the R1CS.
artifacts = setup(r1cs; rng=MersenneTwister(42), prepare_vk=true)

# 3) Prove.
proof = prove(artifacts.pk, artifacts.qap, witness; rng=MersenneTwister(1337))

# 4) Extract public inputs (without the leading constant 1) and verify.
public_inputs = public_inputs_from_witness(r1cs, witness)
ok_unprepared = verify(artifacts.vk, public_inputs, proof)
ok_prepared = verify(artifacts.pvk, public_inputs, proof)

# 5) Mutate one public input and verify rejection.
bad_inputs = copy(public_inputs)
bad_inputs[2] += one(eltype(bad_inputs))
bad_ok = verify(artifacts.vk, bad_inputs, proof)

(ok_unprepared, ok_prepared, bad_ok)
```

Expected result:

- `(true, true, false)`

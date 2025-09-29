
# GrothExamples

Educational scripts that walk through the Groth16 pipeline using the packages in
this monorepo.

## Available scripts

- `test_r1cs_qap.jl` — Builds a multiplication circuit, prints the evaluation
  domain/coset details, and shows that coset- and dense-derived `h(x)` match.
- `multiplication_proof.jl` — Runs the end-to-end Groth16 setup, prover, and
  verifier for the same circuit, printing public inputs and proof elements.

## Running examples

```julia
julia --project=GrothExamples src/test_r1cs_qap.jl
julia --project=GrothExamples src/multiplication_proof.jl
```

These scripts assume `GrothAlgebra`, `GrothCurves`, and `GrothProofs` are
checked out locally (see the repository README for the recommended `Pkg.develop`
workflow).

## Further reading

- RareSkills mapping: `docs/RareSkills_Groth16_Map.md`
- Package reference: `docs/PACKAGE_REFERENCE.md`
- Coset/FFT notes: scripts print parity information so the tutorial output is in
  sync with the implementation.


# GrothExamples

Notebook-first examples for the Groth.jl workspace.

## Available notebook

- `src/r1cs_qap_pluto.jl` — toy R1CS -> QAP walkthrough over `F_17` using `AbstractAlgebra`.
- `src/r1cs_qap_groth_pluto.jl` — the same walkthrough using `GrothAlgebra` and `GrothProofs`.

## Running

```bash
julia --project=GrothExamples -e 'using Pkg; Pkg.instantiate()'
julia --project=GrothExamples -e 'using Pluto; Pluto.run()'
```

Then open `GrothExamples/src/r1cs_qap_pluto.jl` from the Pluto file picker.


# GrothProofs

Constraint systems and Groth16 prover/verifier built on top of GrothAlgebra and
GrothCurves.

## Contents

- `src/R1CS.jl` — Rank-1 constraint system representation, satisfaction checks,
  and example circuits.
- `src/QAP.jl` — R1CS → QAP conversion with coefficient-first interpolation
  feeding the coset FFT pipeline.
- `src/Groth16.jl` — CRS setup, prover, verifier, and prepared verifier path.
- `test/runtests.jl`, `test/random_circuits.jl` — Regression suites covering
  multiple example circuits, randomized fixtures, and prepared-path negatives.

## Usage

```bash
# canonical monorepo setup + test commands
julia --project=. -e 'using Pkg; Pkg.instantiate(workspace=true)'
julia --project=. scripts/test_all.jl

# package-scoped alternative
julia --project=GrothProofs -e 'using Pkg; Pkg.test()'
```

```julia
using GrothProofs
r1cs = create_r1cs_example_sum_of_products()
qap = r1cs_to_qap(r1cs)
```

GrothProofs expects the pairing engine from GrothCurves and the algebra from
GrothAlgebra; in the monorepo this wiring is handled by the root workspace.

## Further reading

- Package reference: `docs/PACKAGE_REFERENCE.md` (GrothProofs section)
- Arkworks comparison: `docs/Implementation_vs_Arkworks.md`
- Roadmap status: `docs/ROADMAP.md`
- Tutorials: `GrothExamples/`

The top-level README outlines the preferred workflow and benchmarking commands.

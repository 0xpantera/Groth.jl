
# GrothAlgebra

Core algebraic primitives used by the Groth16 stack: prime fields, polynomial
arithmetic, and generic group/MSM utilities.

## Contents

- `src/FiniteFields.jl` — BigInt-backed prime fields (BN254, secp256k1, and
  `GaloisField{p}`) with normalized arithmetic and conversions.
- `src/Polynomial.jl` — Polynomial type with evaluation, interpolation,
  derivatives, and FFT scaffolding.
- `src/Group.jl` — Generic group API with scalar multiplication, w-NAF helpers,
  and multi-scalar multiplication utilities.
- `test/` — Unit tests covering field operations, polynomial identities, and
  group helpers.

## Usage

```julia
julia --project=GrothAlgebra -e 'using Pkg; Pkg.test()'   # run tests

julia --project -e 'using Pkg; Pkg.develop("GrothAlgebra")'
using GrothAlgebra
F = bn254_fq(5)
Polynomial([F(1), F(2), F(3)])  # quick REPL smoke test

GF7 = galois_field(7)
one_half = GF7(1) / GF7(2)
one_third = GF7(1) / GF7(3)
@assert one_half + one_third == GF7(5) / GF7(6)
```

## Further reading

- Package reference: `docs/PACKAGE_REFERENCE.md` (GrothAlgebra section)
- Arkworks comparison: `docs/Implementation_vs_Arkworks.md`
- Project roadmap: `docs/ROADMAP.md`

GrothAlgebra is developed alongside the other packages in this monorepo; see the
root README for the full workflow.


# GrothCurves

BN254 curve implementation, extension towers, and pairing engine machinery used
by Groth16.

## Contents

- `src/BN254Fp2.jl`, `src/BN254Fp6_3over2.jl`, `src/BN254Fp12_2over6.jl` —
  Fp2/Fp6/Fp12 tower with optimized multiplication, squaring, and Frobenius maps.
- `src/BN254Curve.jl` — G1/G2 Jacobian curve operations, generators, and
  subgroup checks.
- `src/BN254MillerLoop.jl`, `src/BN254FinalExp.jl`, `src/BN254Pairing.jl` —
  Optimal ate pairing pipeline and pairing engine wrapper (`BN254Engine`).
- `test/test_pairing_engine_interface.jl` — reusable bilinearity and batch
  pairing tests for any `AbstractPairingEngine` implementation.

## Usage

```julia
julia --project=GrothCurves -e 'using Pkg; Pkg.test()'   # run curve/pairing tests

julia --project -e 'using Pkg; Pkg.develop("GrothCurves")'
using GrothCurves
pairing(BN254_ENGINE, G1_GENERATOR, G2_GENERATOR)  # quick check
```

When adding new curves or engines:

1. Implement `miller_loop`, `final_exponentiation`, and (optionally)
   `pairing_batch` for the new engine type.
2. Extend the interface tests in `test/test_pairing_engine_interface.jl`.
3. Document the addition in `docs/PACKAGE_REFERENCE.md`.

## Further reading

- Package reference: `docs/PACKAGE_REFERENCE.md` (GrothCurves section)
- Arkworks comparison: `docs/Implementation_vs_Arkworks.md`
- Roadmap milestones: `docs/ROADMAP.md`

See the repository README for development workflow and related packages.

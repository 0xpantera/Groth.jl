# Implementation Notes

This page centralises repository-wide commentary that used to live in
`PACKAGE_REFERENCE.md` and the pairing-oriented status note. For package-specific
API details, see the dedicated chapters in this manual.

```@contents
Pages = ["implementation-notes.md"]
Depth = 2
```

## Repository Snapshot

- **GrothAlgebra** – prime-field arithmetic (BN254, secp256k1), polynomial
  utilities, and generic MSM helpers. BN254 `Fq` / `Fr` now run on a
  Montgomery backend, and the current follow-through work is focused on the
  remaining `BigInt` escape hatches and prover-shaped specialization.
- **GrothCurves** – BN254 tower fields, Jacobian G1/G2 arithmetic, and the
  optimal ate pairing pipeline. Sparse Fp12 placement and Frobenius corrections
  follow standard BN references; the `ProjectivePoint` abstraction keeps the
  pairing engine generic.
- **GrothProofs** – R1CS/QAP conversion plus Groth16 setup/prove/verify.
  Coset FFT is the default; the dense quotient remains solely as an assertion.
  The latest prover work has shifted from broad backend replacement to the
  remaining hot-path specialization tracked in the roadmap.
- **GrothExamples / benchmarks** – notebook-first tutorials (including
  AbstractAlgebra and package-native R1CS → QAP walkthroughs), alongside JSON/PNG
  benchmark artefacts capturing prover hot paths.

## Pairing-Oriented Overview

- The extension tower (`BN254Fp2`, `BN254Fp6_3over2`, `BN254Fp12_2over6`) uses
  Karatsuba-like decompositions and precomputed non-residue constants to shave
  multiplications.
- `BN254Curve.jl` keeps G1/G2 in Jacobian coordinates, exposing mixed additions
  against affine points plus batch normalisation for CRS vectors.
- `BN254MillerLoop.jl` implements the D-twist line placement with sparse Fp12
  multiplication, while `BN254FinalExp.jl` splits the easy/hard exponentiation
  and reuses Frobenius shortcuts. `BN254Pairing.jl` wraps everything under the
  `BN254Engine` abstraction so additional curves can slot in later.

## Proof System Overview

- `R1CS.jl` defines reusable fixtures (multiplication, sum-of-products,
  affine-product, square-offset) and powers the randomised circuit generator.
- `QAP.jl` records constraint points, builds power-of-two domains, and currently
  recovers coefficients via barycentric interpolation before padding.
- `Groth16.jl` wires the setup, prover, and verifier paths, including the
  prepared verifier that batches pairings (`pairing_batch`) before the single
  final exponentiation.

## Follow-ups from the Roadmap

- Replace the remaining `BigInt`-based inversion path in the BN254 backend.
- Specialize final exponentiation and other extension-field hot paths further.
- Rebaseline `prove_full` after each high-leverage specialization stage before
  deciding when Stage 9 parallelism is worth doing.

## References

- [Implementation vs Arkworks](@ref implementation-vs-arkworks) compares our
  choices against arkworks (`ark-ff`, `ark-poly`, `ark-groth16`).
- [RareSkills ↔ Groth.jl Map](@ref rareskills-map) links textbook chapters to the
  corresponding Julia modules.

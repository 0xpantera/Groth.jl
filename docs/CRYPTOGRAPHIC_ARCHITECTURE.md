# Groth.jl Cryptographic Architecture Map

This is an educational reference for understanding how the Groth.jl packages
fit together. It is not an API specification; prefer the package source,
tests, `ROADMAP.md`, and the Documenter pages under `docs/src/` when a precise
implementation detail matters.

Current status: May 2026. BN254 is implemented end to end. A second pairing
engine, such as BLS12-381, is still future work.

## How To Read This Map

Groth.jl is layered from algebraic primitives up to the Groth16 protocol:

```text
Groth16 setup/prove/verify
    ^
    | uses
R1CS and QAP
    ^
    | use
polynomials, FFT domains, fields
    ^
    | use
finite-field arithmetic

Pairing verification path:

Groth16 verifier
    ^
    | uses
pairing engine abstraction
    ^
    | implemented today by
BN254 tower fields, G1/G2 arithmetic, Miller loop, final exponentiation
```

The important split is:

- generic ZK infrastructure: fields, polynomials, FFTs, group algorithms, R1CS,
  and QAP machinery
- curve-specific cryptography: BN254 field constants, extension tower,
  G1/G2 equations, subgroup checks, GLV helpers, Miller loop, final
  exponentiation
- Groth16-specific protocol logic: CRS shape, proof elements, QAP query
  construction, prover equation, verifier equation

## Package Layers

### GrothAlgebra

`GrothAlgebra` owns reusable algebraic machinery:

- `FiniteFields.jl`
  - `BN254Fq` and `BN254Fr` use a fixed-width Montgomery backend.
  - Generic `GaloisField{p}` and other teaching-friendly fields remain
    available.
  - Canonical conversion through `BigInt` is still supported at API boundaries,
    but hot BN254 paths should avoid unnecessary `BigInt` round trips.
- `Polynomial.jl`
  - polynomial construction, evaluation, derivatives, division, interpolation
  - FFT/IFFT helpers and cached evaluation domains
  - coset-domain helpers used by the Groth16 quotient computation
- `Group.jl`
  - generic scalar multiplication
  - Pippenger-style variable-base MSM with small-input fallback
  - fixed-base tables and batch multiplication
  - BN254 scalar fast paths, including limb-native GLV decomposition where it
    is safe for subgroup-owned points

These pieces are intended to be reusable outside Groth16.

### GrothCurves

`GrothCurves` owns the currently implemented curve and pairing engine:

- `BN254Fp2.jl`, `BN254Fp6_3over2.jl`, `BN254Fp12_2over6.jl`
  - BN254 extension tower
  - concrete field-member storage
  - specialized nonresidue helpers and Frobenius constants
- `BN254Curve.jl`
  - G1 and G2 Jacobian point arithmetic
  - mixed additions, doubling, subgroup checks, batch normalization
  - explicit subgroup-only GLV helpers for paths that already own subgroup
    points
- `BN254MillerLoop.jl`, `BN254FinalExp.jl`, `BN254Pairing.jl`
  - optimal ate pairing for BN254
  - sparse line placement
  - cyclotomic final-exponentiation specialization
  - prepared G2 line coefficients for fixed verifier points
- pairing engine abstraction
  - `BN254Engine` is the live engine.
  - The abstraction is designed to admit another engine later, but BN254 is the
    only implemented engine today.

The verifier must keep arbitrary input validation conservative. Fast subgroup
helpers belong on construction-owned points or points that have already passed
the relevant subgroup discipline.

### GrothProofs

`GrothProofs` owns the proof-system layer:

- `R1CS.jl`
  - constraint representation
  - witness shape checks
  - small reusable circuit fixtures
- `QAP.jl`
  - R1CS to QAP conversion
  - arkworks-shaped domain population:
    - constraint rows first
    - public-input selector rows next
    - explicit zero padding to the next power of two
  - coset quotient path for the prover
- `Groth16.jl`
  - `setup_full`
  - `prove_full`
  - `verify_full`
  - prepared verifier flow: `prepare_verifying_key`, `prepare_inputs`,
    `verify_with_prepared`

Public inputs follow the arkworks convention: callers pass public inputs
without a leading `1`; the verifier key carries the constant `IC[1]` term.

Security note: the setup/prover APIs accept an RNG and default to
`Random.GLOBAL_RNG`. That is convenient for research and tests, but production
use should provide a CSPRNG-backed source of field randomness.

### GrothExamples

`GrothExamples` is notebook-first:

- `src/r1cs_qap_pluto.jl` is the more pedagogical AbstractAlgebra-oriented
  walkthrough.
- `src/r1cs_qap_groth_pluto.jl` shows the same ideas through Groth.jl packages.

The examples are for learning and exploration, not package test aggregation.

### Benchmarks And Docs

The benchmark workflow writes raw run artifacts under `benchmarks/artifacts/`.
Durable summaries and selected plots live under `docs/src/assets/` and are
referenced from:

- `docs/src/benchmarks.md`
- `README.md`

The active roadmap is root `ROADMAP.md`. The old staged Montgomery migration
log is archived at:

- `docs/history/bn254-montgomery-backend-roadmap.md`

The concept-to-code map is:

- `docs/src/rareskills-map.md`
- `docs/src/textbook-to-optimized.md`

Those pages map `zk-book` / RareSkills Groth16 concepts to Groth.jl modules and
explain where optimized implementation paths intentionally diverge from the
direct textbook shape.

## Reusable Vs Curve-Specific Vs Groth16-Specific

### Reusable Infrastructure

These pieces should largely survive a second curve or another pairing-based
protocol:

| Area | Groth.jl location | Notes |
| --- | --- | --- |
| finite-field interface | `GrothAlgebra/src/FiniteFields.jl` | reusable interface; concrete fast backends are field-specific |
| polynomials and FFT domains | `GrothAlgebra/src/Polynomial.jl` | generic over the field type |
| generic group algorithms | `GrothAlgebra/src/Group.jl` | reusable when group operations and identities are implemented |
| R1CS representation | `GrothProofs/src/R1CS.jl` | generic over the scalar field |
| QAP conversion shape | `GrothProofs/src/QAP.jl` | generic over the scalar field, but tied to R1CS/QAP-style systems |

### Curve-Specific Code

These pieces must be reimplemented or specialized for another pairing-friendly
curve:

| Area | BN254 location | Why curve-specific |
| --- | --- | --- |
| base/scalar fields | `BN254Fq`, `BN254Fr` | different moduli and backend constants |
| extension tower | `BN254Fp2/Fp6/Fp12` | different nonresidue constants and Frobenius coefficients |
| G1/G2 equations | `BN254Curve.jl` | different curve coefficients, generators, cofactors |
| subgroup checks | `BN254Curve.jl` | different subgroup/cofactor structure |
| GLV helpers | BN254 scalar/curve helpers | endomorphism parameters are curve-specific |
| Miller loop | `BN254MillerLoop.jl` | loop parameter and line placement are curve-specific |
| final exponentiation | `BN254FinalExp.jl` | exponent decomposition and addition chains are curve-specific |

### Groth16-Specific Code

These pieces are specific to Groth16:

| Area | Groth.jl location | Notes |
| --- | --- | --- |
| proving and verification keys | `GrothProofs/src/Groth16.jl` | CRS shape is Groth16-specific |
| proof tuple | `GrothProofs/src/Groth16.jl` | proof is `(A in G1, B in G2, C in G1)` |
| H/L query construction | `GrothProofs/src/Groth16.jl` | prover-specific use of QAP polynomials |
| verification equation | `GrothProofs/src/Groth16.jl` | pairing product equation for Groth16 |
| prepared verifier | `GrothProofs/src/Groth16.jl` | optimized representation of Groth16 verifier inputs |

## BN254 Flow

### Setup

```text
R1CS
  -> QAP over BN254Fr
  -> sample toxic waste scalars
  -> evaluate QAP queries at tau
  -> build proving key and verification key
  -> optionally prepare verification key
```

The current setup path uses measured BN254 scalar dispatch for G1 query
generation and a fixed-window batch path for G2 query generation. Construction
owned fixed G2 elements may use the explicit subgroup-only GLV helper.

### Proving

```text
proving key + witness
  -> validate witness shape and constraints
  -> compute A and B
  -> compute H through the coset quotient path
  -> combine H and private L contributions into one G1 MSM
  -> output A, B, C
```

The current prover uses a subgroup-owned G1 GLV-MSM helper for the combined H/L
MSM. Dense quotient computation remains useful as a test/debug parity check,
not as the production prover path.

### Verification

```text
verification key + public inputs + proof
  -> enforce public input length
  -> accumulate IC[1] + sum(public_i * IC[i + 1])
  -> check pairing product equation
```

The prepared verifier caches fixed G2 line coefficients and batches pairing
work where possible.

## What A Second Curve Would Require

BLS12-381 remains a useful future validation target, but it should be treated
as new curve engineering, not as a mechanical copy of BN254.

Expected new or heavily specialized pieces:

- base field and scalar field types
- extension tower constants and Frobenius coefficients
- G1/G2 curve definitions, generators, cofactors, subgroup checks
- curve-specific scalar multiplication and possible endomorphism support
- Miller loop with BLS12-381 loop parameter and twist-specific line placement
- final exponentiation for the BLS12-381 parameterization
- engine wrapper and pairing interface tests
- test vectors from independent references

Expected reusable pieces:

- polynomial and FFT infrastructure, once instantiated over the scalar field
- generic group APIs and fallback algorithms
- R1CS and QAP representation
- most high-level Groth16 API shape, assuming the proving and verification keys
  are parameterized by the new engine's G1/G2/GT types
- notebook and benchmark patterns

The main risk is not file count. The main risk is mathematical correctness in
curve-specific details: subgroup checks, cofactor clearing, twist placement,
Frobenius coefficients, and final-exponentiation addition chains.

## Concept Map

For studying from the local `zk-book/` checkout, use
`docs/src/rareskills-map.md` as the maintained map from concepts to code. The
high-level correspondence is:

| Concept | Groth.jl package |
| --- | --- |
| finite fields and modular arithmetic | `GrothAlgebra` |
| group theory and elliptic curve arithmetic | `GrothAlgebra`, `GrothCurves` |
| bilinear pairings | `GrothCurves` |
| R1CS | `GrothProofs` |
| Lagrange interpolation, QAP, R1CS to QAP | `GrothAlgebra`, `GrothProofs` |
| trusted setup and QAP evaluation | `GrothProofs` |
| Groth16 proving and verification | `GrothProofs` |
| worked examples | `GrothExamples` |

## Keeping This File Useful

Update this note when one of the following changes:

- a new package is added or removed
- a second curve engine lands
- the Groth16 API shape changes
- the benchmarked prover/setup paths switch algorithms
- the concept map moves or is renamed

Keep detailed API documentation in `docs/src/` and keep task-specific
execution history in `execplans/`. This file should stay educational and
architectural.

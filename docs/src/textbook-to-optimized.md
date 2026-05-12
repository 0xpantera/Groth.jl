# [From Textbook To Optimized Code](@id textbook-to-optimized)

Groth.jl started as a learning project: a way to turn RareSkills ZK Book
Groth16 concepts into concrete Julia code. The project still serves that role,
but it has also become a performance-oriented BN254 and Groth16 implementation.

That means the source code no longer reads like a direct transcription of the
book. The algebraic objects are the same, but the implementation often uses
representations and algorithms chosen for correctness checks, benchmarked hot
paths, and future extensibility.

```@contents
Pages = ["textbook-to-optimized.md"]
Depth = 2
```

## Two Roles

Groth.jl has two complementary roles:

- **Learning artifact:** it connects finite fields, polynomials, elliptic
  curves, pairings, R1CS, QAPs, trusted setup, proving, and verification to
  inspectable Julia modules.
- **Research implementation:** it uses Montgomery fields, cached FFT domains,
  coset quotient computation, Pippenger-style MSM, GLV-specialized subgroup
  paths, prepared verifier data, and measured benchmark fixtures.

The learning role explains why the modules are deliberately explicit. The
implementation role explains why many internal algorithms are no longer the
naive forms introduced first in a textbook.

## Reading Strategy

When mapping a concept from the book to the current codebase, use three layers:

1. **Concept:** the mathematical object or proof-system idea.
2. **Reference implementation:** the simple version that makes the concept
   visible.
3. **Optimized implementation:** the current production path inside Groth.jl.

For example, "multiply many CRS points by many scalars" is the concept. A
naive loop over scalar multiplication is the reference implementation. The
current prover uses a combined H/L G1 MSM with BN254 subgroup-owned GLV support
because that is the measured hot path.

## Concept Drift Map

| Book concept | Direct mental model | Current Groth.jl shape |
| --- | --- | --- |
| finite-field arithmetic | integers modulo a prime | `BN254Fq` / `BN254Fr` use fixed-width Montgomery representation |
| polynomial evaluation | plug `x` into coefficients | Horner evaluation plus FFT/IFFT domain workflows |
| roots of unity | multiplicative subgroup | cached evaluation domains with explicit coset metadata |
| R1CS | three sparse constraint vectors | typed constraints, witness shape checks, reusable fixtures |
| QAP | interpolate constraint columns | arkworks-shaped full-domain population before IFFT |
| quotient polynomial `H` | divide dense polynomials by `t(x)` | coset quotient path, with dense parity kept for tests/debugging |
| scalar multiplication | double-and-add | w-NAF, fixed-base tables, Pippenger-style MSM, BN254 GLV paths |
| pairing | bilinear map `e(P, Q)` | Miller loop, sparse line placement, final exponentiation |
| verifier equation | direct product of pairings | prepared verifier and batched pairing machinery |
| setup randomness | sample toxic waste scalars | RNG-injected API; convenient defaults for research, CSPRNG needed for production use |

## Examples

### Fields

The book-level model is:

```text
a + b mod p
a * b mod p
a^-1 mod p
```

The current BN254 implementation keeps the same field semantics, but the hot
representation is Montgomery limb arithmetic. That makes the code less
pedagogical at the storage level, while preserving the mathematical API:
addition, multiplication, inversion, exponentiation, and canonical conversion.

Use the simple mental model to understand correctness. Use the Montgomery model
to understand performance.

### Polynomials And Domains

The book-level model for QAPs is interpolation over constraint points, followed
by a divisibility check:

```text
P(x) = A(x)B(x) - C(x)
t(x) divides P(x)
H(x) = P(x) / t(x)
```

Groth.jl still enforces that algebra, but the prover path is shaped around FFT
domains and cosets. Constraint rows, public-input selector rows, and zero
padding are all placed explicitly before IFFT. The quotient is computed through
the coset path because that is the useful proving path; dense quotient
computation remains available as a parity check in tests and debugging helpers.

### MSM

The book-level model is:

```text
sum_i scalar_i * point_i
```

A naive implementation would loop over scalar multiplication. Groth.jl uses
more specialized algorithms because MSM dominates proving and setup work:

- Pippenger-style variable-base MSM
- fixed-base tables for repeated base points
- BN254 scalar dispatch
- explicit subgroup-only GLV helpers for points that Groth.jl construction
  already knows are in the correct subgroup

The important correctness boundary is subgroup ownership. The verifier keeps
generic, conservative paths for arbitrary input; construction-owned prover and
setup points can use faster subgroup-specialized paths.

### Pairings

The book-level model is a bilinear map:

```text
e(aP, bQ) = e(P, Q)^(ab)
```

The implementation is split into:

- Miller loop
- line evaluations
- sparse Fp12 multiplication
- final exponentiation
- prepared G2 line coefficients for fixed verifier points

This split is not a different concept. It is how the bilinear map is computed
efficiently and checked against pairing invariants.

### Groth16 Verification

The book-level model is the Groth16 pairing equation. Groth.jl keeps the same
equation, but the verifier has two shapes:

- direct verification
- prepared verification

The prepared path precomputes fixed verification-key data and batches pairings
where possible. This changes the operational shape, not the verifier equation.

## What Stays Conceptually Close

These parts remain good places to study the underlying concepts:

- `GrothProofs/src/R1CS.jl`
- `GrothProofs/src/QAP.jl`
- `GrothExamples/src/r1cs_qap_pluto.jl`
- `GrothExamples/src/r1cs_qap_groth_pluto.jl`
- the tests around QAP divisibility, coset parity, and full Groth16 examples

The optimized code is still worth reading, but these files are usually better
starting points when the goal is conceptual orientation.

## What Is Performance-Shaped

Expect more implementation distance from the book in:

- `GrothAlgebra/src/FiniteFields.jl`
- `GrothAlgebra/src/Group.jl`
- `GrothCurves/src/BN254Fp*.jl`
- `GrothCurves/src/BN254MillerLoop.jl`
- `GrothCurves/src/BN254FinalExp.jl`
- `GrothProofs/src/Groth16.jl` prover/setup hot paths

These files contain the same math, but they are organized around representation
choices, allocation control, subgroup safety boundaries, and benchmarked hot
paths.

## Related Maps

- [RareSkills ZK Book ↔ Groth.jl Map](@ref rareskills-map) maps conceptual
  chapters to modules.
- [Groth.jl Cryptographic Architecture Map](@ref architecture-map) separates
  reusable ZK infrastructure, curve-specific primitives, and Groth16-specific
  protocol code.
- [Implementation vs Arkworks](@ref implementation-vs-arkworks) compares
  Groth.jl's current implementation shape to arkworks.

Use this page when the code feels "too optimized" to match the book directly.
The goal is to preserve the conceptual path without forcing the implementation
to stay naive.

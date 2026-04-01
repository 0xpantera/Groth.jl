# BN254 Montgomery Backend Roadmap

This roadmap scopes the work needed to move Groth.jl from the current
`BigInt`-backed BN254 primitive layer toward a fixed-width Montgomery backend.
Each stage is intended to be small enough to become its own ExecPlan.

The main motivation is clear from the current primitive comparisons:

- against `py_ecc`, the remaining gap is concentrated in G2-heavy work
- against arkworks, the gap is much larger and shows up even on single
  primitive operations
- that points to the arithmetic representation itself, not just prover wiring

The current low-level reference points are:

- [FiniteFields.jl](./GrothAlgebra/src/FiniteFields.jl)
- [BN254Fp2.jl](./GrothCurves/src/BN254Fp2.jl)
- [BN254Fp6_3over2.jl](./GrothCurves/src/BN254Fp6_3over2.jl)
- [BN254Fp12_2over6.jl](./GrothCurves/src/BN254Fp12_2over6.jl)
- [BN254Curve.jl](./GrothCurves/src/BN254Curve.jl)
- [BN254MillerLoop.jl](./GrothCurves/src/BN254MillerLoop.jl)
- [BN254FinalExp.jl](./GrothCurves/src/BN254FinalExp.jl)

The key arkworks references are:

- [bn254 Fq](../ark-works/algebra/curves/bn254/src/fields/fq.rs)
- [bn254 Fq2](../ark-works/algebra/curves/bn254/src/fields/fq2.rs)
- [bn254 Fq12](../ark-works/algebra/curves/bn254/src/fields/fq12.rs)
- [bn254 G2](../ark-works/algebra/curves/bn254/src/curves/g2.rs)
- [generic prime fields](../ark-works/algebra/ff/src/fields/models/fp/mod.rs)
- [generic Fp2 model](../ark-works/algebra/ff/src/fields/models/fp2.rs)
- [generic variable-base MSM](../ark-works/algebra/ec/src/scalar_mul/variable_base/mod.rs)

## Main Goal

Replace the current BN254 primitive stack with a fixed-width Montgomery-backed
implementation that materially narrows the performance gap to arkworks while
preserving all algebraic and cryptographic invariants.

## Guiding Rules

- measure before changing representation code
- keep the current `BigInt` backend as the oracle during migration
- change one layer at a time: field, tower, curve, pairing, prover
- preserve public semantics until an intentional API change is documented
- do not mix backend replacement with unrelated prover optimizations

## Stage 0: Baseline And Oracle Harness

Goal:
Create the measurement and correctness harness that all later stages depend on.

Scope:

- keep the current BN254 primitive benchmarks current
- add direct benchmarks for `Fq`, `Fr`, `Fp2`, `Fp6`, `Fp12`, G1, G2, Miller
  loop, final exponentiation, and variable-base MSM
- keep the `py_ecc` and arkworks primitive comparisons reproducible
- define oracle checks against the current `BigInt` implementation

Deliverables:

- benchmark fixtures and artifact workflow for primitive comparisons
- cross-check tests between current and future backend results
- a short benchmark methodology note for repeatable local runs

Exit Criteria:

- primitive timing baselines are reproducible
- oracle comparisons exist for field, extension-field, curve, and pairing ops
- later stages can prove correctness without relying on output inspection alone

## Stage 1: Backend Interface Extraction

Goal:
Separate BN254 arithmetic semantics from the current `BigInt` storage model.

Scope:

- define the minimum internal interface needed by BN254 `Fq` and `Fr`
- isolate conversions, normalization, equality, serialization, and arithmetic
- keep the current backend as one implementation of that interface

Deliverables:

- a small backend boundary inside `GrothAlgebra`
- explicit compatibility shims for existing callers
- documentation of what is backend-specific versus BN254-specific

Exit Criteria:

- `BN254Fq` and `BN254Fr` semantics are no longer tied directly to `BigInt`
- the `BigInt` path still passes the full test suite
- the next stage can add a Montgomery backend without immediately touching
  curves or pairings

## Stage 2: BN254 Prime Field Montgomery Core

Goal:
Implement a fixed-width Montgomery backend for BN254 base and scalar fields.

Scope:

- implement limb-based representation for `Fq`
- implement limb-based representation for `Fr`
- support canonical decode/encode
- support add, sub, neg, double, mul, square, inversion, exponentiation, and
  conversion to canonical integers

Deliverables:

- a Montgomery core using fixed-width limbs
- `BN254FqMont` and `BN254FrMont` style internal types or equivalent
- oracle tests against the current `BigInt` backend
- primitive benchmarks for `Fq` and `Fr`

Exit Criteria:

- `Fq` and `Fr` operations match the oracle on deterministic and randomized
  tests
- base-field benchmarks improve materially over the current `BigInt` path
- no curve or pairing code has been migrated yet

## Stage 3: Polynomial And FFT Integration On `Fr`

Goal:
Move `BN254Fr`-based polynomial and domain code onto the new scalar-field
backend.

Scope:

- adapt polynomial arithmetic
- adapt roots of unity and domain setup
- adapt FFT and inverse FFT paths
- keep QAP and R1CS semantics unchanged

Deliverables:

- polynomial/domain compatibility with the new `Fr`
- benchmark deltas for FFT and polynomial work
- cross-checks for interpolation, evaluation, and vanishing-polynomial
  invariants

Exit Criteria:

- polynomial tests and QAP tests pass unchanged
- FFT/domain benchmarks remain correct and reproducible
- no silent degree, domain-size, or truncation regressions

## Stage 4: Extension Tower Migration

Goal:
Rebuild `Fp2`, `Fp6`, and `Fp12` on top of the Montgomery base field.

Scope:

- replace wrapper-heavy tower arithmetic with fixed-width field operations
- encode BN254-specific nonresidue choices explicitly
- add specialized Frobenius and conjugation paths
- preserve the current tower layout and math unless a measured improvement
  requires a documented change

Deliverables:

- Montgomery-backed `Fp2`, `Fp6`, and `Fp12`
- microbenchmarks for extension operations
- oracle and property checks against the current tower implementation

Exit Criteria:

- `Fp2`, `Fp6`, and `Fp12` operations match the oracle
- Frobenius-based identities hold
- extension-field benchmarks improve materially, especially on `Fp2` and `Fp12`

## Stage 5: G1 And G2 Curve Arithmetic Migration

Goal:
Move curve points and hot group operations onto the new field backend.

Scope:

- migrate G1 add, double, mixed-add, normalize, and scalar multiplication
- migrate G2 add, double, mixed-add, normalize, and scalar multiplication
- preserve subgroup behavior and zero-point conventions

Deliverables:

- curve arithmetic running on the new `Fq`/`Fp2`
- direct G1 and G2 primitive benchmarks
- subgroup, on-curve, and generator validation coverage

Exit Criteria:

- deterministic scalar multiplication outputs match the oracle
- G1 and G2 primitive benchmarks move materially in the right direction
- current serialization and subgroup expectations remain valid

## Stage 6: Pairing Engine Migration

Goal:
Move Miller loop and final exponentiation onto the new tower and curve backend.

Scope:

- migrate line-function arithmetic
- migrate Frobenius-based correction steps
- migrate final exponentiation easy and hard parts
- preserve bilinearity and non-degeneracy

Deliverables:

- Montgomery-backed Miller loop and final exponentiation
- pairing primitive benchmarks
- focused pairing invariants and reference checks

Exit Criteria:

- pairings match the oracle on deterministic vectors
- bilinearity and non-degeneracy tests pass
- primitive pairing benchmarks improve materially over the current path

## Stage 7: MSM And Scalar-Mul Specialization

Goal:
Exploit the new backend with better curve-level multiplication strategies.

Scope:

- retune single-point scalar multiplication on the new fields
- improve variable-base MSM on G1 and G2
- evaluate GLV-style acceleration where justified
- revisit fixed-base tables if measurements support it

Deliverables:

- post-backend scalar-mul benchmarks
- post-backend MSM benchmarks
- implementation notes for which specializations are worth keeping

Exit Criteria:

- MSM is no longer bottlenecked by the old arithmetic representation
- `prove_full` prover buckets tied to MSM improve materially
- optimizations remain mathematically transparent

## Stage 8: `prove_full` Integration And Regression Baseline

Goal:
Feed the new backend into the Groth16 prover and re-establish the baseline.

Scope:

- wire the Montgomery-backed primitives into the `prove_full` path
- re-run the staged `prove_full` benchmark and profiling workflow
- compare against previous `prove_full` artifacts and primitive deltas

Deliverables:

- a new end-to-end `prove_full` baseline
- profiler output tied back to the new primitive costs
- explicit notes on what the backend rewrite fixed and what remains

Exit Criteria:

- proof outputs remain correct
- total `prove_full` time drops materially
- hotspot ranking is updated using real measurements, not assumptions

## Stage 9: Parallelism, SIMD, And Optional Accelerators

Goal:
Use the new backend to unlock higher-level acceleration that was not worth doing
on the old representation.

Scope:

- threaded MSM
- task-level `prove_full` parallelism
- SIMD-oriented tuning where the limb layout supports it
- optional feature-gated accelerator work such as GPU MSM

Deliverables:

- parallel benchmark variants
- scalability notes by core count and workload size
- clear feature gates for non-portable acceleration

Exit Criteria:

- parallel speedups are measurable and reproducible
- the single-thread backend remains the correctness reference
- feature-gated acceleration does not complicate the default correctness path

## Portability Beyond BN254

Short answer:
yes, the backend design should be reusable for other curves, but not as a single
hardcoded BN254 type.

The right target is:

- one generic Montgomery prime-field engine
- per-field configuration for modulus, generator, limb count, and two-adicity
- per-extension configuration for nonresidues and Frobenius coefficients
- per-curve configuration for coefficients, generators, subgroup logic, and any
  endomorphism-based optimizations

This is exactly the broad arkworks strategy:

- prime fields use generic `Fp<N>` plus `MontBackend<Config, limbs>` in
  [fp/mod.rs](../ark-works/algebra/ff/src/fields/models/fp/mod.rs)
- each concrete field provides a small config, e.g. BN254 `Fq` in
  [fq.rs](../ark-works/algebra/curves/bn254/src/fields/fq.rs)
- tower fields are generic over configs, e.g. `Fp2Config` in
  [fp2.rs](../ark-works/algebra/ff/src/fields/models/fp2.rs)
- each curve then plugs in its own field configs and curve configs, e.g. BN254
  G2 in [g2.rs](../ark-works/algebra/curves/bn254/src/curves/g2.rs)

That is also why arkworks supports many curves:

- the Montgomery machinery is shared
- the curve-specific constants and formulas are localized
- the limb width changes with the field size, such as `Fp256`, `Fp384`,
  `Fp768`, and `Fp832`

Implication for Groth.jl:

- the first implementation should be BN254-focused to keep scope under control
- the architecture should avoid baking BN254 assumptions into the Montgomery
  core
- once BN254 is stable, the same core could support curves such as secp256k1,
  BLS12-381 base/scalar fields, or other pairing-friendly families with
  additional tower and curve-specific work

## Validation Rules

Every stage must preserve:

- exact finite-field results in canonical representation
- exact extension-field results and Frobenius identities
- exact G1/G2 scalar multiplication outputs on deterministic fixtures
- pairing bilinearity and non-degeneracy
- Groth16 proof correctness

Required checks after each implementation stage:

- `julia --project=GrothAlgebra -e 'using Pkg; Pkg.test()'`
- `julia --project=GrothCurves -e 'using Pkg; Pkg.test()'`
- `julia --project=GrothProofs -e 'using Pkg; Pkg.test()'`
- `julia --project=. scripts/test_all.jl`
- stage-appropriate primitive and prover benchmarks

## Out Of Scope For The Backend Track

- new proof-system features
- unrelated prover API redesign
- premature GPU work before the single-thread backend is stable
- switching curves mid-migration
- claiming performance wins without benchmark artifacts

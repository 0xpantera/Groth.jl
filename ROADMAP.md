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

## Current Status

Completed stages:

- Stage 0: baseline and oracle harness
- Stage 1: backend interface extraction
- Stage 2: BN254 prime-field Montgomery core
- Stage 3: polynomial and FFT integration on `Fr`
- Stage 4: extension tower migration
- Stage 5: G1 and G2 curve arithmetic migration
- Stage 6: pairing engine migration
- Stage 7: MSM and scalar-mul specialization
- Stage 7A: BN254 GLV scalar multiplication
- Stage 8: `prove_full` integration and regression baseline
- Stage 8A: BN254Fr-native prover scalar plumbing

Next concrete work:

- target the remaining measured `BigInt` / GMP work that still appears after
  Stage 8A, starting with limb-native inversion and final-exponentiation
  specialization
- then re-evaluate the next highest-value prover bottleneck from the updated
  `prove_full` profile instead of assuming the next win is only in MSM

Stage 8 remains the end-to-end `prove_full` integration and regression
baseline stage.

## Stage 0: Baseline And Oracle Harness

Goal:
Create the measurement and correctness harness that all later stages depend on.

Status:
Completed on 2026-04-01.

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

Notes:

- Stage 0 added deterministic primitive benchmark groups for `BN254Fq`,
  `BN254Fr`, `Fp2`, `Fp6`, `Fp12`, direct G1/G2 scalar multiplication, Miller
  loop, and final exponentiation.
- Stage 0 also added deterministic oracle tests for field, tower, scalar-mul,
  Miller-loop, and pairing outputs so later backend stages could preserve exact
  semantics while changing representation.
- The first quick Stage 0 artifact is
  `benchmarks/artifacts/2026-04-01_122211`.

## Stage 1: Backend Interface Extraction

Goal:
Separate BN254 arithmetic semantics from the current `BigInt` storage model.

Status:
Completed on 2026-04-01.

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

Notes:

- Stage 1 introduced an explicit internal backend hook layer in
  `GrothAlgebra/src/FiniteFields.jl`, including canonical conversion and
  normalized construction boundaries.
- Production callers that previously read `.value` directly were moved onto the
  abstraction boundary so the Montgomery backend could land without rewriting
  every higher layer at once.

## Stage 2: BN254 Prime Field Montgomery Core

Goal:
Implement a fixed-width Montgomery backend for BN254 base and scalar fields.

Status:
Completed on 2026-04-01 for `BN254Fq` and `BN254Fr`.

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

Notes:

- `BN254Fq` and `BN254Fr` now store 4 `UInt64` limbs in Montgomery form behind
  the Stage 1 backend hooks.
- Canonical conversion still flows through `convert(BigInt, x)`.
- The current inversion path decodes to canonical integers and uses `invmod`,
  while add/sub/mul/square stay on the fixed-width backend.
- Relative to the Stage 0 quick baseline (`2026-04-01_122211`), the
  `2026-04-01_140814` artifact shows:
  - `Fq` add/sub about `40-48x` faster
  - `Fq` mul/square about `14-15x` faster
  - `Fr` add/sub about `35-36x` faster
  - `Fr` mul/square about `15x` faster
  - G1 scalar about `1.39x` faster
  - G2 scalar about `3.50x` faster
  - Miller loop about `3.82x` faster
  - Full pairing about `5.20x` faster

## Stage 3: Polynomial And FFT Integration On `Fr`

Goal:
Move `BN254Fr`-based polynomial and domain code onto the new scalar-field
backend.

Status:
Completed on 2026-04-01 for the `EvaluationDomain`, FFT/IFFT helpers, and QAP
coset recovery path.

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

Notes:

- `EvaluationDomain` now caches stage roots and bit-reversal permutations so
  repeated NTT calls stop rebuilding the same metadata.
- `Polynomial.jl` now exposes internal in-place `fft!` / `ifft!` helpers and
  uses them to avoid extra copies in interpolation and FFT multiplication.
- `QAP.compute_h_polynomial` and the benchmarking helper path now reuse the
  in-place inverse FFT on the coset quotient buffer instead of copying again.
- Stage 3 introduced a dedicated `bn254_polynomials` benchmark family and plot
  so polynomial/FFT work can be measured directly without relying on
  `prove_full`.
- The first Stage 3 artifact is `benchmarks/artifacts/2026-04-01_142526`,
  with medians:
  - `EvaluationDomain(32)`: `7.261 μs`
  - `fft(32)`: `26.116 μs`
  - `ifft(32)`: `33.129 μs`
  - `interpolate_fft(16 -> 32)`: `28.576 μs`
  - `fft_polynomial_multiply(16x16)`: `87.894 μs`
  - `compute_h_polynomial`: `82.718 μs`

## Stage 4: Extension Tower Migration

Goal:
Rebuild `Fp2`, `Fp6`, and `Fp12` on top of the Montgomery base field.

Status:
Completed on 2026-04-01 for the BN254 `Fp2` / `Fp6` / `Fp12` tower, concrete
storage, and specialized nonresidue helpers.

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

Notes:

- `Fp2Element`, `Fp6Element`, and `Fp12Element` now store concrete field
  members instead of `SVector` wrappers.
- `Fp6` now uses a dedicated `mul_fp2_by_nonresidue` helper for `xi = 9 + u`,
  and `Fp12` uses a dedicated `mul_fp6_by_nonresidue` helper for `v = (0, 1,
  0)`.
- The first Stage 4 artifact is `benchmarks/artifacts/2026-04-01_145350`.
- Relative to the Stage 3 baseline, it moved:
  - `Fp6 mul` from `2.324 μs` to `0.533 μs`
  - `Fp6 square` from `11.919 μs` to `0.432 μs`
  - `Fp12 mul` from `19.527 μs` to `1.721 μs`
  - `Fp12 inv` from `47.426 μs` to `6.995 μs`
  - `G2 scalar` from `446.993 μs` to `258.031 μs`
  - full `pairing` from `10.409 ms` to `3.843 ms`

## Stage 5: G1 And G2 Curve Arithmetic Migration

Goal:
Move curve points and hot group operations onto the new field backend.

Status:
Completed on 2026-04-01 for the shared Jacobian kernels, direct curve-kernel
benchmarks, and lower-allocation `batch_to_affine!` paths.

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

Notes:

- Stage 5 added a dedicated `bn254_curve_kernels` benchmark family plus the
  `stage5` benchmark profile so curve arithmetic and normalization can be
  remeasured without running the full suite.
- The first Stage 5 artifact is `benchmarks/artifacts/2026-04-01_154740`.
- Relative to the Stage 4 artifact `2026-04-01_145350`, it moved:
  - `G1 scalar` from `190.685 μs` to `105.483 μs`
  - `G2 scalar` from `258.031 μs` to `226.369 μs`
  - full `pairing` from `3.843 ms` to `3.892 ms`
- The new direct curve-kernel medians are:
  - `G1 double`: `4.630 μs`
  - `G1 add`: `4.755 μs`
  - `G1 to_affine`: `2.361 μs`
  - `G2 double`: `8.544 μs`
  - `G2 add`: `7.687 μs`
  - `G2 to_affine`: `4.133 μs`
- For batch normalization at `N = 128`, the same artifact records:
  - `G1 batch_to_affine!`: `25.262 μs` vs per-point `386.911 μs`
  - `G2 batch_to_affine!`: `87.669 μs` vs per-point `658.845 μs`

## Stage 6: Pairing Engine Migration

Goal:
Move Miller loop and final exponentiation onto the new tower and curve backend.

Status:
Completed on 2026-04-01.

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

Notes:

- Stage 6 added direct `pairing_substeps` benchmarks for `doubling_step`,
  `addition_step`, `evaluate_line`, Frobenius helpers, `exp_by_u`,
  `final_exponentiation_easy`, and `final_exponentiation_hard`.
- The first Stage 6 substep artifact is
  `benchmarks/artifacts/2026-04-01_162018`.
- The final focused Stage 6 artifact is
  `benchmarks/artifacts/2026-04-01_163914`, with the direct Stage 5 comparison
  artifact at `benchmarks/artifacts/2026-04-01_164003`.
- Relative to the Stage 5 baseline `2026-04-01_154740`, Stage 6 improved:
  - `pairing`: `3.892 ms -> 3.107 ms`
  - `miller_loop`: `2.416 ms -> 1.782 ms`
  - `final_exponentiation`: `1.407 ms -> 1.227 ms`

## Stage 7: MSM And Scalar-Mul Specialization

Goal:
Exploit the new backend with better curve-level multiplication strategies.

Status:
Completed on 2026-04-01.

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

Notes:

- Stage 7 added a dedicated `stage7` benchmark profile and a
  `scalar_msm_tuning` benchmark group for explicit scalar-window and
  forced-window MSM comparisons.
- The first Stage 7 full baseline artifact is
  `benchmarks/artifacts/2026-04-01_171130`.
- A synthetic-only G2 small-size retune initially regressed the real
  `prove_full` `B_query_g2` bucket. Direct benchmarking of the actual
  `generated_24_constraints` `pk.B_query_g2` workload showed the best forced
  G2 window at that prover-relevant size was `w = 2`, not the synthetic
  `N = 32` winner `w = 4`.
- The final Stage 7 artifact is `benchmarks/artifacts/2026-04-01_173419`, with
  focused `prove_full` confirmation at `benchmarks/artifacts/2026-04-01_174156`.
- Relative to the Stage 7 baseline `2026-04-01_171130`, Stage 7 improved the
  main `generated_24_constraints` fixture to:
  - `end_to_end`: `30.025 ms -> 29.593 ms`
  - `msm_a_g1`: `3.605 ms -> 2.940 ms`
  - `msm_b_g1`: `2.227 ms -> 2.148 ms`
  - `msm_b_g2`: `4.385 ms -> 4.436 ms`
  - `h_msm`: `7.271 ms -> 7.147 ms`
  - `l_msm`: `3.955 ms -> 3.831 ms`

## Stage 7A: BN254 GLV Scalar Multiplication

Goal:
Add a measured GLV scalar-multiplication path before the Stage 8 prover
re-baseline.

Status:
Completed on 2026-04-01.

Scope:

- implement BN254 GLV decomposition and endomorphism helpers for `G1` and `G2`
- add explicit correctness checks for decomposition and `φ(P) = [λ]P`
- benchmark GLV against the existing scalar path across multiple scalar sizes
- promote GLV into the default path only where the win is measured and the
  semantics remain transparent

Deliverables:

- internal `glv_scalar_mul`, `glv_scalar_decomposition`, and endomorphism
  helpers in `BN254Curve.jl`
- deterministic tests for BN254 GLV decomposition and subgroup-eigenvalue
  checks
- a dedicated `glv_scalar_tuning` benchmark group and `stage7a` profile

Exit Criteria:

- GLV matches the existing scalar path on deterministic subgroup fixtures
- benchmark data identifies whether G1, G2, or both should use GLV by default
- the kept default policy is documented with explicit benchmark evidence

Notes:

- Stage 7A added a GLV helper for both `G1` and `G2`, but only `G1` now uses
  GLV by default.
- The final G1 policy is hybrid:
  - keep w-NAF below `192` bits
  - switch to GLV at `192` bits and above
- `G2` keeps the existing public scalar path for now even though GLV wins on
  full-width subgroup scalars, because `G2Point` does not encode subgroup
  membership and the GLV eigenvalue relation is only guaranteed on the prime
  subgroup.
- The first Stage 7A artifact is
  `benchmarks/artifacts/2026-04-01_182238`, with medians:
  - G1 `bits_128`: default `0.436 ms`, w-NAF `0.440 ms`, GLV `0.612 ms`
  - G1 `bits_192`: default `0.639 ms`, w-NAF `0.654 ms`, GLV `0.648 ms`
  - G1 `bits_254`: default `0.635 ms`, w-NAF `0.967 ms`, GLV `0.640 ms`
  - G2 `bits_128`: default `1.413 ms`, GLV `1.422 ms`
  - G2 `bits_192`: default `1.470 ms`, GLV `1.533 ms`
  - G2 `bits_254`: default `2.249 ms`, GLV `1.435 ms`

## Remaining Gap To Arkworks

The current BN254 primitive gap to arkworks is no longer dominated by the old
`BigInt` representation alone. The biggest remaining causes are now:

- field inversion still escaping to `BigInt invmod`
- final exponentiation still using generic multiplication/squaring in places
  where cyclotomic-specialized paths should be used more aggressively
- extension-field hot paths still being more value-oriented and less in-place
  than the arkworks equivalents
- G2 does have an internal GLV path now, but the current `G2Point` type does
  not encode subgroup membership, so that acceleration is not yet safe as the
  unconditional public default
- MSM still being more generic than arkworks’ specialized variable-base stack

Current priority order for narrowing that gap:

1. Stage 8 `prove_full` re-baseline on the new backend
2. limb-native inversion
3. final-exponentiation specialization
4. more in-place extension-field helpers
5. subgroup-aware G2 typing or another safe way to expose G2 GLV by default
6. deeper MSM specialization

This gap-narrowing track is intentionally placed before Stage 8 so the next
`prove_full` baseline benefits from the higher-value primitive improvements
first.

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

Notes:

- The first Stage 8 baseline artifact is
  `benchmarks/artifacts/2026-04-01_220953`.
- The benchmark fixtures now prove once, assert `verify_full`, and store
  deterministic proof points under `_semantic` before timing begins.
- Relative to the last pre-Stage-8 prover artifact
  (`benchmarks/artifacts/2026-04-01_174156`):
  - `sum_of_products_small` improved from `9.219 ms` to `8.248 ms`
    (`-10.5%`)
  - `generated_24_constraints` stayed effectively flat from `30.045 ms` to
    `30.088 ms` (`+0.1%`)
  - `final_c` improved materially on both fixtures, but the larger fixture saw
    mild regressions in `msm_b_g2`, `h_msm`, and `l_msm`, which offset that
    gain
- The Stage 8 baseline therefore updated the hotspot ranking successfully, but
  it did not yet satisfy the “total `prove_full` time drops materially” exit
  criterion on the primary larger fixture.

### Stage 8A Follow-Through: BN254Fr-Native Prover Scalar Plumbing

Goal:
Eliminate the `BN254Fr -> BigInt` scalar bounce from the prover and setup hot
paths, then measure whether that plumbing change moves the real `prove_full`
baseline.

Status:
Completed on 2026-04-02.

Deliverables:

- direct `BN254Fr` scalar support in `scalar_mul`, `multi_scalar_mul`, and
  fixed-base batch multiplication
- prover/setup/benchmark wiring that uses `BN254Fr` scalars directly instead of
  converting through `BigInt`
- a dedicated scalar-plumbing benchmark comparing `BigInt` and native
  `BN254Fr` prover workloads
- a refreshed `prove_full` baseline and profile tied to the new scalar path

Notes:

- The dedicated scalar-plumbing comparison lives in
  `benchmarks/artifacts/2026-04-01_223814`.
- The full Stage 8A re-baseline artifact is
  `benchmarks/artifacts/2026-04-01_223859`.
- Relative to the first Stage 8 baseline
  (`benchmarks/artifacts/2026-04-01_220953`):
  - `generated_24_constraints` `prove_full` improved from `30.088 ms` to
    `28.873 ms` (`-4.0%`)
  - `msm_b_g2` improved from `4.714 ms` to `4.334 ms` (`-8.1%`)
  - `h_msm` improved from `7.805 ms` to `7.036 ms` (`-9.9%`)
  - `l_msm` improved from `4.210 ms` to `3.856 ms` (`-8.4%`)
  - `final_c` improved from `3.001 ms` to `2.804 ms` (`-6.6%`)
  - the tiny continuity fixture regressed slightly (`8.248 ms -> 8.614 ms`),
    which is acceptable because the larger deterministic fixture is the primary
    optimization target
- The Stage 8A `prove_full` profile no longer contains `canonical_bigint` or
  `limbs_to_bigint` in the main `generated_24_constraints` prover dump, which
  confirms that the hot prover scalar conversions were removed successfully.
- `BigInt` / GMP activity still exists elsewhere in the backend, so the scalar
  plumbing win is not the end of the gap-to-arkworks story.

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

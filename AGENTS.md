# Groth.jl AGENTS.md

## Purpose

Define how agents must operate when modifying Groth.jl, a mathematically sensitive zero-knowledge system implemented as a Julia workspace monorepo.

Groth.jl encodes algebraic and cryptographic constructions in code. Correctness depends on both:
- mathematical validity
- sound software architecture

---

## Primary Objective

Build a **correct, inspectable, and high-performance Groth16 implementation in Julia**, optimized for:

- research and experimentation
- REPL and Pluto workflows
- clarity of algebraic structure
- idiomatic Julia design
- eventual performance optimization (after correctness)

---

## Workspace Model (CRITICAL)

Groth.jl is a **Julia Pkg workspace monorepo**.

- The root directory defines the **canonical workspace environment**
- Each subdirectory (except docs, benchmarks, scripts, CI) is a **separate Julia package**
- Packages include:
  - GrothAlgebra
  - GrothCurves
  - GrothProofs
  - GrothExamples
- The workspace uses the root-environment / root-lockfile model ("Option A").
- The canonical lockfile is the root `Manifest.toml`.
- Do not treat package-local manifests or tool-local manifests as the source of truth unless the user explicitly changes this policy.
- `GrothExamples` is notebook-first (Pluto-focused) and is excluded from `Pkg.test()` aggregation.
- Do not treat `GrothExamples` like a normal package test target unless the user explicitly asks for that behavior to change.

### Canonical Setup

julia --project=. -e 'using Pkg; Pkg.instantiate(workspace=true)'

### Canonical Full Validation

julia --project=. scripts/test_all.jl

This runs the full test suite across all packages.

### Package-Scoped Validation

Allowed only when intentionally scoped:

julia --project=GrothAlgebra -e 'using Pkg; Pkg.test()'

Do not assume per-package testing replaces full workspace validation.

---

## Session Bootstrap (MANDATORY)

Before doing substantial work:

1. Read this file
2. Read workspace `PLANS.md`
3. Read `ROADMAP.md`
4. Read relevant `docs/src/` pages
5. Inspect relevant source files

Do not modify code without context.

---

## Required Workflow

### 1. Understand the Math First

Before modifying code:

- identify the algebraic object (field, polynomial, domain, curve, etc.)
- identify invariants
- identify expected degree / structure
- identify domain assumptions

If you cannot explain the math, do not proceed.

---

### 2. Planning

For non-trivial work:

- create or update an ExecPlan in `execplans/`

Plans must define:

- mathematical invariants
- validation strategy
- expected behavior

Skip only for trivial edits.

---

### 3. Implementation Loop

Always follow:

1. Hypothesis (mathematical expectation)
2. Implementation
3. Verification (invariants)
4. Comparison (if applicable)

Never jump directly to coding.

---

## Mathematical Invariants (CRITICAL)

These must NEVER be violated.

### Polynomial Arithmetic

- deg(A + B) ≤ max(deg(A), deg(B))
- deg(A * B) = deg(A) + deg(B)
- division reduces degree only if exact
- interpolation of n points → degree ≤ n-1

---

### QAP (Groth16 Core)

- t(x) = vanishing polynomial over the evaluation domain
- deg(t) = domain size
- P(x) must be divisible by t(x)
- quotient polynomial degree must match expectation

---

### Domain Rules

- domain size must match or exceed constraint count
- padding must be explicit, never implicit
- FFT domain must not silently truncate
- coset evaluation must preserve domain structure

---

## Common Failure Modes

If something is wrong, suspect:

- degree mismatches
- incorrect vanishing polynomial
- domain size off-by-one
- incorrect coset evaluation
- incorrect interpolation points

If something “almost works”, it is likely a domain or degree issue.

---

## Debugging Rules

When debugging:

1. Check domain size
2. Check polynomial degrees
3. Check interpolation inputs
4. Check vanishing polynomial construction
5. Verify divisibility explicitly

Do not guess. Derive or test.

---

## Validation Requirements

For any change, run:

Pkg.test("GrothAlgebra")
Pkg.test("GrothProofs")

AND:

julia --project=. scripts/test_all.jl

Additionally:

- add targeted tests for new logic
- test small fields (e.g., F7, F11)
- test edge cases
- test known failing seeds

---

## Property-Based Testing

Encouraged for:

- polynomial identities
- field arithmetic
- interpolation correctness
- divisibility properties

---

## Cross-Checking

When uncertain:

- compare with arkworks
- compare with py_ecc
- derive manually in small fields

Never rely on a single implementation.

---

## Benchmarking

If behavior or performance changes:

- julia --project=. benchmarks/run.jl
- julia --project=. benchmarks/plot.jl

Update benchmark documentation and compare before/after.

Do not optimize before correctness is established.

---

## Code Style

- 4-space indentation
- lowercase_with_underscores for functions
- docstrings on exported methods
- prefer multiple dispatch over conditionals
- align with existing APIs

---

## Documentation

Update when behavior changes:

- docs/
- examples/
- README

Documentation must reflect actual behavior.

---

## Definition of Done

A task is complete only if:

- tests pass
- invariants hold
- results match expected mathematics
- validation strategy is satisfied
- benchmarks updated (if relevant)
- documentation updated (if needed)

---

## Anti-Patterns (DO NOT DO)

- silent domain resizing or truncation
- assuming degrees without verifying
- modifying math-heavy code without validation
- fixing tests without understanding the math
- copying reference implementations blindly

---

## Final Rule

Groth.jl is a mathematically grounded Julia codebase: mathematical correctness is non-negotiable, and architecture should make that correctness easier to express, verify, and extend.

If you do not understand why something works mathematically, do not modify the implementation until you do.

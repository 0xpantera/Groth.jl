# Groth.jl

```@meta
CurrentModule = GrothProofs
```

Groth.jl is a Julia research platform for BN254 primitives and Groth16. The
current codebase includes a Montgomery-backed BN254 field layer, extension
tower, G1/G2 arithmetic, optimal ate pairing, and an end-to-end Groth16
pipeline, with benchmarks and profiles tied directly to the active roadmap.

The project has two explicit roles. It remains a learning bridge from
RareSkills / `zk-book` concepts to inspectable Julia code, and it is also a
performance-oriented implementation whose hot paths intentionally depart from
the simplest textbook algorithms. See [From Textbook To Optimized Code](@ref
textbook-to-optimized) for that distinction.

## Choose Your Path

### Learn the concepts

Start with [RareSkills / zk-book ↔ Groth.jl Map](@ref rareskills-map), then use
[Architecture Map](@ref architecture-map) to separate reusable ZK primitives,
curve-specific code, and Groth16-specific protocol logic.

### Run Groth16 end to end

Use [Getting Started](@ref) for workspace setup, then walk through
[Groth16 End-to-End](@ref) to build a circuit, create a proof, and verify it.

### Understand the optimized implementation

[From Textbook To Optimized Code](@ref textbook-to-optimized) explains where the
current code intentionally diverges from the direct textbook shape. For an
arkworks comparison, see [Implementation vs Arkworks](@ref
implementation-vs-arkworks).

### Inspect the packages

The package pages explain ownership and examples:

- [Groth Algebra](@ref)
- [Groth Curves](@ref)
- [Groth Proofs](@ref)

For exported symbols and docstrings, use [API Reference](@ref api-reference).

### Review benchmarks

[Benchmarks](@ref) covers the benchmark workflow and external primitive
comparisons. [Benchmark Snapshots](@ref benchmark-snapshots) keeps the longer
historical prover/setup optimization notes.

## Project Shape

```text
GrothAlgebra  -> finite fields, polynomials, FFTs, group/MSM utilities
GrothCurves   -> BN254 tower fields, G1/G2 arithmetic, pairing engine
GrothProofs   -> R1CS, QAP, Groth16 setup/prove/verify
GrothExamples -> Pluto notebooks and walkthroughs
benchmarks    -> BenchmarkTools profiles, plots, tracked summaries
docs          -> learning maps, implementation notes, benchmark narratives
```

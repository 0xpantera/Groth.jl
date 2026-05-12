# Groth.jl

```@meta
CurrentModule = GrothProofs
```

Welcome to the evolving documentation for Groth.jl. These pages will unify the
algebra, curve, pairing, and Groth16 tooling that lives across the repository.

Groth.jl is a Julia research platform for BN254 primitives and Groth16. The
current codebase includes a Montgomery-backed BN254 field layer, extension
tower, G1/G2 arithmetic, optimal ate pairing, and an end-to-end Groth16
pipeline, with benchmarks and profiles tied directly to the active roadmap.

The project has two explicit roles. It remains a learning bridge from
RareSkills / `zk-book` concepts to inspectable Julia code, and it is also a
performance-oriented implementation whose hot paths intentionally depart from
the simplest textbook algorithms. See [From Textbook To Optimized Code](@ref
textbook-to-optimized) for that distinction.

```@contents
Pages = ["index.md"]
Depth = 2
```

```@index
```

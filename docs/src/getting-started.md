# Getting Started

```@meta
CurrentModule = GrothProofs
DocTestSetup = :(using GrothProofs, GrothAlgebra, GrothCurves)
```

This page documents the canonical workspace-first workflow for Groth.jl.

## 1. Instantiate the workspace

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(workspace=true)'
```

This resolves dependencies for all workspace members from the repository root.

## 2. Run monorepo tests

```bash
julia --project=. scripts/test_all.jl
```

This runs package tests for:

- `GrothAlgebra`
- `GrothCurves`
- `GrothCrypto`
- `GrothProofs`

`GrothExamples` is notebook-first and intentionally excluded from `Pkg.test()`
aggregation.

## 3. Run package-scoped tests (optional)

```bash
julia --project=GrothAlgebra -e 'using Pkg; Pkg.test()'
julia --project=GrothCurves -e 'using Pkg; Pkg.test()'
julia --project=GrothProofs -e 'using Pkg; Pkg.test()'
```

Use these when iterating in one package and you do not need the full repo sweep.

## 4. Open tutorials in Pluto

```bash
julia --project=GrothExamples -e 'using Pkg; Pkg.instantiate()'
julia --project=GrothExamples -e 'using Pluto; Pluto.run()'
```

Open:

- `GrothExamples/src/r1cs_qap_pluto.jl` (AbstractAlgebra derivation)
- `GrothExamples/src/r1cs_qap_groth_pluto.jl` (Groth package-native derivation)

## 5. Build docs and run benchmarks

```bash
julia --project=. docs/make.jl
julia --project=. benchmarks/run.jl
julia --project=. benchmarks/plot.jl
```

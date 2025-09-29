
# Groth.jl Workflow Checklist

Use this guide when working inside the Groth.jl repository. For overall workspace
context (e.g., ark-works/zk-book references), see the parent `AGENTS.md`.

## 1. Review current context
- `docs/ROADMAP.md` — active priorities, completed milestones.
- `docs/PACKAGE_REFERENCE.md` — per-package summaries and implementation notes.
- `docs/Implementation_vs_Arkworks.md` — how our implementation aligns with arkworks.

## 2. Branch & setup
- Create feature branches under Groth.jl (e.g., `feature/…`).
- Ensure local packages are devved:
  ```julia
  julia --project -e 'using Pkg; Pkg.develop("GrothAlgebra"); Pkg.develop("GrothCurves"); Pkg.develop("GrothProofs")'
  ```

## 3. Testing & benchmarks
- Run tests for packages you touch:
  ```julia
  julia --project=GrothAlgebra -e 'using Pkg; Pkg.test()'
  julia --project=GrothProofs  -e 'using Pkg; Pkg.test()'
  ```
- Benchmarks live in `benchmarks/run.jl` / `benchmarks/plot.jl`. Regenerate plots
  when behaviour changes and note results in `benchmarks/README.md`.

## 4. Coding guidelines
- Julia style: 4-space indent, `lowercase_with_underscores` functions, docstrings
  on exported methods, dispatch-friendly signatures matching existing APIs.
- Keep package READMEs/docs in sync when user-visible output changes (e.g.,
  update `GrothExamples/` scripts, `docs/ROADMAP.md`).

## 5. Commits & PRs
- Commit messages: concise imperatives scoped to the area (e.g., `groth16: align coset domain`).
- PR descriptions: list commands run (`Pkg.test`, benchmarks) and link the
  relevant roadmap item.
- Update docs/examples as part of the same PR when behaviour shifts.

Stay focused on Groth.jl; other sibling directories are references only.

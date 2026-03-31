# Benchmarks

The `benchmarks/` project captures runtime and throughput measurements for the
hottest Groth16 paths. This page summarises how to regenerate the data and shows
the latest artefacts.

## Running the Suite

1. Instantiate the root workspace:

   ```bash
   julia --project=. -e 'using Pkg; Pkg.instantiate(workspace=true)'
   ```

2. Run the benchmark harness (artifacts land under `benchmarks/artifacts/<run_id>/`):

   ```bash
   julia --project=. benchmarks/run.jl
   ```

3. Regenerate plots (latest run by default, or pass a run id / JSON file):

   ```bash
   julia --project=. benchmarks/plot.jl
   julia --project=. benchmarks/plot.jl 2026-03-05_144036
   ```

4. Compare two snapshots and flag regressions (20% threshold by default):

   ```bash
   julia --project=. benchmarks/compare.jl 2026-03-05_144036 2026-03-06_091500
   julia --project=. benchmarks/compare.jl 2025-09-29_121914 2026-03-05_144036 10
   ```

5. Profile the `prove_full` path on the deterministic benchmark fixtures:

   ```bash
   julia --project=. benchmarks/profile_prove_full.jl
   julia --project=. benchmarks/profile_prove_full.jl --run-id=2026-03-05_144036 --fixture=sum_of_products_small --repetitions=75
   ```

6. Generate a one-shot markdown report (run -> plot -> compare latest-1 vs latest):

   ```bash
   julia --project=. benchmarks/report.jl
   julia --project=. benchmarks/report.jl --skip-run --threshold=10
   ```

The harness writes a timestamped JSON (raw statistics) and PNG charts covering
MSM, pairing, normalisation, Groth16 end-to-end timings, and `prove_full`
fixture breakdowns. The profiling script writes text profiler dumps under the
same artifact tree, but profiling remains a separate workflow from reproducible
timing baselines.

## Latest Snapshot (2025‑09‑29)

```@example
using JSON
json_path = joinpath(@__DIR__, "assets", "results_2025-09-29_121914.json")
results = JSON.parsefile(json_path)
keys(results)
```

Each entry contains per-benchmark medians, deviations, and configuration
metadata (threading, window sizes, curve parameters). Refer to
`benchmarks/results_2025-09-23_204214_env.md` for the environment capture that
accompanied the latest run.

## Plots

![Pairing throughput](assets/pairing.png)

![Groth16 end-to-end](assets/groth16.png)

![MSM G1 timings](assets/msm_g1.png)

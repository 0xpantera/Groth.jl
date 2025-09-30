# Benchmarks

The `benchmarks/` project captures runtime and throughput measurements for the
hottest Groth16 paths. This page summarises how to regenerate the data and shows
the latest artefacts.

## Running the Suite

1. Activate the benchmarks project:

   ```bash
   julia --project=benchmarks -e 'using Pkg; Pkg.instantiate()'
   ```

2. Run the benchmark harness (JSON results land under `benchmarks/`):

   ```bash
   julia --project=benchmarks benchmarks/run.jl
   ```

3. Regenerate plots from an existing JSON snapshot:

   ```bash
   julia --project=benchmarks benchmarks/plot.jl benchmarks/results_2025-09-29_121914.json
   ```

The harness writes a timestamped JSON (raw statistics) and PNG charts covering
MSM, pairing, normalisation, and Groth16 end-to-end timings.

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

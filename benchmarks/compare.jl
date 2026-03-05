#!/usr/bin/env julia

using JSON
using Printf

const ARTIFACTS_DIR = joinpath(@__DIR__, "artifacts")

is_run_id(s::AbstractString) = occursin(r"^\d{4}-\d{2}-\d{2}_\d{6}$", s)

function usage()
    println("Usage: julia --project=. benchmarks/compare.jl <baseline> <candidate> [threshold_pct]")
    println("Inputs can be either a JSON path or a run id under benchmarks/artifacts/.")
    println("Example: julia --project=. benchmarks/compare.jl 2026-03-05_144036 2026-03-06_091500 20")
end

function resolve_input(input::String)
    direct = isabspath(input) ? input : (isfile(input) ? input : joinpath(@__DIR__, input))
    if isfile(direct)
        return direct
    end

    run_id = basename(normpath(input))
    if is_run_id(run_id)
        candidate = joinpath(ARTIFACTS_DIR, run_id, "results", "benchmark_results.json")
        isfile(candidate) || error("Run id '$run_id' exists but results file is missing at $candidate")
        return candidate
    end

    error("Could not resolve '$input' to a results file or run id.")
end

function load_json(path::String)
    open(path) do io
        return JSON.parse(io)
    end
end

function is_trial_entry(node)
    return node isa AbstractDict && haskey(node, "median_seconds") && node["median_seconds"] isa Number
end

function collect_metrics!(out::Dict{String,Float64}, node, prefix::String = "")
    if is_trial_entry(node)
        out[prefix] = Float64(node["median_seconds"])
        return
    end
    if !(node isa AbstractDict)
        return
    end
    for key in sort!(collect(keys(node)))
        key_str = String(key)
        startswith(key_str, "_") && continue
        child = node[key]
        next_prefix = isempty(prefix) ? key_str : string(prefix, ".", key_str)
        collect_metrics!(out, child, next_prefix)
    end
end

function print_meta(label::String, obj)
    meta = obj isa AbstractDict ? get(obj, "_meta", nothing) : nothing
    if !(meta isa AbstractDict)
        println("[$label] no _meta block found")
        return
    end
    version = get(meta, "julia_version", "unknown")
    cpu = get(meta, "cpu_name", "unknown")
    t = get(meta, "threadpools", Dict{String,Any}())
    td = t isa AbstractDict ? get(t, "default", "?") : "?"
    ti = t isa AbstractDict ? get(t, "interactive", "?") : "?"
    commit = get(meta, "git_commit", "unknown")
    println("[$label] julia=$version threads(default=$td, interactive=$ti) cpu=$cpu commit=$commit")
end

function main()
    if length(ARGS) < 2 || length(ARGS) > 3
        usage()
        exit(1)
    end

    baseline_input = ARGS[1]
    candidate_input = ARGS[2]
    threshold_pct = length(ARGS) == 3 ? parse(Float64, ARGS[3]) : 20.0

    baseline_path = resolve_input(baseline_input)
    candidate_path = resolve_input(candidate_input)
    baseline = load_json(baseline_path)
    candidate = load_json(candidate_path)

    baseline_metrics = Dict{String,Float64}()
    candidate_metrics = Dict{String,Float64}()
    collect_metrics!(baseline_metrics, baseline)
    collect_metrics!(candidate_metrics, candidate)

    shared = intersect(Set(keys(baseline_metrics)), Set(keys(candidate_metrics)))
    missing_in_candidate = setdiff(Set(keys(baseline_metrics)), Set(keys(candidate_metrics)))
    new_in_candidate = setdiff(Set(keys(candidate_metrics)), Set(keys(baseline_metrics)))

    println("Comparing median_seconds")
    println("baseline:  $baseline_path")
    println("candidate: $candidate_path")
    print_meta("baseline", baseline)
    print_meta("candidate", candidate)
    println("shared metrics: $(length(shared))")
    !isempty(missing_in_candidate) && println("missing in candidate: $(length(missing_in_candidate))")
    !isempty(new_in_candidate) && println("new in candidate: $(length(new_in_candidate))")
    println()

    rows = Vector{Tuple{String,Float64,Float64,Float64}}()
    for metric in shared
        b = baseline_metrics[metric]
        c = candidate_metrics[metric]
        delta = (c / b - 1.0) * 100.0
        push!(rows, (metric, b, c, delta))
    end
    sort!(rows, by = x -> abs(x[4]), rev = true)

    println("Top deltas by absolute % change:")
    for (metric, b, c, delta) in Iterators.take(rows, min(15, length(rows)))
        trend = delta > 0 ? "slower" : "faster"
        @printf("  %-45s %9.6f -> %9.6f sec  (%+7.2f%%, %s)\n", metric, b, c, delta, trend)
    end
    println()

    regressions = filter(r -> r[4] >= threshold_pct, rows)
    println("Regression threshold: +$(threshold_pct)%")
    if isempty(regressions)
        println("No regressions above threshold.")
        return
    end

    println("Regressions:")
    for (metric, b, c, delta) in regressions
        @printf("  %-45s %9.6f -> %9.6f sec  (%+7.2f%%)\n", metric, b, c, delta)
    end
end

main()

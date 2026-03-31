#!/usr/bin/env julia

using Dates
using JSON
using StatsPlots
using Printf

const BENCH_DIR = @__DIR__
const ROOT_DIR = normpath(joinpath(BENCH_DIR, ".."))
const ARTIFACTS_DIR = joinpath(BENCH_DIR, "artifacts")
const RUN_ID_RE = r"^\d{4}-\d{2}-\d{2}_\d{6}$"

is_run_id(s::AbstractString) = occursin(RUN_ID_RE, s)

function usage()
    println("Usage: julia --project=. benchmarks/report.jl [--skip-run] [--threshold=<pct>]")
end

function latest_run_ids()
    isdir(ARTIFACTS_DIR) || return String[]
    runs = filter(name -> is_run_id(name) && isdir(joinpath(ARTIFACTS_DIR, name)), readdir(ARTIFACTS_DIR))
    return sort(runs)
end

function run_julia(script::String, args::AbstractVector{<:AbstractString} = String[]; capture::Bool = false)
    jcmd = Base.julia_cmd()
    argv = vcat(collect(jcmd.exec), ["--project=.", script], String.(args))
    cmd = Cmd(argv)
    if capture
        return cd(ROOT_DIR) do
            read(cmd, String)
        end
    end
    cd(ROOT_DIR) do
        run(cmd)
    end
    return nothing
end

function load_results(run_id::String)
    results_path = joinpath(ARTIFACTS_DIR, run_id, "results", "benchmark_results.json")
    isfile(results_path) || error("Results missing for run $run_id at $results_path")
    return JSON.parsefile(results_path)
end

function load_meta(run_id::String)
    obj = load_results(run_id)
    meta = get(obj, "_meta", Dict{String,Any}())
    return meta isa AbstractDict ? Dict{String,Any}(meta) : Dict{String,Any}()
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

format_seconds(x::Real) = @sprintf("%.6f", x)
format_pct(x::Real) = @sprintf("%+.2f%%", x)

function metric_rows(baseline::Dict{String,Float64}, candidate::Dict{String,Float64})
    shared = intersect(Set(keys(baseline)), Set(keys(candidate)))
    rows = Vector{NamedTuple{(:metric, :baseline, :candidate, :delta_pct),Tuple{String,Float64,Float64,Float64}}}()
    for metric in shared
        b = baseline[metric]
        c = candidate[metric]
        delta = (c / b - 1.0) * 100.0
        push!(rows, (metric = metric, baseline = b, candidate = c, delta_pct = delta))
    end
    sort!(rows, by = r -> abs(r.delta_pct), rev = true)
    return rows
end

function meta_value(meta::Dict{String,Any}, key::String; missing::String = "missing")
    return haskey(meta, key) ? string(meta[key]) : missing
end

function threadpool_text(meta::Dict{String,Any})
    t = get(meta, "threadpools", nothing)
    if !(t isa AbstractDict)
        return "_missing (run predates metadata capture)_"
    end
    d = get(t, "default", "?")
    i = get(t, "interactive", "?")
    return "default=$d, interactive=$i"
end

function render_markdown_table(headers::Vector{String}, rows::Vector{Vector{String}})
    out = String[]
    push!(out, "| " * join(headers, " | ") * " |")
    push!(out, "| " * join(fill("---", length(headers)), " | ") * " |")
    for row in rows
        push!(out, "| " * join(row, " | ") * " |")
    end
    return out
end

function rows_for_prefix(all_rows, prefix::String)
    filt = filter(r -> startswith(r.metric, prefix * "."), all_rows)
    sort!(filt, by = r -> r.metric)
    return filt
end

function build_summary_table(rows)
    regressions = count(r -> r.delta_pct > 20.0, rows)
    improvements = count(r -> r.delta_pct < -20.0, rows)
    stable = count(r -> abs(r.delta_pct) <= 20.0, rows)
    return [["Regressions > +20%", string(regressions)],
            ["Improvements < -20%", string(improvements)],
            ["Stable (<= 20% abs delta)", string(stable)],
            ["Total shared metrics", string(length(rows))]]
end

function plot_key_group(rows, title::String, outfile::String)
    isempty(rows) && return
    labels = [replace(split(r.metric, ".")[end], "_" => " ") for r in rows]
    base = [r.baseline * 1e3 for r in rows]
    cand = [r.candidate * 1e3 for r in rows]
    x = collect(1:length(labels))
    bar(x .- 0.18, base, bar_width = 0.35, label = "baseline", legend = :topleft,
        title = title, ylabel = "median time (ms)", xticks = (x, labels), xrotation = 30)
    bar!(x .+ 0.18, cand, bar_width = 0.35, label = "candidate")
    png(outfile)
end

function plot_top_deltas(rows, outfile::String)
    top = filter(r -> abs(r.delta_pct) > 0.0, rows)
    isempty(top) && return
    top = top[1:min(end, 12)]
    labels = [r.metric for r in top]
    deltas = [r.delta_pct for r in top]
    colors = [d > 0 ? :firebrick : :seagreen for d in deltas]
    bar(labels, deltas, legend = false, title = "Top absolute deltas (%)", ylabel = "delta %", xrotation = 50, color = colors)
    hline!([0.0], color = :black, linestyle = :dash)
    png(outfile)
end

function parse_args()
    skip_run = false
    threshold = 20.0
    for arg in ARGS
        if arg == "--skip-run"
            skip_run = true
        elseif startswith(arg, "--threshold=")
            threshold = parse(Float64, split(arg, "=", limit = 2)[2])
        elseif arg in ("-h", "--help")
            usage()
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    return skip_run, threshold
end

function write_report(path::String, baseline::String, candidate::String, threshold::Float64)
    baseline_obj = load_results(baseline)
    candidate_obj = load_results(candidate)
    baseline_metrics = Dict{String,Float64}()
    candidate_metrics = Dict{String,Float64}()
    collect_metrics!(baseline_metrics, baseline_obj)
    collect_metrics!(candidate_metrics, candidate_obj)
    rows = metric_rows(baseline_metrics, candidate_metrics)

    cand_meta = load_meta(candidate)
    baseline_meta = load_meta(baseline)
    missing_meta_note = "_missing (run predates metadata capture)_"
    report_dir = dirname(path)
    mkpath(report_dir)

    groth_rows = rows_for_prefix(rows, "groth16")
    prove_rows = rows_for_prefix(rows, "prove_full")
    pairing_rows = rows_for_prefix(rows, "pairing")
    regressions = filter(r -> r.delta_pct >= threshold, rows)
    improvements = filter(r -> r.delta_pct <= -threshold, rows)
    sort!(regressions, by = r -> r.delta_pct, rev = true)
    sort!(improvements, by = r -> r.delta_pct)

    groth_plot = joinpath(report_dir, "groth16_compare.png")
    pairing_plot = joinpath(report_dir, "pairing_compare.png")
    delta_plot = joinpath(report_dir, "delta_overview.png")
    !isempty(groth_rows) && plot_key_group(groth_rows, "Groth16 Pipeline: Baseline vs Candidate", groth_plot)
    !isempty(pairing_rows) && plot_key_group(pairing_rows[1:min(end, 10)], "Pairing Metrics: Baseline vs Candidate", pairing_plot)
    !isempty(rows) && plot_top_deltas(rows, delta_plot)

    lines = String[]
    push!(lines, "# Benchmark Report")
    push!(lines, "")
    push!(lines, "- Generated (UTC): `$(Dates.now(Dates.UTC))`")
    push!(lines, "- Baseline run: `$baseline`")
    push!(lines, "- Candidate run: `$candidate`")
    push!(lines, "")
    push!(lines, "## Summary")
    push!(lines, "")
    append!(lines, render_markdown_table(["Metric", "Value"], build_summary_table(rows)))
    push!(lines, "")
    push!(lines, "## Run Metadata")
    push!(lines, "")
    meta_rows = Vector{Vector{String}}()
    push!(meta_rows, ["Julia", meta_value(baseline_meta, "julia_version"; missing = missing_meta_note), meta_value(cand_meta, "julia_version"; missing = missing_meta_note)])
    push!(meta_rows, ["Threads", threadpool_text(baseline_meta), threadpool_text(cand_meta)])
    push!(meta_rows, ["CPU", meta_value(baseline_meta, "cpu_name"; missing = missing_meta_note), meta_value(cand_meta, "cpu_name"; missing = missing_meta_note)])
    push!(meta_rows, ["Hostname", meta_value(baseline_meta, "hostname"; missing = missing_meta_note), meta_value(cand_meta, "hostname"; missing = missing_meta_note)])
    push!(meta_rows, ["Git commit", meta_value(baseline_meta, "git_commit"; missing = missing_meta_note), meta_value(cand_meta, "git_commit"; missing = missing_meta_note)])
    push!(meta_rows, ["Git branch", meta_value(baseline_meta, "git_branch"; missing = missing_meta_note), meta_value(cand_meta, "git_branch"; missing = missing_meta_note)])
    append!(lines, render_markdown_table(["Field", "Baseline", "Candidate"], meta_rows))
    push!(lines, "")
    push!(lines, "## Comparison Plots")
    push!(lines, "")
    push!(lines, "![Groth16 Comparison](./groth16_compare.png)")
    push!(lines, "")
    push!(lines, "![Pairing Comparison](./pairing_compare.png)")
    push!(lines, "")
    if isfile(delta_plot)
        push!(lines, "![Top Delta Overview](./delta_overview.png)")
        push!(lines, "")
    end
    push!(lines, "## Groth16 Table")
    push!(lines, "")
    groth_table = [[r.metric, format_seconds(r.baseline), format_seconds(r.candidate), format_pct(r.delta_pct)] for r in groth_rows]
    append!(lines, render_markdown_table(["Metric", "Baseline (s)", "Candidate (s)", "Delta"], groth_table))
    push!(lines, "")
    if !isempty(prove_rows)
        push!(lines, "## prove_full Breakdown")
        push!(lines, "")
        prove_table = [[r.metric, format_seconds(r.baseline), format_seconds(r.candidate), format_pct(r.delta_pct)] for r in prove_rows]
        append!(lines, render_markdown_table(["Metric", "Baseline (s)", "Candidate (s)", "Delta"], prove_table))
        push!(lines, "")
    end
    push!(lines, "## Regressions (>= +$(threshold)%)")
    push!(lines, "")
    reg_rows = isempty(regressions) ? [["None", "-", "-", "-"]] :
        [[r.metric, format_seconds(r.baseline), format_seconds(r.candidate), format_pct(r.delta_pct)] for r in regressions]
    append!(lines, render_markdown_table(["Metric", "Baseline (s)", "Candidate (s)", "Delta"], reg_rows))
    push!(lines, "")
    push!(lines, "## Improvements (<= -$(threshold)%)")
    push!(lines, "")
    imp_rows = isempty(improvements) ? [["None", "-", "-", "-"]] :
        [[r.metric, format_seconds(r.baseline), format_seconds(r.candidate), format_pct(r.delta_pct)] for r in improvements[1:min(end, 20)]]
    append!(lines, render_markdown_table(["Metric", "Baseline (s)", "Candidate (s)", "Delta"], imp_rows))
    push!(lines, "")
    push!(lines, "## Artifacts")
    push!(lines, "")
    push!(lines, "- Results: `benchmarks/artifacts/$candidate/results/benchmark_results.json`")
    push!(lines, "- Plots: `benchmarks/artifacts/$candidate/plots/`")
    push!(lines, "- Report assets: `benchmarks/artifacts/$candidate/report/`")
    open(path, "w") do io
        write(io, join(lines, "\n"))
        write(io, "\n")
    end
end

function main()
    skip_run, threshold = parse_args()

    if !skip_run
        println("Running benchmarks/run.jl ...")
        run_julia("benchmarks/run.jl")
    end

    runs = latest_run_ids()
    isempty(runs) && error("No artifact runs found in benchmarks/artifacts.")

    candidate = runs[end]
    baseline = if length(runs) >= 2
        runs[end - 1]
    else
        error("Need at least two runs to compare; only found $candidate.")
    end

    println("Generating plots for $candidate ...")
    run_julia("benchmarks/plot.jl", [candidate])

    println("Building report for $baseline -> $candidate ...")

    report_dir = joinpath(ARTIFACTS_DIR, candidate, "report")
    mkpath(report_dir)
    report_path = joinpath(report_dir, "benchmark_report.md")
    write_report(report_path, baseline, candidate, threshold)
    println("Report written to $report_path")
end

main()

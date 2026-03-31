#!/usr/bin/env julia

using JSON
using StatsPlots
using Printf
using Dates

const ARTIFACTS_DIR = joinpath(@__DIR__, "artifacts")
const PROVE_FULL_PHASE_ORDER = [
    "end_to_end",
    "witness_to_scalars",
    "msm_a_g1",
    "msm_b_g1",
    "msm_b_g2",
    "compute_h_total",
    "h_poly_assembly",
    "h_dense_quotient",
    "h_coset_fft",
    "h_parity_assert",
    "h_msm",
    "l_msm",
    "final_c",
]

is_run_id(s::AbstractString) = occursin(r"^\d{4}-\d{2}-\d{2}_\d{6}$", s)

function latest_artifact_run_id()
    isdir(ARTIFACTS_DIR) || return nothing
    runs = filter(name -> is_run_id(name) && isdir(joinpath(ARTIFACTS_DIR, name)), readdir(ARTIFACTS_DIR))
    isempty(runs) && return nothing
    return sort(runs)[end]
end

function latest_legacy_results_path()
    files = filter(f -> occursin(r"^results_\d{4}-\d{2}-\d{2}_\d{6}\.json$", f), readdir(@__DIR__))
    isempty(files) && return nothing
    return joinpath(@__DIR__, sort(files)[end])
end

function resolve_input(arg::Union{Nothing,String})
    if arg === nothing
        run_id = latest_artifact_run_id()
        if run_id !== nothing
            return joinpath(ARTIFACTS_DIR, run_id, "results", "benchmark_results.json"), run_id
        end
        legacy = latest_legacy_results_path()
        legacy === nothing && error("No benchmark results found in artifacts/ or legacy benchmarks/results_*.json files.")
        return legacy, nothing
    end

    candidate = isabspath(arg) ? arg : (isfile(arg) ? arg : joinpath(@__DIR__, arg))
    if isfile(candidate)
        return candidate, nothing
    end

    run_id = basename(normpath(arg))
    if is_run_id(run_id)
        run_path = joinpath(ARTIFACTS_DIR, run_id, "results", "benchmark_results.json")
        isfile(run_path) || error("Run id '$run_id' exists but results file is missing at $run_path")
        return run_path, run_id
    end

    error("Could not resolve input '$arg' to a results file or run id.")
end

function materialize_results(results_path::String, run_id_hint::Union{Nothing,String})
    if run_id_hint !== nothing
        run_id = run_id_hint
        canonical_results = joinpath(ARTIFACTS_DIR, run_id, "results", "benchmark_results.json")
        return canonical_results, joinpath(ARTIFACTS_DIR, run_id), false
    end

    run_id = Dates.format(Dates.now(Dates.UTC), "yyyy-mm-dd_HHMMSS")
    artifact_run_dir = joinpath(ARTIFACTS_DIR, run_id)
    results_dir = joinpath(artifact_run_dir, "results")
    plots_dir = joinpath(artifact_run_dir, "plots")
    meta_dir = joinpath(artifact_run_dir, "meta")
    mkpath(results_dir)
    mkpath(plots_dir)
    mkpath(meta_dir)
    canonical_results = joinpath(results_dir, "benchmark_results.json")
    if abspath(results_path) != abspath(canonical_results)
        cp(results_path, canonical_results; force = true)
    end
    return canonical_results, artifact_run_dir, true
end

function load_results(path::String)
    open(path) do io
        return JSON.parse(io)
    end
end

median_ms(entry) = entry["median_seconds"] * 1e3

function plot_group(results::AbstractDict, out_dir::String, key::Symbol, title::String, labels::Vector{String}, outfile::String)
    group = get(results, String(key), nothing)
    group === nothing && return
    Ns = sort(parse.(Int, collect(keys(group))))
    data = [median_ms(group[string(N)][labels[j]]) for j in eachindex(labels), N in Ns]
    groupedbar(string.(Ns), permutedims(data), bar_position=:dodge,
        legend=:topleft, title=title, xlabel="N", ylabel="median time (ms)")
    png(joinpath(out_dir, outfile))
end

function plot_categorical_group(results::AbstractDict, out_dir::String, key::Symbol, categories::Vector{String}, labels::Vector{String}, title::String, outfile::String)
    group = get(results, String(key), nothing)
    group === nothing && return
    present = [category for category in categories if haskey(group, category)]
    isempty(present) && return
    data = [median_ms(group[category][labels[j]]) for j in eachindex(labels), category in present]
    groupedbar(present, permutedims(data), bar_position=:dodge,
        legend=:topleft, title=title, xlabel="", ylabel="median time (ms)")
    png(joinpath(out_dir, outfile))
end

function plot_single_ops(results::AbstractDict, out_dir::String, key::Symbol, order::Vector{String}, title::String, outfile::String)
    group = get(results, String(key), nothing)
    group === nothing && return
    labels = [label for label in order if haskey(group, label)]
    values = [median_ms(group[label]) for label in labels]
    pretty = [replace(label, "_" => " ") for label in labels]
    bar(pretty, values; legend=false, title=title, ylabel="median time (ms)", xticks=:auto, xrotation=30)
    png(joinpath(out_dir, outfile))
end

function plot_prove_full_breakdowns(results::AbstractDict, out_dir::String)
    group = get(results, "prove_full", nothing)
    group === nothing && return
    for fixture_name in sort!(collect(keys(group)))
        startswith(fixture_name, "_") && continue
        fixture = group[fixture_name]
        fixture isa AbstractDict || continue
        labels = [label for label in PROVE_FULL_PHASE_ORDER if haskey(fixture, label)]
        isempty(labels) && continue
        values = [median_ms(fixture[label]) for label in labels]
        pretty = [replace(label, "_" => " ") for label in labels]
        bar(pretty, values; legend=false, title="prove_full breakdown: $(fixture_name)", ylabel="median time (ms)", xrotation=45)
        png(joinpath(out_dir, "prove_full_$(fixture_name).png"))
    end
end

function main()
    resolved_input, run_id_hint = resolve_input(length(ARGS) > 0 ? ARGS[1] : nothing)
    canonical_results, artifact_run_dir, copied = materialize_results(resolved_input, run_id_hint)
    out_dir = joinpath(artifact_run_dir, "plots")
    mkpath(out_dir)

    println("Loading ", canonical_results)
    res = load_results(canonical_results)

    # Fixed-base / MSM / normalization
    plot_group(res, out_dir, :fixed_g1, "Fixed-base G1 (median)", ["naive", "batch"], "fixed_g1.png")
    plot_group(res, out_dir, :fixed_g2, "Fixed-base G2 (median)", ["naive", "batch"], "fixed_g2.png")
    plot_group(res, out_dir, :msm_g1, "MSM G1 (median)", ["naive", "msm"], "msm_g1.png")
    plot_group(res, out_dir, :msm_g2, "MSM G2 (median)", ["naive", "msm"], "msm_g2.png")
    plot_group(res, out_dir, :norm_g1, "Batch norm G1 (median)", ["each", "batch"], "norm_g1.png")
    plot_group(res, out_dir, :norm_g2, "Batch norm G2 (median)", ["each", "batch"], "norm_g2.png")
    plot_categorical_group(res, out_dir, :py_ecc_scalar, ["g1", "g2"], ["groth_jl", "py_ecc"],
        "BN254 scalar multiplication: Groth.jl vs py_ecc", "py_ecc_scalar.png")
    plot_group(res, out_dir, :py_ecc_naive_accum_g1, "Naive variable-base accumulation G1: Groth.jl vs py_ecc", ["groth_jl", "py_ecc"], "py_ecc_naive_accum_g1.png")
    plot_group(res, out_dir, :py_ecc_naive_accum_g2, "Naive variable-base accumulation G2: Groth.jl vs py_ecc", ["groth_jl", "py_ecc"], "py_ecc_naive_accum_g2.png")
    plot_categorical_group(res, out_dir, :py_ecc_pairing, ["single"], ["groth_jl", "py_ecc"],
        "BN254 pairing: Groth.jl vs py_ecc", "py_ecc_pairing.png")

    # Pairing summaries
    plot_group(res, out_dir, :pairing, "Pairing sequential vs batch", ["sequential", "batch"], "pairing.png")
    plot_single_ops(res, out_dir, :pairing_single,
        ["pairing", "miller_loop", "final_exponentiation"],
        "Pairing micro-operations", "pairing_ops.png")

    # Groth16 end-to-end
    plot_single_ops(res, out_dir, :groth16,
        ["r1cs_to_qap", "setup", "prove", "verify_full", "prepare_vk", "prepare_inputs", "verify_prepared"],
        "Groth16 pipeline", "groth16.png")
    plot_prove_full_breakdowns(res, out_dir)

    copied && println("Copied source results into artifact run directory.")
    println("Saved plots to ", out_dir)
end

main()

#!/usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__)

using JSON
using StatsPlots
using Printf
using Dates

function latest_results()
    files = filter(f->occursin(r"^results_\d{4}-\d{2}-\d{2}_\d{6}\.json$", f), readdir(@__DIR__))
    isempty(files) && error("No results_*.json files found in benchmarks/")
    sort(files)[end]
end

function load_results(path::String)
    open(joinpath(@__DIR__, path)) do io
        return JSON.parse(io)
    end
end

median_ms(entry) = entry["median_seconds"] * 1e3

function plot_group(results::Dict, key::Symbol, title::String, labels::Vector{String}, outfile::String)
    group = get(results, String(key), nothing)
    group === nothing && return
    Ns = sort(parse.(Int, collect(keys(group))))
    data = [median_ms(group[string(N)][labels[j]]) for j in eachindex(labels), N in Ns]
    groupedbar(string.(Ns), permutedims(data), bar_position=:dodge,
        legend=:topleft, title=title, xlabel="N", ylabel="median time (ms)", label=labels)
    png(joinpath(@__DIR__, outfile))
end

function plot_single_ops(results::Dict, key::Symbol, order::Vector{String}, title::String, outfile::String)
    group = get(results, String(key), nothing)
    group === nothing && return
    labels = [label for label in order if haskey(group, label)]
    values = [median_ms(group[label]) for label in labels]
    pretty = [replace(label, "_" => " ") for label in labels]
    bar(pretty, values; legend=false, title=title, ylabel="median time (ms)", xticks=:auto, xrotation=30)
    png(joinpath(@__DIR__, outfile))
end

function main()
    file = length(ARGS) > 0 ? ARGS[1] : latest_results()
    println("Loading ", file)
    res = load_results(file)

    # Fixed-base / MSM / normalization
    plot_group(res, :fixed_g1, "Fixed-base G1 (median)", ["naive","batch"], "fixed_g1.png")
    plot_group(res, :fixed_g2, "Fixed-base G2 (median)", ["naive","batch"], "fixed_g2.png")
    plot_group(res, :msm_g1, "MSM G1 (median)", ["naive","msm"], "msm_g1.png")
    plot_group(res, :msm_g2, "MSM G2 (median)", ["naive","msm"], "msm_g2.png")
    plot_group(res, :norm_g1, "Batch norm G1 (median)", ["each","batch"], "norm_g1.png")
    plot_group(res, :norm_g2, "Batch norm G2 (median)", ["each","batch"], "norm_g2.png")

    # Pairing summaries
    plot_group(res, :pairing, "Pairing sequential vs batch", ["sequential","batch"], "pairing.png")
    plot_single_ops(res, :pairing_single,
        ["pairing", "miller_loop", "final_exponentiation"],
        "Pairing micro-operations", "pairing_ops.png")

    # Groth16 end-to-end
    plot_single_ops(res, :groth16,
        ["r1cs_to_qap", "setup", "prove", "verify_full", "prepare_vk", "prepare_inputs", "verify_prepared"],
        "Groth16 pipeline", "groth16.png")

    println("Saved plots to benchmarks/*.png")
end

main()

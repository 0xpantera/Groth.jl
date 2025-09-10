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

function parse_ns(s::String)
    # TrialEstimate string, take the first number and unit
    # e.g., "TrialEstimate(939.764 ms)" -> 939.764, unit ms
    m = match(r"TrialEstimate\(([^\)]+)\)", s)
    v = m === nothing ? NaN : parse(Float64, split(m.captures[1])[1])
    u = m === nothing ? "" : split(m.captures[1])[2]
    # Convert to milliseconds
    factor = u == "ns" ? 1e-6 : u == "μs" ? 1e-3 : u == "ms" ? 1.0 : u == "s" ? 1e3 : 1.0
    return v * factor
end

function plot_group(results::Dict, key::Symbol, title::String, labels::Vector{String}, outfile::String)
    keystr = String(key)
    Ns = sort(parse.(Int, collect(keys(results[keystr]))))
    data = [parse_ns(results[keystr][string(N)][labels[j]]["med"]) for j in eachindex(labels), N in Ns]
    groupedbar(string.(Ns), permutedims(data), bar_position=:dodge,
        legend=:topleft, title=title, xlabel="N", ylabel="median time (ms)", label=labels)
    png(joinpath(@__DIR__, outfile))
end

function main()
    file = length(ARGS) > 0 ? ARGS[1] : latest_results()
    println("Loading ", file)
    res = load_results(file)

    # Fixed-base G1/G2
    plot_group(res, :fixed_g1, "Fixed-base G1 (median)", ["naive","batch"], "fixed_g1.png")
    plot_group(res, :fixed_g2, "Fixed-base G2 (median)", ["naive","batch"], "fixed_g2.png")
    plot_group(res, :msm_g1, "MSM G1 (median)", ["naive","msm"], "msm_g1.png")
    plot_group(res, :msm_g2, "MSM G2 (median)", ["naive","msm"], "msm_g2.png")
    plot_group(res, :norm_g1, "Batch norm G1 (median)", ["each","batch"], "norm_g1.png")
    plot_group(res, :norm_g2, "Batch norm G2 (median)", ["each","batch"], "norm_g2.png")
    println("Saved plots to benchmarks/*.png")
end

main()

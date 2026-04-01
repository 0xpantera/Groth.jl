#!/usr/bin/env julia

using BenchmarkTools
using Random
using JSON
using Dates
using Sockets
using Statistics
using GrothAlgebra
using GrothCurves
using GrothProofs

include("prove_full_common.jl")

const Fr = BN254Fr
const r = BN254_ORDER_R
const ENGINE = BN254_ENGINE
const WORKSPACE_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PY_ECC_ROOT = joinpath(WORKSPACE_ROOT, "py_ecc")
const PY_ECC_SCRIPT = joinpath(@__DIR__, "py_ecc_bn254_bench.py")
const PY_ECC_ACCUM_SIZES = (32, 128)
const BENCHMARK_GROUP_ORDER = [
    :bn254_primitives,
    :fixed_base,
    :variable_msm,
    :batch_norm,
    :pairing_throughput,
    :pairing_micro,
    :py_ecc_primitives,
    :groth16,
    :prove_full,
]
const BENCHMARK_PROFILES = Dict(
    "full" => copy(BENCHMARK_GROUP_ORDER),
    "quick" => [:bn254_primitives, :pairing_micro],
    "stage0" => [:bn254_primitives, :pairing_micro],
    "primitives" => [:bn254_primitives, :pairing_micro, :py_ecc_primitives],
)

# -----------------------------------------------------------------------------
# Sampling helpers
# -----------------------------------------------------------------------------

function rand_scalar(R::AbstractRNG)
    # Sample in [0, r)
    return BigInt(rand(R, 0:r-1))
end

function rand_scalar_nonzero(R::AbstractRNG)
    # Sample in [1, r-1] to avoid the point at infinity in pairings
    return BigInt(rand(R, 1:r-1))
end

function genscalars(R::AbstractRNG, N::Int)
    return [rand_scalar(R) for _ in 1:N]
end

function genpoints_fixed(g, scalars::Vector{BigInt})
    # Naive scalar multiplication loop
    return [scalar_mul(g, s) for s in scalars]
end

det_py_ecc_scalar(tag::Integer) = mod(BigInt(tag)^3 + 19 * BigInt(tag) + 7, r - 1) + 1

function py_ecc_accum_inputs(gen, base_offset::Int, scalar_offset::Int, size::Int)
    bases = [scalar_mul(gen, det_py_ecc_scalar(base_offset + i)) for i in 1:size]
    scalars = [det_py_ecc_scalar(scalar_offset + i) for i in 1:size]
    return bases, scalars
end

function serialize_bn254(x::BN254Fq)
    return string(x.value)
end

function serialize_bn254(x::BN254Fr)
    return string(x.value)
end

function serialize_bn254(x::Fp2Element)
    return Any[serialize_bn254(x[1]), serialize_bn254(x[2])]
end

function serialize_bn254(x::Fp6Element)
    return Any[serialize_bn254(x[1]), serialize_bn254(x[2]), serialize_bn254(x[3])]
end

function serialize_bn254(x::Fp12Element)
    return Any[serialize_bn254(x[1]), serialize_bn254(x[2])]
end

function serialize_affine(p::G1Point)
    x, y = to_affine(p)
    return Any[serialize_bn254(x), serialize_bn254(y)]
end

function serialize_affine(p::G2Point)
    x, y = to_affine(p)
    return Any[serialize_bn254(x), serialize_bn254(y)]
end

# -----------------------------------------------------------------------------
# Trial helpers for consistent JSON output
# -----------------------------------------------------------------------------

function estimate_to_seconds(est)
    # Work with the pretty string so units match printed output
    s = replace(string(est), "\u03bc" => "u")  # normalize micro symbol
    m = match(r"TrialEstimate\(([0-9\.eE+-]+) ([a-z]+)\)", s)
    m === nothing && error("Unable to parse TrialEstimate string: $s")
    value = parse(Float64, m.captures[1])
    unit = m.captures[2]
    factor = unit == "s" ? 1.0 :
             unit == "ms" ? 1e-3 :
             unit == "us" ? 1e-6 :
             unit == "ns" ? 1e-9 :
             error("Unsupported time unit $unit in TrialEstimate")
    return value * factor
end

function trial_summary(tr::BenchmarkTools.Trial)
    min_est = minimum(tr)
    med_est = median(tr)
    Dict(
        "min_pretty" => string(min_est),
        "median_pretty" => string(med_est),
        "min_seconds" => estimate_to_seconds(min_est),
        "median_seconds" => estimate_to_seconds(med_est),
        "memory_bytes" => memory(tr),
    )
end

function record_result!(results::Dict{Symbol,Any}, group::Symbol, nkey::String, label::String, tr::BenchmarkTools.Trial)
    group_dict = get!(results, group, Dict{String,Any}())
    entry = get!(group_dict, nkey, Dict{String,Any}())
    entry[label] = trial_summary(tr)
end

function record_simple!(results::Dict{Symbol,Any}, group::Symbol, label::String, tr::BenchmarkTools.Trial)
    group_dict = get!(results, group, Dict{String,Any}())
    group_dict[label] = trial_summary(tr)
end

function record_external_result!(results::Dict{Symbol,Any}, group::Symbol, nkey::String, label::String, summary::AbstractDict)
    group_dict = get!(results, group, Dict{String,Any}())
    entry = get!(group_dict, nkey, Dict{String,Any}())
    entry[label] = Dict{String,Any}(pairs(summary))
end

function record_fixture_trial!(results::Dict{Symbol,Any}, group::Symbol, fixture_name::String, label::String, tr::BenchmarkTools.Trial)
    group_dict = get!(results, group, Dict{String,Any}())
    fixture_dict = get!(group_dict, fixture_name, Dict{String,Any}())
    fixture_dict[label] = trial_summary(tr)
end

function record_fixture_meta!(results::Dict{Symbol,Any}, group::Symbol, fixture_name::String, meta::Dict{String,Any})
    group_dict = get!(results, group, Dict{String,Any}())
    fixture_dict = get!(group_dict, fixture_name, Dict{String,Any}())
    fixture_dict["_fixture"] = meta
end

function record_semantic!(results::Dict{Symbol,Any}, group::Symbol, label::String, value)
    semantic_root = get!(results, :_semantic, Dict{String,Any}())
    group_dict = get!(semantic_root, String(group), Dict{String,Any}())
    group_dict[label] = value
end

function print_stats(label, tr::BenchmarkTools.Trial)
    tmin = minimum(tr)
    tmed = median(tr)
    tmean = mean(tr)
    bytes = memory(tr)
    println(rpad(label, 40), " min=$(tmin) med=$(tmed) mean=$(tmean) mem=$(bytes)B")
end

function print_external_stats(label::AbstractString, summary::AbstractDict)
    println(rpad(label, 40), " min=$(summary["min_pretty"]) med=$(summary["median_pretty"]) mem=$(summary["memory_bytes"])")
end

function try_readchomp(cmd::Cmd)
    try
        return strip(readchomp(cmd))
    catch
        return "unknown"
    end
end

function run_py_ecc_bench(meta::Dict{String,Any})
    python = Sys.which("python3")
    if python === nothing
        note = Dict{String,Any}(
            "status" => "skipped",
            "reason" => "python3 not found on PATH",
        )
        meta["py_ecc"] = note
        println("\nSkipping py_ecc primitive comparison: ", note["reason"])
        return nothing
    end

    if !isdir(PY_ECC_ROOT)
        note = Dict{String,Any}(
            "status" => "skipped",
            "reason" => "workspace sibling py_ecc checkout not found",
            "expected_root" => PY_ECC_ROOT,
        )
        meta["py_ecc"] = note
        println("\nSkipping py_ecc primitive comparison: ", note["reason"])
        return nothing
    end

    output = read(Cmd([python, PY_ECC_SCRIPT]), String)
    payload = JSON.parse(output)
    py_meta = payload["meta"]
    py_meta["status"] = "measured"
    py_meta["python_executable"] = python
    meta["py_ecc"] = py_meta
    return payload
end

function run_metadata()
    return Dict{String,Any}(
        "timestamp_utc" => string(Dates.now(Dates.UTC)),
        "julia_version" => string(VERSION),
        "threadpools" => Dict(
            "default" => Threads.nthreads(:default),
            "interactive" => Threads.nthreads(:interactive),
        ),
        "cpu_name" => Sys.CPU_NAME,
        "kernel" => string(Sys.KERNEL),
        "machine" => Sys.MACHINE,
        "hostname" => gethostname(),
        "git_commit" => try_readchomp(`git rev-parse HEAD`),
        "git_branch" => try_readchomp(`git rev-parse --abbrev-ref HEAD`),
    )
end

function print_available_profiles()
    println("Available benchmark profiles:")
    for name in sort!(collect(keys(BENCHMARK_PROFILES)))
        groups = join(string.(BENCHMARK_PROFILES[name]), ", ")
        println("  ", rpad(name, 10), " -> ", groups)
    end
    println("\nAvailable benchmark groups:")
    for group in BENCHMARK_GROUP_ORDER
        println("  ", group)
    end
end

function parse_group_list(spec::AbstractString)
    isempty(strip(spec)) && error("Benchmark group list cannot be empty")
    requested = Symbol[]
    seen = Set{Symbol}()
    valid = Set(BENCHMARK_GROUP_ORDER)
    for token in split(spec, ",")
        name = Symbol(strip(token))
        name in valid || error("Unknown benchmark group '$token'. Use --list-profiles to see valid groups.")
        name in seen && continue
        push!(requested, name)
        push!(seen, name)
    end
    return requested
end

function parse_cli_args(args::Vector{String})
    profile = "full"
    selected_groups = nothing
    list_profiles = false

    for arg in args
        if arg == "--list-profiles"
            list_profiles = true
        elseif startswith(arg, "--profile=")
            profile = split(arg, "=", limit = 2)[2]
        elseif startswith(arg, "--groups=")
            selected_groups = parse_group_list(split(arg, "=", limit = 2)[2])
        else
            error("Unrecognized argument '$arg'. Supported flags: --profile=<name>, --groups=a,b,c, --list-profiles")
        end
    end

    if list_profiles
        return (list_profiles = true, profile = profile, selected_groups = Symbol[])
    end

    if selected_groups === nothing
        haskey(BENCHMARK_PROFILES, profile) || error("Unknown benchmark profile '$profile'. Use --list-profiles to inspect available profiles.")
        selected_groups = copy(BENCHMARK_PROFILES[profile])
    else
        profile = "custom"
    end

    return (list_profiles = false, profile = profile, selected_groups = selected_groups)
end

# -----------------------------------------------------------------------------
# Benchmark families
# -----------------------------------------------------------------------------

function bench_fixed_base(results)
    println("\n== Fixed-base precompute: table build and batch_mul (G1) ==")
    g1 = g1_generator()
    rng = MersenneTwister(42)
    for N in (32, 128, 512)
        scalars = genscalars(rng, N)
        println("N = ", N)
        # Warmup
        _ = build_fixed_table(g1)
        _ = genpoints_fixed(g1, scalars)
        # Separate table build vs batch_mul
        tr_build = @benchmark build_fixed_table($g1) seconds = 1 samples = 10
        tab = build_fixed_table(g1)
        tr_naive = @benchmark genpoints_fixed($g1, $scalars) seconds = 1 samples = 10
        tr_batch = @benchmark batch_mul($tab, $scalars) seconds = 1 samples = 10
        print_stats("G1 build", tr_build)
        print_stats("G1 naive", tr_naive)
        print_stats("G1 batch", tr_batch)
        record_result!(results, :fixed_g1, string(N), "build", tr_build)
        record_result!(results, :fixed_g1, string(N), "naive", tr_naive)
        record_result!(results, :fixed_g1, string(N), "batch", tr_batch)
    end

    println("\n== Fixed-base precompute: table build and batch_mul (G2) ==")
    g2 = g2_generator()
    rng = MersenneTwister(123)
    for N in (32, 128, 512)
        scalars = genscalars(rng, N)
        println("N = ", N)
        _ = build_fixed_table(g2)
        _ = genpoints_fixed(g2, scalars)
        tr_build = @benchmark build_fixed_table($g2) seconds = 1 samples = 10
        tab = build_fixed_table(g2)
        tr_naive = @benchmark genpoints_fixed($g2, $scalars) seconds = 1 samples = 10
        tr_batch = @benchmark batch_mul($tab, $scalars) seconds = 1 samples = 10
        print_stats("G2 build", tr_build)
        print_stats("G2 naive", tr_naive)
        print_stats("G2 batch", tr_batch)
        record_result!(results, :fixed_g2, string(N), "build", tr_build)
        record_result!(results, :fixed_g2, string(N), "naive", tr_naive)
        record_result!(results, :fixed_g2, string(N), "batch", tr_batch)
    end
end

function naive_msm(bases, scalars)
    acc = zero(bases[1])
    @inbounds for i in eachindex(bases)
        s = scalars[i]
        iszero(s) && continue
        acc += scalar_mul(bases[i], s)
    end
    return acc
end

function gen_random_bases_g1(rng::AbstractRNG, N::Int)
    g = g1_generator()
    scal = genscalars(rng, N)
    return [scalar_mul(g, s) for s in scal]
end

function gen_random_bases_g2(rng::AbstractRNG, N::Int)
    g = g2_generator()
    scal = genscalars(rng, N)
    return [scalar_mul(g, s) for s in scal]
end

function bench_variable_msm(results)
    println("\n== Variable-base MSM: multi_scalar_mul vs naive (G1) ==")
    rng = MersenneTwister(7)
    for N in (32, 128, 512)
        bases = gen_random_bases_g1(rng, N)
        scalars = genscalars(rng, N)
        println("N = ", N)
        _ = naive_msm(bases, scalars)
        _ = GrothAlgebra.multi_scalar_mul(bases, scalars)
        tr_naive = @benchmark naive_msm($bases, $scalars) seconds = 1 samples = 10
        tr_msm = @benchmark GrothAlgebra.multi_scalar_mul($bases, $scalars) seconds = 1 samples = 10
        print_stats("G1 naive", tr_naive)
        print_stats("G1 MSM", tr_msm)
        record_result!(results, :msm_g1, string(N), "naive", tr_naive)
        record_result!(results, :msm_g1, string(N), "msm", tr_msm)
    end

    println("\n== Variable-base MSM: multi_scalar_mul vs naive (G2) ==")
    rng = MersenneTwister(9)
    for N in (32, 128, 512)
        bases = gen_random_bases_g2(rng, N)
        scalars = genscalars(rng, N)
        println("N = ", N)
        _ = naive_msm(bases, scalars)
        _ = GrothAlgebra.multi_scalar_mul(bases, scalars)
        tr_naive = @benchmark naive_msm($bases, $scalars) seconds = 1 samples = 10
        tr_msm = @benchmark GrothAlgebra.multi_scalar_mul($bases, $scalars) seconds = 1 samples = 10
        print_stats("G2 naive", tr_naive)
        print_stats("G2 MSM", tr_msm)
        record_result!(results, :msm_g2, string(N), "naive", tr_naive)
        record_result!(results, :msm_g2, string(N), "msm", tr_msm)
    end
end

function bench_batch_norm(results)
    println("\n== Batch normalization: batch_to_affine! vs per-point to_affine (G1) ==")
    g = g1_generator()
    rng = MersenneTwister(11)
    for N in (32, 128, 512)
        scalars = genscalars(rng, N)
        proj_pts = [scalar_mul(g, s) for s in scalars]
        println("N = ", N)
        _ = GrothCurves.batch_to_affine!(copy(proj_pts))
        _ = begin
            tmp = copy(proj_pts)
            for i in eachindex(tmp)
                x, y = to_affine(tmp[i])
                tmp[i] = G1Point(x, y, one(BN254Fq))
            end
        end
        tr_batch = @benchmark begin
            tmp = copy($proj_pts)
            GrothCurves.batch_to_affine!(tmp)
        end seconds = 1 samples = 10
        tr_each = @benchmark begin
            tmp = copy($proj_pts)
            for i in eachindex(tmp)
                x, y = to_affine(tmp[i])
                tmp[i] = G1Point(x, y, one(BN254Fq))
            end
        end seconds = 1 samples = 10
        print_stats("G1 batch_norm", tr_batch)
        print_stats("G1 per_point", tr_each)
        record_result!(results, :norm_g1, string(N), "batch", tr_batch)
        record_result!(results, :norm_g1, string(N), "each", tr_each)
    end

    println("\n== Batch normalization: batch_to_affine! vs per-point to_affine (G2) ==")
    g = g2_generator()
    rng = MersenneTwister(13)
    for N in (32, 128, 512)
        scalars = genscalars(rng, N)
        proj_pts = [scalar_mul(g, s) for s in scalars]
        println("N = ", N)
        _ = GrothCurves.batch_to_affine!(copy(proj_pts))
        _ = begin
            tmp = copy(proj_pts)
            for i in eachindex(tmp)
                x, y = to_affine(tmp[i])
                tmp[i] = G2Point(x, y, one(Fp2Element))
            end
        end
        tr_batch = @benchmark begin
            tmp = copy($proj_pts)
            GrothCurves.batch_to_affine!(tmp)
        end seconds = 1 samples = 10
        tr_each = @benchmark begin
            tmp = copy($proj_pts)
            for i in eachindex(tmp)
                x, y = to_affine(tmp[i])
                tmp[i] = G2Point(x, y, one(Fp2Element))
            end
        end seconds = 1 samples = 10
        print_stats("G2 batch_norm", tr_batch)
        print_stats("G2 per_point", tr_each)
        record_result!(results, :norm_g2, string(N), "batch", tr_batch)
        record_result!(results, :norm_g2, string(N), "each", tr_each)
    end
end

function bench_pairing_throughput(results)
    println("\n== Pairing engine: sequential vs batch ==")
    g1 = g1_generator()
    g2 = g2_generator()
    for N in (1, 4, 16)
        rng = MersenneTwister(200 + N)
        p_vec = [scalar_mul(g1, rand_scalar_nonzero(rng)) for _ in 1:N]
        q_vec = [scalar_mul(g2, rand_scalar_nonzero(rng)) for _ in 1:N]
        println("Batch size N = ", N)
        for (P, Q) in zip(p_vec, q_vec)
            pairing(ENGINE, P, Q)
        end
        pairing_batch(ENGINE, p_vec, q_vec)
        tr_seq = @benchmark begin
            acc = one(GTElement)
            @inbounds for i in eachindex($p_vec)
                acc *= pairing($ENGINE, $p_vec[i], $q_vec[i])
            end
            acc
        end seconds = 2 samples = 8
        tr_batch = @benchmark pairing_batch($ENGINE, $p_vec, $q_vec) seconds = 2 samples = 8
        print_stats("Pairing sequential (N=$(N))", tr_seq)
        print_stats("Pairing batch (N=$(N))", tr_batch)
        record_result!(results, :pairing, string(N), "sequential", tr_seq)
        record_result!(results, :pairing, string(N), "batch", tr_batch)
    end
end

function bench_pairing_micro(results)
    println("\n== Pairing engine micro-ops ==")
    inputs = bn254_stage0_inputs()
    P1 = inputs.pair_p
    Q1 = inputs.pair_q
    f_pre = miller_loop(ENGINE, P1, Q1)
    final_pre = final_exponentiation(ENGINE, f_pre)
    final_pre == pairing(ENGINE, P1, Q1) || error("Pairing micro-op benchmark setup mismatch")
    record_semantic!(results, :pairing_single, "miller_loop", serialize_bn254(f_pre))
    record_semantic!(results, :pairing_single, "pairing", serialize_bn254(final_pre))
    tr_single = @benchmark pairing($ENGINE, $P1, $Q1) seconds = 2 samples = 10
    tr_miller = @benchmark miller_loop($ENGINE, $P1, $Q1) seconds = 2 samples = 10
    tr_final = @benchmark final_exponentiation($ENGINE, $f_pre) seconds = 2 samples = 10
    print_stats("Pairing single", tr_single)
    print_stats("Miller loop", tr_miller)
    print_stats("Final exponentiation", tr_final)
    record_simple!(results, :pairing_single, "pairing", tr_single)
    record_simple!(results, :pairing_single, "miller_loop", tr_miller)
    record_simple!(results, :pairing_single, "final_exponentiation", tr_final)
end

function bench_py_ecc_primitives(results, meta)
    payload = run_py_ecc_bench(meta)
    payload === nothing && return

    py_results = payload["results"]
    py_semantic = payload["semantic"]

    println("\n== BN254 primitives: Groth.jl vs py_ecc ==")
    println("Using local py_ecc checkout at ", meta["py_ecc"]["py_ecc_root"])

    g1_base = scalar_mul(g1_generator(), det_py_ecc_scalar(101))
    g2_base = scalar_mul(g2_generator(), det_py_ecc_scalar(102))
    g1_scalar = det_py_ecc_scalar(201)
    g2_scalar = det_py_ecc_scalar(202)

    g1_expected = serialize_affine(scalar_mul(g1_base, g1_scalar))
    g2_expected = serialize_affine(scalar_mul(g2_base, g2_scalar))
    g1_expected == py_semantic["scalar_mul"]["g1"] || error("py_ecc G1 scalar multiplication result mismatch")
    g2_expected == py_semantic["scalar_mul"]["g2"] || error("py_ecc G2 scalar multiplication result mismatch")

    _ = scalar_mul(g1_base, g1_scalar)
    _ = scalar_mul(g2_base, g2_scalar)
    tr_g1 = @benchmark scalar_mul($g1_base, $g1_scalar) seconds = 1 samples = 10
    tr_g2 = @benchmark scalar_mul($g2_base, $g2_scalar) seconds = 1 samples = 10
    print_stats("Groth.jl G1 scalar", tr_g1)
    print_external_stats("py_ecc G1 scalar", py_results["scalar_mul"]["g1"])
    print_stats("Groth.jl G2 scalar", tr_g2)
    print_external_stats("py_ecc G2 scalar", py_results["scalar_mul"]["g2"])
    record_result!(results, :py_ecc_scalar, "g1", "groth_jl", tr_g1)
    record_external_result!(results, :py_ecc_scalar, "g1", "py_ecc", py_results["scalar_mul"]["g1"])
    record_result!(results, :py_ecc_scalar, "g2", "groth_jl", tr_g2)
    record_external_result!(results, :py_ecc_scalar, "g2", "py_ecc", py_results["scalar_mul"]["g2"])

    for N in PY_ECC_ACCUM_SIZES
        g1_bases, g1_scalars = py_ecc_accum_inputs(g1_generator(), 1000 + N, 2000 + N, N)
        g2_bases, g2_scalars = py_ecc_accum_inputs(g2_generator(), 3000 + N, 4000 + N, N)

        g1_expected_acc = serialize_affine(naive_msm(g1_bases, g1_scalars))
        g2_expected_acc = serialize_affine(naive_msm(g2_bases, g2_scalars))
        g1_expected_acc == py_semantic["naive_accum_g1"][string(N)] || error("py_ecc G1 naive accumulation result mismatch for N=$(N)")
        g2_expected_acc == py_semantic["naive_accum_g2"][string(N)] || error("py_ecc G2 naive accumulation result mismatch for N=$(N)")

        _ = naive_msm(g1_bases, g1_scalars)
        _ = naive_msm(g2_bases, g2_scalars)
        tr_g1_acc = @benchmark naive_msm($g1_bases, $g1_scalars) seconds = 1 samples = 10
        tr_g2_acc = @benchmark naive_msm($g2_bases, $g2_scalars) seconds = 1 samples = 10
        print_stats("Groth.jl G1 naive accum N=$(N)", tr_g1_acc)
        print_external_stats("py_ecc G1 naive accum N=$(N)", py_results["naive_accum_g1"][string(N)])
        print_stats("Groth.jl G2 naive accum N=$(N)", tr_g2_acc)
        print_external_stats("py_ecc G2 naive accum N=$(N)", py_results["naive_accum_g2"][string(N)])
        record_result!(results, :py_ecc_naive_accum_g1, string(N), "groth_jl", tr_g1_acc)
        record_external_result!(results, :py_ecc_naive_accum_g1, string(N), "py_ecc", py_results["naive_accum_g1"][string(N)])
        record_result!(results, :py_ecc_naive_accum_g2, string(N), "groth_jl", tr_g2_acc)
        record_external_result!(results, :py_ecc_naive_accum_g2, string(N), "py_ecc", py_results["naive_accum_g2"][string(N)])
    end

    pair_p = scalar_mul(g1_generator(), det_py_ecc_scalar(301))
    pair_q = scalar_mul(g2_generator(), det_py_ecc_scalar(302))
    pair_base = pairing(ENGINE, pair_p, pair_q)
    pair_double_g1 = pairing(ENGINE, scalar_mul(pair_p, BigInt(2)), pair_q)
    pair_double_g2 = pairing(ENGINE, pair_p, scalar_mul(pair_q, BigInt(2)))
    pair_double_g1 == pair_double_g2 || error("Groth.jl pairing bilinearity check failed for py_ecc benchmark workload")
    !isone(pair_base) || error("Groth.jl pairing benchmark workload is degenerate")
    Bool(py_semantic["pairing_checks"]["double_bilinear"]) || error("py_ecc pairing bilinearity check failed")
    Bool(py_semantic["pairing_checks"]["non_degenerate"]) || error("py_ecc pairing benchmark workload is degenerate")
    meta["py_ecc"]["pairing_validation"] = "validated via per-implementation bilinearity/non-degeneracy checks; direct coefficient equality is not enforced across implementations"

    _ = pairing(ENGINE, pair_p, pair_q)
    tr_pair = @benchmark pairing($ENGINE, $pair_p, $pair_q) seconds = 2 samples = 8
    print_stats("Groth.jl pairing", tr_pair)
    print_external_stats("py_ecc pairing", py_results["pairing"]["single"])
    record_result!(results, :py_ecc_pairing, "single", "groth_jl", tr_pair)
    record_external_result!(results, :py_ecc_pairing, "single", "py_ecc", py_results["pairing"]["single"])
end

function bn254_stage0_inputs()
    big(s::AbstractString) = parse(BigInt, s)

    fq_a = bn254_fq(big("1234567890123456789012345678901234567890"))
    fq_b = bn254_fq(big("987654321098765432109876543210987654321"))

    fr_a = bn254_fr(big("1357913579135791357913579135791357913579"))
    fr_b = bn254_fr(big("2468024680246802468024680246802468024680"))

    fp2_a = Fp2Element(fq_a, fq_b)
    fp2_b = Fp2Element(
        bn254_fq(big("3141592653589793238462643383279502884197")),
        bn254_fq(big("2718281828459045235360287471352662497757")),
    )

    fp6_a = Fp6Element(fp2_a, fp2_b, Fp2Element(bn254_fq(42), bn254_fq(99)))
    fp6_b = Fp6Element(Fp2Element(bn254_fq(777), bn254_fq(333)), fp2_a, fp2_b)

    fp12_a = Fp12Element(fp6_a, fp6_b)
    fp12_b = Fp12Element(fp6_b, fp6_a)

    return (
        fq_a = fq_a,
        fq_b = fq_b,
        fr_a = fr_a,
        fr_b = fr_b,
        fp2_a = fp2_a,
        fp2_b = fp2_b,
        fp6_a = fp6_a,
        fp6_b = fp6_b,
        fp12_a = fp12_a,
        fp12_b = fp12_b,
        g1_scalar = BigInt(1234567),
        g2_scalar = BigInt(7654321),
        pair_p = scalar_mul(g1_generator(), BigInt(3333333)),
        pair_q = scalar_mul(g2_generator(), BigInt(4444444)),
    )
end

function bench_bn254_primitives(results)
    println("\n== BN254 primitive fields, tower, and scalar multiplication ==")
    inputs = bn254_stage0_inputs()

    fq_mul = inputs.fq_a * inputs.fq_b
    fq_inv = inv(inputs.fq_a)
    fr_mul = inputs.fr_a * inputs.fr_b
    fr_inv = inv(inputs.fr_a)
    fp2_mul = inputs.fp2_a * inputs.fp2_b
    fp2_inv = inv(inputs.fp2_a)
    fp6_mul = inputs.fp6_a * inputs.fp6_b
    fp6_inv = inv(inputs.fp6_a)
    fp12_mul = inputs.fp12_a * inputs.fp12_b
    fp12_inv = inv(inputs.fp12_a)
    g1_gen = g1_generator()
    g2_gen = g2_generator()
    g1_scalar_point = scalar_mul(g1_gen, inputs.g1_scalar)
    g2_scalar_point = scalar_mul(g2_gen, inputs.g2_scalar)

    record_semantic!(results, :bn254_fq, "mul", serialize_bn254(fq_mul))
    record_semantic!(results, :bn254_fq, "inv", serialize_bn254(fq_inv))
    record_semantic!(results, :bn254_fr, "mul", serialize_bn254(fr_mul))
    record_semantic!(results, :bn254_fr, "inv", serialize_bn254(fr_inv))
    record_semantic!(results, :bn254_fp2, "mul", serialize_bn254(fp2_mul))
    record_semantic!(results, :bn254_fp2, "inv", serialize_bn254(fp2_inv))
    record_semantic!(results, :bn254_fp6, "mul", serialize_bn254(fp6_mul))
    record_semantic!(results, :bn254_fp6, "inv", serialize_bn254(fp6_inv))
    record_semantic!(results, :bn254_fp12, "mul", serialize_bn254(fp12_mul))
    record_semantic!(results, :bn254_fp12, "inv", serialize_bn254(fp12_inv))
    record_semantic!(results, :bn254_scalar_mul, "g1", serialize_affine(g1_scalar_point))
    record_semantic!(results, :bn254_scalar_mul, "g2", serialize_affine(g2_scalar_point))

    println("BN254Fq")
    _ = inputs.fq_a + inputs.fq_b
    _ = inputs.fq_a - inputs.fq_b
    _ = inputs.fq_a * inputs.fq_b
    _ = inputs.fq_a * inputs.fq_a
    _ = inv(inputs.fq_a)
    tr_fq_add = @benchmark $((inputs.fq_a)) + $((inputs.fq_b)) seconds = 1 samples = 10
    tr_fq_sub = @benchmark $((inputs.fq_a)) - $((inputs.fq_b)) seconds = 1 samples = 10
    tr_fq_mul = @benchmark $((inputs.fq_a)) * $((inputs.fq_b)) seconds = 1 samples = 10
    tr_fq_square = @benchmark $((inputs.fq_a)) * $((inputs.fq_a)) seconds = 1 samples = 10
    tr_fq_inv = @benchmark inv($((inputs.fq_a))) seconds = 1 samples = 10
    print_stats("BN254Fq add", tr_fq_add)
    print_stats("BN254Fq sub", tr_fq_sub)
    print_stats("BN254Fq mul", tr_fq_mul)
    print_stats("BN254Fq square", tr_fq_square)
    print_stats("BN254Fq inv", tr_fq_inv)
    record_simple!(results, :bn254_fq, "add", tr_fq_add)
    record_simple!(results, :bn254_fq, "sub", tr_fq_sub)
    record_simple!(results, :bn254_fq, "mul", tr_fq_mul)
    record_simple!(results, :bn254_fq, "square", tr_fq_square)
    record_simple!(results, :bn254_fq, "inv", tr_fq_inv)

    println("\nBN254Fr")
    _ = inputs.fr_a + inputs.fr_b
    _ = inputs.fr_a - inputs.fr_b
    _ = inputs.fr_a * inputs.fr_b
    _ = inputs.fr_a * inputs.fr_a
    _ = inv(inputs.fr_a)
    tr_fr_add = @benchmark $((inputs.fr_a)) + $((inputs.fr_b)) seconds = 1 samples = 10
    tr_fr_sub = @benchmark $((inputs.fr_a)) - $((inputs.fr_b)) seconds = 1 samples = 10
    tr_fr_mul = @benchmark $((inputs.fr_a)) * $((inputs.fr_b)) seconds = 1 samples = 10
    tr_fr_square = @benchmark $((inputs.fr_a)) * $((inputs.fr_a)) seconds = 1 samples = 10
    tr_fr_inv = @benchmark inv($((inputs.fr_a))) seconds = 1 samples = 10
    print_stats("BN254Fr add", tr_fr_add)
    print_stats("BN254Fr sub", tr_fr_sub)
    print_stats("BN254Fr mul", tr_fr_mul)
    print_stats("BN254Fr square", tr_fr_square)
    print_stats("BN254Fr inv", tr_fr_inv)
    record_simple!(results, :bn254_fr, "add", tr_fr_add)
    record_simple!(results, :bn254_fr, "sub", tr_fr_sub)
    record_simple!(results, :bn254_fr, "mul", tr_fr_mul)
    record_simple!(results, :bn254_fr, "square", tr_fr_square)
    record_simple!(results, :bn254_fr, "inv", tr_fr_inv)

    println("\nFp2")
    _ = inputs.fp2_a + inputs.fp2_b
    _ = inputs.fp2_a * inputs.fp2_b
    _ = inputs.fp2_a * inputs.fp2_a
    _ = inv(inputs.fp2_a)
    tr_fp2_add = @benchmark $((inputs.fp2_a)) + $((inputs.fp2_b)) seconds = 1 samples = 10
    tr_fp2_mul = @benchmark $((inputs.fp2_a)) * $((inputs.fp2_b)) seconds = 1 samples = 10
    tr_fp2_square = @benchmark $((inputs.fp2_a)) * $((inputs.fp2_a)) seconds = 1 samples = 10
    tr_fp2_inv = @benchmark inv($((inputs.fp2_a))) seconds = 1 samples = 10
    print_stats("Fp2 add", tr_fp2_add)
    print_stats("Fp2 mul", tr_fp2_mul)
    print_stats("Fp2 square", tr_fp2_square)
    print_stats("Fp2 inv", tr_fp2_inv)
    record_simple!(results, :bn254_fp2, "add", tr_fp2_add)
    record_simple!(results, :bn254_fp2, "mul", tr_fp2_mul)
    record_simple!(results, :bn254_fp2, "square", tr_fp2_square)
    record_simple!(results, :bn254_fp2, "inv", tr_fp2_inv)

    println("\nFp6")
    _ = inputs.fp6_a + inputs.fp6_b
    _ = inputs.fp6_a * inputs.fp6_b
    _ = square(inputs.fp6_a)
    _ = inv(inputs.fp6_a)
    tr_fp6_add = @benchmark $((inputs.fp6_a)) + $((inputs.fp6_b)) seconds = 1 samples = 10
    tr_fp6_mul = @benchmark $((inputs.fp6_a)) * $((inputs.fp6_b)) seconds = 1 samples = 10
    tr_fp6_square = @benchmark square($((inputs.fp6_a))) seconds = 1 samples = 10
    tr_fp6_inv = @benchmark inv($((inputs.fp6_a))) seconds = 1 samples = 10
    print_stats("Fp6 add", tr_fp6_add)
    print_stats("Fp6 mul", tr_fp6_mul)
    print_stats("Fp6 square", tr_fp6_square)
    print_stats("Fp6 inv", tr_fp6_inv)
    record_simple!(results, :bn254_fp6, "add", tr_fp6_add)
    record_simple!(results, :bn254_fp6, "mul", tr_fp6_mul)
    record_simple!(results, :bn254_fp6, "square", tr_fp6_square)
    record_simple!(results, :bn254_fp6, "inv", tr_fp6_inv)

    println("\nFp12")
    _ = inputs.fp12_a + inputs.fp12_b
    _ = inputs.fp12_a * inputs.fp12_b
    _ = square(inputs.fp12_a)
    _ = inv(inputs.fp12_a)
    tr_fp12_add = @benchmark $((inputs.fp12_a)) + $((inputs.fp12_b)) seconds = 1 samples = 10
    tr_fp12_mul = @benchmark $((inputs.fp12_a)) * $((inputs.fp12_b)) seconds = 1 samples = 10
    tr_fp12_square = @benchmark square($((inputs.fp12_a))) seconds = 1 samples = 10
    tr_fp12_inv = @benchmark inv($((inputs.fp12_a))) seconds = 1 samples = 10
    print_stats("Fp12 add", tr_fp12_add)
    print_stats("Fp12 mul", tr_fp12_mul)
    print_stats("Fp12 square", tr_fp12_square)
    print_stats("Fp12 inv", tr_fp12_inv)
    record_simple!(results, :bn254_fp12, "add", tr_fp12_add)
    record_simple!(results, :bn254_fp12, "mul", tr_fp12_mul)
    record_simple!(results, :bn254_fp12, "square", tr_fp12_square)
    record_simple!(results, :bn254_fp12, "inv", tr_fp12_inv)

    println("\nCurve scalar multiplication")
    _ = scalar_mul(g1_gen, inputs.g1_scalar)
    _ = scalar_mul(g2_gen, inputs.g2_scalar)
    tr_g1_scalar = @benchmark scalar_mul($g1_gen, $(inputs.g1_scalar)) seconds = 1 samples = 10
    tr_g2_scalar = @benchmark scalar_mul($g2_gen, $(inputs.g2_scalar)) seconds = 1 samples = 10
    print_stats("G1 scalar", tr_g1_scalar)
    print_stats("G2 scalar", tr_g2_scalar)
    record_simple!(results, :bn254_scalar_mul, "g1", tr_g1_scalar)
    record_simple!(results, :bn254_scalar_mul, "g2", tr_g2_scalar)
end

function bench_groth16(results)
    println("\n== Groth16 end-to-end pipeline ==")
    r1cs = create_r1cs_example_sum_of_products()
    tr_r1cs_to_qap = @benchmark r1cs_to_qap($r1cs) seconds = 2 samples = 10
    print_stats("R1CS -> QAP", tr_r1cs_to_qap)
    record_simple!(results, :groth16, "r1cs_to_qap", tr_r1cs_to_qap)

    qap = r1cs_to_qap(r1cs)
    witness = create_witness_sum_of_products(3, 5, 7, 11)
    public_inputs = r1cs.num_public > 1 ? witness.values[2:r1cs.num_public] : eltype(witness.values)[]

    tr_setup = @benchmark setup_full($qap; rng=MersenneTwister(42)) seconds = 5 samples = 5
    print_stats("Groth16 setup", tr_setup)
    record_simple!(results, :groth16, "setup", tr_setup)

    keypair = setup_full(qap; rng=MersenneTwister(42))
    proof = prove_full(keypair.pk, qap, witness; rng=MersenneTwister(1337))

    tr_prove = @benchmark prove_full($keypair.pk, $qap, $witness; rng=MersenneTwister(1337)) seconds = 5 samples = 5
    print_stats("Groth16 prove", tr_prove)
    record_simple!(results, :groth16, "prove", tr_prove)

    tr_verify = @benchmark verify_full($keypair.vk, $proof, $public_inputs) seconds = 5 samples = 15
    print_stats("Groth16 verify", tr_verify)
    record_simple!(results, :groth16, "verify_full", tr_verify)

    tr_prepare_vk = @benchmark prepare_verifying_key($keypair.vk) seconds = 2 samples = 10
    print_stats("Groth16 prepare_vk", tr_prepare_vk)
    record_simple!(results, :groth16, "prepare_vk", tr_prepare_vk)

    pvk = prepare_verifying_key(keypair.vk)
    tr_prepare_inputs = @benchmark prepare_inputs($pvk, $public_inputs) seconds = 2 samples = 10
    print_stats("Groth16 prepare_inputs", tr_prepare_inputs)
    record_simple!(results, :groth16, "prepare_inputs", tr_prepare_inputs)

    prepared_inputs = prepare_inputs(pvk, public_inputs)
    tr_verify_prepared = @benchmark verify_with_prepared($pvk, $proof, $prepared_inputs) seconds = 5 samples = 15
    print_stats("Groth16 verify prepared", tr_verify_prepared)
    record_simple!(results, :groth16, "verify_prepared", tr_verify_prepared)
end

function bench_prove_full(results)
    println("\n== prove_full fixture baseline and phase breakdown ==")
    fixtures = default_prove_full_fixtures()

    for fixture in fixtures
        name = fixture.name
        desc = fixture.description
        qap = fixture.qap
        witness = fixture.witness
        pk = fixture.keypair.pk
        prove_seed = fixture.prove_seed

        println("Fixture: $(name)")
        println("  ", desc)
        println("  constraints=$(fixture.r1cs.num_constraints) vars=$(fixture.r1cs.num_vars) public=$(fixture.r1cs.num_public) domain=$(qap.domain.size)")

        record_fixture_meta!(results, :prove_full, name, fixture_metadata(fixture))

        _ = prove_full(pk, qap, witness; rng=MersenneTwister(prove_seed))
        tr_end = @benchmark prove_full($pk, $qap, $witness; rng=MersenneTwister($prove_seed)) seconds = 3 samples = 6
        print_stats("prove_full[$(name)] end", tr_end)
        record_fixture_trial!(results, :prove_full, name, "end_to_end", tr_end)

        _ = witness_to_scalars(witness)
        tr_scalars = @benchmark witness_to_scalars($witness) seconds = 1 samples = 10
        print_stats("prove_full[$(name)] scalars", tr_scalars)
        record_fixture_trial!(results, :prove_full, name, "witness_to_scalars", tr_scalars)

        scalars = witness_to_scalars(witness)
        _ = GrothAlgebra.multi_scalar_mul(pk.A_query_g1, scalars)
        _ = GrothAlgebra.multi_scalar_mul(pk.B_query_g1, scalars)
        _ = GrothAlgebra.multi_scalar_mul(pk.B_query_g2, scalars)
        tr_msm_a = @benchmark GrothAlgebra.multi_scalar_mul($pk.A_query_g1, $scalars) seconds = 1 samples = 10
        tr_msm_b1 = @benchmark GrothAlgebra.multi_scalar_mul($pk.B_query_g1, $scalars) seconds = 1 samples = 10
        tr_msm_b2 = @benchmark GrothAlgebra.multi_scalar_mul($pk.B_query_g2, $scalars) seconds = 1 samples = 10
        print_stats("prove_full[$(name)] A_msm", tr_msm_a)
        print_stats("prove_full[$(name)] B1_msm", tr_msm_b1)
        print_stats("prove_full[$(name)] B2_msm", tr_msm_b2)
        record_fixture_trial!(results, :prove_full, name, "msm_a_g1", tr_msm_a)
        record_fixture_trial!(results, :prove_full, name, "msm_b_g1", tr_msm_b1)
        record_fixture_trial!(results, :prove_full, name, "msm_b_g2", tr_msm_b2)

        _ = compute_h_polynomial(qap, witness)
        tr_h_total = @benchmark compute_h_polynomial($qap, $witness) seconds = 2 samples = 8
        print_stats("prove_full[$(name)] h_total", tr_h_total)
        record_fixture_trial!(results, :prove_full, name, "compute_h_total", tr_h_total)

        _ = build_combined_polynomials(qap, witness)
        tr_h_assemble = @benchmark build_combined_polynomials($qap, $witness) seconds = 1 samples = 10
        print_stats("prove_full[$(name)] h_assemble", tr_h_assemble)
        record_fixture_trial!(results, :prove_full, name, "h_poly_assembly", tr_h_assemble)

        u_poly, v_poly, w_poly = build_combined_polynomials(qap, witness)
        _ = compute_h_dense_path(qap, u_poly, v_poly, w_poly)
        tr_h_dense = @benchmark compute_h_dense_path($qap, $u_poly, $v_poly, $w_poly) seconds = 1 samples = 10
        print_stats("prove_full[$(name)] h_dense", tr_h_dense)
        record_fixture_trial!(results, :prove_full, name, "h_dense_quotient", tr_h_dense)

        dense_result = compute_h_dense_path(qap, u_poly, v_poly, w_poly)
        dense_len = length(dense_result.coeffs)
        _ = compute_h_coset_path(qap, u_poly, v_poly, w_poly, dense_len)
        tr_h_coset = @benchmark compute_h_coset_path($qap, $u_poly, $v_poly, $w_poly, $dense_len) seconds = 1 samples = 10
        print_stats("prove_full[$(name)] h_coset", tr_h_coset)
        record_fixture_trial!(results, :prove_full, name, "h_coset_fft", tr_h_coset)

        coset_result = compute_h_coset_path(qap, u_poly, v_poly, w_poly, dense_len)
        _ = assert_h_parity(dense_result, coset_result)
        tr_h_parity = @benchmark assert_h_parity($dense_result, $coset_result) seconds = 1 samples = 10
        print_stats("prove_full[$(name)] h_parity", tr_h_parity)
        record_fixture_trial!(results, :prove_full, name, "h_parity_assert", tr_h_parity)

        h_poly = assert_h_parity(dense_result, coset_result)
        _ = h_msm(pk, h_poly)
        tr_h_msm = @benchmark h_msm($pk, $h_poly) seconds = 1 samples = 10
        print_stats("prove_full[$(name)] H_msm", tr_h_msm)
        record_fixture_trial!(results, :prove_full, name, "h_msm", tr_h_msm)

        _ = l_msm(pk, witness)
        tr_l_msm = @benchmark l_msm($pk, $witness) seconds = 1 samples = 10
        print_stats("prove_full[$(name)] L_msm", tr_l_msm)
        record_fixture_trial!(results, :prove_full, name, "l_msm", tr_l_msm)

        r_rand, s_rand = sample_prove_randomizers(fixture)
        accum = prove_query_accumulators(pk, scalars)
        ab_terms = assemble_ab_terms(pk, accum.A_acc_g1, accum.B_acc_g1, accum.B_acc_g2, r_rand, s_rand)
        H = h_msm(pk, h_poly)
        L = l_msm(pk, witness)
        _ = assemble_c(pk, ab_terms.A, ab_terms.B1_g1, H, L, r_rand, s_rand)
        tr_final_c = @benchmark assemble_c($pk, $(ab_terms.A), $(ab_terms.B1_g1), $H, $L, $r_rand, $s_rand) seconds = 1 samples = 10
        print_stats("prove_full[$(name)] C_final", tr_final_c)
        record_fixture_trial!(results, :prove_full, name, "final_c", tr_final_c)
    end
end

# -----------------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------------

function main()
    cli = parse_cli_args(ARGS)
    if cli.list_profiles
        print_available_profiles()
        return
    end

    println("GrothBenchmarks — Julia $(VERSION)")
    println("Benchmark profile: ", cli.profile)
    println("Selected groups: ", join(string.(cli.selected_groups), ", "))
    run_id = Dates.format(Dates.now(Dates.UTC), "yyyy-mm-dd_HHMMSS")
    artifact_dir = joinpath(@__DIR__, "artifacts", run_id)
    results_dir = joinpath(artifact_dir, "results")
    plots_dir = joinpath(artifact_dir, "plots")
    meta_dir = joinpath(artifact_dir, "meta")
    mkpath(results_dir)
    mkpath(plots_dir)
    mkpath(meta_dir)

    results = Dict{Symbol,Any}()
    meta = run_metadata()
    meta["run_id"] = run_id
    meta["artifact_dir"] = artifact_dir
    meta["benchmark_profile"] = cli.profile
    meta["benchmark_groups"] = string.(cli.selected_groups)
    meta["partial_run"] = cli.profile != "full" || cli.selected_groups != BENCHMARK_PROFILES["full"]
    results[:_meta] = meta
    results[:_semantic] = Dict{String,Any}()
    selected = Set(cli.selected_groups)

    :bn254_primitives in selected && bench_bn254_primitives(results)
    :fixed_base in selected && bench_fixed_base(results)
    :variable_msm in selected && bench_variable_msm(results)
    :batch_norm in selected && bench_batch_norm(results)
    :pairing_throughput in selected && bench_pairing_throughput(results)
    :pairing_micro in selected && bench_pairing_micro(results)
    :py_ecc_primitives in selected && bench_py_ecc_primitives(results, meta)
    :groth16 in selected && bench_groth16(results)
    :prove_full in selected && bench_prove_full(results)

    results_out = joinpath(results_dir, "benchmark_results.json")
    open(results_out, "w") do io
        JSON.print(io, results)
    end
    meta_out = joinpath(meta_dir, "run_meta.json")
    open(meta_out, "w") do io
        JSON.print(io, meta)
    end
    println("\nSaved results to ", results_out)
    println("Saved metadata to ", meta_out)
    println("Plots should be written to ", plots_dir, " via benchmarks/plot.jl")
end

main()

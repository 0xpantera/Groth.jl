#!/usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__)

using BenchmarkTools
using Random
using JSON
using Dates
using GrothAlgebra
using GrothCurves
using GrothProofs

const Fr = BN254ScalarField
const r = BN254_ORDER_R
const ENGINE = BN254_ENGINE

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
    factor = unit == "s"  ? 1.0 :
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

function print_stats(label, tr::BenchmarkTools.Trial)
    tmin = minimum(tr)
    tmed = median(tr)
    tmean = mean(tr)
    bytes = memory(tr)
    println(rpad(label, 40), " min=$(tmin) med=$(tmed) mean=$(tmean) mem=$(bytes)B")
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
        tr_build = @benchmark build_fixed_table($g1) seconds=1 samples=10
        tab = build_fixed_table(g1)
        tr_naive = @benchmark genpoints_fixed($g1, $scalars) seconds=1 samples=10
        tr_batch = @benchmark batch_mul($tab, $scalars) seconds=1 samples=10
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
        tr_build = @benchmark build_fixed_table($g2) seconds=1 samples=10
        tab = build_fixed_table(g2)
        tr_naive = @benchmark genpoints_fixed($g2, $scalars) seconds=1 samples=10
        tr_batch = @benchmark batch_mul($tab, $scalars) seconds=1 samples=10
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
        tr_naive = @benchmark naive_msm($bases, $scalars) seconds=1 samples=10
        tr_msm = @benchmark GrothAlgebra.multi_scalar_mul($bases, $scalars) seconds=1 samples=10
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
        tr_naive = @benchmark naive_msm($bases, $scalars) seconds=1 samples=10
        tr_msm = @benchmark GrothAlgebra.multi_scalar_mul($bases, $scalars) seconds=1 samples=10
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
                tmp[i] = G1Point(x, y, one(BN254Field))
            end
        end
        tr_batch = @benchmark begin
            tmp = copy($proj_pts)
            GrothCurves.batch_to_affine!(tmp)
        end seconds=1 samples=10
        tr_each = @benchmark begin
            tmp = copy($proj_pts)
            for i in eachindex(tmp)
                x, y = to_affine(tmp[i])
                tmp[i] = G1Point(x, y, one(BN254Field))
            end
        end seconds=1 samples=10
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
        end seconds=1 samples=10
        tr_each = @benchmark begin
            tmp = copy($proj_pts)
            for i in eachindex(tmp)
                x, y = to_affine(tmp[i])
                tmp[i] = G2Point(x, y, one(Fp2Element))
            end
        end seconds=1 samples=10
        print_stats("G2 batch_norm", tr_batch)
        print_stats("G2 per_point", tr_each)
        record_result!(results, :norm_g2, string(N), "batch", tr_batch)
        record_result!(results, :norm_g2, string(N), "each", tr_each)
    end
end

function bench_pairings(results)
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
        end seconds=2 samples=8
        tr_batch = @benchmark pairing_batch($ENGINE, $p_vec, $q_vec) seconds=2 samples=8
        print_stats("Pairing sequential (N=$(N))", tr_seq)
        print_stats("Pairing batch (N=$(N))", tr_batch)
        record_result!(results, :pairing, string(N), "sequential", tr_seq)
        record_result!(results, :pairing, string(N), "batch", tr_batch)
    end

    println("\n== Pairing engine micro-ops ==")
    P1 = scalar_mul(g1, BigInt(5))
    Q1 = scalar_mul(g2, BigInt(7))
    f_pre = miller_loop(ENGINE, P1, Q1)
    tr_single = @benchmark pairing($ENGINE, $P1, $Q1) seconds=2 samples=10
    tr_miller = @benchmark miller_loop($ENGINE, $P1, $Q1) seconds=2 samples=10
    tr_final = @benchmark final_exponentiation($ENGINE, $f_pre) seconds=2 samples=10
    print_stats("Pairing single", tr_single)
    print_stats("Miller loop", tr_miller)
    print_stats("Final exponentiation", tr_final)
    record_simple!(results, :pairing_single, "pairing", tr_single)
    record_simple!(results, :pairing_single, "miller_loop", tr_miller)
    record_simple!(results, :pairing_single, "final_exponentiation", tr_final)
end

function bench_groth16(results)
    println("\n== Groth16 end-to-end pipeline ==")
    r1cs = create_r1cs_example_sum_of_products()
    tr_r1cs_to_qap = @benchmark r1cs_to_qap($r1cs) seconds=2 samples=10
    print_stats("R1CS -> QAP", tr_r1cs_to_qap)
    record_simple!(results, :groth16, "r1cs_to_qap", tr_r1cs_to_qap)

    qap = r1cs_to_qap(r1cs)
    witness = create_witness_sum_of_products(3, 5, 7, 11)
    public_inputs = witness.values[1:r1cs.num_public]

    tr_setup = @benchmark setup_full($qap; rng=MersenneTwister(42)) seconds=5 samples=5
    print_stats("Groth16 setup", tr_setup)
    record_simple!(results, :groth16, "setup", tr_setup)

    keypair = setup_full(qap; rng=MersenneTwister(42))
    proof = prove_full(keypair.pk, qap, witness; rng=MersenneTwister(1337))

    tr_prove = @benchmark prove_full($keypair.pk, $qap, $witness; rng=MersenneTwister(1337)) seconds=5 samples=5
    print_stats("Groth16 prove", tr_prove)
    record_simple!(results, :groth16, "prove", tr_prove)

    tr_verify = @benchmark verify_full($keypair.vk, $proof, $public_inputs) seconds=5 samples=15
    print_stats("Groth16 verify", tr_verify)
    record_simple!(results, :groth16, "verify_full", tr_verify)

    tr_prepare_vk = @benchmark prepare_verifying_key($keypair.vk) seconds=2 samples=10
    print_stats("Groth16 prepare_vk", tr_prepare_vk)
    record_simple!(results, :groth16, "prepare_vk", tr_prepare_vk)

    pvk = prepare_verifying_key(keypair.vk)
    tr_prepare_inputs = @benchmark prepare_inputs($pvk, $public_inputs) seconds=2 samples=10
    print_stats("Groth16 prepare_inputs", tr_prepare_inputs)
    record_simple!(results, :groth16, "prepare_inputs", tr_prepare_inputs)

    prepared_inputs = prepare_inputs(pvk, public_inputs)
    tr_verify_prepared = @benchmark verify_with_prepared($pvk, $proof, $prepared_inputs) seconds=5 samples=15
    print_stats("Groth16 verify prepared", tr_verify_prepared)
    record_simple!(results, :groth16, "verify_prepared", tr_verify_prepared)
end

# -----------------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------------

function main()
    println("GrothBenchmarks — Julia $(VERSION)")
    results = Dict{Symbol,Any}()
    bench_fixed_base(results)
    bench_variable_msm(results)
    bench_batch_norm(results)
    bench_pairings(results)
    bench_groth16(results)

    ts = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS")
    out = joinpath(@__DIR__, "results_" * ts * ".json")
    open(out, "w") do io
        JSON.print(io, results)
    end
    println("\nSaved results to ", out)
end

main()

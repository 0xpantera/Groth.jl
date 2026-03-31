#!/usr/bin/env julia

using BenchmarkTools
using Random
using JSON
using Dates
using Sockets
using GrothAlgebra
using GrothCurves
using GrothProofs

include("prove_full_common.jl")

const Fr = BN254Fr
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

function print_stats(label, tr::BenchmarkTools.Trial)
    tmin = minimum(tr)
    tmed = median(tr)
    tmean = mean(tr)
    bytes = memory(tr)
    println(rpad(label, 40), " min=$(tmin) med=$(tmed) mean=$(tmean) mem=$(bytes)B")
end

function try_readchomp(cmd::Cmd)
    try
        return strip(readchomp(cmd))
    catch
        return "unknown"
    end
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
        end seconds = 2 samples = 8
        tr_batch = @benchmark pairing_batch($ENGINE, $p_vec, $q_vec) seconds = 2 samples = 8
        print_stats("Pairing sequential (N=$(N))", tr_seq)
        print_stats("Pairing batch (N=$(N))", tr_batch)
        record_result!(results, :pairing, string(N), "sequential", tr_seq)
        record_result!(results, :pairing, string(N), "batch", tr_batch)
    end

    println("\n== Pairing engine micro-ops ==")
    P1 = scalar_mul(g1, BigInt(5))
    Q1 = scalar_mul(g2, BigInt(7))
    f_pre = miller_loop(ENGINE, P1, Q1)
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
    println("GrothBenchmarks — Julia $(VERSION)")
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
    results[:_meta] = meta
    bench_fixed_base(results)
    bench_variable_msm(results)
    bench_batch_norm(results)
    bench_pairings(results)
    bench_groth16(results)
    bench_prove_full(results)

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

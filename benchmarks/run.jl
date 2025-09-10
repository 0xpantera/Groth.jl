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

function rand_scalar(R::AbstractRNG)
    # Sample in [0, r)
    return BigInt(rand(R, 0:r-1))
end

function genscalars(R::AbstractRNG, N::Int)
    return [rand_scalar(R) for _ in 1:N]
end

function genpoints_fixed(g, scalars::Vector{BigInt})
    # Naive scalar_mul
    return [scalar_mul(g, s) for s in scalars]
end

function print_stats(label, tr::BenchmarkTools.Trial)
    tmin = minimum(tr)
    tmed = median(tr)
    tmean = mean(tr)
    bytes = memory(tr)
    println(rpad(label, 40), " min=$(tmin) med=$(tmed) mean=$(tmean) mem=$(bytes)B")
end

function bench_fixed_base(results)
    println("\n== Fixed-base precompute: table build and batch_mul (G1) ==")
    g1 = g1_generator()
    R = MersenneTwister(42)
    for N in (32, 128, 512)
        scalars = genscalars(R, N)
        println("N = ", N)
        # Warmup
        _ = build_fixed_table(g1)
        _ = genpoints_fixed(g1, scalars)
        # Separate table build vs batch_mul
        tr_build = @benchmark build_fixed_table($g1) seconds=1 samples=10
        tab = build_fixed_table(g1)
        tr_naive = @benchmark genpoints_fixed($g1, $scalars) seconds=1 samples=10
        tr_batch = @benchmark batch_mul($tab, $scalars) seconds=1 samples=10
        print_stats("G1 build", tr_build); print_stats("G1 naive", tr_naive); print_stats("G1 batch", tr_batch)
        results[:fixed_g1][string(N)] = Dict("build"=>Dict("min"=>string(minimum(tr_build)),"med"=>string(median(tr_build))),
                                             "naive"=>Dict("min"=>string(minimum(tr_naive)),"med"=>string(median(tr_naive))),
                                             "batch"=>Dict("min"=>string(minimum(tr_batch)),"med"=>string(median(tr_batch))))
    end

    println("\n== Fixed-base precompute: table build and batch_mul (G2) ==")
    g2 = g2_generator()
    R = MersenneTwister(123)
    for N in (32, 128, 512)
        scalars = genscalars(R, N)
        println("N = ", N)
        _ = build_fixed_table(g2)
        _ = genpoints_fixed(g2, scalars)
        tr_build = @benchmark build_fixed_table($g2) seconds=1 samples=10
        tab = build_fixed_table(g2)
        tr_naive = @benchmark genpoints_fixed($g2, $scalars) seconds=1 samples=10
        tr_batch = @benchmark batch_mul($tab, $scalars) seconds=1 samples=10
        print_stats("G2 build", tr_build); print_stats("G2 naive", tr_naive); print_stats("G2 batch", tr_batch)
        results[:fixed_g2][string(N)] = Dict("build"=>Dict("min"=>string(minimum(tr_build)),"med"=>string(median(tr_build))),
                                             "naive"=>Dict("min"=>string(minimum(tr_naive)),"med"=>string(median(tr_naive))),
                                             "batch"=>Dict("min"=>string(minimum(tr_batch)),"med"=>string(median(tr_batch))))
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

function gen_random_bases_g1(R::AbstractRNG, N::Int)
    g = g1_generator()
    scal = genscalars(R, N)
    return [scalar_mul(g, s) for s in scal]
end

function gen_random_bases_g2(R::AbstractRNG, N::Int)
    g = g2_generator()
    scal = genscalars(R, N)
    return [scalar_mul(g, s) for s in scal]
end

function bench_variable_msm(results)
    println("\n== Variable-base MSM: multi_scalar_mul vs naive (G1) ==")
    R = MersenneTwister(7)
    for N in (32, 128, 512)
        bases = gen_random_bases_g1(R, N)
        scalars = genscalars(R, N)
        println("N = ", N)
        # Warmup
        _ = naive_msm(bases, scalars)
        _ = GrothAlgebra.multi_scalar_mul(bases, scalars)
        tr_naive = @benchmark naive_msm($bases, $scalars) seconds=1 samples=10
        tr_msm = @benchmark GrothAlgebra.multi_scalar_mul($bases, $scalars) seconds=1 samples=10
        print_stats("G1 naive", tr_naive); print_stats("G1 MSM", tr_msm)
        results[:msm_g1][string(N)] = Dict("naive"=>Dict("min"=>string(minimum(tr_naive)),"med"=>string(median(tr_naive))),
                                           "msm"=>Dict("min"=>string(minimum(tr_msm)),"med"=>string(median(tr_msm))))
    end

    println("\n== Variable-base MSM: multi_scalar_mul vs naive (G2) ==")
    R = MersenneTwister(9)
    for N in (32, 128, 512)
        bases = gen_random_bases_g2(R, N)
        scalars = genscalars(R, N)
        println("N = ", N)
        _ = naive_msm(bases, scalars)
        _ = GrothAlgebra.multi_scalar_mul(bases, scalars)
        tr_naive = @benchmark naive_msm($bases, $scalars) seconds=1 samples=10
        tr_msm = @benchmark GrothAlgebra.multi_scalar_mul($bases, $scalars) seconds=1 samples=10
        print_stats("G2 naive", tr_naive); print_stats("G2 MSM", tr_msm)
        results[:msm_g2][string(N)] = Dict("naive"=>Dict("min"=>string(minimum(tr_naive)),"med"=>string(median(tr_naive))),
                                           "msm"=>Dict("min"=>string(minimum(tr_msm)),"med"=>string(median(tr_msm))))
    end
end

function bench_batch_norm(results)
    println("\n== Batch normalization: batch_to_affine! vs per-point to_affine (G1) ==")
    g = g1_generator()
    R = MersenneTwister(11)
    for N in (32, 128, 512)
        scalars = genscalars(R, N)
        pts = [scalar_mul(g, s) for s in scalars]
        pts2 = copy(pts)
        println("N = ", N)
        _ = GrothCurves.batch_to_affine!(copy(pts))
        _ = begin
            tmp = copy(pts2)
            for i in eachindex(tmp); x,y = to_affine(tmp[i]); tmp[i] = G1Point(x,y,one(BN254Field)); end
        end
        tr_batch = @benchmark GrothCurves.batch_to_affine!($pts) seconds=1 samples=10
        tr_each = @benchmark begin
            tmp = copy($pts2)
            for i in eachindex(tmp)
                x,y = to_affine(tmp[i])
                tmp[i] = G1Point(x,y,one(BN254Field))
            end
        end seconds=1 samples=10
        print_stats("G1 batch_norm", tr_batch); print_stats("G1 per_point", tr_each)
        results[:norm_g1][string(N)] = Dict("batch"=>Dict("min"=>string(minimum(tr_batch)),"med"=>string(median(tr_batch))),
                                            "each"=>Dict("min"=>string(minimum(tr_each)),"med"=>string(median(tr_each))))
    end

    println("\n== Batch normalization: batch_to_affine! vs per-point to_affine (G2) ==")
    g = g2_generator()
    R = MersenneTwister(13)
    for N in (32, 128, 512)
        scalars = genscalars(R, N)
        pts = [scalar_mul(g, s) for s in scalars]
        pts2 = copy(pts)
        println("N = ", N)
        _ = GrothCurves.batch_to_affine!(copy(pts))
        _ = begin
            tmp = copy(pts2)
            for i in eachindex(tmp); x,y = to_affine(tmp[i]); tmp[i] = G2Point(x,y,one(Fp2Element)); end
        end
        tr_batch = @benchmark GrothCurves.batch_to_affine!($pts) seconds=1 samples=10
        tr_each = @benchmark begin
            tmp = copy($pts2)
            for i in eachindex(tmp)
                x,y = to_affine(tmp[i])
                tmp[i] = G2Point(x,y,one(Fp2Element))
            end
        end seconds=1 samples=10
        print_stats("G2 batch_norm", tr_batch); print_stats("G2 per_point", tr_each)
        results[:norm_g2][string(N)] = Dict("batch"=>Dict("min"=>string(minimum(tr_batch)),"med"=>string(median(tr_batch))),
                                            "each"=>Dict("min"=>string(minimum(tr_each)),"med"=>string(median(tr_each))))
    end
end

function main()
    println("GrothBenchmarks — Julia $(VERSION)")
    results = Dict(
        :fixed_g1=>Dict(), :fixed_g2=>Dict(),
        :msm_g1=>Dict(), :msm_g2=>Dict(),
        :norm_g1=>Dict(), :norm_g2=>Dict(),
    )
    bench_fixed_base(results)
    bench_variable_msm(results)
    bench_batch_norm(results)
    # TODO: add end-to-end setup_full / prove_full / verify_full timings
    # Save summary JSON
    ts = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS")
    out = joinpath(@__DIR__, "results_"*ts*".json")
    open(out, "w") do io
        JSON.print(io, results)
    end
    println("\nSaved results to ", out)
end

main()

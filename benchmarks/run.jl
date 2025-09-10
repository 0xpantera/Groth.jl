#!/usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__)

using BenchmarkTools
using Random
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

function bench_fixed_base()
    println("\n== Fixed-base precompute: batch_mul vs scalar_mul (G1) ==")
    g1 = g1_generator()
    R = MersenneTwister(42)
    for N in (32, 128, 512)
        scalars = genscalars(R, N)
        @show N
        @btime genpoints_fixed($g1, $scalars);
        tab = build_fixed_table(g1)
        @btime batch_mul($tab, $scalars);
    end

    println("\n== Fixed-base precompute: batch_mul vs scalar_mul (G2) ==")
    g2 = g2_generator()
    R = MersenneTwister(123)
    for N in (32, 128, 512)
        scalars = genscalars(R, N)
        @show N
        @btime genpoints_fixed($g2, $scalars);
        tab = build_fixed_table(g2)
        @btime batch_mul($tab, $scalars);
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

function bench_variable_msm()
    println("\n== Variable-base MSM: multi_scalar_mul vs naive (G1) ==")
    R = MersenneTwister(7)
    for N in (32, 128, 512)
        bases = gen_random_bases_g1(R, N)
        scalars = genscalars(R, N)
        @show N
        @btime naive_msm($bases, $scalars);
        @btime GrothAlgebra.multi_scalar_mul($bases, $scalars);
    end

    println("\n== Variable-base MSM: multi_scalar_mul vs naive (G2) ==")
    R = MersenneTwister(9)
    for N in (32, 128, 512)
        bases = gen_random_bases_g2(R, N)
        scalars = genscalars(R, N)
        @show N
        @btime naive_msm($bases, $scalars);
        @btime GrothAlgebra.multi_scalar_mul($bases, $scalars);
    end
end

function bench_batch_norm()
    println("\n== Batch normalization: batch_to_affine! vs per-point to_affine (G1) ==")
    g = g1_generator()
    R = MersenneTwister(11)
    for N in (32, 128, 512)
        scalars = genscalars(R, N)
        pts = [scalar_mul(g, s) for s in scalars]
        pts2 = copy(pts)
        @show N
        @btime GrothCurves.batch_to_affine!($pts);
        @btime begin
            for i in eachindex($pts2)
                x,y = to_affine($pts2[i])
                $pts2[i] = G1Point(x,y,one(BN254Field))
            end
        end
    end

    println("\n== Batch normalization: batch_to_affine! vs per-point to_affine (G2) ==")
    g = g2_generator()
    R = MersenneTwister(13)
    for N in (32, 128, 512)
        scalars = genscalars(R, N)
        pts = [scalar_mul(g, s) for s in scalars]
        pts2 = copy(pts)
        @show N
        @btime GrothCurves.batch_to_affine!($pts);
        @btime begin
            for i in eachindex($pts2)
                x,y = to_affine($pts2[i])
                $pts2[i] = G2Point(x,y,one(Fp2Element))
            end
        end
    end
end

function main()
    println("GrothBenchmarks — Julia $(VERSION)")
    bench_fixed_base()
    bench_variable_msm()
    bench_batch_norm()
end

main()


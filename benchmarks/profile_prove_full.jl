#!/usr/bin/env julia

using Dates
using JSON
using Profile
using Random
using Sockets
using GrothAlgebra
using GrothCurves
using GrothProofs

include("prove_full_common.jl")

const ARTIFACTS_DIR = joinpath(@__DIR__, "artifacts")
const RUN_ID_RE = r"^\d{4}-\d{2}-\d{2}_\d{6}$"

is_run_id(s::AbstractString) = occursin(RUN_ID_RE, s)

function usage()
    println("Usage: julia --project=. benchmarks/profile_prove_full.jl [--run-id=<id>] [--fixture=<name>] [--repetitions=<n>]")
    println("If --run-id is omitted, a fresh artifact directory is created under benchmarks/artifacts/.")
end

function try_readchomp(cmd::Cmd)
    try
        return strip(readchomp(cmd))
    catch
        return "unknown"
    end
end

function profile_metadata(run_id::String, repetitions::Int, fixture_names::Vector{String})
    return Dict{String,Any}(
        "timestamp_utc" => string(Dates.now(Dates.UTC)),
        "run_id" => run_id,
        "repetitions" => repetitions,
        "fixtures" => fixture_names,
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

function parse_args()
    run_id = nothing
    fixture_name = nothing
    repetitions = 50
    for arg in ARGS
        if startswith(arg, "--run-id=")
            run_id = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--fixture=")
            fixture_name = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--repetitions=")
            repetitions = parse(Int, split(arg, "=", limit=2)[2])
        elseif arg in ("-h", "--help")
            usage()
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    repetitions > 0 || error("repetitions must be positive")
    return run_id, fixture_name, repetitions
end

function resolve_artifact_dir(run_id_arg)
    if run_id_arg === nothing
        run_id = Dates.format(Dates.now(Dates.UTC), "yyyy-mm-dd_HHMMSS")
        artifact_dir = joinpath(ARTIFACTS_DIR, run_id)
        mkpath(artifact_dir)
        return run_id, artifact_dir
    end

    run_id = String(run_id_arg)
    is_run_id(run_id) || error("Invalid run id: $run_id")
    artifact_dir = joinpath(ARTIFACTS_DIR, run_id)
    mkpath(artifact_dir)
    return run_id, artifact_dir
end

function write_profile(io, workload_name::String, fixture, repetitions::Int)
    println(io, "# workload: ", workload_name)
    println(io, "# fixture: ", fixture.name)
    println(io, "# description: ", fixture.description)
    println(io, "# repetitions: ", repetitions)
    println(io, "# constraints: ", fixture.r1cs.num_constraints)
    println(io, "# vars: ", fixture.r1cs.num_vars)
    println(io, "# public: ", fixture.r1cs.num_public)
    println(io, "# domain_size: ", fixture.qap.domain.size)
    println(io)
    println(io, "## flat")
    Profile.print(io; format=:flat, sortedby=:count, maxdepth=24)
    println(io)
    println(io, "## tree")
    Profile.print(io; format=:tree, sortedby=:count, maxdepth=24)
end

function capture_profile(out_path::String, workload_name::String, fixture, workload::Function; repetitions::Int)
    workload()
    GC.gc()
    Profile.clear()
    @profile begin
        for _ in 1:repetitions
            workload()
        end
    end
    open(out_path, "w") do io
        write_profile(io, workload_name, fixture, repetitions)
    end
    Profile.clear()
end

function query_msm_bundle(fixture)
    scalars = witness_to_scalars(fixture.witness)
    return prove_query_accumulators(fixture.keypair.pk, scalars)
end

function main()
    run_id_arg, fixture_name, repetitions = parse_args()
    run_id, artifact_dir = resolve_artifact_dir(run_id_arg)
    profiles_dir = joinpath(artifact_dir, "profiles")
    meta_dir = joinpath(artifact_dir, "meta")
    mkpath(profiles_dir)
    mkpath(meta_dir)

    fixtures = default_prove_full_fixtures()
    if fixture_name !== nothing
        fixtures = filter(f -> f.name == fixture_name, fixtures)
        isempty(fixtures) && error("Unknown fixture '$fixture_name'")
    end

    meta = profile_metadata(run_id, repetitions, [fixture.name for fixture in fixtures])
    meta["fixture_metadata"] = Dict(fixture.name => fixture_metadata(fixture) for fixture in fixtures)
    open(joinpath(meta_dir, "prove_full_profile_meta.json"), "w") do io
        JSON.print(io, meta)
    end

    for fixture in fixtures
        println("Profiling fixture $(fixture.name)")
        capture_profile(
            joinpath(profiles_dir, "prove_full_$(fixture.name).txt"),
            "prove_full",
            fixture,
            () -> prove_full(fixture.keypair.pk, fixture.qap, fixture.witness; rng=MersenneTwister(fixture.prove_seed));
            repetitions = repetitions,
        )
        capture_profile(
            joinpath(profiles_dir, "compute_h_$(fixture.name).txt"),
            "compute_h_polynomial",
            fixture,
            () -> compute_h_polynomial(fixture.qap, fixture.witness);
            repetitions = repetitions,
        )
        capture_profile(
            joinpath(profiles_dir, "query_msms_$(fixture.name).txt"),
            "query_msms",
            fixture,
            () -> query_msm_bundle(fixture);
            repetitions = repetitions,
        )
    end

    println("Profile artifacts written to ", profiles_dir)
end

main()

#!/usr/bin/env julia

using Test

const PACKAGE_DIRS = [
    "GrothAlgebra",
    "GrothCurves",
    "GrothCrypto",
    "GrothProofs",
]

const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))

struct TestResult
    package::String
    ok::Bool
    error::Union{Nothing,String}
end

function ensure_package_load_paths!()
    for package_dir in PACKAGE_DIRS
        package_path = joinpath(ROOT_DIR, package_dir)
        package_path in LOAD_PATH || push!(LOAD_PATH, package_path)
    end
end

function make_test_module(package_dir::String)::Module
    module_name = Symbol("CITests_", package_dir)
    mod = Module(module_name)
    Core.eval(mod, :(using Test))
    Core.eval(mod, :(include(path::AbstractString) = Base.include($mod, path)))
    return mod
end

function run_package_tests(package_dir::String)::TestResult
    test_path = joinpath(ROOT_DIR, package_dir, "test", "runtests.jl")
    println("==> Testing ", package_dir)
    try
        Base.include(make_test_module(package_dir), test_path)
        println("    PASS: ", package_dir)
        return TestResult(package_dir, true, nothing)
    catch err
        msg = sprint(showerror, err, catch_backtrace())
        println("    FAIL: ", package_dir)
        return TestResult(package_dir, false, msg)
    end
end

function main()
    ensure_package_load_paths!()

    results = TestResult[]
    for package_dir in PACKAGE_DIRS
        push!(results, run_package_tests(package_dir))
        println()
    end

    passed = count(r -> r.ok, results)
    failed = length(results) - passed

    println("==== CI Test Summary ====")
    for result in results
        status = result.ok ? "PASS" : "FAIL"
        println(status, "  ", result.package)
    end
    println("=========================")
    println("Passed: ", passed, " / ", length(results))
    println("Failed: ", failed)

    if failed > 0
        println("\nFailure details:")
        for result in results
            if !result.ok
                println("---- ", result.package, " ----")
                println(result.error)
                println()
            end
        end
        exit(1)
    end
end

main()

#!/usr/bin/env julia

using Pkg

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

function run_package_tests(package_dir::String)::TestResult
    package_path = joinpath(ROOT_DIR, package_dir)
    println("==> Testing ", package_dir)
    try
        Pkg.activate(package_path)
        Pkg.test()
        println("    PASS: ", package_dir)
        return TestResult(package_dir, true, nothing)
    catch err
        msg = sprint(showerror, err, catch_backtrace())
        println("    FAIL: ", package_dir)
        return TestResult(package_dir, false, msg)
    end
end

function main()
    results = TestResult[]
    for package_dir in PACKAGE_DIRS
        push!(results, run_package_tests(package_dir))
        println()
    end

    passed = count(r -> r.ok, results)
    failed = length(results) - passed

    println("==== Test Summary ====")
    for result in results
        status = result.ok ? "PASS" : "FAIL"
        println(status, "  ", result.package)
    end
    println("======================")
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

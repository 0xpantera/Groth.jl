using GrothExamples
using Test

@testset "GrothExamples.jl" begin
    @testset "R1CS QAP Example" begin
        # Test that the example runs without error
        @test_nowarn demonstrate_r1cs_qap()
    end
end
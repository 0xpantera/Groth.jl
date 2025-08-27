using GrothCurves
using Test

@testset "GrothCurves.jl" begin
    include("test_bn254_curve.jl")
    include("test_extension_fields.jl")
    include("test_miller_loop.jl")
end

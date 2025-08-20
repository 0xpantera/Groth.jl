using GrothAlgebra
using Test

@testset "GrothAlgebra.jl" begin
    
    # Test finite fields
    include("test_finite_fields.jl")
    
    # Test group operations
    include("test_groups.jl")
    
    # Test polynomial operations  
    include("test_polynomials.jl")
    
end
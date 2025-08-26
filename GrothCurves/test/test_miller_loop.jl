using Test
using GrothCurves
using GrothAlgebra

@testset "Miller Loop Operations" begin
    
    @testset "Line coefficient construction" begin
        # Test that we can create LineCoeffs
        c0 = Fp2Element(1, 2)
        c1 = Fp2Element(3, 4)
        c2 = Fp2Element(5, 6)
        
        coeffs = LineCoeffs(c0, c1, c2)
        @test coeffs.c0 == c0
        @test coeffs.c1 == c1
        @test coeffs.c2 == c2
    end
    
    @testset "Doubling step" begin
        # Test with a simple point
        P = G1Point(1, 2)  # Point on G1
        Q = G2Point(
            Fp2Element(1, 0),
            Fp2Element(2, 0)
        )  # Point on G2
        
        # Compute doubling step
        T_doubled, line_coeffs = doubling_step(Q, P)
        
        # Check that we got a G2Point back
        @test isa(T_doubled, G2Point)
        @test isa(line_coeffs, LineCoeffs)
        
        # Test with point at infinity
        Q_inf = zero(G2Point)
        T_inf, coeffs_inf = doubling_step(Q_inf, P)
        @test iszero(T_inf)
    end
    
    @testset "Addition step" begin
        # Create test points
        P = G1Point(1, 2)
        Q1 = G2Point(Fp2Element(1, 0), Fp2Element(2, 0))
        Q2 = G2Point(Fp2Element(3, 0), Fp2Element(4, 0))
        
        # Test addition
        T_sum, line_coeffs = addition_step(Q1, Q2, P)
        @test isa(T_sum, G2Point)
        @test isa(line_coeffs, LineCoeffs)
        
        # Test with identity
        Q_inf = zero(G2Point)
        T_id, coeffs_id = addition_step(Q1, Q_inf, P)
        @test T_id == Q1
    end
    
    @testset "Line evaluation" begin
        # Test that evaluate_line produces an Fp12Element
        coeffs = LineCoeffs(
            Fp2Element(1, 2),
            Fp2Element(3, 4),
            Fp2Element(5, 6)
        )
        
        f = evaluate_line(coeffs)
        @test isa(f, Fp12Element)
    end
    
    @testset "Miller loop basic" begin
        # Test with identity elements
        P = zero(G1Point)
        Q = zero(G2Point)
        
        f = miller_loop(P, Q)
        @test isa(f, Fp12Element)
        @test isone(f)  # Miller loop of identity should be 1
        
        # Test with generators (even if they're not the actual generators yet)
        P_gen = G1Point(1, 2)
        Q_gen = G2Point(Fp2Element(1, 0), Fp2Element(2, 0))
        
        f_gen = miller_loop(P_gen, Q_gen)
        @test isa(f_gen, Fp12Element)
        @test !iszero(f_gen)  # Should be non-zero
    end
    
    @testset "Miller loop properties" begin
        # Test that miller_loop(0, Q) = 1
        P_zero = zero(G1Point)
        Q = G2Point(Fp2Element(1, 0), Fp2Element(2, 0))
        
        f = miller_loop(P_zero, Q)
        @test isone(f)
        
        # Test that miller_loop(P, 0) = 1
        P = G1Point(1, 2)
        Q_zero = zero(G2Point)
        
        f = miller_loop(P, Q_zero)
        @test isone(f)
    end
end
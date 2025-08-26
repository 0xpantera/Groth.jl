using Test
using GrothCurves
using GrothAlgebra

# Import additional functions we need
import GrothCurves: is_on_curve, real, imag, double

@testset "Miller Loop with Real Test Vectors" begin
    
    @testset "Generator points" begin
        # Verify generators are on their respective curves
        g1 = g1_generator()
        g2 = g2_generator()
        
        @test is_on_curve(g1)
        @test is_on_curve(g2)
        
        # Check coordinates match expected values
        g1_aff = to_affine(g1)
        @test g1_aff[1] == bn254_field(1)
        @test g1_aff[2] == bn254_field(2)
        
        g2_aff = to_affine(g2)
        x0_expected = bn254_field(parse(BigInt, "10857046999023057135944570762232829481370756359578518086990519993285655852781"))
        x1_expected = bn254_field(parse(BigInt, "11559732032986387107991004021392285783925812861821192530917403151452391805634"))
        @test real(g2_aff[1]) == x0_expected
        @test imag(g2_aff[1]) == x1_expected
    end
    
    @testset "Miller loop with generators" begin
        g1 = g1_generator()
        g2 = g2_generator()
        
        # Miller loop of generators should be non-trivial
        f = miller_loop(g1, g2)
        @test !iszero(f)
        @test !isone(f)
        
        # Test identity properties
        @test isone(miller_loop(zero(G1Point), g2))
        @test isone(miller_loop(g1, zero(G2Point)))
    end
    
    @testset "Miller loop linearity (without final exp)" begin
        # Note: Miller loop without final exponentiation doesn't have
        # full bilinearity, but has some related properties
        
        g1 = g1_generator()
        g2 = g2_generator()
        
        # Test with small scalars
        # [2]g1
        g1_2 = double(g1)
        # [2]g2
        g2_2 = double(g2)
        
        f1 = miller_loop(g1, g2)
        f2 = miller_loop(g1_2, g2)
        f3 = miller_loop(g1, g2_2)
        f4 = miller_loop(g1_2, g2_2)
        
        # These should all be different (without final exp)
        @test f1 != f2
        @test f1 != f3
        @test f2 != f4
        
        # But none should be zero or one
        @test !iszero(f1) && !isone(f1)
        @test !iszero(f2) && !isone(f2)
        @test !iszero(f3) && !isone(f3)
        @test !iszero(f4) && !isone(f4)
    end
    
    @testset "Line function consistency" begin
        g1 = g1_generator()
        g2 = g2_generator()
        
        # Test doubling step
        T_doubled, coeffs_double = doubling_step(g2, g1)
        @test !iszero(T_doubled)
        @test T_doubled != g2  # Should be different from original
        
        # Verify 2*g2 via doubling
        g2_doubled_direct = double(g2)
        @test T_doubled == g2_doubled_direct
        
        # Test addition step
        g2_2 = double(g2)
        T_sum, coeffs_add = addition_step(g2, g2_2, g1)
        @test !iszero(T_sum)
        @test T_sum != g2
        @test T_sum != g2_2
        
        # Line evaluations should be non-trivial
        f_double = evaluate_line(coeffs_double)
        f_add = evaluate_line(coeffs_add)
        @test !iszero(f_double)
        @test !iszero(f_add)
        @test f_double != f_add
    end
    
    @testset "Miller loop edge cases" begin
        g1 = g1_generator()
        g2 = g2_generator()
        
        # Test with negated points
        neg_g1 = -g1
        neg_g2 = -g2
        
        f1 = miller_loop(g1, g2)
        f2 = miller_loop(neg_g1, g2)
        f3 = miller_loop(g1, neg_g2)
        f4 = miller_loop(neg_g1, neg_g2)
        
        # These should all be different without final exp
        @test f1 != f2
        @test f1 != f3
        @test f1 != f4
        
        # But none should be trivial
        for f in [f1, f2, f3, f4]
            @test !iszero(f) && !isone(f)
        end
    end
    
    @testset "Ate loop parameter verification" begin
        # Verify our loop count is correct
        u = BigInt(4965661367192848881)
        ate_loop = 6 * u + 2
        @test ate_loop == BigInt(29793968203157093288)
        
        # Check it matches our constant
        @test GrothCurves.ATE_LOOP_COUNT == ate_loop
        @test GrothCurves.BN254_U == u
        @test GrothCurves.BN254_U_IS_NEGATIVE == false
    end
end
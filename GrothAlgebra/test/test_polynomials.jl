"""
Tests for Polynomial operations
"""

using GrothAlgebra
using Test

@testset "Polynomial Tests" begin
    # Use BN254Field for testing
    F = BN254Field
    
    @testset "Polynomial Construction" begin
        # Test basic construction
        coeffs = [bn254_field(1), bn254_field(2), bn254_field(3)]
        p = Polynomial(coeffs)
        @test length(p.coeffs) == 3
        
        # Test construction from integers
        q = Polynomial([1, 2, 3], F)
        @test q.coeffs == coeffs
        
        # Test normalization (removing leading zeros)
        coeffs_with_zeros = [bn254_field(1), bn254_field(2), bn254_field(0)]
        r = Polynomial(coeffs_with_zeros)
        @test length(r.coeffs) == 2
        
        # Test zero polynomial
        zero_poly = Polynomial([bn254_field(0)])
        @test is_zero(zero_poly)
    end
    
    @testset "Polynomial Properties" begin
        # Test degree
        p = Polynomial([bn254_field(1), bn254_field(2), bn254_field(3)])
        @test degree(p) == 2
        
        zero_poly = zero(Polynomial{F})
        @test degree(zero_poly) == -1
        
        # Test leading coefficient
        @test leading_coefficient(p) == bn254_field(3)
        
        # Test predicates
        @test is_zero(zero_poly)
        @test !is_zero(p)
        
        constant_poly = Polynomial([bn254_field(5)])
        @test is_constant(constant_poly)
        @test !is_constant(p)
        
        monic_poly = Polynomial([bn254_field(2), bn254_field(1)])
        @test is_monic(monic_poly)
        @test !is_monic(p)
    end
    
    @testset "Polynomial Special Elements" begin
        # Test zero and one polynomials
        zero_poly = zero(Polynomial{F})
        one_poly = one(Polynomial{F})
        
        @test is_zero(zero_poly)
        @test degree(zero_poly) == -1
        @test isone(leading_coefficient(one_poly))
        @test degree(one_poly) == 0
        
        # Test constant polynomial
        c = bn254_field(7)
        const_poly = constant_polynomial(c)
        @test degree(const_poly) == 0
        @test leading_coefficient(const_poly) == c
        
        # Test monomial
        mono = monomial(3, F)
        @test degree(mono) == 3
        @test leading_coefficient(mono) == one(F)
        
        @test_throws ArgumentError monomial(-1, F)
    end
    
    @testset "Polynomial Equality" begin
        p1 = Polynomial([bn254_field(1), bn254_field(2)])
        p2 = Polynomial([bn254_field(1), bn254_field(2)])
        p3 = Polynomial([bn254_field(1), bn254_field(3)])
        
        @test p1 == p2
        @test p1 != p3
        @test isequal(p1, p2)
        @test !isequal(p1, p3)
    end
    
    @testset "Polynomial Arithmetic" begin
        p1 = Polynomial([bn254_field(1), bn254_field(2)])      # 1 + 2x
        p2 = Polynomial([bn254_field(3), bn254_field(4)])      # 3 + 4x
        
        # Test addition
        sum_poly = p1 + p2
        expected_sum = Polynomial([bn254_field(4), bn254_field(6)])  # 4 + 6x
        @test sum_poly == expected_sum
        
        # Test subtraction
        diff_poly = p1 - p2
        expected_diff = Polynomial([bn254_field(-2), bn254_field(-2)])  # -2 - 2x
        @test diff_poly == expected_diff
        
        # Test negation
        neg_poly = -p1
        expected_neg = Polynomial([bn254_field(-1), bn254_field(-2)])
        @test neg_poly == expected_neg
        
        # Test multiplication
        prod_poly = p1 * p2
        # (1 + 2x)(3 + 4x) = 3 + 4x + 6x + 8x^2 = 3 + 10x + 8x^2
        expected_prod = Polynomial([bn254_field(3), bn254_field(10), bn254_field(8)])
        @test prod_poly == expected_prod
        
        # Test scalar multiplication
        scalar = bn254_field(3)
        scalar_prod = p1 * scalar
        expected_scalar = Polynomial([bn254_field(3), bn254_field(6)])
        @test scalar_prod == expected_scalar
        @test scalar * p1 == scalar_prod
    end
    
    @testset "Polynomial Exponentiation" begin
        p = Polynomial([bn254_field(1), bn254_field(1)])  # 1 + x
        
        # Test basic exponentiation
        @test p^0 == one(Polynomial{F})
        @test p^1 == p
        
        p_squared = p^2
        # (1 + x)^2 = 1 + 2x + x^2
        expected = Polynomial([bn254_field(1), bn254_field(2), bn254_field(1)])
        @test p_squared == expected
        
        # Test invalid negative exponent
        @test_throws MethodError p^(-1)
    end
    
    @testset "Polynomial Evaluation" begin
        p = Polynomial([bn254_field(1), bn254_field(2), bn254_field(3)])  # 1 + 2x + 3x^2
        
        # Test evaluation at specific points
        @test evaluate(p, bn254_field(0)) == bn254_field(1)
        @test evaluate(p, bn254_field(1)) == bn254_field(6)  # 1 + 2 + 3
        @test evaluate(p, bn254_field(2)) == bn254_field(17)  # 1 + 4 + 12
        
        # Test callable syntax
        @test p(bn254_field(1)) == bn254_field(6)
        
        # Test evaluation at multiple points
        points = [bn254_field(0), bn254_field(1), bn254_field(2)]
        values = evaluate(p, points)
        expected_values = [bn254_field(1), bn254_field(6), bn254_field(17)]
        @test values == expected_values
    end
    
    @testset "Polynomial Interpolation" begin
        # Test simple linear interpolation
        points = [bn254_field(0), bn254_field(1)]
        values = [bn254_field(1), bn254_field(3)]
        
        p = interpolate(points, values)
        
        # Verify that the polynomial passes through the points
        @test evaluate(p, points[1]) == values[1]
        @test evaluate(p, points[2]) == values[2]
        
        # Test with single point
        single_p = interpolate([bn254_field(1)], [bn254_field(5)])
        @test is_constant(single_p)
        @test evaluate(single_p, bn254_field(0)) == bn254_field(5)
        
        # Test error cases
        @test_throws ArgumentError interpolate(F[], F[])
        @test_throws ArgumentError interpolate([bn254_field(1)], [bn254_field(1), bn254_field(2)])
    end
    
    @testset "Polynomial Derivative" begin
        # Test derivative of 3 + 2x + x^2 should be 2 + 2x
        p = Polynomial([bn254_field(3), bn254_field(2), bn254_field(1)])
        dp = derivative(p)
        expected = Polynomial([bn254_field(2), bn254_field(2)])
        @test dp == expected
        
        # Test derivative of constant
        const_p = Polynomial([bn254_field(5)])
        @test is_zero(derivative(const_p))
        
        # Test derivative of linear
        linear_p = Polynomial([bn254_field(1), bn254_field(2)])
        linear_dp = derivative(linear_p)
        @test linear_dp == Polynomial([bn254_field(2)])
    end
end
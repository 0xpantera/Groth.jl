"""
Tests for Group operations
"""

using GrothAlgebra
using Test

@testset "Group Tests" begin
    # Create a simple test group element type for testing
    struct TestCurve <: AbstractCurve end
    
    struct TestPoint <: GroupElem{TestCurve}
        x::Int
        y::Int
    end

    struct TestCurveBig <: AbstractCurve end

    struct TestPointBig <: GroupElem{TestCurveBig}
        x::BigInt
        y::BigInt
    end
    
    # Implement required group operations for testing
    Base.zero(::Type{TestPoint}) = TestPoint(0, 0)
    Base.iszero(p::TestPoint) = p.x == 0 && p.y == 0
    Base.:+(p::TestPoint, q::TestPoint) = TestPoint(p.x + q.x, p.y + q.y)
    Base.:-(p::TestPoint) = TestPoint(-p.x, -p.y)
    Base.:(==)(p::TestPoint, q::TestPoint) = p.x == q.x && p.y == q.y
    Base.zero(::Type{TestPointBig}) = TestPointBig(0, 0)
    Base.iszero(p::TestPointBig) = p.x == 0 && p.y == 0
    Base.:+(p::TestPointBig, q::TestPointBig) = TestPointBig(p.x + q.x, p.y + q.y)
    Base.:-(p::TestPointBig) = TestPointBig(-p.x, -p.y)
    Base.:(==)(p::TestPointBig, q::TestPointBig) = p.x == q.x && p.y == q.y

    function naive_msm(points, scalars)
        acc = zero(points[1])
        @inbounds for i in eachindex(points, scalars)
            acc += scalar_mul(points[i], scalars[i])
        end
        return acc
    end
    
    @testset "Group Element Basic Operations" begin
        P = TestPoint(1, 2)
        Q = TestPoint(3, 4)
        O = zero(TestPoint)
        
        # Test zero element
        @test iszero(O)
        @test !iszero(P)
        @test zero(P) == O
        
        # Test addition
        R = P + Q
        @test R == TestPoint(4, 6)
        
        # Test subtraction
        @test P - Q == TestPoint(-2, -2)
        
        # Test negation
        @test -P == TestPoint(-1, -2)
        
        # Test identity operations
        @test group_identity(P) == O
        
        # Test inverse (same as negation for additive groups)
        @test inv(P) == -P
    end
    
    @testset "Scalar Multiplication" begin
        P = TestPoint(2, 3)
        O = zero(TestPoint)
        
        # Test basic scalar multiplication
        @test scalar_mul(P, 0) == O
        @test scalar_mul(P, 1) == P
        @test scalar_mul(P, 2) == TestPoint(4, 6)
        @test scalar_mul(P, 3) == TestPoint(6, 9)
        
        # Test negative scalars
        @test scalar_mul(P, -1) == TestPoint(-2, -3)
        @test scalar_mul(P, -2) == TestPoint(-4, -6)
        
        # Test operator overloads
        @test 3 * P == TestPoint(6, 9)
        @test P * 3 == TestPoint(6, 9)
    end
    
    @testset "Multi-scalar Multiplication" begin
        P1 = TestPoint(1, 2)
        P2 = TestPoint(3, 4)
        P3 = TestPoint(5, 6)
        
        points = [P1, P2, P3]
        scalars = [2, 3, 1]
        
        # Test multi-scalar multiplication
        result = multi_scalar_mul(points, scalars)
        expected = 2 * P1 + 3 * P2 + 1 * P3
        @test result == expected
        
        # Test with empty vectors
        @test_throws ArgumentError multi_scalar_mul(TestPoint[], Int[])
        
        # Test with mismatched lengths
        @test_throws ArgumentError multi_scalar_mul([P1], [1, 2])
        
        # Test with single point
        @test multi_scalar_mul([P1], [5]) == 5 * P1

        # Zero and negative scalars should preserve the same algebraic result
        @test multi_scalar_mul(points, [0, -2, 5]) == naive_msm(points, [0, -2, 5])
    end

    @testset "Multi-scalar Multiplication BigInt scalars" begin
        P1 = TestPointBig(1, 2)
        P2 = TestPointBig(3, 4)
        P3 = TestPointBig(5, 6)
        points = [P1, P2, P3]
        scalars = [BigInt(1) << 60, BigInt(1) << 55, -(BigInt(1) << 54)]

        result = multi_scalar_mul(points, scalars)
        expected = scalar_mul(P1, scalars[1]) + scalar_mul(P2, scalars[2]) + scalar_mul(P3, scalars[3])
        @test result == expected
    end

    @testset "Multi-scalar Multiplication matches naive sums across sizes" begin
        for count in (2, 5, 8, 16, 33)
            points = [TestPoint(mod(7 * i, 19) - 9, mod(11 * i, 23) - 11) for i in 1:count]
            scalars = [mod(13 * i, 41) - 20 for i in 1:count]
            @test multi_scalar_mul(points, scalars) == naive_msm(points, scalars)
        end

        big_points = [TestPointBig(BigInt(i), BigInt(2i)) for i in 1:12]
        big_scalars = [BigInt((-1)^i) * (BigInt(1) << (50 + i)) for i in 1:12]
        @test multi_scalar_mul(big_points, big_scalars) == naive_msm(big_points, big_scalars)
    end

    @testset "Forced Pippenger windows preserve MSM semantics" begin
        points = [TestPoint(mod(5 * i, 17) - 8, mod(9 * i, 29) - 14) for i in 1:12]
        scalars = [mod(11 * i, 53) - 26 for i in 1:12]
        expected = naive_msm(points, scalars)
        max_bits = maximum(GrothAlgebra._bit_length, scalars)

        for w in 2:6
            @test GrothAlgebra._multi_scalar_mul_pippenger(points, scalars, max_bits, w) == expected
        end

        big_points = [TestPointBig(BigInt(i), BigInt(3i)) for i in 1:10]
        big_scalars = [BigInt((-1)^i) * (BigInt(1) << (40 + i)) for i in 1:10]
        big_expected = naive_msm(big_points, big_scalars)
        big_bits = maximum(GrothAlgebra._bit_length, big_scalars)

        for w in 2:6
            @test GrothAlgebra._multi_scalar_mul_pippenger(big_points, big_scalars, big_bits, w) == big_expected
        end
    end
    
    @testset "w-NAF Encoding" begin
        # Test basic w-NAF encoding
        @test wnaf_encode(0) == [0]
        @test wnaf_encode(1) == [1]
        @test wnaf_encode(2) == [0, 1]
        @test wnaf_encode(3) == [3]
        
        # Test window sizes
        naf4 = wnaf_encode(15, 4)
        @test length(naf4) >= 1
        
        # Test negative numbers
        naf_neg = wnaf_encode(-5)
        naf_pos = wnaf_encode(5)
        @test naf_neg == -naf_pos
        
        # Test invalid window size
        @test_throws ArgumentError wnaf_encode(5, 1)
    end
    
    @testset "w-NAF Scalar Multiplication" begin
        P = TestPoint(2, 3)
        
        # Test basic w-NAF scalar multiplication
        @test scalar_mul_wnaf(P, 0) == zero(P)
        @test scalar_mul_wnaf(P, 1) == P
        @test scalar_mul_wnaf(P, -1) == -P
        
        # Test that w-NAF gives same result as regular scalar multiplication
        for k in [5, 10, 15, 23, 42]
            @test scalar_mul_wnaf(P, k) == scalar_mul(P, k)
            @test scalar_mul_wnaf(P, -k) == scalar_mul(P, -k)
        end
    end
    
    @testset "Group Utility Functions" begin
        P = TestPoint(1, 2)
        O = zero(TestPoint)
        
        # Test doubling and tripling
        @test double(P) == P + P
        @test triple(P) == P + P + P
        
        # Test isfinite (non-zero check for our simple group)
        @test isfinite(P)
        @test !isfinite(O)
        
        # Test order computation for simple cases
        @test order(O) == 1
    end
end

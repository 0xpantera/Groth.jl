"""
Tests for the FiniteFields implementation
"""

using GrothAlgebra
using Test

@testset "FiniteFields Tests" begin

    @testset "BN254Fq Tests" begin
        @testset "Field Element Construction" begin
            # Test basic construction
            a = bn254_fq(7)
            @test convert(BigInt, a) == BigInt(7)
            @test prime(BN254Fq) == parse(BigInt, "21888242871839275222246405745257275088696311157297823662689037894645226208583")

            # Test construction from BigInt
            b = BN254Fq(BigInt(12))
            @test convert(BigInt, b) == BigInt(12)

            # Test normalization
            p = prime(BN254Fq)
            large_val = p + BigInt(5)
            c = BN254Fq(large_val)
            @test convert(BigInt, c) == BigInt(5)

            # Test negative values (should normalize correctly)
            d = BN254Fq(-5)
            expected = p - BigInt(5)
            @test convert(BigInt, d) == expected
        end

        @testset "Field Element Equality" begin
            a = bn254_fq(7)
            b = bn254_fq(7)
            c = bn254_fq(6)

            @test a == b
            @test a != c
            @test isequal(a, b)
            @test !isequal(a, c)
        end

        @testset "Zero and One Elements" begin
            zero_elem = zero(BN254Fq)
            one_elem = one(BN254Fq)

            @test iszero(zero_elem)
            @test !iszero(one_elem)
            @test isone(one_elem)
            @test !isone(zero_elem)

            @test is_zero(zero_elem)
            @test is_one(one_elem)
            @test is_unity(one_elem)

            # Test that zero and one work with instances
            a = bn254_fq(5)
            @test zero(a) == zero_elem
            @test one(a) == one_elem
        end

        @testset "Field Addition" begin
            # Test basic addition
            a = bn254_fq(7)
            b = bn254_fq(12)
            c = bn254_fq(19)

            @test a + b == c

            # Test addition with wrapping
            p = prime(BN254Fq)
            large_a = bn254_fq(p - BigInt(1))
            small_b = bn254_fq(2)
            expected = bn254_fq(1)

            @test large_a + small_b == expected

            # Test commutativity
            @test a + b == b + a

            # Test associativity
            d = bn254_fq(3)
            @test (a + b) + d == a + (b + d)

            # Test identity
            @test a + zero(a) == a
        end

        @testset "Field Subtraction" begin
            # Test basic subtraction
            a = bn254_fq(19)
            b = bn254_fq(12)
            c = bn254_fq(7)

            @test a - b == c

            # Test subtraction with wrapping
            small_a = bn254_fq(5)
            large_b = bn254_fq(10)
            p = prime(BN254Fq)
            expected = bn254_fq(p - BigInt(5))

            @test small_a - large_b == expected

            # Test negation
            @test -a == zero(a) - a

            # Test identity
            @test a - zero(a) == a
            @test a - a == zero(a)
        end

        @testset "Field Multiplication" begin
            # Test basic multiplication
            a = bn254_fq(3)
            b = bn254_fq(4)
            c = bn254_fq(12)

            @test a * b == c

            # Test commutativity
            @test a * b == b * a

            # Test associativity
            d = bn254_fq(5)
            @test (a * b) * d == a * (b * d)

            # Test identity
            @test a * one(a) == a

            # Test zero
            @test a * zero(a) == zero(a)

            # Test scalar multiplication
            n = 7
            expected = bn254_fq(21)  # 3 * 7
            @test a * n == expected
            @test n * a == expected
        end

        @testset "Field Division and Inverse" begin
            a = bn254_fq(15)
            b = bn254_fq(3)
            c = bn254_fq(5)

            # Test division
            @test a / b == c

            # Test inverse
            b_inv = inv(b)
            @test b * b_inv == one(b)

            # Test division as multiplication by inverse
            @test a / b == a * inv(b)

            # Test that division by zero throws
            @test_throws DivideError inv(zero(a))
        end

        @testset "Field Exponentiation" begin
            a = bn254_fq(3)

            # Test basic exponentiation
            @test a^0 == one(a)
            @test a^1 == a
            @test a^2 == a * a
            @test a^3 == a * a * a

            # Test negative exponents
            @test a^(-1) == inv(a)
            @test a^(-2) == inv(a * a)
        end

        @testset "Field Display" begin
            a = bn254_fq(42)
            str = string(a)
            @test occursin("BN254", str)
            @test occursin("42", str)
        end
    end

    @testset "Secp256k1Field Tests" begin
        @testset "Basic Operations" begin
            a = secp256k1_field(5)
            b = secp256k1_field(7)

            @test a + b == secp256k1_field(12)
            @test a * b == secp256k1_field(35)
            @test b - a == secp256k1_field(2)

            # Test field is different from BN254
            @test prime(Secp256k1Field) != prime(BN254Fq)
        end
    end

    @testset "GaloisField Tests" begin
        GF7 = galois_field(7)
        one_half = GF7(1) / GF7(2)
        one_third = GF7(1) / GF7(3)
        five_over_six = GF7(5) / GF7(6)

        @test one_half + one_third == five_over_six
        @test convert(BigInt, GF7(9)) == BigInt(2)
        @test string(GF7(4)) == "GF(7)(4)"

        @test_throws ArgumentError galois_field(1)
        @test_throws ArgumentError galois_field(9)
    end

    @testset "Field Axioms" begin
        # Test field axioms with random values
        for _ in 1:20
            a = bn254_fq(rand(1:1000))
            b = bn254_fq(rand(1:1000))
            c = bn254_fq(rand(1:1000))

            # Commutativity
            @test a + b == b + a
            @test a * b == b * a

            # Associativity
            @test (a + b) + c == a + (b + c)
            @test (a * b) * c == a * (b * c)

            # Distributivity
            @test a * (b + c) == a * b + a * c

            # Identity
            @test a + zero(a) == a
            @test a * one(a) == a

            # Inverse
            if !iszero(b)
                @test b * inv(b) == one(b)
                @test a / b == a * inv(b)
            end
        end
    end
end

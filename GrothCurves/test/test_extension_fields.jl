using Test
using GrothCurves
using GrothAlgebra

@testset "Extension Field Operations" begin
    
    @testset "Fp6 Operations" begin
        @testset "Basic operations" begin
            # Zero and one
            @test iszero(zero(Fp6Element))
            @test isone(one(Fp6Element))
            @test !iszero(one(Fp6Element))
            @test !isone(zero(Fp6Element))
            
            # Construction
            a = Fp6Element(2, 3, 4)
            b = Fp6Element(5, 7, 11)
            
            # Addition
            c = a + b
            @test c[1] == Fp2Element(7)
            @test c[2] == Fp2Element(10)
            @test c[3] == Fp2Element(15)
            
            # Subtraction
            d = b - a
            @test d[1] == Fp2Element(3)
            @test d[2] == Fp2Element(4)
            @test d[3] == Fp2Element(7)
            
            # Negation
            neg_a = -a
            @test a + neg_a == zero(Fp6Element)
        end
        
        @testset "Multiplication" begin
            # Test v³ = 9 + u reduction
            v = Fp6Element(0, 1, 0)  # v
            v2 = Fp6Element(0, 0, 1)  # v²
            v3 = v * v2  # Should equal 9 + u
            
            # v³ should equal ξ = 9 + u
            ξ = Fp2Element(9, 1)
            @test v3[1] == ξ
            @test iszero(v3[2])
            @test iszero(v3[3])
            
            # Identity
            a = Fp6Element(2, 3, 4)
            @test a * one(Fp6Element) == a
            @test one(Fp6Element) * a == a
            
            # Zero
            @test a * zero(Fp6Element) == zero(Fp6Element)
        end
        
        @testset "Squaring" begin
            a = Fp6Element(2, 3, 4)
            @test square(a) == a * a
            @test a^Val(2) == a * a
        end
        
        @testset "Inverse" begin
            a = Fp6Element(2, 3, 4)
            a_inv = inv(a)
            @test a * a_inv == one(Fp6Element)
            @test a_inv * a == one(Fp6Element)
            
            # Test division
            b = Fp6Element(5, 7, 11)
            c = b / a
            @test c * a == b
            
            # Zero has no inverse
            @test_throws DivideError inv(zero(Fp6Element))
        end
        
        @testset "Exponentiation" begin
            a = Fp6Element(2, 3, 4)
            
            @test a^0 == one(Fp6Element)
            @test a^1 == a
            @test a^2 == a * a
            @test a^3 == a * a * a
            @test a^(-1) == inv(a)
        end
    end
    
    @testset "Fp12 Operations" begin
        @testset "Basic operations" begin
            # Zero and one
            @test iszero(zero(Fp12Element))
            @test isone(one(Fp12Element))
            
            # Construction
            a0 = Fp6Element(1, 2, 3)
            a1 = Fp6Element(4, 5, 6)
            a = Fp12Element(a0, a1)
            
            b0 = Fp6Element(7, 8, 9)
            b1 = Fp6Element(10, 11, 12)
            b = Fp12Element(b0, b1)
            
            # Addition
            c = a + b
            @test c[1] == a0 + b0
            @test c[2] == a1 + b1
            
            # Subtraction
            d = b - a
            @test d[1] == b0 - a0
            @test d[2] == b1 - a1
            
            # Negation
            neg_a = -a
            @test a + neg_a == zero(Fp12Element)
        end
        
        @testset "Multiplication" begin
            # Test w² = v reduction
            w = Fp12Element(zero(Fp6Element), one(Fp6Element))  # w
            w2 = w * w  # Should equal v
            
            # w² should equal v = (0, 1, 0) in Fp6
            v = Fp6Element(0, 1, 0)
            @test w2[1] == v
            @test iszero(w2[2])
            
            # Identity
            a0 = Fp6Element(1, 2, 3)
            a1 = Fp6Element(4, 5, 6)
            a = Fp12Element(a0, a1)
            
            @test a * one(Fp12Element) == a
            @test one(Fp12Element) * a == a
            
            # Zero
            @test a * zero(Fp12Element) == zero(Fp12Element)
        end
        
        @testset "Squaring" begin
            a0 = Fp6Element(1, 2, 3)
            a1 = Fp6Element(4, 5, 6)
            a = Fp12Element(a0, a1)
            
            @test square(a) == a * a
            @test a^Val(2) == a * a
        end
        
        @testset "Conjugate" begin
            a0 = Fp6Element(1, 2, 3)
            a1 = Fp6Element(4, 5, 6)
            a = Fp12Element(a0, a1)
            
            conj_a = conjugate(a)
            @test conj_a[1] == a0
            @test conj_a[2] == -a1
            
            # Conjugate of conjugate is original
            @test conjugate(conjugate(a)) == a
        end
        
        @testset "Inverse" begin
            a0 = Fp6Element(1, 2, 3)
            a1 = Fp6Element(4, 5, 6)
            a = Fp12Element(a0, a1)
            
            a_inv = inv(a)
            @test a * a_inv == one(Fp12Element)
            @test a_inv * a == one(Fp12Element)
            
            # Test division
            b0 = Fp6Element(7, 8, 9)
            b1 = Fp6Element(10, 11, 12)
            b = Fp12Element(b0, b1)
            
            c = b / a
            @test c * a == b
            
            # Zero has no inverse
            @test_throws DivideError inv(zero(Fp12Element))
        end
        
        @testset "Exponentiation" begin
            a0 = Fp6Element(1, 2, 3)
            a1 = Fp6Element(4, 5, 6)
            a = Fp12Element(a0, a1)
            
            @test a^0 == one(Fp12Element)
            @test a^1 == a
            @test a^2 == a * a
            @test a^3 == a * a * a
            @test a^(-1) == inv(a)
        end
    end
    
    @testset "Field Tower Consistency" begin
        # Test that embeddings work correctly
        
        # Fp2 embeds in Fp6
        a = Fp2Element(3, 5)
        a_in_fp6 = Fp6Element(a, zero(Fp2Element), zero(Fp2Element))
        b = Fp2Element(7, 11)
        b_in_fp6 = Fp6Element(b, zero(Fp2Element), zero(Fp2Element))
        
        @test (a_in_fp6 * b_in_fp6)[1] == a * b
        @test iszero((a_in_fp6 * b_in_fp6)[2])
        @test iszero((a_in_fp6 * b_in_fp6)[3])
        
        # Fp6 embeds in Fp12
        c = Fp6Element(1, 2, 3)
        c_in_fp12 = Fp12Element(c, zero(Fp6Element))
        d = Fp6Element(4, 5, 6)
        d_in_fp12 = Fp12Element(d, zero(Fp6Element))
        
        @test (c_in_fp12 * d_in_fp12)[1] == c * d
        @test iszero((c_in_fp12 * d_in_fp12)[2])
    end
end
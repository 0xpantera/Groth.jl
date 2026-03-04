using Test
using GrothCurves
using GrothAlgebra

# Explicitly import all needed functions
using GrothCurves: is_on_curve, double, to_affine, g1_generator, g2_generator, 
                   bn254_fq, x_coord, y_coord, z_coord, real, imag,
                   G1Point, G2Point, BN254Fq, Fp2Element,
                   conjugate, norm, frobenius

@testset "BN254 Curve Operations" begin
    
    @testset "G1 Point Operations" begin
        @testset "Point creation and identity" begin
            # Test identity element
            O = zero(G1Point)
            @test iszero(O)
            @test iszero(z_coord(O))
            
            # Test generator point
            g1 = g1_generator()
            @test !iszero(g1)
            @test is_on_curve(g1)
            
            # Test affine conversion
            x_aff, y_aff = to_affine(g1)
            @test x_aff == bn254_fq(1)
            @test y_aff == bn254_fq(2)
        end
        
        @testset "Point addition" begin
            g1 = g1_generator()
            O = zero(G1Point)
            
            # Identity tests
            @test g1 + O == g1
            @test O + g1 == g1
            @test O + O == O
            
            # Addition test
            p2 = g1 + g1
            @test !iszero(p2)
            @test is_on_curve(p2)
            
            # Should be same as doubling
            @test p2 == double(g1)
            
            # Commutativity
            p3 = double(g1)
            p4 = g1 + p2
            @test p4 + g1 == g1 + p4
            
            # Associativity
            @test (g1 + g1) + g1 == g1 + (g1 + g1)
        end
        
        @testset "Point doubling" begin
            g1 = g1_generator()
            
            # Double of identity is identity
            @test double(zero(G1Point)) == zero(G1Point)
            
            # Double of generator
            p2 = double(g1)
            @test !iszero(p2)
            @test is_on_curve(p2)
            
            # Double should equal addition with self
            @test p2 == g1 + g1
            
            # Test specific coordinates (known values for BN254)
            x2, y2 = to_affine(p2)
            # These are known values for 2*G1 on BN254
            expected_x = bn254_fq(parse(BigInt, "1368015179489954701390400359078579693043519447331113978918064868415326638035"))
            expected_y = bn254_fq(parse(BigInt, "9918110051302171585080402603319702774565515993150576347155970296011118125764"))
            @test x2 == expected_x
            @test y2 == expected_y
        end
        
        @testset "Point negation" begin
            g1 = g1_generator()
            
            # Negation test
            neg_g1 = -g1
            @test is_on_curve(neg_g1)
            
            # Adding point and its negation gives identity
            @test g1 + neg_g1 == zero(G1Point)
            @test neg_g1 + g1 == zero(G1Point)
            
            # Double negation
            @test -(-g1) == g1
            
            # Negation of identity
            @test -zero(G1Point) == zero(G1Point)
        end
        
        @testset "Scalar multiplication" begin
            g1 = g1_generator()
            
            # Multiplication by 0
            @test g1 * 0 == zero(G1Point)
            @test 0 * g1 == zero(G1Point)
            
            # Multiplication by 1
            @test g1 * 1 == g1
            @test 1 * g1 == g1
            
            # Multiplication by 2
            @test g1 * 2 == double(g1)
            @test 2 * g1 == double(g1)
            
            # Multiplication by 3
            @test g1 * 3 == g1 + g1 + g1
            @test 3 * g1 == g1 + g1 + g1
            
            # Multiplication by negative
            @test g1 * (-1) == -g1
            @test (-1) * g1 == -g1
            
            # Distributivity
            @test g1 * 5 == g1 * 2 + g1 * 3
            @test g1 * 10 == (g1 * 5) * 2
        end
        
        @testset "Subtraction" begin
            g1 = g1_generator()
            p2 = double(g1)
            
            # Basic subtraction
            @test p2 - g1 == g1
            @test p2 - p2 == zero(G1Point)
            @test g1 - zero(G1Point) == g1
            @test zero(G1Point) - g1 == -g1
        end
    end
    
    @testset "G2 Point Operations" begin
        @testset "Point creation and identity" begin
            # Test identity element
            O = zero(G2Point)
            @test iszero(O)
            @test iszero(z_coord(O))
            
            # Test generator point
            g2 = g2_generator()
            @test !iszero(g2)
            # Note: is_on_curve not implemented for G2 yet
        end
        
        @testset "Point addition" begin
            g2 = g2_generator()
            O = zero(G2Point)
            
            # Identity tests
            @test g2 + O == g2
            @test O + g2 == g2
            @test O + O == O
            
            # Addition test
            p2 = g2 + g2
            @test !iszero(p2)
            @test p2 == double(g2)
            
            # Commutativity
            p3 = double(g2)
            p4 = g2 + p2
            @test p4 + g2 == g2 + p4
            
            # Associativity
            @test (g2 + g2) + g2 == g2 + (g2 + g2)
        end
        
        @testset "Point doubling" begin
            g2 = g2_generator()
            
            # Double of identity is identity
            @test double(zero(G2Point)) == zero(G2Point)
            
            # Double of generator
            p2 = double(g2)
            @test !iszero(p2)
            @test p2 == g2 + g2
        end
        
        @testset "Point negation" begin
            g2 = g2_generator()
            
            # Negation test
            neg_g2 = -g2
            
            # Adding point and its negation gives identity
            @test g2 + neg_g2 == zero(G2Point)
            @test neg_g2 + g2 == zero(G2Point)
            
            # Double negation
            @test -(-g2) == g2
            
            # Negation of identity
            @test -zero(G2Point) == zero(G2Point)
        end
        
        @testset "Scalar multiplication" begin
            g2 = g2_generator()
            
            # Multiplication by 0
            @test g2 * 0 == zero(G2Point)
            
            # Multiplication by 1
            @test g2 * 1 == g2
            
            # Multiplication by 2
            @test g2 * 2 == double(g2)
            
            # Multiplication by 3
            @test g2 * 3 == g2 + g2 + g2
            
            # Multiplication by negative
            @test g2 * (-1) == -g2
            
            # Distributivity
            @test g2 * 5 == g2 * 2 + g2 * 3
        end
    end
    
    @testset "Fp2 Field Operations" begin
        @testset "Basic operations" begin
            a = Fp2Element(2, 3)
            b = Fp2Element(5, 7)
            
            # Addition
            c = a + b
            @test real(c) == bn254_fq(7)
            @test imag(c) == bn254_fq(10)
            
            # Subtraction
            d = b - a
            @test real(d) == bn254_fq(3)
            @test imag(d) == bn254_fq(4)
            
            # Multiplication (using u² = -1)
            # (2 + 3u) * (5 + 7u) = 10 + 14u + 15u + 21u²
            #                     = 10 + 29u - 21
            #                     = -11 + 29u
            e = a * b
            @test real(e) == bn254_fq(-11)
            @test imag(e) == bn254_fq(29)
            
            # Conjugate
            conj_a = conjugate(a)
            @test real(conj_a) == bn254_fq(2)
            @test imag(conj_a) == bn254_fq(-3)
            
            # Norm
            n = norm(a)
            @test n == bn254_fq(13)  # 2² + 3² = 4 + 9 = 13
        end
        
        @testset "Inverse and division" begin
            a = Fp2Element(2, 3)
            
            # Inverse
            a_inv = inv(a)
            product = a * a_inv
            @test isone(product)
            
            # Division
            b = Fp2Element(5, 7)
            c = b / a
            @test c * a == b
        end
        
        @testset "Exponentiation" begin
            a = Fp2Element(2, 3)
            
            # Square
            a2 = a^2
            @test a2 == a * a
            
            # Cube
            a3 = a^3
            @test a3 == a * a * a
            
            # Identity
            @test a^0 == one(Fp2Element)
            @test a^1 == a
            
            # Negative exponent
            @test a^(-1) == inv(a)
        end
    end
end